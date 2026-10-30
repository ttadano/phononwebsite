# On-Demand Eigenvectors

The current website format stores full eigenvectors for every q-point:

- `eigenvalues`: `Nq x 3N`
- `vectors`: `Nq x 3N x N x 3 x 2`

That is convenient for rendering, but it is also the dominant storage cost for
large PhononDB-derived datasets.

## What Phonopy Actually Uses

After tracing phonopy, the minimal path for a q-point eigensystem is:

1. build a `Phonopy` object from `phonon.yaml` + `FORCE_SETS` and optional `BORN`
2. produce force constants
3. construct a `DynamicalMatrix`
4. for a requested q-point, build `D(q)`
5. diagonalize the Hermitian complex matrix with `numpy.linalg.eigh`

For band structures, `phonopy/phonon/band_structure.py` mainly does path
bookkeeping and repeated calls to the same two core steps:

- `phonopy.harmonic.dynamical_matrix.DynamicalMatrix._run_py_dynamical_matrix`
- `numpy.linalg.eigh`

## Minimal Runtime Payload

We do not need to port phonopy's primitive-cell construction to JavaScript if
we serialize the resolved bookkeeping once during preprocessing. The minimal
runtime payload is:

- primitive masses: `masses[i]`
- compact force constants: `force_constants_compact[i, k, 3, 3]`
- dense shortest vectors in primitive reduced coordinates: `shortest_vectors`
- multiplicity/address table: `multiplicity[k, i] = [count, offset]`
- supercell-to-primitive map in primitive indexing: `s2pp_map[k]`

That is enough to rebuild the same matrix that phonopy uses:

```text
D_ij(q) = sum_k sum_l Phi(i, k) * exp(2 pi i q . r_kl) / sqrt(m_i m_j) / m_ki
```

with the same loop structure as phonopy's Python fallback.

## Why This Is Attractive

- We can keep the existing precomputed `eigenvalues` array for the dispersion plot.
- We only compute eigenvectors when the user clicks a mode.
- The browser-side code becomes a small matrix builder plus a Hermitian eigensolver.
- The heavy, calculator-specific parsing stays offline in Python.

## What Still Needs To Be Stored

The upstream PhononDB archives already contain what phonopy starts from:

- `phonon.yaml`
- `FORCE_SETS`
- optional `BORN`

From those, preprocessing can emit a compact runtime JSON block instead of full
eigenvectors. The new helper script for that is:

```bash
python -m phononweb.scripts.export_phonopy_runtime phonon.yaml FORCE_SETS --born BORN
```

## Practical Caveats

- NAC / LO-TO splitting is the one part that is not in the first minimal port.
  Non-polar materials are straightforward. Polar materials will need the NAC
  correction terms and q-direction handling added afterwards.
- We still need a JS dependency for Hermitian complex diagonalization.
- For stability, the runtime path should recompute eigenvectors only. Frequencies
  can continue coming from the precomputed website JSON.

## Suggested Integration Path

1. Extend the PhononDB preparation step to emit a `dynamical_matrix` block.
2. Teach `PhononJson` to load either precomputed `vectors` or runtime data.
3. On mode selection, compute eigenvectors for the selected q-point if missing.
4. Add NAC support as a second phase.

## Sanity Check

`python/phononweb/tests/test_runtime_dynamical_matrix.py` now verifies that the
exported minimal payload is sufficient to rebuild the same dynamical matrix as
phonopy for a sample q-point.
