"""Helpers to export the minimal data needed to rebuild phonopy dynamical matrices."""

from __future__ import annotations

import copy

import numpy as np

from phononweb.lattice import red_car


def _get_attr_or_call(obj, attr_name, getter_name):
    if hasattr(obj, getter_name):
        return getattr(obj, getter_name)()
    return getattr(obj, attr_name)


def get_primitive_cell(phonon):
    """Return the phonopy primitive cell across supported API variants."""
    return _get_attr_or_call(phonon, "primitive", "get_primitive")


def get_supercell(phonon):
    """Return the phonopy supercell across supported API variants."""
    return _get_attr_or_call(phonon, "supercell", "get_supercell")


def get_force_constants(phonon):
    """Return force constants across supported API variants."""
    if hasattr(phonon, "force_constants"):
        return phonon.force_constants
    if hasattr(phonon, "get_force_constants"):
        return phonon.get_force_constants()
    raise AttributeError("Phonopy object does not expose force constants.")


def _get_dense_smallest_vectors(primitive):
    svecs, multi = primitive.get_smallest_vectors()
    if getattr(primitive, "store_dense_svecs", True):
        return np.array(svecs, dtype="double", order="C"), np.array(
            multi, dtype="int64", order="C"
        )

    from phonopy.structure.cells import sparse_to_dense_svecs

    dense_svecs, dense_multi = sparse_to_dense_svecs(svecs, multi)
    return np.array(dense_svecs, dtype="double", order="C"), np.array(
        dense_multi, dtype="int64", order="C"
    )


def get_primitive_structure_payload(phonon):
    """Serialize the primitive-cell structural data needed by the website."""
    primitive = get_primitive_cell(phonon)
    lattice = np.array(_get_attr_or_call(primitive, "cell", "get_cell"), dtype=float)
    atom_pos_red = np.array(
        _get_attr_or_call(primitive, "scaled_positions", "get_scaled_positions"),
        dtype=float,
    )
    atom_numbers = np.array(
        _get_attr_or_call(primitive, "numbers", "get_atomic_numbers"), dtype=int
    )
    masses = np.array(_get_attr_or_call(primitive, "masses", "get_masses"), dtype=float)
    symbols = list(_get_attr_or_call(primitive, "symbols", "get_chemical_symbols"))

    return {
        "natoms": int(len(primitive)),
        "lattice": lattice,
        "atom_pos_red": atom_pos_red,
        "atom_pos_car": red_car(atom_pos_red, lattice),
        "atom_numbers": atom_numbers,
        "atom_types": symbols,
        "masses": masses,
    }


def prepare_phonon_for_runtime_export(phonon, symmetrize_force_constants=True):
    """Return a phonopy object prepared for runtime export.

    Runtime payloads benefit from phonopy's stronger force-constant
    symmetrization. Work on a copy so the caller's phonon object keeps its
    original state for any other export path.
    """

    if not symmetrize_force_constants or not hasattr(phonon, "symmetrize_force_constants"):
        return phonon

    runtime_phonon = copy.deepcopy(phonon)
    runtime_phonon.symmetrize_force_constants(show_drift=False)
    return runtime_phonon


def get_runtime_dynamical_matrix_payload(phonon):
    """Serialize the minimal phonopy data needed to rebuild D(q).

    The returned payload matches the data consumed by
    ``phonopy.harmonic.dynamical_matrix.DynamicalMatrix._run_py_dynamical_matrix``
    after the primitive/supercell bookkeeping has already been resolved offline.
    """

    primitive = get_primitive_cell(phonon)
    supercell = get_supercell(phonon)
    force_constants = np.array(get_force_constants(phonon), dtype="double", order="C")
    primitive_lattice = np.array(
        _get_attr_or_call(primitive, "cell", "get_cell"), dtype="double", order="C"
    )
    positions_car = np.array(
        _get_attr_or_call(primitive, "positions", "get_positions"),
        dtype="double",
        order="C",
    )

    p2s_map = np.array(primitive.p2s_map, dtype="int64")
    s2p_map = np.array(primitive.s2p_map, dtype="int64")
    p2p_map = primitive.p2p_map
    s2pp_map = np.array(
        [p2p_map[s2p_map[i]] for i in range(len(s2p_map))], dtype="int64"
    )

    shortest_vectors, multiplicity = _get_dense_smallest_vectors(primitive)

    nac_params = None
    if hasattr(phonon, "nac_params"):
        nac_params = phonon.nac_params
    elif hasattr(phonon, "get_nac_params"):
        nac_params = phonon.get_nac_params()

    nac_payload = None
    dm = None
    if nac_params is not None and hasattr(phonon, "_set_dynamical_matrix"):
        phonon._set_dynamical_matrix()
        dm = getattr(phonon, "dynamical_matrix", None) or getattr(phonon, "_dynamical_matrix", None)

    if dm is not None and dm.__class__.__name__ == "DynamicalMatrixGL":
        dm.make_Gonze_nac_dataset()
        gonze_fc, dd_q0, _g_cutoff, g_list, lambda_value = dm.Gonze_nac_dataset
        force_constants_source = np.array(gonze_fc, dtype="double", order="C")
        nac_payload = {
            "method": "gonze",
            "born": np.array(nac_params["born"], dtype="double", order="C"),
            "dielectric": np.array(nac_params["dielectric"], dtype="double", order="C"),
            "nac_factor": float(dm.nac_factor),
            "positions_car": positions_car,
            "g_list": np.array(g_list, dtype="double", order="C"),
            "dd_q0": {
                "real": np.array(dd_q0.real, dtype="double", order="C"),
                "imag": np.array(dd_q0.imag, dtype="double", order="C"),
            },
            "lambda": float(lambda_value),
            "q_direction_tolerance": float(dm.Q_DIRECTION_TOLERANCE),
        }
    else:
        force_constants_source = force_constants
        if nac_params is not None:
            nac_payload = {
                "method": "wang",
                "born": np.array(nac_params["born"], dtype="double", order="C"),
                "dielectric": np.array(nac_params["dielectric"], dtype="double", order="C"),
                "factor": float(nac_params["factor"]),
            }

    if force_constants_source.shape[0] == force_constants_source.shape[1]:
        force_constants_compact = force_constants_source[p2s_map]
    else:
        force_constants_compact = force_constants_source

    payload = {
        "format": "phonopy-dynamical-matrix-v1",
        "acoustic_sum_rule": "off",
        "primitive_natoms": int(len(primitive)),
        "supercell_natoms": int(len(supercell)),
        "primitive_lattice": primitive_lattice,
        "frequency_conversion_factor": float(
            getattr(phonon, "unit_conversion_factor", None)
            or getattr(phonon, "factor", None)
            or 1.0
        ),
        "masses": np.array(primitive.masses, dtype="double", order="C"),
        "force_constants_compact": np.array(
            force_constants_compact, dtype="double", order="C"
        ),
        "shortest_vectors": shortest_vectors,
        "multiplicity": multiplicity,
        "s2pp_map": s2pp_map,
    }
    if nac_payload is not None:
        payload["nac"] = nac_payload
    return payload


def build_dynamical_matrix_from_payload(payload, qpoint, q_direction=None):
    """Reference NumPy implementation of the runtime dynamical-matrix builder."""
    q = np.array(qpoint, dtype="double")
    masses = np.array(payload["masses"], dtype="double")
    force_constants = np.array(payload["force_constants_compact"], dtype="double")
    shortest_vectors = np.array(payload["shortest_vectors"], dtype="double")
    multiplicity = np.array(payload["multiplicity"], dtype="int64")
    s2pp_map = np.array(payload["s2pp_map"], dtype="int64")

    natoms = len(masses)
    dm = np.zeros((natoms * 3, natoms * 3), dtype="complex128", order="C")

    if payload.get("nac", {}).get("method") == "gonze":
        rec_lat = np.linalg.inv(np.array(payload["primitive_lattice"], dtype="double"))
        q_cart = rec_lat @ q
        q_dir_cart = None
        if q_direction is not None:
            q_dir_cart = rec_lat @ np.array(q_direction, dtype="double")
        born = np.array(payload["nac"]["born"], dtype="double")
        dielectric = np.array(payload["nac"]["dielectric"], dtype="double")
        positions_car = np.array(payload["nac"]["positions_car"], dtype="double")
        g_list = np.array(payload["nac"]["g_list"], dtype="double")
        dd_q0 = np.array(payload["nac"]["dd_q0"]["real"], dtype="double") + 1j * np.array(
            payload["nac"]["dd_q0"]["imag"], dtype="double"
        )
        lambda_value = float(payload["nac"]["lambda"])
        tolerance = float(payload["nac"]["q_direction_tolerance"])

        dd_in = np.zeros((natoms, natoms, 3, 3), dtype="complex128")
        l2 = 4 * lambda_value * lambda_value
        for g in g_list:
            q_k = g + q_cart
            norm = np.linalg.norm(q_k)
            if norm < tolerance:
                if q_dir_cart is None:
                    kk = np.zeros((3, 3), dtype="double")
                else:
                    dielectric_part = q_dir_cart @ dielectric @ q_dir_cart
                    kk = np.outer(q_dir_cart, q_dir_cart) / dielectric_part
            else:
                dielectric_part = q_k @ dielectric @ q_k
                kk = np.outer(q_k, q_k) / dielectric_part * np.exp(-dielectric_part / l2)

            for i in range(natoms):
                for j in range(natoms):
                    phase = 2j * np.pi * np.dot(positions_car[i] - positions_car[j], g)
                    dd_in[i, j] += kk * np.exp(phase)

        dd = np.zeros((natoms, natoms, 3, 3), dtype="complex128")
        for i in range(natoms):
            for j in range(natoms):
                for alpha in range(3):
                    for beta in range(3):
                        value = 0j
                        for alpha_prime in range(3):
                            for beta_prime in range(3):
                                value += (
                                    dd_in[i, j, alpha_prime, beta_prime]
                                    * born[i, alpha_prime, alpha]
                                    * born[j, beta_prime, beta]
                                )
                        dd[i, j, alpha, beta] = value

        for i in range(natoms):
            dd[i, i] -= dd_q0[i]

        dd *= float(payload["nac"]["nac_factor"])
        for i in range(natoms):
            for j in range(natoms):
                dm[(i * 3) : (i * 3 + 3), (j * 3) : (j * 3 + 3)] += dd[i, j] / np.sqrt(
                    masses[i] * masses[j]
                )

    for i in range(natoms):
        for j in range(natoms):
            sqrt_mm = np.sqrt(masses[i] * masses[j])
            block = np.zeros((3, 3), dtype="complex128", order="C")
            for super_index, primitive_index in enumerate(s2pp_map):
                if primitive_index != j:
                    continue
                multiplicity_count, address = multiplicity[super_index, i]
                svecs = shortest_vectors[address : address + multiplicity_count]
                phase_factor = np.exp(2j * np.pi * np.dot(svecs, q)).sum()
                block += (
                    force_constants[i, super_index]
                    * phase_factor
                    / sqrt_mm
                    / multiplicity_count
                )
            dm[(i * 3) : (i * 3 + 3), (j * 3) : (j * 3 + 3)] += block

    return (dm + dm.conj().T) / 2
