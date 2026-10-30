from pathlib import Path
import json
import shutil
import sys

import numpy as np
import pytest

pytest.importorskip('h5py')

from phononweb.alamodephonon import AlamodePhonon
from phononweb.scripts.read_alamode_phonon import main as read_alamode_phonon_main

# PbTe band structure (G-X, X-G, G-L; 6 points per segment) computed with ALAMODE 2.0.0
# anphon (PRINTEVEC = 1); pbte.bands is the band structure written by the same run.
FIXTURE_ROOT = Path(__file__).resolve().parents[3] / 'test' / 'fixtures' / 'alamode'
FIXTURE_EVEC = FIXTURE_ROOT / 'pbte.band.evec.h5'
FIXTURE_BANDS = FIXTURE_ROOT / 'pbte.bands'


def test_alamode_reads_evec_h5_with_absolute_path():
    phonon = AlamodePhonon(str(FIXTURE_EVEC), 'PbTe')

    assert phonon.natoms == 2
    assert phonon.nqpoints == 18
    assert phonon.nphons == 6
    assert phonon.atom_types == ['Pb', 'Te']
    assert phonon.atom_numbers == [82, 52]
    assert phonon.chemical_formula == 'PbTe'
    np.testing.assert_allclose(phonon.pos, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], atol=1e-8)
    # fcc lattice in angstrom
    np.testing.assert_allclose(np.diag(phonon.cell), 0.0, atol=1e-8)
    assert 3.0 < phonon.cell[0][1] < 3.5
    assert phonon.eigenvectors.shape == (18, 6, 2, 3, 2)


def test_alamode_frequencies_match_anphon_bands_file():
    phonon = AlamodePhonon(str(FIXTURE_EVEC), 'PbTe')
    bands = np.loadtxt(FIXTURE_BANDS)

    np.testing.assert_allclose(phonon.eigenvalues, bands[:, 1:], atol=1e-4)


def test_alamode_acoustic_modes_at_gamma_are_rigid_translations():
    # after the 1/sqrt(m) scaling all atoms move together in the acoustic modes at Gamma
    phonon = AlamodePhonon(str(FIXTURE_EVEC), 'PbTe')

    for mode in range(3):
        displacements = phonon.eigenvectors[0, mode]
        np.testing.assert_allclose(displacements[0], displacements[1], atol=1e-6)


def test_alamode_band_range():
    phonon = AlamodePhonon(str(FIXTURE_EVEC), 'PbTe', band_range='4 6')
    full = AlamodePhonon(str(FIXTURE_EVEC), 'PbTe')

    assert phonon.nphons == 3
    assert phonon.eigenvectors.shape == (18, 3, 2, 3, 2)
    np.testing.assert_allclose(phonon.eigenvalues, full.eigenvalues[:, 3:6])


def test_alamode_reads_legacy_hdf5_extension(tmp_path):
    # ALAMODE < 2.0 writes the same data to PREFIX.band.evec.hdf5
    legacy = tmp_path / 'pbte.band.evec.hdf5'
    shutil.copy(FIXTURE_EVEC, legacy)

    phonon = AlamodePhonon(legacy.name, 'PbTe', folder=str(tmp_path))

    assert phonon.nqpoints == 18


def test_alamode_missing_file_raises(tmp_path):
    with pytest.raises(ValueError):
        AlamodePhonon(str(tmp_path / 'missing.band.evec.h5'), 'PbTe')


def test_read_alamode_phonon_main_writes_json(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', [
        'read_alamode_phonon',
        str(FIXTURE_EVEC),
        'PbTe',
        '-w',
    ])

    read_alamode_phonon_main()

    data = json.loads((tmp_path / 'PbTe.json').read_text())
    assert data['name'] == 'PbTe'
    assert data['natoms'] == 2
    assert data['atom_types'] == ['Pb', 'Te']
    assert len(data['qpoints']) == 18
    assert len(data['eigenvalues']) == 18
    np.testing.assert_allclose(data['eigenvalues'], np.loadtxt(FIXTURE_BANDS)[:, 1:], atol=1e-4)
