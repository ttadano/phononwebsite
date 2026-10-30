from pathlib import Path

import numpy as np
import pytest

phonopy = pytest.importorskip("phonopy")

from phononweb.phonopyphonon import PhonopyPhonon
from phononweb.runtime_dynamical_matrix import (
    build_dynamical_matrix_from_payload,
    get_runtime_dynamical_matrix_payload,
)


def test_runtime_payload_rebuilds_phonopy_dynamical_matrix():
    fixture_root = (
        Path(__file__).resolve().parents[3]
        / "test"
        / "fixtures"
        / "phonondb"
        / "2015"
        / "mp-149"
        / "gruneisen-00"
    )
    phonon_yaml = fixture_root / "phonon.yaml"
    force_sets = fixture_root / "FORCE_SETS"

    phonopy_phonon = PhonopyPhonon.from_files(str(phonon_yaml), str(force_sets))
    phonon = phonopy_phonon.phonon
    payload = get_runtime_dynamical_matrix_payload(phonon)

    assert payload["format"] == "phonopy-dynamical-matrix-v1"
    assert payload["force_constants_compact"].shape[0] == payload["primitive_natoms"]
    assert payload["force_constants_compact"].shape[1] == payload["supercell_natoms"]
    assert payload["multiplicity"].shape[:2] == (
        payload["supercell_natoms"],
        payload["primitive_natoms"],
    )

    qpoint = np.array([0.25, 0.0, 0.25], dtype=float)
    expected = phonon.get_dynamical_matrix_at_q(qpoint)
    rebuilt = build_dynamical_matrix_from_payload(payload, qpoint)

    np.testing.assert_allclose(rebuilt, expected, atol=1e-10)


def test_runtime_payload_rebuilds_phonopy_gonze_dynamical_matrix():
    archive_path = Path('/Users/henriquemiranda/phonondb/phonondb2017/mp-1000.tar.lzma')
    if not archive_path.exists():
        pytest.skip('local PhononDB archive with NAC not available')

    import tarfile
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(archive_path, mode='r:xz') as tf:
            tf.extractall(tmp)
        fixture_root = Path(tmp) / 'mp-1000'
        phonon_yaml = fixture_root / 'phonon.yaml'
        force_sets = fixture_root / 'FORCE_SETS'
        born = fixture_root / 'BORN'

        phonopy_phonon = PhonopyPhonon.from_files(str(phonon_yaml), str(force_sets), nac_filename=str(born))
        phonon = phonopy_phonon.phonon
        payload = get_runtime_dynamical_matrix_payload(phonon)

        assert payload['nac']['method'] == 'gonze'
        qpoint = np.array([0.25, 0.0, 0.25], dtype=float)
        expected = phonon.get_dynamical_matrix_at_q(qpoint)
        rebuilt = build_dynamical_matrix_from_payload(payload, qpoint)

        np.testing.assert_allclose(rebuilt, expected, atol=1e-8)
