from pathlib import Path

import pytest

phonopy = pytest.importorskip('phonopy')

from phononweb.phonopyphonon import PhonopyPhonon


def test_phonopy_from_files_and_band_yaml(tmp_path):
    fixture_root = Path(__file__).resolve().parents[3] / 'test' / 'fixtures' / 'phonondb' / '2015' / 'mp-149' / 'gruneisen-00'
    phonon_yaml = fixture_root / 'phonon.yaml'
    force_sets = fixture_root / 'FORCE_SETS'

    phonon = PhonopyPhonon.from_files(str(phonon_yaml), str(force_sets))
    phonon.set_bandstructure_seekpath()
    phonon.get_bandstructure()

    outfile = tmp_path / 'band.yaml'
    phonon.write_band_yaml(filename=str(outfile))

    assert outfile.exists()
    assert outfile.stat().st_size > 0
