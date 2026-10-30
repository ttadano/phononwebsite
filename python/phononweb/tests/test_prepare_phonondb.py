from pathlib import Path
import gzip
import json

import pytest

phonopy = pytest.importorskip('phonopy')
seekpath = pytest.importorskip('seekpath')

from phononweb.scripts.prepare_phonondb import prepare_archive


def test_prepare_archive_writes_compressed_internal_json(tmp_path):
    fixture_root = (
        Path(__file__).resolve().parents[3]
        / 'test'
        / 'fixtures'
        / 'phonondb'
        / '2015'
        / 'mp-149'
        / 'gruneisen-00'
    )

    archive_path = tmp_path / 'mp-149.tar.lzma'
    import tarfile
    with tarfile.open(archive_path, mode='w:xz') as tf:
        tf.add(fixture_root, arcname='mp-149')

    payload, output_path = prepare_archive(
        archive_path,
        output_dir=tmp_path,
        repetitions=[3, 3, 3],
        band_points=21,
        name_mode='mp-id',
    )

    assert output_path.exists()
    assert output_path.suffixes[-2:] == ['.json', '.gz']
    assert payload['formula'] == 'Si2'
    assert payload['repetitions'] == [3, 3, 3]
    assert len(payload['line_breaks']) > 0
    assert len(payload['qpoints']) == len(payload['distances'])
    assert len(payload['qpoints']) == len(payload['eigenvalues'])
    assert 'vectors_compressed' in payload
    assert 'vectors' not in payload

    with gzip.open(output_path, 'rt', encoding='utf-8') as handle:
        written = json.load(handle)

    assert written['formula'] == 'Si2'
    assert written['name'] == 'Si2'
    assert written['vectors_compressed']['format'] == 'q11-int16-base64'
    assert written['atom_pos_car'][1] != written['atom_pos_red'][1]


def test_prepare_archive_runtime_mode_writes_dynamical_matrix_payload(tmp_path):
    fixture_root = (
        Path(__file__).resolve().parents[3]
        / 'test'
        / 'fixtures'
        / 'phonondb'
        / '2015'
        / 'mp-149'
        / 'gruneisen-00'
    )

    archive_path = tmp_path / 'mp-149.tar.lzma'
    import tarfile
    with tarfile.open(archive_path, mode='w:xz') as tf:
        tf.add(fixture_root, arcname='mp-149')

    payload, output_path = prepare_archive(
        archive_path,
        output_dir=tmp_path,
        repetitions=[3, 3, 3],
        band_points=21,
        name_mode='mp-id',
        vector_format='runtime',
    )

    assert output_path.exists()
    assert 'vectors' not in payload
    assert 'vectors_compressed' not in payload
    assert 'eigenvalues' not in payload
    assert len(payload['qpoints']) == len(payload['distances'])
    assert len(payload['qpoints']) > 0
    assert len(payload['line_breaks']) > 0
    assert payload['mode_amplitude_convention'] == 'avg-mass-normalized'
    assert payload['average_mass'] > 0
    assert payload['dynamical_matrix']['format'] == 'phonopy-dynamical-matrix-v1'
    assert payload['dynamical_matrix']['acoustic_sum_rule'] == 'off'
    assert len(payload['dynamical_matrix']['masses']) == payload['natoms']
