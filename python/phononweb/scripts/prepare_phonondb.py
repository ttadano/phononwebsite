#!/usr/bin/env python
# Copyright (c) 2026, Henrique Miranda
# All rights reserved.
#
# This file is part of the phononwebsite project
#
"""Prepare ready-to-serve compressed JSON files from raw PhononDB archives."""

import argparse
import base64
import concurrent.futures
import gzip
import json
import sys
import tarfile
import tempfile
import time
import zlib
from pathlib import Path

import numpy as np
import yaml

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phononweb.jsonencoder import JsonEncoder
from phononweb.lattice import red_car
from phononweb.phonopyphonon import PhonopyPhonon
from phononweb.runtime_dynamical_matrix import (
    get_runtime_dynamical_matrix_payload,
    prepare_phonon_for_runtime_export,
)
from phononweb.units import atomic_numbers
from phononweb.units import atomic_mass
from phononweb.utils import estimate_band_connection

THZ_TO_CM1 = 33.35641


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert raw PhononDB tar.lzma archives into the compressed JSON format "
            "used by the phonon website."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more .tar.lzma archives or directories containing such archives.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory where the generated .json.gz files will be written.",
    )
    parser.add_argument(
        "--band-points",
        type=int,
        default=15,
        help="Number of points per high-symmetry segment.",
    )
    parser.add_argument(
        "--repetitions",
        default="3,3,3",
        help="Default repetitions stored in the generated website JSON.",
    )
    parser.add_argument(
        "--name-from",
        choices=["mp-id", "formula"],
        default="mp-id",
        help="How to name the generated output files and material entries.",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Optional JSON file receiving a materials manifest for the generated files.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of archives to convert after sorting the inputs.",
    )
    parser.add_argument(
        "--max-atoms",
        type=int,
        default=None,
        help="Optional maximum number of atoms in the primitive cell; larger archives are skipped.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip archives whose output .json.gz already exists.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of archives to convert in parallel. Use 1 for serial execution.",
    )
    parser.add_argument(
        "--vector-decimals",
        type=int,
        default=4,
        help="Decimal digits kept for eigenvector components in the generated JSON. Use a negative value to disable rounding.",
    )
    parser.add_argument(
        "--vector-format",
        choices=["json", "q11-int16-base64", "runtime"],
        default="q11-int16-base64",
        help="How eigenvectors are stored in the generated JSON payload.",
    )
    return parser.parse_args()


def expand_input_paths(raw_inputs):
    paths = []
    for item in raw_inputs:
        path = Path(item).expanduser().resolve()
        if path.is_dir():
            paths.extend(sorted(path.glob("*.tar.lzma")))
        elif path.is_file():
            paths.append(path)
    return paths


def parse_repetitions(text):
    chunks = [chunk.strip() for chunk in text.replace(",", " ").split() if chunk.strip()]
    if len(chunks) != 3:
        raise ValueError("Repetitions must contain exactly three integers, e.g. 3,3,3")
    return [int(value) for value in chunks]


def extract_archive(archive_path: Path, workdir: Path) -> Path:
    with tarfile.open(archive_path, mode="r:xz") as tf:
        tf.extractall(workdir)

    entries = [entry for entry in workdir.iterdir() if entry.is_dir()]
    if len(entries) == 1:
        return entries[0]
    return workdir


def load_band_yaml(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_gzip_json(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def choose_output_stem(archive_path: Path, data, mode: str):
    if mode == "formula":
        formula = data.get("formula") or data.get("name")
        if formula:
            return formula.replace(" ", "")
    return archive_path.stem.replace(".tar", "")


def get_archive_primitive_natoms(archive_path: Path):
    with tempfile.TemporaryDirectory(prefix="phonondb-scan-") as tmp:
        extracted_root = extract_archive(archive_path, Path(tmp))
        phonon_yaml = extracted_root / "phonon.yaml"
        obj = load_band_yaml(phonon_yaml)
        return len(obj["points"])


def build_highsym_and_breaks(phonon_entries):
    highsym = []
    line_breaks = []
    current_start = 0

    if phonon_entries and phonon_entries[0].get("label"):
        highsym.append([0, phonon_entries[0]["label"]])

    for index in range(1, len(phonon_entries)):
        current = phonon_entries[index]
        previous = phonon_entries[index - 1]
        current_label = current.get("label")
        previous_label = previous.get("label")

        if current_label and previous_label and current_label != previous_label:
            highsym.append([index - 1, previous_label])
            line_breaks.append([current_start, index])
            current_start = index
            highsym.append([index, current_label])

    if phonon_entries:
        final_label = phonon_entries[-1].get("label")
        if final_label:
            if not highsym or highsym[-1][0] != len(phonon_entries) - 1:
                highsym.append([len(phonon_entries) - 1, final_label])

    if phonon_entries:
        line_breaks.append([current_start, len(phonon_entries)])

    return highsym, line_breaks


def build_path_metadata(phonopy_phonon):
    bands = getattr(phonopy_phonon, "bands", None)
    labels = getattr(phonopy_phonon, "labels", None)
    if not bands:
        return [], []

    line_breaks = []
    highsym_qpts = []
    start_index = 0

    if labels and labels[0]:
        highsym_qpts.append([0, labels[0]])

    for segment_index, branch in enumerate(bands):
        branch_length = len(branch)
        end_index = start_index + branch_length
        line_breaks.append([start_index, end_index])

        if labels and segment_index + 1 < len(labels) and labels[segment_index + 1]:
            highsym_qpts.append([end_index - 1, labels[segment_index + 1]])

        start_index = end_index

    return highsym_qpts, line_breaks


def build_formula(atom_numbers_list):
    counts = {}
    for atomic_number in atom_numbers_list:
        counts[atomic_number] = counts.get(atomic_number, 0) + 1

    items = sorted(
        ((atomic_number, counts[atomic_number]) for atomic_number in counts),
        key=lambda item: atomic_numbers_inverse_sort_key(item[0]),
    )
    formula = []
    for atomic_number, count in items:
        symbol = atomic_number_to_symbol(atomic_number)
        formula.append(symbol)
        if count > 1:
            formula.append(str(count))
    return "".join(formula)


def atomic_number_to_symbol(atomic_number):
    for symbol, value in atomic_numbers.items():
        if value == atomic_number:
            return symbol
    raise KeyError(f"Unknown atomic number {atomic_number}")


def atomic_numbers_inverse_sort_key(atomic_number):
    symbol = atomic_number_to_symbol(atomic_number)
    if symbol == "H":
        return (0, symbol)
    if symbol == "C":
        return (1, symbol)
    return (2, symbol)


def get_structure_metadata(phonopy_phonon):
    primitive = (
        phonopy_phonon.phonon.get_primitive()
        if hasattr(phonopy_phonon.phonon, "get_primitive")
        else phonopy_phonon.phonon.primitive
    )
    lattice = primitive.get_cell() if hasattr(primitive, "get_cell") else primitive.cell
    atom_pos_red = (
        primitive.get_scaled_positions()
        if hasattr(primitive, "get_scaled_positions")
        else primitive.scaled_positions
    )
    atom_numbers_list = (
        primitive.get_atomic_numbers()
        if hasattr(primitive, "get_atomic_numbers")
        else primitive.numbers
    )
    atom_types = [atomic_number_to_symbol(int(number)) for number in atom_numbers_list]
    atom_pos_car = red_car(atom_pos_red, lattice)
    formula = build_formula(atom_numbers_list)
    return {
        "lattice": lattice,
        "atom_pos_red": atom_pos_red,
        "atom_pos_car": atom_pos_car,
        "atom_numbers": [int(number) for number in atom_numbers_list],
        "atom_types": atom_types,
        "masses": [float(mass) for mass in primitive.masses],
        "formula": formula,
    }


def apply_average_mass_normalization(payload):
    atom_numbers_list = payload.get("atom_numbers") or []
    if not atom_numbers_list:
        return payload

    masses = np.array([float(atomic_mass[int(number)]) for number in atom_numbers_list], dtype=float)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        return payload

    average_mass = float(np.mean(masses))
    scale = np.sqrt(average_mass / masses).reshape(1, 1, len(masses), 1, 1)
    vectors = np.array(payload["vectors"], dtype=float)
    payload["vectors"] = vectors * scale
    payload["average_mass"] = average_mass
    payload["mode_amplitude_convention"] = "avg-mass-normalized"
    return payload


def add_average_mass_metadata(payload):
    atom_numbers_list = payload.get("atom_numbers") or []
    if not atom_numbers_list:
        return payload

    masses = np.array([float(atomic_mass[int(number)]) for number in atom_numbers_list], dtype=float)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        return payload

    payload["average_mass"] = float(np.mean(masses))
    payload["mode_amplitude_convention"] = "avg-mass-normalized"
    return payload


def convert_seekpath_bands_to_runtime_json(
    phonopy_phonon,
    structure,
    repetitions,
    archive_path,
    name_mode,
):
    bands = getattr(phonopy_phonon, "bands", None) or []
    highsym_qpts, line_breaks = build_path_metadata(phonopy_phonon)
    reciprocal_metric = np.linalg.inv(np.array(structure["lattice"], dtype=float)).T

    qpoints = []
    distances = []
    cumulative_distance = 0.0

    for branch in bands:
        branch = np.array(branch, dtype=float)
        if len(branch) == 0:
            continue
        previous_qpoint = branch[0].copy()
        for branch_index, qpoint in enumerate(branch):
            if branch_index > 0:
                cumulative_distance += float(
                    np.linalg.norm(np.dot(qpoint - previous_qpoint, reciprocal_metric))
                )
                previous_qpoint = qpoint.copy()
            qpoints.append([float(value) for value in qpoint])
            distances.append(cumulative_distance)

    payload = {
        "name": choose_output_stem(archive_path, structure, name_mode),
        "natoms": len(structure["atom_types"]),
        "lattice": structure["lattice"],
        "atom_types": structure["atom_types"],
        "atom_numbers": structure["atom_numbers"],
        "formula": structure["formula"],
        "qpoints": qpoints,
        "repetitions": repetitions,
        "atom_pos_red": structure["atom_pos_red"],
        "atom_pos_car": structure["atom_pos_car"],
        "masses": structure["masses"],
        "distances": distances,
        "highsym_qpts": highsym_qpts,
        "line_breaks": line_breaks,
    }
    add_average_mass_metadata(payload)
    return payload


def convert_band_yaml_to_site_json(
    band_yaml,
    structure,
    repetitions,
    archive_path,
    name_mode,
    path_metadata=None,
):
    phonon_entries = band_yaml["phonon"]
    natoms = len(structure["atom_types"])
    nmodes = natoms * 3

    data = {
        "name": choose_output_stem(archive_path, structure, name_mode),
        "natoms": natoms,
        "lattice": structure["lattice"],
        "atom_types": structure["atom_types"],
        "atom_numbers": structure["atom_numbers"],
        "formula": structure["formula"],
        "qpoints": [],
        "repetitions": repetitions,
        "atom_pos_red": structure["atom_pos_red"],
        "atom_pos_car": structure["atom_pos_car"],
        "masses": structure["masses"],
        "eigenvalues": [],
        "distances": [],
        "highsym_qpts": [],
        "line_breaks": [],
        "vectors": [],
    }

    for qpoint in phonon_entries:
        data["qpoints"].append(qpoint["q-position"])
        data["distances"].append(qpoint["distance"])

        eig = []
        vec = []
        bands = qpoint["band"]
        for mode_index in range(nmodes):
            band = bands[mode_index]
            eig.append(band["frequency"] * THZ_TO_CM1)
            mode_vectors = []
            for atom_index in range(natoms):
                components = band["eigenvector"][atom_index]
                mode_vectors.append([
                    [float(components[0][0]), float(components[0][1])],
                    [float(components[1][0]), float(components[1][1])],
                    [float(components[2][0]), float(components[2][1])],
                ])
            vec.append(mode_vectors)
        data["eigenvalues"].append(eig)
        data["vectors"].append(vec)

    if path_metadata is None:
        highsym, line_breaks = build_highsym_and_breaks(phonon_entries)
    else:
        highsym, line_breaks = path_metadata
    data["highsym_qpts"] = highsym
    data["line_breaks"] = line_breaks

    return data


def reorder_payload_band_connection(payload):
    eigenvalues = np.array(payload["eigenvalues"], dtype=float)
    vectors = np.array(payload["vectors"], dtype=float)

    if eigenvalues.ndim != 2 or vectors.ndim != 5:
        return payload

    nqpoints, nphons = eigenvalues.shape
    complex_vectors = vectors[..., 0] + 1j * vectors[..., 1]
    complex_vectors = complex_vectors.reshape(nqpoints, nphons, nphons)

    ordered_eigenvalues = np.empty_like(eigenvalues)
    ordered_vectors = np.empty_like(vectors)
    band_order = list(range(nphons))
    ordered_eigenvalues[0] = eigenvalues[0]
    ordered_vectors[0] = vectors[0]

    for qpoint_index in range(1, nqpoints):
        band_order = estimate_band_connection(
            complex_vectors[qpoint_index - 1].T,
            complex_vectors[qpoint_index].T,
            band_order,
        )
        ordered_eigenvalues[qpoint_index] = eigenvalues[qpoint_index, band_order]
        ordered_vectors[qpoint_index] = vectors[qpoint_index, band_order]

    payload["eigenvalues"] = ordered_eigenvalues
    payload["vectors"] = ordered_vectors
    return payload


def quantize_payload(payload, vector_decimals):
    if vector_decimals is None or vector_decimals < 0:
        return payload

    payload["vectors"] = np.round(np.array(payload["vectors"], dtype=float), decimals=vector_decimals)
    return payload


def encode_vectors_q11_int16_base64(payload):
    vectors = np.array(payload["vectors"], dtype=np.float32)
    quantized = np.clip(np.rint(vectors * 2047.0), -2047, 2047).astype(np.int16)
    payload["vectors_compressed"] = {
        "format": "q11-int16-base64",
        "scale": 2047,
        "shape": list(quantized.shape),
        "data": base64.b64encode(quantized.tobytes(order="C")).decode("ascii"),
    }
    del payload["vectors"]
    return payload


def write_gzip_json(path: Path, payload):
    encoded = json.dumps(payload, cls=JsonEncoder, separators=(",", ":")).encode("utf-8")
    with gzip.open(path, "wb", compresslevel=9) as handle:
        handle.write(encoded)


def build_manifest_entry(payload, output_path: Path):
    stem = output_path.stem.replace(".json", "")
    material_id = stem[3:] if stem.startswith("mp-") else stem
    return {
        "id": material_id,
        "name": payload["formula"],
        "file": output_path.name,
        "type": "json",
        "link": f"https://materialsproject.org/materials/mp-{material_id}" if stem.startswith("mp-") else None,
    }


def prepare_archive(
    archive_path: Path,
    output_dir: Path,
    repetitions,
    band_points: int,
    name_mode: str,
    vector_decimals: int = 4,
    vector_format: str = "q11-int16-base64",
):
    with tempfile.TemporaryDirectory(prefix="phonondb-prepare-") as tmp:
        tmpdir = Path(tmp)
        extracted_root = extract_archive(archive_path, tmpdir)

        phonon_yaml = extracted_root / "phonon.yaml"
        force_sets = extracted_root / "FORCE_SETS"
        born = extracted_root / "BORN"
        if not phonon_yaml.exists() or not force_sets.exists():
            raise FileNotFoundError(f"{archive_path} is missing phonon.yaml or FORCE_SETS")

        phonopy_phonon = PhonopyPhonon.from_files(
            str(phonon_yaml),
            str(force_sets),
            nac_filename=str(born) if born.exists() else None,
        )
        phonopy_phonon.set_bandstructure_seekpath_points(band_points=band_points)
        structure = get_structure_metadata(phonopy_phonon)

        if vector_format == "runtime":
            runtime_phonon = prepare_phonon_for_runtime_export(phonopy_phonon.phonon)
            payload = convert_seekpath_bands_to_runtime_json(
                phonopy_phonon,
                structure=structure,
                repetitions=repetitions,
                archive_path=archive_path,
                name_mode=name_mode,
            )
            payload["dynamical_matrix"] = get_runtime_dynamical_matrix_payload(runtime_phonon)
        else:
            phonopy_phonon.get_bandstructure(is_eigenvectors=True, is_band_connection=True)
            band_yaml_path = tmpdir / "band.yaml"
            phonopy_phonon.write_band_yaml(filename=str(band_yaml_path))
            band_yaml = load_band_yaml(band_yaml_path)
            path_metadata = build_path_metadata(phonopy_phonon)

            payload = convert_band_yaml_to_site_json(
                band_yaml,
                structure=structure,
                repetitions=repetitions,
                archive_path=archive_path,
                name_mode=name_mode,
                path_metadata=path_metadata,
            )
            reorder_payload_band_connection(payload)
            apply_average_mass_normalization(payload)
            if vector_format == "json":
                quantize_payload(payload, vector_decimals)
            elif vector_format == "q11-int16-base64":
                encode_vectors_q11_int16_base64(payload)

        output_stem = choose_output_stem(archive_path, payload, name_mode)
        payload["name"] = payload["formula"]
        output_path = output_dir / f"{output_stem}.json.gz"
        write_gzip_json(output_path, payload)
        return payload, output_path


def prepare_archive_task(task):
    archive_path = Path(task["archive_path"])
    output_dir = Path(task["output_dir"])
    repetitions = task["repetitions"]
    band_points = task["band_points"]
    name_mode = task["name_mode"]
    vector_decimals = task["vector_decimals"]
    vector_format = task["vector_format"]
    started_at = time.perf_counter()
    payload, output_path = prepare_archive(
        archive_path,
        output_dir=output_dir,
        repetitions=repetitions,
        band_points=band_points,
        name_mode=name_mode,
        vector_decimals=vector_decimals,
        vector_format=vector_format,
    )
    elapsed = time.perf_counter() - started_at
    return {
        "index": task["index"],
        "total": task["total"],
        "archive_name": archive_path.name,
        "payload": payload,
        "output_path": str(output_path),
        "elapsed": elapsed,
    }


def main():
    args = parse_args()
    repetitions = parse_repetitions(args.repetitions)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    archives = expand_input_paths(args.inputs)
    if not archives:
        raise SystemExit("No .tar.lzma archives were found in the provided inputs.")
    if args.limit is not None:
        archives = archives[: max(0, args.limit)]

    manifest = []
    total_archives = len(archives)
    pending_tasks = []
    for index, archive in enumerate(archives, start=1):
        if args.max_atoms is not None:
            natoms = get_archive_primitive_natoms(archive)
            if natoms > args.max_atoms:
                print(
                    f"[{index}/{total_archives}] Skipping {archive.name}: {natoms} atoms > max {args.max_atoms}",
                    flush=True,
                )
                continue

        if args.skip_existing and args.name_from == "mp-id":
            output_stem = archive.stem.replace(".tar", "")
            output_path = output_dir / f"{output_stem}.json.gz"
            if output_path.exists():
                try:
                    payload = load_gzip_json(output_path)
                except (OSError, EOFError, json.JSONDecodeError, zlib.error) as exc:
                    print(
                        f"[{index}/{total_archives}] Existing output is unreadable, rebuilding {archive.name}: "
                        f"{output_path} ({exc})",
                        flush=True,
                    )
                else:
                    manifest.append(build_manifest_entry(payload, output_path))
                    print(
                        f"[{index}/{total_archives}] Skipping {archive.name}: {output_path.name} already exists",
                        flush=True,
                    )
                    continue

        print(f"[{index}/{total_archives}] Queued {archive.name}...", flush=True)
        pending_tasks.append({
            "index": index,
            "total": total_archives,
            "archive_path": str(archive),
            "output_dir": str(output_dir),
            "repetitions": repetitions,
            "band_points": args.band_points,
            "name_mode": args.name_from,
            "vector_decimals": args.vector_decimals,
            "vector_format": args.vector_format,
        })

    jobs = max(1, int(args.jobs))
    if jobs == 1:
        for task in pending_tasks:
            print(f"[{task['index']}/{task['total']}] Processing {Path(task['archive_path']).name}...", flush=True)
            result = prepare_archive_task(task)
            output_path = Path(result["output_path"])
            manifest.append(build_manifest_entry(result["payload"], output_path))
            print(
                f"[{result['index']}/{result['total']}] Wrote {output_path} in {result['elapsed']:.1f}s",
                flush=True,
            )
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
            future_map = {}
            for task in pending_tasks:
                print(
                    f"[{task['index']}/{task['total']}] Processing {Path(task['archive_path']).name} "
                    f"(worker pool, jobs={jobs})...",
                    flush=True,
                )
                future = executor.submit(prepare_archive_task, task)
                future_map[future] = task

            results = []
            for future in concurrent.futures.as_completed(future_map):
                result = future.result()
                output_path = Path(result["output_path"])
                results.append((result["index"], build_manifest_entry(result["payload"], output_path)))
                print(
                    f"[{result['index']}/{result['total']}] Wrote {output_path} in {result['elapsed']:.1f}s",
                    flush=True,
                )

            for _, manifest_entry in sorted(results, key=lambda item: item[0]):
                manifest.append(manifest_entry)

    if args.manifest:
        manifest_path = Path(args.manifest).expanduser().resolve()
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with manifest_path.open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)
        print(f"Wrote manifest {manifest_path}")


if __name__ == "__main__":
    main()
