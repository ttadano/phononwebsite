#!/usr/bin/env python
"""Export the minimal phonopy runtime payload needed for on-demand eigenvectors."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phononweb.jsonencoder import JsonEncoder
from phononweb.phonopyphonon import PhonopyPhonon
from phononweb.runtime_dynamical_matrix import (
    get_primitive_structure_payload,
    get_runtime_dynamical_matrix_payload,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Export the minimal phonopy data needed to rebuild dynamical matrices "
            "and eigenvectors at arbitrary q-points."
        )
    )
    parser.add_argument("phonon_yaml", help="Path to phonopy's phonon.yaml file.")
    parser.add_argument("force_sets", help="Path to the FORCE_SETS file.")
    parser.add_argument(
        "--born",
        default=None,
        help="Optional BORN file for NAC-enabled materials.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output JSON file. Use '-' to print to stdout.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    phonopy_phonon = PhonopyPhonon.from_files(
        args.phonon_yaml,
        args.force_sets,
        nac_filename=args.born,
    )
    phonon = phonopy_phonon.phonon
    payload = {
        "structure": get_primitive_structure_payload(phonon),
        "dynamical_matrix": get_runtime_dynamical_matrix_payload(phonon),
    }

    if args.pretty:
        text = json.dumps(payload, cls=JsonEncoder, indent=2)
    else:
        text = json.dumps(payload, cls=JsonEncoder, separators=(",", ":"))

    if args.output == "-":
        print(text)
        return

    Path(args.output).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
