#!/usr/bin/env python3
# Copyright (c) 2026, Henrique Miranda
# All rights reserved.
#
# This file is part of the phononwebsite project
#
"""Render index.html from README.md using pandoc and the homepage template."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


PLACEHOLDER = "PANDOC"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate index.html from README.md and ref_index.html."
    )
    parser.add_argument(
        "--readme",
        default="README.md",
        help="Markdown source used for the homepage body.",
    )
    parser.add_argument(
        "--template",
        default="ref_index.html",
        help="HTML template containing the PANDOC placeholder.",
    )
    parser.add_argument(
        "--output",
        default="index.html",
        help="Output homepage path.",
    )
    return parser.parse_args()


def render_markdown(readme_path: Path) -> str:
    command = [
        "pandoc",
        "--columns",
        "10000",
        "--email-obfuscation=javascript",
        str(readme_path),
        "-f",
        "markdown",
        "-t",
        "html",
    ]
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def main():
    args = parse_args()
    readme_path = Path(args.readme).resolve()
    template_path = Path(args.template).resolve()
    output_path = Path(args.output).resolve()

    rendered = render_markdown(readme_path)
    template = template_path.read_text(encoding="utf-8")
    if PLACEHOLDER not in template:
        raise RuntimeError(f"Template {template_path} does not contain {PLACEHOLDER!r}.")

    output = template.replace(PLACEHOLDER, rendered)
    output_path.write_text(output, encoding="utf-8")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
