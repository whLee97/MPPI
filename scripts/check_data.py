#!/usr/bin/env python3
"""Verify the ten local raw survey inputs against the reference hashes."""
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
args = parser.parse_args()
root = args.root.resolve()
expected = json.loads((root / "data/input_checksums.json").read_text())
failed = False
for name, digest in expected.items():
    path = root / "data/raw_dxa" / name
    if not path.is_file():
        print("MISSING", name)
        failed = True
        continue
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    matches = h.hexdigest() == digest
    print("OK" if matches else "MISMATCH", name)
    failed = failed or not matches
raise SystemExit(1 if failed else 0)
