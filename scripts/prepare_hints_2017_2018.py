#!/usr/bin/env python3
"""Prepare HINTS 5 Cycle 1/2 for the 2017-2018 DXA threshold application."""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw_dxa"
OUT_DIR = ROOT / "data" / "processed_dxa"


HINTS_FILES = [
    {
        "year": 2017,
        "cycle": "HINTS 5 Cycle 1",
        "zip": "hints5_cycle1_2017_stata.zip",
        "dta": "hints5_cycle1_public.dta",
    },
    {
        "year": 2018,
        "cycle": "HINTS 5 Cycle 2",
        "zip": "hints5_cycle2_2018_stata.zip",
        "dta": "hints5_cycle2_public.dta",
    },
]


def numeric_clean(x: pd.Series, valid_min: float, valid_max: float) -> pd.Series:
    out = pd.to_numeric(x, errors="coerce")
    out = out.mask(out < valid_min)
    out = out.mask(out > valid_max)
    out = out.mask(out.isin([-9, -7, -5, -4, -1]))
    return out


def read_hints_file(spec: dict[str, object], raw_dir: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    zip_path = raw_dir / str(spec["zip"])
    if not zip_path.exists():
        raise FileNotFoundError(f"Missing HINTS zip: {zip_path}")

    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        dta_name = str(spec["dta"])
        if dta_name not in names:
            dta_matches = [name for name in names if name.lower().endswith(".dta")]
            if len(dta_matches) != 1:
                raise FileNotFoundError(f"Could not identify a single .dta in {zip_path}")
            dta_name = dta_matches[0]
        with zf.open(dta_name) as fh:
            dat = pd.read_stata(fh, convert_categoricals=False)

    dat.columns = [col.lower() for col in dat.columns]
    required = {"age", "genderc", "bmi"}
    missing = sorted(required.difference(dat.columns))
    if missing:
        raise ValueError(f"{zip_path.name} is missing required columns: {missing}")

    age = numeric_clean(dat["age"], 40, 59)
    bmi = numeric_clean(dat["bmi"], 10, 80)
    sex = pd.to_numeric(dat["genderc"], errors="coerce")
    male = pd.Series(np.where(sex == 1, 1, np.where(sex == 2, 0, np.nan)), index=dat.index)

    out = pd.DataFrame(
        {
            "Y": np.nan,
            "x_age": age,
            "x_male": male,
            "x_bmi": bmi,
            "source": "HINTS",
            "year": int(spec["year"]),
            "cycle": str(spec["cycle"]),
        }
    )
    out = out.dropna(subset=["x_age", "x_male", "x_bmi"]).copy()

    diagnostics = {
        "cycle": str(spec["cycle"]),
        "year": int(spec["year"]),
        "raw_n": int(len(dat)),
        "complete_40_59_n": int(len(out)),
        "male_40_59": int((out["x_male"] == 1).sum()),
        "female_40_59": int((out["x_male"] == 0).sum()),
        "bmi_mean": float(out["x_bmi"].mean()),
        "bmi_min": float(out["x_bmi"].min()),
        "bmi_max": float(out["x_bmi"].max()),
    }
    return out, diagnostics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()

    root = args.root.resolve()
    raw_dir = root / "data" / "raw_dxa"
    out_dir = root / "data" / "processed_dxa"
    out_dir.mkdir(parents=True, exist_ok=True)

    frames: list[pd.DataFrame] = []
    diagnostics: list[dict[str, object]] = []
    for spec in HINTS_FILES:
        frame, diag = read_hints_file(spec, raw_dir)
        frames.append(frame)
        diagnostics.append(diag)

    hints = pd.concat(frames, ignore_index=True)
    hints_path = out_dir / "source_hints_2017_2018.csv"
    diag_path = out_dir / "source_hints_2017_2018_diagnostics.json"
    hints.to_csv(hints_path, index=False)

    total_diag = {
        "dataset": "HINTS source",
        "time_span": "2017+2018",
        "total_40_59": int(len(hints)),
        "male_40_59": int((hints["x_male"] == 1).sum()),
        "female_40_59": int((hints["x_male"] == 0).sum()),
        "cycles": diagnostics,
    }
    diag_path.write_text(json.dumps(total_diag, indent=2), encoding="utf-8")

    print(json.dumps(total_diag, indent=2))
    print(f"Wrote {hints_path}")
    print(f"Wrote {diag_path}")


if __name__ == "__main__":
    main()
