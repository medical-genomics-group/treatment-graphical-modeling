#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
import sys

import polars as pl


KEY = ["chr", "snp", "ref"]
NULLS = ["nan", "NaN", "NA", ""]


def compute_effective_n(setup_dir):
    sds_path = setup_dir / "merged_mxp_sds.tsv"
    cors_path = setup_dir / "merged_mxp_cors.tsv"

    sds = pl.scan_csv(str(sds_path), separator="\t", null_values=NULLS)
    cors = pl.scan_csv(str(cors_path), separator="\t", null_values=NULLS)

    cols = [c for c in sds.collect_schema().names() if c not in KEY]
    if not cols:
        return {}

    sds_num = sds.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in cols])
    cors_num = cors.with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in cols])

    r_med = cors_num.select([pl.col(c).median().alias(c) for c in cols])
    s_med = sds_num.select([pl.col(c).median().alias(c) for c in cols])

    n_df = r_med.join(s_med, how="cross", suffix="_s").select(
        [(((1 - pl.col(c) ** 2) / pl.col(f"{c}_s")) ** 2).alias(c) for c in cols]
    )

    row = n_df.collect().row(0, named=True)

    out = {}
    for k, v in row.items():
        if v is None or v != v: 
            out[k] = None
        else:
            out[k] = float(v)
    return out


def iter_setup_dirs(root):
    for p in sorted(root.iterdir()):
        if p.is_dir():
            yield p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        type=Path,
        default=Path("path/to/mxp/main"),
        help="Main mxp folder containing setup subfolders",
    )
    ap.add_argument(
        "--out-name",
        default="effective_n.json",
        help="Output JSON filename to write inside each setup folder",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output JSON if present",
    )
    args = ap.parse_args()

    root = args.root
    if not root.exists():
        print(f"ERROR: root does not exist: {root}", file=sys.stderr)
        return 2

    n_ok = 0
    n_skip = 0
    n_fail = 0

    for setup_dir in iter_setup_dirs(root):
        sds_path = setup_dir / "merged_mxp_sds.tsv"
        cors_path = setup_dir / "merged_mxp_cors.tsv"
        out_path = setup_dir / args.out_name

        if not (sds_path.exists() and cors_path.exists()):
            print(f"SKIP {setup_dir.name}: missing merged_mxp_sds.tsv or merged_mxp_cors.tsv")
            n_skip += 1
            continue

        if out_path.exists() and not args.overwrite:
            print(f"SKIP {setup_dir.name}: {out_path.name} exists (use --overwrite)")
            n_skip += 1
            continue

        try:
            out = compute_effective_n(setup_dir)
            out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
            print(f"OK   {setup_dir.name}: wrote {out_path}")
            n_ok += 1
        except Exception as e:
            print(f"FAIL {setup_dir.name}: {e}", file=sys.stderr)
            n_fail += 1

    print(f"\nDone. OK={n_ok}  SKIP={n_skip}  FAIL={n_fail}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
