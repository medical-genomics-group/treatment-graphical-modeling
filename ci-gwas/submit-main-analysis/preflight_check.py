#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
import pandas as pd


def upper_triangular_ok(df, atol=1e-12):
    a = df.to_numpy(dtype=float)
    # compare to upper triangle; allow NaNs to match NaNs
    return np.allclose(a, np.triu(a), atol=atol, rtol=0.0, equal_nan=True)


def read_time_index(path, n):

    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] == 1:
        v = df.iloc[:, 0].to_numpy()
    else:
        raise ValueError
  
    if len(v) < n:
        raise ValueError(f"time index too short: {len(v)} < {n}")
    return v[:n]


def check_inputs_for_run(gen_setup, pxp_path, mxp_dir, time_index_path, out_time_table):
    ok = True
    print(f"Checking inputs for {gen_setup}")


    pxp = pd.read_csv(pxp_path, sep="\t", index_col=0)
    print(f"PXP shape: {pxp.shape}")

    if upper_triangular_ok(pxp):
        print("OK: PXP is upper triangular")
    else:
        print("ERROR: PXP is not upper triangular")
        ok = False

    rows = pxp.index.astype(str).tolist()
    cols = pxp.columns.astype(str).tolist()

    if rows == cols:
        print("OK: PXP row and column names match and are in the same order.")
    else:
        print("ERROR: PXP row and column names DO NOT match.")
        ok = False

    expected_mxp_cols = ["chr", "snp", "ref"] + cols

    for c in range(1, 23):
        cors_path = os.path.join(mxp_dir, f"c{c}_cors.tsv")
        sds_path  = os.path.join(mxp_dir, f"c{c}_sds.tsv")

        if not os.path.exists(cors_path):
            print(f"ERROR: missing {cors_path}")
            ok = False
            continue
        if not os.path.exists(sds_path):
            print(f"ERROR: missing {sds_path}")
            ok = False
            continue

        mxp_cors = pd.read_csv(cors_path, sep="\t")
        mxp_sds  = pd.read_csv(sds_path,  sep="\t")

        cors_cols = mxp_cors.columns.astype(str).tolist()
        sds_cols  = mxp_sds.columns.astype(str).tolist()

        if cors_cols == sds_cols:
            print(f"OK: chr{c} mxp cors and sds have the same columns")
        else:
            print(f"ERROR: chr{c} mxp cors and sds columns differ")
            ok = False

        if cors_cols == expected_mxp_cols:
            print(f"OK: chr{c} mxp compatible with pxp")
        else:
            print(f"ERROR: chr{c} mxp not matching pxp")
            # tiny helpful debug
            print(f"  expected first cols: {expected_mxp_cols[:6]}")
            print(f"  observed first cols: {cors_cols[:6]}")
            ok = False

    if not os.path.exists(time_index_path):
        print(f"ERROR: missing time index file: {time_index_path}")
        ok = False
        time_idx = None
    else:
        try:
            time_idx = read_time_index(time_index_path, n=len(cols))
            print(f"OK: time index loaded ({len(time_idx)} entries)")
        except Exception as e:
            print(f"ERROR: failed reading time index: {e}")
            ok = False
            time_idx = None

    if time_idx is not None:
        out_df = pd.DataFrame(
            {
                "gen_setup_name": [gen_setup] * len(cols),
                "phenotype_name": cols,
                "time_index": time_idx,
            }
        )
        out_dir = os.path.dirname(out_time_table)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        out_df.to_csv(out_time_table, sep="\t", index=False)
        print(f"Wrote time index table: {out_time_table}")

    if not ok:
        print("Preflight FAILED.")
        sys.exit(1)

    print("Preflight OK.")
    sys.exit(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen-setup", required=True)
    ap.add_argument("--pxp", required=True, help="path to cor_pxp.tsv")
    ap.add_argument("--mxp-dir", required=True, help="dir containing c{chr}_cors.tsv and c{chr}_sds.tsv")
    ap.add_argument("--time-index", required=True, help="path to time index file for this GEN_SETUP")
    ap.add_argument("--out-time-table", required=True, help="output TSV: gen_setup_name phenotype_name time_index")
    args = ap.parse_args()

    check_inputs_for_run(
        gen_setup=args.gen_setup,
        pxp_path=args.pxp,
        mxp_dir=args.mxp_dir,
        time_index_path=args.time_index,
        out_time_table=args.out_time_table,
    )


if __name__ == "__main__":
    main()
