#!/usr/bin/env python3
import os
import sys
import argparse
import polars as pl

KINDS = ["cors", "sds"]
CHROMS = range(1, 23)

def fast_count_lines(path, buf_size = 16 * 1024 * 1024):
    n = 0
    with open(path, "rb") as f:
        while True:
            b = f.read(buf_size)
            if not b:
                break
            n += b.count(b"\n")
    return n

def list_subdirs(root):
    out = []
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if os.path.isdir(p):
            out.append(p)
    return out

def merge_kind(outdir, kind):
    
    p1 = os.path.join(outdir, f"c1_{kind}.tsv")
    if not os.path.exists(p1):
        raise FileNotFoundError(f"Missing file: {p1}")

    df1 = pl.read_csv(p1, separator="\t", null_values=["nan", "NA", ""])
    cols = df1.columns

    dfs = [df1]
    total_rows = max(0, fast_count_lines(p1) - 1)

    for c in range(2, 23):
        p = os.path.join(outdir, f"c{c}_{kind}.tsv")
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing file: {p}")

        d = pl.read_csv(
            p,
            separator="\t",
            skip_rows=1,
            has_header=False,
            null_values=["nan"],
        )
        d.columns = cols
        dfs.append(d)

        total_rows += max(0, fast_count_lines(p) - 1)

    merged = pl.concat(dfs, how="vertical")
    out_path = os.path.join(outdir, f"merged_mxp_{kind}.tsv")
    merged.write_csv(out_path, separator="\t", null_value = "nan")

    return total_rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mxp-root", required=True, help="Path to mxp/main")
    ap.add_argument("--bim", default=None, help="Path to autosomes.bim for row-count check")
    ap.add_argument("--task-id", type=int, default=None, help="1-based task index (optional; SLURM_ARRAY_TASK_ID otherwise)")
    ap.add_argument("--only", action="append",default=[],help="Only process these subdir")

    args = ap.parse_args()

    subdirs = list_subdirs(args.mxp_root)

    if args.only:
        wanted = set(args.only)
        subdirs = [p for p in subdirs if os.path.basename(p) in wanted]

    task_id = args.task_id
    if task_id is None:
        task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    if task_id <= 0 or task_id > len(subdirs):
        raise SystemExit(f"Bad task id {task_id}. Must be 1..{len(subdirs)}")

    outdir = subdirs[task_id - 1]
    print(f"[task {task_id}/{len(subdirs)}] merging: {outdir}")

    bim_rows = fast_count_lines(args.bim) if args.bim else None

    rows_cors = merge_kind(outdir, "cors")
    rows_sds  = merge_kind(outdir, "sds")

    if bim_rows is not None:
        if rows_cors != bim_rows:
            raise SystemExit(f"{outdir}: merged_mxp_cors.tsv has {rows_cors} rows, bim has {bim_rows}")
        if rows_sds != bim_rows:
            raise SystemExit(f"{outdir}: merged_mxp_sds.tsv has {rows_sds} rows, bim has {bim_rows}")

    print(f"[OK] {outdir} (cors={rows_cors}, sds={rows_sds})")

if __name__ == "__main__":
    main()
