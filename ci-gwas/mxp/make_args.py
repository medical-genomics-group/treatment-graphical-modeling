#!/usr/bin/env python3
import os, re

def main():
    phenodir = "path/to/phenotype/files"
    outroot  = "path/to/mxp/output"

    runs = set()
    for fn in os.listdir(phenodir):
        if not fn.endswith(".tsv"):
            continue
        m = re.match(r"(.+?)__.+\.tsv$", fn)
        if m:
            runs.add(m.group(1))
    runs = sorted(runs)

    os.makedirs("log", exist_ok=True)
    n = 0
    with open("args.txt", "w") as f:
        for run in runs:
            outdir = os.path.join(outroot, run)
            for c in range(1, 23):
                f.write(f"{phenodir} {run} {c} {outdir}\n")
                n += 1
    print(f"[done] discovered {len(runs)} runs; wrote {n} lines to args.txt")

if __name__ == "__main__":
    main()
