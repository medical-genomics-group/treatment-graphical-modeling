#!/usr/bin/env python
"""Calculate marker-marker correlations within blocks"""
import sys
import subprocess
import pandas as pd

def main():
    # chromosome number
    chrom = sys.argv[1]
    # path to output of ci-gwas.py block
    blockfile = sys.argv[2]
    # path to filestem of plink bedfiles
    bfiles = sys.argv[3]
    outdir = sys.argv[4]
    # path to plink binary
    plinkpath = sys.argv[5]

    subprocess.run(f"mkdir -p {outdir}", shell=True, check=True)

    bimfile = f"{bfiles}.bim"
    blocks = pd.read_csv(
        blockfile, sep="\t", header=None, names=["chr", "first", "last"]
    )
    bim = pd.read_csv(
        filepath_or_buffer=bimfile,
        sep="\t",
        names=["chr", "rsid", "X1", "bp", "X2", "X3"],
        header=None,
    )

    block_ix = 0

    for first, last in zip(blocks["first"], blocks["last"]):
        block_ix += 1
        print(f"At block_ix: {block_ix}")
        first_rsid = bim[bim["chr"] == int(chrom)].iloc[int(first), 1]
        last_rsid = bim[bim["chr"] == int(chrom)].iloc[int(last), 1]
        subprocess.run(
            " ".join(
                [
                    plinkpath,
                    "--r triangle bin4",
                    "--memory 32000",
                    "--threads 4",
                    f"--bfile {bfiles}",
                    f"--from {first_rsid}",
                    f"--to {last_rsid}",
                    "--allow-no-sex",
                    f"--out {outdir}/c{chrom}_b{block_ix}",
                ]
            ),
            shell=True,
            check=True,
        )


if __name__ == "__main__":
    main()
