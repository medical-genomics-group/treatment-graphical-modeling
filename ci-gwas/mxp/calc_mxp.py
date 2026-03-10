#!/usr/bin/env python3

import os, sys
import numpy as np
import pandas as pd
from bed_reader import open_bed
from bicorr import biserial_correlation

def read_bed_entries(path, sample_ix, start, stop):
    with open_bed(path) as bed:
        return bed.read(index=np.s_[sample_ix, start:stop])

def get_iids(fam_path):
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None)
    return fam.iloc[:, 1].astype(str).values

def read_pxp_order(pxp_root, run_name):
    cor_path = os.path.join(pxp_root, run_name, "cor_pxp.tsv")
    df = pd.read_csv(cor_path, sep="\t", index_col=0, nrows=1)
    return list(df.columns)

def main():
    phenodir = sys.argv[1]
    run      = sys.argv[2]
    chrom    = sys.argv[3]
    outdir   = sys.argv[4]

    bed_dir = "path/to/bedfiles"
    bed_path = f"{bed_dir}/c{chrom}.bed"
    bim_path = f"{bed_dir}/c{chrom}.bim"
    fam_path = f"{bed_dir}/c{chrom}.fam"

    pxp_root = "path/to/pxp"

    traits = read_pxp_order(pxp_root, run)
    print(f"[mxp] run={run} chr={chrom} traits={len(traits)}")

    fam_iids = get_iids(fam_path)  # str
    sample_ix = np.arange(len(fam_iids))
    iids = fam_iids

    T = pd.DataFrame(index=iids)
    missing_files = []
    zero_overlap = []

    for t in traits:
        f = os.path.join(phenodir, f"{run}__{t}.tsv")
        if not os.path.exists(f):
            T[t] = np.nan
            missing_files.append(t)
            continue
        df = pd.read_csv(f, sep="\t", dtype={"IID": str})
        s = pd.to_numeric(df[t], errors="coerce")
        s.index = df["IID"].astype(str)
        s = s.reindex(iids)
        if s.notna().sum() == 0:
            zero_overlap.append(t)
        T[t] = s

    if missing_files:
        print(f"[warn] {run}: {len(missing_files)} trait files missing (e.g. {missing_files[:5]})")
    if zero_overlap:
        print(f"[warn] {run}: {len(zero_overlap)} traits have 0 overlapping non-NA with FAM (e.g. {zero_overlap[:5]})")

    T.reset_index(inplace=True)
    T.rename(columns={"index": "IID"}, inplace=True)

    bim = pd.read_csv(bim_path, sep=r"\s+", header=None)
    out_cor = bim.iloc[:, [0, 1, 4]].copy()
    out_cor.columns = ["chr", "snp", "ref"]
    out_ses = out_cor.copy()
    for t in traits:
        out_cor[t] = np.nan
        out_ses[t] = np.nan

    num_snps = len(bim)
    step = 10000

    for start in range(0, num_snps, step):
        stop = min(start + step, num_snps)
        print(f"[chr{chrom}] {start}/{num_snps} ({start/num_snps:.2%})", flush=True)
        G = read_bed_entries(bed_path, sample_ix, start, stop)

        for local_ix in range(stop - start):
            g = G[:, local_ix]

            for t in (tt for tt in traits if not tt.endswith("ADJ")):
                r, s = biserial_correlation(g, T[t].values)
                out_cor.iat[start + local_ix, out_cor.columns.get_loc(t)] = r
                out_ses.iat[start + local_ix, out_ses.columns.get_loc(t)] = s

            for t in (tt for tt in traits if tt.endswith("ADJ")):
                y = T[t].values
                ok = ~np.isnan(g) & ~np.isnan(y)
                if ok.sum() < 3:
                    continue
                x = g[ok]; yy = y[ok]
                if np.std(x) == 0 or np.std(yy) == 0:
                    continue
                r = np.corrcoef(x, yy)[0, 1]
                n = len(x)
                s = (1 - r**2) / np.sqrt(max(n - 3, 1))
                out_cor.iat[start + local_ix, out_cor.columns.get_loc(t)] = r
                out_ses.iat[start + local_ix, out_ses.columns.get_loc(t)] = s

    os.makedirs(outdir, exist_ok=True)
    out_cor.to_csv(os.path.join(outdir, f"c{chrom}_cors.tsv"), sep="\t", index=False, na_rep="nan")
    out_ses.to_csv(os.path.join(outdir, f"c{chrom}_sds.tsv"),  sep="\t", index=False, na_rep="nan")
    print(f"[done] {run} chr{chrom}")

if __name__ == "__main__":
    main()
