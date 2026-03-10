#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import subprocess
from scipy.io import mmread


def load_bim(bim_path):
    return pd.read_csv(
        filepath_or_buffer=bim_path,
        sep="\t",
        names=["chrm", "rsid", "X1", "bp", "X2", "X3"],
        header=None,
    )


def load_traits(traits_path, n_phen):
    if traits_path is None or not os.path.exists(traits_path):
        return [f"phenotype_{i}" for i in range(n_phen)]
    with open(traits_path, "r") as fin:
        traits = [l.strip() for l in fin if l.strip()]
    if len(traits) >= n_phen:
        return traits[:n_phen]
    return [f"phenotype_{i}" for i in range(n_phen)]


def selected_markers(
    bim_path,
    corr_path,
    adj_path,
    ixs_path,
    n_phen,
    traits_path=None,
    bim_df=None,
):
    if bim_df is None:
        bim_df = load_bim(bim_path)

    rs_ids = bim_df["rsid"].values
    chrs = bim_df["chrm"].values
    on_chr_positions = bim_df["bp"].values

    adj = mmread(adj_path).toarray()
    corr = mmread(corr_path).toarray()
    glob_ixs = np.fromfile(ixs_path, dtype=np.int32)

    if adj.shape[0] != adj.shape[1]:
        raise ValueError(f"adj matrix not square: {adj.shape}")

    n_total = adj.shape[0]
    n_markers = glob_ixs.shape[0]

    if n_total != n_phen + n_markers:
        raise ValueError(
            f"dimension mismatch: adj={n_total}, "
            f"pxp phen={n_phen}, ixs markers={n_markers}"
        )

    p_names = load_traits(traits_path, n_phen)
    assoc_markers = []

    for pix in range(n_phen):
        row = adj[pix, n_phen:]          
        mask = row != 0
        if not np.any(mask):
            continue

        bim_lines = glob_ixs[mask]       
        corrs = corr[pix, n_phen:][mask]

        for bim_line, c in zip(bim_lines, corrs):
            assoc_markers.append(
                {
                    "phenotype": p_names[pix],
                    "rsID": rs_ids[bim_line],
                    "bim_line_ix": int(bim_line),
                    "chr": int(chrs[bim_line]),
                    "bp": int(on_chr_positions[bim_line]),
                    "corr": float(c),
                }
            )

    return pd.DataFrame(assoc_markers)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: calc_mxm_merged.py CUSK_DIR PXP_PATH\n"
            "  CUSK_DIR  directory with merged_blocks_* files\n"
            "  PXP_PATH  path to cor_pxp.tsv used in the run"
        )

    cusk_dir = sys.argv[1]
    pxp_path = sys.argv[2]

    pxp_df = pd.read_csv(pxp_path, sep="\t", index_col=0)
    n_phen = pxp_df.shape[1]

    bim_path = "path/to/bedfiles/autosomes.bim"
    bfile = "path/to/bedfiles/autosomes"

    corr_path = os.path.join(cusk_dir, "merged_blocks_scm.mtx")
    adj_path = os.path.join(cusk_dir, "merged_blocks_sam.mtx")
    ixs_path = os.path.join(cusk_dir, "merged_blocks.ixs")
    traits_path = os.path.join(cusk_dir, "traits.txt")

    plinkpath = "path/to/plink"
    outfile_prefix = os.path.join(cusk_dir, "merged_blocks.mxm")

    sel_df = selected_markers(
        bim_path=bim_path,
        corr_path=corr_path,
        adj_path=adj_path,
        ixs_path=ixs_path,
        n_phen=n_phen,
        traits_path=traits_path,
    )

    sel_df.to_csv(
        os.path.join(cusk_dir, "merged_blocks_selected_markers.tsv"),
        sep="\t",
        na_rep="nan",
        index=False,
    )

    rsids_path = os.path.join(cusk_dir, "merged_blocks.rsids")

    if sel_df.empty:
        open(rsids_path, "w").close()
        return

    rsids = list(sel_df.sort_values("bim_line_ix")["rsID"].unique())
    with open(rsids_path, "w") as fout:
        for rsid in rsids:
            fout.write(f"{rsid}\n")

    cmd = " ".join(
        [
            plinkpath,
            "--r triangle bin4",
            "--memory 10000",
            "--threads 4",
            f"--bfile {bfile}",
            "--allow-no-sex",
            f"--extract {rsids_path}",
            f"--out {outfile_prefix}",
        ]
    )
    subprocess.run(cmd, shell=True, check=True)


if __name__ == "__main__":
    main()
