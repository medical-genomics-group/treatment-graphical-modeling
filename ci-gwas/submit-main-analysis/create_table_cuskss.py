#!/usr/bin/env python3

import sys
import json
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy import stats


def load_bim(bim_path):
    return pd.read_csv(
        filepath_or_buffer=bim_path,
        sep="\t",
        names=["chrm", "rsid", "X1", "bp", "X2", "X3"],
        header=None,
    )


def load_traits_from_pxp(pxp_path):
    df = pd.read_csv(pxp_path, sep="\t", index_col=0)
    return df.index.to_list()


def selected_markers(
    bim_path,
    corr_path,
    adj_path,
    ixs_path,
    minz_path,
    traits,
    bim_df=None,
):
    num_phen = len(traits)
    if bim_df is None:
        bim_df = load_bim(bim_path)

    rs_ids = bim_df["rsid"].values
    chrs = bim_df["chrm"].values
    on_chr_positions = bim_df["bp"].values

    adj = mmread(adj_path).toarray()
    corr = mmread(corr_path).toarray()
    minz = mmread(minz_path).toarray()
    glob_ixs = np.fromfile(ixs_path, dtype=np.int32)

    assoc_markers = []
    for pix in np.arange(0, num_phen):
        mask = adj[pix, num_phen:]
        bim_lines = glob_ixs[np.where(mask)]
        corrs = corr[pix, num_phen:][np.where(mask)]
        minzs = minz[pix, num_phen:][np.where(mask)]
        for bim_line, c, z in zip(bim_lines, corrs, minzs):
            assoc_markers.append(
                {
                    "phenotype": traits[pix],
                    "rsID": rs_ids[bim_line],
                    "bim_line_ix": bim_line,
                    "chr": chrs[bim_line],
                    "bp": on_chr_positions[bim_line],
                    "corr": c,
                    "min_z": z,
                }
            )

    return pd.DataFrame(assoc_markers)


def FDR(x):
    """
    Copied from p.adjust function in R
    source: https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v: v[1], reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v: v[1])]
    q = sum([1.0 / i for i in range(1, len(x) + 1)])
    l = [q * len(x) / i * x[j] for i, j in zip(reversed(range(1, len(x) + 1)), o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l


def main():
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: postprocess_cusk.py <cusk_dir> <pxp_path>"
        )

    cusk_dir = sys.argv[1]
    pxp_path = sys.argv[2]

    bim_path = "path/to/bedfiles/autosomes.bim"
    corr_path = f"{cusk_dir}/cuskss_merged_scm.mtx"
    adj_path = f"{cusk_dir}/cuskss_merged_sam.mtx"
    minz_path = f"{cusk_dir}/cuskss_merged_szm.mtx"
    ixs_path = f"{cusk_dir}/cuskss_merged.ixs"

    traits = load_traits_from_pxp(pxp_path)

    sel_df = selected_markers(
        bim_path=bim_path,
        corr_path=corr_path,
        adj_path=adj_path,
        ixs_path=ixs_path,
        minz_path=minz_path,
        traits=traits,
    )

    setup = Path(pxp_path).parent.name

    n_dir = Path("path/to/mxp/main")
    
    p = n_dir / setup / "effective_n.json"
    with p.open("r") as f:
        effn = json.load(f)

    sel_df["n"] = sel_df["phenotype"].map(effn)

    t = sel_df["min_z"] * np.sqrt((sel_df["n"] - 3) / (1 - sel_df["min_z"]))
    sel_df["p_value"] = 2 * (1 - stats.norm.cdf(abs(t)))

    sel_df["p_fdr"] = FDR(sel_df["p_value"].to_numpy())

    sel_df.to_csv(
        f"{cusk_dir}/cuskss_selected_markers.tsv",
        sep="\t",
        na_rep="nan",
        index=False,
    )

    rsids = list(sel_df.sort_values("bim_line_ix")["rsID"].unique())
    with open(f"{cusk_dir}/cuskss.rsids", "w") as fout:
        for rsid in rsids:
            fout.write(f"{rsid}\n")


if __name__ == "__main__":
    main()
