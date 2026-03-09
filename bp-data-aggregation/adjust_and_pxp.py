#!/usr/bin/env python3

import os
import re
import warnings
import numpy as np
import pandas as pd
from bicorr import biserial_correlation, tetrachoric_correlation
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression



BIN_BASES = [
    'ACE_inhibitor','angiotensin_receptor_blocker',
    'calcium_channel_blocker','beta_blocker','diuretic',
    'statin','ever_cvd','prior_cvd'
]

BIN_GLOBAL_UNSUFFIXED = {"ever_cvd"}

CONT_GLOBAL_UNSUFFIXED = {"afs", "als", "yob", "YOB"}

CONT_BASES = [
    'SBP_pre','DBP_pre','SBP_post','DBP_post','LDL_pre','LDL_post',
    'measno','age','afs','als','yob','YOB'
]

PHENO_FILES = [
    "bp_pre_post_1to60_with_copying.tsv",
    "bp_pre_post_1to30.tsv",
    "bp_pre_post_1to30_with_copying.tsv",
    "bp_pre_post_1to60.tsv",
    "bp_pre_post_1to60_age4_with_statins.tsv",
    "bp_pre_post_1to60_age5_no_statins.tsv",
    "bp_pre_post_1to60_age5_shifted_with_statins.tsv",
    "bp_pre_post_1to60_age5_with_statins.tsv",
    "bp_pre_post_1to60_no_statins.tsv",
    "bp_pre_post_1to60_on_drugs.tsv",
    "bp_pre_post_1to60_rel_up_to_2nd.tsv",
    "bp_pre_post_1to60_rel_up_to_3rd.tsv",
    "ldl_pre_post_1to30.tsv",
    "ldl_pre_post_1to60.tsv",
    "ldl_pre_post_1to60_age5.tsv",
    "ldl_pre_post_1to60_on_drugs.tsv",
    "ldl_pre_post_1to60_statins_only.tsv",
]

def read_pheno(path):
    df = pd.read_csv(path, sep="\t")
    df["eid"] = df["eid"].astype(str)
    df["FID"] = df["eid"]
    df["IID"] = df["eid"]
    df = df.drop(columns=["eid"])
    cols = ["FID", "IID"] + [c for c in df.columns if c not in ("FID","IID")]
    return df.loc[:, cols]

def read_covars(cov_path):
    name_path = cov_path.replace(".cov", ".Namecov")
    names = [line.strip() for line in open(name_path, "r") if line.strip() != ""]
    cv = pd.read_csv(cov_path, sep=r"\s+", header=None, engine="python")
    cv = cv.iloc[:, :len(names)]
    cv.columns = names[:cv.shape[1]]
    cv["FID"] = cv["FID"].astype(str)
    cv["IID"] = cv["IID"].astype(str)
    cov_cols = [c for c in cv.columns if c not in ("FID","IID")]
    if cov_cols:
        cv[cov_cols] = cv[cov_cols].apply(pd.to_numeric, errors="coerce")
    cols = ["FID","IID"] + cov_cols
    return cv.loc[:, cols]

def group_suffix(name):
    m = re.match(r"^(.+?)_([0-9]+)$", str(name))
    if m:
        return m.group(1), m.group(2)
    return name, None

def find_groups_from_df_and_name(df, stem):
    if "age" not in stem:
        return []
    groups = []
    for c in df.columns:
        m = re.match(r"^age_([0-9]+)$", str(c))
        if m:
            groups.append(m.group(1))
    return sorted(set(groups), key=lambda x: int(x))

def actual_trait_columns(ph, stem):
    groups = find_groups_from_df_and_name(ph, stem)
    
    if not groups:
        
        bin_cols = [b for b in BIN_BASES if b in ph.columns]
        cont_cols = [b for b in CONT_BASES if b in ph.columns]
        
    else:
        bin_cols = []
        
        for b in BIN_BASES:  
            if b in BIN_GLOBAL_UNSUFFIXED:
                if b in ph.columns:
                    bin_cols.append(b)
            else:
                for g in groups:
                    col = f"{b}_{g}"
                    if col in ph.columns:
                        bin_cols.append(col)

        cont_cols = []
        
        for b in CONT_BASES:
            if b in CONT_GLOBAL_UNSUFFIXED:
                if b in ph.columns:
                    cont_cols.append(b)
            else:
                for g in groups:
                    col = f"{b}_{g}"
                    if col in ph.columns:
                        cont_cols.append(col)
                        
    return bin_cols, cont_cols

def design_matrix(df, covars):
    use = [c for c in covars if c in df.columns]
    cols = [np.ones(len(df))]
    for c in use:
        cols.append(pd.to_numeric(df[c], errors="coerce").to_numpy())
    if len(cols) == 1:
        return np.ones((len(df), 1))
    return np.column_stack(cols)

def ols_residual(y, X):
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)

    mask = ~np.isnan(y)

    out = np.full_like(y, np.nan, dtype=float)
    
    X_use = X[mask, :]
    y_use = y[mask]

    lm = LinearRegression(fit_intercept=False)
    lm.fit(X_use, y_use)
    y_pred = lm.predict(X_use)

    out[mask] = y_use - y_pred
    return out

def zscore(x):
    x = np.asarray(x, dtype=float)

    mu = np.nanmean(x)
    sd = np.nanstd(x, ddof=1)

    return (x - mu) / sd

def safe_tetrachoric(a, b, name_a, name_b):
    """ if a correlation can't compute due to lack of overlap it is set to zero,
    if it can't compute for cvd (as people can't gen undiagnosed it is set to 1 """
    
    x = pd.to_numeric(a, errors="coerce").to_numpy()
    y = pd.to_numeric(b, errors="coerce").to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]; y = y[valid]

    pat = r"^prior_cvd_([0-9]+)$"
    ma = re.match(pat, str(name_a))
    mb = re.match(pat, str(name_b))
    prior_pair = (ma is not None) and (mb is not None) and (ma.group(1) != mb.group(1))

    n11 = np.sum((x == 1) & (y == 1))
    n10 = np.sum((x == 1) & (y == 0))
    n01 = np.sum((x == 0) & (y == 1))
    n00 = np.sum((x == 0) & (y == 0))

    degenerate = (n11 == 0) or (n10 == 0) or (n01 == 0) or (n00 == 0)

    if degenerate:
        if prior_pair:
            print(f"[warn] tetrachoric {name_a} x {name_b}: degenerate prior_cvd age groups, set r=1, se=0")
            return 1.0, 0.0
        else:
            print(f"[warn] tetrachoric {name_a} x {name_b}: degenerate 2x2 (n11={n11}, n10={n10}, n01={n01}, n00={n00}) r=0, se=0")
            return 0.0, 0.0
            
    else:
        return tetrachoric_correlation(x, y)

def safe_biserial(cont, b, name_c, name_b):
    x = pd.to_numeric(cont, errors="coerce").to_numpy()
    y = pd.to_numeric(b, errors="coerce").to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]; y = y[valid]
    has0 = np.any(y == 0); has1 = np.any(y == 1)
    if (not has0) or (not has1):
        print(f"[warn] biserial {name_c}×{name_b}: missing class (has0={has0}, has1={has1})  r=0, se=0")
        return 0.0, 0.0
    if np.nanstd(x, ddof=1) == 0:
        print(f"[warn] biserial {name_c}×{name_b}: zero SD in continuous r=0, se=0")
        return 0.0, 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        return biserial_correlation(x, y)

def adjust_file(pheno_path, covar_path, fam):
    stem = os.path.basename(pheno_path)
    if stem.lower().endswith(".tsv"):
        stem = stem[:-4]
    print(f"[adjust] {stem}")

    ph = read_pheno(pheno_path)
    cv = read_covars(covar_path)

    ph = ph[ph["IID"].isin(fam["IID"])]
    cv = cv[cv["IID"].isin(ph["IID"])]

    merged = ph.merge(cv, on=["FID", "IID"], how="left")

    general_covs = [c for c in cv.columns if c not in ("FID", "IID")]

    bin_cols, cont_cols = actual_trait_columns(ph, stem)

    groups = find_groups_from_df_and_name(ph, stem)
    is_age_grouped = len(groups) > 0

    out = ph[["FID", "IID"] + bin_cols].copy()

    for t in cont_cols:
        base_covs = list(general_covs)
        extra_covs = []

        base, grp = group_suffix(t)

        if is_age_grouped:
            if grp is not None:    
                if f"ASCEN_{grp}" in ph.columns:
                    extra_covs.append(f"ASCEN_{grp}")
                if f"GP_{grp}" in ph.columns:
                    extra_covs.append(f"GP_{grp}")  
            else:
                pass
        else:
            
            if "GP" in ph.columns:
                extra_covs.append("GP")
            if "ASCEN" in ph.columns:
                extra_covs.append("ASCEN")

        covs = base_covs + extra_covs

        X = design_matrix(merged, covs)
        y = pd.to_numeric(ph[t], errors="coerce").astype(float).values
        out[t + "_ADJ"] = zscore(ols_residual(y, X))

    adj_path = os.path.join(os.path.dirname(pheno_path), f"{stem}_ADJ.tsv")
    out.to_csv(adj_path, sep="\t", index=False, na_rep="nan")
    
    return adj_path, stem

def pearson_se(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan, np.nan
    xv, yv = x[mask], y[mask]
    if xv.std(ddof=1) == 0 or yv.std(ddof=1) == 0:
        return np.nan, np.nan
    r = np.corrcoef(xv, yv)[0,1]
    n = len(xv)
    se = (1 - r**2) / np.sqrt(max(n - 3, 1))
    return float(r), float(se)

def pxp_from_adjusted(adj_path, out_pxp_root, stem):
    df = pd.read_csv(adj_path, sep="\t")
    cols = [c for c in df.columns if c not in ("FID","IID")]
    cont_adj = [c for c in cols if c.endswith("ADJ")]
    binary = [c for c in cols if c not in cont_adj]
    traits = sorted(binary + cont_adj)
    cor = pd.DataFrame(np.zeros((len(traits), len(traits))), index=traits, columns=traits, dtype=float)
    se = pd.DataFrame(np.zeros_like(cor.values), index=traits, columns=traits, dtype=float)
    for i, ti in enumerate(traits):
        for j, tj in enumerate(traits):
            if i == j:
                r, s = 1.0, 0.0
            else:
                ic = ti.endswith("ADJ")
                jc = tj.endswith("ADJ")
                if ic and jc:
                    x = pd.to_numeric(df[ti], errors="coerce").astype(float).values
                    y = pd.to_numeric(df[tj], errors="coerce").astype(float).values
                    r, s = pearson_se(x, y)
                elif (not ic) and (not jc):
                    r, s = safe_tetrachoric(df[ti], df[tj], ti, tj)
                elif ic:
                    r, s = safe_biserial(df[ti], df[tj], ti, tj)
                else:
                    r, s = safe_biserial(df[tj], df[ti], tj, ti)
            cor.iloc[i, j] = r
            se.iloc[i, j] = s
    m = cor.values.copy()
    svals = se.values.copy()
    ii = np.tril_indices_from(m, k=0)
    m[ii] = 0.0
    svals[ii] = 0.0
    cor_nan = int(np.isnan(m).sum()); cor_inf = int(np.isinf(m).sum())
    se_nan  = int(np.isnan(svals).sum()); se_inf  = int(np.isinf(svals).sum())
    print(f"[pxp_scan] {stem}: cor nan={cor_nan} inf={cor_inf}; se nan={se_nan} inf={se_inf}")
    outdir = os.path.join(out_pxp_root, stem)
    os.makedirs(outdir, exist_ok=True)
    pd.DataFrame(m, index=traits, columns=traits).to_csv(os.path.join(outdir, "cor_pxp.tsv"), sep="\t", na_rep="nan")
    pd.DataFrame(svals, index=traits, columns=traits).to_csv(os.path.join(outdir, "ses_pxp.tsv"), sep="\t", na_rep="nan")
    print(f"[pxp] {stem} -> {os.path.join(outdir, 'cor_pxp.tsv')}")
    try:
        fig_w = min(max(6.0, 0.30 * len(traits)), 24.0)
        fig_h = min(max(5.0, 0.30 * len(traits)), 24.0)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(m, vmin=-1.0, vmax=1.0, cmap="bwr", interpolation="nearest")
        ax.set_xticks(range(len(traits)))
        ax.set_yticks(range(len(traits)))
        ax.set_xticklabels(traits, rotation=90, fontsize=7)
        ax.set_yticklabels(traits, fontsize=7)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=8)
        ax.set_title(f"P×P correlations: {stem}", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "cor_pxp_heatmap.png"), dpi=200)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] heatmap {stem}: {e}")

def export_adj_individual_phenos(base_dir, out_dir, fam_path):
    print("[fam] reading skeleton")
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None, engine="python")
    fam = fam.rename(columns={0: "FID", 1: "IID"})
    fam["FID"] = fam["FID"].astype(str)
    fam["IID"] = fam["IID"].astype(str)
    fam = fam[["FID","IID"]]
    os.makedirs(out_dir, exist_ok=True)
    for name in PHENO_FILES:
        stem = name[:-4] if name.endswith(".tsv") else name
        adj_path = os.path.join(base_dir, f"{stem}_ADJ.tsv")
        if not os.path.exists(adj_path):
            continue
        print(f"[export] {stem}")
        df = pd.read_csv(adj_path, sep="\t")
        df["FID"] = df["FID"].astype(str)
        df["IID"] = df["IID"].astype(str)
        trait_cols = [c for c in df.columns if c not in ("FID","IID")]
        for t in trait_cols:
            tmp = df[["FID","IID", t]]
            merged = fam.merge(tmp, on=["FID","IID"], how="left")
            out_path = os.path.join(out_dir, f"{stem}__{t}.tsv")
            merged.to_csv(out_path, sep="\t", index=False, na_rep="NA")

def process_all(base_dir, covar_file, out_pxp_root, fam_path=""):
    fam = None
    if fam_path:
        print("[fam] reading for adjustment")
        fam = pd.read_csv(fam_path, sep=r"\s+", header=None, engine="python")
        fam = fam.rename(columns={0: "FID", 1: "IID"})
        fam["FID"] = fam["FID"].astype(str)
        fam["IID"] = fam["IID"].astype(str)
        fam = fam[["FID","IID"]]
    rows = []
    for name in PHENO_FILES:
        ph = os.path.join(base_dir, name)
        if not os.path.exists(ph):
            print(f"[skip] {name} (missing)")
            continue
        adj_path, stem = adjust_file(ph, covar_file, fam=fam)
        pxp_from_adjusted(adj_path, out_pxp_root, stem)
        rows.append({"input": ph, "adjusted": adj_path, "pxp_dir": os.path.join(out_pxp_root, stem)})
    man = pd.DataFrame(rows)
    if rows:
        os.makedirs(out_pxp_root, exist_ok=True)
        man.to_csv(os.path.join(out_pxp_root, "manifest.tsv"), sep="\t", index=False)
        print("[done] manifest written")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="CI-GWAS single-use prep")
    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--covar-file", required=True)
    ap.add_argument("--out-pxp-root", required=True)
    ap.add_argument("--fam-path", default="")
    ap.add_argument("--out-ind-dir", default="")
    args = ap.parse_args()
    process_all(args.base_dir, args.covar_file, args.out_pxp_root, args.fam_path)
    if args.fam_path and args.out_ind_dir:
        export_adj_individual_phenos(args.base_dir, args.out_ind_dir, args.fam_path)
