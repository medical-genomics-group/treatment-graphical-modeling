#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

CI_GWAS_ROOT = Path('path/to/ci-gwas')
MXP_ROOT = CI_GWAS_ROOT / 'mxp' / 'pruned'

COV_TOKENS = ['afs', 'als', 'yob', 'YOB', 'age', 'measno']
DRUG_TOKENS = [
    'ACE_inhibitor',
    'angiotensin_receptor_blocker',
    'calcium_channel_blocker',
    'beta_blocker',
    'diuretic',
    'statin',
]


def is_covariate(col):
    return any(tok in col for tok in COV_TOKENS)


def is_drug(col):
    return any(tok in col for tok in DRUG_TOKENS)


def split_all_bp_pre_post_mxp():

    for entry in sorted(MXP_ROOT.iterdir()):
        if not entry.is_dir():
            continue

        name = entry.name
        if not name.startswith("bp_pre_post"):
            continue

        src_folder = MXP_ROOT / name
        suffix = name[len("bp_"):]
        sbp_out = "sbp_" + suffix
        dbp_out = "dbp_" + suffix

        sbp_dir = MXP_ROOT / sbp_out
        dbp_dir = MXP_ROOT / dbp_out
        sbp_dir.mkdir(parents=True, exist_ok=True)
        dbp_dir.mkdir(parents=True, exist_ok=True)

        for ch in range(1, 23):
            cors_p = src_folder / f"c{ch}_cors.tsv"
            if not cors_p.exists():
                continue
            sds_p = src_folder / f"c{ch}_sds.tsv"

            dfc = pd.read_csv(cors_p, sep="\t")
            cols = list(dfc.columns)
            meta_cols = cols[:3]     
            pheno_cols = cols[3:]

            sbp_pheno = [c for c in pheno_cols if "DBP" not in c]
            dbp_pheno = [c for c in pheno_cols if "SBP" not in c]

            sbp_cols = meta_cols + sbp_pheno
            dbp_cols = meta_cols + dbp_pheno

            dfc[sbp_cols].to_csv(sbp_dir / cors_p.name, sep="\t", index=False, na_rep="nan")
            dfc[dbp_cols].to_csv(dbp_dir / cors_p.name, sep="\t", index=False, na_rep="nan")

            dfs = pd.read_csv(sds_p, sep="\t")
            dfs[sbp_cols].to_csv(sbp_dir / sds_p.name, sep="\t", index=False, na_rep="nan")
            dfs[dbp_cols].to_csv(dbp_dir / sds_p.name, sep="\t", index=False, na_rep="nan")


def make_cvd_and_drugs_mxp():

    for base in (
        "sbp_pre_post_1to60",
        "dbp_pre_post_1to60",
        "dbp_pre_post_1to60_on_drugs",
        "sbp_pre_post_1to60_on_drugs",
        "ldl_pre_post_1to60"
    ):
        src_dir = MXP_ROOT / base
        if not src_dir.exists():
            continue

        ever_dir = MXP_ROOT / (base + "_with_ever_cvd")
        prior_dir = MXP_ROOT / (base + "_with_prior_cvd")
        ever_dir.mkdir(parents=True, exist_ok=True)
        prior_dir.mkdir(parents=True, exist_ok=True)

        for ch in range(1, 23):
            cors_p = src_dir / f"c{ch}_cors.tsv"
            if not cors_p.exists():
                continue
            sds_p = src_dir / f"c{ch}_sds.tsv"

            dfc = pd.read_csv(cors_p, sep="\t")
            cols = list(dfc.columns)

            cols_ever = [c for c in cols if c != "prior_cvd"]
            dfc[cols_ever].to_csv(ever_dir / cors_p.name, sep="\t", index=False, na_rep="nan")

            cols_prior = [c for c in cols if c != "ever_cvd"]
            dfc[cols_prior].to_csv(prior_dir / cors_p.name, sep="\t", index=False, na_rep="nan")

            dfs = pd.read_csv(sds_p, sep="\t")
            dfs[cols_ever].to_csv(ever_dir / sds_p.name, sep="\t", index=False, na_rep="nan")
            dfs[cols_prior].to_csv(prior_dir / sds_p.name, sep="\t", index=False, na_rep="nan")

    for base in (
        "sbp_pre_post_1to60_age5_with_statins",
        "dbp_pre_post_1to60_age5_with_statins",
    ):
        src_dir = MXP_ROOT / base
        if not src_dir.exists():
            continue

        out_name = base.replace("_age5_with_statins", "_age5_with_prior_cvd")
        out_dir = MXP_ROOT / out_name
        out_dir.mkdir(parents=True, exist_ok=True)

        for ch in range(1, 23):
            cors_p = src_dir / f"c{ch}_cors.tsv"
            if not cors_p.exists():
                continue
            sds_p = src_dir / f"c{ch}_sds.tsv"

            dfc = pd.read_csv(cors_p, sep="\t")
            dfc.to_csv(out_dir / cors_p.name, sep="\t", index=False)

            dfs = pd.read_csv(sds_p, sep="\t")
            dfs.to_csv(out_dir / sds_p.name, sep="\t", index=False)

    src = "bp_pre_post_1to60"
    src_dir = MXP_ROOT / src
    if not src_dir.exists():
        return

    for out_name, mode in [
        ("prior_cvd", "prior"),
        ("ever_cvd", "ever"),
        ("ever_cvd_with_drugs", "ever_drugs"),
    ]:
        out_dir = MXP_ROOT / out_name
        out_dir.mkdir(parents=True, exist_ok=True)

        for ch in range(1, 23):
            cors_p = src_dir / f"c{ch}_cors.tsv"
            if not cors_p.exists():
                continue
            sds_p = src_dir / f"c{ch}_sds.tsv"

            dfc = pd.read_csv(cors_p, sep="\t")
            cols_all = list(dfc.columns)
            meta_cols = cols_all[:3]
            pheno_cols = cols_all[3:]

            if mode == "prior":
                kept_pheno = [
                    c for c in pheno_cols
                    if c == "prior_cvd" or is_covariate(c)
                ]
            elif mode == "ever":
                kept_pheno = [
                    c for c in pheno_cols
                    if c == "ever_cvd" or is_covariate(c)
                ]
            else: 
                kept_pheno = [
                    c for c in pheno_cols
                    if (c == "ever_cvd") or is_covariate(c) or is_drug(c)
                ]

            if not kept_pheno:
                continue

            keep_cols = meta_cols + kept_pheno

            dfc[keep_cols].to_csv(out_dir / cors_p.name, sep="\t", index=False, na_rep="nan")

            dfs = pd.read_csv(sds_p, sep="\t")
            dfs[keep_cols].to_csv(out_dir / sds_p.name, sep="\t", index=False, na_rep="nan")


def make_no_cvd_versions_mxp():

    prefixes = ("dbp_pre_post_1to60",
        "sbp_pre_post_1to60", "dbp_pre_post_1to60_age5_with_statins", "sbp_pre_post_1to60_age5_with_statins")

    for entry in sorted(MXP_ROOT.iterdir()):
        if not entry.is_dir():
            continue

        name = entry.name
        if not name.startswith(prefixes):
            continue
        if name.endswith("_no_cvd"):
            continue

        src_dir = MXP_ROOT / name
        out_name = name + "_no_cvd"
        out_dir = MXP_ROOT / out_name
        out_dir.mkdir(parents=True, exist_ok=True)

        for ch in range(1, 23):
            cors_p = src_dir / f"c{ch}_cors.tsv"
            if not cors_p.exists():
                continue
            sds_p = src_dir / f"c{ch}_sds.tsv"

            dfc = pd.read_csv(cors_p, sep="\t")
            drop_labels = [lab for lab in (
                "ever_cvd",
                "prior_cvd",
                "cvd",
                "cvd_1",
                "cvd_2",
                "cvd_3",
                "cvd_4",
                "cvd_5",
                "prior_cvd_1",
                "prior_cvd_2",
                "prior_cvd_3",
                "prior_cvd_4",
                "prior_cvd_5"
            ) if lab in dfc.columns]
            if drop_labels:
                keep_cols = [c for c in dfc.columns if c not in drop_labels]
            else:
                keep_cols = list(dfc.columns)

            dfc[keep_cols].to_csv(out_dir / cors_p.name, sep="\t", index=False, na_rep="nan")

            dfs = pd.read_csv(sds_p, sep="\t")
            dfs[keep_cols].to_csv(out_dir / sds_p.name, sep="\t", index=False, na_rep="nan")

        print(f"finished {cors_p.name}")
