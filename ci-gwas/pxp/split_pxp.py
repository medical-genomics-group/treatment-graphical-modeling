#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

CI_GWAS_ROOT = Path('path/to/ci-gwas')
PXP_ROOT = CI_GWAS_ROOT / 'pxp'

COV_TOKENS = ['afs', 'als', 'yob', 'YOB', 'age', 'measno']
DRUG_TOKENS = [
    'ACE_inhibitor',
    'angiotensin_receptor_blocker',
    'calcium_channel_blocker',
    'beta_blocker',
    'diuretic',
    'statin',
]


def read_pxp(dataset):
    f = PXP_ROOT / dataset
    cor = pd.read_csv(f / 'cor_pxp.tsv', sep='\t', index_col=0)
    ses = pd.read_csv(f / 'ses_pxp.tsv', sep='\t', index_col=0)
    return cor, ses


def split_all_bp_pre_post():
    
    created = []
    
    for entry in sorted(PXP_ROOT.iterdir()):
        
        if not entry.is_dir():
            continue
        
        name = entry.name
        
        if not name.startswith('bp_pre_post'):
            continue

        cor, ses = read_pxp(name)
        vars_all = list(cor.columns)

        sbp_keep = [v for v in vars_all if 'DBP' not in v]
        dbp_keep = [v for v in vars_all if 'SBP' not in v]

        suffix = name[len('bp_'):]
        sbp_out = 'sbp_' + suffix
        dbp_out = 'dbp_' + suffix

        # SBP-only pxp
        sbp_cor = cor.loc[sbp_keep, sbp_keep]
        sbp_ses = ses.loc[sbp_keep, sbp_keep]
        sbp_dir = PXP_ROOT / sbp_out
        sbp_dir.mkdir(parents=True, exist_ok=True)
        sbp_cor.to_csv(sbp_dir / 'cor_pxp.tsv', sep='\t')
        sbp_ses.to_csv(sbp_dir / 'ses_pxp.tsv', sep='\t')

        # DBP-only pxp
        dbp_cor = cor.loc[dbp_keep, dbp_keep]
        dbp_ses = ses.loc[dbp_keep, dbp_keep]
        dbp_dir = PXP_ROOT / dbp_out
        dbp_dir.mkdir(parents=True, exist_ok=True)
        dbp_cor.to_csv(dbp_dir / 'cor_pxp.tsv', sep='\t')
        dbp_ses.to_csv(dbp_dir / 'ses_pxp.tsv', sep='\t')

        created.append((name, sbp_out, dbp_out))

    return created


def make_cvd_and_drugs():

    """ makes extra datasets pooled with ever or prior cvd
    then age-split with ever cvd then only ever and prior cvd then ever cvd with drugs """
    
    for base in (
        "sbp_pre_post_1to60",
        "dbp_pre_post_1to60",
        "dbp_pre_post_1to60_on_drugs",
        "sbp_pre_post_1to60_on_drugs",
        "ldl_pre_post_1to60"
    ):
        cor, ses = read_pxp(base)

        cor_ever = cor.drop(index="prior_cvd", columns="prior_cvd", errors="ignore")
        ses_ever = ses.drop(index="prior_cvd", columns="prior_cvd", errors="ignore")
        out_ever = base + "_with_ever_cvd"
        out_dir = PXP_ROOT / out_ever
        out_dir.mkdir(parents=True, exist_ok=True)
        cor_ever.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses_ever.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

        cor_prior = cor.drop(index="ever_cvd", columns="ever_cvd", errors="ignore")
        ses_prior = ses.drop(index="ever_cvd", columns="ever_cvd", errors="ignore")
        out_prior = base + "_with_prior_cvd"
        out_dir = PXP_ROOT / out_prior
        out_dir.mkdir(parents=True, exist_ok=True)
        cor_prior.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses_prior.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

    for base in (
        "sbp_pre_post_1to60_age5_with_statins",
        "dbp_pre_post_1to60_age5_with_statins",
    ):
        cor, ses = read_pxp(base)
        out_name = base.replace("_age5_with_statins", "_age5_with_prior_cvd")
        out_dir = PXP_ROOT / out_name
        out_dir.mkdir(parents=True, exist_ok=True)
        cor.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

    src = "bp_pre_post_1to60"
    cor_base, ses_base = read_pxp(src)
    cols = list(cor_base.columns)

    def is_covariate(col):
        return any(tok in col for tok in COV_TOKENS)

    def is_drug(col):
        return any(tok in col for tok in DRUG_TOKENS)

    keep_prior = [
        c for c in cols if c == "prior_cvd" or is_covariate(c)
    ]

    #prior cvd + covariates
    if keep_prior:
        cor = cor_base.loc[keep_prior, keep_prior]
        ses = ses_base.loc[keep_prior, keep_prior]
        out_dir = PXP_ROOT / "prior_cvd"
        out_dir.mkdir(parents=True, exist_ok=True)
        cor.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

    # ever_cvd + covariates
    keep_ever = [
        c for c in cols if c == "ever_cvd" or is_covariate(c)
    ]
    if keep_ever:
        cor = cor_base.loc[keep_ever, keep_ever]
        ses = ses_base.loc[keep_ever, keep_ever]
        out_dir = PXP_ROOT / "ever_cvd"
        out_dir.mkdir(parents=True, exist_ok=True)
        cor.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

    # ever_cvd + drugs + covariates
    keep_ever_drugs = [
        c
        for c in cols
        if (c == "ever_cvd") or is_covariate(c) or is_drug(c)
    ]
    if keep_ever_drugs:
        cor = cor_base.loc[keep_ever_drugs, keep_ever_drugs]
        ses = ses_base.loc[keep_ever_drugs, keep_ever_drugs]
        out_dir = PXP_ROOT / "ever_cvd_with_drugs"
        out_dir.mkdir(parents=True, exist_ok=True)
        cor.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses.to_csv(out_dir / "ses_pxp.tsv", sep="\t")


def make_no_cvd_versions():

    prefixes = ("bp", "sbp", "dbp", "ldl")

    for entry in sorted(PXP_ROOT.iterdir()):
        if not entry.is_dir():
            continue

        name = entry.name
        if not name.startswith(prefixes):
            continue
        if name.endswith("_no_cvd"):
            continue

        cor, ses = read_pxp(name)

        drop_labels = []
        for label in ["ever_cvd", "prior_cvd", "prior_cvd_1", "prior_cvd_2", "prior_cvd_3", "prior_cvd_4", "prior_cvd_5"]:
            if label in cor.columns:
                drop_labels.append(label)

        if not drop_labels:
            continue

        cor_nc = cor.drop(index=drop_labels, columns=drop_labels, errors="ignore")
        ses_nc = ses.drop(index=drop_labels, columns=drop_labels, errors="ignore")

        out_name = name + "_no_cvd"
        out_dir = PXP_ROOT / out_name
        out_dir.mkdir(parents=True, exist_ok=True)
        cor_nc.to_csv(out_dir / "cor_pxp.tsv", sep="\t")
        ses_nc.to_csv(out_dir / "ses_pxp.tsv", sep="\t")

