#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import ast
import numpy as np

drugs = [
    "ACE_inhibitor",
    "angiotensin_receptor_blocker",
    "calcium_channel_blocker",
    "beta_blocker",
    "diuretic",
    "statin",
]

def ph_path(out_root, run_id):
    return str(out_root / "regenie_phenos" / f"{run_id}.phen")

def cv_path(out_root, run_id):
    return str(out_root / "regenie_covars" / f"{run_id}.cov")

def base_covars_for_scope(scope):
    if scope in {"pooled", "interaction"}:
        return ["measno", "afs", "als", "yob", "age", "GP", "ASCEN"]

    if scope.startswith("age5g"):
        g = scope.replace("age5g", "")
        return [
            "afs", "als", "yob",
            f"GP_{g}", f"ASCEN_{g}", f"measno_{g}", f"age_{g}",
        ]


def covars_for_row(row):
    covars = []
    scope = row["scope"]

    covars.extend(base_covars_for_scope(scope))

    covars_in_model = row["covars_in_model"]

    def _age_suffix(col):
        if scope.startswith("age5g"):
            g = scope.replace("age5g", "")
            return f"{col}_{g}"
        return col

    if "+classes" in covars_in_model:
        covars.extend([_age_suffix(d) for d in drugs])

    if "+pre" in covars_in_model:
        trait = row["trait"]
        if trait in {"SBP", "DBP", "LDL"}:
            covars.append(_age_suffix(f"{trait}_pre"))

    if covars_in_model.startswith("base+others+"):
        pheno = row["phenotype"]

        if scope.startswith("age5g"):
            g = scope.replace("age5g", "")
            covars.extend([f"{d}_{g}" for d in drugs if f"{d}_{g}" != pheno])
        else:
            covars.extend([d for d in drugs if d != pheno])

        tail = covars_in_model.split("base+others+", 1)[1]
        covars.append(_age_suffix(tail))

    if scope == "interaction":
        covars.append(row["interaction_var"])

    return list(dict.fromkeys(covars))


def _add_bp_ldl_rows(add, scope, bp_tsv, ldl_tsv, g=None):
    suf = f"_{g}" if g is not None else ""

    for trait in ["SBP", "DBP"]:
        add(f"{scope}__{trait}_pre__base", scope, bp_tsv, trait, f"{trait}_pre{suf}", "base")
        add(f"{scope}__{trait}_post__base+classes+pre", scope, bp_tsv, trait, f"{trait}_post{suf}", "base+classes+pre")

        for cls in drugs:
            add(
                f"{scope}__CLASS_{cls}__base+others+{trait}_pre",
                scope, bp_tsv, f"CLASS_{cls}", f"{cls}{suf}", f"base+others+{trait}_pre",
            )

        for covars_in_model in ["base", "base+pre", "base+classes", "base+classes+pre"]:
            add(
                f"{scope}__{trait}_delta__{covars_in_model}",
                scope, bp_tsv, trait, f"{trait}_delta{suf}", covars_in_model,
            )

        for name in ["logdelta", "prop"]:
            for covars_in_model in ["base", "base+classes"]:
                add(
                    f"{scope}__{trait}_{name}__{covars_in_model}",
                    scope, bp_tsv, trait, f"{trait}_{name}{suf}", covars_in_model,
                )

    add(f"{scope}__LDL_pre__base", scope, ldl_tsv, "LDL", f"LDL_pre{suf}", "base")
    add(f"{scope}__LDL_post__base+classes+pre", scope, ldl_tsv, "LDL", f"LDL_post{suf}", "base+classes+pre")

    for cls in drugs:
        add(
            f"{scope}__CLASS_{cls}__base+others+LDL_pre",
            scope, ldl_tsv, f"CLASS_{cls}", f"{cls}{suf}", "base+others+LDL_pre",
        )

    for covars_in_model in ["base", "base+pre", "base+classes", "base+classes+pre"]:
        add(
            f"{scope}__LDL_delta__{covars_in_model}",
            scope, ldl_tsv, "LDL", f"LDL_delta{suf}", covars_in_model,
        )

    for name in ["logdelta", "prop"]:
        for covars_in_model in ["base", "base+classes"]:
            add(
                f"{scope}__LDL_{name}__{covars_in_model}",
                scope, ldl_tsv, "LDL", f"LDL_{name}{suf}", covars_in_model,
            )


def make_normal_inputs(out_root, bp_pooled, ldl_pooled, bp_age5, ldl_age5):
    rows = []

    def add(run_id, scope, input_tsv, trait, phenotype, covars_in_model):
        rows.append({
            "run_id": run_id,
            "scope": scope,
            "input_tsv": input_tsv,
            "trait": trait,
            "phenotype": phenotype,
            "covars_in_model": covars_in_model,
            "pheno_path": ph_path(out_root, run_id),
            "covar_path": cv_path(out_root, run_id),
        })

    _add_bp_ldl_rows(add, "pooled", bp_pooled, ldl_pooled)
    add("pooled__ever_cvd__base", "pooled", bp_pooled, "CVD", "ever_cvd", "base")
    add("pooled__ever_cvd__base+classes", "pooled", bp_pooled, "CVD", "ever_cvd", "base+classes")
    add("pooled__prior_cvd__base", "pooled", bp_pooled, "CVD", "prior_cvd", "base")

    for g in [1, 2, 3, 4, 5]:
        scope = f"age5g{g}"
        _add_bp_ldl_rows(add, scope, bp_age5, ldl_age5, g=g)
        add(f"{scope}__prior_cvd__base", scope, bp_age5, "CVD", f"prior_cvd_{g}", "base")

    df = pd.DataFrame(
        rows,
        columns=[
            "run_id",
            "scope",
            "input_tsv",
            "trait",
            "phenotype",
            "covars_in_model",
            "pheno_path",
            "covar_path",
        ],
    )
    df["covars"] = df.apply(covars_for_row, axis=1)
    return df


def make_int_inputs(out_root, bp_pooled, ldl_pooled):
    rows = []

    def add(run_id, input_tsv, trait, phenotype, covars_in_model, interaction_var):
        rows.append({
            "run_id": run_id,
            "scope": "interaction",
            "input_tsv": input_tsv,
            "trait": trait,
            "phenotype": phenotype,
            "covars_in_model": covars_in_model,
            "interaction_var": interaction_var,
            "pheno_path": ph_path(out_root, run_id),
            "covar_path": cv_path(out_root, run_id),
        })

    for trait, tsv in [("SBP", bp_pooled), ("DBP", bp_pooled), ("LDL", ldl_pooled)]:
        for cls in drugs:
            for covars_in_model in ["base+classes", "base+classes+pre"]:
                add(
                    f"intx__{trait}_delta__INT_{cls}__{covars_in_model}",
                    tsv,
                    trait,
                    f"{trait}_delta",
                    covars_in_model,
                    cls,
                )

            for name in ["logdelta", "prop"]:
                add(
                    f"intx__{trait}_{name}__INT_{cls}__base+classes",
                    tsv,
                    trait,
                    f"{trait}_{name}",
                    "base+classes",
                    cls,
                )

    df = pd.DataFrame(
        rows,
        columns=[
            "run_id",
            "scope",
            "input_tsv",
            "trait",
            "phenotype",
            "covars_in_model",
            "interaction_var",
            "pheno_path",
            "covar_path",
        ],
    )
    df["covars"] = df.apply(covars_for_row, axis=1)
    return df


# the above part prepared tables with input description the part below prepares the actual inputs

def compute_pheno(df, phenotype, scope):

    if phenotype in df.columns:
        return df[phenotype]

    def _split_suffix(ph, tag):
        if ph.endswith(tag):
            return ph[: -len(tag)], None
        if tag in ph:
            # e.g. SBP_delta_1
            parts = ph.split(tag)
            if len(parts) == 2 and parts[1].startswith("_"):
                base = parts[0]
                g = parts[1][1:]
                return base, g
        return None, None

    base, g = _split_suffix(phenotype, "_delta")
    if base is not None:
        if g is not None:
            return df[f"{base}_post_{g}"] - df[f"{base}_pre_{g}"]
        return df[f"{base}_post"] - df[f"{base}_pre"]

    base, g = _split_suffix(phenotype, "_prop")
    if base is not None:
        if g is not None:
            return df[f"{base}_post_{g}"] / df[f"{base}_pre_{g}"]
        return df[f"{base}_post"] / df[f"{base}_pre"]

    base, g = _split_suffix(phenotype, "_logdelta")
    if base is not None:
        if g is not None:
            return np.log(df[f"{base}_post_{g}"]) - np.log(df[f"{base}_pre_{g}"])
        return np.log(df[f"{base}_post"]) - np.log(df[f"{base}_pre"])

    return df[phenotype]


def prepare_regenie_inputs_from_manifest(
    manifest_tsv,
    out_root,
    ukb_covs,
    ukb_namecovs,
    fam
):

    out_root = Path(out_root)

    manifest = pd.read_csv(manifest_tsv, sep=",")
    manifest["covars"] = manifest["covars"].apply(lambda x: list(ast.literal_eval(x)))

    fam_ids = pd.read_csv(
        fam,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
    )

    namecov_path = Path(ukb_namecovs)
    names = [ln.strip() for ln in namecov_path.read_text().splitlines() if ln.strip()]

    ukb = pd.read_csv(
        ukb_covs,
        sep=r"\s+",
        header=None,
        names=names,
    )
    ukb = ukb[["FID", "IID"] + [c for c in ukb.columns if c not in {"FID", "IID"}]]

    pheno_cache = {}

    for _, row in manifest.iterrows():
        run_id = row["run_id"]
        input_tsv = row["input_tsv"]
        phenotype = row["phenotype"]
        covars_extra = row["covars"]

        if input_tsv not in pheno_cache:
            pheno_cache[input_tsv] = pd.read_csv(input_tsv, sep="\t")

        df_in = pheno_cache[input_tsv]

        scope = row["scope"]
        y = compute_pheno(df_in, phenotype, scope)

        ph = pd.DataFrame({"FID": df_in["eid"], "IID": df_in["eid"], phenotype: y})
        ph = fam_ids.merge(ph, on=["FID", "IID"], how="left")
        ph.to_csv(out_root / "regenie_phenos" / f"{run_id}.phen", sep="\t", index=False, na_rep = "NA")

        extra = df_in[["eid"] + covars_extra].rename(columns={"eid": "FID"}).assign(IID=lambda x: x["FID"])[["FID","IID"] + covars_extra]

        cov = fam_ids.merge(ukb, on=["FID", "IID"], how="left")
        cov = cov.merge(extra, on=["FID", "IID"], how="left")
        cov.to_csv(out_root / "regenie_covars" / f"{run_id}.cov", sep="\t", index=False, na_rep = "NA")
