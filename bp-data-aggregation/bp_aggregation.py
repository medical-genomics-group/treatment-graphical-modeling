import os
import numpy as np
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM

DRUG_CLASSES = ""
YOB_CSV      = ""
CVD_TSV      = ""
AFS_ALS_TSV  = ""

BASE_BP_COLS  = {"eid","date","age","systolic","diastolic","GP","ASCEN"}
BASE_LDL_COLS = {"eid","date","age","ldl","GP","ASCEN"}
EXTRA_NON_DRUG = {
    "row_id","yob","cvd","afs","als","measno","median_age",
    "drug","drug_cumsum","days_to_previous","date_cumsum","starts_with_0", "prior_cvd", "ever_cvd"
}

AGE_MAP_5_DEFAULT = {
    1: (-np.inf, 49.10922669),
    2: (49.10922669, 55.07205019),
    3: (55.07205019, 60.07187883),
    4: (60.07187883, 64.00752076),
    5: (64.00752076, np.inf),
}

CLASS_ORDER = [
    "ACE_inhibitor",
    "angiotensin_receptor_blocker",
    "calcium_channel_blocker",
    "beta_blocker",
    "diuretic",
    "statin",
]

def detect_substance_cols(df, base_cols):
    non_drug = set(base_cols) | EXTRA_NON_DRUG
    return [c for c in df.columns if c not in non_drug]

def binarize_by_window(df_raw, lo, hi, base_cols):
    df = df_raw.copy()
    df["eid"]  = df["eid"].astype(str)
    df["date"] = pd.to_numeric(df["date"], errors="coerce")

    subs = detect_substance_cols(df, base_cols)
    if subs:
        df[subs] = df[subs].apply(pd.to_numeric, errors="coerce")
        bin_block = (df[subs].ge(lo) & df[subs].le(hi)).astype("Int8")
        out = pd.concat([df[list(base_cols)], bin_block], axis=1)
    else:
        out = df[list(base_cols)].copy()

    for c in ("GP","ASCEN"):
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype("Int8")
    if "age" in out.columns:
        out["age"] = pd.to_numeric(out["age"], errors="coerce")
    return out, subs

def map_to_classes(df_bin, substance_names):
    drug = pd.read_csv(DRUG_CLASSES, sep="\t", dtype=str)
    drug["substance"] = drug["substance"].astype(str)
    drug["class"]     = drug["class"].astype(str)
    drug["class"] = drug["class"].replace({
        "Angiotensin_II_receptor_blocker": "angiotensin_receptor_blocker",
        "beta_blockers": "beta_blocker",
        "diuretics": "diuretic"
    })

    present = set(substance_names)
    class_to_subs = {
        cl: sorted(set(drug.loc[drug["class"] == cl, "substance"]).intersection(present))
        for cl in drug["class"].unique()
    }

    out = df_bin.copy()
    cls_cols = []
    for cl, subs in class_to_subs.items():
        if not subs:
            continue
        out[cl] = (out[subs].sum(axis=1) > 0).astype("Int8")
        cls_cols.append(cl)

    drop_cols = [c for c in substance_names if c in out.columns]
    if drop_cols:
        out.drop(columns=drop_cols, inplace=True, errors="ignore")
    return out, cls_cols

def add_yob(df):
    yob = pd.read_csv(YOB_CSV, dtype={"eid": str})
    df["eid"] = df["eid"].astype(str)
    df["yob"] = df["eid"].map(yob.set_index("eid")["YOB"])


def filter_by_year(df):
    df["yob"] = pd.to_numeric(df["yob"], errors="coerce")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    keep = (df["yob"] + df["age"]) <= 2021
    df.drop(index=df.index[~keep], inplace=True)

def add_yob_and_filter_year(df):
    add_yob(df)
    filter_by_year(df)

def add_ever_cvd(df):
    cvd = pd.read_csv(CVD_TSV, sep="\t", dtype={"eid": str})
    cvd_eids = set(cvd["eid"].unique())
    df["ever_cvd"] = df["eid"].astype(str).isin(cvd_eids).astype("Int8")

def add_prior_cvd(df):
    cvd = pd.read_csv(CVD_TSV, sep="\t", dtype={"eid": str})
    cvd = cvd.rename(columns={"age": "cvd_age"})
    
    tmp = df.merge(cvd, on="eid", how="left")
    prior = (
        tmp["cvd_age"] < tmp["afs"]
    ).astype("Int8")
    df["prior_cvd"] = prior

def add_afs_als(df):
    phen = pd.read_csv(AFS_ALS_TSV, sep="\t", dtype={"eid": str})
    phen = phen.groupby("eid", as_index=False)["age"].agg(afs="min", als="max")
    return df.merge(phen, on="eid", how="left")

def filter_bp_plausibility(df):
    if "systolic" not in df or "diastolic" not in df:
        return
    keep = (
        (df["age"] > 0) &
        (((df["systolic"]  < 375) & (df["systolic"]  > 40)) | df["systolic"].isna()) &
        (((df["diastolic"] < 365) & (df["diastolic"] > 30)) | df["diastolic"].isna())
    )
    df.drop(index=df.index[~keep], inplace=True)

def filter_two_years_after_therapy(df, class_cols):
    if not class_cols:
        return
    d = df.sort_values(["eid","date"]).copy()
    d["drug"] = (d[class_cols].sum(axis=1) > 0)
    d["drug_cumsum"] = d.groupby("eid")["drug"].cumsum()
    d["days_to_previous"] = d.groupby("eid")["date"].diff()
    d.loc[d["drug"] == True, "days_to_previous"] = 0
    d["starts_with_0"] = (d.groupby(["eid","drug_cumsum"])["drug"].transform("sum") == 0)
    d["date_cumsum"] = d.groupby(["eid","drug_cumsum"])["days_to_previous"].cumsum()
    remove = ((d["date_cumsum"] < 730) & (~d["starts_with_0"]) & (d["drug"] == False))
    df.drop(index=d.index[remove], inplace=True)

def yob_map_df():
    y = pd.read_csv(YOB_CSV, dtype={"eid": str})
    y = y.rename(columns={"YOB": "yob"})[["eid", "yob"]]
    y["yob"] = pd.to_numeric(y["yob"])
    return y

def afs_als_df():
    phen = pd.read_csv(AFS_ALS_TSV, sep="\t", dtype={"eid": str})
    phen = (
        phen
        .groupby("eid", as_index=False)["age"]
        .agg(afs="min", als="max")
    )
    return phen[["eid", "afs", "als"]]

def attach_meta_yob_afs_als(df_wide, age_grouped=False): #this will add afs, als, yob to all 
    meta = yob_map_df().merge(afs_als_df(), on="eid", how="outer")
    if not age_grouped:
        add_ever_cvd(meta)
        add_prior_cvd(meta)
    out = df_wide.merge(meta, on="eid", how="left")
    return out

def aggregate_pre_post(df, class_cols, measures, pre_names_map, post_names_map):

    df = df.copy()
    any_on = df[class_cols].sum(axis=1) > 0

    pre = (
        df.loc[~any_on]
          .groupby("eid", as_index=False)[measures]
          .mean()
          .rename(columns=pre_names_map)
    )
    post = (
        df.loc[any_on]
          .groupby("eid", as_index=False)[measures]
          .mean()
          .rename(columns=post_names_map)
    )
    out = pre.merge(post, on="eid", how="outer")

    measno = df.groupby("eid").size().reset_index(name="measno")
    age_med = df.groupby("eid")["age"].median().reset_index(name="age")
    out = out.merge(measno, on="eid", how="left")
    out = out.merge(age_med, on="eid", how="left")

    bin_cols = list(class_cols) + ["GP", "ASCEN"]
    bin_agg = (
        df.groupby("eid")[bin_cols]
          .max()
          .astype("Int8")
          .reset_index()
    )
    out = out.merge(bin_agg, on="eid", how="left")
    for c in bin_cols:
        out[c] = out[c].astype("Int8")

    out = attach_meta_yob_afs_als(out)
    return out

def aggregate_pre_post_by_age(df, class_cols, age_map, measures, pre_names_map, post_names_map):

    df = df.copy()
    add_age_group_column(df, age_map)

    cvd = pd.read_csv(CVD_TSV, sep="\t", dtype={"eid": str})
    cvd = cvd.rename(columns={"age": "cvd_age"})
    df = df.merge(cvd, on="eid", how="left")

    any_on = df[class_cols].sum(axis=1) > 0

    group_lo = {g: age_map[g][0] for g in age_map}
    df["group_lo"] = df["age_group"].map(group_lo)
    df["prior_cvd"] = (
        (df["cvd_age"] < df["group_lo"]) &
        df["cvd_age"].notna() &
        df["group_lo"].notna()
    ).astype("Int8")

    pre = (
        df.loc[~any_on]
          .groupby(["eid", "age_group"], as_index=False)[measures]
          .mean()
          .rename(columns=pre_names_map)
    )
    post = (
        df.loc[any_on]
          .groupby(["eid", "age_group"], as_index=False)[measures]
          .mean()
          .rename(columns=post_names_map)
    )
    out = pre.merge(post, on=["eid", "age_group"], how="outer")

    measno = df.groupby(["eid", "age_group"]).size().reset_index(name="measno")
    age_med = df.groupby(["eid", "age_group"])["age"].median().reset_index(name="age")
    out = out.merge(measno, on=["eid", "age_group"], how="left")
    out = out.merge(age_med, on=["eid", "age_group"], how="left")

    bin_cols = list(class_cols) + ["GP", "ASCEN"]
    bin_agg = (
        df.groupby(["eid", "age_group"])[bin_cols]
          .max()
          .astype("Int8")
          .reset_index()
    )
    out = out.merge(bin_agg, on=["eid", "age_group"], how="left")
    for c in bin_cols:
        out[c] = out[c].astype("Int8")

    prior_agg = (
        df.groupby(["eid", "age_group"])["prior_cvd"]
          .max()
          .reset_index()
    )
    out = out.merge(prior_agg, on=["eid", "age_group"], how="left")

    out["age_group"] = out["age_group"].astype(int)
    out_wide = out.pivot(index="eid", columns="age_group")
    out_wide.columns = [f"{c}_{int(g)}" for (c, g) in out_wide.columns]
    out_wide = out_wide.reset_index()

    out_wide = attach_meta_yob_afs_als(out_wide, age_grouped = True) #this will add afs, als and yob only
    return out_wide

def reorder_for_bp(out, class_cols):
    ordered = (
        ["eid","SBP_pre","DBP_pre","SBP_post","DBP_post"] +
        [c for c in CLASS_ORDER if c in class_cols] +
        ["measno","ever_cvd","prior_cvd","afs","als","yob","age","GP","ASCEN"]

    )
    return out[[c for c in ordered if c in out.columns]]

def reorder_for_ldl(out, class_cols):
    middle = [c for c in ["age","GP","ASCEN"] if c in out.columns]
    ordered = (
        ["eid","LDL_pre"] + middle + ["LDL_post"] +
        [c for c in CLASS_ORDER if c in class_cols] +
        ["measno","ever_cvd","prior_cvd","afs","als","yob"]
    )
    return out[[c for c in ordered if c in out.columns]]

def copy_pre_into_post_for_never_on(agg_df, class_cols):
    df = agg_df.copy()
    if not class_cols:
        return df
    never_on = (df[class_cols].sum(axis=1) == 0)

    if "SBP_pre" in df and "SBP_post" in df:
        mask = never_on & df["SBP_post"].isna() & df["SBP_pre"].notna()
        df.loc[mask, "SBP_post"] = df.loc[mask, "SBP_pre"]
    if "DBP_pre" in df and "DBP_post" in df:
        mask = never_on & df["DBP_post"].isna() & df["DBP_pre"].notna()
        df.loc[mask, "DBP_post"] = df.loc[mask, "DBP_pre"]
    return df

def write_bp_aggregates(raw_tsv):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t")
    df_raw["eid"] = df_raw["eid"].astype(str)

    outputs = []
    for (lo, hi) in [(1,30),(1,60)]:
        df_bin, subs = binarize_by_window(df_raw, lo, hi, base_cols=BASE_BP_COLS)
        df_cls, cls_cols = map_to_classes(df_bin, subs)

        add_yob_and_filter_year(df_cls)
        filter_bp_plausibility(df_cls)
        filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

        measures = ["systolic","diastolic"]
        pre_map  = {"systolic":"SBP_pre","diastolic":"DBP_pre"}
        post_map = {"systolic":"SBP_post","diastolic":"DBP_post"}

        agg = aggregate_pre_post(df_cls, cls_cols, measures, pre_map, post_map)
        agg = reorder_for_bp(agg, cls_cols)

        path_all = os.path.join(outdir, f"bp_pre_post_{lo}to{hi}.tsv")
        agg.to_csv(path_all, sep="\t", index=False, na_rep="NA")
        outputs.append(path_all)

        agg_copy = copy_pre_into_post_for_never_on(agg, cls_cols)
        path_copy = os.path.join(outdir, f"bp_pre_post_{lo}to{hi}_with_copying.tsv")
        agg_copy.to_csv(path_copy, sep="\t", index=False, na_rep="NA")
        outputs.append(path_copy)

        if (lo, hi) == (1,60):
            on_mask = (agg[[c for c in cls_cols]].sum(axis=1) > 0)
            on_only = agg.loc[on_mask].copy()
            path_on = os.path.join(outdir, f"bp_pre_post_{lo}to{hi}_on_drugs.tsv")
            on_only.to_csv(path_on, sep="\t", index=False, na_rep="NA")
            outputs.append(path_on)

        both = agg["SBP_pre"].notna() & agg["SBP_post"].notna()
        only_pre  = agg["SBP_pre"].notna() & ~agg["SBP_post"].notna()
        only_post = ~agg["SBP_pre"].notna() & agg["SBP_post"].notna()
        print(f"[BP {lo}–{hi}] wrote: {path_all}  &  {path_copy}")
        print(f"  N={len(agg):,} | both={both.sum():,} | only pre={only_pre.sum():,} | only post={only_post.sum():,}")
        if (lo, hi) == (1,60):
            print(f"  on-drugs subset: {path_on} (N={len(on_only):,})")
    return outputs

def write_ldl_aggregates(raw_tsv):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t")
    df_raw["eid"] = df_raw["eid"].astype(str)

    outputs = []
    for (lo, hi) in [(1,30),(1,60)]:
        df_bin, subs = binarize_by_window(df_raw, lo, hi, base_cols=BASE_LDL_COLS)
        df_cls, cls_cols = map_to_classes(df_bin, subs)

        if "age" in df_cls.columns:
            add_yob_and_filter_year(df_cls)
        filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

        measures = ["ldl"]
        pre_map  = {"ldl":"LDL_pre"}
        post_map = {"ldl":"LDL_post"}

        agg = aggregate_pre_post(df_cls, cls_cols, measures, pre_map, post_map)
        agg = reorder_for_ldl(agg, cls_cols)

        path_all = os.path.join(outdir, f"ldl_pre_post_{lo}to{hi}.tsv")
        agg.to_csv(path_all, sep="\t", index=False, na_rep="NA")
        outputs.append(path_all)

        if (lo,hi) == (1,60):
            on_mask = (agg[[c for c in cls_cols]].sum(axis=1) > 0)
            on_only = agg.loc[on_mask].copy()
            path_on = os.path.join(outdir, f"ldl_pre_post_{lo}to{hi}_on_drugs.tsv")
            on_only.to_csv(path_on, sep="\t", index=False, na_rep="NA")
            outputs.append(path_on)

        both      = agg["LDL_pre"].notna() & agg["LDL_post"].notna()
        only_pre  = agg["LDL_pre"].notna() & ~agg["LDL_post"].notna()
        only_post = ~agg["LDL_pre"].notna() & agg["LDL_post"].notna()
        print(f"[LDL {lo}–{hi}] wrote: {path_all}")
        print(f"  N={len(agg):,} | both={both.sum():,} | only pre={only_pre.sum():,} | only post={only_post.sum():,}")
    return outputs

def write_bp_aggregates_nostatins_1to60(raw_tsv):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t").assign(eid=lambda d: d["eid"].astype(str))
    df_bin, subs = binarize_by_window(df_raw, 1, 60, base_cols=BASE_BP_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)

    cls_cols = [c for c in cls_cols if c.lower() != "statin"]
    if "statin" in df_cls.columns:
        df_cls.drop(columns=["statin"], inplace=True, errors="ignore")

    add_yob_and_filter_year(df_cls)
    filter_bp_plausibility(df_cls)
    filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

    measures = ["systolic","diastolic"]
    pre_map  = {"systolic":"SBP_pre","diastolic":"DBP_pre"}
    post_map = {"systolic":"SBP_post","diastolic":"DBP_post"}

    agg = aggregate_pre_post(df_cls, cls_cols, measures, pre_map, post_map)
    agg = reorder_for_bp(agg, cls_cols)

    path = os.path.join(outdir, "bp_pre_post_1to60_no_statins.tsv")
    agg.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(f"[BP 1–60 no_statins] wrote: {path}  | options: no_statins=True; age_split=None; copying=False")
    return path

def write_ldl_aggregates_statins_only_1to60(raw_tsv):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t").assign(eid=lambda d: d["eid"].astype(str))
    df_bin, subs = binarize_by_window(df_raw, 1, 60, base_cols=BASE_LDL_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)

    keep = [c for c in cls_cols if c.lower() == "statin"]
    drop = [c for c in cls_cols if c.lower() != "statin"]
    if drop:
        df_cls.drop(columns=drop, inplace=True, errors="ignore")
    cls_cols = keep

    add_yob_and_filter_year(df_cls)
    filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

    measures = ["ldl"]
    pre_map  = {"ldl":"LDL_pre"}
    post_map = {"ldl":"LDL_post"}

    agg = aggregate_pre_post(df_cls, cls_cols, measures, pre_map, post_map)
    agg = reorder_for_ldl(agg, cls_cols)

    path = os.path.join(outdir, "ldl_pre_post_1to60_statins_only.tsv")
    agg.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(f"[LDL 1–60 statins_only] wrote: {path}  | options: keep_only_statins=True; age_split=None; copying=False")
    return path

def add_age_group_column(df, age_map):
    ages = pd.to_numeric(df["age"], errors="coerce")
    groups = np.zeros(len(df), dtype=int)
    for g, (lo, hi) in age_map.items():
        groups[(ages >= lo) & (ages < hi)] = g
    df["age_group"] = groups

def write_bp_age_split_1to60(raw_tsv, age_map, include_statins=True, tag=None):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t").assign(eid=lambda d: d["eid"].astype(str))
    df_bin, subs = binarize_by_window(df_raw, 1, 60, base_cols=BASE_BP_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)

    if not include_statins:
        cls_cols = [c for c in cls_cols if c.lower() != "statin"]
        if "statin" in df_cls.columns:
            df_cls.drop(columns=["statin"], inplace=True, errors="ignore")

    add_yob_and_filter_year(df_cls)
    filter_bp_plausibility(df_cls)
    filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

    measures = ["systolic","diastolic"]
    pre_map  = {"systolic":"SBP_pre","diastolic":"DBP_pre"}
    post_map = {"systolic":"SBP_post","diastolic":"DBP_post"}

    agg = aggregate_pre_post_by_age(df_cls, cls_cols, age_map, measures, pre_map, post_map)

    stat_tag = "with_statins" if include_statins else "no_statins"
    tag = (tag or "age5")
    path = os.path.join(outdir, f"bp_pre_post_1to60_{tag}_{stat_tag}.tsv")
    agg.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(f"[BP 1–60 age_split] wrote: {path}  | options: age_split={tag}; {stat_tag}; copying=False")
    return path

def write_bp_age_split_1to60_four_groups(raw_tsv, age_map4, include_statins=True):
    return write_bp_age_split_1to60(raw_tsv, age_map=age_map4, include_statins=include_statins, tag="age4")

def write_bp_age_split_1to60_shifted_edges(raw_tsv, age_map_shifted, include_statins=True):
    return write_bp_age_split_1to60(raw_tsv, age_map=age_map_shifted, include_statins=include_statins, tag="age5_shifted")

def write_ldl_age_split_1to60(raw_tsv, age_map):
    outdir = os.path.dirname(raw_tsv)
    df_raw = pd.read_csv(raw_tsv, sep="\t").assign(eid=lambda d: d["eid"].astype(str))
    df_bin, subs = binarize_by_window(df_raw, 1, 60, base_cols=BASE_LDL_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)

    if "age" in df_cls.columns:
        add_yob_and_filter_year(df_cls)
    filter_two_years_after_therapy(df_cls, class_cols=cls_cols)

    measures = ["ldl"]
    pre_map  = {"ldl":"LDL_pre"}
    post_map = {"ldl":"LDL_post"}

    agg = aggregate_pre_post_by_age(df_cls, cls_cols, age_map, measures, pre_map, post_map)
    path = os.path.join(outdir, f"ldl_pre_post_1to60_age5.tsv")
    agg.to_csv(path, sep="\t", index=False, na_rep="NA")
    print(f"[LDL 1–60 age_split] wrote: {path}  | options: age_split=age5; copying=False")
    return path

def make_age_map_from_quantiles(df_with_age, k):
    a = pd.to_numeric(df_with_age["age"], errors="coerce").dropna()
    if a.empty:
        raise ValueError("No non-missing ages to compute quantiles.")
    qs = np.linspace(0, 1, k+1)
    edges = list(np.unique(np.quantile(a, qs)))
    edges[0] = -np.inf
    edges[-1] = np.inf
    age_map = {i+1: (edges[i], edges[i+1]) for i in range(len(edges)-1)}
    return age_map

def expand_middle_groups_age5(age_map_5, years=1.0):
    bounds = [age_map_5[1][0], age_map_5[1][1], age_map_5[2][1], age_map_5[3][1], age_map_5[4][1], age_map_5[5][1]]
    b0, b1, b2, b3, b4, b5 = bounds
    b1p = b1 - years
    b2p = b2 - years
    b3p = b3 + years
    b4p = b4 + years
    new = {
        1: (b0, b1p),
        2: (b1p, b2p),
        3: (b2p, b3p),
        4: (b3p, b4p),
        5: (b4p, b5)
    }
    print("new age map expanded:"+str(new))
    return new

def make_age_map(df):
    wdf = df.loc[:, ("eid", "age")]
    wdf.drop_duplicates(inplace=True)

def age_group_size(offsets, df):
    age_boundaries = [-np.inf] + list(np.cumsum(offsets)) + [np.inf]
    num_groups = len(offsets) + 1
    age_map = {i + 1: (age_boundaries[i], age_boundaries[i + 1]) for i in range(num_groups)}
    add_age_group_column(df, age_map)
    return df.groupby("age_group")["eid"].nunique()

def age_group_size_sd(offsets, df):
    return age_group_size(offsets, df).std()

def age_group_size_jac(offsets, df, delta=1):
    start_sd = age_group_size_sd(offsets, df)
    res = np.zeros_like(offsets)
    for i in range(len(offsets)):
        oc = offsets.copy()
        oc[i] += delta
        res[i] = (age_group_size_sd(oc, df) - start_sd) / delta
    return res

def gd(x0, grad_fn, step=0.0001, n_iter=1000, args=(), verbose=False, eval_fn=None):
    x = x0.copy().astype("float64")
    for iter_ix in range(n_iter):
        if verbose and eval_fn is not None:
            print(f"{iter_ix}: \t {x} \t {eval_fn(x, args)}")
        x -= grad_fn(x, args) * step
    return x

def gd_discrete(x0, grad_fn, step=1, n_iter=1000, args=(), verbose=False, eval_fn=None):
    x = x0.copy().astype("float64")
    for iter_ix in range(n_iter):
        if verbose and eval_fn is not None:
            print(f"{iter_ix}: \t {x} \t {eval_fn(x, args)}")
        x -= np.sign(grad_fn(x, args))
    return x

# filters for unrelated

def filter_aggregated_by_id(id_list_tsv, aggregated_tsv, out_suffix, eid_col = "eid"):
    df = pd.read_csv(aggregated_tsv, sep="\t")
    df[eid_col] = df[eid_col].astype(str)
    ids_df = pd.read_csv(id_list_tsv, sep="\t", dtype=str)
    eid_set = set(ids_df["0"].tolist())
    filtered = df[df[eid_col].isin(eid_set)].copy()
    base, ext = os.path.splitext(aggregated_tsv)
    out_tsv = f"{base}"+out_suffix
    filtered.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(f"filtered {len(filtered):,}/{len(df):,} rows using {len(eid_set):,} IDs → {out_tsv}")
    return out_tsv

# mixed models

def drop_constant(d, cols):
    return [c for c in cols if d[c].dropna().nunique() >= 2]

def build_model_frame_bp(df_raw, lo=1, hi=60):
    df_bin, subs = binarize_by_window(df_raw, lo, hi, base_cols=BASE_BP_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)
    add_yob_and_filter_year(df_cls)
    filter_bp_plausibility(df_cls)
    mdf = df_cls[["eid","systolic","diastolic","age","GP","ASCEN"] + cls_cols].copy()
    mdf["eid"] = mdf["eid"].astype(str)
    cls_cols = drop_constant(mdf, cls_cols)
    return mdf, cls_cols

def build_model_frame_ldl(df_raw, lo=1, hi=60):
    df_bin, subs = binarize_by_window(df_raw, lo, hi, base_cols=BASE_LDL_COLS)
    df_cls, cls_cols = map_to_classes(df_bin, subs)
    if "age" in df_cls.columns:
        add_yob_and_filter_year(df_cls)
    keep = ["eid","ldl","age","GP","ASCEN"] + cls_cols
    keep = [c for c in keep if c in df_cls.columns]
    mdf = df_cls[keep].copy()
    mdf["eid"] = mdf["eid"].astype(str)
    cls_cols = drop_constant(mdf, cls_cols)
    return mdf, cls_cols

def fit_mixed(long_df, outcome, cls_cols):
    rhs = [c for c in (cls_cols + ["GP","ASCEN","age"]) if c in long_df.columns]
    formula = f"{outcome} ~ " + " + ".join(rhs)

    d = long_df.dropna(subset=[outcome]).copy()
    if "age" in d.columns:
        d = d.dropna(subset=["age"])

    for c in cls_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0).astype("Int8")
    for c in ["GP","ASCEN"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0).astype("Int8")

    md = MixedLM.from_formula(formula, groups="eid", data=d, re_formula="1")
    m  = md.fit(method="lbfgs", reml=True, maxiter=300, disp=False)

    fe = m.fe_params
    se = m.bse_fe
    z  = fe / se
    p  = m.pvalues.loc[fe.index] if hasattr(m.pvalues, "loc") else pd.Series(m.pvalues, index=fe.index)
    ci = m.conf_int()
    ci = ci.loc[fe.index] if hasattr(ci, "loc") else pd.DataFrame(ci, index=fe.index)

    out = (pd.DataFrame({
            "term": fe.index,
            "beta": fe.values,
            "se":   se.values,
            "z":    z.values,
            "p_value": p.values,
            "ci_low":  ci.iloc[:,0].values,
            "ci_high": ci.iloc[:,1].values
         })
         .query("term != 'Intercept'")
         .reset_index(drop=True))

    class_rows = out[out["term"].isin(cls_cols)].copy()
    class_rows["outcome"] = outcome.upper()
    class_rows["formula"] = formula
    class_rows["n_rows"]  = int(m.nobs)
    class_rows["n_eid"]   = int(d["eid"].nunique())
    class_rows["converged"] = bool(m.converged)
    return class_rows, m

def run_mixed_bp_1to60(raw_tsv, save_prefix="mixed_bp_1to60"):
    df_raw = pd.read_csv(raw_tsv, sep="\t")
    df_raw["eid"] = df_raw["eid"].astype(str)
    mdf, cls_cols = build_model_frame_bp(df_raw, lo=1, hi=60)

    rows, models = [], {}
    for outcome in ["systolic","diastolic"]:
        coef_df, mdl = fit_mixed(mdf, outcome, cls_cols)
        rows.append(coef_df)
        models[outcome.upper()] = mdl
        print(f"[BP 1–60] {outcome.upper()}: rows={coef_df['n_rows'].iloc[0]}, eids={coef_df['n_eid'].iloc[0]}, converged={coef_df['converged'].iloc[0]}")

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(f"{save_prefix}_coef.tsv", sep="\t", index=False)
    return out, models

def run_mixed_ldl_1to60(raw_tsv, save_prefix="mixed_ldl_1to60"):
    df_raw = pd.read_csv(raw_tsv, sep="\t")
    df_raw["eid"] = df_raw["eid"].astype(str)
    mdf, cls_cols = build_model_frame_ldl(df_raw, lo=1, hi=60)

    coef_df, mdl = fit_mixed(mdf, outcome="ldl", cls_cols=cls_cols)
    coef_df.to_csv(f"{save_prefix}_coef.tsv", sep="\t", index=False)
    print(f"[LDL 1–60] rows={coef_df['n_rows'].iloc[0]}, eids={coef_df['n_eid'].iloc[0]}, converged={coef_df['converged'].iloc[0]}")
    return coef_df, mdl

# annotate ldl with age

def annotate_ldl_age(in_tsv, yob_csv, out_tsv=None, eid_col="eid", date_col="date",
                     yob_col_in_csv="YOB", out_col="age"):
    df_ldl = pd.read_csv(in_tsv, sep="\t")
    df_ldl[eid_col] = df_ldl[eid_col].astype(str)

    yob = pd.read_csv(yob_csv, dtype={eid_col: str})
    yob = yob[[eid_col, yob_col_in_csv]].copy()
    yob[yob_col_in_csv] = pd.to_numeric(yob[yob_col_in_csv], errors="coerce")

    df = df_ldl.merge(yob, on=eid_col, how="left")
    date_num = pd.to_numeric(df[date_col], errors="coerce")

    meas_year = date_num / 364.24
    age = meas_year - df[yob_col_in_csv]
    age = np.floor(age + 0.5)
    df[out_col] = age
    df.loc[(df[out_col] < 0) | (df[out_col] > 120), out_col] = np.nan
    df.drop(columns=[yob_col_in_csv], inplace=True)

    if out_tsv is None:
        out_tsv = in_tsv[:-4] + "_with_age.tsv" if in_tsv.lower().endswith(".tsv") else in_tsv + "_with_age.tsv"
    df.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(f"Wrote LDL with age → {out_tsv}  (n={len(df):,})")
    return out_tsv
