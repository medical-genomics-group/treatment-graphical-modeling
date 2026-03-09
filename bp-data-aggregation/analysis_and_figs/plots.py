import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


COL = {
    "ascen": "#9DB4A5",
    "gp": "#5B7D8D",
    "mix": "#7FA39C",
    "other": "#3F6270",
    "a_window": "#8E5D9F",
    "gray": "0.25",
}

FIGSIZE = (5, 3)
DPI = 120
KDE_SMOOTH_MULT = 5.8
BASE_COLS = {"eid", "date", "age", "systolic", "diastolic", "GP", "ASCEN"}

def detect_drug_cols(df):
    return [c for c in df.columns if c not in BASE_COLS]


def compute_days_series(df):
    """All (measurement, drug) day offsets in [0, 90] across all drug columns."""
    drug_cols = detect_drug_cols(df)
    if not drug_cols:
        return pd.Series(dtype="int64")

    vals = df[drug_cols].apply(pd.to_numeric, errors="coerce")
    vals = vals.where(vals.ge(0) & vals.le(90))
    s = vals.stack().dropna().astype(int)
    return s


def add_flags_1_60(df):
    """Adds on_drug flags (1..60) + derived source/pattern flags."""
    out = df.copy()
    drug_cols = detect_drug_cols(out)

    if drug_cols:
        vals = out[drug_cols].apply(pd.to_numeric, errors="coerce")
        out["on_drug"] = (vals.ge(1) & vals.le(60)).any(axis=1)
    else:
        out["on_drug"] = False

    out["src"] = np.where(out["ASCEN"] == 1, "ascen", np.where(out["GP"] == 1, "gp", "other"))
    out["on_gp"] = out["on_drug"] & (out["GP"] == 1)

    out["off_any"] = ~out["on_drug"]
    out["off_ascen"] = out["off_any"] & (out["ASCEN"] == 1)
    out["off_ascen_only"] = out["off_ascen"] & (out["GP"] == 0)
    out["off_gp"] = out["off_any"] & (out["GP"] == 1)
    out["off_other"] = out["off_any"] & (out["ASCEN"] == 0) & (out["GP"] == 0)
    return out


def per_eid_summary(df):
    g = df.groupby("eid", sort=False)
    base = g.agg(
        n_meas=("eid", "size"),
        any_on=("on_drug", "any"),
        all_on=("on_drug", "all"),
        has_gp=("GP", "any"),
        has_ascen=("ASCEN", "any"),
        on_gp_cnt=("on_gp", "sum"),
        off_cnt=("off_any", "sum"),
        off_ascen_cnt=("off_ascen", "sum"),
        off_ascen_only_cnt=("off_ascen_only", "sum"),
        off_gp_cnt=("off_gp", "sum"),
        off_other_cnt=("off_other", "sum"),
    )
    base["any_off"] = ~base["all_on"]
    base["both_on_off"] = base["any_on"] & base["any_off"]
    base["both_sources"] = base["has_gp"] & base["has_ascen"]

    single_src = df.loc[g["eid"].transform("size") == 1].set_index("eid")["src"]
    return base.join(single_src.rename("single_src"), how="left")


def _nice_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _annotate_stacked(ax, x, bottoms, heights, labels, fontsize=9):
    for xi, btm, h, text in zip(x, bottoms, heights, labels):
        if h > 0:
            ax.text(xi, btm + h / 2, str(text), ha="center", va="center", fontsize=fontsize)


def _kde_curve(x, gridsize=512, smooth_mult=KDE_SMOOTH_MULT):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        xs = np.linspace((np.nanmin(x) - 1) if x.size else 0, (np.nanmax(x) + 1) if x.size else 1, 100)
        return xs, np.zeros_like(xs)

    xmin, xmax = np.percentile(x, [0.5, 99.5])
    rng = xmax - xmin if xmax > xmin else (x.std() or 1.0)
    xs = np.linspace(xmin - 0.05 * rng, xmax + 0.05 * rng, gridsize)

    try:
        from scipy.stats import gaussian_kde

        kde = gaussian_kde(x, bw_method="scott")
        try:
            kde.set_bandwidth(bw_method=kde.factor * smooth_mult)
        except Exception:
            kde = gaussian_kde(x, bw_method=smooth_mult)
        return xs, kde(xs)
    except Exception:
        hist, edges = np.histogram(x, bins=int(np.sqrt(max(10, x.size))) * 2, density=True)
        centers = (edges[:-1] + edges[1:]) / 2.0
        sigma_bins = 2.8 * smooth_mult
        radius = int(max(3, np.ceil(3 * sigma_bins)))
        t = np.arange(-radius, radius + 1)
        kern = np.exp(-(t**2) / (2 * sigma_bins**2))
        kern /= kern.sum()
        smooth = np.convolve(hist, kern, mode="same")
        ys = np.interp(xs, centers, smooth, left=0, right=0)
        return xs, ys

def _plot_density_two_with_medians(ax, arr_ascen, arr_gp, title_prefix, metric_label):
    xs1, ys1 = _kde_curve(arr_ascen)
    xs2, ys2 = _kde_curve(arr_gp)

    ax.plot(xs1, ys1, label="ascen", lw=2.0, color=COL["ascen"])
    ax.plot(xs2, ys2, label="gp", lw=2.0, color=COL["gp"])

    ax.set_ylim(bottom=0)
    if ax.get_ylim()[1] < 0.006:
        ax.set_ylim(top=0.006)

    xmin, xmax = ax.get_xlim()
    xrng = max(xmax - xmin, 1e-9)

    med_ascen = np.nanmedian(arr_ascen)
    med_gp = np.nanmedian(arr_gp)

    if np.isfinite(med_ascen):
        ax.axvline(med_ascen, color=COL["ascen"], lw=1.5, ls="--")
        ax.text(med_ascen - 0.02 * xrng, 0.006, f"median (ascen): {med_ascen:.3f}",
                color=COL["ascen"], ha="right", va="bottom", fontsize=9)
    if np.isfinite(med_gp):
        ax.axvline(med_gp, color=COL["gp"], lw=1.5, ls="--")
        ax.text(med_gp - 0.02 * xrng, 0.002, f"median (gp): {med_gp:.3f}",
                color=COL["gp"], ha="right", va="bottom", fontsize=9)

    ax.set_title(f"{title_prefix}. {metric_label} density by source")
    ax.set_xlabel(metric_label)
    ax.set_ylabel("density")
    ax.legend(frameon=False)
    _nice_axes(ax)

def plot_A_days_since_hist(df, save_path="plot_a_days_since.svg", x_mark=48.71):
    s = compute_days_series(df)
    if s.empty:
        print("no (measurement, drug) pairs in [0, 90] days.")
        return

    counts = s.value_counts().reindex(range(0, 91), fill_value=0).sort_index()
    frac = counts / counts.sum()

    colors = [
        COL["a_window"] if (1 <= d <= 60) else COL["mix"] 
        for d in counts.index
    ]

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    ax.bar(counts.index.values, frac.values, color=colors, edgecolor="white", width=0.98)
    ax.set_xlim(-0.5, 90.5)
    ax.set_xlabel("days since prescription")
    ax.set_ylabel("fraction of measurement–drug pairs")
    ax.set_title("A. BP measurement to prescription distance")
    ax.grid(axis="y", alpha=0.2, linestyle=":", linewidth=0.7)
    _nice_axes(ax)

    from matplotlib.patches import Patch
    ax.legend(
        handles=[
            Patch(facecolor=COL["a_window"], edgecolor="white", label="days 1–60"),
            Patch(facecolor=COL["mix"], edgecolor="white", label="day 0 & 61–90"),
        ],
        frameon=False,
    )

    ax.axvline(x_mark, color=COL["gray"], lw=1.8)
    y_top = ax.get_ylim()[1] * 0.95
    ax.text(x_mark + 1.0, y_top, f"{x_mark:.2f}", color=COL["gray"], ha="left", va="top", fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, format="svg")
    plt.show()

    drug_cols = detect_drug_cols(df)
    if drug_cols:
        vals = df[drug_cols].apply(pd.to_numeric, errors="coerce")
        vals = vals.where(vals.ge(0) & vals.le(90))

        row_has_zero = vals.eq(0).any(axis=1)
        row_has_1_60 = (vals.ge(1) & vals.le(60)).any(axis=1)

        denom = int(row_has_zero.sum())
        numer = int((row_has_zero & row_has_1_60).sum())
        pct = (100.0 * numer / denom) if denom else float("nan")

        print(f"[a] day-0 measurement rows with any 1–60d prescription on the same eid/date: "
              f"{numer}/{denom} ({pct:.1f}%)")
    else:
        print("[a] no drug columns detected — cannot compute day-0 overlap metric.")


def plot_B_hist_n_measurements(df1_60, save_path="plot_b_n_measurements.svg"):
    summ = per_eid_summary(df1_60)
    n = summ["n_meas"].to_numpy()

    bin_vals = np.arange(1, 51)
    counts_1_50 = np.array([(n == k).sum() for k in bin_vals])
    count_over = int((n > 50).sum())

    x = np.arange(1, 52) 
    heights = np.concatenate([counts_1_50, [count_over]])

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    ax.bar(x, heights, width=0.9, color=COL["ascen"], edgecolor="white")
    ax.set_xlim(0, 52)
    ax.set_ylabel("no of participants")
    ax.set_xlabel("no of measurements")
    ax.set_title("B. BP measurements per participant")

    tick_pos = [1, 5, 10, 20, 30, 40, 51]
    tick_lbl = ["1", "5", "10", "20", "30", "40", ">50"]
    ax.set_xticks(tick_pos, tick_lbl)

    med = float(np.mean(n)) if n.size else float("nan")
    if np.isfinite(med):
        xline = 51 if med > 50 else med
        ax.axvline(xline, color=COL["gray"], lw=1.8, ls="--")
        y_top = ax.get_ylim()[1] * 0.95
        ax.text(xline + 0.5, y_top, f"median={med:.1f}", color=COL["gray"], ha="left", va="top", fontsize=9)

    _nice_axes(ax)
    plt.tight_layout()
    plt.savefig(save_path, format="svg")
    plt.show()

    print(f"[b] individuals: {len(n)} | mean={n.mean():.2f} | median={np.median(n):.0f} | max={n.max()} | >50: {count_over}")


def plot_C_two_bars_by_source(df1_60, save_path="plot_c_by_source.svg"):
    summ = per_eid_summary(df1_60)

    is_single = (summ["n_meas"] == 1)
    s_src = summ.loc[is_single, "single_src"].fillna("other")
    s_counts = {"gp": int((s_src == "gp").sum()), "ascen": int((s_src == "ascen").sum()), "mix": 0, "other": int((s_src == "other").sum())}

    is_multi = (summ["n_meas"] >= 2)
    only_gp = is_multi & summ["has_gp"] & ~summ["has_ascen"]
    only_ascen = is_multi & ~summ["has_gp"] & summ["has_ascen"]
    mix = is_multi & summ["has_gp"] & summ["has_ascen"]
    none_src = is_multi & ~summ["has_gp"] & ~summ["has_ascen"]

    m_counts = {"gp": int(only_gp.sum()), "ascen": int(only_ascen.sum()), "mix": int(mix.sum()), "other": int(none_src.sum())}

    cats = ["single measurement", "two or more measurements"]
    x = np.arange(len(cats))
    width = 0.65

    seg_order = ["gp", "ascen", "mix", "other"]
    seg_colors = {"gp": COL["gp"], "ascen": COL["ascen"], "mix": COL["mix"], "other": COL["other"]}

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    bottoms = np.zeros_like(x, dtype=float)

    for key in seg_order:
        heights = np.array([s_counts[key], m_counts[key]], dtype=float)
        if key == "other" and heights.sum() == 0:
            continue
        ax.bar(x, heights, width, bottom=bottoms, label=key, color=seg_colors[key], edgecolor="white")
        _annotate_stacked(ax, x, bottoms, heights, [int(h) for h in heights])
        bottoms += heights

    ax.set_xticks(x, cats)
    ax.set_ylabel("no of participants")
    ax.set_title("C. Participants with 2 or more BP measurements stratified by source")
    ax.legend(frameon=False, ncols=4)
    _nice_axes(ax)

    plt.tight_layout()
    plt.savefig(save_path, format="svg")
    plt.show()


def plot_D_patterns_for_two_plus(df1_60, save_path="plot_d_patterns.svg"):
    summ = per_eid_summary(df1_60)

    is_two = (summ["n_meas"] == 2)
    is_3plus = (summ["n_meas"] >= 3)

    cond_two_target = is_two & (summ["on_gp_cnt"] >= 1) & (summ["off_ascen_cnt"] >= 1)
    cond_two_other = is_two & ~cond_two_target

    cond_3plus_target = is_3plus & (summ["off_cnt"] >= 1) & (summ["off_cnt"] == summ["off_ascen_only_cnt"])
    cond_3plus_other = is_3plus & ~cond_3plus_target

    counts = [int(cond_two_target.sum()), int(cond_two_other.sum()), int(cond_3plus_target.sum()), int(cond_3plus_other.sum())]
    labels = [
        "2 measures: gp on-drug and assesment centre off-drug",
        "2 measures: other pattern",
        "3+ measures: the only off-drug is assesment centre",
        "3+ measures: other pattern",
    ]
    colors = [COL["ascen"], COL["gp"], COL["mix"], COL["other"]]

    x = np.arange(4)
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    ax.bar(x, counts, color=colors, edgecolor="white")
    for xi, h in zip(x, counts):
        if h > 0:
            ax.text(xi, h + max(1, 0.01 * h), str(h), ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x, labels, rotation=15, ha="right")
    ax.set_ylabel("no of participants")
    ax.set_title("D. Measurement source patterns")
    _nice_axes(ax)

    plt.tight_layout()
    plt.savefig(save_path, format="svg")
    plt.show()

    print("[d] counts (in order):", dict(zip(labels, counts)))


def _plot_density(df1_60, col, title_prefix, metric_label, save_path):
    a = df1_60.loc[df1_60["ASCEN"] == 1, col].dropna().to_numpy()
    g = df1_60.loc[df1_60["ASCEN"] == 0, col].dropna().to_numpy()

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    _plot_density_two_with_medians(ax, a, g, title_prefix, metric_label)
    plt.tight_layout()
    plt.savefig(save_path, format="svg")
    plt.show()

    print(f"[{title_prefix}] n (ascen) {metric_label} = {len(a)} | n (gp) = {len(g)}")


def plot_E_density_sbp(df1_60, save_path="plot_e_density_sbp.svg"):
    _plot_density(df1_60, "systolic", "e", "SBP", save_path)


def plot_F_density_dbp(df1_60, save_path="plot_f_density_dbp.svg"):
    _plot_density(df1_60, "diastolic", "f", "DBP", save_path)

def make_all_plots(df, outdir="."):
    os.makedirs(outdir, exist_ok=True)

    plot_A_days_since_hist(df, save_path=os.path.join(outdir, "plot_a_days_since.svg"))
    df1_60 = add_flags_1_60(df)

    plot_B_hist_n_measurements(df1_60, save_path=os.path.join(outdir, "plot_b_n_measurements.svg"))
    plot_C_two_bars_by_source(df1_60, save_path=os.path.join(outdir, "plot_c_by_source.svg"))
    plot_D_patterns_for_two_plus(df1_60, save_path=os.path.join(outdir, "plot_d_patterns.svg"))
    plot_E_density_sbp(df1_60, save_path=os.path.join(outdir, "plot_e_density_sbp.svg"))
    plot_F_density_dbp(df1_60, save_path=os.path.join(outdir, "plot_f_density_dbp.svg"))
