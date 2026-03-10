# postproc_utils.py

from dataclasses import dataclass
from pathlib import Path
import pandas as pd


AUTOSOMES = list(range(1, 23))


def read_pairs(pairs_tsv):
    df = pd.read_csv(pairs_tsv, sep="\t")
    return df


@dataclass(frozen=True)
class ChrCheck:
    run_id: str
    chr: int
    chr_dir: Path
    has_chr_dir: bool
    has_step2: bool
    has_clumps: bool


def _has_any(path, pattern):
    return next(path.glob(pattern), None) is not None


def check_run(
    step2_dir,
    run_id,
    chrs=AUTOSOMES,
    step2_glob="*step2*.regenie*",
    clumps_glob="*clump*.clumps*",
):
    step2_dir = Path(step2_dir)
    out = []

    run_dir = step2_dir / run_id
    for c in chrs:
        chr_dir = run_dir / f"chr{c}"
        has_chr_dir = chr_dir.is_dir()
        has_step2 = has_chr_dir and _has_any(chr_dir, step2_glob)
        has_clumps = has_chr_dir and _has_any(chr_dir, clumps_glob)

        out.append(
            ChrCheck(
                run_id=run_id,
                chr=c,
                chr_dir=chr_dir,
                has_chr_dir=has_chr_dir,
                has_step2=has_step2,
                has_clumps=has_clumps,
            )
        )
    return out


def check_all_runs(
    pairs_df,
    step2_dir,
    chrs=AUTOSOMES,
):
    rows = []
    for run_id in pairs_df["run_id"].astype(str).tolist():
        for rec in check_run(step2_dir=step2_dir, run_id=run_id, chrs=chrs):
            rows.append(
                {
                    "run_id": rec.run_id,
                    "chr": rec.chr,
                    "chr_dir": str(rec.chr_dir),
                    "has_chr_dir": rec.has_chr_dir,
                    "has_step2": rec.has_step2,
                    "has_clumps": rec.has_clumps,
                    "missing_step2": rec.has_chr_dir and (not rec.has_step2),
                    "missing_clumps": rec.has_chr_dir and (not rec.has_clumps),
                    "missing_chr_dir": not rec.has_chr_dir,
                }
            )
    return pd.DataFrame(rows)


def summarize_missing(check_df):
    """
    Returns:
      - per_run: one row per run_id with counts
      - missing_rows: only rows where something is missing
    """
    missing_rows = check_df[
        check_df["missing_chr_dir"] | check_df["missing_step2"] | check_df["missing_clumps"]
    ].copy()

    per_run = (
        check_df.groupby("run_id", as_index=False)
        .agg(
            n_chr=("chr", "count"),
            n_missing_chr_dir=("missing_chr_dir", "sum"),
            n_missing_step2=("missing_step2", "sum"),
            n_missing_clumps=("missing_clumps", "sum"),
        )
        .sort_values(
            ["n_missing_chr_dir", "n_missing_step2", "n_missing_clumps", "run_id"],
            ascending=[False, False, False, True],
        )
    )

    return per_run, missing_rows
