#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd

ALWAYS_ONE_TOKS = {"yob", "afs", "als"}

TOK_FIRST = {"pre", "age", "measno"}

re_group_suffix = re.compile(r'_(\d+)$')

def strip_group_suffix(t):
    return re_group_suffix.sub("", t)

def group_id(t):
    m = re.search(r'_(\d+)(?:_[A-Za-z]+)?$', t)
    return int(m.group(1)) if m else None

def has_token(name, tok):
    return re.search(rf'(?:^|_){re.escape(tok)}(?:_|$)', name) is not None

def is_post(name):
    b = strip_group_suffix(name).lower()
    return has_token(b, "post")

def is_first_block(name):
    b = strip_group_suffix(name).lower()
    if any(has_token(b, tok) for tok in TOK_FIRST):
        return True
    return strip_group_suffix(name).endswith("ADJ")

def is_always_one(name):
    b = strip_group_suffix(name).lower()
    return any(has_token(b, tok) for tok in ALWAYS_ONE_TOKS)

def is_prior_cvd(name):
    b = strip_group_suffix(name).lower()
    return has_token(b, "prior_cvd")

def is_ever_cvd(name):
    b = strip_group_suffix(name).lower()
    return has_token(b, "ever_cvd")

def read_traits_from_pxp(cor_path):
    cols = list(pd.read_csv(cor_path, sep="\t", nrows=0).columns)
    return [c.strip() for c in cols if c and not c.startswith("Unnamed:")]

def assign_time_indices(traits):
    has_groups = any(group_id(t) is not None for t in traits)
    idxs = []

    for t in traits:
        g = group_id(t)

        if is_prior_cvd(t):
            cat = 1
            
        elif is_ever_cvd(t):
            cat = 2
            
        elif is_always_one(t):
            idxs.append(1)
            continue
        
        elif is_post(t):
            cat = 3
            
        elif is_first_block(t):
            cat = 1
            
        else:
            cat = 2

        if has_groups and g is not None:
            idx = 3 * (g - 1) + cat
        else:
            idx = cat

        idxs.append(idx)

    return idxs

def find_pxp_dirs(root):
    root = os.path.abspath(root)
    dirs = []
    if not os.path.isdir(root):
        return dirs
    for name in sorted(os.listdir(root)):
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        cor_path = os.path.join(p, "cor_pxp.tsv")
        if os.path.exists(cor_path):
            dirs.append(p)
    return dirs

def main():
    if len(sys.argv) not in (2, 3):
        print(f"Usage: {sys.argv[0]} PXP_ROOT_DIR [OUTDIR_TIME_IND]", file=sys.stderr)
        sys.exit(1)

    pxp_root = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) == 3 else os.path.join(pxp_root, "time_indices")
    os.makedirs(outdir, exist_ok=True)

    pxp_dirs = find_pxp_dirs(pxp_root)
    if not pxp_dirs:
        print(f"[ERROR] No pxp subfolders with cor_pxp.tsv found under: {pxp_root}", file=sys.stderr)
        sys.exit(2)

    for pxp_dir in pxp_dirs:
        cor_path = os.path.join(pxp_dir, "cor_pxp.tsv")

        traits = read_traits_from_pxp(cor_path)
        if not traits:
            print(f"[WARN] No traits found in {cor_path} — skipping", file=sys.stderr)
            continue

        idxs = assign_time_indices(traits)

        base = os.path.basename(pxp_dir.rstrip("/"))
        out_path = os.path.join(outdir, f"{base}_time_index.txt")
        with open(out_path, "w") as fout:
            for ix in idxs:
                fout.write(f"{ix}\n")

        print(f"# {base}")
        for t, ix in zip(traits, idxs):
            print(f"{t}\t{ix}")

if __name__ == "__main__":
    main()
