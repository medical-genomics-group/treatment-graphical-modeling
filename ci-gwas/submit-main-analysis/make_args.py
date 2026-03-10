#!/usr/bin/env python3
import os
import sys

BLOCK_DIR = "path/to/blockfiles/"
OUT = "args-2.txt"

GEN_SETUPS = [
    "sbp_pre_post_1to60_rel_up_to_2nd_no_cvd",
]


def num_lines(path):
    with open(path, "r") as f:
        return sum(1 for _ in f)


def main():
    if not GEN_SETUPS:
        print("ERROR: GEN_SETUPS is empty", file=sys.stderr)
        return 1

    total = 0
    with open(OUT, "w") as fout:
        for cix in range(1, 23):  # chromosomes 1..22
            blockfile = os.path.join(BLOCK_DIR, f"c{cix}_m11000.blocks")
            if not os.path.exists(blockfile):
                print(f"WARNING: missing {blockfile}", file=sys.stderr)
                continue

            nblocks = num_lines(blockfile)
            for bix in range(1, nblocks + 1):
                for gs in GEN_SETUPS:
                    fout.write(f"{cix} {bix} {gs}\n")
                    total += 1

    print(f"wrote {total} lines to {OUT}")
    print(f"gen_setups: {len(GEN_SETUPS)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
