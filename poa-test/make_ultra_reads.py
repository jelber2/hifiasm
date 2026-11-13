#!/usr/bin/env python3
"""Build 1 Mbp ultra-long ONT-like test reads from a reference.

The reference is embedded in random sequence so the resulting read is much
longer than any repeat, then a forward and a reverse-complement copy are
written. Useful for checking that ultra-long reads are trimmed to the locus
instead of being passed to POA whole.

Usage:
  ./make_ultra_reads.py --ref ref.fa --outdir data/ont_vntr_30x \
      --length 1000000 --seed 42
"""

import argparse
import random


def read_fasta(path):
    parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="reference FASTA")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--length", type=int, default=1000000,
                    help="read length in bp [1000000]")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref = read_fasta(args.ref)
    if len(ref) >= args.length:
        raise SystemExit("reference is longer than the requested read length")

    pad = args.length - len(ref)
    seq = ("".join(rng.choice("ACGT") for _ in range(pad // 2)) + ref +
           "".join(rng.choice("ACGT") for _ in range(pad - pad // 2)))
    comp = str.maketrans("ACGT", "TGCA")
    rc = seq.translate(comp)[::-1]

    fwd = f"{args.outdir}/ultra_fwd.fa"
    rev = f"{args.outdir}/ultra_rev.fa"
    with open(fwd, "w") as fh:
        fh.write(f">ultra_fwd len={len(seq)}\n{seq}\n")
    with open(rev, "w") as fh:
        fh.write(f">ultra_rev len={len(rc)}\n{rc}\n")
    print(f"wrote {fwd} and {rev} ({len(seq)} bp each)")


if __name__ == "__main__":
    main()
