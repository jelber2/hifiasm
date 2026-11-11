#!/usr/bin/env python3
"""Simulate an ONT-like read set over a reference carrying several VNTR/STR loci.

The reference is a single ~250 kb haplotype with four tandem-repeat loci that
differ in unit length (6, 21, 37, 61 bp) and copy number. Reads are emitted with
a length profile and error profile loosely matching ONT R9.4/R10 simplex data:

  * per-read error rate ~ N(5%, 1%) clipped to [2%, 10%]
  * substitutions ~3.5% per base
  * deletions ~3.5% per base, elevated in homopolymer runs (polymerase slippage)
  * insertions ~3.5% per base, also homopolymer-biased (usually a duplication of
    the current base, which is what ONT tends to produce)

Usage:
    python3 sim_ont_vntr.py --outdir data/ont_vntr --coverage 45 --seed 7
"""

import argparse
import gzip
import json
import os
import random


COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.41):
    return "".join(
        rng.choices("ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])[0]
        for _ in range(n)
    )


def make_imperfect_copies(rng, unit, copies, sub_rate=0.02):
    parts = []
    for _ in range(copies):
        u = list(unit)
        for i in range(len(u)):
            if rng.random() < sub_rate:
                u[i] = rng.choice([b for b in "ACGT" if b != u[i]])
        parts.append("".join(u))
    return "".join(parts)


def build_reference(rng, loci_spec):
    """Interleave unique sequence with tandem arrays at the requested positions."""
    pieces = []
    cursor = 0
    truth = []
    for spec in loci_spec:
        gap = spec["pos"] - cursor
        pieces.append(rand_seq(rng, gap))
        unit = spec["unit"]
        repeats = make_imperfect_copies(rng, unit, spec["copies"])
        start = sum(len(p) for p in pieces)
        pieces.append(repeats)
        end = start + len(repeats)
        truth.append(
            {
                "name": spec["name"],
                "start": start,
                "end": end,
                "unit": unit,
                "unit_len": len(unit),
                "copies": spec["copies"],
                "length": len(repeats),
            }
        )
        cursor = spec["pos"] + len(repeats)
    pieces.append(rand_seq(rng, 8000))
    return "".join(pieces), truth


def simulate_reads(rng, ref, coverage, n50, max_len, err_rate, hp_del, hp_ins):
    """Emit (name, sequence) pairs with ONT-ish length and error profiles."""
    target_bases = int(len(ref) * coverage)
    total = 0
    reads = []
    i = 0
    # gamma shape/scale tuned so the length N50 lands near `n50`
    shape, scale = 2.0, n50 / 3.0
    while total < target_bases:
        length = int(rng.gammavariate(shape, scale))
        length = max(1200, min(max_len, length))
        if length >= len(ref) - 1000:
            length = len(ref) - 1000
        # allow overhang so the reference ends are covered like the middle
        start = rng.randrange(-length + 1, len(ref))
        lo = max(0, start)
        hi = min(len(ref), start + length)
        if hi - lo < 600:
            continue
        sub = ref[lo:hi]
        if rng.random() < 0.5:
            sub = revcomp(sub)
        seq = introduce_errors(rng, sub, err_rate, hp_del, hp_ins)
        name = f"ont_read_{i:06d}"
        reads.append((name, seq))
        total += len(seq)
        i += 1
    return reads


def introduce_errors(rng, s, base_err=0.09, hp_del=2.4, hp_ins=2.0):
    """Apply ONT-like substitutions/indels, biased inside homopolymer runs."""
    err_rate = max(0.02, min(0.14, rng.gauss(base_err, 0.012)))
    # ONT R9.4/R10 simplex: roughly equal substitutions/insertions/deletions,
    # with deletions dominating in homopolymers (slippage).
    p_sub = err_rate * 0.34
    p_del = err_rate * 0.36
    p_ins = err_rate * 0.30
    out = []
    n = len(s)
    i = 0
    while i < n:
        # homopolymer context: how far the current run extends
        j = i + 1
        while j < n and s[j] == s[i]:
            j += 1
        run = j - i
        p_del_here = p_del * (hp_del if run >= 4 else 1.0)
        if rng.random() < p_del_here:
            # drop this base (slippage); in long runs occasionally drop one more
            i += 1
            if run >= 5 and rng.random() < 0.25:
                i += 1
            continue

        base = s[i]
        if rng.random() < p_sub:
            base = rng.choice([b for b in "ACGT" if b != s[i]])
        out.append(base)

        p_ins_here = p_ins * (hp_ins if run >= 4 else 1.0)
        if rng.random() < p_ins_here:
            out.append(rng.choice([s[i], base, rng.choice("ACGT")]))
        i += 1
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/ont_vntr")
    ap.add_argument("--coverage", type=float, default=45.0)
    ap.add_argument("--n50", type=int, default=12000)
    ap.add_argument("--err-rate", type=float, default=0.045,
                    help="mean per-read ONT error rate (default 0.045 = ONT R10)")
    ap.add_argument("--hp-del", type=float, default=2.4,
                    help="deletion multiplier inside homopolymer runs")
    ap.add_argument("--hp-ins", type=float, default=2.0,
                    help="insertion multiplier inside homopolymer runs")
    ap.add_argument("--panel", choices=["short", "long", "hp", "xl", "xlong"], default="short")
    ap.add_argument("--max-len", type=int, default=60000)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # VNTR units are arbitrary but must be non-periodic, so the locus is a true
    # VNTR rather than a longer STR.
    loci_spec = [
        {"name": "VNTR_A", "pos": 30000, "unit": "GGTGTGGGACAGGTGTGGGAC", "copies": 42},
        {"name": "VNTR_B", "pos": 95000,
         "unit": "CCTAGGCAGTCCATGGACTTCAGGTCAGTCCATGGACT", "copies": 27},
        {"name": "VNTR_C", "pos": 165000,
         "unit": "TTCCAGGCTAGGTACCGATGCTAGGCATCCGATTCGAGGCCTAGGTACCGA",
         "copies": 60},
        {"name": "STR_D", "pos": 230000, "unit": "CAG", "copies": 300},
        {"name": "VNTR_E", "pos": 290000,
         "unit": "AAAAGAAGAGGGAAAGGGAGAAAAGAGGG", "copies": 45},
    ]
    if args.panel == "long":
        loci_spec = [
            {"name": "VNTR_A", "pos": 30000, "unit": "GGTGTGGGACAGGTGTGGGAC", "copies": 60},
            {"name": "VNTR_B", "pos": 95000,
             "unit": "CCTAGGCAGTCCATGGACTTCAGGTCAGTCCATGGACT", "copies": 40},
            {"name": "VNTR_C", "pos": 165000,
             "unit": "TTCCAGGCTAGGTACCGATGCTAGGCATCCGATTCGAGGCCTAGGTACCGA",
             "copies": 100},
            {"name": "STR_D", "pos": 230000, "unit": "CAG", "copies": 300},
            {"name": "VNTR_E", "pos": 290000,
             "unit": "AAAAGAAGAGGGAAAGGGAGAAAAGAGGG", "copies": 80},
        ]
    elif args.panel == "xl":
        loci_spec = [
            {"name": "VNTR_A", "pos": 30000, "unit": "GGTGTGGGACAGGTGTGGGAC", "copies": 42},
            {"name": "VNTR_B", "pos": 95000,
             "unit": "CCTAGGCAGTCCATGGACTTCAGGTCAGTCCATGGACT", "copies": 40},
            {"name": "VNTR_C", "pos": 165000,
             "unit": "TTCCAGGCTAGGTACCGATGCTAGGCATCCGATTCGAGGCCTAGGTACCGA",
             "copies": 160},
            {"name": "STR_D", "pos": 240000, "unit": "CAG", "copies": 300},
            {"name": "VNTR_E", "pos": 300000,
             "unit": "AAAAGAAGAGGGAAAGGGAGAAAAGAGGG", "copies": 170},
            {"name": "VNTR_F", "pos": 360000, "unit": "AAG", "copies": 1200},
        ]
    elif args.panel == "xlong":
        loci_spec = [
            {"name": "VNTR_A", "pos": 30000, "unit": "GGTGTGGGACAGGTGTGGGAC", "copies": 42},
            {"name": "VNTR_B", "pos": 95000,
             "unit": "CCTAGGCAGTCCATGGACTTCAGGTCAGTCCATGGACT", "copies": 40},
            {"name": "VNTR_C", "pos": 165000,
             "unit": "TTCCAGGCTAGGTACCGATGCTAGGCATCCGATTCGAGGCCTAGGTACCGA",
             "copies": 300},
            {"name": "STR_D", "pos": 250000, "unit": "CAG", "copies": 300},
            {"name": "VNTR_E", "pos": 310000,
             "unit": "AAAAGAAGAGGGAAAGGGAGAAAAGAGGG", "copies": 250},
            {"name": "VNTR_F", "pos": 380000, "unit": "AAG", "copies": 1200},
        ]
    elif args.panel == "hp":
        loci_spec = [
            {"name": "VNTR_A", "pos": 30000, "unit": "GGTGTGGGACAGGTGTGGGAC", "copies": 42},
            {"name": "VNTR_B", "pos": 95000,
             "unit": "CCTAGGCAGTCCATGGACTTCAGGTCAGTCCATGGACT", "copies": 40},
            {"name": "VNTR_C", "pos": 165000,
             "unit": "TTCCAGGCTAGGTACCGATGCTAGGCATCCGATTCGAGGCCTAGGTACCGA",
             "copies": 80},
            {"name": "STR_D", "pos": 230000, "unit": "CAG", "copies": 300},
            {"name": "VNTR_E", "pos": 290000,
             "unit": "AAAAAAGAGGGAAAAGGGGAAAAAGAGGGGAAAGGG", "copies": 55},
            {"name": "VNTR_F", "pos": 350000, "unit": "AAG", "copies": 1200},
        ]

    ref, truth = build_reference(rng, loci_spec)
    ref_name = "sim_haplotype"
    with open(os.path.join(args.outdir, "ref.fa"), "w") as fh:
        fh.write(f">{ref_name}\n")
        for i in range(0, len(ref), 80):
            fh.write(ref[i : i + 80] + "\n")

    with open(os.path.join(args.outdir, "vntr_truth.bed"), "w") as fh:
        for t in truth:
            fh.write(
                f"{ref_name}\t{t['start']}\t{t['end']}\t{t['name']}\t"
                f"{t['unit_len']}\t{t['copies']}\n"
            )

    reads = simulate_reads(rng, ref, args.coverage, args.n50, args.max_len,
                           args.err_rate, args.hp_del, args.hp_ins)
    qtab = [chr(33 + q) for q in range(0, 42)]
    with gzip.open(os.path.join(args.outdir, "reads.fq.gz"), "wt") as fh:
        for name, seq in reads:
            # ONT-ish per-base quality in Q8..Q25
            quals = "".join(rng.choice(qtab[8:26]) for _ in range(len(seq)))
            fh.write(f"@{name}\n{seq}\n+\n{quals}\n")

    stats = {
        "reference": ref_name,
        "reference_length": len(ref),
        "coverage": sum(len(s) for _, s in reads) / len(ref),
        "n_reads": len(reads),
        "n50": n50_of([len(s) for _, s in reads]),
        "mean_read_len": sum(len(s) for _, s in reads) / len(reads),
        "loci": truth,
        "seed": args.seed,
    }
    with open(os.path.join(args.outdir, "truth.json"), "w") as fh:
        json.dump(stats, fh, indent=2)

    print(json.dumps({k: v for k, v in stats.items() if k != "loci"}, indent=2))
    for t in truth:
        print(
            f"  {t['name']}: unit={t['unit_len']}bp copies={t['copies']} "
            f"len={t['length']}bp at {t['start']}-{t['end']}"
        )


def n50_of(lengths):
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    acc = 0
    for l in lengths:
        acc += l
        if acc >= half:
            return l
    return 0


if __name__ == "__main__":
    main()
