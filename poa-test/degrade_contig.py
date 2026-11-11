#!/usr/bin/env python3
"""Deliberately remove repeat copies from an assembled contig GFA.

This produces a controlled "assembler got the copy number wrong" input so the
poa-consensus polishing stage can be tested in isolation:

    ./polish_gfa degraded reads.fq.gz

Truth loci are mapped onto the contig with a minimap2 PAF (contig -> reference)
and `--delete-frac` of each array's length is removed, rounded to whole repeat
units (at least 2 units).

Usage:
    python3 degrade_contig.py --gfa asm.bp.p_ctg.gfa --paf asm.paf \
        --truth truth.json --out degraded.p_ctg.gfa --delete-frac 0.08
"""

import argparse
import json
import sys

from eval import parse_cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--delete-frac", type=float, default=0.08)
    args = ap.parse_args()

    loci = json.load(open(args.truth))["loci"]

    seqs = {}
    order = []
    for line in open(args.gfa):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            seqs[f[1]] = f[2]
            order.append(f[1])

    # contig coordinate range for each locus, via the best PAF alignment
    spans = {i: None for i in range(len(loci))}
    for line in open(args.paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        cs = None
        for tag in f[12:]:
            if tag.startswith("cs:Z:"):
                cs = tag[5:]
                break
        if cs is None:
            continue
        per = parse_cs(cs, int(f[7]), int(f[8]), int(f[2]), loci)
        for i, loc in enumerate(loci):
            p = per[i]
            if p["q_at_start"] is None or p["q_at_end"] is None:
                continue
            qlen_total = int(f[1])
            qs, qe = p["q_at_start"], p["q_at_end"]
            if f[4] == "-":
                # PAF/cs query coordinates run along the reverse-complemented
                # query; convert back to coordinates on the stored contig.
                qs, qe = qlen_total - qe, qlen_total - qs
            if qe <= qs:
                continue
            if spans[i] is None or (qe - qs) > (spans[i][1] - spans[i][0]):
                spans[i] = (qs, qe, f[0])

    edits = []
    for i, loc in enumerate(loci):
        if spans[i] is None:
            print(f"{loc['name']}: not on any contig; skipped")
            continue
        q0, q1, qname = spans[i]
        cur_len = q1 - q0
        unit = loc["unit_len"]
        n_del = max(2, round(cur_len * args.delete_frac / unit))
        del_bp = n_del * unit
        if del_bp <= 0 or cur_len - del_bp < unit:
            continue
        # remove from the middle of the array so both flanks stay intact
        at = q0 + (cur_len - del_bp) // 2
        seqs[qname] = seqs[qname][:at] + seqs[qname][at + del_bp:]
        edits.append((loc["name"], qname, q0, at, del_bp, n_del))
        print(f"{loc['name']}: {qname}:{q0}-{q1} delete {n_del} x {unit} bp "
              f"= {del_bp} bp at contig {at}")

    with open(args.out, "w") as fh:
        for n in order:
            s = seqs[n]
            fh.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
