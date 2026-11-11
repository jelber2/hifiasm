#!/usr/bin/env python3
"""Evaluate an assembly against the simulated VNTR reference.

Uses a PAF produced with `minimap2 -cx asm5 --cs=long` and walks the cs tag to
report, per VNTR locus, the assembled length, base identity and copy-number
error. Also reports whole-assembly summary statistics.

Usage:
    python3 eval.py --ref ref.fa --truth truth.json --paf asm.paf [--json out.json]
"""

import argparse
import json
import sys


def parse_cs(cs, tstart, tend, qstart, loci):
    """Walk a cs tag, attributing matches/errors to reference loci.

    Returns per-locus counters and the query coordinate at each locus boundary.
    """
    per = {
        i: {"match": 0, "mismatch": 0, "ins": 0, "del": 0,
            "q_at_start": None, "q_at_end": None}
        for i in range(len(loci))
    }
    t = tstart
    q = qstart
    i = 0
    n = len(cs)

    def note(t0, t1, kind, qlen=0):
        for li, loc in enumerate(loci):
            lo = max(t0, loc["start"])
            hi = min(t1, loc["end"])
            if hi > lo:
                per[li][kind] += hi - lo

    while i < n:
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            L = int(cs[i + 1 : j])
            note(t, t + L, "match")
            t += L
            q += L
            i = j
        elif c == "=":
            # --cs=long: '=' followed by the matching sequence
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            L = j - i - 1
            note(t, t + L, "match")
            t += L
            q += L
            i = j
        elif c == "*":
            b1, b2 = cs[i + 1], cs[i + 2]
            note(t, t + 1, "mismatch")
            t += 1
            q += 1
            i += 3
        elif c == "+":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            L = j - i - 1
            for li, loc in enumerate(loci):
                if loc["start"] <= t <= loc["end"]:
                    per[li]["ins"] += L
            q += L
            i = j
        elif c == "-":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            L = j - i - 1
            note(t, t + L, "del")
            t += L
            i = j
        elif c == "~":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            # long skip (splice); not expected in asm5
            t += 0
            i = j
        else:
            i += 1

    # query coordinate at each locus boundary: walk again, recording q(t)
    t = tstart
    q = qstart
    i = 0
    bounds = {li: {} for li in range(len(loci))}
    while i < n:
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            L = int(cs[i + 1 : j])
            for li, loc in enumerate(loci):
                if t <= loc["start"] <= t + L and "s" not in bounds[li]:
                    bounds[li]["s"] = q + (loc["start"] - t)
                if t <= loc["end"] <= t + L and "e" not in bounds[li]:
                    bounds[li]["e"] = q + (loc["end"] - t)
            t += L
            q += L
            i = j
        elif c == "=":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            L = j - i - 1
            for li, loc in enumerate(loci):
                if t <= loc["start"] <= t + L and "s" not in bounds[li]:
                    bounds[li]["s"] = q + (loc["start"] - t)
                if t <= loc["end"] <= t + L and "e" not in bounds[li]:
                    bounds[li]["e"] = q + (loc["end"] - t)
            t += L
            q += L
            i = j
        elif c == "*":
            for li, loc in enumerate(loci):
                if t <= loc["start"] <= t + 1 and "s" not in bounds[li]:
                    bounds[li]["s"] = q + (loc["start"] - t)
                if t <= loc["end"] <= t + 1 and "e" not in bounds[li]:
                    bounds[li]["e"] = q + (loc["end"] - t)
            t += 1
            q += 1
            i += 3
        elif c == "+":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            q += j - i - 1
            i = j
        elif c == "-":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            t += j - i - 1
            i = j
        elif c == "~":
            j = i + 1
            while j < n and cs[j] not in "=*+-~":
                j += 1
            i = j
        else:
            i += 1

    for li, loc in enumerate(loci):
        if "s" in bounds[li]:
            per[li]["q_at_start"] = bounds[li]["s"]
        elif loc["start"] <= tend:
            per[li]["q_at_start"] = qstart + max(0, loc["start"] - tstart)
        if "e" in bounds[li]:
            per[li]["q_at_end"] = bounds[li]["e"]
        elif loc["end"] <= tend:
            per[li]["q_at_end"] = qstart + max(0, loc["end"] - tstart)
    return per


def load_loci(truth_fn):
    with open(truth_fn) as fh:
        return json.load(fh)["loci"]


def load_fasta(fn):
    seqs = {}
    name = None
    parts = []
    for line in open(fn):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


_COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def flank_lengths(ref_seqs, asm_seqs, loci):
    """Repeat length between the unique flanks, measured by exact matching.

    Independent of alignment artifacts: the distance between the bases right
    before the array and right after it. Probe lengths are shortened until an
    exact hit is found, since the contig itself carries a few sequencing errors.
    """
    ref = list(ref_seqs.values())[0]
    out = []
    for loc in loci:
        best = None
        for probe in (40, 30, 25, 20, 16):
            lp = ref[max(0, loc["start"] - probe):loc["start"]]
            rp = ref[loc["end"]:loc["end"] + probe]
            if len(lp) < probe or len(rp) < probe:
                continue
            for name, s in asm_seqs.items():
                for strand_s in (s, revcomp(s)):
                    a = strand_s.find(lp)
                    b = strand_s.find(rp)
                    if a >= 0 and b >= 0 and b > a:
                        best = (b - (a + len(lp)), name)
                        break
                if best is not None:
                    break
            if best is not None:
                break
        out.append(best)
    return out


def contig_stats(fasta):
    lengths = []
    name = None
    length = 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths.append(length)
                name = line[1:].strip()
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths.append(length)
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    acc = 0
    n50 = 0
    for l in lengths:
        acc += l
        if acc >= half:
            n50 = l
            break
    return {"n_contigs": len(lengths), "total": total, "n50": n50,
            "longest": lengths[0] if lengths else 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--ref-fasta", default=None,
                    help="reference FASTA; enables flank-anchored repeat lengths")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    loci = load_loci(args.truth)
    aligned = {li: None for li in range(len(loci))}

    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            tstart, tend = int(f[7]), int(f[8])
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            if cs is None:
                continue
            per = parse_cs(cs, tstart, tend, int(f[2]), loci)
            for li, loc in enumerate(loci):
                ov = max(0, min(tend, loc["end"]) - max(tstart, loc["start"]))
                if ov <= 0:
                    continue
                score = per[li]["match"] + per[li]["mismatch"] + per[li]["del"]
                if aligned[li] is None or score > aligned[li]["score"]:
                    aligned[li] = {
                        "score": score,
                        "qname": f[0],
                        "per": per[li],
                        "ref_span": ov,
                    }

    flank_len = None
    if args.ref_fasta and args.fasta:
        ref_seqs = load_fasta(args.ref_fasta)
        asm_seqs = load_fasta(args.fasta)
        flank_len = flank_lengths(ref_seqs, asm_seqs, loci)

    rows = []
    for li, loc in enumerate(loci):
        a = aligned[li]
        if a is None:
            fl = flank_len[li] if flank_len is not None else None
            rows.append({
                "name": loc["name"],
                "unit_len": loc["unit_len"],
                "truth_len": loc["length"],
                "asm_len": fl[0] if fl else None,
                "len_delta": (fl[0] - loc["length"]) if fl else None,
                "status": "no_paf_alignment",
            })
            continue
        p = a["per"]
        denom = p["match"] + p["mismatch"] + p["ins"] + p["del"]
        ident = p["match"] / denom if denom else 0.0
        if p["q_at_start"] is not None and p["q_at_end"] is not None:
            asm_len = p["q_at_end"] - p["q_at_start"]
        else:
            asm_len = None
        if flank_len is not None:
            fl = flank_len[li]
            if fl is not None:
                asm_len = fl[0]
        rows.append({
            "name": loc["name"],
            "unit_len": loc["unit_len"],
            "truth_len": loc["length"],
            "truth_copies": loc["copies"],
            "asm_len": asm_len,
            "len_delta": (asm_len - loc["length"]) if asm_len is not None else None,
            "est_copies": (round(asm_len / loc["unit_len"], 1)
                           if asm_len is not None else None),
            "identity": round(ident, 5),
            "errors": denom - p["match"],
            "qname": a["qname"],
            "status": "mapped",
        })

    # whole-assembly identity over mapped bases
    tot_m = tot_e = 0
    with open(args.paf) as fh:
        for line in fh:
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
            per = parse_cs(cs, int(f[7]), int(f[8]), int(f[2]), [])
            # no loci -> per empty; count directly instead
            tot_m += int(f[9])
            tot_e += int(f[10]) - int(f[9])

    if args.fasta:
        stats = contig_stats(args.fasta)
        print(f"Contigs: {stats['n_contigs']}  total: {stats['total']}  "
              f"N50: {stats['n50']}  longest: {stats['longest']}")
    else:
        stats = {}
    print(f"{'locus':10s} {'unit':>4s} {'truth':>6s} {'asm':>6s} {'delta':>6s} "
          f"{'ident':>7s} {'errors':>6s}  status")
    for r in rows:
        if r["status"] != "mapped":
            al = r.get("asm_len")
            dl = r.get("len_delta")
            print(f"{r['name']:10s} {r['unit_len']:4d} {r['truth_len']:6d} "
                  f"{(al if al is not None else '-'):>6} {(dl if dl is not None else '-'):>6} "
                  f"{'-':>7s} {'-':>6s}  {r['status']}")
        else:
            print(f"{r['name']:10s} {r['unit_len']:4d} {r['truth_len']:6d} "
                  f"{r['asm_len']:6d} {r['len_delta']:+6d} {r['identity']*100:6.2f}% "
                  f"{r['errors']:6d}  {r['status']}")

    result = {"loci": rows, "aligned_bases": tot_m, "aligned_errors": tot_e,
              "aligned_identity": (tot_m / (tot_m + tot_e)) if tot_m + tot_e else 0.0,
              "assembly": stats}
    print(f"\nMapped bases: {tot_m}  errors: {tot_e}  identity: "
          f"{result['aligned_identity']*100:.3f}%")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
