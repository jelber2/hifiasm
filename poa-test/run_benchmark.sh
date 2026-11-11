#!/usr/bin/env bash
# Reproduce the ONT VNTR benchmark for hifiasm-fast +/- poa-consensus polishing.
#
# For each coverage setting:
#   1. simulate a reference with VNTR/STR arrays and ONT-like reads
#   2. assemble with plain hifiasm --ont
#   3. assemble with hifiasm --ont --poa-polish
#   4. map both assemblies to the truth and report per-locus copy-number error
#
# Finally, a controlled repair test: take the (correct) 30x assembly, remove 8%
# of every array to imitate an assembler that lost copies, and run only the
# polish stage (poa-test/polish_gfa).
#
# Usage: ./run_benchmark.sh [--force]
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
DATA="$HERE/data"
PY=${PYTHON:-python3}
THREADS=${THREADS:-8}
FORCE=${1:-}

if [[ ! -x "$ROOT/hifiasm" ]]; then
    echo "building hifiasm (and the poa-ffi static library)..." >&2
    make -C "$ROOT" -j"$THREADS" >/dev/null
fi
if [[ ! -x "$HERE/polish_gfa" ]]; then
    make -C "$ROOT" poa-polish-gfa >/dev/null
fi

extract_gfa_fa() { # $1 = gfa, $2 = fasta
    awk '/^S\t/{print ">"$2; print $3}' "$1" > "$2"
}

simulate() { # name coverage panel seed
    local name=$1 cov=$2 panel=$3 seed=$4
    local d="$DATA/$name"
    if [[ -f "$d/truth.json" && "$FORCE" != "--force" ]]; then
        echo "  (reusing $d)" >&2
        return
    fi
    "$PY" "$HERE/sim_ont_vntr.py" --outdir "$d" --coverage "$cov" \
        --err-rate 0.045 --hp-del 3 --hp-ins 2 --n50 12000 --seed "$seed" \
        --panel "$panel" >/dev/null
}

run_pair() { # name coverage panel seed
    local name=$1 cov=$2 panel=$3 seed=$4
    local d="$DATA/$name"
    simulate "$name" "$cov" "$panel" "$seed"
    echo
    echo "=================================================================="
    echo " dataset $name  (coverage ${cov}x, panel $panel, ONT R10-like)"
    echo "=================================================================="
    (
        cd "$d"
        echo "-- assembling (plain --ont)..." >&2
        "$ROOT/hifiasm" --ont -t "$THREADS" -o base reads.fq.gz > base.log 2>&1
        extract_gfa_fa base.bp.p_ctg.gfa base.ctg.fa
        minimap2 -cx asm5 --cs=long ref.fa base.ctg.fa > base.paf 2>/dev/null

        echo "-- assembling (--ont --poa-polish)..." >&2
        "$ROOT/hifiasm" --ont --poa-polish -t "$THREADS" -o poa reads.fq.gz > poa.log 2>&1
        extract_gfa_fa poa.bp.p_ctg.poa.gfa poa.ctg.fa
        minimap2 -cx asm5 --cs=long ref.fa poa.ctg.fa > poa.paf 2>/dev/null

        echo
        echo "--- ORIGINAL hifiasm-fast ---"
        "$PY" "$HERE/eval.py" --ref ref.fa --truth truth.json --paf base.paf \
            --fasta base.ctg.fa --ref-fasta ref.fa --json base.eval.json
        echo
        echo "--- WITH poa-consensus polishing ---"
        "$PY" "$HERE/eval.py" --ref ref.fa --truth truth.json --paf poa.paf \
            --fasta poa.ctg.fa --ref-fasta ref.fa --json poa.eval.json
        echo
        echo "polish log:"
        grep "poa_polish" poa.log | grep -v "hap1\|hap2" | sed 's/^/  /'
    )
}

run_controlled_repair() {
    local d="$DATA/ont_vntr_30x"
    echo
    echo "=================================================================="
    echo " controlled repair: 30x assembly with 8% of every array deleted"
    echo "=================================================================="
    (
        cd "$d"
        "$PY" "$HERE/degrade_contig.py" --gfa base.bp.p_ctg.gfa --paf base.paf \
            --truth truth.json --out degraded.p_ctg.gfa --delete-frac 0.08
        extract_gfa_fa degraded.p_ctg.gfa degraded.ctg.fa
        minimap2 -cx asm5 --cs=long ref.fa degraded.ctg.fa > degraded.paf 2>/dev/null
        echo
        echo "--- damaged contig ---"
        "$PY" "$HERE/eval.py" --ref ref.fa --truth truth.json --paf degraded.paf \
            --fasta degraded.ctg.fa --ref-fasta ref.fa
        echo
        echo "-- running only the polish stage..." >&2
        "$HERE/polish_gfa" degraded reads.fq.gz > polish.log 2>&1
        minimap2 -cx asm5 --cs=long ref.fa degraded.p_ctg.poa.fa > repaired.paf 2>/dev/null
        echo "--- after poa-consensus polishing ---"
        "$PY" "$HERE/eval.py" --ref ref.fa --truth truth.json --paf repaired.paf \
            --fasta degraded.p_ctg.poa.fa --ref-fasta ref.fa
        echo
        column -t degraded.poa.tsv | sed 's/^/  /'
    )
}

run_pair ont_vntr_30x 30 xl   5
run_pair ont_vntr_20x 20 xl   5
run_pair ont_vntr_10x 10 long 5
run_controlled_repair

echo
echo "done. raw data + GFA/PAF/JSON under poa-test/data/."
