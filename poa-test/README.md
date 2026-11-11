# poa-consensus post-assembly VNTR polishing

Experimental integration of [poa-consensus](https://github.com/Psy-Fer/poa-consensus)
(v0.7.0, pinned as a git submodule) into hifiasm-fast, following the author's
suggestion:

> after the contig assembly step, cut out repeat regions with flanks, feed them
> to the single allele mode in poa-consensus with the reads that made them, you
> could get a better sequence, then insert it back in using the flanks as an
> anchor

## What was added

| File | Purpose |
|---|---|
| `poa-consensus/` | submodule, pinned to `f90b3dd` (v0.7.0) |
| `poa-ffi/` | tiny Rust crate exposing `consensus_flanked()` through a C ABI (`libpoa_ffi.a`), so the C++ assembler can call it in-process |
| `poa_polish.{h,cpp}` | the post-assembly pass: VNTR detection, read recovery, POA call, guarded splice-back, GFA/FASTA/TSV output |
| `main.cpp`, `CommandLines.*`, `Makefile` | wiring, CLI flags, build rules |
| `poa-test/` | ONT VNTR simulator, evaluation scripts, benchmark harness, data |

## How the pass works

1. **Walk every contig** in the GFA hifiasm just wrote
   (`*.p_ctg.gfa`, `*.bp.p_ctg.gfa`, alternate and trio hap files).
2. **Find tandem repeats** from the contig sequence: a single linear pass with a
   packed 11-mer table (`last occurrence` gives the unit length), then tolerant
   extension that allows a small phase drift, so a 1–2 bp indel inside an ONT
   array does not truncate the call. Periods 2–200 bp, arrays ≥ 60 bp and
   ≥ 2.5 copies.
3. **Recover the reads that made the contig — without loading the read set.**
   Every k-mer (17 bp) of the four flank edges (left-start, left-end,
   right-start, right-end; 60 bp each by default, see `--poa-edge-len`) is
   indexed per locus, in both strand orientations. The input is then streamed
   **once** per GFA and (by default) every 2nd k-mer of each read is looked up
   (`--poa-seed-step`). A read is recruited only if it hits all four sides in
   the right order, i.e. it reaches into the distal end of the left flank,
   through the repeat, and past the distal end of the right flank. Only these
   cut-down reads are retained (at most 128 per locus); everything else is
   discarded as it streams past. Read orientation is handled by indexing both
   the flank and its reverse complement. A read must also reach a seed near the
   distal end of each flank and be at least 70% as long as the contig array, so
   short reads that only dip into one flank cannot masquerade as spanning
   reads.
4. **Cut each read** to `[left flank | repeat | right flank]` and call
   `poa_consensus::consensus_flanked(reads, left, right, &PoaConfig::default())`
   — the author's single-allele flanked entry point, with homopolymer
   compression/gap costs and tandem-gap costs on (the crate defaults for ONT).
5. **Splice back between the flanks** only when the new repeat is credible:
   * the crate actually anchored both flanks in its consensus (checked against
     the flank cores, so an unsliced fallback is rejected),
   * the new length is within 0.4–2.5× the old length,
   * k-mer containment of the new repeat in the recruited reads is not worse
     than the old repeat's,
   * the contig's current length is a clear outlier against the median repeat
     length measured in the recruited reads (median absolute deviation), and
     the consensus moves it closer. This prevents chasing consensus noise on an
     already-correct array,
   * at least `--poa-min-reads` reads span the locus.

The original GFA is never modified. Outputs:

```
<prefix>.*.poa.gfa     patched contig graph (S sequences, A offsets shifted)
<prefix>.*.poa.fa      patched contigs as FASTA
<prefix>.poa.tsv       one row per locus: period, copies, lengths, reads,
                       k-mer support before/after, status
```

`A`-line offsets are shifted by the splice deltas so downstream read-to-contig
annotations stay consistent.

## Build

Requires `cargo`/`rustc` ≥ 1.85 (the crate is edition 2024) in addition to the
usual C++ toolchain and zlib. `make` builds the Rust static library
automatically and links it into `hifiasm`:

```bash
git submodule update --init --recursive
make -j
```

The standalone polish driver is built with:

```bash
make poa-polish-gfa
poa-test/polish_gfa --help
```

CI note: `.github/workflows/ci.yaml` must check out submodules
(`submodules: recursive` in `actions/checkout`) or the Rust crate will be
missing.

## Usage

```bash
hifiasm --ont --poa-polish -t32 -o sample sample.fq.gz
# -> sample.bp.p_ctg.poa.gfa / .poa.fa / sample.poa.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `--poa-polish` | off | enable the post-assembly pass |
| `--poa-flank INT` | 120 | flank length (bp) used as the anchor on each side |
| `--poa-min-reads INT` | 6 | minimum spanning reads to rewrite a locus |
| `--poa-max-period INT` | 200 | longest repeat unit considered |
| `--poa-max-candidate-bases INT` | 4,000,000,000 | cap on the cut-down reads kept in memory (repeat+flanks, counted once per locus they span); `0` = unlimited |
| `--poa-seed-step INT` | 2 | check every Nth read k-mer; 1 is most sensitive, 4–8 faster but may miss heavily errored flanks |
| `--poa-edge-len INT` | 60 | bp at each flank end used for seeding/anchoring (20–120); larger is more sensitive and uses more index memory |

For example, maximum sensitivity on a 256 GiB machine:

```bash
hifiasm --ont --poa-polish \
    --poa-max-candidate-bases 8000000000 \
    --poa-seed-step 1 \
    --poa-edge-len 120 \
    -t32 -o sample sample.fq.gz
```

Only `POA_POLISH_DEBUG_DIR` remains an environment variable (it dumps the cut
reads and consensus per locus).

Memory now scales with the number of repeat loci and the reads that span them,
not with the input read set. On the 30× benchmark (371 kb genome, 6 arrays),
peak RSS is ~109 MB; feeding the same reads 8× (0.09 Gbp) raises it only to
~139 MB because the per-locus read cap kicks in. The index costs roughly
350 seed k-mers per locus at the default 60 bp edge (≈15–18 KB/locus in the
current hash representation); `--poa-edge-len 35` roughly halves that, while
`--poa-edge-len 120` roughly doubles it. A human genome with ~1M detected arrays
therefore uses tens of GB, independent of coverage.

## Test data

`sim_ont_vntr.py` builds a single-haplotype reference with four to six tandem
arrays and emits ONT-like simplex reads (fastq):

* read lengths ~ Gamma(k=2, scale=N50/3), default N50 ≈ 10–12 kb, max 60 kb;
* per-read error rate ~ N(4.5%, 1.2%), clipped to [2%, 14%];
* substitutions/insertions/deletions ≈ 34/30/36% of the error budget,
  with deletions and insertions elevated 3× / 2× inside homopolymer runs
  (polymerase slippage);
* uniform coverage including the reference ends (reads may overhang).

Panels (`--panel`): unit lengths 3, 21, 29, 38, 51 bp, arrays 0.9–15.3 kb,
including homopolymer-rich units.

`eval.py` reports, per locus:

* **copy number** measured by exact flank anchoring (distance between the bases
  immediately before and after the array, strand-aware, probe lengths
  40→16 bp), which is independent of alignment artefacts;
* **identity/errors** inside the locus from `minimap2 -cx asm5 --cs=long`.

`degrade_contig.py` deletes a fraction of every array from a contig, producing a
controlled "assembler lost copies" input for the polish stage.

## Results (this repository, 8 threads)

Reproduce with `./run_benchmark.sh` (~1 min). Reads: 4.5% ONT-R10-like error;
truth arrays: A 21×42, B 38×40, C 51×160 (8.2 kb), D CAG×300, E 29×170 (4.9 kb),
F AAG×1200.

### 30× coverage (1 contig, 371,606 bp)

| locus | truth | original | +poa-polish |
|---|---:|---:|---:|
| VNTR_A | 882 | 882 (100%) | 882 |
| VNTR_B | 1520 | 1520 (100%) | 1520 |
| VNTR_C | 8160 | 8160 (100%) | 8160 |
| STR_D | 900 | 900 (100%) | 900 |
| VNTR_E | 4930 | 4930 (100%) | 4930 |
| VNTR_F | 3600 | 3600 (100%) | 3600 |

Mapped identity 99.997% in both. The pass found 6/6 loci (12, 20, 7, 25, 7, 10
spanning reads) and changed **nothing** (18 calls across primary/hap1/hap2).

### 20× coverage

| locus | truth | original | +poa-polish |
|---|---:|---:|---:|
| VNTR_A | 882 | 882 (100%) | 882 |
| VNTR_B | 1520 | 1520 (100%) | 1520 |
| VNTR_C | 8160 | **8165 (+5)** | 8165 (only 3 spanning reads, below floor) |
| STR_D | 900 | 900 (100%) | 900 |
| VNTR_E | 4930 | 4930 (100%) | 4930 |
| VNTR_F | 3600 | 3600 (100%) | 3600 |

Only the 8 kb array is slightly off, and only 3 reads span it at 20×, so the
pass correctly declines to touch it.

### 10× coverage (assembly fragments: 4 contigs, 192 kb)

| locus | truth | original | +poa-polish |
|---|---:|---:|---:|
| VNTR_A | 1260 | **1243 (−17)** | 1243 (4 reads, below floor) |
| VNTR_B | 1520 | **1414 (−106)** | 1414 (2 reads) |
| VNTR_C | 5100 | 5100 (100%) | 5100 (2 reads) |
| STR_D | 900 | **897 (−3)** | 897 (8 reads, within read noise) |

At 10× hifiasm does lose copies, but only 2–8 reads are recruited per array;
below the `--poa-min-reads` floor the pass makes no change. This is the fundamental constraint of the
author's flanked approach: it needs reads that span the array **and** both
flanks.

### Controlled repair (validates the mechanism)

Starting from the correct 30× contig, 8% of every array was deleted
(63, 114, 663, 72, 406, 288 bp) and only the polish stage was run with the same
reads:

| locus | truth | damaged | after polish |
|---|---:|---:|---:|
| VNTR_A | 882 | 819 (−63) | **882** |
| VNTR_B | 1520 | 1406 (−114) | **1520** |
| VNTR_C | 8160 | 7497 (−663) | **8155** |
| STR_D | 900 | 828 (−72) | **900** |
| VNTR_E | 4930 | 4524 (−406) | **4926** |
| VNTR_F | 3600 | 3312 (−288) | **3593** |

Locus identity 91.8–92.9% before, 99.75–100% after; all six arrays restored to
within 7 bp (three exactly, three within 4–7 bp).

## Interpretation

On these **idealised simulated** ONT data, hifiasm-fast's own `--ont` pipeline
already reconstructs VNTR copy numbers exactly once coverage is sufficient
(≥20×, and even a 15 kb array at 30×). The poa-consensus pass is therefore a
no-op at normal coverage, and it cannot rescue low-coverage assemblies because
the flanked single-allele method requires spanning reads. The author's predicted
benefit did not materialise here.

Where it can be expected to help in real data, and what would be worth testing
next:

* ONT runs with systematic homopolymer/stutter bias that survival of read
  correction does not remove (simulated per-base randomness is too kind to
  hifiasm);
* heterozygous VNTRs, where the pass currently keeps the allele nearest the
  median read length but could be given a haplotype-restricted read set from
  trio A-lines;
* arrays close to read length, where hifiasm's unitig graph is at its most
  uncertain but enough reads still span;
* arrays longer than the reads, which need `bridged_consensus` (reads covering
  only the two ends) rather than `consensus_flanked`.

Caveats: single-allele mode only; loci are processed sequentially; the TR
detector is heuristic (2–200 bp units); the seed index costs ~15–18 KB per
detected locus (so a highly repetitive genome with ~1M microsatellites will use
tens of GB); and the Rust static library adds a build-time dependency.

## Files

```
poa-ffi/                 Rust C-ABI shim (staticlib)
poa_polish.cpp/.h        hifiasm post-assembly pass
poa-test/sim_ont_vntr.py ONT VNTR read simulator
poa-test/eval.py         copy-number + identity evaluation
poa-test/degrade_contig.py  controlled copy-number damage
poa-test/polish_gfa.cpp  standalone polish driver
poa-test/run_benchmark.sh  end-to-end benchmark
poa-test/data/           30x / 20x / 10x datasets and outputs
```
