/**
 * poa_polish.cpp -- post-assembly VNTR/tandem-repeat polishing for hifiasm-fast.
 *
 * Implements the workflow suggested by the poa-consensus author:
 *
 *   1. After the contig graph is written, walk every contig in the GFA.
 *   2. Find tandem-repeat / VNTR loci from the contig sequence itself
 *      (period 2..max_period, >= ~2.5 copies).
 *   3. Recover the reads that were assembled into that contig from the GFA
 *      `A` lines; keep the ones that contain both flanking regions (located by
 *      exact k-mer hits, so read orientation does not matter).
 *   4. Cut each read down to [left flank | repeat | right flank] and ask
 *      poa-consensus' single-allele flanked mode for a fresh repeat consensus.
 *   5. If the new repeat is better supported by the reads, splice it back
 *      between the same flank anchors.
 *
 * The original GFA is never modified: patched assemblies are written to
 * <prefix>.poa.gfa / <prefix>.poa.fa, together with a .poa.tsv report.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <zlib.h>

#include "kseq.h"
#include "CommandLines.h"
#include "poa_polish.h"
#include "poa-ffi/poa_ffi.h"

KSEQ_INIT(gzFile, gzread)

namespace {

const int TR_KMER = 11;          // k for period detection in the contig
const int HIT_KMER = 15;         // k for locating flanks inside reads
const int MAX_READS_PER_LOCUS = 128;
const double MIN_COPIES = 2.5;
const int MIN_TR_LEN = 60;

struct TrRegion {
    long long start;
    long long end;
    int period;
};

struct ReadSet {
    std::unordered_map<std::string, std::string> seq;
    std::vector<std::string> order;   // deterministic scan order
    size_t total_bases;
};

struct Segment {
    std::string name, seq, tags;
    std::vector<std::string> a_lines;
    std::vector<std::string> read_names;
    std::unordered_set<std::string> seen_reads;
};

struct LocusStat {
    std::string contig;
    long long start, end;
    int period;
    long long old_len, new_len;
    int n_reads;
    double cont_old, cont_new;
    std::string status;
};

/* ------------------------------------------------------------------------- */
/* small helpers                                                              */
/* ------------------------------------------------------------------------- */

inline int base_code(char c)
{
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

std::string revcomp(const std::string &s)
{
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
        case 'A': case 'a': r[i] = 'T'; break;
        case 'C': case 'c': r[i] = 'G'; break;
        case 'G': case 'g': r[i] = 'C'; break;
        case 'T': case 't': r[i] = 'A'; break;
        default:  r[i] = 'N'; break;
        }
    }
    return r;
}

bool has_suffix(const std::string &s, const std::string &suf)
{
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

bool file_exists(const std::string &fn)
{
    FILE *fp = fopen(fn.c_str(), "r");
    if (!fp) return false;
    fclose(fp);
    return true;
}

std::string strip_gfa_suffix(const std::string &fn)
{
    if (has_suffix(fn, ".gfa")) return fn.substr(0, fn.size() - 4);
    return fn;
}

void split_tabs(const std::string &s, std::vector<std::string> &out)
{
    out.clear();
    size_t p = 0;
    for (;;) {
        size_t q = s.find('\t', p);
        if (q == std::string::npos) {
            out.push_back(s.substr(p));
            break;
        }
        out.push_back(s.substr(p, q - p));
        p = q + 1;
    }
}

/**
 * Tolerant periodicity test at position i for period p, allowing the local
 * phase to drift by a couple of bases. A 1-2 bp indel inside a long VNTR is
 * common in an ONT contig and otherwise terminates repeat extension early,
 * leaving the tail of the array undetected.
 */
bool periodic_at(const std::string &s, long long i, int p, int &p_used)
{
    long long L = (long long)s.size();
    double best = 0.0;
    int best_p = p;
    for (int d = -3; d <= 3; ++d) {
        int pp = p + d;
        if (pp < 2 || i + pp >= L) continue;
        int m = 0, n = 0;
        for (int t = 0; t < 12 && i + t + pp < L; ++t, ++n)
            m += (s[i + t] == s[i + t + pp]);
        if (n < 6) continue;
        double f = (double)m / (double)n;
        if (f > best) { best = f; best_p = pp; }
    }
    p_used = best_p;
    return best >= 0.60;
}

/** Detect tandem repeats in a contig. Single pass with a packed-k-mer table. */
std::vector<TrRegion> detect_trs(const std::string &s, int max_period)
{
    std::vector<TrRegion> regions;
    const long long L = (long long)s.size();
    if (L < MIN_TR_LEN) return regions;

    const size_t TBL = (size_t)1 << (2 * TR_KMER);
    static std::vector<int> last;
    if (last.size() != TBL) last.assign(TBL, -1);
    else std::fill(last.begin(), last.end(), -1);

    const uint32_t mask = (uint32_t)(TBL - 1);
    uint32_t code = 0;
    int run = 0;                    // consecutive valid bases
    long long next_allowed = 0;

    for (long long i = 0; i < L; ++i) {
        int bc = base_code(s[i]);
        if (bc < 0) { run = 0; code = 0; continue; }
        code = ((code << 2) | (uint32_t)bc) & mask;
        ++run;
        if (run < TR_KMER) continue;

        int prev = last[code];
        last[code] = (int)i;
        if (prev < 0 || i < next_allowed) continue;

        int p = (int)(i - prev);
        if (p < 2 || p > max_period) continue;

        /* Compare at least two full units. The window must therefore scale
         * with p; a fixed 100 bp cap silently made every period > 50 bp
         * undetectable. */
        int W = std::min(2 * p + 50, 500);
        int m = 0, n = 0;
        for (int t = 0; t < W && i + t + p < L; ++t, ++n)
            m += (s[i + t] == s[i + t + p]);
        if (n < std::min(2 * p, W - 5) || m * 100 < 72 * n) continue;

        long long a = i;
        int pa = p;
        while (a > 0) {
            int np = pa;
            if (!periodic_at(s, a - 1, pa, np)) break;
            pa = np;
            --a;
        }
        long long b = i;
        int pb = p;
        while (b + p < L) {
            int np = pb;
            if (!periodic_at(s, b, pb, np)) break;
            pb = np;
            ++b;
        }
        long long e = std::min<long long>(L, b + pb);

        if (e - a >= MIN_TR_LEN && (e - a) >= (long long)(MIN_COPIES * p)) {
            regions.push_back(TrRegion{a, e, p});
            next_allowed = e;
        }
    }
    /* Collapse overlapping and near-adjacent calls: a 19 bp sub-period inside a
     * 38 bp unit, or the same array re-detected a few bases later, must not be
     * patched twice (the second splice would use stale coordinates). */
    std::sort(regions.begin(), regions.end(), [](const TrRegion &x, const TrRegion &y) {
        if (x.start != y.start) return x.start < y.start;
        return (x.end - x.start) > (y.end - y.start);
    });
    std::vector<TrRegion> merged;
    for (const TrRegion &r : regions) {
        if (!merged.empty() && r.start <= merged.back().end) {
            TrRegion &m = merged.back();
            long long ml = m.end - m.start, rl = r.end - r.start;
            if (std::abs(r.period - m.period) <= 2) {
                m.start = std::min(m.start, r.start);
                m.end = std::max(m.end, r.end);
                if (rl > ml) m.period = r.period;
            } else if (rl > ml) {
                m = r;
            }
        } else if (!merged.empty() &&
                   r.start - merged.back().end <= 2 * std::max(r.period, merged.back().period) + 30 &&
                   std::abs(r.period - merged.back().period) <= 2) {
            merged.back().end = r.end;
        } else {
            merged.push_back(r);
        }
    }
    return merged;
}

/**
 * Locate a repeat inside a read, requiring the FULL flanks on both sides.
 *
 * `l_start` / `l_end` are k-mer sets from the first / last 35 bp of the left
 * flank; `r_start` / `r_end` likewise for the right flank. A read only counts
 * when it reaches into the distal end of both flanks, so reads that start
 * inside the repeat (and would not let poa-consensus anchor the flank) are
 * rejected. `repeat_span` is the read's own repeat length between the flank
 * cores, used to keep a homogeneous allele set.
 */
bool locate_flanks(const std::string &read, int K,
                   const std::unordered_set<uint64_t> &l_start,
                   const std::unordered_set<uint64_t> &l_end,
                   const std::unordered_set<uint64_t> &r_start,
                   const std::unordered_set<uint64_t> &r_end,
                   long long min_span,
                   long long &lpos, long long &rend, long long &repeat_span)
{
    const long long n = (long long)read.size();
    if (n < K) return false;
    const uint64_t mask = ((1ULL << (2 * K)) - 1);
    uint64_t code = 0;
    int run = 0;
    long long ls = -1;           /* first distal-left hit */
    long long le = -1;           /* end of proximal-left hit */
    long long rs = -1;           /* start of proximal-right hit */
    long long re = -1;           /* end of distal-right hit */

    for (long long i = 0; i < n; ++i) {
        int b = base_code(read[i]);
        if (b < 0) { code = 0; run = 0; continue; }
        code = ((code << 2) | (uint64_t)b) & mask;
        ++run;
        if (run < K) continue;
        long long pos = i - (K - 1);
        if (l_start.count(code) && ls < 0) ls = pos;
        if (l_end.count(code)) le = pos + K;
        if (r_start.count(code) && rs < 0) rs = pos;
        if (r_end.count(code)) re = pos + K;
    }
    if (ls < 0 || le < 0 || rs < 0 || re < 0) return false;
    if (!(ls < le && le <= rs && rs < re)) return false;
    repeat_span = rs - le;
    if (repeat_span < min_span) return false;
    lpos = ls;
    rend = re;
    return true;
}

/** Distinct k-mer set of a sequence (non-canonical, ACGT only). */
std::unordered_set<uint64_t> kmer_set(const std::string &s, int K)
{
    std::unordered_set<uint64_t> out;
    const long long n = (long long)s.size();
    if (n < K) return out;
    const uint64_t mask = ((1ULL << (2 * K)) - 1);
    uint64_t code = 0;
    int run = 0;
    for (long long i = 0; i < n; ++i) {
        int b = base_code(s[i]);
        if (b < 0) { run = 0; code = 0; continue; }
        code = ((code << 2) | (uint64_t)b) & mask;
        if (++run >= K) out.insert(code);
    }
    return out;
}

/** Fraction of `seq`'s distinct k-mers present in `kmers`; -1 if too short. */
double containment(const std::string &seq, const std::unordered_set<uint64_t> &kmers, int K)
{
    if ((long long)seq.size() < K) return -1.0;
    std::unordered_set<uint64_t> mine = kmer_set(seq, K);
    if (mine.empty()) return -1.0;
    size_t hit = 0;
    for (uint64_t k : mine)
        if (kmers.count(k)) ++hit;
    return (double)hit / (double)mine.size();
}

/* ------------------------------------------------------------------------- */
/* reads                                                                      */
/* ------------------------------------------------------------------------- */

ReadSet load_reads(const hifiasm_opt_t *opt)
{
    ReadSet rs;
    rs.total_bases = 0;
    long long limit = 4000000000LL;
    if (const char *e = getenv("POA_POLISH_MAX_READ_BASES")) limit = atoll(e);

    for (int f = 0; f < opt->num_reads; ++f) {
        gzFile fp = gzopen(opt->read_file_names[f], "r");
        if (!fp) {
            fprintf(stderr, "[M::poa_polish] warning: cannot open %s\n", opt->read_file_names[f]);
            continue;
        }
        kseq_t *ks = kseq_init(fp);
        int l;
        while ((l = kseq_read(ks)) >= 0) {
            if (rs.total_bases > (size_t)limit) {
                fprintf(stderr,
                        "[M::poa_polish] read set exceeds %lld bases; raise "
                        "POA_POLISH_MAX_READ_BASES or disable --poa-polish\n", limit);
                kseq_destroy(ks);
                gzclose(fp);
                rs.seq.clear();
                rs.order.clear();
                return rs;
            }
            std::string nm(ks->name.s);
            if (rs.seq.emplace(nm, std::string(ks->seq.s, ks->seq.l)).second)
                rs.order.push_back(nm);
            rs.total_bases += ks->seq.l;
        }
        kseq_destroy(ks);
        gzclose(fp);
    }
    fprintf(stderr, "[M::poa_polish] loaded %zu reads (%.2f Mbp)\n",
            rs.seq.size(), rs.total_bases / 1e6);
    return rs;
}

/* ------------------------------------------------------------------------- */
/* one locus                                                                  */
/* ------------------------------------------------------------------------- */

struct CutRead {
    std::string seq;
};

bool locate_in_either_orientation(const std::string &read, int K,
                                  const std::unordered_set<uint64_t> &l_start,
                                  const std::unordered_set<uint64_t> &l_end,
                                  const std::unordered_set<uint64_t> &r_start,
                                  const std::unordered_set<uint64_t> &r_end,
                                  long long min_span,
                                  std::string &out,
                                  long long &lpos, long long &rend,
                                  long long &repeat_span)
{
    if (locate_flanks(read, K, l_start, l_end, r_start, r_end, min_span,
                      lpos, rend, repeat_span)) {
        out = read;
        return true;
    }
    std::string rc = revcomp(read);
    if (locate_flanks(rc, K, l_start, l_end, r_start, r_end, min_span,
                      lpos, rend, repeat_span)) {
        out.swap(rc);
        return true;
    }
    return false;
}

/**
 * Try to polish one repeat. Returns true when `seq` was modified.
 */
bool polish_locus(std::string &seq,
                  const TrRegion &tr,
                  const std::vector<std::string> &contig_reads,
                  const ReadSet &reads,
                  const hifiasm_opt_t *opt,
                  LocusStat &stat)
{
    const long long L = (long long)seq.size();
    int flank = opt->poa_flank_len > 0 ? opt->poa_flank_len : 120;
    if (flank < 40) flank = 40;
    flank = (int)std::min<long long>(flank, std::min<long long>(tr.start, L - tr.end));
    if (flank < 40) {
        stat.status = "skip:flank_too_short";
        return false;
    }

    const std::string left = seq.substr(tr.start - flank, flank);
    const std::string right = seq.substr(tr.end, flank);

    const int edge = std::min(35, flank);
    std::unordered_set<uint64_t> l_start = kmer_set(left.substr(0, edge), HIT_KMER);
    std::unordered_set<uint64_t> l_end = kmer_set(left.substr(flank - edge, edge), HIT_KMER);
    std::unordered_set<uint64_t> r_start = kmer_set(right.substr(0, edge), HIT_KMER);
    std::unordered_set<uint64_t> r_end = kmer_set(right.substr(flank - edge, edge), HIT_KMER);
    if (l_start.size() < 5 || l_end.size() < 5 ||
        r_start.size() < 5 || r_end.size() < 5) {
        stat.status = "skip:low_complexity_flank";
        return false;
    }

    const long long region_len = tr.end - tr.start;
    const long long min_span = (long long)(region_len * 0.5);
    const int pad = 10;

    std::vector<std::string> cut_reads;
    std::vector<long long> spans;
    cut_reads.reserve(std::min<size_t>(reads.order.size(), MAX_READS_PER_LOCUS));

    /* Prefer the reads listed in the GFA `A` lines (they are literally the ones
     * the assembler used at this contig), then scan the remaining input reads.
     * In hifiasm's ONT mode the A-line set is only a thin tiling, so the full
     * scan is what actually supplies depth. */
    std::unordered_set<std::string> seen;
    auto try_read = [&](const std::string &rn) {
        if (cut_reads.size() >= (size_t)MAX_READS_PER_LOCUS) return;
        if (!seen.insert(rn).second) return;
        auto it = reads.seq.find(rn);
        if (it == reads.seq.end()) return;
        std::string oriented;
        long long lpos = 0, rend = 0, span = 0;
        if (!locate_in_either_orientation(it->second, HIT_KMER,
                                          l_start, l_end, r_start, r_end,
                                          min_span, oriented, lpos, rend, span))
            return;
        long long a = std::max<long long>(0, lpos - pad);
        long long b = std::min<long long>((long long)oriented.size(), rend + pad);
        if (b - a < region_len * 0.6) return;
        cut_reads.push_back(oriented.substr(a, b - a));
        spans.push_back(span);
    };
    for (const std::string &rn : contig_reads) try_read(rn);
    for (const std::string &rn : reads.order) try_read(rn);

    /* Keep one allele's worth of reads: if the read set is a mixture (diploid
     * locus), cluster around the median repeat length instead of blending two
     * arrays into a longer mosaic. */
    if (spans.size() >= 3) {
        std::vector<long long> sorted = spans;
        std::sort(sorted.begin(), sorted.end());
        long long med = sorted[sorted.size() / 2];
        std::vector<std::string> keep_r;
        std::vector<long long> keep_s;
        for (size_t i = 0; i < cut_reads.size(); ++i) {
            if (spans[i] * 10 >= med * 7 && spans[i] * 10 <= med * 14) {
                keep_r.push_back(cut_reads[i]);
                keep_s.push_back(spans[i]);
            }
        }
        cut_reads.swap(keep_r);
        spans.swap(keep_s);
    }

    stat.n_reads = (int)cut_reads.size();
    long long span_med = region_len, span_tol = 40;
    if (!spans.empty()) {
        std::vector<long long> sorted = spans;
        std::sort(sorted.begin(), sorted.end());
        span_med = sorted[sorted.size() / 2];
        std::vector<long long> dev;
        dev.reserve(sorted.size());
        for (long long x : sorted) dev.push_back(std::llabs(x - span_med));
        std::sort(dev.begin(), dev.end());
        long long mad = dev[dev.size() / 2];
        span_tol = std::max<long long>(40, std::max(3 * mad, (long long)(0.02 * span_med)));
    }
    const int min_reads = std::max(3, opt->poa_min_reads);
    if ((int)cut_reads.size() < min_reads) {
        stat.status = "skip:insufficient_reads";
        return false;
    }

    std::unordered_set<uint64_t> read_kmers;
    for (const std::string &r : cut_reads) {
        std::unordered_set<uint64_t> s = kmer_set(r, HIT_KMER);
        read_kmers.insert(s.begin(), s.end());
    }
    const std::string old_repeat = seq.substr(tr.start, region_len);
    stat.cont_old = containment(old_repeat, read_kmers, HIT_KMER);

    std::vector<const char *> ptrs;
    std::vector<size_t> lens;
    ptrs.reserve(cut_reads.size());
    lens.reserve(cut_reads.size());
    for (const std::string &r : cut_reads) {
        ptrs.push_back(r.data());
        lens.push_back(r.size());
    }

    const char *dbg_dir = getenv("POA_POLISH_DEBUG_DIR");
    if (dbg_dir && *dbg_dir) {
        char fn[1024];
        snprintf(fn, sizeof(fn), "%s/locus_%s_%lld_%lld.reads.fa", dbg_dir,
                 stat.contig.c_str(), tr.start, tr.end);
        FILE *fp = fopen(fn, "w");
        if (fp) {
            for (size_t i = 0; i < cut_reads.size(); ++i)
                fprintf(fp, ">r%zu\n%s\n", i, cut_reads[i].c_str());
            fclose(fp);
        }
    }

    PoaOpts o;
    o.band_width = 0;                       /* crate default (50 + adaptive) */
    o.min_reads = (size_t)min_reads;
    o.no_homopolymer = 0;                   /* ONT: keep HPC + tandem gaps */
    o.max_read_len = 0;

    PoaResult *res = poa_consensus_flanked_c(ptrs.data(), lens.data(), ptrs.size(),
                                             left.data(), left.size(),
                                             right.data(), right.size(), &o);
    if (!res || res->status != POA_OK || !res->repeat) {
        stat.status = (res && res->error) ? (std::string("skip:") + res->error)
                                          : std::string("skip:poa_error");
        if (res) poa_result_free(res);
        return false;
    }

    const std::string new_repeat(res->repeat, res->repeat_len);
    if (dbg_dir && *dbg_dir) {
        char fn[1024];
        snprintf(fn, sizeof(fn), "%s/locus_%s_%lld_%lld.cons.fa", dbg_dir,
                 stat.contig.c_str(), tr.start, tr.end);
        FILE *fp = fopen(fn, "w");
        if (fp) {
            fprintf(fp, ">repeat_%zu\n%s\n", new_repeat.size(), new_repeat.c_str());
            if (res->whole)
                fprintf(fp, ">whole_%zu\n%s\n", res->whole_len, res->whole);
            fclose(fp);
        }
    }
    /* Both flanks must actually be locatable in the full consensus: if the
     * crate could not anchor, it returns the whole (flank-carrying) sequence
     * instead of a slice, and splicing that would duplicate or delete sequence
     * at the junctions. */
    bool anchored = false;
    if (res->whole && res->whole_len >= res->repeat_len + (size_t)(1.2 * flank)) {
        std::string whole(res->whole, res->whole_len);
        for (int probe : {20, 16, 12}) {
            if (probe > flank) continue;
            size_t pl = whole.find(left.substr(flank - probe, probe));
            size_t pr = whole.find(right.substr(0, probe));
            if (pl != std::string::npos && pr != std::string::npos && pr > pl) {
                anchored = true;
                break;
            }
        }
    }
    poa_result_free(res);

    if (!anchored) {
        stat.status = "skip:flank_not_anchored";
        return false;
    }

    double ratio = region_len > 0 ? (double)new_repeat.size() / (double)region_len : 0.0;
    if (new_repeat.empty() || ratio < 0.4 || ratio > 2.5) {
        stat.status = "skip:length_outlier";
        return false;
    }

    stat.cont_new = containment(new_repeat, read_kmers, HIT_KMER);
    if (stat.cont_new >= 0 && stat.cont_old >= 0 && stat.cont_new + 0.05 < stat.cont_old) {
        stat.status = "skip:lower_support";
        return false;
    }
    const long long delta = (long long)new_repeat.size() - region_len;
    if (new_repeat == old_repeat ||
        (std::llabs(delta) <= 2 && stat.cont_new <= stat.cont_old + 0.005)) {
        stat.status = "unchanged";
        return false;
    }

    /* Do not chase consensus noise: only rewrite when the current array length
     * is inconsistent with the length distribution of the reads that span it
     * (a real copy-number error), and the new length is closer to that median. */
    long long old_err = std::llabs(region_len - span_med);
    long long new_err = std::llabs((long long)new_repeat.size() - span_med);
    if (old_err <= span_tol || new_err >= old_err) {
        stat.status = "unchanged_read_length";
        return false;
    }

    seq = seq.substr(0, tr.start) + new_repeat + seq.substr(tr.end);
    stat.old_len = region_len;
    stat.new_len = (long long)new_repeat.size();
    stat.status = "patched";
    return true;
}

/* ------------------------------------------------------------------------- */
/* GFA walking                                                                */
/* ------------------------------------------------------------------------- */

std::string make_out_name(const std::string &gfa, const char *ext)
{
    return strip_gfa_suffix(gfa) + ext;
}

bool polish_one_gfa(const std::string &gfa_name,
                    const ReadSet &reads,
                    const hifiasm_opt_t *opt,
                    FILE *tsv,
                    int &n_loci,
                    int &n_patched,
                    long long &bases_added)
{
    FILE *in = fopen(gfa_name.c_str(), "r");
    if (!in) return false;

    std::string out_gfa_name = make_out_name(gfa_name, ".poa.gfa");
    std::string out_fa_name = make_out_name(gfa_name, ".poa.fa");
    FILE *out = fopen(out_gfa_name.c_str(), "w");
    FILE *fa = fopen(out_fa_name.c_str(), "w");
    if (!out || !fa) {
        fprintf(stderr, "[M::poa_polish] cannot write %s\n", out_gfa_name.c_str());
        if (out) fclose(out);
        if (fa) fclose(fa);
        fclose(in);
        return false;
    }

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    Segment cur;
    bool have_seg = false;
    int n_seg = 0, n_seg_patched = 0;
    std::vector<std::string> tk;

    auto flush_segment = [&]() {
        if (!have_seg) return;
        ++n_seg;
        bool patched = false;
        std::vector<std::pair<long long, long long>> shifts; // (region end, delta)

        if (cur.seq != "*" && cur.seq.size() >= (size_t)MIN_TR_LEN && !reads.order.empty()) {
            std::vector<TrRegion> trs = detect_trs(cur.seq, opt->poa_max_period);
            /* right-to-left keeps earlier coordinates valid */
            for (int k = (int)trs.size() - 1; k >= 0; --k) {
                LocusStat st;
                st.contig = cur.name;
                st.start = trs[k].start;
                st.end = trs[k].end;
                st.period = trs[k].period;
                st.old_len = st.new_len = trs[k].end - trs[k].start;
                st.n_reads = 0;
                st.cont_old = st.cont_new = -1.0;
                ++n_loci;
                if (polish_locus(cur.seq, trs[k], cur.read_names, reads, opt, st)) {
                    patched = true;
                    ++n_patched;
                    bases_added += st.new_len - st.old_len;
                    shifts.push_back({trs[k].end, st.new_len - st.old_len});
                }
                fprintf(stderr,
                        "[M::poa_polish] %s %lld-%lld period=%d copies=%.1f reads=%d "
                        "support=%.3f->%.3f %s\n",
                        st.contig.c_str(), st.start, st.end, st.period,
                        (double)(trs[k].end - trs[k].start) / trs[k].period,
                        st.n_reads, st.cont_old, st.cont_new, st.status.c_str());
                if (tsv)
                    fprintf(tsv,
                            "%s\t%s\t%lld\t%lld\t%d\t%.2f\t%lld\t%lld\t%d\t%.4f\t%.4f\t%s\n",
                            gfa_name.c_str(), st.contig.c_str(), st.start, st.end,
                            st.period,
                            (double)(trs[k].end - trs[k].start) / trs[k].period,
                            st.old_len, st.new_len, st.n_reads,
                            st.cont_old, st.cont_new, st.status.c_str());
            }
            if (patched) ++n_seg_patched;
        }

        /* S line: update LN:i: but preserve every other tag */
        fprintf(out, "S\t%s\t%s", cur.name.c_str(), cur.seq.c_str());
        split_tabs(cur.tags, tk);
        std::string rebuilt;
        bool had_ln = false;
        for (const std::string &t : tk) {
            if (t.rfind("LN:i:", 0) == 0) { had_ln = true; continue; }
            if (t.empty()) continue;
            if (!rebuilt.empty()) rebuilt += '\t';
            rebuilt += t;
        }
        if (!had_ln && cur.seq != "*") {
            char buf[64];
            snprintf(buf, sizeof(buf), "LN:i:%lld", (long long)cur.seq.size());
            if (!rebuilt.empty()) rebuilt += '\t';
            rebuilt += buf;
        }
        if (!rebuilt.empty()) fprintf(out, "\t%s", rebuilt.c_str());
        fprintf(out, "\n");

        /* Original `A` lines, with contig offsets shifted by the splices that
         * happened upstream of them (otherwise every annotation after a
         * lengthened/shortened array would point at the wrong base). */
        for (const std::string &a : cur.a_lines) {
            if (shifts.empty()) {
                fprintf(out, "%s\n", a.c_str());
                continue;
            }
            split_tabs(a, tk);
            if (tk.size() >= 3) {
                char *endp = NULL;
                long long off = strtoll(tk[2].c_str(), &endp, 10);
                if (endp && *endp == '\0') {
                    long long shift = 0;
                    for (const auto &sh : shifts)
                        if (off >= sh.first) shift += sh.second;
                    char buf[32];
                    snprintf(buf, sizeof(buf), "%lld", off + shift);
                    tk[2] = buf;
                }
            }
            std::string rebuilt;
            for (size_t ti = 0; ti < tk.size(); ++ti) {
                if (ti) rebuilt += '\t';
                rebuilt += tk[ti];
            }
            fprintf(out, "%s\n", rebuilt.c_str());
        }

        if (cur.seq != "*") {
            fprintf(fa, ">%s\n", cur.name.c_str());
            for (size_t i = 0; i < cur.seq.size(); i += 60)
                fprintf(fa, "%.*s\n", (int)std::min<size_t>(60, cur.seq.size() - i),
                        cur.seq.c_str() + i);
        }

        cur = Segment();
        have_seg = false;
    };

    while ((len = getline(&line, &cap, in)) > 0) {
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len == 0) continue;

        if (line[0] == 'S' && line[1] == '\t') {
            flush_segment();
            std::string s(line, len);
            size_t t1 = s.find('\t');
            size_t t2 = (t1 == std::string::npos) ? std::string::npos : s.find('\t', t1 + 1);
            size_t t3 = (t2 == std::string::npos) ? std::string::npos : s.find('\t', t2 + 1);
            if (t2 == std::string::npos) continue;
            cur = Segment();
            cur.name = s.substr(t1 + 1, t2 - t1 - 1);
            if (t3 == std::string::npos) {
                cur.seq = s.substr(t2 + 1);
                cur.tags.clear();
            } else {
                cur.seq = s.substr(t2 + 1, t3 - t2 - 1);
                cur.tags = s.substr(t3 + 1);
            }
            have_seg = true;
        } else if (line[0] == 'A' && line[1] == '\t') {
            if (!have_seg) continue;
            split_tabs(std::string(line, len), tk);
            if (tk.size() >= 5 && tk[4] != "FAKE" && tk[4] != "*") {
                if (cur.seen_reads.insert(tk[4]).second) cur.read_names.push_back(tk[4]);
            }
            cur.a_lines.push_back(std::string(line, len));
        } else {
            flush_segment();
            fprintf(out, "%s\n", line);
        }
    }
    flush_segment();

    free(line);
    fclose(in);
    fclose(out);
    fclose(fa);

    fprintf(stderr, "[M::poa_polish] %s: %d contig(s), %d polished\n",
            gfa_name.c_str(), n_seg, n_seg_patched);
    return true;
}

} // namespace

/* ------------------------------------------------------------------------- */

int poa_polish_run(const hifiasm_opt_t *opt)
{
    if (!opt || !opt->is_poa_polish) return 0;
    if (!opt->output_file_name || !opt->read_file_names || opt->num_reads == 0) {
        fprintf(stderr, "[M::poa_polish] no reads/output; skipping\n");
        return 1;
    }

    fprintf(stderr, "[M::poa_polish] post-assembly tandem-repeat polishing (%s)\n",
            poa_version());

    ReadSet reads = load_reads(opt);
    if (reads.seq.empty()) {
        fprintf(stderr, "[M::poa_polish] no reads loaded; skipping\n");
        return 1;
    }

    std::vector<std::string> candidates;
    const std::string out = opt->output_file_name;
    const char *suffixes[] = {".p_ctg.gfa", ".a_ctg.gfa",
                              ".bp.p_ctg.gfa", ".bp.hap1.p_ctg.gfa",
                              ".bp.hap2.p_ctg.gfa",
                              ".hap1.p_ctg.gfa", ".hap2.p_ctg.gfa",
                              ".stg.p_ctg.gfa"};
    for (const char *suf : suffixes) {
        std::string fn = out + suf;
        if (file_exists(fn)) candidates.push_back(fn);
    }
    if (candidates.empty()) {
        fprintf(stderr, "[M::poa_polish] no contig GFA found for prefix %s; skipping\n",
                out.c_str());
        return 1;
    }

    std::string tsv_name = out + ".poa.tsv";
    FILE *tsv = fopen(tsv_name.c_str(), "w");
    if (tsv)
        fprintf(tsv, "#gfa\tcontig\tstart\tend\tperiod\tcopies\told_len\tnew_len\treads\t"
                     "support_old\tsupport_new\tstatus\n");

    int n_loci = 0, n_patched = 0;
    long long bases_added = 0;
    for (const std::string &fn : candidates)
        polish_one_gfa(fn, reads, opt, tsv, n_loci, n_patched, bases_added);

    if (tsv) fclose(tsv);

    fprintf(stderr,
            "[M::poa_polish] done: %d locus/loci found, %d patched, %+lld bp; report: %s\n",
            n_loci, n_patched, bases_added, tsv_name.c_str());
    return 0;
}
