/**
 * poa_polish.cpp -- post-assembly VNTR/tandem-repeat polishing for hifiasm-fast.
 *
 * Implements the workflow suggested by the poa-consensus author:
 *
 *   1. After the contig graph is written, walk every contig in the GFA.
 *   2. Find tandem-repeat / VNTR loci from the contig sequence itself
 *      (period 2..max_period, >= ~2.5 copies).
 *   3. Recover the reads that span each locus plus BOTH flanks. This is done in
 *      a single streaming pass over the input reads using a minimizer index
 *      built from the locus flanks only -- the read set is never held in
 *      memory, only the (small) set of reads that actually span a repeat.
 *   4. Cut each recruited read down to [left flank | repeat | right flank] and
 *      ask poa-consensus' single-allele flanked mode for a fresh repeat
 *      consensus.
 *   5. If the new repeat is better supported by the reads, splice it back
 *      between the same flank anchors.
 *
 * Memory therefore scales with the number of repeat loci and the reads that
 * span them, not with the size of the input read set. The original GFA is
 * never modified: patched assemblies are written to <prefix>.poa.gfa /
 * <prefix>.poa.fa, together with a .poa.tsv report.
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#include "kseq.h"
#include "CommandLines.h"
#include "poa_polish.h"
#include "poa-ffi/poa_ffi.h"

/*
 * One reader for plain files, gzip via zlib, and gzip via an external
 * rapidgzip process (parallel decompression of a single gzip stream).
 */
struct InputStream {
    int kind;      /* 0 = gzFile, 1 = FILE* (plain file or pipe) */
    gzFile gz;
    FILE *fp;
    pid_t pid;
};

static int input_read(InputStream *in, unsigned char *buf, int len)
{
    if (in->kind == 0) return gzread(in->gz, buf, len);
    return (int)fread(buf, 1, len, in->fp);
}

#define KSEQ_INIT_BUF(type_t, __read, __bufsize) \
    KSTREAM_INIT(type_t, __read, __bufsize) \
    __KSEQ_TYPE(type_t) \
    __KSEQ_BASIC(static klib_unused, type_t) \
    __KSEQ_READ(static klib_unused)
KSEQ_INIT_BUF(InputStream*, input_read, 1048576)

namespace {

const int TR_KMER = 11;            // k for period detection in the contig
const int RECRUIT_K = 17;          // seed k-mer for read recruitment
/* Check every Nth k-mer of a read. 2 keeps ~10 seeds per flank edge and
 * recovers nearly all spanning reads; raise it with --poa-seed-step for speed
 * at the cost of sensitivity. */
int g_seed_step = 2;
int g_edge_len = 60;               // bp at each flank end used for anchoring
const int HIT_KMER = 15;           // k for the post-hoc k-mer support check
const int MAX_READS_PER_LOCUS = 128;
const double MIN_COPIES = 2.5;
const int MIN_TR_LEN = 60;

struct TrRegion {
    long long start;
    long long end;
    int period;
};

/**
 * One repeat locus plus the reads recruited for it. Only these reads are held
 * in memory; everything else is discarded as it streams past.
 */
struct Locus {
    std::string contig;
    long long start, end;
    int period;
    int flank;
    std::string left, right;
    std::vector<std::string> cut_reads;
    std::vector<long long> spans;
};

static double now_s()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct LocusStat {
    std::string contig;
    long long start, end;
    int period;
    long long old_len, new_len;
    int n_reads;
    double cont_old, cont_new;
    long long span_med = 0, span_min = 0, span_max = 0;
    std::string status;
};

struct Mini {
    uint64_t kmer;
    int pos;
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

/* ------------------------------------------------------------------------- */
/* tandem-repeat detection                                                    */
/* ------------------------------------------------------------------------- */

/** Tolerant periodicity test at position i for period p (phase may drift). */
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
std::vector<TrRegion> detect_trs(const std::string &s, int max_period, int min_len)
{
    std::vector<TrRegion> regions;
    const long long L = (long long)s.size();
    if (L < min_len) return regions;

    const size_t TBL = (size_t)1 << (2 * TR_KMER);
    static std::vector<int> last;
    if (last.size() != TBL) last.assign(TBL, -1);
    else std::fill(last.begin(), last.end(), -1);

    const uint32_t mask = (uint32_t)(TBL - 1);
    uint32_t code = 0;
    int run = 0;
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

        /* Compare at least two full units: the window must scale with p or
         * every period > 50 bp becomes undetectable. */
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

        if (e - a >= min_len && (e - a) >= (long long)(MIN_COPIES * p)) {
            regions.push_back(TrRegion{a, e, p});
            next_allowed = e;
        }
    }

    /* Collapse overlapping and near-adjacent calls (a sub-period inside a
     * longer unit, or the same array re-detected). */
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

/* ------------------------------------------------------------------------- */
/* minimizers + seed index                                                    */
/* ------------------------------------------------------------------------- */

/**
 * kmer -> list of (locus, strand, flank side). Strand is 0 for the contig
 * orientation and 1 when the indexed k-mer is the reverse complement of the
 * flank edge; side is 0..3 = left-start, left-end, right-start, right-end.
 */
struct SeedHit {
    uint64_t packed;   // locus<<15 | strand<<14 | side<<12 | offset
    uint32_t next;
};

struct SeedIndex {
    std::unordered_map<uint64_t, uint32_t> head;
    std::vector<SeedHit> hits;

    void clear()
    {
        head.clear();
        hits.clear();
    }

    void add(uint64_t kmer, uint64_t locus, int strand, int side, int offset)
    {
        uint32_t &h = head.emplace(kmer, UINT32_MAX).first->second;
        uint32_t idx = (uint32_t)hits.size();
        hits.push_back(SeedHit{(locus << 15) | ((uint64_t)strand << 14) |
                               ((uint64_t)side << 12) | (uint64_t)offset, h});
        h = idx;
    }
};

std::string flank_edge(const Locus &l, int side)
{
    const int F = l.flank;
    const int E = std::min(g_edge_len, F);
    switch (side) {
    case 0: return l.left.substr(0, E);
    case 1: return l.left.substr(F - E, E);
    case 2: return l.right.substr(0, E);
    default: return l.right.substr(F - E, E);
    }
}

/**
 * Index every k-mer of the flank edges (not just their minimizers): the read
 * side is sampled with minimizers, and a minimizer of the *read* is only
 * guaranteed to be one of the flank's k-mers if the whole edge is indexed
 * (minimizers are not closed under context extension).
 */
bool low_complexity_kmer(uint64_t code, int k)
{
    int prev = -1, run = 1, maxrun = 1;
    unsigned seen = 0;
    int distinct = 0;
    for (int j = 0; j < k; ++j) {
        int b = (int)((code >> (2 * (k - 1 - j))) & 3);
        if (!(seen & (1u << b))) { seen |= 1u << b; ++distinct; }
        if (b == prev) {
            ++run;
            if (run > maxrun) maxrun = run;
        } else {
            run = 1;
            prev = b;
        }
    }
    return distinct <= 2 || maxrun >= 10;
}

void build_index(const std::vector<Locus> &loci, SeedIndex &idx,
                 size_t begin, size_t end)
{
    idx.clear();
    size_t n = end > begin ? end - begin : 0;
    size_t seeds_per_locus = 8 * (size_t)std::max(1, g_edge_len - RECRUIT_K + 1);
    idx.head.reserve(n * seeds_per_locus + 16);
    idx.hits.reserve(n * seeds_per_locus + 64);
    for (size_t i = begin; i < end; ++i) {
        for (int side = 0; side < 4; ++side) {
            std::string edge = flank_edge(loci[i], side);
            for (int strand = 0; strand < 2; ++strand) {
                const std::string s = (strand == 0) ? edge : revcomp(edge);
                std::unordered_set<uint64_t> seen;
                for (int p = 0; p + RECRUIT_K <= (int)s.size(); ++p) {
                    uint64_t code = 0;
                    bool ok = true;
                    for (int j = 0; j < RECRUIT_K; ++j) {
                        int b = base_code(s[p + j]);
                        if (b < 0) { ok = false; break; }
                        code = (code << 2) | (uint64_t)b;
                    }
                    /* keep the first offset for a repeated k-mer within an
                     * edge; skip non-discriminative low-complexity seeds */
                    if (ok && !low_complexity_kmer(code, RECRUIT_K) &&
                        seen.insert(code).second)
                        idx.add(code, i, strand, side, p);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/* read recruitment (producer / worker pipeline)                              */
/* ------------------------------------------------------------------------- */

struct FoundHit {
    uint32_t locus;
    int pos;
    int offset;
    uint8_t strand;
    uint8_t side;
};

struct CutRecord {
    uint64_t key;
    uint32_t locus;
    uint8_t strand;
    std::string seq;
    long long span;
};

struct RecruitCtx {
    std::vector<FoundHit> found;
    std::string scratch;
};

struct RecruitStats {
    long long reads;
    long long bases;
    long long stored_reads;
    long long stored_bases;
    bool capped;
};

/**
 * Scan a read against the locus seed index and append one CutRecord per locus
 * the read spans with both full flanks. Thread-safe: all mutable state is in
 * `ctx`/`out` (one per worker).
 */
void recruit_read(const std::string &read, uint64_t key,
                  const std::vector<Locus> &loci,
                  const SeedIndex &idx,
                  RecruitCtx &ctx,
                  std::vector<CutRecord> &out,
                  size_t cap_bases,
                  std::atomic<size_t> *soft_bases,
                  std::atomic<bool> *soft_capped)
{
    const long long n = (long long)read.size();
    std::vector<FoundHit> &found = ctx.found;
    std::string &scratch = ctx.scratch;
    found.clear();

    const uint64_t mask = (1ULL << (2 * RECRUIT_K)) - 1;
    uint64_t code = 0;
    int run = 0;
    for (long long p = 0; p < n; ++p) {
        int b = base_code(read[p]);
        if (b < 0) { code = 0; run = 0; continue; }
        code = ((code << 2) | (uint64_t)b) & mask;
        ++run;
        if (run < RECRUIT_K) continue;
        int start = (int)(p - RECRUIT_K + 1);
        if (g_seed_step > 1 && (start % g_seed_step)) continue;
        auto it = idx.head.find(code);
        if (it == idx.head.end()) continue;
        int chain = 0;
        for (uint32_t h = it->second; h != UINT32_MAX; h = idx.hits[h].next) {
            if (++chain > 512) break;  /* non-discriminative k-mer */
            uint64_t q = idx.hits[h].packed;
            FoundHit f;
            f.locus = (uint32_t)(q >> 15);
            f.strand = (uint8_t)((q >> 14) & 1);
            f.side = (uint8_t)((q >> 12) & 3);
            f.offset = (int)(q & 0xFFF);
            f.pos = start;
            found.push_back(f);
            if (found.size() > 500000) { found.clear(); return; }
        }
    }
    if (found.empty()) return;

    std::sort(found.begin(), found.end(), [](const FoundHit &a, const FoundHit &b) {
        if (a.locus != b.locus) return a.locus < b.locus;
        if (a.strand != b.strand) return a.strand < b.strand;
        return a.pos < b.pos;
    });

    bool have_rc = false;
    size_t i = 0;
    while (i < found.size()) {
        size_t j = i + 1;
        while (j < found.size() && found[j].locus == found[i].locus &&
               found[j].strand == found[i].strand)
            ++j;

        const int strand = found[i].strand;
        const uint32_t lid = found[i].locus;
        const Locus &loc = loci[lid];
        const long long region_len = loc.end - loc.start;
        const int E = std::min(g_edge_len, loc.flank);
        std::vector<long long> ls_v, le_v, rs_v, re_v;
        int l_start_min_off = INT_MAX;
        int r_end_max_off = -1;
        int n_l_end = 0, n_r_start = 0;

        for (size_t t = i; t < j; ++t) {
            int pos = found[t].pos;
            int off = found[t].offset;
            if (strand == 1) {
                pos = (int)(n - pos - RECRUIT_K);
                off = E - off - RECRUIT_K;
            }
            switch (found[t].side) {
            case 0:
                ls_v.push_back(pos - off);
                l_start_min_off = std::min(l_start_min_off, off);
                break;
            case 1:
                le_v.push_back(pos + E - off);
                ++n_l_end;
                break;
            case 2:
                rs_v.push_back(pos - off);
                ++n_r_start;
                break;
            default:
                re_v.push_back(pos + E - off);
                r_end_max_off = std::max(r_end_max_off, off);
                break;
            }
        }
        i = j;

        if (ls_v.empty() || re_v.empty() || le_v.empty() || rs_v.empty()) continue;
        if (n < (long long)(0.7 * region_len)) continue;
        if (l_start_min_off > 20 || r_end_max_off < E - RECRUIT_K - 20) continue;
        if (n_l_end < 1 || n_r_start < 1) continue;

        auto median = [](std::vector<long long> &v) {
            std::sort(v.begin(), v.end());
            return v[v.size() / 2];
        };
        const long long ls = median(ls_v);
        const long long le = median(le_v);
        const long long rs = median(rs_v);
        const long long re = median(re_v);
        if (!(ls < le && le <= rs && rs < re)) continue;
        long long span = rs - le;

        if (strand == 1 && !have_rc) { scratch = revcomp(read); have_rc = true; }
        const std::string &rd = (strand == 0) ? read : scratch;

        int lprobe = -1, loff = -1, rprobe = -1, roff = -1;
        for (int off : {loc.flank - 16, loc.flank - 24, loc.flank - 32, loc.flank - 40}) {
            if (off < 0 || off + 16 > (int)loc.left.size()) continue;
            size_t p = rd.find(loc.left.substr(off, 16));
            if (p != std::string::npos) { lprobe = (int)p; loff = off; break; }
        }
        for (int off : {0, 8, 16, 24}) {
            if (off + 16 > (int)loc.right.size()) continue;
            size_t p = rd.rfind(loc.right.substr(off, 16));
            if (p != std::string::npos) { rprobe = (int)p; roff = off; break; }
        }

        long long a = ls - 10;
        long long b = re + 10;
        if (lprobe >= 0 && rprobe > lprobe) {
            const long long pspan = (long long)(rprobe - roff) -
                                    ((long long)lprobe - loff + loc.flank);
            if (pspan >= (long long)(0.2 * region_len)) {
                span = pspan;
                a = (long long)lprobe - loff - 10;
                b = (long long)rprobe + loc.flank - roff + 10;
            }
        }
        if (span > (long long)(5.0 * region_len) + 2000) continue;
        if (a < 0) a = 0;
        if (b > n) b = n;
        if (b - a < span + 2 * (long long)loc.flank - 40) continue;
        if (b - a < 100) continue;

        const size_t len = (size_t)(b - a);
        if (cap_bases != std::numeric_limits<size_t>::max()) {
            size_t cur = soft_bases->load(std::memory_order_relaxed);
            if (cur + len > cap_bases) {
                soft_capped->store(true, std::memory_order_relaxed);
                return;
            }
            soft_bases->fetch_add(len, std::memory_order_relaxed);
        }
        out.push_back(CutRecord{key, lid, (uint8_t)strand,
                                rd.substr(a, len), span});
    }
}

/** One batch of reads handed from a producer to a worker. */
struct ReadBatch {
    std::vector<std::string> reads;
    std::vector<uint64_t> keys;
    size_t bases = 0;
};

struct BatchQueue {
    std::deque<ReadBatch> q;
    std::mutex mtx;
    std::condition_variable not_empty, not_full;
    std::vector<CutRecord> results;
    std::mutex results_mtx;
    std::atomic<long long> reads{0};
    std::atomic<long long> bases{0};
    std::atomic<size_t> soft_bases{0};
    std::atomic<bool> soft_capped{false};
    std::atomic<int> next_file{0};
    int producers_left = 0;
    bool done = false;
    size_t max_batches = 16;
    std::mutex log_mtx;
};


/* ------------------------------------------------------------------------- */
/* input backends: plain file, zlib, rapidgzip pipe                          */
/* ------------------------------------------------------------------------- */

bool path_is_gzip(const char *path)
{
    size_t n = strlen(path);
    if (n >= 3 && !strcmp(path + n - 3, ".gz")) return true;
    if (n >= 4 && !strcmp(path + n - 4, ".bgz")) return true;
    FILE *fp = fopen(path, "rb");
    if (!fp) return false;
    unsigned char m[2] = {0, 0};
    size_t got = fread(m, 1, 2, fp);
    fclose(fp);
    return got == 2 && m[0] == 0x1f && m[1] == 0x8b;
}

bool command_available(const char *name)
{
    const char *path = getenv("PATH");
    if (!path || !*path) return false;
    std::string all(path);
    size_t pos = 0;
    for (;;) {
        size_t q = all.find(':', pos);
        std::string dir = all.substr(pos, q == std::string::npos ? std::string::npos
                                                                : q - pos);
        if (dir.empty()) dir = ".";
        std::string fn = dir + "/" + name;
        if (access(fn.c_str(), X_OK) == 0) return true;
        if (q == std::string::npos) break;
        pos = q + 1;
    }
    return false;
}

/**
 * Open a read file. `decompressor`: 0 = auto, 1 = zlib, 2 = rapidgzip.
 * rapidgzip is spawned as `rapidgzip -d -c -P N <file>` and read through a
 * pipe; its output is ordinary FASTQ, so the parser is unchanged.
 */
InputStream *open_input(const char *path, int decompressor, int gzip_threads)
{
    InputStream *in = new InputStream;
    in->kind = 0;
    in->gz = NULL;
    in->fp = NULL;
    in->pid = -1;

    if (!path_is_gzip(path)) {
        in->fp = fopen(path, "rb");
        if (!in->fp) { delete in; return NULL; }
        in->kind = 1;
        return in;
    }

    bool want_rapidgzip = (decompressor == 2) ||
                          (decompressor == 0 && command_available("rapidgzip"));
    if (want_rapidgzip) {
        int fds[2];
        if (pipe(fds) == 0) {
            pid_t pid = fork();
            if (pid == 0) {
                close(fds[0]);
                if (dup2(fds[1], STDOUT_FILENO) < 0) _exit(127);
                close(fds[1]);
                char pbuf[32];
                snprintf(pbuf, sizeof(pbuf), "%d", gzip_threads > 0 ? gzip_threads : 0);
                execlp("rapidgzip", "rapidgzip", "-d", "-c", "-P", pbuf, path,
                       (char *)NULL);
                _exit(127);
            }
            if (pid > 0) {
                close(fds[1]);
                FILE *fp = fdopen(fds[0], "r");
                if (fp) {
                    in->kind = 1;
                    in->fp = fp;
                    in->pid = pid;
                    return in;
                }
                close(fds[0]);
                int st = 0;
                waitpid(pid, &st, 0);
            } else {
                close(fds[0]);
                close(fds[1]);
            }
        }
        if (decompressor == 2)
            fprintf(stderr, "[M::poa_polish] warning: could not start rapidgzip; "
                            "falling back to zlib\n");
    }

    in->gz = gzopen(path, "rb");
    if (!in->gz) { delete in; return NULL; }
    gzbuffer(in->gz, 1 << 20);
    in->kind = 0;
    return in;
}

void close_input(InputStream *in)
{
    if (!in) return;
    if (in->kind == 0) {
        if (in->gz) gzclose(in->gz);
    } else {
        if (in->fp) fclose(in->fp);
        if (in->pid > 0) {
            int st = 0;
            waitpid(in->pid, &st, 0);
            if (WIFEXITED(st) && WEXITSTATUS(st) != 0)
                fprintf(stderr, "[M::poa_polish] warning: rapidgzip exited with status %d\n",
                        WEXITSTATUS(st));
        }
    }
    delete in;
}

struct PipelineCfg {
    const std::vector<Locus> *loci;
    const SeedIndex *idx;
    size_t cap_bases;
    size_t batch_reads;
    size_t batch_bases;
    int decompressor;
    int gzip_threads;
};


void push_batch(BatchQueue &bq, ReadBatch &batch)
{
    if (batch.reads.empty()) return;
    std::unique_lock<std::mutex> lk(bq.mtx);
    bq.not_full.wait(lk, [&] { return bq.q.size() < bq.max_batches; });
    bq.q.push_back(std::move(batch));
    batch = ReadBatch();
    bq.not_empty.notify_one();
}

void producer(const hifiasm_opt_t *opt, BatchQueue &bq, const PipelineCfg &cfg)
{
    for (;;) {
        if (bq.soft_capped.load(std::memory_order_relaxed)) break;
        int f = bq.next_file.fetch_add(1);
        if (f >= opt->num_reads) break;

        InputStream *in = open_input(opt->read_file_names[f], cfg.decompressor,
                                     cfg.gzip_threads);
        if (!in) {
            fprintf(stderr, "[M::poa_polish] warning: cannot open %s\n",
                    opt->read_file_names[f]);
            continue;
        }
        kseq_t *ks = kseq_init(in);
        ReadBatch batch;
        uint64_t read_idx = 0;
        int l;
        while ((l = kseq_read(ks)) >= 0) {
            if (bq.soft_capped.load(std::memory_order_relaxed)) break;
            batch.keys.push_back(((uint64_t)f << 40) | read_idx++);
            batch.reads.emplace_back(ks->seq.s, ks->seq.l);
            batch.bases += (size_t)l;
            long long nr = bq.reads.fetch_add(1) + 1;
            bq.bases.fetch_add(l);
            if (nr % 1000000 == 0) {
                std::lock_guard<std::mutex> lk(bq.log_mtx);
                fprintf(stderr,
                        "[M::poa_polish] streamed %lld reads (%.2f Gbp), kept %.2f Mbp\n",
                        nr, bq.bases.load() / 1e9, bq.soft_bases.load() / 1e6);
            }
            if (batch.reads.size() >= cfg.batch_reads ||
                batch.bases >= cfg.batch_bases)
                push_batch(bq, batch);
        }
        kseq_destroy(ks);
        close_input(in);
        push_batch(bq, batch);
    }

    std::unique_lock<std::mutex> lk(bq.mtx);
    if (--bq.producers_left == 0) {
        bq.done = true;
        bq.not_empty.notify_all();
    }
}

void worker(BatchQueue &bq, const PipelineCfg &cfg)
{
    RecruitCtx ctx;
    std::vector<CutRecord> local;
    for (;;) {
        ReadBatch batch;
        {
            std::unique_lock<std::mutex> lk(bq.mtx);
            bq.not_empty.wait(lk, [&] { return !bq.q.empty() || bq.done; });
            if (bq.q.empty() && bq.done) break;
            batch = std::move(bq.q.front());
            bq.q.pop_front();
            bq.not_full.notify_one();
        }
        for (size_t i = 0; i < batch.reads.size(); ++i)
            recruit_read(batch.reads[i], batch.keys[i], *cfg.loci, *cfg.idx, ctx,
                         local, cfg.cap_bases, &bq.soft_bases, &bq.soft_capped);
        if (!local.empty()) {
            std::lock_guard<std::mutex> lk(bq.results_mtx);
            for (CutRecord &r : local) bq.results.push_back(std::move(r));
            local.clear();
        }
    }
}

/** Merge per-worker records in read order so results match a serial pass. */
void merge_records(std::vector<CutRecord> &recs, std::vector<Locus> &loci,
                   size_t cap_bases, RecruitStats &st)
{
    std::sort(recs.begin(), recs.end(), [](const CutRecord &a, const CutRecord &b) {
        if (a.key != b.key) return a.key < b.key;
        if (a.locus != b.locus) return a.locus < b.locus;
        return a.strand < b.strand;
    });
    std::vector<uint32_t> counts(loci.size(), 0);
    size_t total = 0;
    for (CutRecord &r : recs) {
        if (counts[r.locus] >= (uint32_t)MAX_READS_PER_LOCUS) continue;
        if (cap_bases != std::numeric_limits<size_t>::max() &&
            total + r.seq.size() > cap_bases) {
            st.capped = true;
            break;
        }
        total += r.seq.size();
        ++counts[r.locus];
        st.stored_bases += (long long)r.seq.size();
        ++st.stored_reads;
        loci[r.locus].cut_reads.push_back(std::move(r.seq));
        loci[r.locus].spans.push_back(r.span);
    }
}

RecruitStats stream_recruit(const hifiasm_opt_t *opt, std::vector<Locus> &loci,
                            const SeedIndex &idx, size_t cap_total,
                            size_t already_stored)
{
    RecruitStats st{0, 0, 0, 0, false};
    size_t cap_bases = std::numeric_limits<size_t>::max();
    if (cap_total != std::numeric_limits<size_t>::max()) {
        if (already_stored >= cap_total) {
            st.capped = true;
            return st;
        }
        cap_bases = cap_total - already_stored;
    }

    int nthreads = opt->thread_num > 0 ? opt->thread_num : 1;

    int nprod = std::min(opt->num_reads, 4);
    if (nprod < 1) nprod = 1;

    bool have_rapidgzip = command_available("rapidgzip");
    int decompressor = opt->poa_decompressor; /* 0 auto, 1 zlib, 2 rapidgzip */
    if (decompressor == 2 && !have_rapidgzip) {
        fprintf(stderr, "[M::poa_polish] warning: rapidgzip not found in PATH; "
                        "using zlib\n");
        decompressor = 1;
    }
    if (decompressor == 0) decompressor = have_rapidgzip ? 2 : 1;
    int gzip_threads = opt->poa_gzip_threads > 0
                           ? opt->poa_gzip_threads
                           : std::max(1, nthreads / (2 * nprod));

    PipelineCfg cfg;
    cfg.loci = &loci;
    cfg.idx = &idx;
    cfg.cap_bases = cap_bases;
    cfg.batch_reads = 256;
    cfg.batch_bases = 2 * 1024 * 1024;
    cfg.decompressor = decompressor;
    cfg.gzip_threads = gzip_threads;

    static bool logged_decompressor = false;
    if (!logged_decompressor) {
        std::string dmsg = decompressor == 2
            ? "rapidgzip (-P " + std::to_string(gzip_threads) + " per file; plain files bypass it)"
            : "zlib (plain files bypass it)";
        fprintf(stderr, "[M::poa_polish] decompressor: %s\n", dmsg.c_str());
        logged_decompressor = true;
    }

    BatchQueue bq;
    bq.max_batches = (size_t)std::max(8, 4 * nthreads);
    bq.producers_left = nprod;

    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) workers.emplace_back(worker, std::ref(bq), std::cref(cfg));

    std::vector<std::thread> prods;
    prods.reserve(nprod);
    for (int t = 0; t < nprod; ++t)
        prods.emplace_back(producer, opt, std::ref(bq), std::cref(cfg));

    for (std::thread &th : prods) th.join();
    for (std::thread &th : workers) th.join();

    st.reads = bq.reads.load();
    st.bases = bq.bases.load();
    merge_records(bq.results, loci, cap_bases, st);
    st.capped = st.capped || bq.soft_capped.load();

    if (st.capped) {
        double capv = (cap_total != std::numeric_limits<size_t>::max()
                           ? (double)cap_total : 0.0);
        if (capv >= 1e9)
            fprintf(stderr,
                    "[M::poa_polish] candidate read cap reached (%.2f Gbp); "
                    "raise --poa-max-candidate-bases to scan the remaining reads\n",
                    capv / 1e9);
        else
            fprintf(stderr,
                    "[M::poa_polish] candidate read cap reached (%.2f Mbp); "
                    "raise --poa-max-candidate-bases to scan the remaining reads\n",
                    capv / 1e6);
    }
    return st;
}

/* ------------------------------------------------------------------------- */
/* k-mer support                                                              */
/* ------------------------------------------------------------------------- */

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
/* one locus                                                                  */
/* ------------------------------------------------------------------------- */

struct PatchResult {
    bool patch = false;
    double seconds = 0.0;
    std::string repeat;
    LocusStat stat;
};

/**
 * Try to polish one repeat using the reads already recruited for it.
 * Pure with respect to `seq`: returns the replacement repeat and a status, so
 * many loci can be computed in parallel and spliced afterwards.
 */
PatchResult compute_locus(const std::string &seq,
                          const TrRegion &tr,
                          const Locus &loc,
                          const hifiasm_opt_t *opt)
{
    PatchResult out;
    LocusStat &stat = out.stat;
    stat.contig = loc.contig;
    stat.start = tr.start;
    stat.end = tr.end;
    stat.period = tr.period;
    stat.old_len = stat.new_len = tr.end - tr.start;
    stat.cont_old = stat.cont_new = -1.0;
    const int flank = loc.flank;
    const long long region_len = tr.end - tr.start;
    stat.n_reads = (int)loc.cut_reads.size();

    if (opt->poa_max_locus_len > 0 && region_len > opt->poa_max_locus_len) {
        stat.status = "skip:array_too_long";
        return out;
    }

    stat.span_med = region_len;
    const int min_reads = std::max(3, opt->poa_min_reads);
    if ((int)loc.cut_reads.size() < min_reads) {
        stat.status = "skip:insufficient_reads";
        return out;
    }

    /* Median read-array length + tolerance, used to avoid chasing consensus
     * noise on an already-correct array. */
    long long span_med = region_len, span_tol = 40;
    if (!loc.spans.empty()) {
        std::vector<long long> sorted = loc.spans;
        std::sort(sorted.begin(), sorted.end());
        span_med = sorted[sorted.size() / 2];
        stat.span_min = sorted.front();
        stat.span_max = sorted.back();
        std::vector<long long> dev;
        dev.reserve(sorted.size());
        for (long long x : sorted) dev.push_back(std::llabs(x - span_med));
        std::sort(dev.begin(), dev.end());
        long long mad = dev[dev.size() / 2];
        span_tol = std::max<long long>(40, std::max(3 * mad, (long long)(0.02 * span_med)));
    }

    /* The splice guard below refuses to touch a locus whose length already
     * agrees with the spanning reads. Check that *before* running POA: on a
     * correct assembly this skips essentially every locus, and POA on long
     * tandem arrays is the expensive part. */
    stat.span_med = span_med;
    if (std::llabs(region_len - span_med) <= span_tol) {
        stat.status = "unchanged_read_length";
        return out;
    }

    std::unordered_set<uint64_t> read_kmers;
    for (const std::string &r : loc.cut_reads) {
        std::unordered_set<uint64_t> s = kmer_set(r, HIT_KMER);
        read_kmers.insert(s.begin(), s.end());
    }
    const std::string old_repeat = seq.substr(tr.start, region_len);
    stat.cont_old = containment(old_repeat, read_kmers, HIT_KMER);

    /* Bound the DP work for very long arrays: pass at most ~50 kb of reads
     * per locus (evenly sampled), with at least --poa-min-reads. */
    const size_t POA_BUDGET = 65536;
    size_t n_rec = loc.cut_reads.size();
    long long total_bases = 0;
    for (const std::string &r : loc.cut_reads) total_bases += (long long)r.size();
    size_t n_use = n_rec;
    if (total_bases > (long long)POA_BUDGET) {
        n_use = (size_t)((double)n_rec * (double)POA_BUDGET / (double)total_bases);
        if (n_use < (size_t)min_reads) n_use = (size_t)min_reads;
        if (n_use > n_rec) n_use = n_rec;
    }
    std::vector<const char *> ptrs;
    std::vector<size_t> lens;
    ptrs.reserve(n_use);
    lens.reserve(n_use);
    for (size_t j = 0; j < n_use; ++j) {
        size_t idx = (size_t)((double)j * (double)n_rec / (double)n_use);
        if (idx >= n_rec) idx = n_rec - 1;
        ptrs.push_back(loc.cut_reads[idx].data());
        lens.push_back(loc.cut_reads[idx].size());
    }

    const char *dbg_dir = getenv("POA_POLISH_DEBUG_DIR");
    if (dbg_dir && *dbg_dir) {
        char fn[1024];
        snprintf(fn, sizeof(fn), "%s/locus_%s_%lld_%lld.reads.fa", dbg_dir,
                 stat.contig.c_str(), tr.start, tr.end);
        FILE *fp = fopen(fn, "w");
        if (fp) {
            for (size_t i = 0; i < loc.cut_reads.size(); ++i)
                fprintf(fp, ">r%zu\n%s\n", i, loc.cut_reads[i].c_str());
            fclose(fp);
        }
    }

    PoaOpts o;
    o.band_width = 0;                       /* crate default (50 + adaptive) */
    o.min_reads = (size_t)min_reads;
    o.no_homopolymer = 0;                   /* ONT: keep HPC + tandem gaps */
    o.max_read_len = 0;

    PoaResult *res = poa_consensus_flanked_c(ptrs.data(), lens.data(), ptrs.size(),
                                             loc.left.data(), loc.left.size(),
                                             loc.right.data(), loc.right.size(), &o);
    if (!res || res->status != POA_OK || !res->repeat) {
        stat.status = (res && res->error) ? (std::string("skip:") + res->error)
                                          : std::string("skip:poa_error");
        if (res) poa_result_free(res);
        return out;
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
            size_t pl = whole.find(loc.left.substr(flank - probe, probe));
            size_t pr = whole.find(loc.right.substr(0, probe));
            if (pl != std::string::npos && pr != std::string::npos && pr > pl) {
                anchored = true;
                break;
            }
        }
    }
    poa_result_free(res);

    if (!anchored) {
        stat.status = "skip:flank_not_anchored";
        return out;
    }

    double ratio = region_len > 0 ? (double)new_repeat.size() / (double)region_len : 0.0;
    if (new_repeat.empty() || ratio < 0.4 || ratio > 2.5) {
        stat.status = "skip:length_outlier";
        return out;
    }

    stat.cont_new = containment(new_repeat, read_kmers, HIT_KMER);
    if (stat.cont_new >= 0 && stat.cont_old >= 0 && stat.cont_new + 0.05 < stat.cont_old) {
        stat.status = "skip:lower_support";
        return out;
    }

    const long long delta = (long long)new_repeat.size() - region_len;
    if (new_repeat == old_repeat ||
        (std::llabs(delta) <= 2 && stat.cont_new <= stat.cont_old + 0.005)) {
        stat.status = "unchanged";
        return out;
    }

    long long old_err = std::llabs(region_len - span_med);
    long long new_err = std::llabs((long long)new_repeat.size() - span_med);
    if (old_err <= span_tol || new_err >= old_err) {
        stat.status = "unchanged_read_length";
        return out;
    }

    out.patch = true;
    out.repeat = new_repeat;
    stat.old_len = region_len;
    stat.new_len = (long long)new_repeat.size();
    stat.status = "patched";
    return out;
}

/* ------------------------------------------------------------------------- */
/* GFA: pass 1 (detect loci), read streaming, pass 2 (patch/write)            */
/* ------------------------------------------------------------------------- */

int locus_flank(long long L, const TrRegion &tr, const hifiasm_opt_t *opt)
{
    int flank = opt->poa_flank_len > 0 ? opt->poa_flank_len : 120;
    if (flank < 40) flank = 40;
    flank = (int)std::min<long long>(flank, std::min<long long>(tr.start, L - tr.end));
    return flank;
}

/** Parse the first three fields of an S line and keep the tag remainder. */
bool parse_s_line(const std::string &line, std::string &name, std::string &seq,
                  std::string &tags)
{
    size_t t1 = line.find('\t');
    size_t t2 = (t1 == std::string::npos) ? std::string::npos : line.find('\t', t1 + 1);
    if (t2 == std::string::npos) return false;
    size_t t3 = line.find('\t', t2 + 1);
    name = line.substr(t1 + 1, t2 - t1 - 1);
    if (t3 == std::string::npos) {
        seq = line.substr(t2 + 1);
        tags.clear();
    } else {
        seq = line.substr(t2 + 1, t3 - t2 - 1);
        tags = line.substr(t3 + 1);
    }
    return true;
}

/** Pass 1: collect repeat loci (contig coordinates + flanks) from a GFA. */
std::vector<Locus> collect_loci(const std::string &gfa_name, const hifiasm_opt_t *opt)
{
    std::vector<Locus> loci;
    FILE *in = fopen(gfa_name.c_str(), "r");
    if (!in) return loci;

    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    std::string name, seq, tags;
    while ((len = getline(&line, &cap, in)) > 0) {
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len < 3 || line[0] != 'S' || line[1] != '\t') continue;
        if (!parse_s_line(std::string(line, len), name, seq, tags)) continue;
        if (seq == "*" || seq.size() < (size_t)opt->poa_min_array) continue;

        std::vector<TrRegion> trs = detect_trs(seq, opt->poa_max_period,
                                              opt->poa_min_array);
        for (const TrRegion &tr : trs) {
            if (opt->poa_max_locus_len > 0 &&
                tr.end - tr.start > opt->poa_max_locus_len) continue;
            int flank = locus_flank((long long)seq.size(), tr, opt);
            if (flank < 40) continue;
            Locus loc;
            loc.contig = name;
            loc.start = tr.start;
            loc.end = tr.end;
            loc.period = tr.period;
            loc.flank = flank;
            loc.left = seq.substr(tr.start - flank, flank);
            loc.right = seq.substr(tr.end, flank);
            loci.push_back(std::move(loc));
        }
    }
    free(line);
    fclose(in);
    return loci;
}

std::string make_out_name(const std::string &gfa, const char *ext)
{
    return strip_gfa_suffix(gfa) + ext;
}

bool polish_one_gfa(const std::string &gfa_name,
                    const hifiasm_opt_t *opt,
                    FILE *tsv,
                    int &n_loci,
                    int &n_patched,
                    long long &bases_added)
{
    const double t_start = now_s();

    /* ---- pass 1: find loci and their flanks ---------------------------- */
    std::vector<Locus> loci = collect_loci(gfa_name, opt);
    const double t_detect = now_s();
    if (loci.empty()) {
        fprintf(stderr, "[M::poa_polish] %s: no tandem repeats found\n", gfa_name.c_str());
    }

    /* ---- read recruitment, in index chunks ------------------------------ */
    RecruitStats st{0, 0, 0, 0, false};
    double t_index_total = 0.0, t_stream_total = 0.0;
    size_t nchunks_done = 0;
    double t_index_done = t_detect, t_stream = t_detect;
    if (!loci.empty()) {
        const size_t seeds_per_locus =
            8 * (size_t)std::max(1, g_edge_len - RECRUIT_K + 1);
        size_t chunk = opt->poa_index_chunk > 0
                           ? (size_t)opt->poa_index_chunk
                           : loci.size();
        if (chunk == 0) chunk = loci.size();
        size_t nchunks = (loci.size() + chunk - 1) / chunk;
        double full_gb = (double)loci.size() * (double)seeds_per_locus * 72.0 / 1e9;
        double chunk_gb = (double)std::min(chunk, loci.size()) *
                          (double)seeds_per_locus * 72.0 / 1e9;
        if (loci.size() > 500000 || full_gb > 8.0) {
            fprintf(stderr,
                    "[M::poa_polish] warning: %zu loci; full seed index ~%.1f GB "
                    "(edge %d bp); %zu chunk(s) of <=%zu loci -> ~%.1f GB per chunk, "
                    "%zu read pass(es)\n",
                    loci.size(), full_gb, g_edge_len, nchunks,
                    std::min(chunk, loci.size()), chunk_gb, nchunks);
            if (nchunks == 1)
                fprintf(stderr,
                        "[M::poa_polish] consider --poa-min-array 150/200 or "
                        "--poa-index-chunk 250000 to bound memory\n");
        }

        size_t cap_total = std::numeric_limits<size_t>::max();
        if (opt->poa_max_cand_bases > 0)
            cap_total = (size_t)opt->poa_max_cand_bases;
        size_t stored = 0;

        for (size_t c = 0; c < nchunks; ++c) {
            size_t begin = c * chunk;
            size_t end = std::min(begin + chunk, loci.size());
            SeedIndex idx;
            double ti = now_s();
            build_index(loci, idx, begin, end);
            t_index_total += now_s() - ti;

            fprintf(stderr,
                    "[M::poa_polish] %s: chunk %zu/%zu: loci %zu-%zu (%zu loci, "
                    "%zu seed k-mers); streaming reads...\n",
                    gfa_name.c_str(), c + 1, nchunks, begin, end, end - begin,
                    idx.head.size());
            double ts = now_s();
            RecruitStats part = stream_recruit(opt, loci, idx, cap_total, stored);
            double dt = now_s() - ts;
            t_stream_total += dt;
            stored += (size_t)part.stored_bases;
            st = part;
            ++nchunks_done;

            size_t kept = 0;
            for (size_t i = begin; i < end; ++i) kept += loci[i].cut_reads.size();
            fprintf(stderr,
                    "[M::poa_polish] %s: chunk %zu/%zu: %lld reads, kept %zu "
                    "(%.2f Mbp), %.2fs (%.1f Mbp/s)\n",
                    gfa_name.c_str(), c + 1, nchunks, part.reads, kept,
                    part.stored_bases / 1e6, dt,
                    dt > 0.001 ? part.bases / 1e6 / dt : 0.0);
            if (part.capped) {
                fprintf(stderr,
                        "[M::poa_polish] candidate cap reached; skipping remaining "
                        "%zu chunk(s)\n", nchunks - c - 1);
                break;
            }
        }

        t_index_done = t_detect + t_index_total;
        t_stream = t_index_done + t_stream_total;
        size_t kept_all = 0;
        for (const Locus &l : loci) kept_all += l.cut_reads.size();
        fprintf(stderr,
                "[M::poa_polish] %s: %zu loci, %zu chunk(s), kept %zu reads "
                "(%.2f Mbp); detect %.2fs, index %.2fs, stream %.2fs "
                "(%.1f Mbp/s aggregate)\n",
                gfa_name.c_str(), loci.size(), nchunks_done, kept_all,
                (double)stored / 1e6, t_detect - t_start, t_index_total,
                t_stream_total,
                t_stream_total > 0.001 && nchunks_done > 0
                    ? (double)st.bases * (double)nchunks_done / 1e6 / t_stream_total
                    : 0.0);
    }

    /* ---- pass 2: patch and write --------------------------------------- */
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
    std::string name, seq, tags;
    std::vector<std::string> a_lines;
    std::vector<std::string> tk;
    size_t next_locus = 0;
    int n_seg = 0, n_seg_patched = 0;
    double t_poa = 0.0;
    bool have_seg = false;
    std::vector<TrRegion> cur_trs;
    std::vector<int> cur_ids;

    auto flush_segment = [&]() {
        if (!have_seg) return;
        ++n_seg;
        bool patched = false;
        std::vector<std::pair<long long, long long>> shifts; // (region end, delta)

        /* POA calls are independent per locus: run them in parallel, then
         * splice right-to-left so earlier coordinates stay valid. */
        const int L = (int)cur_trs.size();
        std::vector<PatchResult> results(L);
        std::vector<char> computed(L, 0);
        int nthreads = opt->thread_num > 0 ? opt->thread_num : 1;
        const double p0 = now_s();
        if (nthreads > 1 && L > 1) {
            std::atomic<size_t> next{0};
            auto worker = [&]() {
                for (;;) {
                    size_t k = next.fetch_add(1);
                    if (k >= (size_t)L) break;
                    if (cur_ids[k] < 0) continue;
                    const double q0 = now_s();
                    results[k] = compute_locus(seq, cur_trs[k], loci[cur_ids[k]], opt);
                    results[k].seconds = now_s() - q0;
                    computed[k] = 1;
                }
            };
            std::vector<std::thread> pool;
            pool.reserve(nthreads - 1);
            for (int t = 1; t < nthreads; ++t) pool.emplace_back(worker);
            worker();
            for (std::thread &th : pool) th.join();
        } else {
            for (int k = 0; k < L; ++k) {
                if (cur_ids[k] < 0) continue;
                const double q0 = now_s();
                results[k] = compute_locus(seq, cur_trs[k], loci[cur_ids[k]], opt);
                results[k].seconds = now_s() - q0;
                computed[k] = 1;
            }
        }
        t_poa += now_s() - p0;

        for (int k = L - 1; k >= 0; --k) {
            if (!computed[k]) continue;
            LocusStat &stt = results[k].stat;
            ++n_loci;
            if (results[k].patch) {
                patched = true;
                ++n_patched;
                bases_added += stt.new_len - stt.old_len;
                shifts.push_back({cur_trs[k].end, stt.new_len - stt.old_len});
                seq = seq.substr(0, cur_trs[k].start) + results[k].repeat +
                      seq.substr(cur_trs[k].end);
            }
            fprintf(stderr,
                    "[M::poa_polish] %s %lld-%lld period=%d copies=%.1f reads=%d "
                    "support=%.3f->%.3f med=%lld [%lld..%lld] %s (%.2fs)\n",
                    stt.contig.c_str(), stt.start, stt.end, stt.period,
                    (double)(cur_trs[k].end - cur_trs[k].start) / cur_trs[k].period,
                    stt.n_reads, stt.cont_old, stt.cont_new, stt.span_med,
                    stt.span_min, stt.span_max, stt.status.c_str(),
                    results[k].seconds);
            if (tsv)
                fprintf(tsv,
                        "%s\t%s\t%lld\t%lld\t%d\t%.2f\t%lld\t%lld\t%d\t%.4f\t%.4f\t%s\n",
                        gfa_name.c_str(), stt.contig.c_str(), stt.start, stt.end,
                        stt.period,
                        (double)(cur_trs[k].end - cur_trs[k].start) / cur_trs[k].period,
                        stt.old_len, stt.new_len, stt.n_reads,
                        stt.cont_old, stt.cont_new, stt.status.c_str());
        }
        if (patched) ++n_seg_patched;

        /* S line: update LN:i: but preserve every other tag */
        fprintf(out, "S\t%s\t%s", name.c_str(), seq.c_str());
        split_tabs(tags, tk);
        std::string rebuilt;
        bool had_ln = false;
        for (const std::string &t : tk) {
            if (t.rfind("LN:i:", 0) == 0) { had_ln = true; continue; }
            if (t.empty()) continue;
            if (!rebuilt.empty()) rebuilt += '\t';
            rebuilt += t;
        }
        if (!had_ln && seq != "*") {
            char buf[64];
            snprintf(buf, sizeof(buf), "LN:i:%lld", (long long)seq.size());
            if (!rebuilt.empty()) rebuilt += '\t';
            rebuilt += buf;
        }
        if (!rebuilt.empty()) fprintf(out, "\t%s", rebuilt.c_str());
        fprintf(out, "\n");

        /* A lines, with contig offsets shifted by upstream splices */
        for (const std::string &a : a_lines) {
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

        if (seq != "*") {
            fprintf(fa, ">%s\n", name.c_str());
            for (size_t i = 0; i < seq.size(); i += 60)
                fprintf(fa, "%.*s\n", (int)std::min<size_t>(60, seq.size() - i),
                        seq.c_str() + i);
        }

        a_lines.clear();
        cur_trs.clear();
        cur_ids.clear();
        have_seg = false;
    };

    while ((len = getline(&line, &cap, in)) > 0) {
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len == 0) continue;

        if (line[0] == 'S' && line[1] == '\t') {
            flush_segment();
            std::string s(line, len);
            if (!parse_s_line(s, name, seq, tags)) continue;
            have_seg = true;
            a_lines.clear();
            cur_trs.clear();
            cur_ids.clear();
            if (seq != "*") {
                cur_trs = detect_trs(seq, opt->poa_max_period,
                                     opt->poa_min_array);
                for (const TrRegion &tr : cur_trs) {
                    if (opt->poa_max_locus_len > 0 &&
                        tr.end - tr.start > opt->poa_max_locus_len) {
                        fprintf(stderr,
                                "[M::poa_polish] %s %lld-%lld array is %lld bp "
                                "(> --poa-max-locus-len %d); skipped\n",
                                name.c_str(), tr.start, tr.end,
                                tr.end - tr.start, opt->poa_max_locus_len);
                        cur_ids.push_back(-1);
                        continue;
                    }
                    int flank = locus_flank((long long)seq.size(), tr, opt);
                    if (flank < 40) { cur_ids.push_back(-1); continue; }
                    if (next_locus >= loci.size()) {
                        fprintf(stderr,
                                "[M::poa_polish] warning: locus mismatch in %s at %s:%lld\n",
                                gfa_name.c_str(), name.c_str(), tr.start);
                        cur_ids.push_back(-1);
                        continue;
                    }
                    const Locus &loc = loci[next_locus];
                    if (loc.contig != name || loc.start != tr.start || loc.end != tr.end) {
                        fprintf(stderr,
                                "[M::poa_polish] warning: locus mismatch in %s at %s:%lld\n",
                                gfa_name.c_str(), name.c_str(), tr.start);
                        cur_ids.push_back(-1);
                        ++next_locus; // keep the sequences aligned
                        continue;
                    }
                    cur_ids.push_back((int)next_locus);
                    ++next_locus;
                }
            }
        } else if (line[0] == 'A' && line[1] == '\t') {
            if (have_seg) a_lines.push_back(std::string(line, len));
        } else {
            flush_segment();
            fprintf(out, "%s\n", line);
        }
    }
    flush_segment();

    if (next_locus != loci.size()) {
        fprintf(stderr,
                "[M::poa_polish] warning: %s: used %zu of %zu loci\n",
                gfa_name.c_str(), next_locus, loci.size());
    }

    free(line);
    fclose(in);
    fclose(out);
    fclose(fa);

    const double t_end = now_s();
    fprintf(stderr, "[M::poa_polish] %s: %d contig(s), %d polished; timing: "
            "detect %.2fs, index %.2fs, stream %.2fs, patch %.2fs (poa %.2fs), total %.2fs\n",
            gfa_name.c_str(), n_seg, n_seg_patched, t_detect - t_start,
            t_index_done - t_detect, t_stream - t_index_done, t_end - t_stream,
            t_poa, t_end - t_start);
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

    if (opt->poa_seed_step >= 1 && opt->poa_seed_step <= 64) {
        g_seed_step = opt->poa_seed_step;
    } else {
        fprintf(stderr,
                "[M::poa_polish] warning: --poa-seed-step %d out of range [1,64]; using 2\n",
                opt->poa_seed_step);
        g_seed_step = 2;
    }
    if (opt->poa_edge_len >= 20 && opt->poa_edge_len <= 120) {
        g_edge_len = opt->poa_edge_len;
    } else {
        fprintf(stderr,
                "[M::poa_polish] warning: --poa-edge-len %d out of range [20,120]; using 60\n",
                opt->poa_edge_len);
        g_edge_len = 60;
    }
    if (opt->poa_max_cand_bases < 0) {
        fprintf(stderr,
                "[M::poa_polish] warning: --poa-max-candidate-bases %lld is negative; "
                "treating it as unlimited\n", opt->poa_max_cand_bases);
    }
    char capbuf[64];
    if (opt->poa_max_cand_bases > 0)
        snprintf(capbuf, sizeof(capbuf), "%.2f Gbp", opt->poa_max_cand_bases / 1e9);
    else
        snprintf(capbuf, sizeof(capbuf), "unlimited");
    fprintf(stderr,
            "[M::poa_polish] post-assembly tandem-repeat polishing (%s); streaming "
            "recruitment (flank %d bp, edge %d bp, seed step %d, candidate cap %s)\n",
            poa_version(), opt->poa_flank_len, g_edge_len, g_seed_step, capbuf);

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
        polish_one_gfa(fn, opt, tsv, n_loci, n_patched, bases_added);

    if (tsv) fclose(tsv);

    fprintf(stderr,
            "[M::poa_polish] done: %d locus/loci found, %d patched, %+lld bp; report: %s\n",
            n_loci, n_patched, bases_added, tsv_name.c_str());
    return 0;
}
