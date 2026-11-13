/**
 * polish_gfa -- run only hifiasm-fast's poa-consensus post-assembly pass on an
 * existing contig GFA. Useful for testing or re-polishing an assembly without
 * re-running hifiasm.
 *
 * The driver mirrors the hifiasm options (--poa-*) and also accepts shorter
 * aliases (--edge-len, --seed-step, ...). Run `polish_gfa --help` for details.
 */

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include "CommandLines.h"
#include "poa_polish.h"
#include "poa-ffi/poa_ffi.h"

static const char *kProg = "polish_gfa";

static void print_help(FILE *fh)
{
    fprintf(fh,
"polish_gfa -- run only hifiasm-fast's poa-consensus VNTR polishing pass\n"
"\n"
"Usage:\n"
"  %s [options] <prefix> <reads.fq[.gz]> [reads2.fq.gz ...]\n"
"  %s [options] -o <prefix> <reads.fq[.gz]> [...]\n"
"\n"
"<prefix> is the assembly prefix used for the hifiasm run (for example\n"
"'sample' for 'sample.bp.p_ctg.gfa'), not a GFA filename. The driver looks\n"
"for every contig graph hifiasm may have written with that prefix:\n"
"  <prefix>.p_ctg.gfa        <prefix>.a_ctg.gfa\n"
"  <prefix>.bp.p_ctg.gfa     <prefix>.bp.hap1.p_ctg.gfa\n"
"  <prefix>.bp.hap2.p_ctg.gfa\n"
"  <prefix>.hap1.p_ctg.gfa   <prefix>.hap2.p_ctg.gfa\n"
"  <prefix>.stg.p_ctg.gfa\n"
"\n"
"For each tandem-repeat locus it streams the reads once, rebuilds the array\n"
"with poa-consensus' single-allele flanked mode, and splices it back between\n"
"the same flank anchors when the reads support it. The input GFA is left\n"
"untouched; outputs are:\n"
"  <prefix>.poa.tsv              one row per locus, with status\n"
"  <gfa-prefix>.poa.gfa          patched contig graph\n"
"  <gfa-prefix>.poa.fa           patched contigs as FASTA\n"
"\n"
"Options:\n"
"  -o, --output-prefix PREFIX    assembly prefix (alternative to the first\n"
"                                positional argument)\n"
"  -t, --threads INT             POA worker threads, 1-64 [all available CPUs]\n"
"  -f, --flank INT               flank length used to anchor a repeat, >=40 [120]\n"
"      --min-reads INT           minimum spanning reads to rewrite a locus\n"
"                                (effective minimum 3) [6]\n"
"      --max-period INT          longest repeat unit considered [200]\n"
"      --edge-len INT            flank edge length used for seeding,\n"
"                                20-120 [60]\n"
"      --seed-step INT           check every Nth read k-mer when recruiting\n"
"                                reads, 1-64 [2]; 1 is most sensitive\n"
"      --max-candidate-bases INT cap on recruited read bases held in memory;\n"
"                                0 = unlimited [4000000000]\n"
"      --max-locus-len INT       skip arrays longer than this many bp;\n"
"                                0 = no limit [100000]\n"
"      --min-array INT           minimum array length considered [60];\n"
"                                use 150-200 for a WGS run (fewer seed loci)\n"
"      --index-chunk INT         build the seed index for this many loci at a\n"
"                                time; 0 = all at once [1000000]\n"
"      --decompressor STR        read stream backend: auto, zlib or rapidgzip\n"
"                                [auto; rapidgzip if found in PATH]\n"
"      --gzip-threads INT        rapidgzip threads per input file; 0 = auto\n"
"                                [auto: (threads / 2) / files]\n"
"  -h, --help                    show this help and exit\n"
"  -V, --version                 show version information and exit\n"
"\n"
"  The same options are also accepted with the hifiasm spelling:\n"
"  --poa-flank, --poa-min-reads, --poa-max-period, --poa-edge-len,\n"
"  --poa-seed-step, --poa-max-candidate-bases.\n"
"\n"
"Environment:\n"
"  POA_POLISH_DEBUG_DIR=DIR      dump the cut reads and consensus for every\n"
"                                locus (one .reads.fa and .cons.fa per locus)\n"
"\n"
"Examples:\n"
"  %s asm reads.fq.gz\n"
"  %s -o asm reads.fq.gz reads2.fq.gz\n"
"  %s asm reads.fq.gz --edge-len 120 --seed-step 1\n"
"  %s asm reads.fq.gz --min-reads 4 --max-candidate-bases 8000000000\n"
"\n"
"Defaults are the same as hifiasm's --poa-* options. Reads may be FASTQ or\n"
"FASTA, optionally gzip-compressed.\n"
"\n"
"Exit status: 0 success, 1 runtime failure, 2 usage error.\n"
"\n"
"To assemble and polish in one step, use:\n"
"  hifiasm --ont --poa-polish -o asm reads.fq.gz\n",
            kProg, kProg, kProg, kProg, kProg, kProg);
}

static void short_usage(FILE *fh)
{
    fprintf(fh,
            "Usage: %s [options] <prefix> <reads.fq[.gz]> [...]\n"
            "Try '%s --help' for the full option list.\n",
            kProg, kProg);
}

static int usage_error(const char *fmt, const char *arg = NULL)
{
    fprintf(stderr, "%s: ", kProg);
    fprintf(stderr, fmt, arg);
    fprintf(stderr, "\n");
    short_usage(stderr);
    return 2;
}

static bool parse_int(const char *s, long long &out)
{
    if (!s || !*s) return false;
    errno = 0;
    char *end = NULL;
    long long v = strtoll(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0') return false;
    out = v;
    return true;
}

int main(int argc, char **argv)
{
    const char *prefix = NULL;
    std::vector<char *> positional;

    hifiasm_opt_t opt;
    memset(&opt, 0, sizeof(opt));
    opt.is_poa_polish = 1;
    opt.poa_flank_len = 120;
    opt.poa_min_reads = 6;
    opt.poa_max_period = 200;
    opt.poa_max_cand_bases = 4000000000LL;
    opt.poa_seed_step = 2;
    opt.poa_edge_len = 60;
    opt.poa_max_locus_len = 100000;
    opt.poa_min_array = 60;
    opt.poa_index_chunk = 1000000;
    opt.poa_decompressor = 0;
    opt.poa_gzip_threads = 0;
    unsigned hw = std::thread::hardware_concurrency();
    opt.thread_num = hw ? (int)hw : 1;
    if (opt.thread_num > 64) opt.thread_num = 64;

    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];

        if (!strcmp(arg, "-h") || !strcmp(arg, "--help")) {
            print_help(stdout);
            return 0;
        }
        if (!strcmp(arg, "-V") || !strcmp(arg, "--version")) {
            printf("%s (hifiasm-fast poa-consensus driver)\n", kProg);
            printf("%s\n", poa_version());
            return 0;
        }
        if (!strcmp(arg, "--")) {
            for (++i; i < argc; ++i) positional.push_back(argv[i]);
            break;
        }

        if (arg[0] == '-' && arg[1] != '\0') {
            const char *val = NULL;
            const char *name = arg;

            auto take_value = [&](const char *what) -> const char * {
                if (i + 1 >= argc) {
                    fprintf(stderr, "%s: option '%s' requires a value\n", kProg, what);
                    short_usage(stderr);
                    exit(2);
                }
                return argv[++i];
            };

            if (!strcmp(name, "-o") || !strcmp(name, "--output-prefix") ||
                !strcmp(name, "--prefix")) {
                prefix = take_value(name);
                continue;
            }
            if (!strcmp(name, "-t") || !strcmp(name, "--threads")) {
                const char *tv = take_value(name);
                long long v = 0;
                if (!parse_int(tv, v) || v < 1 || v > 64) {
                    fprintf(stderr,
                            "%s: --threads must be an integer between 1 and 64 "
                            "(got '%s')\n", kProg, tv);
                    return 2;
                }
                opt.thread_num = (int)v;
                continue;
            }

            bool is_flank = !strcmp(name, "-f") || !strcmp(name, "--flank") ||
                            !strcmp(name, "--poa-flank");
            bool is_min_reads = !strcmp(name, "--min-reads") ||
                                !strcmp(name, "--poa-min-reads");
            bool is_max_period = !strcmp(name, "--max-period") ||
                                 !strcmp(name, "--poa-max-period");
            bool is_edge = !strcmp(name, "--edge-len") ||
                           !strcmp(name, "--poa-edge-len");
            bool is_step = !strcmp(name, "--seed-step") ||
                           !strcmp(name, "--poa-seed-step");
            bool is_cand = !strcmp(name, "--max-candidate-bases") ||
                           !strcmp(name, "--poa-max-candidate-bases");
            bool is_maxlocus = !strcmp(name, "--max-locus-len") ||
                               !strcmp(name, "--poa-max-locus-len");
            bool is_minarray = !strcmp(name, "--min-array") ||
                               !strcmp(name, "--poa-min-array");
            bool is_chunk = !strcmp(name, "--index-chunk") ||
                            !strcmp(name, "--poa-index-chunk");
            bool is_decomp = !strcmp(name, "--decompressor") ||
                             !strcmp(name, "--poa-decompressor");
            bool is_gzthreads = !strcmp(name, "--gzip-threads") ||
                                !strcmp(name, "--poa-gzip-threads");

            if (!is_flank && !is_min_reads && !is_max_period && !is_edge &&
                !is_step && !is_cand && !is_maxlocus && !is_decomp &&
                !is_gzthreads && !is_minarray && !is_chunk) {
                fprintf(stderr, "%s: unknown option '%s'\n", kProg, name);
                short_usage(stderr);
                return 2;
            }

            val = take_value(name);

            if (is_decomp) {
                if (!strcmp(val, "auto")) opt.poa_decompressor = 0;
                else if (!strcmp(val, "zlib")) opt.poa_decompressor = 1;
                else if (!strcmp(val, "rapidgzip")) opt.poa_decompressor = 2;
                else {
                    fprintf(stderr, "%s: --decompressor must be auto, zlib or "
                                    "rapidgzip (got '%s')\n", kProg, val);
                    return 2;
                }
                continue;
            }

            long long v = 0;
            if (!parse_int(val, v)) {
                fprintf(stderr, "%s: option '%s' expects an integer, got '%s'\n",
                        kProg, name, val);
                return 2;
            }

            if (is_flank) {
                if (v < 40) return usage_error("--flank must be >= 40 (got %s)",
                                               val);
                opt.poa_flank_len = (int)v;
            } else if (is_min_reads) {
                if (v < 1) return usage_error("--min-reads must be >= 1 (got %s)",
                                              val);
                opt.poa_min_reads = (int)v;
            } else if (is_max_period) {
                if (v < 2) return usage_error("--max-period must be >= 2 (got %s)",
                                              val);
                opt.poa_max_period = (int)v;
            } else if (is_edge) {
                if (v < 20 || v > 120)
                    return usage_error("--edge-len must be between 20 and 120 "
                                       "(got %s)", val);
                opt.poa_edge_len = (int)v;
            } else if (is_step) {
                if (v < 1 || v > 64)
                    return usage_error("--seed-step must be between 1 and 64 "
                                       "(got %s)", val);
                opt.poa_seed_step = (int)v;
            } else if (is_cand) {
                if (v < 0)
                    return usage_error("--max-candidate-bases must be >= 0 "
                                       "(got %s)", val);
                opt.poa_max_cand_bases = v;
            } else if (is_minarray) {
                if (v < 40)
                    return usage_error("--min-array must be >= 40 (got %s)", val);
                opt.poa_min_array = (int)v;
            } else if (is_chunk) {
                if (v < 0)
                    return usage_error("--index-chunk must be >= 0 (got %s)", val);
                opt.poa_index_chunk = (int)v;
            } else if (is_maxlocus) {
                if (v < 0)
                    return usage_error("--max-locus-len must be >= 0 (got %s)", val);
                opt.poa_max_locus_len = (int)v;
            } else {
                if (v < 0)
                    return usage_error("--gzip-threads must be >= 0 (got %s)", val);
                opt.poa_gzip_threads = (int)v;
            }
            continue;
        }

        positional.push_back(argv[i]);
    }

    if (!prefix) {
        if (positional.empty()) {
            fprintf(stderr, "%s: missing assembly prefix and read files\n", kProg);
            short_usage(stderr);
            return 2;
        }
        prefix = positional.front();
        positional.erase(positional.begin());
    }
    if (positional.empty()) {
        fprintf(stderr, "%s: missing read file(s)\n", kProg);
        short_usage(stderr);
        return 2;
    }
    if (!*prefix) {
        fprintf(stderr, "%s: assembly prefix must not be empty\n", kProg);
        return 2;
    }
    size_t plen = strlen(prefix);
    bool looks_like_gfa =
        (plen >= 4 && !strcmp(prefix + plen - 4, ".gfa")) ||
        (plen >= 7 && !strcmp(prefix + plen - 7, ".gfa.gz"));
    if (looks_like_gfa) {
        fprintf(stderr,
                "%s: '%s' looks like a GFA filename, but an assembly prefix is "
                "expected\n"
                "  e.g. pass 'asm' for 'asm.bp.p_ctg.gfa'\n",
                kProg, prefix);
        return 2;
    }

    /* Fail early with a readable message instead of producing empty outputs. */
    static const char *kGfaSuffixes[] = {
        ".p_ctg.gfa", ".a_ctg.gfa", ".bp.p_ctg.gfa",
        ".bp.hap1.p_ctg.gfa", ".bp.hap2.p_ctg.gfa",
        ".hap1.p_ctg.gfa", ".hap2.p_ctg.gfa", ".stg.p_ctg.gfa"
    };
    bool have_gfa = false;
    std::string expected;
    for (const char *suf : kGfaSuffixes) {
        std::string fn = std::string(prefix) + suf;
        FILE *fp = fopen(fn.c_str(), "r");
        if (fp) {
            have_gfa = true;
            fclose(fp);
            break;
        }
        if (expected.empty()) expected = std::string(prefix) + ".bp.p_ctg.gfa";
    }
    if (!have_gfa) {
        fprintf(stderr,
                "%s: no contig GFA found for prefix '%s'\n"
                "  expected something like %s (from `hifiasm -o %s ...`)\n",
                kProg, prefix, expected.c_str(), prefix);
        return 1;
    }

    int n_readable = 0;
    for (char *r : positional) {
        FILE *fp = fopen(r, "rb");
        if (fp) {
            ++n_readable;
            fclose(fp);
        } else {
            fprintf(stderr, "%s: warning: cannot open read file '%s'\n", kProg, r);
        }
    }
    if (n_readable == 0) {
        fprintf(stderr, "%s: no readable read files\n", kProg);
        return 1;
    }

    opt.output_file_name = (char *)prefix;
    opt.num_reads = (int)positional.size();
    opt.read_file_names = positional.data();

#ifdef POA_POLISH_DRIVER_VERBOSE
    fprintf(stderr, "%s: prefix=%s reads=%d\n", kProg, prefix, opt.num_reads);
#endif

    return poa_polish_run(&opt);
}
