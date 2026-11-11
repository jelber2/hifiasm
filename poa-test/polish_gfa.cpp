/**
 * polish_gfa -- run only hifiasm-fast's poa-consensus post-assembly pass on an
 * existing contig GFA. Useful for testing the polish stage independently.
 *
 * Usage: ./polish_gfa <prefix> <reads.fq[.gz]> [...]
 *          [--edge-len N] [--seed-step N] [--max-candidate-bases N]
 *          [--min-reads N] [--flank N]
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>

#include "CommandLines.h"
#include "poa_polish.h"

int main(int argc, char **argv)
{
    const char *prefix = NULL;
    std::vector<char *> reads;

    hifiasm_opt_t opt;
    memset(&opt, 0, sizeof(opt));
    opt.is_poa_polish = 1;
    opt.poa_flank_len = 120;
    opt.poa_min_reads = 6;
    opt.poa_max_period = 200;
    opt.poa_max_cand_bases = 4000000000LL;
    opt.poa_seed_step = 2;
    opt.poa_edge_len = 60;
    opt.thread_num = 1;

    for (int i = 1; i < argc; ++i) {
        auto need = [&](const char *name) -> const char * {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s requires a value\n", name);
                exit(1);
            }
            return argv[++i];
        };
        if (!strcmp(argv[i], "--edge-len")) opt.poa_edge_len = atol(need("--edge-len"));
        else if (!strcmp(argv[i], "--seed-step")) opt.poa_seed_step = atol(need("--seed-step"));
        else if (!strcmp(argv[i], "--max-candidate-bases"))
            opt.poa_max_cand_bases = atoll(need("--max-candidate-bases"));
        else if (!strcmp(argv[i], "--min-reads")) opt.poa_min_reads = atol(need("--min-reads"));
        else if (!strcmp(argv[i], "--flank")) opt.poa_flank_len = atol(need("--flank"));
        else if (argv[i][0] == '-') {
            fprintf(stderr, "unknown option %s\n", argv[i]);
            return 1;
        } else if (!prefix) {
            prefix = argv[i];
        } else {
            reads.push_back(argv[i]);
        }
    }

    if (!prefix || reads.empty()) {
        fprintf(stderr,
                "Usage: %s <prefix> <reads.fq[.gz]> [...] "
                "[--edge-len N] [--seed-step N] [--max-candidate-bases N] "
                "[--min-reads N] [--flank N]\n",
                argv[0]);
        return 1;
    }

    opt.output_file_name = (char *)prefix;
    opt.num_reads = (int)reads.size();
    opt.read_file_names = reads.data();

    return poa_polish_run(&opt);
}
