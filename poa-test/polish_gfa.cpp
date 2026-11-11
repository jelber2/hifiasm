/**
 * polish_gfa -- run only hifiasm-fast's poa-consensus post-assembly pass on an
 * existing contig GFA. Useful for testing the polish stage independently.
 *
 * Usage: ./polish_gfa <prefix> <reads.fq[.gz]> [...]
 *
 * Reads <prefix>.p_ctg.gfa / .bp.p_ctg.gfa / ... and writes
 * <prefix>....poa.gfa / .poa.fa / <prefix>.poa.tsv.
 */

#include <cstdio>
#include <cstring>

#include "CommandLines.h"
#include "poa_polish.h"

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <prefix> <reads.fq[.gz]> [...]\n", argv[0]);
        return 1;
    }

    hifiasm_opt_t opt;
    memset(&opt, 0, sizeof(opt));
    opt.is_poa_polish = 1;
    opt.poa_flank_len = 120;
    opt.poa_min_reads = 6;
    opt.poa_max_period = 200;
    opt.thread_num = 1;
    opt.output_file_name = argv[1];
    opt.num_reads = argc - 2;
    opt.read_file_names = argv + 2;

    return poa_polish_run(&opt);
}
