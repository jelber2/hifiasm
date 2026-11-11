#ifndef __POA_POLISH__
#define __POA_POLISH__

#include "CommandLines.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Post-assembly tandem-repeat polishing.
 *
 * Runs after the contig graph has been written. For every tandem repeat longer
 * than a couple of copies, the reads that span the locus plus both flanks are
 * recovered by streaming the input once against a seed index built from the
 * locus flanks (the read set itself is never held in memory), cut down to
 * [repeat + both flanks], and fed to poa-consensus' single-allele flanked mode.
 * The resulting repeat is spliced back between the same flank anchors when it
 * is better supported by the reads.
 *
 * Writes <prefix>.poa.gfa / .poa.fa next to the assembly and a .poa.tsv report.
 * Returns 0 on success, non-zero if no assembly could be polished.
 */
int poa_polish_run(const hifiasm_opt_t *opt);

#ifdef __cplusplus
}
#endif

#endif
