#ifndef POA_FFI_H
#define POA_FFI_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Status codes mirrored from Rust. */
#define POA_OK          0
#define POA_ERR_DEPTH   1
#define POA_ERR_FLANK   2
#define POA_ERR_OTHER   3

typedef struct PoaOpts {
    size_t band_width;      /* 0 => crate default (50) */
    size_t min_reads;       /* 0 => crate default (3) */
    int    no_homopolymer;  /* 1 => disable homopolymer/tandem candidates */
    size_t max_read_len;    /* 0 => no limit */
} PoaOpts;

typedef struct PoaResult {
    char  *repeat;      /* consensus between flanks, NUL-terminated */
    size_t repeat_len;
    char  *whole;       /* full consensus including flanks, NUL-terminated */
    size_t whole_len;
    size_t n_reads;
    int    status;      /* POA_OK / POA_ERR_* */
    char  *error;       /* NUL-terminated error text, or NULL */
} PoaResult;

PoaResult *poa_consensus_flanked_c(const char *const *read_ptrs,
                                   const size_t *read_lens,
                                   size_t n_reads,
                                   const char *left, size_t left_len,
                                   const char *right, size_t right_len,
                                   const PoaOpts *opts);

void poa_result_free(PoaResult *r);
const char *poa_version(void);
const char *poa_result_error(const PoaResult *r);
int poa_is_dna(const char *seq, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* POA_FFI_H */
