//! Minimal C ABI around `poa_consensus::consensus_flanked`.
//!
//! hifiasm-fast calls this once per tandem-repeat locus, passing the reads that
//! were assembled into the contig, cut down to the repeat plus both flanks.
//! Both the flank-sliced repeat and the whole (flank-included) consensus are
//! returned so the caller can verify the flanks were actually anchored before
//! splicing the result back into the contig.

use std::ffi::{c_char, CStr, CString};
use std::slice;
use std::ptr;

use poa_consensus::config::PoaConfig;
use poa_consensus::{consensus_flanked, PoaError};

/// Status codes returned in `poa_result_t::status`.
pub const POA_OK: i32 = 0;
pub const POA_ERR_DEPTH: i32 = 1;
pub const POA_ERR_FLANK: i32 = 2;
pub const POA_ERR_OTHER: i32 = 3;

#[repr(C)]
pub struct PoaOpts {
    /// Band half-width; 0 uses the crate default (50).
    pub band_width: usize,
    /// Minimum reads; 0 uses the crate default (3).
    pub min_reads: usize,
    /// 1 disables homopolymer compression and homopolymer/tandem gap costs.
    pub no_homopolymer: i32,
    /// Reads longer than this are skipped; 0 means no limit.
    pub max_read_len: usize,
}

#[repr(C)]
pub struct PoaResult {
    /// Consensus BETWEEN the flanks (flanks excluded), NUL-terminated.
    pub repeat: *mut c_char,
    pub repeat_len: usize,
    /// Full consensus, flanks included, NUL-terminated.
    pub whole: *mut c_char,
    pub whole_len: usize,
    pub n_reads: usize,
    /// POA_OK on success, POA_ERR_* otherwise.
    pub status: i32,
    /// NUL-terminated error text, or NULL.
    pub error: *mut c_char,
}

fn cstr_opt(p: *const c_char, len: usize) -> Option<&'static [u8]> {
    if p.is_null() {
        return None;
    }
    // SAFETY: the caller guarantees the pointer is valid for `len` bytes.
    Some(unsafe { slice::from_raw_parts(p as *const u8, len) })
}

unsafe fn into_c_string(bytes: &[u8]) -> (*mut c_char, usize) {
    match CString::new(bytes.to_vec()) {
        Ok(c) => {
            let len = c.as_bytes().len();
            (c.into_raw(), len)
        }
        Err(_) => {
            let empty = CString::new("").unwrap();
            (empty.into_raw(), 0)
        }
    }
}

unsafe fn error_c_string(msg: &str) -> *mut c_char {
    CString::new(msg.replace('\0', " "))
        .map(|c| c.into_raw())
        .unwrap_or(ptr::null_mut())
}

/// Build a flanked single-allele consensus of the repeat between `left` and
/// `right`, from `n_reads` reads given as (pointer, length) pairs.
///
/// Returns a heap-allocated `PoaResult` which must be released with
/// `poa_result_free`.
#[no_mangle]
pub unsafe extern "C" fn poa_consensus_flanked_c(
    read_ptrs: *const *const c_char,
    read_lens: *const usize,
    n_reads: usize,
    left: *const c_char,
    left_len: usize,
    right: *const c_char,
    right_len: usize,
    opts: *const PoaOpts,
) -> *mut PoaResult {
    let mut out = Box::new(PoaResult {
        repeat: ptr::null_mut(),
        repeat_len: 0,
        whole: ptr::null_mut(),
        whole_len: 0,
        n_reads: 0,
        status: POA_ERR_OTHER,
        error: ptr::null_mut(),
    });

    let Some(left_b) = cstr_opt(left, left_len) else {
        out.error = error_c_string("null left flank");
        return Box::into_raw(out);
    };
    let Some(right_b) = cstr_opt(right, right_len) else {
        out.error = error_c_string("null right flank");
        return Box::into_raw(out);
    };

    let mut config = PoaConfig::default();
    let mut max_read_len = 0usize;
    if !opts.is_null() {
        let o = &*opts;
        if o.band_width > 0 {
            config.band_width = o.band_width;
        }
        if o.min_reads > 0 {
            config.min_reads = o.min_reads;
        }
        if o.no_homopolymer != 0 {
            config.homopolymer_min_run = 0;
            config.homopolymer_gap = None;
            config.tandem_gap = None;
        }
        max_read_len = o.max_read_len;
    }

    if read_ptrs.is_null() || n_reads == 0 {
        out.error = error_c_string("no reads");
        return Box::into_raw(out);
    }

    let ptrs = slice::from_raw_parts(read_ptrs, n_reads);
    let lens = slice::from_raw_parts(read_lens, n_reads);
    let mut reads: Vec<&[u8]> = Vec::with_capacity(n_reads);
    for i in 0..n_reads {
        if ptrs[i].is_null() || lens[i] == 0 {
            continue;
        }
        if max_read_len > 0 && lens[i] > max_read_len {
            continue;
        }
        reads.push(slice::from_raw_parts(ptrs[i] as *const u8, lens[i]));
    }
    if reads.is_empty() {
        out.error = error_c_string("all reads filtered out");
        return Box::into_raw(out);
    }
    out.n_reads = reads.len();

    // Build the whole-region consensus once, then slice it to the repeat in C++
    // (the crate exposes both the sliced and the unsliced forms through the
    // same call, so call the internal entry point via consensus_flanked and
    // recover the whole sequence when the flanks were not found).
    //
    // `consensus_flanked` returns the repeat only; to get the whole sequence we
    // additionally run the plain consensus on the same (flank-carrying) reads.
    // The whole consensus is needed as a flank-anchoring sanity check.
    match consensus_flanked(&reads, left_b, right_b, &config) {
        Ok(c) => {
            let (rp, rl) = into_c_string(&c.sequence);
            out.repeat = rp;
            out.repeat_len = rl;
            out.status = POA_OK;
            // Plain consensus of the same reads = full region including flanks.
            match poa_consensus::consensus(&reads, 0, &config) {
                Ok(w) => {
                    let (wp, wl) = into_c_string(&w.sequence);
                    out.whole = wp;
                    out.whole_len = wl;
                }
                Err(_) => { /* sanity check unavailable; caller treats whole as absent */ }
            }
            Box::into_raw(out)
        }
        Err(e) => {
            out.status = match e {
                PoaError::InsufficientDepth { .. } => POA_ERR_DEPTH,
                _ => POA_ERR_OTHER,
            };
            out.error = error_c_string(&e.to_string());
            Box::into_raw(out)
        }
    }
}

/// Release a result previously returned by `poa_consensus_flanked_c`.
#[no_mangle]
pub unsafe extern "C" fn poa_result_free(r: *mut PoaResult) {
    if r.is_null() {
        return;
    }
    let boxed = Box::from_raw(r);
    if !boxed.repeat.is_null() {
        drop(CString::from_raw(boxed.repeat));
    }
    if !boxed.whole.is_null() {
        drop(CString::from_raw(boxed.whole));
    }
    if !boxed.error.is_null() {
        drop(CString::from_raw(boxed.error));
    }
}

/// Version string of the linked poa-consensus crate (for logging).
#[no_mangle]
pub extern "C" fn poa_version() -> *const c_char {
    static VERSION: &[u8] = b"poa-consensus staticlib via poa-ffi 0.1.0\0";
    VERSION.as_ptr() as *const c_char
}

/// Unused helper kept for symmetry with the crate API; returns 1 if `seq` is
/// valid ASCII DNA (A/C/G/T/N, any case).
#[no_mangle]
pub unsafe extern "C" fn poa_is_dna(seq: *const c_char, len: usize) -> i32 {
    if seq.is_null() {
        return 0;
    }
    let s = slice::from_raw_parts(seq as *const u8, len);
    if s.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n')) {
        1
    } else {
        0
    }
}

/// Expose the NUL-terminated error text of a result (NULL when none).
#[no_mangle]
pub unsafe extern "C" fn poa_result_error(r: *const PoaResult) -> *const c_char {
    if r.is_null() {
        return ptr::null();
    }
    let e = (*r).error;
    if e.is_null() {
        return ptr::null();
    }
    CStr::from_ptr(e).as_ptr()
}
