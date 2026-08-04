//! One rule, shared by every dev tool: should this run prove the reference FASTA still
//! matches its `.fai`?
//!
//! **Why this exists.** The check reads the whole FASTA — a *fixed* cost, about 11 seconds
//! of CPU on GRCh38, the same whether the tool then walks one contig or all of them. On a
//! real run that rounds to nothing. On a short one it is most of the run: the generic walk
//! over chr21 takes under two seconds, so an unconditional check makes the tool spend ~85 %
//! of its life proving something the operator already knows. That is the cost this knob
//! exists to remove, and **only** in the tools — the library default stays
//! [`ReferenceCheck::VerifyAgainstIndex`].
//!
//! **The default here is still to check.** These tools produce the dumps used as the
//! byte-identity gate for performance work, and a reference that silently disagreed with its
//! index would corrupt that gate without raising an error anywhere. So skipping is opt-in,
//! per run, by the person who knows they are re-running a fixture they trust.
//!
//! ```text
//! PVC_TRUST_REFERENCE_INDEX=1 ./target/release/examples/ng_generic_walk_probe …
//! ```
//!
//! The value is read the way the other knobs in these tools are read — `1` on, `0` off,
//! anything else is a usage error rather than a silent default, because a misspelled value
//! that quietly meant "check" would show up only as an unexplained 11 seconds.

use pop_var_caller::ng::reference_info::ReferenceCheck;

/// The environment variable that turns the check off. Named for what it grants, not for what
/// it saves: the operator is asserting the index is trustworthy, and the time is the
/// consequence.
pub const TRUST_REFERENCE_INDEX_VAR: &str = "PVC_TRUST_REFERENCE_INDEX";

/// Read the knob. Absent → check (the safe default). `1` → skip. `0` → check, said out loud.
///
/// # Errors
///
/// Any other value, so a typo is a usage error rather than a silent fallback.
pub fn reference_check_from_env() -> Result<ReferenceCheck, String> {
    match std::env::var(TRUST_REFERENCE_INDEX_VAR) {
        Err(_) => Ok(ReferenceCheck::VerifyAgainstIndex),
        Ok(value) => match value.as_str() {
            "1" => Ok(ReferenceCheck::TrustIndexWithoutChecking),
            "0" => Ok(ReferenceCheck::VerifyAgainstIndex),
            other => Err(format!(
                "{TRUST_REFERENCE_INDEX_VAR}={other:?} is not 1 or 0 — \
                 set it to 1 to skip proving the reference matches its .fai, \
                 or leave it unset to check"
            )),
        },
    }
}

/// How to say which mode a run used, for the tool's own output.
///
/// **A timing from a skipped-check run is not comparable to one from a checked run**, and
/// nothing else in a tool's output distinguishes them — so a tool that prints timings prints
/// this too, and a pasted result stays attributable.
///
/// `dead_code` is allowed because this file is included by `#[path]` into several tools and
/// only the ones that report timings call this: the dump tools deliberately do **not**, since
/// their stdout is compared byte-for-byte against a stored baseline and an extra header line
/// would break exactly the check they exist to provide.
#[allow(dead_code)]
pub fn reference_check_label(check: ReferenceCheck) -> &'static str {
    match check {
        ReferenceCheck::VerifyAgainstIndex => "verified_against_fai",
        ReferenceCheck::TrustIndexWithoutChecking => "trusted_unverified",
    }
}
