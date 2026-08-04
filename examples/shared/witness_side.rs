//! Which border of a locus a read's witness is anchored against — **the derivation the STR
//! dump tools label their `read_witness` column with**, in one body instead of three.
//!
//! Not an example itself: cargo discovers `examples/*.rs` and `examples/*/main.rs`, so a plain
//! file in a subdirectory is only compiled where an example asks for it by
//! `#[path = "shared/witness_side.rs"] mod witness_side;`.
//!
//! # Why it is shared, and why only the derivation is
//!
//! Three STR dumps carried this match with a **byte-identical seven-line comment** and had
//! already drifted apart in their string literals: `ng_ssr_loci_dump` emitted
//! `partial:left` / `partial:right` where `ng_ssr_cohort_stutter` and `ng_ssr_aligner_bakeoff`
//! emitted `partial_left` / `partial_right` — while all three said `partial:interior`, so two of
//! them mixed both separators inside one function. A consumer grepping `partial:` got one tool's
//! sides and another tool's interiors (Milestone C review, F6).
//!
//! A comment is not a mechanism; one body is. **The strings stay with the tools** — a dump's
//! output is its own contract and must not move because a sibling's did — so what is shared is
//! the *question* ("which border does this witness hold?") and not the answer's spelling. D4
//! then decided the drift rather than inheriting it: every tool now spells the colon form, which
//! is the one that was internally consistent.
//!
//! **A fourth site remains, and it cannot use this** (Milestone D naming review, which found the
//! "one body instead of three" above understating what is left). `ng_ssr_divergent_reads` derives
//! the side from a `(ReadWitness, Vec<u8>)` pair with **no locus length in hand**, so it cannot
//! call `witness_side` at all; it reads "not flush left" as "flush right", which is sound on the
//! STR path — a partial there always anchors exactly one border — and mislabels an interior run,
//! as its own comment says. Its spellings are `partialL` / `partialR`, a fifth and sixth. Folding
//! it in means threading a `LocusLen` to its call site, which is a change to that tool and not to
//! this one.
//!
//! # What it deliberately cannot answer
//!
//! Whether a **partial** witness measured the locus. `Complete` answers that for the complete
//! case and this enum carries it through; what it cannot carry is the difference between a
//! partial that is flush at both borders because it holds a hole and one that is flush at both
//! because an STR reach counted in read bases saturated — both land in `Left`, since `(true, _)`
//! matched first. **That is what [`WitnessSide::BothBorders`] fixed** (owner, 2026-07-31): the
//! case now has a name, so a tool has to spell it and a consumer has to handle it. `Complete`
//! remains the only "the length is pinned" test. See the note on `ReadWitness` itself.

use pop_var_caller::ng::locus_generation::{LocusLen, ReadWitness};

/// The border constraint a witness carries, once the side is a **derivation** rather than a
/// variant: a run flush with the left border is a prefix constraint on the allele, one flush
/// with the right border a suffix, one flush with neither is interior — and one flush with
/// **both** that is still not `Complete` is the case below, which used to be reported as a
/// prefix.
///
/// `Interior` is unreachable from the STR path, which anchors a border or yields no observation
/// at all — it is named because the generic path can mint one, and because a label that cannot
/// say "neither" would quietly report it as a prefix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WitnessSide {
    /// The read reached both borders and witnessed every position between them — **the only
    /// case that measured the locus.**
    Complete,
    /// Flush with the left border only — a prefix constraint on the allele.
    Left,
    /// Flush with the right border only — a suffix constraint.
    Right,
    /// **Flush with both borders and still not a measurement** (owner, 2026-07-31).
    ///
    /// Two different reads land here, and neither pinned the allele's length:
    ///
    /// - **A repeat read that ran out.** It anchored **one** flank, and its repeat — counted in
    ///   read bases — reached or passed the reference tract's length, so laying it down from the
    ///   anchored border covers the tract end to end. It covered every reference position and
    ///   still did not measure the allele, because the allele can run on past what the read
    ///   showed. On chromosome 1 of a tomato sample this is **2,530 of 6,216 partial
    ///   observations, 41 %** — and before this case existed, every one of them was labelled a
    ///   left-edge prefix, including the reads anchored on the *right*.
    /// - **A read blind in the middle.** It reached both borders with a hole between them, which
    ///   is what the generic path mints for a spliced read across a widened record.
    ///
    /// Keeping them together is deliberate: what a consumer must not do is read either as a
    /// measurement, and `Complete` is the test that separates measurement from lower bound.
    /// Telling the two apart *within* this case is `positions_covered` against the locus length,
    /// which the label does not need.
    BothBorders,
    /// Flush with neither.
    Interior,
}

/// Classify `witness` against the locus it was measured on.
///
/// Destructured rather than guarded on `_`, so a future `ReadWitness` variant is a compile error
/// here — the guard form is exactly what let the compiler stop forcing these sites to be
/// revisited during the rename. **The border pair is exhaustive too, for the same reason**: it
/// was `(true, _)` until the both-borders case was given a name, and that wildcard is precisely
/// what swallowed it.
pub fn witness_side(witness: &ReadWitness, locus_len: LocusLen) -> WitnessSide {
    match witness {
        ReadWitness::Complete => WitnessSide::Complete,
        run @ ReadWitness::Partial { .. } => {
            match (run.is_flush_left(), run.is_flush_right(locus_len)) {
                (true, true) => WitnessSide::BothBorders,
                (true, false) => WitnessSide::Left,
                (false, true) => WitnessSide::Right,
                (false, false) => WitnessSide::Interior,
            }
        }
    }
}
