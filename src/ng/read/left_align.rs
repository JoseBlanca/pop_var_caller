//! Step 2's v1 implementation — **pass-through + left-alignment**, producing production's
//! `PreparedRead`.
//!
//! A read whose CIGAR carries no indel is built straight from the record: nothing to shift, no
//! reference fetched, no buffer touched. A read that does carry one has its indels rewritten into
//! their leftmost equivalent spelling, so that equivalent variants get an identical one and the
//! reads supporting them pool into a single candidate instead of scattering across several weak
//! ones. Bases and qualities are copied through untouched; only the CIGAR changes.
//!
//! **The shifting is not implemented here.** It is production's `left_align_indels` (itself a port
//! of GATK's `AlignmentUtils.leftAlignIndels`), reached through the
//! [`AlignmentNormalizer`](crate::ng::alignment::AlignmentNormalizer) trait in
//! [`alignment`](crate::ng::alignment) — this module supplies the reference window, the round-trip
//! into and out of an `Alignment`, and the policy of when to fetch at all.
//!
//! Design: `doc/devel/ng/spec/read_preparation.md`, `doc/devel/ng/arch/read_preparation.md`.
