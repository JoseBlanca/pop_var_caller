//! **Stage 1 of the parity oracle: ng's walker *is* production's walker.**
//!
//! One `PreparedRead` stream, two walkers, and the two
//! `Result<PileupRecord, WalkerError>` streams must be equal element for element — plus
//! the `RunSummary`s. That is the whole claim of this plan
//! (`doc/devel/ng/spec/locus_generation_pileup.md` §3, §13.1), and it is what every
//! deliberate divergence in plan 3 will be measured against.
//!
//! # Why "byte-identical" is well defined here
//!
//! ng's copy still emits **production's** [`PileupRecord`], so the comparison is
//! `assert_eq!` on one type rather than a projection. Its hand-written `PartialEq`
//! compares the two `f32`s **by bits** ([pileup_record.rs](crate::pileup_record)), so the
//! `NaN` placeholders `finalise()` writes for `windowed_gc` / `windowed_coverage` compare
//! equal and the comparison is total.
//!
//! # One stream, fed to both
//!
//! Cases are generated as production's `PreparedRead` and converted to ng's through
//! `PreparedRead::from_production`, so the two walkers see the same bytes in the same
//! order by construction. Preparing separately would inject read preparation's uppercase
//! divergence (`read_preparation.md` §6) into a comparison that is about the *walk*.
//! Each walk builds its own `MockFasta` from the same bytes — the type is stateless, so
//! that is equivalent to sharing one; the real-data test lends a single reference to both
//! because *its* accessor is not stateless.
//!
//! # The two sides no longer reach the reference the same way — A0
//!
//! Production's walker takes a `MultiChromRefFetcher`; ng's takes a
//! [`RefSeq`](crate::ng::ref_seq::RefSeq). The **same `MockFasta` value** serves both,
//! through its two impls — production's own, and the canonicalising one in
//! [`mock_reference`](super::mock_reference) — so the bytes are shared by construction
//! rather than by two fixtures agreeing. Where they could still part company is
//! canonicalisation, which is the identity only because the generator draws from
//! `ACGTN`: `both_sides_of_the_differential_are_served_the_same_bytes` pins that, so a
//! future generator change introducing a lower-case or ambiguity-coded base fails there
//! and not as "the walkers disagree".
//!
//! A0 also changes ng's `WalkerError::Fasta` to carry a `RefSeqError` where production's
//! carries a `ChromRefFetchError`, so the two error streams can no longer be compared by
//! `Debug` alone. `render_*_error` normalises **that one variant** and leaves the other
//! eight rendered verbatim; `both_walkers_report_the_same_error_on_the_same_malformed_input`
//! carries a fixture that reaches it.
//!
//! # It has been shown to fail — B2, re-run 2026-07-29
//!
//! A differential that has only ever passed is a claim. Each behaviour the plan names was
//! mutated **in ng's copy**, one at a time, and `ng_holds_the_same_evidence_as_production_on_complete_reads` was
//! required to fail; then the mutation was reverted and the differential re-run green.
//! Every one died inside the **first nineteen cases of the first seed**, against a default
//! of 400 cases × 4 seeds.
//!
//! | # | behaviour | mutation applied to ng's copy | first divergence |
//! |---|---|---|---|
//! | 1 | mate-overlap reconciliation | early `return` from `genome_walk::resolve_mate_overlap_at_pos` | seed 0 case 2, item 6 |
//! | 2 | adaptor masking | `cigar_cursor::base_in_adaptor` always `false` | seed 0 case 0 — 27 records against production's 25 |
//! | 3 | record widening | `open_record::widen` extends only `alleles[0]`, not every bucket | seed 0 case 18, item 17 |
//! | 4 | the subtract-then-add re-fold | the `subtract_contribution` half dropped, so a re-folding read double-counts | seed 0 case 0, item 16 |
//! | 5 | the column depth cap | the `truncate(cap)` removed, counter left incrementing | seed 0 case 0, item 10 |
//! | 6 | **the panic *cause*** | ng's copy of the reachable `debug_assert!` replaced by an unrelated `panic!` | seed 0 case 15 — "the two walkers did not stop the same way" |
//!
//! Two of these are worth keeping in mind rather than skimming.
//!
//! **Mutation 5** leaves `column_depth_truncations` incrementing, so a differential that
//! compared only the `RunSummary` — a plausible way to write this harness — would have
//! passed it. The records caught it.
//!
//! **Mutation 6 is here because a review demonstrated it passing.** The first version of
//! this harness stored `catch_unwind(..).is_err()`, so "both stopped" counted as
//! agreement; replacing ng's copy of production's reachable `debug_assert!` with a
//! semantically unrelated `panic!` left the *entire* parity module green, including the
//! test whose own doc claimed to check that ng reaches the same precondition.
//! `open_record.rs` alone carries eight distinct `debug_assert!`s. `WalkOutcome` now
//! carries the panic **message**.
//!
//! # Run at scale — B3, re-run 2026-07-29 after the review. **Zero divergences everywhere.**
//!
//! | run | scale | result |
//! |---|---|---|
//! | synthetic, release, `PVC_PARITY_CASES=5000` | 20,000 cases | **1,010,515 records, 0 divergences** |
//! | synthetic, debug, `PVC_PARITY_CASES=2500` | 10,000 cases | **505,979 records, 0 divergences**, **0 panics** |
//! | GIAB HG002 10×, `chr1:1000000-1400000` | targeted TR bundle | **4,600 records, 0 divergences** |
//! | GIAB HG002 300×, `chr1:100000000-120000000` | 20 Mb | **137,591 records, 0 divergences** |
//! | tomato CRAM `SRR7279481.p1`, `SL4.0ch01:3406886-3506886` | 100 kb | **96,260 records, 0 divergences** |
//! | tomato CRAM `SRR7279481.p1`, `SL4.0ch01:13806669-15092603` | 1.3 Mb | **198,673 records, 0 divergences** |
//!
//! 437,124 records of real sequencing data across two organisms, a BAM and a CRAM, 10× and
//! 300×; 1.5 M more synthetic.
//!
//! **The debug row used to read "415 cases (4.2%) panicked — in *both* walkers".** That was
//! the production defect this harness found; it is now **fixed in both copies** (see
//! `a_deletion_anchored_before_its_record_contributes_none_of_the_bases_it_deleted`), so the
//! rate is zero and 13,755 more records are compared in debug than before, because those
//! cases now run to completion instead of being truncated at the panic.
//!
//! The panic channel stays, and it is not vestigial: it is what verifies the fix was applied
//! **identically to both copies**. A fix in one and not the other would show up here as a
//! panic-message divergence rather than as a silent difference in bases.
//!
//! # This harness dies in plan 3, by design
//!
//! Plan 3 makes the two walkers differ on purpose — the no-fill haplotype builder,
//! REF-only widening, the read-group split. What survives is narrower: loci where every
//! folded read witnessed the whole footprint must agree forever. **So this is the last
//! moment the baseline can be banked**, and `the_generator_exercises_what_the_port_can_break`
//! is what says the banking is worth anything: a differential that has never been shown
//! to fail is a claim, not evidence.

use std::rc::Rc;
use std::sync::Arc;

use super::tests::MockFasta;
use super::{PreparedRead, WalkerConfig};
use crate::ng::read::PLACEHOLDER_READ_GROUP;
// Aliased, so which walker a call reaches is legible at the call site rather than carried
// by a `super::` — this file is the one place both are in scope at once.
use super::super::SampleLocusObservations;
use crate::pileup::walker::{
    CigarOp, MateRole as ProductionMateRole, PreparedRead as ProductionPreparedRead,
    run as production_run,
};
use crate::pileup_record::{AlleleSupportStats, ChainId, PileupRecord};

/// A deterministic PRNG, so a failure is reproducible from its seed alone.
///
/// SplitMix64 — the same generator the delimiter and left-alignment parity harnesses use,
/// written out rather than pulled from a crate because a parity fixture must be
/// reproducible from the source in front of you, not from a dependency's version.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }

    /// `true` with probability `1 / n`.
    fn one_in(&mut self, n: usize) -> bool {
        self.below(n) == 0
    }

    fn base(&mut self) -> u8 {
        b"ACGTN"[self.below(5)]
    }
}

/// The two fixture contigs. Two rather than one so the chromosome-boundary path — where
/// production flushes the open-record table and the allocator keeps its counter — is on
/// the walk rather than assumed.
const CONTIG_LENGTH: usize = 160;
const CONTIGS: usize = 2;

/// One generated case: a reference, a coordinate-sorted read stream, and the config to
/// walk them under.
struct Case {
    reference: Vec<String>,
    reads: Vec<ProductionPreparedRead>,
    config: WalkerConfig,
    /// Reads carrying an adaptor boundary that actually falls inside their own span —
    /// the only ones for which the G1 filter can silence a base. Counted at generation
    /// because the summary has no counter for it.
    reads_with_live_adaptor_boundary: usize,
}

impl Case {
    fn fasta(&self) -> MockFasta {
        MockFasta::with_chromosomes(
            &self
                .reference
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
        )
    }
}

/// Build one read: a placement, a CIGAR, and bases mostly agreeing with the reference.
///
/// The invariant `PreparedRead::length()` checks — `seq.len() == bq_baq.len()` and the
/// read-consuming CIGAR ops summing to `seq.len()` — is maintained by construction, since
/// a read that violates it is rejected at admission and the walk stops there. What the
/// generator *does* vary is everything the fold reads: op mix, strand, quality, mapping
/// quality, and whether an adaptor boundary silences part of the read.
fn generate_read(
    rng: &mut SplitMix64,
    reference: &[String],
    chrom_id: u32,
    start: u32,
    qname: &str,
    mate_role: ProductionMateRole,
) -> (ProductionPreparedRead, bool) {
    let contig = reference[chrom_id as usize].as_bytes();
    let mut cigar: Vec<CigarOp> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut ref_pos = start as usize - 1; // 0-based cursor into the contig

    // A leading hard clip on some reads: it consumes **neither** reference nor read, so it
    // is the op that probes the cursor's offset-table walk without moving either cursor.
    if rng.one_in(6) {
        cigar.push(CigarOp::HardClip(1 + rng.below(3) as u32));
    }

    // A leading soft clip on some reads: it consumes read bases but no reference, which
    // is the offset arithmetic the cursor is most easily got wrong on.
    if rng.one_in(4) {
        let clip = 1 + rng.below(3);
        cigar.push(CigarOp::SoftClip(clip as u32));
        seq.extend((0..clip).map(|_| rng.base()));
    }

    // Between one and four blocks of match/indel, bounded so the read stays inside the
    // contig with room for a record to widen past it.
    let blocks = 1 + rng.below(4);
    for _ in 0..blocks {
        let remaining = contig.len().saturating_sub(ref_pos);
        if remaining < 12 {
            break;
        }
        let matched = 2 + rng.below(8);
        // `M`, `=` and `X` are the **same op to the walker** — they share the `Match` arm at
        // four sites in the cursor — which is exactly why a copy that treated them
        // differently has to be caught here rather than assumed. `=`/`X` are what minimap2
        // `--eqx` and DRAGEN emit, so this is not an exotic input class.
        cigar.push(match rng.below(4) {
            0 => CigarOp::SeqMatch(matched as u32),
            1 => CigarOp::SeqMismatch(matched as u32),
            _ => CigarOp::Match(matched as u32),
        });
        for offset in 0..matched {
            // Mostly the reference base, so REF alleles dominate as they do in real data;
            // a substitution now and then, so records carry more than one allele.
            let base = if rng.one_in(6) {
                rng.base()
            } else {
                contig[ref_pos + offset]
            };
            seq.push(base);
        }
        ref_pos += matched;

        match rng.below(10) {
            // A deletion — this is what widens a record past the reads that opened it,
            // and the widen path is one of the five behaviours B2 mutates.
            0 | 1 => {
                let deleted = 1 + rng.below(6);
                if contig.len().saturating_sub(ref_pos) > deleted + 4 {
                    cigar.push(CigarOp::Deletion(deleted as u32));
                    ref_pos += deleted;
                }
            }
            // An insertion — footprint of one reference position, several bases.
            2 => {
                let inserted = 1 + rng.below(3);
                cigar.push(CigarOp::Insertion(inserted as u32));
                seq.extend((0..inserted).map(|_| rng.base()));
            }
            // Padding — like a hard clip, it consumes neither axis. Emitted here rather
            // than at an end so it is never the first or last op, where the walker's
            // first/last-op rules would make it uninteresting.
            4 => cigar.push(CigarOp::Padding(1 + rng.below(2) as u32)),
            // A reference skip — emits nothing and lets both flanks emit independently.
            3 => {
                let skipped = 1 + rng.below(4);
                if contig.len().saturating_sub(ref_pos) > skipped + 4 {
                    cigar.push(CigarOp::Skip(skipped as u32));
                    ref_pos += skipped;
                }
            }
            _ => {}
        }
    }

    // A read with no reference-consuming op is rejected at admission (`ZeroRefSpan`),
    // which would end both walks at the same place but stop the case testing anything.
    // The **position** is clamped, not just the base index: a read whose `alignment_end`
    // ran past the contig would be rejected with a fatal `WalkerError::Fasta`, which
    // silently converts a fold test into a rejection-path test.
    if !cigar.iter().any(|op| {
        matches!(
            op,
            CigarOp::Match(_)
                | CigarOp::SeqMatch(_)
                | CigarOp::SeqMismatch(_)
                | CigarOp::Deletion(_)
                | CigarOp::Skip(_)
        )
    }) {
        ref_pos = ref_pos.min(contig.len() - 1);
        cigar.push(CigarOp::Match(1));
        seq.push(contig[ref_pos]);
        ref_pos += 1;
    }
    debug_assert!(
        ref_pos <= contig.len(),
        "a generated read must end inside its contig"
    );

    let alignment_end = ref_pos as u32; // 1-based inclusive == 0-based exclusive
    let is_reverse_strand = rng.one_in(2);

    // The adaptor boundary, and whether it is **live** — i.e. whether the G1 filter will
    // actually silence a base on this read.
    //
    // "Inside the alignment span" is the wrong predicate, and an earlier version used it
    // (worse: it returned `true` unconditionally, so the counter was just `is_some()`).
    // `base_in_adaptor` is consulted **only at Match-emit sites**, so a forward-strand
    // boundary landing in a read's trailing `D`/`N` tail silences nothing. Liveness is
    // therefore computed the way the cursor computes it: over the positions a `Match`,
    // `=` or `X` will actually emit.
    let adaptor_boundary = if rng.one_in(4) {
        let span = alignment_end.saturating_sub(start);
        Some(if span >= 2 {
            start + 1 + (rng.below(span as usize - 1) as u32)
        } else {
            start
        })
    } else {
        None
    };
    let live_boundary = adaptor_boundary.is_some_and(|boundary| {
        let mut pos = start;
        cigar.iter().any(|op| match *op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                let emitted = pos..pos + n;
                pos += n;
                emitted.into_iter().any(|p| {
                    if is_reverse_strand {
                        p <= boundary
                    } else {
                        p >= boundary
                    }
                })
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                pos += n;
                false
            }
            _ => false,
        })
    });

    let read = ProductionPreparedRead {
        chrom_id,
        alignment_start: start,
        alignment_end,
        cigar,
        bq_baq: (0..seq.len()).map(|_| (rng.below(41)) as u8).collect(),
        seq,
        mq_log_err: -(1.0 + rng.below(40) as f64 / 10.0),
        mapq: rng.below(61) as u8,
        is_reverse_strand,
        qname: Arc::from(qname),
        mate_role,
        adaptor_boundary,
    };
    (read, live_boundary)
}

/// Build one case: a reference, a read stream, and a config.
fn generate(rng: &mut SplitMix64) -> Case {
    let reference: Vec<String> = (0..CONTIGS)
        .map(|_| {
            (0..CONTIG_LENGTH)
                .map(|_| b"ACGT"[rng.below(4)] as char)
                .collect()
        })
        .collect();

    let mut reads: Vec<ProductionPreparedRead> = Vec::new();
    let mut reads_with_live_adaptor_boundary = 0usize;
    let read_count = 2 + rng.below(14);

    for index in 0..read_count {
        // Reads cluster near the contig start so columns get deep enough for the depth
        // cap to bite; a scattered stream would leave every column at depth one or two.
        let chrom_id = if CONTIGS > 1 && rng.one_in(6) { 1 } else { 0 };
        // Most reads cluster near the contig start so columns get deep enough for the
        // depth cap to bite. **A minority sit at the far end**, where a record widened by
        // the longest possible deletion produces the fetch most likely to be off by one —
        // without them the last third of every contig is never touched and the bounds
        // guards below never fire even once. `- 6` leaves room for a second mate's offset,
        // so no read is placed past the contig.
        let start = if rng.one_in(8) {
            (CONTIG_LENGTH - 6) as u32 + rng.below(6) as u32
        } else {
            1 + rng.below(CONTIG_LENGTH / 3) as u32
        };

        if rng.one_in(3) {
            // A pair. Placed to overlap, because mate-overlap reconciliation is detected
            // by shared chain id at a *shared position* — two mates that never meet
            // exercise the allocator but not the reconciliation.
            let qname = format!("pair{index}");
            let (first, live_a) = generate_read(
                rng,
                &reference,
                chrom_id,
                start,
                &qname,
                ProductionMateRole::FirstOfPair,
            );
            let overlap_start = start + (rng.below(6) as u32);
            let (second, live_b) = generate_read(
                rng,
                &reference,
                chrom_id,
                overlap_start,
                &qname,
                ProductionMateRole::SecondOfPair,
            );
            reads_with_live_adaptor_boundary += usize::from(live_a) + usize::from(live_b);
            reads.push(first);
            reads.push(second);
        } else {
            let qname = format!("solo{index}");
            let (read, live) = generate_read(
                rng,
                &reference,
                chrom_id,
                start,
                &qname,
                ProductionMateRole::Solo,
            );
            reads_with_live_adaptor_boundary += usize::from(live);
            reads.push(read);
        }
    }

    // The walker's coordinate-order invariant is fatal (`WalkerError::OutOfOrder`), and a
    // run that ends at the first read would test nothing beyond it. A **stable** sort, so
    // the admission order of equal-position reads is a property of the generator rather
    // than of the sort — which matters because the depth cap truncates in admission order.
    reads.sort_by_key(|read| (read.chrom_id, read.alignment_start));

    // Most cases run production's defaults. A third run tiny column caps, because the cap
    // never fires at these depths otherwise — and "the cap is on the walk" is one of the
    // five things B2 mutates.
    let mut config = WalkerConfig::default();
    if rng.one_in(3) {
        config.max_snp_column_depth = 1 + rng.below(4) as u32;
        config.max_indel_column_depth = 1 + rng.below(2) as u32;
    }
    // The window only bites when it is smaller than the gap between a first mate and the
    // next admitted read, and mates are placed within 5 bp of each other — so a window
    // drawn up to 40 almost never evicts. Drawn small, so the path where a pair silently
    // degrades to two solos (different chain id, no reconciliation, a different fold) is
    // actually on the walk; `the_generator_exercises_what_the_port_can_break` asserts it.
    if rng.one_in(8) {
        config.mate_lookup_window = 1 + rng.below(8) as u32;
    }

    Case {
        reference,
        reads,
        config,
        reads_with_live_adaptor_boundary,
    }
}

/// Build a read that **witnesses every reference position it spans** — the class ng's
/// no-fabrication rule does not change, and therefore the class the two walkers must
/// agree on forever.
///
/// Four things are excluded, and each is excluded because it makes a read *silent* at a
/// position its alignment covers — which is precisely the fabrication primitive: "no event
/// → reference base", wider than "outside the read's span" (spec §6).
///
/// - **`N` read bases** — the cursor emits no `Match` for one.
/// - **an adaptor boundary** — the G1 filter silences every base past it.
/// - **a reference skip** — it consumes reference and emits nothing.
/// - **an indel as the first or last op** — the first/last-op rule drops it, leaving its
///   reference footprint unwitnessed.
///
/// Everything else stays: substitutions, insertions, deletions, clips, padding, both
/// strands, quality variation. What makes the *record* class complete is that the read
/// runs from position 1 to the contig's end, so every record the walk opens lies inside
/// its span and every position of every such record carries one of its events.
fn generate_complete_read(
    rng: &mut SplitMix64,
    reference: &[String],
    chrom_id: u32,
    qname: &str,
    mate_role: ProductionMateRole,
) -> ProductionPreparedRead {
    let contig = reference[chrom_id as usize].as_bytes();
    let mut cigar: Vec<CigarOp> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut ref_pos: usize = 0; // 0-based cursor into the contig

    // A leading soft clip on some reads: it consumes read bases but no reference, so it
    // silences nothing. Never an indel here — the first/last-op rule would drop that.
    if rng.one_in(4) {
        let clip = 1 + rng.below(3);
        cigar.push(CigarOp::SoftClip(clip as u32));
        seq.extend((0..clip).map(|_| b"ACGT"[rng.below(4)]));
    }

    while ref_pos < contig.len() {
        let remaining = contig.len() - ref_pos;
        // The final block always runs to the contig's end, so the read spans it whole.
        let matched = if remaining <= 12 {
            remaining
        } else {
            2 + rng.below(8)
        };
        cigar.push(match rng.below(4) {
            0 => CigarOp::SeqMatch(matched as u32),
            1 => CigarOp::SeqMismatch(matched as u32),
            _ => CigarOp::Match(matched as u32),
        });
        for offset in 0..matched {
            // `ACGT` only: an `N` would silence this position, which is the one thing
            // this generator exists to exclude.
            let base = if rng.one_in(6) {
                b"ACGT"[rng.below(4)]
            } else {
                contig[ref_pos + offset]
            };
            seq.push(base);
        }
        ref_pos += matched;

        // An indel only while enough reference is left for a further match block, so it
        // is never the last op.
        if contig.len().saturating_sub(ref_pos) > 12 {
            match rng.below(6) {
                0 | 1 => {
                    let deleted = 1 + rng.below(6);
                    cigar.push(CigarOp::Deletion(deleted as u32));
                    ref_pos += deleted;
                }
                2 => {
                    let inserted = 1 + rng.below(3);
                    cigar.push(CigarOp::Insertion(inserted as u32));
                    seq.extend((0..inserted).map(|_| b"ACGT"[rng.below(4)]));
                }
                // Padding consumes neither axis, so it silences nothing.
                3 => cigar.push(CigarOp::Padding(1 + rng.below(2) as u32)),
                _ => {}
            }
        }
    }

    ProductionPreparedRead {
        chrom_id,
        alignment_start: 1,
        alignment_end: contig.len() as u32,
        cigar,
        bq_baq: (0..seq.len()).map(|_| 20 + rng.below(21) as u8).collect(),
        seq,
        mq_log_err: -3.0 - (rng.below(4) as f64),
        mapq: 20 + rng.below(41) as u8,
        is_reverse_strand: rng.one_in(2),
        qname: Arc::from(qname),
        mate_role,
        // No boundary: the G1 filter would silence bases the alignment spans.
        adaptor_boundary: None,
    }
}

/// A case in which **every folded read witnessed every position of every record it folded
/// into** — the permanent anchor's fixture.
///
/// Depth is kept high enough that the column cap can bite and that mate overlap and
/// re-folds after a widen are on the walk; the mates are placed at the same start rather
/// than offset, because a read that begins later would not span the records before it.
fn generate_complete(rng: &mut SplitMix64) -> Case {
    let reference: Vec<String> = (0..CONTIGS)
        .map(|_| {
            (0..CONTIG_LENGTH)
                .map(|_| b"ACGT"[rng.below(4)] as char)
                .collect()
        })
        .collect();

    let mut reads: Vec<ProductionPreparedRead> = Vec::new();
    let read_count = 2 + rng.below(10);
    for index in 0..read_count {
        let chrom_id = if CONTIGS > 1 && rng.one_in(6) { 1 } else { 0 };
        if rng.one_in(3) {
            let qname = format!("pair{index}");
            for role in [
                ProductionMateRole::FirstOfPair,
                ProductionMateRole::SecondOfPair,
            ] {
                reads.push(generate_complete_read(
                    rng, &reference, chrom_id, &qname, role,
                ));
            }
        } else {
            let qname = format!("solo{index}");
            reads.push(generate_complete_read(
                rng,
                &reference,
                chrom_id,
                &qname,
                ProductionMateRole::Solo,
            ));
        }
    }
    reads.sort_by_key(|read| (read.chrom_id, read.alignment_start));

    let mut config = WalkerConfig::default();
    if rng.one_in(3) {
        config.max_snp_column_depth = 1 + rng.below(4) as u32;
        config.max_indel_column_depth = 1 + rng.below(2) as u32;
    }

    Case {
        reference,
        reads,
        config,
        // No read carries one, by construction — see `generate_complete_read`.
        reads_with_live_adaptor_boundary: 0,
    }
}

/// The eight `RunSummary` counters, **named**.
///
/// Named rather than an `[u64; 8]` beside a parallel `[&str; 8]`, because two arrays that
/// must stay in step are two arrays that can drift — and the failure mode is a divergence
/// report that names the wrong counter, which is worse than no name at all.
#[derive(Debug, PartialEq, Eq)]
struct SummaryCounters {
    reads_admitted: u64,
    records_emitted: u64,
    record_widen_events: u64,
    mate_overlap_positions: u64,
    chain_allocations: u64,
    active_reads_high_water: u64,
    mate_lookup_evictions: u64,
    column_depth_truncations: u64,
}

/// Production's `RunSummary`, named. **Exhaustively destructured, no `..`**: a ninth field
/// on production's summary must stop this compiling, or it would silently leave the parity
/// claim. `RunSummary::merge` destructures for exactly this reason.
fn production_counters(summary: crate::pileup::walker::RunSummary) -> SummaryCounters {
    let crate::pileup::walker::RunSummary {
        reads_admitted,
        records_emitted,
        record_widen_events,
        mate_overlap_positions,
        chain_allocations,
        active_reads_high_water,
        mate_lookup_evictions,
        column_depth_truncations,
    } = summary;
    SummaryCounters {
        reads_admitted,
        records_emitted,
        record_widen_events,
        mate_overlap_positions,
        chain_allocations,
        active_reads_high_water: u64::from(active_reads_high_water),
        mate_lookup_evictions,
        column_depth_truncations,
    }
}

/// ng's copy of `RunSummary`, named. The field list appears twice because the two types are
/// nominally distinct and Rust cannot bound on field access — and writing it twice is what
/// makes *both* exhaustive.
fn ng_counters(summary: super::RunSummary) -> SummaryCounters {
    let super::RunSummary {
        reads_admitted,
        records_emitted,
        record_widen_events,
        mate_overlap_positions,
        chain_allocations,
        active_reads_high_water,
        mate_lookup_evictions,
        column_depth_truncations,
    } = summary;
    SummaryCounters {
        reads_admitted,
        records_emitted,
        record_widen_events,
        mate_overlap_positions,
        chain_allocations,
        active_reads_high_water: u64::from(active_reads_high_water),
        mate_lookup_evictions,
        column_depth_truncations,
    }
}

/// What one walk produced: its record/error stream, its summary, and — if it stopped —
/// **why**.
///
/// The panic channel is not defensiveness. Production's `apply_events_to_ref_into` carries
/// a `debug_assert!` that every event's anchor is at or after the record's own position,
/// and that precondition is **reachable on a legal read stream** — see
/// `both_walkers_panic_on_a_deletion_anchored_before_its_record` for the three-read case
/// and the mechanism. So on those inputs a debug build panics before either walker can
/// finish, and a harness that could not represent that would have to exclude exactly the
/// long-deletion inputs this port exists to get right.
///
/// **It is the panic *message*, not a `bool`.** An earlier version stored
/// `catch_unwind(..).is_err()`, and a review demonstrated the hole: replacing ng's copy of
/// that `debug_assert!` with an unrelated `panic!` left the whole suite green, because
/// "both stopped" is not "both reached the same precondition". `open_record.rs` alone
/// carries eight distinct `debug_assert!`s. The message is compared and the *location* is
/// not — the two copies share the format string verbatim, but ng panics from
/// `src/ng/locus_generation/pileup/` and production from `src/pileup/walker/`.
struct WalkOutcome {
    records: Vec<Result<PileupRecord, String>>,
    /// `None` when the walk panicked — there is no summary to read from a walker that
    /// unwound out from under us.
    summary: Option<SummaryCounters>,
    /// `Some(message)` when the walk panicked.
    panic_message: Option<String>,
}

/// A reference-fetch failure, said in terms **both** sides can state.
///
/// A0 gave ng's `WalkerError::Fasta` a [`RefSeqError`](crate::ng::ref_seq::RefSeqError)
/// source where production's carries a `ChromRefFetchError`. The two are
/// variant-for-variant equivalents but nominally distinct types, so `Debug` renders them
/// differently and a raw `{:?}` comparison would report every fetch failure as a
/// divergence — including the ones where the two walkers agree exactly.
///
/// **The contig identifier is deliberately not part of this.** Production names contigs
/// by name and ng by `ContigId`, and there is nothing to compare; `WalkerError::Fasta`
/// carries `chrom_id` *outside* the source, where it is compared verbatim.
///
/// One collapse is worth naming: an unknown contig reaches production as
/// `Io { NotFound }` and ng as `UnknownContig`, because that is the convention
/// `MockFasta` and its ng-side view use (`mock_reference.rs`). Both land on
/// [`FetchFailure::UnknownContig`], which loses the distinction between "the contig is
/// not in this reference" and "a file went missing mid-run" — neither of which any
/// fixture here can produce by another route.
#[derive(Debug, PartialEq, Eq)]
enum FetchFailure {
    OutOfBounds {
        contig_length: u64,
        start: u64,
        end: u64,
    },
    InvalidStart,
    UnknownContig,
    Io(std::io::ErrorKind, String),
    /// A variant with no counterpart on the other side. Rendered rather than dropped, so
    /// it can never compare equal to anything by accident.
    Unmatched(String),
}

fn production_fetch_failure(error: &crate::fasta::ChromRefFetchError) -> FetchFailure {
    use crate::fasta::ChromRefFetchError as E;
    match error {
        E::OutOfBounds {
            chrom_length,
            start,
            end,
            chrom_name: _,
        } => FetchFailure::OutOfBounds {
            contig_length: u64::from(*chrom_length),
            start: u64::from(*start),
            end: u64::from(*end),
        },
        E::InvalidStart => FetchFailure::InvalidStart,
        E::Io {
            source,
            chrom_name: _,
        } if source.kind() == std::io::ErrorKind::NotFound => FetchFailure::UnknownContig,
        E::Io {
            source,
            chrom_name: _,
        } => FetchFailure::Io(source.kind(), source.to_string()),
        // The streaming fetcher's; no `RefSeq` implementation has an equivalent.
        E::OutOfPattern { .. } => FetchFailure::Unmatched(format!("{error:?}")),
    }
}

fn ng_fetch_failure(error: &crate::ng::ref_seq::RefSeqError) -> FetchFailure {
    use crate::ng::ref_seq::RefSeqError as E;
    match error {
        E::OutOfBounds {
            contig_length,
            start,
            end,
            contig: _,
        } => FetchFailure::OutOfBounds {
            contig_length: *contig_length,
            start: *start,
            end: *end,
        },
        E::InvalidStart => FetchFailure::InvalidStart,
        E::UnknownContig(_) => FetchFailure::UnknownContig,
        E::Io { source, contig: _ } => FetchFailure::Io(source.kind(), source.to_string()),
    }
}

/// Production's `WalkerError`, rendered for comparison against ng's.
///
/// Only `Fasta` is rewritten — the other eight variants are structurally identical in
/// the two enums and their derived `Debug` output prints no module path, so it is the
/// same string on both sides. They are **listed by name rather than matched with `_`**,
/// so a variant added to either enum stops this compiling instead of silently taking the
/// verbatim path.
fn render_production_error(error: &crate::pileup::walker::WalkerError) -> String {
    use crate::pileup::walker::WalkerError as E;
    match error {
        E::Fasta {
            chrom_id,
            start,
            start_plus_len,
            source,
        } => render_fasta_error(
            *chrom_id,
            *start,
            *start_plus_len,
            production_fetch_failure(source),
        ),
        E::OutOfOrder { .. }
        | E::ZeroRefSpan { .. }
        | E::ActiveReadsExhausted { .. }
        | E::ChainIdSpaceExhausted { .. }
        | E::PendingMatesExhausted { .. }
        | E::RecordTooWide { .. }
        | E::Internal { .. }
        | E::MalformedRead { .. } => format!("{error:?}"),
    }
}

/// ng's `WalkerError`, rendered the same way. Written out a second time rather than
/// shared: the two enums are nominally distinct, and writing it twice is what makes
/// *both* matches exhaustive — the same reason `production_counters` and `ng_counters`
/// repeat their field lists.
fn render_ng_error(error: &super::WalkerError) -> String {
    use super::WalkerError as E;
    match error {
        E::Fasta {
            chrom_id,
            start,
            start_plus_len,
            source,
        } => render_fasta_error(*chrom_id, *start, *start_plus_len, ng_fetch_failure(source)),
        E::OutOfOrder { .. }
        | E::ZeroRefSpan { .. }
        | E::ActiveReadsExhausted { .. }
        | E::ChainIdSpaceExhausted { .. }
        | E::PendingMatesExhausted { .. }
        | E::RecordTooWide { .. }
        | E::Internal { .. }
        | E::MalformedRead { .. } => format!("{error:?}"),
    }
}

fn render_fasta_error(
    chrom_id: u32,
    start: u32,
    start_plus_len: u32,
    cause: FetchFailure,
) -> String {
    format!(
        "Fasta {{ chrom_id: {chrom_id}, start: {start}, start_plus_len: {start_plus_len}, cause: {cause:?} }}"
    )
}

/// Render a panic payload as its message, which is the part the two walkers must share.
fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "<non-string panic payload>".to_string()
    }
}

/// Drive one walker to exhaustion, surviving a panic with whatever it emitted first.
///
/// Generic over the walker, so production's and ng's — nominally distinct types with no
/// shared trait — go through **one** body. An earlier version had this block four times and
/// the two inlined copies had already drifted.
///
/// `records` is filled inside the closure and read after the unwind: `Vec::push` is
/// exception-safe, so on a panic it holds exactly the prefix that was emitted, which is
/// what makes the comparison element-wise rather than all-or-nothing.
///
/// Errors are rendered to a `String` rather than compared as values: the two
/// `WalkerError` types are nominally distinct, so they cannot be compared directly. The
/// rendering is mostly `{:?}` — which shows the variant and every field, where `Display`
/// shows only what the format string chose — with `Fasta`'s source normalised, since A0
/// gave the two enums different source types (`render_production_error` /
/// `render_ng_error`). The renderer is a **parameter** rather than a `Debug` bound, so
/// which side is being rendered is a fact at the call site and neither side can quietly
/// fall back to a bound the other does not satisfy.
fn drive_production<W, E>(
    mut walker: W,
    render_error: impl Fn(&E) -> String,
    summary_of: impl FnOnce(&W) -> SummaryCounters,
) -> WalkOutcome
where
    W: Iterator<Item = Result<PileupRecord, E>>,
{
    let mut records = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(item.map_err(|error| render_error(&error)));
        }
        summary = Some(summary_of(&walker));
    }))
    .err()
    .map(panic_message);
    WalkOutcome {
        records,
        summary,
        panic_message,
    }
}

/// The same as [`drive_production`], for a walker that yields ng's own locus type — each record laid
/// back out as a [`PileupRecord`] by [`to_pileup_record`](super::to_pileup_record) so the
/// two streams stay comparable. Separate rather than generic over the item type: the
/// projection is the thing worth seeing at the call site.
fn drive_ng<W, E>(
    mut walker: W,
    render_error: impl Fn(&E) -> String,
    summary_of: impl FnOnce(&W) -> SummaryCounters,
) -> WalkOutcome
where
    W: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    let mut records = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(
                item.map(|locus| super::to_pileup_record(&locus))
                    .map_err(|error| render_error(&error)),
            );
        }
        summary = Some(summary_of(&walker));
    }))
    .err()
    .map(panic_message);
    WalkOutcome {
        records,
        summary,
        panic_message,
    }
}

/// Production's answer.
fn production_walk(case: &Case) -> WalkOutcome {
    let fasta = case.fasta();
    drive_production(
        production_run(case.reads.clone(), &fasta, &case.config),
        render_production_error,
        |walker| production_counters(walker.summary()),
    )
}

/// ng's answer, over the same reads converted to ng's read type.
///
/// The `MockFasta` goes in **by value** where production's went in by reference:
/// `MultiChromRefFetcher` has a blanket impl for `&T` and `RefSeq` has none. The type is
/// stateless and cheap, so this is a call-shape difference and not a fixture difference —
/// and `both_sides_of_the_differential_are_served_the_same_bytes` checks the bytes rather
/// than trusting that.
fn ng_walk(case: &Case) -> WalkOutcome {
    let fasta = case.fasta();
    let reads: Vec<PreparedRead> = case
        .reads
        .iter()
        .cloned()
        // The read group plays no part in what the copy computes — it is carried for
        // plan 3 — so a placeholder is right here, and a *varying* one would be a
        // difference between the two streams rather than a property under test.
        .map(|read| PreparedRead::from_production(read, PLACEHOLDER_READ_GROUP))
        .collect();
    // **ng emits `SampleLocusObservations` from B2 on; the differential still compares
    // `PileupRecord`s.** Laid back out through `to_pileup_record` rather than projecting
    // production forward, which is D1's job and needs the five divergence classes named
    // first. What the projection cannot carry is listed on it, and the two that matter here
    // are named in `record_evidence` (`placed_start`, no longer computed at all) and in
    // `classify_record` (chain ids, now dropped per read).
    drive_ng(super::run(reads, fasta, &case.config), render_ng_error, {
        |walker| ng_counters(walker.summary())
    })
}

/// Project a record onto the surface the two walkers can still be compared on — **two
/// named differences, both of them projections rather than excuses** (spec §3, classes 4
/// and 5).
///
/// 1. **Drop every non-REF bucket no read is folded into.** From A3 `widen` extends only
///    `alleles[0]`, so a read that re-folds after a widen lands in a new bucket and leaves
///    the old one behind at `num_obs == 0`; ng evicts those at the end of every fold, and
///    production keeps them, because its append-to-every-bucket kept the re-fold landing
///    back where it started. Neither side's *evidence* moves — an empty bucket supports no
///    read and carries no chain id.
/// 2. **Sort the non-REF buckets by their bases.** Production's order is bucket-creation
///    order; ng's creation order now changes with eviction, and B2 makes ng sort before
///    emitting anyway, for a reason that outlives this comparison: its rows come from an
///    `AHashMap` whose iteration order is seeded per process.
///
/// `alleles[0]` is left where it is whatever its support: it is the record's reference
/// sequence, production creates it with zero observations by design, and moving it would
/// break the positional REF invariant both sides still rely on.
///
/// Applied to **both** sides, so neither can hide a difference the other does not have: a
/// bucket with `num_obs > 0` survives on either side and still has to match, bases and
/// support and chain ids alike.
fn comparable(record: &PileupRecord) -> PileupRecord {
    let mut out = comparable_exact_q_sum(record);
    for allele in &mut out.alleles {
        allele.support.q_sum = round_q_sum(allele.support.q_sum);
    }
    out
}

/// Q_SUM_GRAIN: the granularity `q_sum` is compared at — see [`round_q_sum`].
const Q_SUM_GRAIN: f64 = 1e9;

/// Round `q_sum` to ~1e-9, because ng changed the **order** its accumulation happens in
/// and nothing else.
///
/// `q_sum` is an `f64` running sum, and `f64` addition is not associative. **Two changes
/// reorder it, and the count below moved sharply when the second landed** — 521 records at
/// A5, 5,368 at B1, over the same 20,000 cases:
///
/// - **A3's eviction.** Production keeps a bucket alive at `num_obs == 0` and keeps
///   accumulating into it, so a read that leaves and returns leaves `+q -q +q` behind; ng
///   evicts the empty bucket and recreates it, so the same read's sum starts from `0.0` and
///   is *exactly* `q`. Production's `-2.999999999999999` against ng's `-3.0` is this.
/// - **B1's per-read re-derivation.** Production's bucket total is accumulated during the
///   walk, with a subtract-then-add on every re-fold; ng's row sums each read's contribution
///   **once**, in `read_id` order. Same addends, different order — and ng's is the more
///   accurate of the two, since nothing cancels.
///
/// Neither is a difference in evidence: no read moved, no base changed, and the two numbers
/// differ in the last representable bits. **It is a named divergence class rather than an
/// absorbed one** — spec §3 lists five and warns that an unlisted one gets triaged as a
/// listed one and contaminates the measurement.
///
/// The grain is nine decimal places on values of order `-3` to `-50`, where the smallest
/// *real* difference is a whole read's `ln` contribution — order 1. So a genuine divergence
/// cannot hide under it, and `float_only_divergences` counts how often it fires so it can be
/// seen rather than assumed.
///
fn round_q_sum(q_sum: f64) -> f64 {
    (q_sum * Q_SUM_GRAIN).round() / Q_SUM_GRAIN
}

/// The projection itself, minus the `q_sum` rounding [`comparable`] adds — so a caller can
/// count how many records agree **only** because of the rounding, and it is shown to be
/// doing work rather than quietly matching nothing.
fn comparable_exact_q_sum(record: &PileupRecord) -> PileupRecord {
    let mut out = record.clone();
    let mut index = 0usize;
    out.alleles.retain(|allele| {
        let keep = index == 0 || allele.support.num_obs > 0;
        index += 1;
        keep
    });
    out.alleles[1..].sort_by(|a, b| a.seq.cmp(&b.seq));
    // **`placed_start` is zeroed on both sides — a third named projection, from B2.** ng
    // stops computing the quantity entirely (spec §6: nothing consumes it, and it is a pure
    // function of the read's start against the anchor), so its side is a structural zero and
    // every record would otherwise diverge on it. Zeroing *both* keeps the comparison
    // symmetric, which is the rule every projection here follows: neither side may be
    // normalised in a way the other is not. See `RecordEvidence` for why this is a
    // deliberate removal rather than the oversight it replaced.
    for allele in &mut out.alleles {
        allele.support.placed_start = 0;
    }
    out
}

/// Records that agree after the `q_sum` rounding and disagree before it.
fn float_only_divergences(
    ours: &[Result<PileupRecord, String>],
    theirs: &[Result<PileupRecord, String>],
) -> usize {
    ours.iter()
        .zip(theirs)
        .filter(|(ours, theirs)| match (ours, theirs) {
            (Ok(ours), Ok(theirs)) => {
                comparable_exact_q_sum(ours) != comparable_exact_q_sum(theirs)
                    && comparable(ours) == comparable(theirs)
            }
            _ => false,
        })
        .count()
}

/// Assert the two walks agree, and return whether they agreed *by both panicking*.
///
/// Shared by the synthetic differential and the real-data one, so the two cannot drift
/// into asserting different things about the same claim.
#[track_caller]
fn assert_same_walk(where_: &str, ours: &WalkOutcome, theirs: &WalkOutcome) -> bool {
    // **The panic comparison runs first**, and it compares the *message*. A verbatim copy
    // must panic verbatim: the claim is that the two reach the *same* precondition on the
    // same input, not merely that both stopped. Ordering matters for the diagnosis, not
    // the verdict — with the length check first, a case where ng panics early and
    // production does not reported "ng emitted 3 stream items, production 24", which names
    // the symptom and hides the cause.
    assert_eq!(
        ours.panic_message, theirs.panic_message,
        "{where_}: the two walkers did not stop the same way",
    );

    assert_eq!(
        ours.records.len(),
        theirs.records.len(),
        "{where_}: ng emitted {} stream items, production {}",
        ours.records.len(),
        theirs.records.len(),
    );
    for (position, (ours, theirs)) in ours.records.iter().zip(theirs.records.iter()).enumerate() {
        // Normalised on both sides — see `comparable` for why an empty bucket is a
        // projection rather than an excuse.
        let project = |item: &Result<PileupRecord, String>| {
            item.as_ref().map(comparable).map_err(String::clone)
        };
        assert_eq!(
            project(ours),
            project(theirs),
            "{where_}: stream item {position} diverged"
        );
    }

    if ours.panic_message.is_some() {
        return true;
    }

    assert_eq!(
        ours.summary
            .as_ref()
            .expect("a walk that did not panic has a summary"),
        theirs
            .summary
            .as_ref()
            .expect("a walk that did not panic has a summary"),
        "{where_}: the RunSummary counters diverged",
    );
    false
}

/// Cases per seed. Small enough to stay a unit test, large enough that every one of the
/// five behaviours B2 mutates fires many times over — which
/// `the_generator_exercises_what_the_port_can_break` asserts rather than assumes.
const CASES_PER_SEED: usize = 400;

/// Cases per seed, overridable by `PVC_PARITY_CASES` so a soak run is one command away —
/// the convention `delimit_parity` set.
///
/// **Use `--profile soak`, not `--release`** (Cargo.toml): `[profile.release]` leaves
/// `debug-assertions` off, so a release soak proves the two walkers diverge in the same
/// places and proves nothing about the invariants asserted along the way — which is most of
/// what this walk asserts. `soak` is release-speed with the assertions and overflow checks
/// armed:
///
/// ```text
/// PVC_PARITY_CASES=5000 cargo test --profile soak --lib ng::locus_generation::pileup::parity
/// ```
///
/// **Host-native.** `scripts/dev.sh` forwards only `CARGO_TARGET_DIR` and `HOME`, so
/// `PVC_PARITY_CASES` never reaches the container: a soak invoked through it silently walks
/// the default case count and finishes in under a second, looking like it worked.
fn cases_per_seed() -> usize {
    std::env::var("PVC_PARITY_CASES")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(CASES_PER_SEED)
}

/// Several seeds, because one generator stream can systematically miss a corner.
const SEEDS: [u64; 4] = [
    0x5EED_0001,
    0xC0FF_EE42,
    0x1234_5678_9ABC_DEF0,
    0xDEAD_BEEF_CAFE,
];

/// How one record's evidence relates to production's.
#[derive(Debug, PartialEq, Eq)]
enum RecordAgreement {
    /// Every supported bucket matches, bases and support and chain ids alike.
    Exact,
    /// The record's reference bytes and every support total match, and only the bases of
    /// some rows differ — production having widened past a read it did not re-fold. See
    /// [`classify_record`].
    EvidenceIntact,
    /// Anything else.
    Divergent,
}

/// **Every** support scalar of a whole record, summed across its buckets — "no evidence was
/// created, lost or double counted", independent of which bucket any of it ended up in.
///
/// # Why this is a struct built by an exhaustive destructure
///
/// It was a six-field tuple read field by field, and it **omitted `placed_start`** — the one
/// field this milestone made fragile, since A1 dropped it from ng's stats and `finalise`
/// reconstructs it from a per-read flag. That is not a hypothetical: injecting a real wrong
/// number (the REF bucket stops counting `placed_start`) left
/// [`ng_holds_the_same_evidence_as_production_on_complete_reads`] — the permanent anchor, the test
/// that *replaced* the retired stage-1 differential — **green**, while quietly moving
/// `EvidenceIntact` from 264 records (0.15 %) to 2,806 (1.58 %) and printing "same support
/// totals" and "Every other record is identical, field for field" as it went. 2,542 wrong
/// records absorbed. Only an inherited test from production's own suite noticed.
///
/// So the sum is taken by **destructuring production's [`AlleleSupportStats`]
/// exhaustively**: a field added to that type stops this file compiling instead of being
/// silently left out of the oracle. Milestone A's whole verification argument rests on this
/// function being total, and it was not.
///
/// `q_sum` is summed as `f64` and rounded **once, at the end**. Rounding each bucket first
/// would make the total depend on how the reads were distributed across buckets: two reads
/// at `-4.835428695` in one bucket sum to `-9.670857391`, while the same two in separate
/// buckets round to `-4835428695` twice and total `-9670857390`. A3 moves reads between
/// buckets by design, so that is exactly the difference this function must not see.
/// **`placed_start` is deliberately absent from this struct, and the difference between
/// that and its earlier absence is the whole point.**
///
/// Until B2 it was missing by *oversight*: the sum was a six-field tuple read field by
/// field, and a real wrong number in `placed_start` went unnoticed while 2,542 records were
/// absorbed into the tolerated class. From B2 ng does not compute the quantity **at all** —
/// spec §6: no model consumes it, it is a pure function of the read's start against the
/// anchor, and a later consumer re-derives it without touching the fold — so there is
/// nothing to compare, and comparing production's value against a structural zero would
/// fail every record.
///
/// The exhaustive destructure below is what keeps the two cases apart. A field added to
/// production's stats still stops this file compiling; `placed_start` is bound and dropped
/// **by name**, at one site, with this comment attached. An oversight cannot look like this.
#[derive(Debug, PartialEq, Eq)]
struct RecordEvidence {
    num_obs: u32,
    fwd: u32,
    placed_left: u32,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    q_sum_rounded: i64,
}

fn record_evidence(record: &PileupRecord) -> RecordEvidence {
    let mut totals = RecordEvidence {
        num_obs: 0,
        fwd: 0,
        placed_left: 0,
        mapq_sum: 0,
        mapq_sum_sq: 0,
        q_sum_rounded: 0,
    };
    let mut q_sum = 0.0_f64;
    for allele in &record.alleles {
        // The exhaustive destructure is the point — see the doc above.
        let AlleleSupportStats {
            num_obs,
            q_sum: bucket_q_sum,
            fwd,
            placed_left,
            // Bound and dropped by name — see `RecordEvidence`. ng stops computing this
            // at B2, so there is nothing on its side to compare against.
            placed_start: _,
            mapq_sum,
            mapq_sum_sq,
        } = allele.support;
        totals.num_obs += num_obs;
        totals.fwd += fwd;
        totals.placed_left += placed_left;
        totals.mapq_sum += mapq_sum;
        totals.mapq_sum_sq += mapq_sum_sq;
        q_sum += bucket_q_sum;
    }
    totals.q_sum_rounded = (q_sum * Q_SUM_GRAIN).round() as i64;
    totals
}

/// Every chain id a record carries, sorted and deduplicated — the *other* half of "no
/// evidence moved", and the half [`record_evidence`] cannot see because chain ids are a set
/// per bucket rather than a scalar.
///
/// **Equality here is the wrong test, and the fixture proved it.** Requiring the two sets to
/// match fails on `seed 0x5eed0001 case 11, record 30`, and correctly so: production folds a
/// read into the REF bucket — `num_obs: 2`, chain ids dropped by the `allele_index == 0`
/// rule — having missed an insertion it never re-folded, while ng emits the nine bases the
/// read actually witnessed as its own row, carrying chain id `6`. That is the defect being
/// fixed showing up in the ids, not evidence going missing.
///
/// So what is asserted is the **subset** direction, which is invariant: ng's REF bucket
/// holds the record's full reference bytes, so any read whose witness is partial or carries
/// an unseen event lands outside it and keeps its id. ng's REF bucket can therefore only
/// ever be a *subset* of production's, and ng's id set only ever a **superset** of
/// production's. An id production has that ng lacks means ng lost a read's identity, which
/// no part of this change is allowed to do.
fn record_chain_ids(record: &PileupRecord) -> Vec<ChainId> {
    let mut ids: Vec<ChainId> = record
        .alleles
        .iter()
        .flat_map(|allele| allele.chain_ids.iter().copied())
        .collect();
    ids.sort_unstable();
    ids.dedup();
    ids
}

/// **The one class of divergence A3 leaves on the complete-reads fixture — named, because
/// an unnamed class gets triaged as a named one** (spec §3).
///
/// After the owner's decision of 2026-07-29, `widen` re-folds every **live** read already
/// in the record, so a read sitting inside its own deletion no longer keeps a bucket — or,
/// more importantly, a `witnessed` extent — pinned to the pre-widen footprint. What
/// remains is the case production **cannot** get right by construction: a record that
/// widens past a read it does not re-fold, which after (b) means a read that has already
/// **expired**. Its cursor is gone, so neither walker can consult its events again;
/// production appends the new reference bases to its bucket anyway, and ng leaves the
/// bucket saying what the read saw.
///
/// The difference runs in **both** directions, which is why this cannot be stated as "ng's
/// bases are production's with a reference tail stripped":
///
/// - production credits the read with bases sequenced after it left the active set;
/// - production also *misses* events the read did witness inside the final footprint,
///   because it folded that read against a narrower window and never went back. An
///   insertion anchored in the widened region is the case that shows it — ng's row is then
///   one base **longer** than production's.
///
/// So what is checked is that **no evidence moved**: the record's own reference bytes are
/// identical, and every support scalar summed across the buckets is identical, so no read
/// was created, lost or double counted — only the bases of some rows differ.
///
/// # What this does not prove, and where it gets sharper
///
/// It does not check *which* rows differ or by how much. The right filter is the spec's
/// own definition of the anchor class — "loci where every folded read witnessed the whole
/// footprint" — which needs `FoldedReadState::witnessed` resolved against the final
/// footprint, i.e. **A4's `coverage_of`**. When that lands, this classifier should be
/// replaced by that predicate and the surviving class should be empty rather than merely
/// evidence-preserving. Recorded here so it is tightened rather than rediscovered (D1).
fn classify_record(ours: &PileupRecord, theirs: &PileupRecord) -> RecordAgreement {
    if comparable(ours) == comparable(theirs) {
        return RecordAgreement::Exact;
    }
    // **The totals come from the originals, not from `comparable`'s output.** That
    // projection rounds each bucket's `q_sum` to 1e-9, and a record's fourteen buckets
    // then carry up to fourteen half-grains of slack — more than the grain itself, so a
    // record whose reads merely sit in different buckets would fail this check on the
    // rounding alone. Summing first and rounding once is the whole point.
    //
    // **The chain ids are checked too, by subset and not by equality** — see
    // `record_chain_ids` for the record that settles which of the two is right. ng may hold
    // an id production dropped into its REF bucket; it may never *lose* one.
    let ng_ids = record_chain_ids(ours);
    if ours.alleles[0].seq == theirs.alleles[0].seq
        && record_evidence(ours) == record_evidence(theirs)
        && record_chain_ids(theirs)
            .iter()
            .all(|id| ng_ids.binary_search(id).is_ok())
    {
        return RecordAgreement::EvidenceIntact;
    }
    RecordAgreement::Divergent
}

/// Set in the child process by the determinism test below, to select its other half.
const DETERMINISM_CHILD_VAR: &str = "PVC_DETERMINISM_CHILD";

/// Set alongside it, so an inherited child marker cannot silently disarm the test — see the
/// guard in the test body.
const DETERMINISM_PARENT_VAR: &str = "PVC_DETERMINISM_PARENT";

/// The name the determinism test spawns itself under. A `const` beside the test rather
/// than a literal at the call site, because a rename that missed one of the two would make
/// the child run *no* test, print no digest, and fail with "the child printed no digest" —
/// which reads like a bug in the walk rather than a stale string.
const DETERMINISM_TEST_PATH: &str =
    "ng::locus_generation::pileup::parity::ng_emits_the_same_bytes_in_a_second_process";

/// Walk a fixed set of cases and reduce the whole emitted stream to one digest.
///
/// Over ng's **own** type, not the `PileupRecord` projection: the claim under test is about
/// what the generator emits, and the projection merges rows back together, which is exactly
/// the axis a hash-order bug would show up on.
fn determinism_digest() -> String {
    determinism_digest_with(|_| {})
}

/// The same, with a hook that perturbs each case's reads — used only by the positive
/// control, which needs the digest to be shown *sensitive* rather than merely stable.
fn determinism_digest_with(perturb: impl Fn(&mut Vec<PreparedRead>)) -> String {
    use std::hash::{Hash, Hasher};

    let mut rng = SplitMix64(0xD37E_2E1D);
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    let mut records = 0usize;
    for _ in 0..32 {
        let case = generate(&mut rng);
        let fasta = case.fasta();
        let reads: Vec<PreparedRead> = case
            .reads
            .iter()
            .cloned()
            .map(|read| PreparedRead::from_production(read, PLACEHOLDER_READ_GROUP))
            .collect();
        let mut reads = reads;
        perturb(&mut reads);
        for item in super::run(reads, fasta, &case.config) {
            match item {
                Ok(locus) => {
                    records += 1;
                    format!("{locus:?}").hash(&mut hasher);
                }
                Err(error) => render_ng_error(&error).hash(&mut hasher),
            }
        }
    }
    format!("{records}:{:016x}", hasher.finish())
}

/// **The same input walked in two separate processes emits the same bytes** — spec §7's
/// "output is a deterministic function of (reference, config, reads)", which until now was
/// claimed and never tested.
///
/// # Why it has to be two processes
///
/// `folded_reads` is an `AHashMap`, and `ahash` seeds itself **per process**. Inside one
/// process that order is arbitrary but *fixed*, so a hash-order dependency is invisible to
/// any number of runs in the same binary — and every other test in this file compares ng
/// against production in that same process, never one ng run against another. Two
/// mechanisms in `refold_live_reads` are unpinnable for exactly this reason; the owner's
/// decision (2026-07-29) was that B2's sort is the guarantee and that this is where it gets
/// proven, rather than writing cross-process tests for mechanisms the sort makes redundant.
///
/// # The canary, which is what stops this passing vacuously
///
/// If `ahash` ever stopped randomising, the digests would match for a reason that has
/// nothing to do with the sort, and this test would go on passing while proving nothing —
/// the failure mode this branch has hit repeatedly. So each child also prints a hash of a
/// fixed string under a fresh `RandomState`, and the parent asserts those **differ**. The
/// two assertions together say: the seed really did change, and the output really did not.
#[test]
fn ng_emits_the_same_bytes_in_a_second_process() {
    // **A child marker inherited from the environment would make this run take the child
    // branch and assert nothing, while still reporting `ok`.** The parent sets the variable
    // itself, so seeing it here at top level means it came from outside — which is a broken
    // invocation, not a child.
    if std::env::var(DETERMINISM_CHILD_VAR).is_ok() && std::env::var(DETERMINISM_PARENT_VAR).is_ok()
    {
        println!("DIGEST {}", determinism_digest());
        println!(
            "CANARY {:016x}",
            ahash::RandomState::new().hash_one("canary")
        );
        return;
    }

    let exe = std::env::current_exe().expect("a test binary knows its own path");
    let run_child = || -> (String, String) {
        let output = std::process::Command::new(&exe)
            .args(["--exact", DETERMINISM_TEST_PATH, "--nocapture"])
            .env(DETERMINISM_CHILD_VAR, "1")
            .env(DETERMINISM_PARENT_VAR, "1")
            .output()
            .expect("the child test binary runs");
        let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
        assert!(
            output.status.success(),
            "the child walk failed:\n{stdout}\n{}",
            String::from_utf8_lossy(&output.stderr),
        );
        let field = |prefix: &str| {
            stdout
                .lines()
                .find_map(|line| line.strip_prefix(prefix))
                .unwrap_or_else(|| {
                    panic!(
                        "the child printed no {prefix}line — has {DETERMINISM_TEST_PATH} \
                         been renamed?\n{stdout}"
                    )
                })
                .to_owned()
        };
        (field("DIGEST "), field("CANARY "))
    };

    let (first_digest, first_canary) = run_child();
    let (second_digest, second_canary) = run_child();

    assert_ne!(
        first_canary, second_canary,
        "the two child processes hashed a fixed string to the same value, so `ahash` is \
         not seeding per process here and this test cannot detect a hash-order dependency \
         at all — it would pass whatever the walk did"
    );
    assert!(
        first_digest.split(':').next().is_some_and(|records| {
            records.parse::<usize>().is_ok_and(|records| records > 1000)
        }),
        "only {first_digest} — too few records for the fixture to exercise the fold's \
         hash-order surface"
    );
    assert_eq!(
        first_digest, second_digest,
        "the same reads, reference and config emitted different bytes in two processes — \
         something in the emission depends on the per-process hash seed, which spec §7 \
         forbids outright"
    );
}

/// **The positive control for the test above** — without it, a digest that stopped hashing
/// the loci would compare equal across processes for a reason that has nothing to do with
/// the walk, and a real hash-order bug could sit underneath it undetected.
///
/// `records > 1000` does not cover this: it counts loci, not sensitivity. Nor does removing
/// a read — that changes the record *count*, so a digest hashing nothing but the count still
/// moves and the control passes while proving nothing. (Checked: it does.) The perturbation
/// has to change **what the records say** while leaving which records exist alone, so this
/// one rewrites every read's MAPQ: same loci, same bases, different `mapq_sum`.
#[test]
fn the_determinism_digest_responds_to_the_evidence() {
    let full = determinism_digest();
    let remapped = determinism_digest_with(|reads| {
        for read in reads.iter_mut() {
            read.mapq = read.mapq.wrapping_add(7).max(1);
        }
    });
    assert_ne!(
        full, remapped,
        "every read's MAPQ changed and the digest did not, so it is not a function of what \
         the walk emitted — `ng_emits_the_same_bytes_in_a_second_process` would compare two \
         constants and pass whatever the walk did"
    );
    assert_eq!(
        full.split(':').next(),
        remapped.split(':').next(),
        "MAPQ must not change which records exist, or this control is testing the record \
         count rather than the evidence"
    );
}

/// **The permanent anchor: on reads that witnessed the whole footprint, the two walkers
/// agree forever.**
///
/// From A2 the two walkers differ *on purpose* — production fills the positions a read did
/// not witness from the reference, and ng does not — so the full stage-1 differential
/// stops being a claim anyone can make. **What survives is narrower and is not a
/// snapshot:** every change in this plan leaves the complete class alone by construction
/// (there are no gaps to fill), so this must stay green through A3, A4, A5 and Milestone
/// B, and it is what replaces the differential rather than merely outliving it.
///
/// The fixture is `generate_complete`, whose reads span their contig end to end and are
/// silent nowhere inside it (`generate_complete_read` names the four exclusions and why
/// each one matters). Everything else the walk does is still on it: substitutions,
/// insertions, deletions that widen records past the reads that opened them, re-folds,
/// mate overlap, both strands, the column cap.
///
/// A failure means ng has changed an answer on the class it promised not to touch; the
/// seed and case index replay the exact input.
#[test]
fn ng_holds_the_same_evidence_as_production_on_complete_reads() {
    let mut compared_records = 0usize;
    let mut widened_records = 0usize;
    let mut widen_stale = 0usize;
    let mut float_only = 0usize;
    let total = SEEDS.len() * cases_per_seed();

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate_complete(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");

            let theirs = production_walk(&case);
            let ours = ng_walk(&case);
            compared_records += ours.records.len();
            widened_records += theirs
                .records
                .iter()
                .filter(|item| item.as_ref().is_ok_and(|r| r.alleles[0].seq.len() > 1))
                .count();
            assert_eq!(
                ours.panic_message, None,
                "{where_}: a complete-reads case must not panic — this fixture excludes \
                 every input that reaches production's reachable precondition"
            );
            assert_eq!(
                ours.panic_message, theirs.panic_message,
                "{where_}: the two walkers did not stop the same way",
            );
            assert_eq!(
                ours.records.len(),
                theirs.records.len(),
                "{where_}: ng emitted {} stream items, production {}",
                ours.records.len(),
                theirs.records.len(),
            );
            for (position, (ours, theirs)) in
                ours.records.iter().zip(theirs.records.iter()).enumerate()
            {
                match (ours, theirs) {
                    (Ok(ours), Ok(theirs)) => match classify_record(ours, theirs) {
                        RecordAgreement::Exact => {}
                        RecordAgreement::EvidenceIntact => widen_stale += 1,
                        RecordAgreement::Divergent => panic!(
                            "{where_}: record {position} diverged on a complete-reads \
                             case, and not as production's widen appending reference \
                             bases to a bucket belonging to an already-expired read — \
                             which is the only class this fixture is allowed to produce.\n  ng   {:?}\n  \
                             prod {:?}",
                            comparable(ours),
                            comparable(theirs),
                        ),
                    },
                    (ours, theirs) => assert_eq!(
                        ours, theirs,
                        "{where_}: stream item {position} diverged on the error channel"
                    ),
                }
            }
            assert_eq!(
                ours.summary.as_ref(),
                theirs.summary.as_ref(),
                "{where_}: the RunSummary counters diverged",
            );
            float_only += float_only_divergences(&ours.records, &theirs.records);
        }
    }

    // The anchor is only worth the ground it covered, and "multi-base records exist" is
    // the part that matters most: a one-base record is `Complete` for every contributor
    // whatever the builder does, so a fixture of only those would hold this test green
    // against any implementation at all.
    assert!(
        compared_records > total * 5,
        "only {compared_records} records compared over {total} cases — the generator has \
         stopped producing walks worth comparing"
    );
    assert!(
        widened_records > total,
        "only {widened_records} multi-base records over {total} cases — without them this \
         anchor cannot tell a filling builder from a witnessing one"
    );
    // The class has to be **present**, or the classifier is a branch nothing takes and a
    // record that started diverging some other way would be reported as this one.
    assert!(
        widen_stale > 0,
        "no record showed production widening past a read it never re-folded — either the \
         generator stopped producing widens, or `classify_record` is matching something \
         it should not"
    );
    assert!(
        widen_stale * 20 < compared_records,
        "{widen_stale} of {compared_records} records fell outside byte-identity — that \
         class is supposed to be the rare corner where a record widens past an expired \
         read, and at this rate it is something else"
    );
    eprintln!(
        "complete-reads differential: {compared_records} records compared over {total} \
         cases, {widened_records} of them multi-base. {widen_stale} ({:.2}%) hold the \
         same evidence — same reference bytes, same support totals — with some rows' \
         bases differing, because production widened past a read it never re-folded. \
         Every other record is identical, field for field. {float_only} agree only after \
         `q_sum` is rounded to 1e-9 — the order the sum accumulates in, from two causes: \
         A3's eviction recreating a bucket the reads return to, and B1 summing each read's \
         contribution once where production accumulates into the bucket with a \
         subtract-then-add per re-fold.",
        100.0 * widen_stale as f64 / compared_records as f64,
    );
}

/// **What A2 is allowed to change, and what it is not** — the comparison both the
/// synthetic census and the real-data one run.
///
/// Records are opened and widened from the **events**, which the no-fabrication rule does
/// not touch; only the allele buckets move. So the two streams must still be the same
/// length, at the same anchors, with the same REF bytes, the same error items and the same
/// `RunSummary`. Anything else is a bug rather than a design.
///
/// Returns `(records compared, records whose allele lists differ)` — the census.
#[track_caller]
fn assert_only_allele_bytes_moved(
    where_: &str,
    ours: &WalkOutcome,
    theirs: &WalkOutcome,
) -> (usize, usize) {
    assert_eq!(
        ours.panic_message, theirs.panic_message,
        "{where_}: the two walkers did not stop the same way",
    );
    assert_eq!(
        ours.records.len(),
        theirs.records.len(),
        "{where_}: ng emitted {} stream items, production {} — A2 moves allele bytes, \
         never whether a record exists",
        ours.records.len(),
        theirs.records.len(),
    );

    let mut records = 0usize;
    let mut diverged = 0usize;
    for (position, (ours, theirs)) in ours.records.iter().zip(theirs.records.iter()).enumerate() {
        records += 1;
        match (ours, theirs) {
            (Ok(ours), Ok(theirs)) => {
                assert_eq!(
                    (ours.chrom_id, ours.pos, &ours.alleles[0].seq),
                    (theirs.chrom_id, theirs.pos, &theirs.alleles[0].seq),
                    "{where_}: record {position} moved its anchor or its REF bytes, which \
                     A2 does not touch",
                );
                if comparable(ours) != comparable(theirs) {
                    diverged += 1;
                }
            }
            (ours, theirs) => assert_eq!(
                ours, theirs,
                "{where_}: stream item {position} diverged on the error channel, which A2 \
                 does not touch",
            ),
        }
    }

    if ours.panic_message.is_none() {
        assert_eq!(
            ours.summary.as_ref(),
            theirs.summary.as_ref(),
            "{where_}: the RunSummary counters diverged — every one of them is driven by \
             the events, which A2 does not touch",
        );
    }
    (records, diverged)
}

/// **The fabrication is gone, and here is how much of it there was.**
///
/// The general fixture — partial reads, adaptor boundaries, `N` bases, ref-skips — is
/// where the two walkers now differ, and this is the census of that difference rather than
/// an assertion that it is absent. Two things are asserted, and they are what makes the
/// census trustworthy:
///
/// - **A2 changes no record's existence and no record's footprint.** Records are opened
///   and widened from the *events*, which this step does not touch; only the allele
///   buckets move. So the two streams must still be the same length, at the same
///   positions, with the same REF bytes and the same `RunSummary`. A divergence there is a
///   bug, not a design.
/// - **Some records must differ.** A run in which nothing moved would mean the fill is
///   still there and every other test in this module is passing for the wrong reason.
///
/// The number reported here is the same quantity D3 measures on real data: how many loci
/// production credits to reads that never sequenced them.
#[test]
fn ng_diverges_from_production_only_where_a_read_did_not_witness() {
    let mut records = 0usize;
    let mut diverged = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");

            let theirs = production_walk(&case);
            let ours = ng_walk(&case);
            let (seen, differing) = assert_only_allele_bytes_moved(&where_, &ours, &theirs);
            records += seen;
            diverged += differing;
        }
    }

    assert!(
        diverged > 0,
        "{records} records compared and not one differed — the fill is still there, and \
         every complete-reads assertion in this module is passing for the wrong reason"
    );
    // **A ceiling as well as a floor.** The floor alone makes this a census that cannot
    // fail upwards: the headline can be driven from 1.9% to 91.8% by a genuine defect and
    // the test still passes, reporting the defect as the measurement. The class is
    // "production fabricated bases here", which is a minority of loci by construction —
    // most positions carry no widened record at all — so a figure in the tens of percent
    // means the walk changed, not that the defect is larger than thought.
    assert!(
        diverged * 10 < records,
        "{diverged} of {records} records ({:.1}%) differ — this census measures production's \
         fabrication, which is a small minority of loci; at this rate something in the walk \
         has moved and the number is no longer the defect's size",
        100.0 * diverged as f64 / records as f64,
    );
    eprintln!(
        "the fabrication census: {diverged} of {records} records ({:.1}%) carried bases \
         production credited to a read that had not witnessed them",
        100.0 * diverged as f64 / records as f64,
    );
}

/// **No record leaves the walker carrying a bucket no read is folded into** — A3's
/// eviction, asserted on the records that are *emitted* rather than on the table inside.
///
/// Nothing checked this, and the gap was structural rather than accidental: A3's own
/// eviction fixtures reach into `OpenPileupRecordTable` while the record is still open, and
/// the differential is **blind by construction**, because [`comparable`] drops unsupported
/// non-REF buckets on *both* sides before comparing — precisely so that ng evicting them and
/// production keeping them is not read as a divergence. So the one projection that makes the
/// two walkers comparable is also the one that hides whether ng evicted anything at all.
///
/// The concrete hole: moving `evict_unsupported_alleles` to *before* the contributor fold
/// loop leaves the whole suite green. Buckets emptied by the widen are still caught, because
/// the widen runs first — but buckets emptied by the **fold loop itself**, when a contributor
/// re-folds into a different bucket at this position, survive to `finalise` and are emitted.
/// That is the accumulation spec §7 says to design for rather than discover, against a
/// `find_allele_index` that is a linear scan with a full byte compare.
///
/// `alleles[0]` is exempt: it is the REF sequence, and production creates it with zero
/// observations by design.
#[test]
fn ng_emits_no_allele_bucket_without_support() {
    let mut records = 0usize;
    let mut multi_allele = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");
            for (position, item) in ng_walk(&case).records.iter().enumerate() {
                let Ok(record) = item else { continue };
                records += 1;
                if record.alleles.len() > 2 {
                    multi_allele += 1;
                }
                for (index, allele) in record.alleles.iter().enumerate().skip(1) {
                    assert!(
                        allele.support.num_obs > 0,
                        "{where_}: emitted record {position} carries bucket {index} \
                         ({:?}) with no supporting read — a re-fold left it behind and the \
                         eviction did not run after it: {record:?}",
                        String::from_utf8_lossy(&allele.seq),
                    );
                }
            }
        }
    }

    // Without a record that has somewhere to move a read *to*, the property is vacuous —
    // a walk of one-allele records satisfies it against an implementation that never
    // evicts anything.
    assert!(
        multi_allele * 100 > records,
        "only {multi_allele} of {records} emitted records carry more than one non-REF \
         bucket — this fixture cannot exercise a read moving between buckets, so the \
         property it asserts is vacuous"
    );
    eprintln!(
        "eviction census: {records} emitted records, {multi_allele} of them with more than \
         one non-REF bucket; none carried an unsupported one"
    );
}

/// **The regression test for the production defect this differential found** (owner:
/// *"we should fix the production bug"*, 2026-07-29).
///
/// `events_overlapping` **does not clip a deletion to the window** — spec §8 records that,
/// and a deletion anchored before a record whose run reaches into it comes back whole. What
/// nobody had recorded is that `apply_events_to_ref_into` then took
/// `offset = anchor.saturating_sub(record_pos)`, and the saturation was **silently wrong
/// twice**: it emitted `ref_seq[0]` — a base this read had *deleted* — and skipped
/// `offset + 1 + deleted_len`, one position too many. In debug a `debug_assert!` caught it;
/// in release it produced wrong allele bytes with no error, at exactly the long-deletion
/// loci this port exists to get right.
///
/// The fixture is the three-read case the differential's generator produced, shrunk:
///
/// 1. `pair3/Second` carries a deletion anchored at **17**, spanning 18–22.
/// 2. At 17 it overlaps its own mate. Mate-overlap reconciliation in the indel regime
///    **collapses the pair to a single observation**, removing the contributor that carried
///    the indel — so **no record opens at 17**.
/// 3. `pair9/First`'s deletion, anchored at **19**, opens a record there — *inside* the
///    footprint of a deletion that never opened one.
/// 4. Where `pair3/Second` matches again (23+), it folds into that record, and the fold is
///    handed its deletion anchored at 17, two positions before `record_pos`.
///
/// The fix computes the skip from **absolute coordinates** —
/// `anchor + 1 + deleted_len - record_pos` — and emits the anchor base only when the anchor
/// is inside the record. For a deletion anchored *within* its record that is arithmetically
/// identical to the old expression, which is why it changed nothing else.
///
/// **This test is the parity claim and the regression together:** the two walkers must
/// agree, and the record must show `pair3/Second` contributing none of the bases it
/// deleted.
#[test]
fn a_deletion_anchored_before_its_record_contributes_none_of_the_bases_it_deleted() {
    fn read(
        qname: &str,
        role: ProductionMateRole,
        start: u32,
        end: u32,
        cigar: Vec<CigarOp>,
        seq_len: usize,
    ) -> ProductionPreparedRead {
        ProductionPreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: end,
            cigar,
            seq: vec![b'A'; seq_len],
            bq_baq: vec![30; seq_len],
            mq_log_err: -3.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from(qname),
            mate_role: role,
            adaptor_boundary: None,
        }
    }

    let case = Case {
        reference: vec!["ACGT".repeat(40), "ACGT".repeat(40)],
        reads: vec![
            // The mate that reconciliation keeps, suppressing the record at 17.
            read(
                "pair3",
                ProductionMateRole::FirstOfPair,
                4,
                22,
                vec![CigarOp::Match(6), CigarOp::Match(8), CigarOp::Match(5)],
                19,
            ),
            // The deletion anchored at 17, spanning 18–22, that opens no record.
            read(
                "pair3",
                ProductionMateRole::SecondOfPair,
                8,
                38,
                vec![
                    CigarOp::Match(8),
                    CigarOp::Match(2),
                    CigarOp::Deletion(5),
                    CigarOp::Match(8),
                    CigarOp::Match(8),
                ],
                26,
            ),
            // The deletion anchored at 19 that opens the record the first one folds into.
            read(
                "pair9",
                ProductionMateRole::FirstOfPair,
                10,
                42,
                vec![
                    CigarOp::Match(2),
                    CigarOp::Match(8),
                    CigarOp::Deletion(6),
                    CigarOp::Match(9),
                    CigarOp::Deletion(6),
                    CigarOp::Match(2),
                ],
                21,
            ),
        ],
        config: WalkerConfig::default(),
        reads_with_live_adaptor_boundary: 0,
    };

    let theirs = production_walk(&case);
    let ours = ng_walk(&case);

    assert_eq!(
        theirs.panic_message, None,
        "the input that used to trip production's debug_assert! must now walk cleanly"
    );
    assert_eq!(
        ours.panic_message, theirs.panic_message,
        "ng's copy carries the same fix, so it must stop — or not — the same way"
    );

    // The record at 19 is the one the defect corrupted: `pair9/First`'s deletion opens it
    // spanning 19..=25, and `pair3/Second` folds in carrying a deletion anchored at 17
    // whose run covers 18–22 — so of this record's seven positions it witnessed only
    // 23, 24, 25.
    let at_19 = |outcome: &WalkOutcome| {
        outcome
            .records
            .iter()
            .filter_map(|item| item.as_ref().ok())
            .find(|record| record.pos == 19)
            .map(|record| {
                record
                    .alleles
                    .iter()
                    .map(|allele| String::from_utf8_lossy(&allele.seq).to_string())
                    .collect::<Vec<_>>()
            })
            .expect("pair9's deletion opens a record at 19")
    };
    let (theirs_at_19, ours_at_19) = (at_19(&theirs), at_19(&ours));
    assert_eq!(
        theirs_at_19[0].len(),
        7,
        "the record spans 19..=25: the anchor plus six deleted bases"
    );

    // **The wrong answer, named.** Every fixture read's `seq` is all `A`, so
    // `pair3/Second`'s honest contribution to this record is its own three bases at 23, 24
    // and 25 — the positions past the end of its own deletion — and nothing before them.
    //
    // Before `5f32a62` the saturated offset made the fold emit `ref_seq[0]` first: the base
    // at **19**, which this read had deleted and never sequenced. So the allele was one
    // base longer and began with a reference base borrowed from a position the read
    // explicitly says is absent. Both spellings are checked, because "the right bases are
    // present" and "the wrong bases are gone" are different claims and a regression could
    // satisfy either alone. **Asserted on both walkers**: this is production's fix as much
    // as ng's, and it is the one claim in this test that A2 does not move.
    let before_the_fix = format!("{}AAA", &theirs_at_19[0][..1]);
    for (whose, alleles) in [("production", &theirs_at_19), ("ng", &ours_at_19)] {
        assert!(
            alleles.iter().any(|allele| allele == "AAA"),
            "{whose}: the read whose deletion covers 19–22 should contribute the three \
             bases it witnessed (AAA), but the record holds {alleles:?}",
        );
        assert!(
            !alleles.contains(&before_the_fix),
            "{whose}: the record still holds {before_the_fix}, the pre-fix spelling — a \
             leading base at 19 that this read had deleted. Records: {alleles:?}",
        );
    }

    // **And this is the witnessed-extent rule, on the fixture the defect handed it.**
    // `pair3/First` matched 4–22, so of the record's seven positions it witnessed four:
    // 19, 20, 21, 22. Production folded it as `AAAA` **plus the reference bases at 23, 24
    // and 25** — three bases it never sequenced. ng emits the four it saw.
    assert!(
        theirs_at_19.iter().any(|allele| allele == "AAAAGTA"),
        "production should fabricate the tail here, or this fixture no longer shows the \
         difference: {theirs_at_19:?}",
    );
    assert!(
        ours_at_19.iter().any(|allele| allele == "AAAA"),
        "ng should carry only the four positions the read witnessed: {ours_at_19:?}",
    );
    assert!(
        !ours_at_19.iter().any(|allele| allele == "AAAAGTA"),
        "ng still carries production's reference-padded spelling: {ours_at_19:?}",
    );
    // Everything else about the record is untouched: same REF bytes, same anchor.
    assert_eq!(
        theirs_at_19[0], ours_at_19[0],
        "the REF bucket is the record's own reference bytes and A2 does not touch it",
    );
}

/// **The differential is worth exactly what its generator reaches, and this is what says
/// what it reaches.**
///
/// A parity test that passes on inputs exercising none of the interesting paths is a claim
/// dressed as evidence — the failure mode this branch has hit in four consecutive
/// milestones. So the five behaviours plan 2's B2 mutates are counted here, over the same
/// generator, and each is required to fire. If a change to the generator stops producing
/// deep columns or overlapping mates, this fails immediately and loudly rather than
/// leaving `ng_holds_the_same_evidence_as_production_on_complete_reads` quietly weaker.
///
/// The counts come from production's own `RunSummary` wherever it has one, so they measure
/// what the *walker* did rather than what the generator intended. Adaptor boundaries have
/// no counter, so they are tallied at generation.
#[test]
fn the_generator_exercises_what_the_port_can_break() {
    let mut widens = 0u64;
    let mut mate_overlaps = 0u64;
    let mut cap_truncations = 0u64;
    let mut adaptor_boundaries = 0usize;
    let mut chain_allocations = 0u64;
    let mut multi_base_records = 0usize;
    let mut multi_allele_records = 0usize;
    let mut mate_evictions = 0u64;
    let mut errors = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        // Scales with `PVC_PARITY_CASES` like the differential does: a soak that widened
        // the comparison but not the coverage floor would report more confidence in the
        // same evidence.
        for _ in 0..cases_per_seed() {
            let case = generate(&mut rng);
            adaptor_boundaries += case.reads_with_live_adaptor_boundary;
            let outcome = production_walk(&case);
            if let Some(summary) = &outcome.summary {
                widens += summary.record_widen_events;
                mate_overlaps += summary.mate_overlap_positions;
                chain_allocations += summary.chain_allocations;
                cap_truncations += summary.column_depth_truncations;
                mate_evictions += summary.mate_lookup_evictions;
            }
            for item in &outcome.records {
                match item {
                    Ok(record) => {
                        if record.alleles[0].seq.len() > 1 {
                            multi_base_records += 1;
                        }
                        if record.alleles.len() > 1 {
                            multi_allele_records += 1;
                        }
                    }
                    Err(_) => errors += 1,
                }
            }
        }
    }

    eprintln!(
        "generator coverage: {widens} widens, {mate_overlaps} mate-overlap positions, \
         {cap_truncations} cap truncations, {adaptor_boundaries} live adaptor boundaries, \
         {chain_allocations} chain allocations, {mate_evictions} mate evictions, \
         {multi_base_records} multi-base records, {multi_allele_records} multi-allele \
         records, {errors} walker errors"
    );

    // Each of these is a behaviour B2 mutates. A zero here means that mutation could not
    // fail the differential, whatever B2 reports.
    assert!(widens > 0, "no record ever widened");
    assert!(mate_overlaps > 0, "mate-overlap reconciliation never fired");
    assert!(cap_truncations > 0, "the column depth cap never fired");
    assert!(
        adaptor_boundaries > 0,
        "no read carried a live adaptor boundary"
    );
    assert!(chain_allocations > 0, "no chain id was ever allocated");
    // The eighth compared counter. Without this it is compared but never shown to be
    // non-trivially populated — and it is the path where a pair silently degrades to two
    // solos, which is a different chain id and a different fold.
    assert!(
        mate_evictions > 0,
        "the mate-lookup window never evicted, so that counter's comparison is vacuous"
    );

    // A widened record with several alleles is where the subtract-then-add re-fold runs:
    // a live read re-folds against the wider window and must move between buckets exactly
    // once, not once per position of the footprint.
    assert!(
        multi_base_records > 0,
        "no record footprint ever exceeded one base, so no re-fold ran"
    );
    assert!(
        multi_allele_records > 0,
        "every record carried a single allele, so nothing distinguishes the buckets"
    );

    // Not a behaviour to exercise — a diagnostic. The generator is built to keep reads
    // in bounds and in order, so a large error count means it has drifted into testing
    // the walker's rejection paths instead of its fold.
    assert!(
        errors * 20 < multi_allele_records.max(1),
        "the generator is producing too many walker errors ({errors}) to be testing the fold"
    );
}

/// **The `Err` half of the stream, which the main generator deliberately never produces.**
///
/// `ng_holds_the_same_evidence_as_production_on_complete_reads` is built to keep reads in bounds and in order — its
/// own coverage test reports zero walker errors — so `assert_same_walk`'s element-wise
/// comparison only ever sees `Ok`, and the `map_err` machinery that exists to compare two
/// nominally distinct `WalkerError` types is dead. Spec §3 states the claim as "the two
/// `Result<PileupRecord, WalkerError>` streams are equal element for element"; without this
/// test, half of that has no input behind it.
///
/// Every `WalkerError` variant is fatal and terminal for the iterator, so these cases
/// cannot live in the main generator without truncating the walks that test the fold. They
/// get their own fixtures, and each fixture is required to actually **reach** its error —
/// otherwise this compares two clean walks and calls it error parity, which is the exact
/// failure it exists to close.
#[test]
fn both_walkers_report_the_same_error_on_the_same_malformed_input() {
    fn read(
        qname: &str,
        start: u32,
        end: u32,
        cigar: Vec<CigarOp>,
        seq_len: usize,
    ) -> ProductionPreparedRead {
        ProductionPreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: end,
            cigar,
            seq: vec![b'A'; seq_len],
            bq_baq: vec![30; seq_len],
            mq_log_err: -3.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from(qname),
            mate_role: ProductionMateRole::Solo,
            adaptor_boundary: None,
        }
    }

    let reference = vec!["ACGT".repeat(40), "ACGT".repeat(40)];
    // `(name, reads, the variant the fixture must reach)`. The third element is what
    // makes each fixture prove its *own* claim: without it a fixture that started failing
    // for some unrelated reason would still satisfy "production emitted an error", and
    // the case would go on passing while testing something else.
    let fixtures: Vec<(&str, Vec<ProductionPreparedRead>, &str)> = vec![
        (
            "out of order",
            vec![
                read("a", 20, 27, vec![CigarOp::Match(8)], 8),
                read("b", 4, 11, vec![CigarOp::Match(8)], 8),
            ],
            "OutOfOrder",
        ),
        (
            // The check is `alignment_end < alignment_start`, not "the CIGAR consumes no
            // reference" — an all-insertion read whose `alignment_end` equals its start
            // sails through, which the fixture's own reach assertion caught.
            "zero reference span",
            vec![read("i", 4, 3, vec![CigarOp::Insertion(4)], 4)],
            "ZeroRefSpan",
        ),
        (
            "cigar consumes more read bases than seq provides",
            vec![read("m", 4, 11, vec![CigarOp::Match(8)], 5)],
            "MalformedRead",
        ),
        (
            "seq and bq of different lengths",
            vec![{
                let mut malformed = read("q", 4, 11, vec![CigarOp::Match(8)], 8);
                malformed.bq_baq.truncate(7);
                malformed
            }],
            "MalformedRead",
        ),
        (
            // The reference is 160 bases and this read's footprint runs past its end, so
            // `open_new` asks for bases that do not exist. **This is the fixture that pays
            // for A0's error normalisation:** `WalkerError::Fasta` is the one variant whose
            // source type differs between the two enums, so without a case reaching it
            // `render_production_error` and `render_ng_error` would be two pieces of dead
            // code agreeing with each other. The main generator clamps every read inside
            // its contig precisely to avoid this path, so it has to be built by hand — and
            // the marker names the normalised *cause*, not just the variant, since it is
            // the cause that the two sides state in different types.
            //
            // The failing fetch is at **161**, not at the read's start: a record is opened
            // per covered position, so the walk fails at the first position past the
            // contig rather than at the first position of the offending read.
            "fetch past the contig end",
            vec![read("far", 156, 175, vec![CigarOp::Match(20)], 20)],
            "Fasta { chrom_id: 0, start: 161, start_plus_len: 162, cause: OutOfBounds { \
             contig_length: 160, start: 161, end: 162 } }",
        ),
    ];

    for (name, reads, expected) in fixtures {
        let case = Case {
            reference: reference.clone(),
            reads,
            config: WalkerConfig::default(),
            reads_with_live_adaptor_boundary: 0,
        };
        let theirs = production_walk(&case);
        let ours = ng_walk(&case);
        let rendered: Vec<&String> = theirs
            .records
            .iter()
            .filter_map(|item| item.as_ref().err())
            .collect();
        assert!(
            rendered.iter().any(|error| error.contains(expected)),
            "{name}: production reached no `{expected}`, so this fixture tests something \
             else — it emitted {rendered:?}"
        );
        assert_same_walk(name, &ours, &theirs);
    }
}

/// **The two sides of the differential are handed the same reference bytes.**
///
/// A0 split how they reach it: production's walker takes `MockFasta` through
/// `MultiChromRefFetcher`, ng's takes the same value through the canonicalising
/// [`RefSeq`](crate::ng::ref_seq::RefSeq) view in
/// [`mock_reference`](super::mock_reference). Those two agree only because the generator
/// draws its reference from `ACGTN`, where canonicalisation is the identity — a fact
/// about the *fixture*, not about either impl.
///
/// So it is checked rather than relied on. Without this, a generator change introducing a
/// lower-case or ambiguity-coded base would surface as `ng_holds_the_same_evidence_as_production_on_complete_reads`
/// failing on the bases inside a record, and be chased into the walk — where there is
/// nothing to find.
#[test]
fn both_sides_of_the_differential_are_served_the_same_bytes() {
    use crate::fasta::MultiChromRefFetcher;
    use crate::ng::ref_seq::RefSeq;
    use crate::ng::types::ContigId;

    let mut compared = 0usize;
    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        // A handful of cases per seed: the reference is drawn from the same alphabet in
        // every one, so this is about covering the *alphabet*, not the case count.
        for index in 0..16 {
            let case = generate(&mut rng);
            let fasta = case.fasta();
            for contig in 0..CONTIGS {
                let length = case.reference[contig].len() as u32;
                let theirs = MultiChromRefFetcher::fetch(&fasta, contig as u32, 1, length)
                    .expect("the whole fixture contig is in range");
                let ours = RefSeq::fetch(&fasta, ContigId(contig as u32), 1, u64::from(length))
                    .expect("the whole fixture contig is in range");
                assert_eq!(
                    ours, theirs,
                    "seed {seed:#x} case {index} contig {contig}: the two views of the \
                     fixture reference disagree, so the differential is comparing two \
                     walks over different bases"
                );
                compared += 1;
            }
        }
    }
    assert!(compared > 0, "no fixture reference was compared");
}

/// **The adaptor counter's real claim, checked end to end.**
///
/// `reads_with_live_adaptor_boundary` is tallied at generation and cannot see the walk, so
/// however carefully it computes liveness it remains a statement about the *input*. The
/// only honest test of "the G1 filter is on this walk" is that **removing it changes the
/// answer** — which is exactly what B2's mutation 2 exploits, made permanent here so the
/// property survives a generator change rather than resting on a hand-run exercise.
#[test]
fn the_adaptor_filter_changes_the_records_the_walk_emits() {
    let mut cases_changed = 0usize;
    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for _ in 0..CASES_PER_SEED {
            let case = generate(&mut rng);
            if case.reads_with_live_adaptor_boundary == 0 {
                continue;
            }
            let without = Case {
                reference: case.reference.clone(),
                reads: case
                    .reads
                    .iter()
                    .cloned()
                    .map(|mut read| {
                        read.adaptor_boundary = None;
                        read
                    })
                    .collect(),
                config: case.config,
                reads_with_live_adaptor_boundary: 0,
            };
            if production_walk(&case).records != production_walk(&without).records {
                cases_changed += 1;
            }
        }
    }
    assert!(
        cases_changed > 0,
        "clearing every adaptor boundary changed no case's records — the generator places \
         boundaries where they silence nothing, so the live-boundary count is measuring the \
         RNG rather than the walk"
    );
}

// ---------------------------------------------------------------------
// B3 — the differential at scale, on real reads
// ---------------------------------------------------------------------

/// **The differential on real alignments** — GIAB HG002 and a tomato CRAM (spec §13.1).
///
/// `#[ignore]`d, because it needs data that is not in the tree: the bundles live in the
/// main repo under `benchmarks/`, and the container mounts `$HOME/genomes` for the
/// references. Driven by environment, so the same test serves both organisms:
///
/// ```text
/// PVC_PARITY_FASTA=$HOME/genomes/h_sapiens/gca_grch38/GCA_….fna \
/// PVC_PARITY_READS=…/benchmarks/ssr_hg002/bam/10x/HG002_TR_v1.0.1_Tier_10x.bam \
/// PVC_PARITY_REGION=chr1:1000000-1200000 \
///   cargo test --release --lib ng_diverges_from_production_on_real_reads_only_where_a_read_did_not_witness \
///     -- --ignored --nocapture
/// ```
///
/// **`--release` is not an optimisation here, it is the point.** Real paired-end data hits
/// production's reachable `debug_assert!` (see
/// `both_walkers_panic_on_a_deletion_anchored_before_its_record`) constantly — overlapping
/// mates carrying deletions are ordinary — so a debug build would abort the walk early and
/// measure almost nothing. Release is also where the walker actually runs.
///
/// # What it compares, and why one fetcher
///
/// Reads are ingested through ng's own step 1 and prepared **once** by ng's
/// `LeftAlignPreparer`; production's stream is that same stream converted down through
/// `PreparedRead::into_production`. Preparing twice would compare two different inputs.
/// A single `WindowedRefSeq` is lent to both walkers for the same reason: identical bytes
/// by construction, so any divergence is the walk's. ng's walker takes it directly (A0);
/// production's reaches it through [`SharedReference`], the local adaptor below.
#[test]
#[ignore = "needs a real BAM/CRAM and reference; see the doc comment for the invocation"]
fn ng_diverges_from_production_on_real_reads_only_where_a_read_did_not_witness() {
    use std::path::PathBuf;

    use crate::ng::read::ReadFilterConfig;
    use crate::ng::read::input::SampleReads;
    use crate::ng::read::left_align::LeftAlignPreparer;
    use crate::ng::read::{ReadPreparer, prepared_read::PreparedRead as NgPreparedRead};
    use crate::ng::ref_seq::WindowedRefSeq;
    use crate::ng::reference_info::{ReferenceInfoCache, read_reference_verifying_or_creating_fai};
    use crate::ng::types::{ContigId, GenomeRegion, Position};
    use std::sync::Arc as StdArc;

    let Ok(fasta) = std::env::var("PVC_PARITY_FASTA") else {
        panic!("set PVC_PARITY_FASTA to a reference FASTA with a sibling .fai");
    };
    let Ok(reads_path) = std::env::var("PVC_PARITY_READS") else {
        panic!("set PVC_PARITY_READS to a coordinate-sorted, indexed BAM or CRAM");
    };
    let fasta = PathBuf::from(fasta);

    // The convenience path, for two reasons: it sets `ReferenceInfo.fasta_path`, without
    // which a CRAM cannot be opened at all (this test serves the tomato CRAMs as well as
    // the HG002 BAMs), and with a `.fai` already present it verifies in the **background**
    // rather than making a 3 GB whole-genome pass before the first read is decoded.
    let cache = StdArc::new(ReferenceInfoCache::new());
    let (reference_info, verification) =
        read_reference_verifying_or_creating_fai(&cache, fasta.clone())
            .expect("the reference is readable and has (or can derive) a .fai");
    let contigs = reference_info.contig_list();

    // `name:start-end`, 1-based inclusive — the same convention `GenomeRegion` uses, so
    // the string a user types and the region the walk covers are the same numbers.
    // Defaults to the first 200 kb of the first contig, which is enough to be a real
    // measurement and small enough to hold two prepared read streams in memory.
    let region_spec = std::env::var("PVC_PARITY_REGION").unwrap_or_default();
    let region = if region_spec.is_empty() {
        let first = &contigs.entries[0];
        GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(first.length.min(200_000)),
        }
    } else {
        let (name, span) = region_spec
            .split_once(':')
            .expect("PVC_PARITY_REGION looks like chr1:1000000-1200000");
        let (start, end) = span
            .split_once('-')
            .expect("PVC_PARITY_REGION looks like chr1:1000000-1200000");
        let index = contigs
            .entries
            .iter()
            .position(|entry| entry.name == name)
            .unwrap_or_else(|| panic!("the reference has no contig named {name}"));
        GenomeRegion {
            contig: ContigId(index as u32),
            start: Position(start.parse().expect("a 1-based start")),
            end: Position(end.parse().expect("an inclusive end")),
        }
    };

    let sample = SampleReads::open_only_sample(
        &[PathBuf::from(&reads_path)],
        &reference_info,
        ReadFilterConfig::default(),
        true,
    )
    .expect("the alignment file opens against this reference");

    let stream = sample
        .reads_in_region(region, || {
            WindowedRefSeq::new(fasta.clone(), contigs.clone())
        })
        .expect("the region query opens");

    // Prepared **once**. The preparer holds its own accessor, separate from the ones the
    // query factory mints — read preparation's rule.
    let preparer = LeftAlignPreparer::with_default_normalizer(WindowedRefSeq::new(
        fasta.clone(),
        contigs.clone(),
    ));
    let mut scratch = <LeftAlignPreparer<WindowedRefSeq> as ReadPreparer>::Scratch::default();

    let mut ng_reads: Vec<NgPreparedRead> = Vec::new();
    for item in stream {
        let read = item.expect("the read stream is readable");
        if let Some(prepared) = preparer
            .prepare_read(read, &mut scratch)
            .expect("preparation does not fail on a well-formed reference")
        {
            ng_reads.push(prepared);
        }
    }
    assert!(
        !ng_reads.is_empty(),
        "no reads in {region:?} of {reads_path} — the region is empty, so this run would \
         prove nothing"
    );
    let prepared_reads = ng_reads.len();
    let production_reads: Vec<ProductionPreparedRead> = ng_reads
        .iter()
        .cloned()
        .map(NgPreparedRead::into_production)
        .collect();

    // One reference, lent to both walks: identical bytes by construction, and one reader
    // rather than two — building `WindowedRefSeq` twice would re-read the `.fai` and open
    // the FASTA again, which is the per-query cost `4bc3ef9` removed on the STR side.
    let reference = SharedReference(Rc::new(WindowedRefSeq::new(fasta.clone(), contigs.clone())));
    let config = WalkerConfig::default();

    // Driven one after the other, not interleaved, so the shared reader's `RefCell` is
    // never borrowed by both walks at once.
    let theirs = drive_production(
        production_run(production_reads, reference.clone(), &config),
        render_production_error,
        |walker| production_counters(walker.summary()),
    );
    let ours = drive_ng(
        super::run(ng_reads, reference, &config),
        render_ng_error,
        |walker| ng_counters(walker.summary()),
    );

    let where_ = format!("{reads_path} {region:?}");
    let ok_records = theirs.records.iter().filter(|item| item.is_ok()).count();
    let first_error = theirs
        .records
        .iter()
        .find_map(|item| item.as_ref().err().cloned());
    // **The census, on real alignments** — the same comparison the synthetic one runs, and
    // the same two claims: every record still exists at the same anchor with the same REF
    // bytes, and the number whose alleles moved is the size of production's defect on real
    // data (spec §13.2, D3's headline). Not `assert_same_walk`: from A2 the two walkers
    // differ on purpose wherever a read did not witness a whole footprint, which on
    // paired-end sequencing is most long-deletion loci.
    let (records_compared, diverged) = assert_only_allele_bytes_moved(&where_, &ours, &theirs);
    let panicked = ours.panic_message.is_some();

    assert!(
        !panicked,
        "{where_}: both walkers panicked, which agrees but measures nothing — re-run with \
         --release, where production's debug_assert! is compiled out"
    );
    // **The floor the synthetic differential has and this one lacked.** Every `WalkerError`
    // is fatal and terminal for the iterator, so a walk that dies on its first read yields
    // one identical `Err` on each side, `assert_same_walk` agrees, and the run prints
    // "1 records compared, zero divergences" — green, and proof of nothing. This is the
    // only evidence in the milestone that the two agree on real data, and it is hand-run,
    // so a green-but-empty run must not look like a real one.
    assert!(
        first_error.is_none(),
        "{where_}: the walk ended in a WalkerError, so the two streams agree only up to \
         where both stopped: {first_error:?}"
    );
    assert!(
        ok_records * 4 > prepared_reads,
        "{where_}: {ok_records} records from {prepared_reads} prepared reads — far too few \
         for a walk that ran to completion, so this run compared a prefix and proves nothing"
    );
    // Joined at the end rather than dropped: a `.fai` that does not describe this FASTA
    // would mean the two walks agreed about the wrong bases, which is a green run that
    // proves nothing.
    if let Some(handle) = verification {
        handle
            .join()
            .expect("the .fai beside the reference describes it");
    }
    eprintln!(
        "real-data differential: {where_} — {records_compared} records compared from \
         {prepared_reads} prepared reads; every anchor, REF sequence and counter identical. \
         {diverged} records ({:.2}%) carried bases production credited to a read that had \
         not witnessed them.",
        100.0 * diverged as f64 / records_compared.max(1) as f64,
    );
}

/// **One reference, two traits — the real-data differential's shared accessor.**
///
/// Production's walker takes a `MultiChromRefFetcher`; ng's takes a
/// [`RefSeq`](crate::ng::ref_seq::RefSeq). Lending **one** accessor to both is what makes
/// "identical bytes by construction" true rather than hopeful, so one value implements
/// both — exactly the shape `MockFasta` has on the synthetic side (`mock_reference.rs`),
/// where production's own impl and ng's meet on one type.
///
/// **This is not `RefSeqFetcher` under a new name.** That type made *ng's walker* speak
/// production's trait, which is what A0 deleted. This makes *production's* walker speak
/// ng's reference, for the length of one `#[ignore]`d differential, and it dies with this
/// file when the two walkers begin to diverge on purpose (A2).
///
/// `Rc`, not `Arc`: [`WindowedRefSeq`](crate::ng::ref_seq::WindowedRefSeq) buffers behind
/// a `RefCell` and is therefore not `Sync`, and this test is single-threaded. The two
/// walks are driven one after the other, so the buffer is never borrowed by both.
struct SharedReference(Rc<crate::ng::ref_seq::WindowedRefSeq>);

impl Clone for SharedReference {
    fn clone(&self) -> Self {
        Self(Rc::clone(&self.0))
    }
}

impl crate::ng::ref_seq::RefSeq for SharedReference {
    fn fetch_into(
        &self,
        contig: crate::ng::types::ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), crate::ng::ref_seq::RefSeqError> {
        self.0.fetch_into(contig, start_1based, length, dst)
    }
}

impl crate::fasta::MultiChromRefFetcher for SharedReference {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, crate::fasta::ChromRefFetchError> {
        crate::ng::ref_seq::RefSeq::fetch(
            self,
            crate::ng::types::ContigId(chrom_id),
            u64::from(start_1based),
            u64::from(length),
        )
        // The failure is reported as I/O carrying ng's own rendering. It is not a
        // faithful variant-for-variant translation and does not need to be: the two sides
        // are compared through `FetchFailure`, which reads the *ng* error on ng's side,
        // and this test refuses any error outright (`first_error.is_none()` above)
        // because a walk that stopped early proves nothing. A divergence here is
        // therefore loud, and it names the real cause.
        .map_err(|error| crate::fasta::ChromRefFetchError::Io {
            chrom_name: format!("chrom_id {chrom_id}"),
            source: std::io::Error::other(format!("{error:?}")),
        })
    }
}
