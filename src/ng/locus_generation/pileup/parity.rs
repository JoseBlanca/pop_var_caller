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
//! REF-only widening, the read-group split. What survives is narrower, and **narrower in a
//! different way than this header used to claim**: not "loci where every folded read
//! witnessed the whole footprint", which class 6 disproves — a read can witness every
//! position and production still be wrong about it — but *every* locus of a fixture on which
//! production fabricates nothing (`generate_uniform_events`; spec §3). **So this is the last
//! moment the baseline can be banked**, and `the_generator_exercises_what_the_port_can_break`
//! is what says the banking is worth anything: a differential that has never been shown
//! to fail is a claim, not evidence.

use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;
use std::sync::Arc;

use super::tests::MockFasta;
use super::{PreparedRead, WalkerConfig};
use crate::ng::read::PLACEHOLDER_READ_GROUP;
use crate::ng::types::{ContigId, GenomeRegion, Position};
// Aliased, so which walker a call reaches is legible at the call site rather than carried
// by a `super::` — this file is the one place both are in scope at once.
use super::super::{LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation};
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

/// A case in which every read spans the contig whole **and shares one event set with every
/// other read on its contig** — the permanent anchor's fixture (D1).
///
/// # Why the events have to be shared, which `generate_complete` did not know
///
/// "Every read witnessed the whole footprint" is not enough to make the two walkers agree,
/// and D1 found the case that proves it. Production's `widen` appends reference bases to
/// **every** bucket and re-folds nothing; ng's re-folds every live read
/// (`refold_live_reads`, A3/option (b) — a function production has no counterpart to). So a
/// record that widens after a read folded into it leaves production holding that read's
/// haplotype computed against a **stale, narrower** footprint: reference bases where the
/// read's own were, and events inside the widened region missed entirely. The read is a
/// complete witness on ng's side and production is still wrong about it. (Both directions
/// occur — the divergence that found this had ng's row one base *longer* than production's,
/// an insertion anchored in the widened region.)
///
/// One shared event set removes the cause rather than filtering the symptom — **not by
/// stopping records widening, which it does not.** Records widen here as often as anywhere
/// else: a pileup opens a record at every covered position, so a deletion anchored at one
/// widens the record already standing there. What the shared event set changes is *who is
/// left stale by it*. Every read in the record carries **its own copy of the widening
/// event**, so every read is re-folded by that event rather than having reference bases
/// appended to a bucket nobody revisits. Production's append-to-every-bucket still runs; it
/// is simply overwritten, for every read, by the read's own subtract-then-add.
///
/// The fixture's defining property is asserted rather than argued — the anchor checks that
/// every read on a contig carries **one CIGAR**, and separately that widened records were
/// reached at all, so neither "the generator stopped sharing events" nor "the widen path
/// went unexercised" can quietly weaken the anchor.
///
/// What still varies, because it is what the anchor is for: the bases at every matched
/// position (so records carry several alleles), qualities, MAPQ, strand, pairing and mate
/// overlap, the column caps, and both contigs.
fn generate_uniform_events(rng: &mut SplitMix64) -> Case {
    let reference: Vec<String> = (0..CONTIGS)
        .map(|_| {
            (0..CONTIG_LENGTH)
                .map(|_| b"ACGT"[rng.below(4)] as char)
                .collect()
        })
        .collect();

    // One template per contig: reads on different contigs never share a record, so their
    // events cannot widen each other's.
    let templates: Vec<ProductionPreparedRead> = (0..CONTIGS)
        .map(|contig| {
            generate_complete_read(
                rng,
                &reference,
                contig as u32,
                "template",
                ProductionMateRole::Solo,
            )
        })
        .collect();

    let mut reads: Vec<ProductionPreparedRead> = Vec::new();
    let read_count = 2 + rng.below(10);
    for index in 0..read_count {
        let chrom_id = if CONTIGS > 1 && rng.one_in(6) { 1 } else { 0 };
        let contig = reference[chrom_id as usize].as_bytes();
        // The template's CIGAR, with this read's own bases over it. The matched positions
        // are walked in the same order the template built them, so the seq stays the length
        // the CIGAR demands — a read whose seq and CIGAR disagree is rejected at admission
        // and the case would test the rejection path instead of the fold.
        fn vary(
            rng: &mut SplitMix64,
            contig: &[u8],
            chrom_id: u32,
            template: &ProductionPreparedRead,
            qname: String,
            role: ProductionMateRole,
        ) -> ProductionPreparedRead {
            let mut seq: Vec<u8> = Vec::with_capacity(template.seq.len());
            let mut ref_pos = 0usize;
            for op in &template.cigar {
                match op {
                    CigarOp::Match(len) | CigarOp::SeqMatch(len) | CigarOp::SeqMismatch(len) => {
                        for offset in 0..*len as usize {
                            seq.push(if rng.one_in(6) {
                                b"ACGT"[rng.below(4)]
                            } else {
                                contig[ref_pos + offset]
                            });
                        }
                        ref_pos += *len as usize;
                    }
                    CigarOp::Insertion(len) | CigarOp::SoftClip(len) => {
                        seq.extend((0..*len).map(|_| b"ACGT"[rng.below(4)]));
                    }
                    CigarOp::Deletion(len) | CigarOp::Skip(len) => ref_pos += *len as usize,
                    CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
                }
            }
            ProductionPreparedRead {
                chrom_id,
                alignment_start: template.alignment_start,
                alignment_end: template.alignment_end,
                cigar: template.cigar.clone(),
                bq_baq: (0..seq.len()).map(|_| 20 + rng.below(21) as u8).collect(),
                seq,
                mq_log_err: -3.0 - (rng.below(4) as f64),
                mapq: 20 + rng.below(41) as u8,
                is_reverse_strand: rng.one_in(2),
                qname: Arc::from(qname.as_str()),
                mate_role: role,
                adaptor_boundary: None,
            }
        }
        let template = &templates[chrom_id as usize];
        if rng.one_in(3) {
            let qname = format!("pair{index}");
            for role in [
                ProductionMateRole::FirstOfPair,
                ProductionMateRole::SecondOfPair,
            ] {
                reads.push(vary(rng, contig, chrom_id, template, qname.clone(), role));
            }
        } else {
            reads.push(vary(
                rng,
                contig,
                chrom_id,
                template,
                format!("solo{index}"),
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
/// claim. `PileupGeneratorCounts::fold_region_walk` destructures for exactly this reason.
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
///
/// **One field is ng's alone and is dropped by name.** `reads_silent_over_footprint` (D2)
/// counts reads that were admitted and never contributed anywhere; production has no
/// counterpart, so there is nothing to compare and comparing it against a structural zero
/// would fail every walk. It is bound and dropped here, at one site, with this comment — the
/// same treatment `placed_start` gets in [`project`], and the reason the destructure is
/// exhaustive is so that treatment has to be *chosen*.
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
        reads_silent_over_footprint: _,
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
    /// **ng's type on both sides, from D1 on.** Production's `PileupRecord`s arrive here
    /// through [`project`]; ng's loci arrive as they are emitted. The comparison surface is
    /// the type ng actually ships, so a field ng carries and production cannot say is a
    /// *named* divergence rather than something the old back-projection quietly merged away.
    records: Vec<Result<SampleLocusObservations, String>>,
    /// What the projection had to drop to produce each stream item — one entry per item, so
    /// the two index together. Class 4 (and class 5) are facts about the *projection*, not
    /// about the emitted locus, so they cannot be read back off `records`: a dropped bucket
    /// leaves no trace in what it was dropped from. Default — nothing dropped — for an error
    /// item and for every item of an ng walk, which projects nothing.
    drops: Vec<ProjectionDrops>,
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
    let mut drops = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(
                item.map(|record| {
                    let (locus, dropped) = project_counting_drops(&record);
                    drops.push(dropped);
                    locus
                })
                .map_err(|error| {
                    drops.push(ProjectionDrops::default());
                    render_error(&error)
                }),
            );
        }
        summary = Some(summary_of(&walker));
    }))
    .err()
    .map(panic_message);
    WalkOutcome {
        records,
        drops,
        summary,
        panic_message,
    }
}

/// The same as [`drive_production`], for a walker that yields ng's own locus type — which
/// needs no projection at all, being the surface both sides are now compared on. Separate
/// rather than generic over the item type: which side projects is the thing worth seeing at
/// the call site, and the asymmetry is the whole of D1.
fn drive_ng<W, E>(
    mut walker: W,
    render_error: impl Fn(&E) -> String,
    summary_of: impl FnOnce(&W) -> SummaryCounters,
) -> WalkOutcome
where
    W: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    let mut records = Vec::new();
    let mut drops = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(item.map_err(|error| render_error(&error)));
            // ng projects nothing, so nothing was dropped on the way here. The vector is
            // kept the same length as `records` so the two index together whichever walk an
            // outcome came from.
            drops.push(ProjectionDrops::default());
        }
        summary = Some(summary_of(&walker));
    }))
    .err()
    .map(panic_message);
    WalkOutcome {
        records,
        drops,
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
    ng_walk_in_groups(case, 1)
}

/// ng's answer with the sample's reads dealt round-robin into `groups` read groups.
///
/// **The read group is the one part of a row's identity production has no way to say**, so
/// spec §3's class 2 — "an allele supported from several read groups is several rows" — is
/// unreachable on a one-group fixture: it would be a class counted zero forever, which is
/// how an unexercised branch passes for a tested one. Dealing the same reads into two groups
/// is what makes the class fire, and it is safe to do here because the group reaches
/// **nothing but the row key**: `open_record.rs` reads `read_group` only to build
/// `ObservationKey` and to order the emitted rows. So the walk is bit-for-bit the walk of
/// `groups == 1`, and the only difference in the output is the split this class names —
/// which is pinned in two places — **not** by a
/// `two_read_groups_split_rows_without_moving_evidence`, which this comment cited until
/// 2026-07-30 and which exists nowhere in the repository. What discharges it: this module's
/// [`evidence_by_bases`] reconciliation, applied at **every** census locus, which groups the
/// evidence by bases and so is blind to how many rows carry it; and the dump tool's
/// `two_read_groups_split_one_allele_into_two_rows_that_sum_back`, which asserts the split
/// rows' `num_obs` sum back to the one-group total on a committed fixture.
fn ng_walk_in_groups(case: &Case, groups: u32) -> WalkOutcome {
    assert!(groups >= 1, "a walk has at least one read group");
    let fasta = case.fasta();
    let reads: Vec<PreparedRead> = case
        .reads
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, read)| {
            let group = if groups == 1 {
                PLACEHOLDER_READ_GROUP
            } else {
                crate::ng::types::ReadGroupId(index as u32 % groups)
            };
            PreparedRead::from_production(read, group)
        })
        .collect();
    drive_ng(super::run(reads, fasta, &case.config), render_ng_error, {
        |walker| ng_counters(walker.summary())
    })
}

/// The read group every projected row carries: production has none, so one is chosen and
/// named rather than defaulted at four call sites.
///
/// It is [`PLACEHOLDER_READ_GROUP`], the group [`ng_walk`] gives every read on a one-group
/// walk — so on that fixture the field compares equal and the group is not a divergence.
/// On a multi-group walk it is what spec §3's class 2 diverges *from*.
const PROJECTED_READ_GROUP: crate::ng::types::ReadGroupId = PLACEHOLDER_READ_GROUP;

/// **Production's record, said in ng's terms — the stage-2 projection** (spec §3).
///
/// This is the direction D1 exists to reverse. Until now the differential ran the other way,
/// laying ng's locus back out as a `PileupRecord` through `to_pileup_record` — and a
/// back-projection can only ever compare what the *older, smaller* type can hold. Every
/// field B2 added (`read_witness`, `read_group`, the two per-record counters, the region's
/// end) was merged or dropped on the way out, which is why Milestone B's review found three
/// live surfaces the suite could not see. Projecting **forward** makes each of them visible:
/// production has no way to say them, so each shows up as a difference with a name.
///
/// # Total, and the drops are the classes
///
/// `PileupRecord` is destructured exhaustively, so a field added to production's type stops
/// this compiling instead of being silently left out of the oracle — the same rule
/// [`RecordEvidence`] earned the hard way. What each field becomes:
///
/// - `chrom_id`, `pos` and `alleles[0].seq.len()` become the [`GenomeRegion`]. The **end**
///   is production's only statement of the footprint's extent, and it is implicit in the REF
///   bucket's length; ng carries it explicitly. Nothing projected *back* could see it.
/// - `alleles[0].seq` becomes `reference_bases` — the REF bucket's bytes *are* the record's
///   reference sequence, which is why ng does not store the two separately.
/// - each allele with `num_obs > 0` becomes one [`SequenceObservation`]. **Class 4:** the ones
///   with no support are dropped, because production creates `alleles[0]` at record open
///   regardless and A3's widen strands empty non-REF buckets, while ng derives rows from
///   reads that folded. `dropped_unsupported` on the returned census says how many, and
///   the caller asserts nothing supported was ever dropped.
/// - `read_witness` is **`Complete` on every projected row. Class 1:** production has no
///   notion of a read that witnessed part of a footprint — that is the defect, and this is
///   where it becomes visible instead of being merged away.
/// - `read_group` is [`PROJECTED_READ_GROUP`] on every projected row. **Class 2:** the
///   group joins ng's row identity and production cannot split a row by it.
/// - `reads_without_observation` and `reads_discarded_by_cap` are **zero. Class 3:**
///   production keeps neither per record — only run-level totals on `RunSummary` — so these
///   are asserted against hand-counted fixtures and against ng's own read accounting, never
///   against production.
/// - `windowed_gc` / `windowed_coverage` are dropped by name: `f32::NAN` placeholders the
///   pileup→`.psp` seam fills from a sliding window the walker cannot compute (spec §3).
///   ng's locus type has no counterpart and deliberately does not gain one here.
/// - `placed_start` is dropped by name from every allele's support: ng stops computing it at
///   B2 (spec §6), so there is nothing on the other side to compare and comparing
///   production's value against a structural zero would fail every locus. The difference
///   between this and the *oversight* it replaced is spelled out on [`RecordEvidence`].
///
/// Rows come out in ng's emission order — [`sort_rows`] — so **class 5** (production's order
/// is bucket-creation order) is normalised here and counted by the census rather than
/// tolerated silently.
fn project(record: &PileupRecord) -> SampleLocusObservations {
    project_counting_drops(record).0
}

/// How much [`project`] had to drop — the number that turns class 4 from an excuse into a
/// measurement.
#[derive(Debug, Default, PartialEq, Eq)]
struct ProjectionDrops {
    /// Alleles dropped for having no supporting read.
    ///
    /// *(There is deliberately no "observations dropped" field beside this: the filter
    /// **is** `num_obs == 0`, so such a counter could only ever read zero — a number that
    /// cannot move is not evidence. The claim "no supported allele was dropped" is made
    /// where it can fail, by [`assert_reads_are_accounted_for`], which reconciles ng's
    /// per-locus observation total against production's.)*
    unsupported: usize,
    /// Chain ids carried by those alleles — a claim that **can** fail. Production's re-fold
    /// subtracts a read's scalars from the bucket it leaves and clears its id
    /// (`refold_after_widen_clears_chain_id_from_old_bucket`); an id stranded on a bucket at
    /// `num_obs == 0` would be a read whose identity this projection silently discarded.
    unsupported_chain_ids: usize,
    /// Whether production's own allele order differed from ng's emission order — class 5,
    /// counted so "ng sorts" is shown to be doing work.
    reordered: bool,
}

fn project_counting_drops(record: &PileupRecord) -> (SampleLocusObservations, ProjectionDrops) {
    // Exhaustive: a field added to production's record stops this compiling.
    let PileupRecord {
        chrom_id,
        pos,
        alleles,
        // Filled at the pileup→`.psp` seam from a window the walker cannot see; ng's locus
        // type has no counterpart (spec §3).
        windowed_gc: _,
        windowed_coverage: _,
    } = record;

    let reference_bases: Box<[u8]> = alleles
        .first()
        .expect("production creates alleles[0] with the record")
        .seq
        .clone()
        .into_boxed_slice();

    let mut drops = ProjectionDrops::default();
    let mut rows: Vec<SequenceObservation> = Vec::with_capacity(alleles.len());
    for allele in alleles {
        // Exhaustive on the support too: `placed_start` is bound and dropped **by name**.
        let AlleleSupportStats {
            num_obs,
            q_sum,
            fwd,
            placed_left,
            placed_start: _,
            mapq_sum,
            mapq_sum_sq,
        } = allele.support;
        if num_obs == 0 {
            drops.unsupported += 1;
            drops.unsupported_chain_ids += allele.chain_ids.len();
            continue;
        }
        let mut chain_ids = allele.chain_ids.clone();
        chain_ids.sort_unstable();
        chain_ids.dedup();
        rows.push(SequenceObservation {
            bases: allele.seq.clone().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: PROJECTED_READ_GROUP,
            num_obs,
            num_fwd: fwd,
            q_sum,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
            chain_ids,
        });
    }
    let before = rows.clone();
    sort_rows(&mut rows);
    drops.reordered = before != rows;

    let locus = SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(*chrom_id),
            // 1-based **inclusive**: the last covered position, not one past it.
            start: Position(u64::from(*pos)),
            end: Position(u64::from(*pos) + reference_bases.len() as u64 - 1),
        },
        reference_bases,
        observations: rows,
        // Class 3 — production keeps neither per record.
        reads_without_observation: 0,
        reads_discarded_by_cap: 0,
        kind: LocusKind::Generic,
    };
    (locus, drops)
}

/// ng's emission order, applied to both sides — spec §3's class 5.
///
/// The `ReadWitness` half of the comparator is **the walk's own**
/// (`open_record::witness_order`, lifted to `pub(super)` for this) rather than a second
/// spelling here: `finalise` sorts `ObservationRow`s and this sorts `SequenceObservation`s, so
/// the loop cannot be shared, but the one piece that could silently drift is.
/// `the_projection_orders_rows_as_the_walk_does` covers the rest, by asserting that sorting
/// an ng locus's rows with this function leaves them where the walk emitted them.
fn sort_rows(rows: &mut [SequenceObservation]) {
    use super::open_record::witness_order;
    rows.sort_by(|a, b| {
        a.bases
            .cmp(&b.bases)
            .then_with(|| witness_order(a.read_witness).cmp(&witness_order(b.read_witness)))
            .then_with(|| a.read_group.0.cmp(&b.read_group.0))
    });
}

/// The surface the two sides are compared on — **three named projections, applied to both**
/// (spec §3).
///
/// 1. `q_sum` compared within [`Q_SUM_TOLERANCE`] — a **relative** tolerance, applied by
///    [`ComparableLocus`]'s `PartialEq` rather than by normalising the value — because ng
///    changed the **order** the sum accumulates in and nothing else. It was a fixed 1e-9
///    rounding until D3; the two constants it named, `Q_SUM_GRAIN` and `round_q_sum`, no
///    longer exist, and the reason the grain had to go is on [`Q_SUM_TOLERANCE`].
/// 2. Rows sorted into ng's emission order (class 5). ng's are already; production's are in
///    bucket-creation order.
/// 3. The two per-record counters zeroed (class 3). Production has no counterpart, so an ng
///    locus that counted a read out would otherwise differ on every comparison that has
///    nothing to do with the rows. They are asserted **directly**, on ng's own locus, by
///    [`assert_reads_are_accounted_for`] — which is a stronger claim than any equality
///    against a structural zero could be.
///
/// Applied to **both** sides, so neither can hide a difference the other does not have.
fn comparable(locus: &SampleLocusObservations) -> ComparableLocus {
    ComparableLocus(comparable_exact_q_sum(locus))
}

/// **The relative tolerance `q_sum` is compared at, and the reason it is relative.**
///
/// `q_sum` is an `f64` running sum and `f64` addition is not associative. Two changes reorder
/// it, and neither moves any evidence:
///
/// - **A3's eviction.** Production keeps a bucket alive at `num_obs == 0` and keeps
///   accumulating into it, so a read that leaves and returns leaves `+q -q +q` behind; ng
///   evicts the empty bucket and recreates it, so the same read's sum starts from `0.0` and is
///   *exactly* `q`. Production's `-2.999999999999999` against ng's `-3.0` is this.
/// - **B1's per-read re-derivation.** Production's bucket total is accumulated during the walk
///   with a subtract-then-add on every re-fold; ng's row sums each read's contribution **once**,
///   in `read_id` order. Same addends, different order — and ng's is the more accurate of the
///   two, since nothing cancels.
///
/// # It was a fixed 1e-9 grain, and real data at 300× broke it
///
/// The original comparison rounded both sides to nine decimal places, on the argument that
/// "the grain is nine decimal places on values of order −3 to −50, where the smallest *real*
/// difference is a whole read's `ln` contribution — order 1". The premise is a statement about
/// **depth**: a locus's `q_sum` is order −3 only when a handful of reads support it. D3's first
/// real-data run — HG002 at **300×**, chr1 — hit a locus with 414 observations and
/// `q_sum ≈ −3360.39`, where 1e-9 *absolute* is 3.4 × 10¹² grains and one reordered
/// accumulation lands a single grain apart:
///
/// ```text
/// ng   … q_sum_rounded: -3360392684715
/// prod … q_sum_rounded: -3360392684716
/// ```
///
/// The census reported it as an **unlisted divergence**, which is exactly what it is supposed
/// to do with a difference it cannot name — and the thing that could not be named was the
/// tolerance, not the walk.
///
/// So the comparison is a **relative** one, and it is a tolerance rather than a rounding.
/// Rounding decides equality by which side of a grain boundary each value falls, so two values
/// one ulp apart can round apart however fine the grain — at millions of loci that happens.
/// A tolerance cannot: `|a − b| ≤ ε · max(1, |a|, |b|)` is true for every pair that close.
///
/// At `ε = 1e-9` the allowance at that locus is 3.4 × 10⁻⁶, and the smallest *real* difference
/// is still a whole read's contribution — order 1. Six orders of magnitude of headroom, at any
/// depth, which is what the fixed grain only had at low ones.
const Q_SUM_TOLERANCE: f64 = 1e-9;

/// Whether two `q_sum`s differ by no more than accumulation order can explain — see
/// [`Q_SUM_TOLERANCE`].
fn q_sum_close(ours: f64, theirs: f64) -> bool {
    let scale = ours.abs().max(theirs.abs()).max(1.0);
    (ours - theirs).abs() <= Q_SUM_TOLERANCE * scale
}

/// A locus compared the way this differential means to compare loci: **every field exactly,
/// except `q_sum`, which is compared within [`Q_SUM_TOLERANCE`]**.
///
/// A newtype because `SampleLocusObservations`' own `PartialEq` is derived and compares `f64`s
/// exactly — which is right for the shared type, and wrong for a differential whose whole
/// subject is two implementations accumulating the same addends in different orders.
struct ComparableLocus(SampleLocusObservations);

impl PartialEq for ComparableLocus {
    fn eq(&self, other: &Self) -> bool {
        let (ours, theirs) = (&self.0, &other.0);
        // Exhaustive destructure on both sides: a field added to the locus type stops this
        // compiling rather than going silently uncompared, which is the same rule
        // `LocusEvidence` and `project` follow.
        let SampleLocusObservations {
            region: our_region,
            reference_bases: our_bases,
            observations: our_rows,
            reads_without_observation: our_without,
            reads_discarded_by_cap: our_capped,
            kind: our_kind,
        } = ours;
        let SampleLocusObservations {
            region: their_region,
            reference_bases: their_bases,
            observations: their_rows,
            reads_without_observation: their_without,
            reads_discarded_by_cap: their_capped,
            kind: their_kind,
        } = theirs;
        our_region == their_region
            && our_bases == their_bases
            && our_without == their_without
            && our_capped == their_capped
            && our_kind == their_kind
            && our_rows.len() == their_rows.len()
            && our_rows.iter().zip(their_rows).all(|(ours, theirs)| {
                let SequenceObservation {
                    bases: our_row_bases,
                    read_witness: our_witness,
                    read_group: our_group,
                    num_obs: our_obs,
                    num_fwd: our_fwd,
                    q_sum: our_q_sum,
                    mapq_sum: our_mapq,
                    mapq_sum_sq: our_mapq_sq,
                    placed_left: our_placed_left,
                    chain_ids: our_ids,
                } = ours;
                our_row_bases == &theirs.bases
                    && our_witness == &theirs.read_witness
                    && our_group == &theirs.read_group
                    && our_obs == &theirs.num_obs
                    && our_fwd == &theirs.num_fwd
                    && our_mapq == &theirs.mapq_sum
                    && our_mapq_sq == &theirs.mapq_sum_sq
                    && our_placed_left == &theirs.placed_left
                    && our_ids == &theirs.chain_ids
                    && q_sum_close(*our_q_sum, theirs.q_sum)
            })
    }
}

impl std::fmt::Debug for ComparableLocus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// The normalisation [`comparable`] wraps, with `q_sum` left exact — so a caller can count the
/// loci that agree **only** within the tolerance, and the tolerance is shown to be doing work
/// rather than quietly matching nothing.
fn comparable_exact_q_sum(locus: &SampleLocusObservations) -> SampleLocusObservations {
    let mut out = locus.clone();
    sort_rows(&mut out.observations);
    // **Class 3, zeroed on both sides.** Production keeps neither counter per record, so its
    // side is a structural zero; zeroing ng's too keeps the comparison symmetric, which is
    // the rule every projection here follows — neither side may be normalised in a way the
    // other is not. The counters are not thereby unchecked: they are asserted on ng's own
    // locus by `assert_reads_are_accounted_for`, which is where a claim about a quantity
    // production cannot state belongs.
    out.reads_without_observation = 0;
    out.reads_discarded_by_cap = 0;
    out
}

/// Loci that agree within the `q_sum` tolerance and disagree on an exact comparison.
fn float_only_divergences(
    ours: &[Result<SampleLocusObservations, String>],
    theirs: &[Result<SampleLocusObservations, String>],
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
        // Normalised on both sides — see `comparable` for the three projections and why
        // each is a named class rather than an excuse.
        let normalise = |item: &Result<SampleLocusObservations, String>| {
            item.as_ref().map(comparable).map_err(String::clone)
        };
        assert_eq!(
            normalise(ours),
            normalise(theirs),
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

/// **The six divergence classes spec §3 names**, as the set present at one locus.
///
/// A locus may be in several at once — a widened record with two read groups and a blind
/// read is in three — so this is a set of flags rather than an enum. The rule D1 enforces is
/// not "at most one class" but **"no divergence outside the six"**: [`classify_locus`]
/// panics when two loci differ and no class is present, which is the whole of "counted, not
/// excused".
///
/// The classes are stated as facts about **ng's** locus and the projection's drops, not as
/// descriptions of the difference. That matters: a class read off the difference itself
/// would be a label chosen after the fact, and the one thing spec §3 warns against is an
/// unlisted class getting triaged as a listed one.
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy)]
struct DivergenceClasses {
    /// **Class 1** — a read that did not witness the whole footprint becomes an `Observed`
    /// row whose bases are what it saw, where production folded a reference-filled
    /// haplotype into a bucket. *This is the deliverable* (spec §13.2).
    partial_witness: bool,
    /// **Class 2** — the read group joins ng's row identity and production cannot state it
    /// at all, so any row whose group is not the one [`project`] had to invent is in this
    /// class. The *split* — one allele carried by two groups becoming two rows — is the
    /// visible consequence and is counted separately
    /// ([`DivergenceCensus::group_split_rows`]), because a class flag that only fired on
    /// the split would leave a locus whose single row merely names a different group
    /// unclassified.
    group_split: bool,
    /// **Class 3** — `reads_without_observation` / `reads_discarded_by_cap` are non-zero.
    /// Production keeps neither per record.
    counters: bool,
    /// **Class 4** — production emitted a bucket with `num_obs == 0` and ng emits no such
    /// row; the projection dropped it.
    unsupported_bucket: bool,
    /// **Class 5** — row order: production's is bucket-creation order, ng's is sorted.
    row_order: bool,
    /// **Class 6 — the one spec §3's table does not list, found by D1 building the anchor.**
    ///
    /// Production's `widen` appends reference bases to every bucket and re-folds nothing;
    /// ng's re-folds every live read (`refold_live_reads`, A3/option (b) — a function
    /// production has no counterpart to) and leaves an expired one saying what it saw. So a
    /// record that widens after a read has folded into it leaves production holding that
    /// read's haplotype computed against a **stale, narrower footprint**, in either
    /// direction: reference bases where the read's own belong, or events inside the widened
    /// region missed entirely.
    ///
    /// **It is not class 1 and must not be filed there.** The read in the case that found
    /// this witnessed every position of the final footprint — ng reports it `Complete` —
    /// and production is wrong about it anyway. Filing it under "a read did not witness the
    /// whole footprint" would put reads production mis-folded into the count of reads
    /// production credited with bases they never sequenced, which is the deliverable, and
    /// spec §13.2 asks for the two numbers **separately**. §3's table listed five until this
    /// class was added to it; it now lists six, and §3 carries the mechanism.
    ///
    /// **Which reads this is, and it is not the ones §13.2 originally named.** The gate below
    /// is `!classes.partial_witness`, so a read that had already **expired** before the widen
    /// is *not* here — it leaves ng an `Observed` row, which makes the locus class 1, whose
    /// triple counts it. This class is the other non-contributor: a read **still live** at the
    /// widening step but silent there. Spec §13.2 records that correction.
    stale_widen: bool,
}

impl DivergenceClasses {
    fn any(self) -> bool {
        self.partial_witness
            || self.group_split
            || self.counters
            || self.unsupported_bucket
            || self.row_order
            || self.stale_widen
    }
}

/// A running count of each class over a census, plus the three numbers spec §13.2 calls the
/// deliverable.
#[derive(Debug, Default)]
struct DivergenceCensus {
    loci: usize,
    /// Loci equal to the projection under [`comparable`].
    exact: usize,
    partial_witness: usize,
    group_split: usize,
    counters: usize,
    unsupported_bucket: usize,
    row_order: usize,
    stale_widen: usize,
    /// Loci where one allele really is carried by two read groups and so becomes two rows —
    /// class 2's *visible* consequence, and the check spec §13 asks for by name. Counted
    /// apart from the class flag because the flag fires on any non-projected group, which is
    /// a weaker event and would let this one go to zero unnoticed.
    group_split_rows: usize,
    /// Loci that agree only once `q_sum` is compared within [`Q_SUM_TOLERANCE`] — the
    /// accumulation-order class. A subset of the **exact** count, since `classify_locus`
    /// decides `exact` with [`comparable`], which applies the tolerance. Bounded above as
    /// well as below: see the ceiling at the census's assertions.
    float_only: usize,
    /// **The deliverable, one:** loci carrying at least one partial witness.
    fabricating_loci: usize,
    /// **The deliverable, two:** reads production credited with bases they never sequenced —
    /// the observations sitting on `Observed` rows.
    fabricated_reads: u64,
    /// **The deliverable, three:** reference bases production credited to those reads —
    /// summed as (footprint − positions witnessed) × reads, which is exactly the count of
    /// bases `apply_events_to_ref_into` used to copy out of `ref_seq` on their behalf.
    fabricated_ref_bases: u64,
    /// **The stale-widen deliverable, two** (spec §13.2's *second* triple): reads production
    /// folded against a footprint it never revisited. `stale_widen` above is the loci; these
    /// two are the reads and the bases, and without them the class this milestone *discovered*
    /// was the one number nobody could size.
    stale_widen_reads: u64,
    /// **The stale-widen deliverable, three:** reference bases `widen` appended on those
    /// reads' behalf — the tail of each production row past the point it stops matching any
    /// ng row, times the reads carrying it.
    stale_widen_ref_bases: u64,
}

impl DivergenceCensus {
    fn record(&mut self, classes: DivergenceClasses, exact: bool) {
        self.loci += 1;
        if exact {
            self.exact += 1;
        }
        self.partial_witness += usize::from(classes.partial_witness);
        self.group_split += usize::from(classes.group_split);
        self.counters += usize::from(classes.counters);
        self.unsupported_bucket += usize::from(classes.unsupported_bucket);
        self.row_order += usize::from(classes.row_order);
        self.stale_widen += usize::from(classes.stale_widen);
    }

    /// Add one locus's contribution to the fabrication measurement (spec §13.2), and to the
    /// tally of alleles genuinely split across read groups.
    fn measure_fabrication(&mut self, locus: &SampleLocusObservations) {
        self.group_split_rows += usize::from(rows_split_by_group(locus));
        let footprint = locus.region.len();
        let mut fabricating = false;
        for row in &locus.observations {
            let ReadWitness::Observed {
                positions_covered, ..
            } = row.read_witness
            else {
                continue;
            };
            fabricating = true;
            self.fabricated_reads += u64::from(row.num_obs);
            self.fabricated_ref_bases +=
                footprint.saturating_sub(u64::from(positions_covered)) * u64::from(row.num_obs);
        }
        self.fabricating_loci += usize::from(fabricating);
    }

    /// Add one locus's contribution to the **stale-widen** measurement — spec §13.2's second
    /// three-number deliverable, which until now was a locus count and nothing else.
    ///
    /// **The prefix logic is [`stale_widen_shape`]'s, deliberately duplicated in shape and not
    /// in code path.** That function decides *whether* a locus is class 6 by asking that every
    /// production row ng does not have is some ng row's bases followed by a reference tail;
    /// this one measures *how long* those tails are. Reusing its `any(...)` would give a
    /// boolean, not a length, so what is shared is the definition of the tail — the bytes past
    /// the longest prefix any ng row agrees with. Where several ng rows explain a production
    /// row, the **longest** shared prefix wins, which charges production the smallest tail
    /// consistent with the class.
    fn measure_stale_widen(
        &mut self,
        classes: DivergenceClasses,
        ours: &SampleLocusObservations,
        theirs: &SampleLocusObservations,
    ) {
        if !classes.stale_widen {
            return;
        }
        let ng_bases: BTreeSet<&[u8]> = ours
            .observations
            .iter()
            .map(|row| row.bases.as_ref())
            .collect();
        for row in &theirs.observations {
            let theirs_bases: &[u8] = &row.bases;
            if ng_bases.contains(theirs_bases) {
                continue;
            }
            let longest_shared = ng_bases
                .iter()
                .filter_map(|ours_bases| {
                    let shared = ours_bases
                        .iter()
                        .zip(theirs_bases)
                        .take_while(|(a, b)| a == b)
                        .count();
                    // The same admissibility test `stale_widen_shape` applies, so a row it
                    // would not have called explained cannot contribute a length here.
                    (theirs_bases.len() >= shared
                        && ours.reference_bases.ends_with(&theirs_bases[shared..]))
                    .then_some(shared)
                })
                .max();
            // `None` cannot happen on a locus `stale_widen_shape` accepted — it returns
            // `false` unless every such row is explained — so this is the invariant, stated.
            let Some(shared) = longest_shared else {
                debug_assert!(
                    false,
                    "a class-6 locus has a production row no ng row explains: {:?}",
                    String::from_utf8_lossy(theirs_bases),
                );
                continue;
            };
            self.stale_widen_reads += u64::from(row.num_obs);
            self.stale_widen_ref_bases +=
                (theirs_bases.len() - shared) as u64 * u64::from(row.num_obs);
        }
    }
}

/// **Every** support scalar of a whole locus, summed across its rows — "no evidence was
/// created, lost or double counted", independent of which bucket any of it ended up in.
///
/// # Why this is a struct built by an exhaustive destructure
///
/// It was a six-field tuple read field by field, and it **omitted `placed_start`** — the one
/// field Milestone A made fragile, since A1 dropped it from ng's stats and `finalise`
/// reconstructed it from a per-read flag. That is not a hypothetical: injecting a real wrong
/// number (the REF bucket stops counting `placed_start`) left the permanent anchor **green**,
/// while quietly moving the tolerated class from 264 records (0.15 %) to 2,806 (1.58 %) and
/// printing "same support totals" and "Every other record is identical, field for field" as
/// it went. 2,542 wrong records absorbed. Only an inherited test from production's own suite
/// noticed.
///
/// So the sum is taken by **destructuring [`SequenceObservation`] exhaustively**: a field added
/// to ng's row type stops this file compiling instead of being silently left out of the
/// oracle. D1 turns the comparison around — the sum is now taken over ng's rows on both
/// sides, production's having been projected into them — and the guard has to turn with it,
/// or a field added to `SequenceObservation` would go uncompared exactly as `placed_start` did.
/// (`AlleleSupportStats` is still destructured exhaustively, in [`project`], which is the
/// one place production's type is now read.)
///
/// `q_sum` is summed as `f64` and compared **once, at the end**, within
/// [`Q_SUM_TOLERANCE`] — never per bucket. The reason survives D3's move from a rounding to a
/// tolerance unchanged, and is why the sum is taken here rather than compared bucket by
/// bucket: a per-bucket comparison makes the verdict depend on how the reads were distributed
/// across buckets, and A3 moves reads between buckets by design. (Under the old grain the same
/// hazard was arithmetic rather than a verdict: two reads at `-4.835428695` in one bucket
/// summed to `-9.670857391`, while the same two in separate buckets rounded to `-4835428695`
/// twice and totalled `-9670857390`.)
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
/// ng's row type stops this file compiling; the two fields this sum deliberately excludes
/// are bound and dropped **by name**, at one site, with this comment attached. An oversight
/// cannot look like this.
#[derive(Debug, Default)]
struct LocusEvidence {
    num_obs: u32,
    fwd: u32,
    placed_left: u32,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    q_sum: f64,
}

/// Every integer field exactly; `q_sum` within [`Q_SUM_TOLERANCE`], for the reason given there
/// — and it is the reason this is a hand-written impl rather than a derive over a scaled
/// integer. A grain decides equality by which side of a boundary each value falls on, so two
/// sums one ulp apart can land in different grains; at 300× depth and millions of loci, they
/// do. That is the failure D3's first real-data run produced.
impl PartialEq for LocusEvidence {
    fn eq(&self, other: &Self) -> bool {
        // Exhaustive destructure: a field added to this struct stops the comparison compiling
        // rather than going uncompared.
        let Self {
            num_obs,
            fwd,
            placed_left,
            mapq_sum,
            mapq_sum_sq,
            q_sum,
        } = self;
        num_obs == &other.num_obs
            && fwd == &other.fwd
            && placed_left == &other.placed_left
            && mapq_sum == &other.mapq_sum
            && mapq_sum_sq == &other.mapq_sum_sq
            && q_sum_close(*q_sum, other.q_sum)
    }
}

/// The evidence on one row, so a caller can sum over whichever grouping it is reconciling.
fn row_evidence(row: &SequenceObservation, totals: &mut LocusEvidence, q_sum: &mut f64) {
    // The exhaustive destructure is the point — see the doc above.
    let SequenceObservation {
        // **Not summed, and each for its own reason.** `bases` is *what* the evidence says
        // rather than how much of it there is — it is the thing class 1 moves, and summing
        // it would be the excuse this census exists to refuse. `read_witness` and
        // `read_group` are the two halves of the row identity production cannot state
        // (classes 1 and 2): including them would make every split locus differ here, where
        // the question this function asks is only "did a read go missing".
        bases: _,
        read_witness: _,
        read_group: _,
        num_obs,
        num_fwd,
        q_sum: row_q_sum,
        mapq_sum,
        mapq_sum_sq,
        placed_left,
        // Chain ids are a *set*, not a scalar; they are checked separately and by subset —
        // see `locus_chain_ids` for the record that settles why subset and not equality.
        chain_ids: _,
    } = row;
    totals.num_obs += num_obs;
    totals.fwd += num_fwd;
    totals.placed_left += placed_left;
    totals.mapq_sum += mapq_sum;
    totals.mapq_sum_sq += mapq_sum_sq;
    *q_sum += row_q_sum;
}

fn locus_evidence(locus: &SampleLocusObservations) -> LocusEvidence {
    let mut totals = LocusEvidence::default();
    let mut q_sum = 0.0_f64;
    for row in &locus.observations {
        row_evidence(row, &mut totals, &mut q_sum);
    }
    totals.q_sum = q_sum;
    totals
}

/// The same evidence, gathered per distinct `bases` — **spec §3's class-2 reconciliation**:
/// "summing a locus's rows by `(bases, coverage)` must reproduce production's per-allele
/// totals".
///
/// Grouped by `bases` alone rather than by `(bases, coverage)`, and the difference is
/// deliberate: coverage is class **1**, which moves the bases too, so a grouping that kept
/// coverage in the key would answer "the rows differ" for the class it is not asking about.
/// What this reconciles is the split — several rows where production has one allele — and
/// the split is invisible to `bases`.
fn evidence_by_bases(locus: &SampleLocusObservations) -> BTreeMap<Vec<u8>, LocusEvidence> {
    let mut per_bases: BTreeMap<Vec<u8>, (LocusEvidence, f64)> = BTreeMap::new();
    for row in &locus.observations {
        let entry = per_bases.entry(row.bases.to_vec()).or_default();
        row_evidence(row, &mut entry.0, &mut entry.1);
    }
    per_bases
        .into_iter()
        .map(|(bases, (mut totals, q_sum))| {
            totals.q_sum = q_sum;
            (bases, totals)
        })
        .collect()
}

/// Every chain id a locus carries, sorted and deduplicated — the *other* half of "no
/// evidence moved", and the half [`locus_evidence`] cannot see because chain ids are a set
/// per row rather than a scalar.
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
fn locus_chain_ids(locus: &SampleLocusObservations) -> Vec<ChainId> {
    let mut ids: Vec<ChainId> = locus
        .observations
        .iter()
        .flat_map(|row| row.chain_ids.iter().copied())
        .collect();
    ids.sort_unstable();
    ids.dedup();
    ids
}

/// **Whether ng's locus differs from the projection, and if so under which of the six
/// classes — with every difference outside them a panic** (spec §3, §13.2).
///
/// This replaces `classify_record`, whose own doc said what it was missing: it checked only
/// that *no evidence moved* — the reference bytes and the support totals — and never *which*
/// rows differed or by how much. It could not, because the surface it compared on was
/// production's type, where a partial witness has no spelling. Projecting forward gives
/// coverage a name, so the class can be read off ng's own rows rather than inferred from a
/// difference.
///
/// # What holds at every locus, class or no class
///
/// These are asserted, not classified. A2 moves allele bytes between rows; it does not move
/// a record's existence, its anchor, its footprint or its reads:
///
/// - the **region** and the **reference bases** are identical, and so is the kind;
/// - **the reads are all accounted for**: ng's observations plus the reads it counted out
///   equal production's observations, exactly (see [`assert_reads_are_accounted_for`]);
/// - **no chain id was lost**: production's ids are a subset of ng's. ng may hold an id
///   production dropped into its REF bucket — see [`locus_chain_ids`] for the record that
///   settles why this is a subset check and not an equality.
///
/// # And the reconciliation that makes class 2 a measurement
///
/// When ng has no partial witness and counted no read out — so classes 1 and 3 are absent —
/// the per-`bases` evidence must be **equal**, field for field. That is spec §3's own
/// wording for class 2 ("summing a locus's rows must reproduce production's per-allele
/// totals") and it is what stops a split from hiding a lost read: several rows summing to
/// the wrong total would fail here where a row count never could.
#[track_caller]
fn classify_locus(
    where_: &str,
    ours: &SampleLocusObservations,
    theirs: &SampleLocusObservations,
    drops: &ProjectionDrops,
) -> (DivergenceClasses, bool) {
    assert_eq!(
        (&ours.region, &ours.reference_bases, &ours.kind),
        (&theirs.region, &theirs.reference_bases, &theirs.kind),
        "{where_}: the locus moved its anchor, its footprint or its reference bytes — none \
         of which the no-fabrication rule touches",
    );
    assert_eq!(
        drops.unsupported_chain_ids, 0,
        "{where_}: the projection dropped an unsupported bucket that still carried a chain \
         id, so a read's identity went missing on production's side of the comparison",
    );
    assert_reads_are_accounted_for(where_, ours, theirs);

    let mut classes = DivergenceClasses {
        partial_witness: ours
            .observations
            .iter()
            .any(|row| row.read_witness != ReadWitness::Complete),
        group_split: ours
            .observations
            .iter()
            .any(|row| row.read_group != PROJECTED_READ_GROUP),
        counters: ours.reads_without_observation > 0 || ours.reads_discarded_by_cap > 0,
        unsupported_bucket: drops.unsupported > 0,
        row_order: drops.reordered,
        // Filled below: unlike the other five, class 6 is not a fact about ng's rows — a
        // read production mis-folded looks, on ng's side, exactly like one it folded
        // correctly. It is recognised by its **shape**, and the shape is checked.
        stale_widen: false,
    };

    // The reconciliation, on the loci where it is meaningful. Class 1 moves the bases
    // themselves, so grouping by them cannot line the two sides up; class 3 means a read ng
    // dropped is still in production's totals; class 6 recomputes a read's bases against a
    // footprint production never revisited.
    let bases_reconcile = evidence_by_bases(ours) == evidence_by_bases(theirs);
    // Computed once and used twice: it is both what excuses the reconciliation below and
    // class 6 itself.
    let stale_widen =
        !classes.partial_witness && !bases_reconcile && stale_widen_shape(where_, ours, theirs);
    if !classes.partial_witness && !classes.counters && !stale_widen {
        assert!(
            bases_reconcile,
            "{where_}: ng's rows do not sum back to production's per-allele totals, and \
             none of a partial witness, a counted-out read or a stale widen explains it\n  \
             ng   {:?}\n  prod {:?}",
            evidence_by_bases(ours),
            evidence_by_bases(theirs),
        );
    }

    // **Chain ids, where the two rules are the same rule.** ng drops a read's id when the
    // read agreed with the reference across *everything it witnessed*; production drops it
    // when the read's bucket is `alleles[0]`. Those coincide exactly when the read witnessed
    // the whole footprint and both walkers put it in a bucket with the same bases — so that
    // is where equality is asserted, and it is a stronger claim than the subset the old
    // record-level comparison could make.
    //
    // Where the bases *do* differ the ids are part of that difference and not a separate
    // fact: a read whose witnessed window happens to match the reference carries no id on
    // ng's side while production, having folded it into a bucket its fill made non-REF,
    // keeps one. An earlier draft asserted "ng's ids are a superset of production's" here
    // and this census disproved it at `seed 0x5eed0001 case 22, locus 37` — the invariant
    // held on the *record* comparison it was written for and does not survive the turn to a
    // per-read rule.
    if !classes.partial_witness && bases_reconcile {
        assert_eq!(
            locus_chain_ids(ours),
            locus_chain_ids(theirs),
            "{where_}: the two walkers put the same reads in the same buckets and disagree \
             about which of them carry a chain id",
        );
    }

    classes.stale_widen = stale_widen;
    let exact = comparable(ours) == comparable(theirs);
    assert!(
        exact || classes.any(),
        "{where_}: ng's locus differs from the projection and none of spec §3's six \
         classes is present — an unlisted divergence, which is the one thing this census \
         may not absorb.\n  ng   {:?}\n  prod {:?}",
        comparable(ours),
        comparable(theirs),
    );
    (classes, exact)
}

/// **Every read production folded is somewhere in ng's locus** — the identity that makes
/// class 3 a claim rather than a shrug.
///
/// A read folds into exactly one bucket per record on production's side, contributing
/// exactly one observation. ng either emits it in a row or removes it and counts it in
/// `reads_without_observation` (A5) — so the two sides balance **exactly**, and any other
/// arithmetic means evidence was created or lost rather than moved.
///
/// `reads_discarded_by_cap` is deliberately *not* in this identity: the cap truncates in the
/// walk, before any record exists, so both walkers lose the same reads and neither counts
/// them in its totals. ng merely says how many, which is the whole of that counter.
#[track_caller]
fn assert_reads_are_accounted_for(
    where_: &str,
    ours: &SampleLocusObservations,
    theirs: &SampleLocusObservations,
) {
    let ours_obs: u32 = ours.observations.iter().map(|row| row.num_obs).sum();
    let theirs_obs: u32 = theirs.observations.iter().map(|row| row.num_obs).sum();
    assert_eq!(
        ours_obs + ours.reads_without_observation,
        theirs_obs,
        "{where_}: ng emitted {ours_obs} observations and counted \
         {} reads out, against production's {theirs_obs} — a read was created or lost, not \
         moved",
        ours.reads_without_observation,
    );
}

/// **Does this divergence have the shape of production's stale widen?** — class 6's
/// recogniser, and the one class that has to be recognised from the difference rather than
/// read off ng's rows.
///
/// That asymmetry is a liability and is treated as one: a class recognised from the
/// difference will absorb any *other* difference of the same shape, which is the failure
/// mode that made the old `EvidenceIntact` classifier too weak to be an anchor. So the shape
/// is made as narrow as the phenomenon allows, and every part of it is required:
///
/// 0. **The bases genuinely differ.** The caller has already established that grouping both
///    sides by `bases` does not reconcile — without that, a locus differing only in its read
///    groups would be filed here, and the class would fire everywhere.
/// 1. **The footprint spans more than one position.** A widen strictly grows a record, and a
///    one-base record has never been widened, so a one-base locus can never be in this class.
/// 2. **No evidence moved.** Every support scalar summed over the locus is identical — no
///    read was created, lost or double counted; the rows are the same reads, said
///    differently.
/// 3. **Every row production holds that ng does not is production's stale fold of an ng
///    row**: some ng row's bases share a prefix with it, and what production has past that
///    prefix is a **suffix of the reference bases** — which is precisely what appending the
///    widened tail from `ref_seq` produces. The empty suffix is allowed and is the other
///    direction of the same defect: production folded against a narrower footprint and so
///    *missed* an event ng saw, leaving its row the shorter of the two.
///
/// Anything failing any of the three is unclassified, and [`classify_locus`] panics on it.
#[track_caller]
fn stale_widen_shape(
    where_: &str,
    ours: &SampleLocusObservations,
    theirs: &SampleLocusObservations,
) -> bool {
    if ours.reference_bases.len() < 2 {
        return false;
    }
    if locus_evidence(ours) != locus_evidence(theirs) {
        return false;
    }
    let ng_bases: BTreeSet<&[u8]> = ours
        .observations
        .iter()
        .map(|row| row.bases.as_ref())
        .collect();
    for row in &theirs.observations {
        let theirs_bases: &[u8] = &row.bases;
        if ng_bases.contains(theirs_bases) {
            continue;
        }
        let explained = ng_bases.iter().any(|ours_bases| {
            let shared = ours_bases
                .iter()
                .zip(theirs_bases)
                .take_while(|(a, b)| a == b)
                .count();
            theirs_bases.len() >= shared && ours.reference_bases.ends_with(&theirs_bases[shared..])
        });
        if !explained {
            eprintln!(
                "{where_}: production's row {:?} is not any ng row's bases plus a reference \
                 tail — the stale-widen shape does not fit, so this is an unlisted class",
                String::from_utf8_lossy(theirs_bases),
            );
            return false;
        }
    }
    true
}

/// Whether any allele at this locus is carried by more than one read group — spec §3's
/// class 2, read off ng's own rows.
///
/// Two rows sharing `(bases, read_witness)` can only differ in the group, the three
/// together being the whole row identity.
fn rows_split_by_group(locus: &SampleLocusObservations) -> bool {
    /// A row's identity **without** its read group: the bases, and the coverage run as the
    /// walk's own total order gives it (`witness_order`).
    type RowIdentityWithoutGroup<'a> = (&'a [u8], (u8, u16, u16));

    let mut seen: BTreeSet<RowIdentityWithoutGroup<'_>> = BTreeSet::new();
    for row in &locus.observations {
        if !seen.insert((
            &row.bases,
            super::open_record::witness_order(row.read_witness),
        )) {
            return true;
        }
    }
    false
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

/// **The permanent anchor: at a locus where every folded read witnessed the whole footprint,
/// ng and production agree exactly — forever.**
///
/// From A2 the two walkers differ *on purpose* — production fills the positions a read did
/// not witness from the reference, and ng does not — so the full stage-1 differential stops
/// being a claim anyone can make. **What survives is narrower and is not a snapshot:** every
/// change in this plan leaves the complete class alone by construction (there are no gaps to
/// fill), so this must stay green through the rest of the plan and everything after it.
///
/// # D1 sharpened the predicate, and that is the point of this step
///
/// The previous version of this test could not state its own class. It compared
/// `PileupRecord`s — a type with no way to say "this read witnessed four of seven
/// positions" — so it ran on *every* locus of a fixture chosen to make partial witnesses
/// rare, and tolerated the ones that occurred anyway under `EvidenceIntact`: same reference
/// bytes, same support totals, *some rows' bases differ*. Its own doc said what was missing:
/// "It does not check which rows differ or by how much… the right filter is the spec's own
/// definition of the anchor class, which needs A4's `witness_of`."
///
/// That filter is now readable straight off the emitted locus. A locus qualifies when every
/// row is [`Complete`](ReadWitness::Complete) and no read was counted out — which *is*
/// "every folded read witnessed the whole footprint" — and on those loci the comparison is
/// **equality**, not evidence-preservation. The tolerated class is gone rather than
/// counted: a locus that once landed in `EvidenceIntact` now fails the predicate and is
/// measured by the census instead, where it belongs.
///
/// # And the predicate needed a second half, which is D1's other finding
///
/// "Every read witnessed the whole footprint" is **not sufficient**. Production's `widen`
/// appends reference bases to every bucket and re-folds nothing, where ng re-folds every
/// live read — so a record that widened after a read folded into it leaves production
/// holding that read's haplotype against a stale footprint, whether or not the read went on
/// to witness the whole thing. The fixture is what supplies the second half:
/// [`generate_uniform_events`] gives every read on a contig the same event set, so **no widen
/// leaves any read stale** — every read in a record carries its own copy of every widening
/// event and is re-folded by it. Records widen here as often as anywhere else, and this test
/// asserts that they do (`widens > 0`); what the shared event set removes is the *staleness*,
/// not the widen. That is the class production genuinely cannot get wrong, and it is the one
/// that must hold forever.
///
/// The loci this leaves out are not unchecked: they are the census's, where the stale widen
/// is a **named** class with its own count
/// (`every_divergence_from_production_is_one_of_the_six_named_classes`).
///
/// *Measured, so nobody mistakes which half is doing the work:* on this fixture the
/// every-read-`Complete` filter currently excludes **nothing** — all 216,203 loci qualify —
/// because uniform events leave no way for a read to be blind over part of a footprint (no
/// adaptor boundary, no `N`, and the column cap removes a read outright rather than in part).
/// So what is actually asserted here is that the two walkers agree at *every* locus of a
/// fixture built to contain no fabrication. The filter is the guard that keeps that true if
/// the fixture ever gains a partial witness; it is not the thing under test.
#[test]
fn ng_agrees_with_production_where_production_fabricated_nothing() {
    let mut compared = 0usize;
    let mut anchored = 0usize;
    let mut anchored_multi_base = 0usize;
    let mut float_only = 0usize;
    let mut widens = 0u64;
    let total = SEEDS.len() * cases_per_seed();

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate_uniform_events(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");

            let theirs = production_walk(&case);
            let ours = ng_walk(&case);
            assert_eq!(
                ours.panic_message, None,
                "{where_}: an anchor case must not panic — this fixture excludes every \
                 input that reaches production's reachable precondition"
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
            // **The fixture's own property, checked.** One CIGAR per contig is what leaves
            // no read holding a stale fold after a widen; if the generator ever stops
            // delivering that, this anchor would silently start tolerating production's
            // stale folds and would no longer be an anchor.
            for contig in 0..CONTIGS as u32 {
                let mut on_contig = case.reads.iter().filter(|read| read.chrom_id == contig);
                let Some(first) = on_contig.next() else {
                    continue;
                };
                for read in on_contig {
                    assert_eq!(
                        (&read.cigar, read.alignment_start),
                        (&first.cigar, first.alignment_start),
                        "{where_}: two reads on contig {contig} carry different event sets, \
                         so one of them can be left folded against a stale footprint and \
                         this fixture is no longer an anchor",
                    );
                }
            }
            widens += theirs
                .summary
                .as_ref()
                .expect("a walk that did not panic has a summary")
                .record_widen_events;
            for (position, (ours, theirs)) in
                ours.records.iter().zip(theirs.records.iter()).enumerate()
            {
                match (ours, theirs) {
                    (Ok(ours), Ok(theirs)) => {
                        compared += 1;
                        if !every_read_witnessed_the_whole_footprint(ours) {
                            continue;
                        }
                        anchored += 1;
                        if ours.reference_bases.len() > 1 {
                            anchored_multi_base += 1;
                        }
                        assert_eq!(
                            comparable(ours),
                            comparable(theirs),
                            "{where_}: locus {position} covers a footprint every one of \
                             its reads witnessed in full, in a walk where no record \
                             widened — so there is nothing for production to have \
                             fabricated and nothing for ng to have withheld",
                        );
                    }
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

    // The anchor is only worth the ground it covered, and **multi-base loci are the part
    // that matters**: at a one-base locus every contributor is `Complete` whatever the
    // builder does, so an anchor made only of those would stay green against a builder that
    // filled every gap it found.
    assert!(
        anchored > total * 5,
        "only {anchored} loci qualified over {total} cases — the generator has stopped \
         producing walks worth anchoring"
    );
    assert!(
        anchored_multi_base > total,
        "only {anchored_multi_base} qualifying loci span more than one base over {total} \
         cases — without them this anchor cannot tell a filling builder from a witnessing \
         one"
    );
    // **And the widen path has to be on the walk.** The class this anchor exists to protect
    // is not "records that never widened" — it is "reads no widen left stale", which is a
    // claim about widened records. A fixture that stopped producing widens would satisfy
    // every assertion above while testing the easy half.
    assert!(
        widens > 0,
        "no record widened over {total} cases — the anchor covers only records that never \
         grew, which is the half of the property that was never in doubt"
    );
    // **A ceiling on `float_only` — the quantity that *explains* the headline is the one that
    // could move silently.** `float_only` counts the loci the `q_sum` tolerance rescued: the
    // two sides differ on an exact comparison and agree within `Q_SUM_TOLERANCE`. Its stated
    // purpose is that the tolerance "is shown to be doing work rather than quietly matching
    // nothing", and until now it was printed and never asserted — so an arithmetic change
    // could take the walk from "agrees" to "agrees only because the tolerance is wide" while
    // this test reported a green anchor.
    //
    // **Measured, by injecting a uniform relative `q_sum` error into every emitted row**
    // (review, 2026-07-30): at 5 × 10⁻¹⁰ every test stayed green and this count went
    // **103 → 215,659 of 216,203**; at 2 × 10⁻⁹ three tests fail on the tolerance itself. So
    // the window where the anchor holds and the explanation is nonsense is real, and this
    // assertion closes it.
    //
    // One in ten, against a measured **103 of 216,203 (0.05 %)** at the default case count and
    // **1,014 of 1,620,856 (0.06 %)** at `PVC_PARITY_CASES=3000` — a ratio stable across a 7.5×
    // fixture, so the bound is a statement about the property and not about the case count.
    // That leaves two orders of headroom below the 99.7 % the injection produces.
    assert!(
        float_only * 10 < compared,
        "{float_only} of {compared} loci agree only within the `q_sum` tolerance — the \
         tolerance is meant to absorb summation order, which touches a fraction of a percent \
         of loci; at this rate the two walkers no longer agree on `q_sum` and the anchor is \
         passing on the tolerance rather than on the arithmetic",
    );
    eprintln!(
        "the anchor: {anchored} of {compared} loci ({:.1}%) had every read witness the whole \
         footprint, {anchored_multi_base} of them spanning more than one base, across \
         {widens} widened records; all identical to the projection, field for field. \
         {float_only} of them agree on `q_sum` only within the 1e-9 relative tolerance — the \
         order the sum accumulates in, from two causes: A3's eviction recreating a bucket the \
         reads return to, and B1 summing each read's contribution once where production \
         accumulates into the bucket with a subtract-then-add per re-fold.",
        100.0 * anchored as f64 / compared as f64,
    );
}

/// The **filter** the anchor applies on top of its fixture: *loci where every folded read
/// witnessed the whole footprint*.
///
/// **This is not the anchor's predicate** — spec §3 used to say it was, and class 6 disproves
/// it. The anchor's guarantee comes from the fixture (`generate_uniform_events`, which leaves
/// no read stale); this filter rides on top as a tripwire for fixture drift.
///
/// **And on that fixture it currently excludes nothing** — measured, 216,203 of 216,203 loci
/// qualify, because uniform events leave no way to be blind over part of a footprint. So
/// "both halves are needed" is *unverified* rather than false: neither half excludes a locus
/// today, and both would start to if the fixture ever gained a partial witness. Each half is
/// stated below so that whichever fires first is legible.
///
/// Every row `Complete` says the reads that
/// **produced** an observation each saw the whole footprint; `reads_without_observation`
/// covers the reads that saw it in pieces and so produced none at all (A5) — those have no
/// row to be `Observed`, and production folded them anyway, with the gaps filled.
///
/// `reads_discarded_by_cap` is deliberately not part of it: the cap acts in the walk, before
/// any record exists, so both walkers lose the same reads and the locus is still one where
/// everything that folded witnessed everything.
fn every_read_witnessed_the_whole_footprint(locus: &SampleLocusObservations) -> bool {
    locus.reads_without_observation == 0
        && locus
            .observations
            .iter()
            .all(|row| row.read_witness == ReadWitness::Complete)
}

/// **What A2 is allowed to change, and what it is not** — the comparison both the
/// synthetic census and the real-data one run.
///
/// Records are opened and widened from the **events**, which the no-fabrication rule does
/// not touch; only the allele buckets move. So the two streams must still be the same
/// length, at the same anchors, with the same REF bytes, the same error items and the same
/// `RunSummary`. Anything else is a bug rather than a design.
///
/// Returns the census: every locus classified, and every difference outside the six classes
/// a panic ([`classify_locus`]).
#[track_caller]
fn assert_only_allele_bytes_moved(
    where_: &str,
    ours: &WalkOutcome,
    theirs: &WalkOutcome,
    census: &mut DivergenceCensus,
) {
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
    assert_eq!(
        theirs.records.len(),
        theirs.drops.len(),
        "{where_}: the projection's drop census is a different length from the stream it \
         came from",
    );

    for (position, ((ours, theirs), drops)) in ours
        .records
        .iter()
        .zip(theirs.records.iter())
        .zip(theirs.drops.iter())
        .enumerate()
    {
        match (ours, theirs) {
            (Ok(ours), Ok(theirs)) => {
                let at = format!("{where_}: locus {position}");
                let (classes, exact) = classify_locus(&at, ours, theirs, drops);
                census.record(classes, exact);
                census.measure_fabrication(ours);
                census.measure_stale_widen(classes, ours, theirs);
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
    census.float_only += float_only_divergences(&ours.records, &theirs.records);
}

/// **Every divergence is one of the six, each is counted, and here is how big the first one
/// is** — spec §3's stage-2 assertion and §13.2's deliverable, in one pass.
///
/// The general fixture — partial reads, adaptor boundaries, `N` bases, ref-skips — is where
/// the two walkers now differ, and this is the census of that difference rather than an
/// assertion that it is absent. What makes the census trustworthy:
///
/// - **A2 changes no locus's existence and no locus's footprint.** Records are opened and
///   widened from the *events*, which this step does not touch; only the rows move. So the
///   two streams must still be the same length, at the same regions, with the same reference
///   bytes and the same `RunSummary`.
/// - **No read is created or lost**, only moved — [`assert_reads_are_accounted_for`], at
///   every locus.
/// - **A divergence outside the six classes is a panic, not a bucket.** That is the whole
///   of "counted, not excused": [`classify_locus`] reads each class off ng's own rows and
///   the projection's drops, and fails on anything left over.
/// - **Every class must fire.** A class counted zero for a whole census is a branch nothing
///   takes, and a later divergence of that kind would be triaged as one of the classes that
///   *did* fire — the exact contamination spec §3 warns about.
///
/// The three numbers reported are the ones D3 re-measures on real data.
///
/// **Six, and the sixth was found by this harness rather than by reading.** Spec §3's table
/// listed five; building the anchor produced a locus none of them explained, and §3 now lists
/// six with the mechanism. It is **not** "production's `widen` re-folds nothing" — production
/// re-folds every *contributor* into an affected record. What it leaves unrevised is the
/// appended reference tail on a read that was **not** a contributor at the widening step. See
/// [`DivergenceClasses::stale_widen`].
/// # Two passes, because one read group and two prove different things
///
/// At **one** read group class 2 cannot fire, so the "no unlisted divergence" panic has its
/// full force: any difference at all has to be explained by one of the other five. At
/// **two**, class 2 fires almost everywhere — the group is part of every row's identity —
/// which exercises the class and the split it exists for but leaves the unlisted-divergence
/// check with nothing to catch. Both passes run over the same cases, and each asserts what
/// it is in a position to assert: the one-group pass that class 2 is *silent*, the two-group
/// pass that it is not.
#[test]
fn every_divergence_from_production_is_one_of_the_six_named_classes() {
    let mut one_group = DivergenceCensus::default();
    let mut two_groups = DivergenceCensus::default();

    // Both generators: the general one reaches the classes that need partial witnesses, the
    // complete-reads one is where a divergence must be the widen or nothing, so running it
    // here is what keeps class 6 honest.
    for generator in [
        generate as fn(&mut SplitMix64) -> Case,
        generate_complete as fn(&mut SplitMix64) -> Case,
    ] {
        for seed in SEEDS {
            let mut rng = SplitMix64(seed);
            for index in 0..cases_per_seed() {
                let case = generator(&mut rng);
                let where_ = format!("seed {seed:#x} case {index}");

                let theirs = production_walk(&case);
                for (groups, census) in [(1, &mut one_group), (2, &mut two_groups)] {
                    // The read group reaches nothing but the row key — see
                    // `ng_walk_in_groups` — so the walk itself is the same walk both times.
                    let ours = ng_walk_in_groups(&case, groups);
                    assert_only_allele_bytes_moved(&where_, &ours, &theirs, census);
                }
            }
        }
    }

    let census = &one_group;
    let diverged = census.loci - census.exact;
    assert!(
        diverged > 0,
        "{} loci compared and not one differed — the fill is still there, and every anchor \
         assertion in this module is passing for the wrong reason",
        census.loci,
    );
    // **Each class, present.** A class counted zero is a branch nothing takes.
    for (class, count) in [
        ("1 (a partial witness)", census.partial_witness),
        ("3 (a per-locus counter)", census.counters),
        ("4 (an unsupported bucket)", census.unsupported_bucket),
        ("5 (row order)", census.row_order),
        ("6 (production's stale widen)", census.stale_widen),
        // Class 2's own two counts, from the pass that can produce them.
        (
            "2 (a read group production cannot say)",
            two_groups.group_split,
        ),
        // Not a class: the consequence class 2 exists for. Spec §13 asks for it by name —
        // "an allele supported by both groups is two rows" — and the class flag fires on the
        // weaker event of a row merely naming a different group, so this would go to zero
        // unnoticed if the split ever stopped happening.
        (
            "2's split, an allele in two rows",
            two_groups.group_split_rows,
        ),
    ] {
        assert!(
            count > 0,
            "class {class} never fired over {} loci — it is a branch nothing takes, and the \
             next divergence of that kind will be triaged as one of the classes that did",
            census.loci,
        );
    }
    // **And silent where it cannot be true.** One read group means every row carries the
    // group `project` invents, so class 2 must count zero — without this, a bug that tagged
    // rows with an arbitrary group would look like the class working.
    assert_eq!(
        census.group_split, 0,
        "class 2 fired on a one-read-group walk, where every row's group is the one the \
         projection uses — so the class is firing on something other than the group",
    );
    // **A ceiling as well as a floor on the deliverable.** The floor alone makes this a
    // census that cannot fail upwards: the headline can be driven from 1.9 % to 91.8 % by a
    // genuine defect and the test still passes, reporting the defect as the measurement.
    // "Production fabricated bases here" is a minority of loci by construction — most
    // positions carry no widened record at all — so a figure in the tens of percent means
    // the walk changed, not that the defect is larger than thought.
    assert!(
        census.fabricating_loci * 10 < census.loci,
        "{} of {} loci ({:.1} %) carry a partial witness — this census measures production's \
         fabrication, which is a small minority of loci; at this rate something in the walk \
         has moved and the number is no longer the defect's size",
        census.fabricating_loci,
        census.loci,
        100.0 * census.fabricating_loci as f64 / census.loci as f64,
    );
    // **The stale-widen triple needs its own floor, and for a sharper reason than symmetry.**
    // The class flag is asserted non-zero below with the other five, so `stale_widen` cannot
    // silently die. These two can: they are read off production's *rows* rather than off the
    // classification, so a change that made every class-6 locus report a zero-length tail —
    // or made `measure_stale_widen` skip every row — would leave the class count intact and
    // the two numbers at zero, which reads as "production mis-folds reads but appends nothing
    // on their behalf". That is not a state the class can be in: a locus is class 6 *because*
    // a production row is an ng row plus a reference tail, so the tail has length.
    // Mutation-verified: making `measure_stale_widen` return without measuring reports
    // "264 loci, 0 reads, 0 bases" and fails here.
    //
    // **Two tighter assertions were tried and rejected as unsound**, recorded so nobody adds
    // them later on the strength of the measured figures (267 reads and 544 bases over 264
    // loci, so both ratios look safe):
    //
    // - `stale_widen_reads >= stale_widen` — one mis-folded read per class-6 locus — fails on a
    //   locus where every production row *does* appear among ng's and `!bases_reconcile` comes
    //   from the counts on matching bases instead. `stale_widen_shape`'s loop body never fires
    //   there, so it returns `true` having measured nothing.
    // - `stale_widen_ref_bases >= stale_widen_reads` — one appended base per mis-folded read —
    //   fails when production's row is a strict **prefix** of an ng row: `shared` is then the
    //   whole production row, the tail is empty, and `ends_with(&[])` accepts it. That is not a
    //   corner case but the shape of D1's original counter-example, where production held eight
    //   bases against ng's nine.
    // **The same ceiling on `float_only` as the anchor carries, and for the same reason.**
    // Measured here at **863 of 256,974 loci (0.34 %)**, and **6,666 of 1,957,023 (0.34 %)** at
    // `PVC_PARITY_CASES=3000` — the same fraction over a 7.5× fixture. The review's uniform
    // 5 × 10⁻¹⁰ relative `q_sum` injection takes it to **253,143** while every class count and
    // the headline hold. `float_only` is a subset of the **exact** count, not of `diverged`:
    // `classify_locus` decides `exact` with `comparable`, which applies the tolerance, so a
    // locus the tolerance rescued is reported as identical. That is what makes it able to
    // move silently — nothing else in this census would change.
    assert!(
        census.float_only * 10 < census.loci,
        "{} of {} loci agree only within the `q_sum` tolerance — the tolerance absorbs \
         summation order, which touches a fraction of a percent of loci; at this rate the \
         census's \"identical to the projection\" is carried by the tolerance rather than by \
         the arithmetic",
        census.float_only,
        census.loci,
    );
    assert!(
        census.stale_widen_reads > 0 && census.stale_widen_ref_bases > 0,
        "class 6 fired on {} loci but the stale-widen deliverable is {} reads / {} bases — a \
         locus is class 6 because production holds a row that is an ng row plus a reference \
         tail, so both numbers are positive whenever the class is",
        census.stale_widen,
        census.stale_widen_reads,
        census.stale_widen_ref_bases,
    );
    eprintln!(
        "the divergence census over {} loci (one read group): {} identical to the \
         projection, {diverged} differing, {} of the identical ones agreeing on `q_sum` only \
         within the 1e-9 relative tolerance.\n  \
         class 1 partial witness {}   class 3 counters {}   class 4 unsupported bucket {}   \
         class 5 row order {}   class 6 stale widen {}\n  \
         at two read groups: class 2 on {} loci, of which {} carry one allele in two rows.\n  \
         the deliverable: production credited {} reads over {} loci with {} reference bases \
         they never sequenced ({:.2} fabricated bases per fabricating locus).\n  \
         the stale-widen deliverable (§13.2's second triple): production mis-folded {} reads \
         over {} loci, appending {} reference bases on their behalf.",
        census.loci,
        census.exact,
        census.float_only,
        census.partial_witness,
        census.counters,
        census.unsupported_bucket,
        census.row_order,
        census.stale_widen,
        two_groups.group_split,
        two_groups.group_split_rows,
        census.fabricated_reads,
        census.fabricating_loci,
        census.fabricated_ref_bases,
        census.fabricated_ref_bases as f64 / census.fabricating_loci.max(1) as f64,
        census.stale_widen_reads,
        census.stale_widen,
        census.stale_widen_ref_bases,
    );
}

/// **The projection lays production's record out the way the walk lays out ng's** — the
/// promise [`sort_rows`] makes, checked rather than commented.
///
/// `finalise` sorts `ObservationRow`s and [`sort_rows`] sorts `SequenceObservation`s, so the two
/// loops cannot be shared even though the `ReadWitness` comparator is (`witness_order`,
/// lifted to `pub(super)` for this). What could still drift is the *rest* of the key — the
/// bases, then the group — so this walks a fixture and asserts that sorting ng's own emitted
/// rows with this function does not move them. If either spelling of the order changes, the
/// projection stops being comparable to the walk and spec §3's class 5 stops being a
/// normalisation.
///
/// The fixture is the general one at two read groups, because a one-group walk cannot
/// distinguish an order that considers the group from one that does not.
#[test]
fn the_projection_orders_rows_as_the_walk_does() {
    let mut loci = 0usize;
    let mut multi_row = 0usize;
    let mut grouped = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");
            for (position, item) in ng_walk_in_groups(&case, 2).records.iter().enumerate() {
                let Ok(locus) = item else { continue };
                loci += 1;
                if locus.observations.len() > 1 {
                    multi_row += 1;
                }
                if rows_split_by_group(locus) {
                    grouped += 1;
                }
                let mut sorted = locus.observations.clone();
                sort_rows(&mut sorted);
                assert_eq!(
                    sorted, locus.observations,
                    "{where_}: locus {position} came out of the walk in an order this \
                     projection would not have produced, so the two orders have drifted",
                );
            }
        }
    }

    assert!(
        multi_row * 10 > loci && grouped > 0,
        "only {multi_row} of {loci} loci carry more than one row and {grouped} split by read \
         group — an order test needs rows to order"
    );
}

/// **Every field of a `PileupRecord` arrives somewhere, or is dropped by name** — the
/// projection's own oracle, on a record built by hand so the expected locus can be written
/// out in full.
///
/// The differential can only ever say "the two agree"; it cannot say the projection is
/// *right*. If `project` mapped `fwd` onto `placed_left`, both sides of a comparison that
/// went through it would be equally wrong and the census would be green — the A-versus-B
/// blindness this branch has already been caught by once.
#[test]
fn the_projection_says_everything_a_record_says() {
    use crate::pileup_record::AlleleObservation;

    let record = PileupRecord {
        chrom_id: 3,
        pos: 940,
        alleles: vec![
            AlleleObservation {
                seq: b"ACGT".to_vec(),
                support: AlleleSupportStats {
                    num_obs: 7,
                    q_sum: -12.5,
                    fwd: 3,
                    placed_left: 2,
                    // Dropped by name: ng stops computing it at B2 (spec §6).
                    placed_start: 5,
                    mapq_sum: 210,
                    mapq_sum_sq: 6_300,
                },
                // Production drops ids from `alleles[0]`; an id here would be a bug on its
                // side, so the REF row carries none and the projection must not invent any.
                chain_ids: Vec::new(),
            },
            // A supported alt that sorts **before** the REF bucket, so production's
            // creation order is genuinely not ng's emission order and class 5 fires on the
            // rows the projection keeps. Its ids are out of order and duplicated — the
            // projection sorts and dedups, as `finalise` does.
            AlleleObservation {
                seq: b"AAGT".to_vec(),
                support: AlleleSupportStats {
                    num_obs: 2,
                    q_sum: -4.0,
                    fwd: 1,
                    placed_left: 1,
                    placed_start: 1,
                    mapq_sum: 100,
                    mapq_sum_sq: 5_000,
                },
                chain_ids: vec![9, 4, 9],
            },
            // Class 4: a bucket A3's widen stranded. Dropped — and it sorts *between* the
            // two kept rows, so a projection that kept it would fail the row list as well
            // as the drop count.
            AlleleObservation {
                seq: b"GCGT".to_vec(),
                support: AlleleSupportStats::default(),
                chain_ids: Vec::new(),
            },
        ],
        windowed_gc: f32::NAN,
        windowed_coverage: f32::NAN,
    };

    let (locus, drops) = project_counting_drops(&record);

    assert_eq!(
        locus.region,
        GenomeRegion {
            contig: ContigId(3),
            start: Position(940),
            // 1-based inclusive over a four-base REF bucket: 940, 941, 942, 943.
            end: Position(943),
        },
        "the region's end is the only statement production makes about the footprint's \
         extent, and it is implicit in the REF bucket's length"
    );
    assert_eq!(&*locus.reference_bases, b"ACGT");
    assert_eq!(locus.kind, LocusKind::Generic);
    assert_eq!(
        (
            locus.reads_without_observation,
            locus.reads_discarded_by_cap
        ),
        (0, 0)
    );
    assert_eq!(
        drops,
        ProjectionDrops {
            unsupported: 1,
            unsupported_chain_ids: 0,
            // The rows the projection *keeps* come out of production as `ACGT`, `AAGT` —
            // not sorted, so class 5 fired on this record. (Measured on the kept rows
            // deliberately: dropping a bucket that happened to be out of place is class 4,
            // and counting it here would make the two classes indistinguishable.)
            reordered: true,
        }
    );
    assert_eq!(
        locus.observations,
        vec![
            SequenceObservation {
                bases: Box::from(&b"AAGT"[..]),
                read_witness: ReadWitness::Complete,
                read_group: PROJECTED_READ_GROUP,
                num_obs: 2,
                num_fwd: 1,
                q_sum: -4.0,
                mapq_sum: 100,
                mapq_sum_sq: 5_000,
                placed_left: 1,
                chain_ids: vec![4, 9],
            },
            SequenceObservation {
                bases: Box::from(&b"ACGT"[..]),
                read_witness: ReadWitness::Complete,
                read_group: PROJECTED_READ_GROUP,
                num_obs: 7,
                num_fwd: 3,
                q_sum: -12.5,
                mapq_sum: 210,
                mapq_sum_sq: 6_300,
                placed_left: 2,
                chain_ids: Vec::new(),
            },
        ],
        "every support scalar must arrive on its own field, `fwd` as `num_fwd` and \
         `placed_left` as `placed_left` — a projection that crossed the two would leave both \
         sides of the differential equally wrong"
    );
    // The entry point the plan names, and the one the differential's other users call.
    assert_eq!(project(&record), locus);
}

/// **Every emitted row is a row some read is folded into** — which is *not* A3's eviction,
/// and D1 is where the difference was found.
///
/// This test used to claim the eviction. It ran on `PileupRecord`s and asserted that no
/// emitted bucket had `num_obs == 0`, and while the walk emitted production's type that was a
/// real check. **B2 made it vacuous:** ng's rows are derived from `folded_reads`, per read, so
/// a bucket no read is folded into produces no row at all and leaves no trace in the output.
/// D1 mutated the code the test named — moving `evict_unsupported_alleles` above the
/// contributor fold loop, which strands every bucket that loop empties — and all 198 tests
/// in this module stayed green. The eviction is now enforced by a `debug_assert!` in
/// `finalise`, where the buckets still exist and where every walk in the suite checks it.
///
/// What is left here is worth keeping and is a different claim: **`num_obs` is never zero on
/// an emitted row**, i.e. the re-derivation cannot mint a row for nothing. That would fail if
/// `observation_rows` ever created rows from the bucket table instead of from the reads —
/// which is precisely the B1 mistake it was written to avoid.
#[test]
fn every_emitted_row_carries_a_read() {
    let mut records = 0usize;
    let mut multi_allele = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");
            for (position, item) in ng_walk(&case).records.iter().enumerate() {
                let Ok(locus) = item else { continue };
                records += 1;
                if locus.observations.len() > 2 {
                    multi_allele += 1;
                }
                for (index, row) in locus.observations.iter().enumerate() {
                    assert!(
                        row.num_obs > 0,
                        "{where_}: emitted locus {position} carries row {index} \
                         ({:?}) with no supporting read: {locus:?}",
                        String::from_utf8_lossy(&row.bases),
                    );
                }
            }
        }
    }

    // Loci where the reads disagreed, so a read has somewhere to move *to* — without them
    // the fixture never exercises a re-fold at all.
    assert!(
        multi_allele * 100 > records,
        "only {multi_allele} of {records} emitted loci carry more than two rows — this \
         fixture cannot exercise a read moving between buckets, so the property it asserts \
         is vacuous"
    );
    eprintln!(
        "{records} emitted loci, {multi_allele} of them with more than two rows; every row \
         carried at least one read"
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
    // **The reference bytes first, then the rows.** On ng's type the locus's own reference
    // sequence is a field rather than `alleles[0]`, so the two are read separately — which
    // is why the "REF bucket" assertions below name `reference_bases`.
    let at_19 = |outcome: &WalkOutcome| {
        outcome
            .records
            .iter()
            .filter_map(|item| item.as_ref().ok())
            .find(|locus| locus.region.start.get() == 19)
            .map(|locus| {
                let reference = String::from_utf8_lossy(&locus.reference_bases).to_string();
                let rows: Vec<String> = locus
                    .observations
                    .iter()
                    .map(|row| String::from_utf8_lossy(&row.bases).to_string())
                    .collect();
                (reference, rows)
            })
            .expect("pair9's deletion opens a record at 19")
    };
    let ((theirs_reference, theirs_at_19), (ours_reference, ours_at_19)) =
        (at_19(&theirs), at_19(&ours));
    assert_eq!(
        theirs_reference.len(),
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
    let before_the_fix = format!("{}AAA", &theirs_reference[..1]);
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
    // Everything else about the record is untouched: same reference bytes, same anchor.
    assert_eq!(
        theirs_reference, ours_reference,
        "the locus's reference bytes are the record's own, and A2 does not touch them",
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
                    Ok(locus) => {
                        if locus.reference_bases.len() > 1 {
                            multi_base_records += 1;
                        }
                        // On the projection a row *is* a supported allele, so "more than
                        // one" counts loci where the reads disagreed — which is what this
                        // number always meant, production's unsupported `alleles[0]` having
                        // been the reason the old spelling needed the `> 1`.
                        if locus.observations.len() > 1 {
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
    let mut census = DivergenceCensus::default();
    assert_only_allele_bytes_moved(&where_, &ours, &theirs, &mut census);
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
    // **The three numbers spec §13.2 calls the deliverable**, on real alignments — the same
    // census the synthetic differential prints, over the same code path, so the two cannot
    // measure subtly different quantities. D3 is what runs this at chromosome scale.
    eprintln!(
        "real-data differential: {where_} — {} loci compared from {prepared_reads} prepared \
         reads; every region, reference sequence and counter identical.\n  \
         class 1 partial witness {}   class 2 group split {}   class 3 counters {}   \
         class 4 unsupported bucket {}   class 5 row order {}   class 6 stale widen {}\n  \
         the deliverable: production credited {} reads over {} loci ({:.2} % of them) with \
         {} reference bases they never sequenced.\n  \
         the stale-widen deliverable (§13.2's second triple): production mis-folded {} reads \
         over {} loci, appending {} reference bases on their behalf.",
        census.loci,
        census.partial_witness,
        census.group_split,
        census.counters,
        census.unsupported_bucket,
        census.row_order,
        census.stale_widen,
        census.fabricated_reads,
        census.fabricating_loci,
        100.0 * census.fabricating_loci as f64 / census.loci.max(1) as f64,
        census.fabricated_ref_bases,
        census.stale_widen_reads,
        census.stale_widen,
        census.stale_widen_ref_bases,
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
