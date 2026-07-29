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
//! that is equivalent to sharing one; the real-data test lends a single `RefSeqFetcher`
//! to both because *its* accessor is not stateless.
//!
//! # It has been shown to fail — B2, re-run 2026-07-29
//!
//! A differential that has only ever passed is a claim. Each behaviour the plan names was
//! mutated **in ng's copy**, one at a time, and `ng_walks_identically_to_production` was
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

use std::sync::Arc;

use super::tests::MockFasta;
use super::{PreparedRead, WalkerConfig};
use crate::ng::read::PLACEHOLDER_READ_GROUP;
// Aliased, so which walker a call reaches is legible at the call site rather than carried
// by a `super::` — this file is the one place both are in scope at once.
use crate::pileup::walker::{
    CigarOp, MateRole as ProductionMateRole, PreparedRead as ProductionPreparedRead,
    run as production_run,
};
use crate::pileup_record::PileupRecord;

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
/// Errors are rendered with `{:?}` rather than `to_string()`: the two `WalkerError` types
/// are nominally distinct so they cannot be compared directly, and `Debug` shows the
/// variant and every field where `Display` shows only what the format string chose — an
/// `io::ErrorKind` inside `WalkerError::Fasta` is invisible through `Display`.
fn drive<W, E>(mut walker: W, summary_of: impl FnOnce(&W) -> SummaryCounters) -> WalkOutcome
where
    W: Iterator<Item = Result<PileupRecord, E>>,
    E: std::fmt::Debug,
{
    let mut records = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(item.map_err(|error| format!("{error:?}")));
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
    drive(
        production_run(case.reads.clone(), &fasta, &case.config),
        |walker| production_counters(walker.summary()),
    )
}

/// ng's answer, over the same reads converted to ng's read type.
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
    drive(super::run(reads, &fasta, &case.config), |walker| {
        ng_counters(walker.summary())
    })
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
        assert_eq!(ours, theirs, "{where_}: stream item {position} diverged");
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
/// the convention `delimit_parity` set: `PVC_PARITY_CASES=20000 cargo test
/// ng_walks_identically_to_production`.
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

/// **The port anchor.** Every record production emits, ng emits identically — and every
/// error, and every counter.
///
/// A failure means the copy has drifted from the walker it was transcribed from; the seed
/// and case index in the message replay the exact input.
#[test]
fn ng_walks_identically_to_production() {
    let mut compared_records = 0usize;
    let mut panicking_cases = 0usize;
    let total = SEEDS.len() * cases_per_seed();

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for index in 0..cases_per_seed() {
            let case = generate(&mut rng);
            let where_ = format!("seed {seed:#x} case {index}");

            let theirs = production_walk(&case);
            let ours = ng_walk(&case);
            compared_records += ours.records.len();
            if assert_same_walk(&where_, &ours, &theirs) {
                panicking_cases += 1;
            }
        }
    }

    // The differential is only worth the ground it covered. A generator change that
    // quietly stopped emitting records would leave every assertion above vacuous.
    assert!(
        compared_records > total * 5,
        "only {compared_records} records compared over {total} cases — the generator has \
         stopped producing walks worth comparing"
    );
    // Reported, not asserted at a threshold. It was 4.2% before the production defect this
    // harness found was fixed; it is 0 now. Deliberately **not** asserted to be zero
    // either: the channel's job is to compare *how* the two walkers stop, not to claim
    // nothing can stop them.
    eprintln!(
        "stage-1 differential: {compared_records} records compared over {total} cases; \
         {panicking_cases} cases ({:.1}%) reached production's reachable debug_assert! and \
         panicked in both walkers",
        100.0 * panicking_cases as f64 / total as f64,
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
    assert_eq!(ours.records, theirs.records, "and the records must agree");

    // The record at 19 is the one the defect corrupted: `pair9/First`'s deletion opens it
    // spanning 19..=25, and `pair3/Second` folds in carrying a deletion anchored at 17
    // whose run covers 18–22 — so of this record's seven positions it witnessed only
    // 23, 24, 25.
    let record = theirs
        .records
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .find(|record| record.pos == 19)
        .expect("pair9's deletion opens a record at 19");
    let reference = &record.alleles[0].seq;
    assert_eq!(
        reference.len(),
        7,
        "the record spans 19..=25: the anchor plus six deleted bases"
    );

    // **The wrong answer, named.** Every fixture read's `seq` is all `A`, so
    // `pair3/Second`'s honest contribution to this record is its own three bases at 23, 24
    // and 25 — the positions past the end of its own deletion — and nothing before them.
    //
    // Before the fix the saturated offset made the fold emit `ref_seq[0]` first: the base
    // at **19**, which this read had deleted and never sequenced. So the allele was one
    // base longer and began with a reference base borrowed from a position the read
    // explicitly says is absent. Both spellings are written out here, because "the right
    // bases are present" and "the wrong bases are gone" are different claims and a
    // regression could satisfy either alone.
    let witnessed: &[u8] = b"AAA";
    let before_the_fix: Vec<u8> = [&reference[..1], witnessed].concat();
    let folded: Vec<String> = record
        .alleles
        .iter()
        .map(|allele| String::from_utf8_lossy(&allele.seq).to_string())
        .collect();

    assert!(
        record.alleles.iter().any(|allele| allele.seq == witnessed),
        "the read whose deletion covers 19–22 should contribute the three bases it \
         witnessed ({}), but the record holds {folded:?}",
        String::from_utf8_lossy(witnessed),
    );
    assert!(
        !record
            .alleles
            .iter()
            .any(|allele| allele.seq == before_the_fix),
        "the record still holds {}, the pre-fix spelling: a leading base at 19 that this \
         read had deleted. Records: {folded:?}",
        String::from_utf8_lossy(&before_the_fix),
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
/// leaving `ng_walks_identically_to_production` quietly weaker.
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
/// `ng_walks_identically_to_production` is built to keep reads in bounds and in order — its
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
    let fixtures: Vec<(&str, Vec<ProductionPreparedRead>)> = vec![
        (
            "out of order",
            vec![
                read("a", 20, 27, vec![CigarOp::Match(8)], 8),
                read("b", 4, 11, vec![CigarOp::Match(8)], 8),
            ],
        ),
        (
            // The check is `alignment_end < alignment_start`, not "the CIGAR consumes no
            // reference" — an all-insertion read whose `alignment_end` equals its start
            // sails through, which the fixture's own reach assertion caught.
            "zero reference span",
            vec![read("i", 4, 3, vec![CigarOp::Insertion(4)], 4)],
        ),
        (
            "cigar consumes more read bases than seq provides",
            vec![read("m", 4, 11, vec![CigarOp::Match(8)], 5)],
        ),
        (
            "seq and bq of different lengths",
            vec![{
                let mut malformed = read("q", 4, 11, vec![CigarOp::Match(8)], 8);
                malformed.bq_baq.truncate(7);
                malformed
            }],
        ),
    ];

    for (name, reads) in fixtures {
        let case = Case {
            reference: reference.clone(),
            reads,
            config: WalkerConfig::default(),
            reads_with_live_adaptor_boundary: 0,
        };
        let theirs = production_walk(&case);
        let ours = ng_walk(&case);
        assert!(
            theirs.records.iter().any(|item| item.is_err()),
            "{name}: production emitted no error, so this fixture tests nothing"
        );
        assert_same_walk(name, &ours, &theirs);
    }
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
///   cargo test --release --lib ng_walks_identically_to_production_on_real_reads \
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
/// A single `RefSeqFetcher` is lent to both walkers for the same reason: identical bytes
/// by construction, so any divergence is the walk's.
#[test]
#[ignore = "needs a real BAM/CRAM and reference; see the doc comment for the invocation"]
fn ng_walks_identically_to_production_on_real_reads() {
    use std::path::PathBuf;

    use super::RefSeqFetcher;
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

    // One fetcher, lent to both walks: identical bytes by construction.
    let fetcher = RefSeqFetcher(WindowedRefSeq::new(fasta.clone(), contigs.clone()));
    let config = WalkerConfig::default();

    let theirs = drive(
        production_run(production_reads, &fetcher, &config),
        |walker| production_counters(walker.summary()),
    );
    let ours = drive(super::run(ng_reads, &fetcher, &config), |walker| {
        ng_counters(walker.summary())
    });

    let where_ = format!("{reads_path} {region:?}");
    let records_compared = theirs.records.len();
    let ok_records = theirs.records.iter().filter(|item| item.is_ok()).count();
    let first_error = theirs
        .records
        .iter()
        .find_map(|item| item.as_ref().err().cloned());
    let panicked = assert_same_walk(&where_, &ours, &theirs);

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
         {prepared_reads} prepared reads, zero divergences"
    );
}
