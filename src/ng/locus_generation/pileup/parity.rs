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
//! The reference is one `MockFasta`, lent to both.
//!
//! # It has been shown to fail — B2, 2026-07-29
//!
//! A differential that has only ever passed is a claim. Each of the five behaviours the
//! plan names was mutated **in ng's copy**, one at a time, and
//! `ng_walks_identically_to_production` was required to fail; then the mutation was
//! reverted and the differential re-run green. All five died inside the **first six
//! cases of the first seed**, against a default of 400 cases × 4 seeds — so the margin
//! is roughly two orders of magnitude, not a coin flip.
//!
//! | # | behaviour | mutation applied to ng's copy | first divergence |
//! |---|---|---|---|
//! | 1 | mate-overlap reconciliation | early `return` from `genome_walk::resolve_mate_overlap_at_pos` | seed 0 case 0, item 17 |
//! | 2 | adaptor masking | `cigar_cursor::base_in_adaptor` always `false` | seed 0 case 0 — 29 records against production's 24 |
//! | 3 | record widening | `open_record::widen` extends only `alleles[0]`, not every bucket | seed 0 case 5, item 5 |
//! | 4 | the subtract-then-add re-fold | the `subtract_contribution` half dropped, so a re-folding read double-counts | seed 0 case 3, item 8 |
//! | 5 | the column depth cap | the `truncate(cap)` removed, counter left incrementing | seed 0 case 1, item 2 |
//!
//! Mutation 5 is the one worth noting: it leaves `column_depth_truncations` incrementing,
//! so a differential that compared only the summary would have passed it. The records
//! caught it.
//!
//! # Run at scale — B3, 2026-07-29. **Zero divergences everywhere.**
//!
//! | run | scale | result |
//! |---|---|---|
//! | synthetic, release, `PVC_PARITY_CASES=5000` | 20,000 cases | **968,852 records, 0 divergences** |
//! | synthetic, debug, `PVC_PARITY_CASES=2500` | 10,000 cases | **469,069 records, 0 divergences**; 457 cases (4.6%) panicked — in *both* walkers |
//! | GIAB HG002 10×, `chr1:1000000-1400000` | targeted TR bundle | **4,600 records, 0 divergences** |
//! | GIAB HG002 300×, `chr1:100000000-120000000` | 20 Mb | **137,591 records, 0 divergences** |
//! | tomato CRAM `SRR7279481.p1`, `SL4.0ch01:3406886-3506886` | 100 kb | **96,260 records, 0 divergences** |
//! | tomato CRAM `SRR7279481.p1`, `SL4.0ch01:13806669-15092603` | 1.3 Mb | **198,673 records, 0 divergences** |
//!
//! 437,124 records of real sequencing data across two organisms, a BAM and a CRAM, 10× and
//! 300×; 1.4 M more synthetic. The debug row is the one that says something the others
//! cannot: at scale the two walkers agree not only on outputs but on **which inputs reach
//! production's reachable `debug_assert!`** — 457 cases, both walkers, every time.
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
use crate::ng::types::ReadGroupId;
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
        cigar.push(CigarOp::Match(matched as u32));
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
    if !cigar.iter().any(|op| {
        matches!(
            op,
            CigarOp::Match(_) | CigarOp::Deletion(_) | CigarOp::Skip(_)
        )
    }) {
        cigar.push(CigarOp::Match(1));
        seq.push(contig[ref_pos.min(contig.len() - 1)]);
        ref_pos += 1;
    }

    let alignment_end = ref_pos as u32; // 1-based inclusive == 0-based exclusive
    let is_reverse_strand = rng.one_in(2);

    // The adaptor boundary. Placed *inside* the read's own span on most of the reads that
    // get one, because a boundary outside it silences nothing and would make the G1
    // filter a no-op the differential could not see.
    let (adaptor_boundary, live_boundary) = if rng.one_in(4) {
        let span = alignment_end.saturating_sub(start);
        if span >= 2 {
            let boundary = start + 1 + (rng.below(span as usize - 1) as u32);
            (Some(boundary), true)
        } else {
            (Some(start), true)
        }
    } else {
        (None, false)
    };

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
        let start = 1 + rng.below(CONTIG_LENGTH / 3) as u32;

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
    if rng.one_in(8) {
        config.mate_lookup_window = 1 + rng.below(40) as u32;
    }

    Case {
        reference,
        reads,
        config,
        reads_with_live_adaptor_boundary,
    }
}

/// What one walk produced: its record/error stream, its summary, and whether it **panicked**.
///
/// The panic flag is not defensiveness. Production's `apply_events_to_ref_into` carries a
/// `debug_assert!` that every event's anchor is at or after the record's own position, and
/// that precondition is **reachable on a legal read stream** — see
/// `both_walkers_panic_on_a_deletion_anchored_before_its_record` for the three-read case
/// and the mechanism. So on those inputs a debug build panics before either walker can
/// finish, and a harness that could not represent that would have to exclude exactly the
/// long-deletion inputs this port exists to get right.
///
/// Representing it instead makes the parity claim *stronger*: the two walkers must agree
/// on **which inputs they panic on**, and on every record emitted before the panic. A
/// verbatim copy must panic verbatim.
struct WalkOutcome {
    records: Vec<Result<PileupRecord, String>>,
    /// `None` when the walk panicked — there is no summary to read from a walker that
    /// unwound out from under us.
    summary: Option<[u64; 8]>,
    panicked: bool,
}

/// Flatten a `RunSummary` — production's or ng's copy, which are structurally identical
/// and nominally distinct — into the array `SUMMARY_FIELDS` names.
macro_rules! summary_array {
    ($summary:expr) => {{
        let summary = $summary;
        [
            summary.reads_admitted,
            summary.records_emitted,
            summary.record_widen_events,
            summary.mate_overlap_positions,
            summary.chain_allocations,
            u64::from(summary.active_reads_high_water),
            summary.mate_lookup_evictions,
            summary.column_depth_truncations,
        ]
    }};
}

/// Production's answer.
fn production_walk(case: &Case) -> WalkOutcome {
    let fasta = case.fasta();
    let mut records = Vec::new();
    let mut summary = None;
    // `records` is filled inside the closure and survives an unwind, so a panicking walk
    // still yields everything it emitted first — which is what makes the comparison
    // element-wise rather than all-or-nothing.
    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut walker = production_run(case.reads.clone(), &fasta, &case.config);
        for item in &mut walker {
            records.push(item.map_err(|error| error.to_string()));
        }
        summary = Some(summary_array!(walker.summary()));
    }))
    .is_err();
    WalkOutcome {
        records,
        summary,
        panicked,
    }
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
        .map(|read| PreparedRead::from_production(read, ReadGroupId(0)))
        .collect();
    let mut records = Vec::new();
    let mut summary = None;
    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut walker = super::run(reads, &fasta, &case.config);
        for item in &mut walker {
            records.push(item.map_err(|error| error.to_string()));
        }
        summary = Some(summary_array!(walker.summary()));
    }))
    .is_err();
    WalkOutcome {
        records,
        summary,
        panicked,
    }
}

/// The `RunSummary` field names, in the order `WalkOutcome::summary` stores them, so a
/// divergence names the counter instead of an index.
const SUMMARY_FIELDS: [&str; 8] = [
    "reads_admitted",
    "records_emitted",
    "record_widen_events",
    "mate_overlap_positions",
    "chain_allocations",
    "active_reads_high_water",
    "mate_lookup_evictions",
    "column_depth_truncations",
];

/// Assert the two walks agree, and return whether they agreed *by both panicking*.
///
/// Shared by the synthetic differential and the real-data one, so the two cannot drift
/// into asserting different things about the same claim.
#[track_caller]
fn assert_same_walk(where_: &str, ours: &WalkOutcome, theirs: &WalkOutcome) -> bool {
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

    // A verbatim copy must panic verbatim: the two must agree on *which* inputs reach
    // production's reachable `debug_assert!`, not merely on the outputs of the inputs
    // that do not.
    assert_eq!(
        ours.panicked, theirs.panicked,
        "{where_}: ng panicked={} but production panicked={}",
        ours.panicked, theirs.panicked,
    );
    if ours.panicked {
        return true;
    }

    let (ours, theirs) = (
        ours.summary
            .expect("a walk that did not panic has a summary"),
        theirs
            .summary
            .expect("a walk that did not panic has a summary"),
    );
    for (field, (ours, theirs)) in ours.iter().zip(theirs.iter()).enumerate() {
        assert_eq!(
            ours, theirs,
            "{where_}: RunSummary::{} diverged",
            SUMMARY_FIELDS[field],
        );
    }
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
    // Reported, not asserted at a threshold: this is the size of the production defect
    // this harness found, and it is expected to be a small minority of cases. Asserting
    // it is non-zero would make the differential fail the day production is fixed.
    eprintln!(
        "stage-1 differential: {compared_records} records compared over {total} cases; \
         {panicking_cases} cases ({:.1}%) reached production's reachable debug_assert! and \
         panicked in both walkers",
        100.0 * panicking_cases as f64 / total as f64,
    );
}

/// **A production defect this differential found, pinned — and the reason `WalkOutcome`
/// carries a panic flag at all.**
///
/// `apply_events_to_ref_into` asserts in debug that every event it is handed is anchored
/// at or after the record's own position
/// ([open_record.rs](../../../pileup/walker/open_record.rs)). Spec §8 already records that
/// `events_overlapping` **does not clip a deletion to the window** — "one anchored before
/// the record can report an anchor below `record_pos`" — but not that production carries a
/// `debug_assert!` requiring the opposite. It does, and it is reachable:
///
/// 1. `pair3/Second` carries a deletion anchored at **17**, spanning 18–22.
/// 2. At position 17 it overlaps its own mate, `pair3/First`. Mate-overlap reconciliation
///    in the indel regime **collapses the pair to a single observation**, and the
///    contributor that carried the indel is the one removed — so no record opens at 17.
/// 3. `pair9/First`'s own deletion, anchored at **19**, then opens a record at 19 — inside
///    the footprint of a deletion that never opened one.
/// 4. Later, where `pair3/Second` matches again (23+), it folds into that record, and
///    `events_overlapping` hands the fold its deletion **anchored at 17**, two positions
///    before `record_pos`.
///
/// In a debug build that panics. **In release, `saturating_sub` clamps the offset to 0**,
/// so the deletion is applied at the record's first base — wrong allele bytes, no error,
/// at exactly the long-deletion loci this port exists to get right.
///
/// Production is frozen, so this is **recorded, not fixed**. What this test asserts is the
/// parity claim: the copy reaches the same precondition on the same input. Its value is
/// that it fails if ng's copy ever stops behaving like production here — including when
/// plan 3 changes the haplotype builder, which is precisely the code that would fix it.
#[test]
#[cfg(debug_assertions)]
fn both_walkers_panic_on_a_deletion_anchored_before_its_record() {
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

    assert!(
        theirs.panicked,
        "production should still reach its own debug_assert! on this input; if it no \
         longer does, production has been fixed and this test is the record of what it \
         used to do"
    );
    assert!(
        ours.panicked,
        "ng's copy is verbatim, so it must reach the same precondition on the same input"
    );
    assert_eq!(
        ours.records, theirs.records,
        "and the records emitted before the panic must agree"
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
    let mut errors = 0usize;

    for seed in SEEDS {
        let mut rng = SplitMix64(seed);
        for _ in 0..CASES_PER_SEED {
            let case = generate(&mut rng);
            adaptor_boundaries += case.reads_with_live_adaptor_boundary;
            let outcome = production_walk(&case);
            if let Some(summary) = outcome.summary {
                widens += summary[2];
                mate_overlaps += summary[3];
                chain_allocations += summary[4];
                cap_truncations += summary[7];
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
         {chain_allocations} chain allocations, {multi_base_records} multi-base records, \
         {multi_allele_records} multi-allele records, {errors} walker errors"
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
    let production_reads: Vec<ProductionPreparedRead> = ng_reads
        .iter()
        .cloned()
        .map(NgPreparedRead::into_production)
        .collect();

    // One fetcher, lent to both walks: identical bytes by construction.
    let fetcher = RefSeqFetcher(WindowedRefSeq::new(fasta.clone(), contigs.clone()));
    let config = WalkerConfig::default();

    let mut theirs = Vec::new();
    let mut their_summary = None;
    let their_panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut walker = production_run(production_reads, &fetcher, &config);
        for item in &mut walker {
            theirs.push(item.map_err(|error| error.to_string()));
        }
        their_summary = Some(summary_array!(walker.summary()));
    }))
    .is_err();

    let mut ours = Vec::new();
    let mut our_summary = None;
    let our_panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut walker = super::run(ng_reads, &fetcher, &config);
        for item in &mut walker {
            ours.push(item.map_err(|error| error.to_string()));
        }
        our_summary = Some(summary_array!(walker.summary()));
    }))
    .is_err();

    let reads = production_reads_len(&theirs);
    let where_ = format!("{reads_path} {region:?}");
    let panicked = assert_same_walk(
        &where_,
        &WalkOutcome {
            records: ours,
            summary: our_summary,
            panicked: our_panic,
        },
        &WalkOutcome {
            records: theirs,
            summary: their_summary,
            panicked: their_panic,
        },
    );
    assert!(
        !panicked,
        "{where_}: both walkers panicked, which agrees but measures nothing — re-run with \
         --release, where production's debug_assert! is compiled out"
    );
    // Joined at the end rather than dropped: a `.fai` that does not describe this FASTA
    // would mean the two walks agreed about the wrong bases, which is a green run that
    // proves nothing.
    if let Some(handle) = verification {
        handle
            .join()
            .expect("the .fai beside the reference describes it");
    }
    eprintln!("real-data differential: {where_} — {reads} records compared, zero divergences");
}

/// The number of stream items, named so the message above reads as what it is.
fn production_reads_len(records: &[Result<PileupRecord, String>]) -> usize {
    records.len()
}
