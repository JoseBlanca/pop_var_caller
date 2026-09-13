//! **What is left of the walker parity oracle once production is gone: ng's walk over a
//! seeded read generator, and production's answers frozen where a comparison survives.**
//!
//! This file began as a differential — one `PreparedRead` stream handed to production's
//! walker and to ng's copy, the two outputs required to agree
//! (`doc/devel/ng/spec/locus_generation_pileup.md` §3, §13.1). Its record, banked before it
//! was cut down, is in this file's history at `d9e7b076`: six mutations of ng's copy each
//! caught inside the first nineteen generated cases, and zero divergences over 437,124
//! records of real reads (GIAB HG002 at 10× and 300×, a tomato CRAM) and 1.5 million
//! synthetic ones.
//!
//! # What was retired, and why (promotion plan step C1, 2026-09-13)
//!
//! - **The two whole-output differentials** — "ng agrees with production where production
//!   fabricated nothing" and "every divergence is one of the six named classes" — had been
//!   `#[ignore]`d since 2026-09-11, when an insertion's record was widened to the ground the
//!   insertion could occupy and the two walkers' streams stopped holding the same records.
//!   The owner's ruling then was that ng's aim is to improve on production, not to match it.
//!   The census, the divergence classes and the projection of production's record into ng's
//!   type existed only for those two, and went with them.
//! - **The real-data differential** asserted the same equality on GIAB and tomato reads and
//!   needed data not in the tree; its figures are in the first paragraph.
//! - **Two tests of the harness's own plumbing** — that the projection maps every field of
//!   production's record, and that both walkers were served the same reference bytes — have
//!   no subject once there is no second walker.
//!
//! # What stays
//!
//! - **The error stream** (`both_walkers_report_the_same_error_on_the_same_malformed_input`)
//!   compares ng's walk against production's answers, written out by hand below from a run of
//!   production's walker at `d9e7b076`. It is the one comparison that still means what it
//!   says: no fixture that gets past admission carries an insertion, so the widened insertion
//!   record that ended the whole-output comparison does not reach it.
//! - **Determinism across processes** and its positive control, which never involved
//!   production.
//! - **The generator's coverage floor and the adaptor-filter control**, which used to read
//!   production's walk and now read ng's: they are claims about what the generator reaches
//!   and ng's walker has every path they count.
//! - **The deletion-anchored-before-its-record regression**, on ng's walk alone. What
//!   production emitted on that fixture is recorded in the test's comments.
//! - **The two ng-only properties** — observations come out in one fixed order, and none is
//!   emitted without a read.

use std::collections::BTreeSet;
use std::sync::Arc;

use super::tests::MockFasta;
use super::{CigarOp, MateRole, PreparedRead, WalkerConfig};
use crate::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation, WitnessedLocusPositions,
};
use crate::ng::read::PLACEHOLDER_READ_GROUP;
use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};

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

/// The two fixture contigs. Two rather than one so the chromosome-boundary path — where the
/// walk flushes the open-record table and the allocator keeps its counter — is on the walk
/// rather than assumed.
const CONTIG_LENGTH: usize = 160;
const CONTIGS: usize = 2;

/// One generated case: a reference, a coordinate-sorted read stream, and the config to
/// walk them under.
struct Case {
    reference: Vec<String>,
    reads: Vec<PreparedRead>,
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
    mate_role: MateRole,
) -> (PreparedRead, bool) {
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

    let read = PreparedRead {
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
        read_group: PLACEHOLDER_READ_GROUP,
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

    let mut reads: Vec<PreparedRead> = Vec::new();
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
                MateRole::FirstOfPair,
            );
            let overlap_start = start + (rng.below(6) as u32);
            let (second, live_b) = generate_read(
                rng,
                &reference,
                chrom_id,
                overlap_start,
                &qname,
                MateRole::SecondOfPair,
            );
            reads_with_live_adaptor_boundary += usize::from(live_a) + usize::from(live_b);
            reads.push(first);
            reads.push(second);
        } else {
            let qname = format!("solo{index}");
            let (read, live) =
                generate_read(rng, &reference, chrom_id, start, &qname, MateRole::Solo);
            reads_with_live_adaptor_boundary += usize::from(live);
            reads.push(read);
        }
    }

    // The walker's coordinate-order invariant is fatal (`WalkerError::OutOfOrder`), and a
    // run that ends at the first read would test nothing beyond it. A **stable** sort, so
    // the admission order of equal-position reads is a property of the generator rather
    // than of the sort — which matters because the depth cap truncates in admission order.
    reads.sort_by_key(|read| (read.chrom_id, read.alignment_start));

    // Most cases run the default limits. A third run tiny column caps, because the cap
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

/// The eight `RunSummary` counters production's walker also kept, **named** — the part of a
/// walk's summary the frozen answers below can be compared against.
///
/// Named rather than an `[u64; 8]` beside a parallel `[&str; 8]`, because two arrays that
/// must stay in step are two arrays that can drift — and the failure mode is a report that
/// names the wrong counter, which is worse than no name at all.
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

/// ng's `RunSummary`, cut down to [`SummaryCounters`].
///
/// **Exhaustively destructured, no `..`**, so a counter added to ng's summary stops this
/// compiling and whoever adds it decides whether it belongs in the comparison. The ones
/// dropped by name below are ng's alone: production had no counterpart, so its frozen answers
/// say nothing about them. `reads_silent_over_footprint` counts reads that were admitted and
/// never contributed anywhere; `reads_with_holed_witness` and `hole_positions` count reads
/// blind in the middle of a record, which production could not represent;
/// `reads_shed_at_admission` counts reads ng refuses at the door once `max_active_reads` are
/// open, where production aborted the walk instead.
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
        reads_shed_at_admission: _,
        reads_silent_over_footprint: _,
        reads_with_holed_witness: _,
        hole_positions: _,
        // **Five more of ng's alone, dropped by name for the same reason as the four
        // above** (the depth-cap change, 2026-08-05). `reads_evicted_at_ceiling` is the
        // other half of `reads_shed_at_admission`: the hold ceiling now gives back a read
        // it is holding rather than only refusing the arrival, and production, which
        // aborts the walk at its ceiling, has no counterpart to either.
        // `positions_short_of_cap` and `short_of_cap_deficit` measure what those two cost
        // the evidence, a question production cannot be asked because it never gets past
        // the ceiling to have an answer. `column_depth_high_water` and
        // `pending_mates_high_water` are peaks ng reports and production does not track.
        reads_evicted_at_ceiling: _,
        positions_short_of_cap: _,
        short_of_cap_deficit: _,
        column_depth_high_water: _,
        pending_mates_high_water: _,
        // The fairness census over capped columns — a measurement of ng's own sampling
        // rule, which is the thing production does not have.
        capped_column_reads_seen: _,
        capped_column_reads_seen_placed_left: _,
        capped_column_reads_kept: _,
        capped_column_reads_kept_placed_left: _,
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
/// **The panic channel is the panic *message*, not a `bool`.** An earlier version stored
/// `catch_unwind(..).is_err()`, and a review demonstrated the hole: replacing ng's copy of a
/// reachable `debug_assert!` with an unrelated `panic!` left the whole suite green, because
/// "it stopped" is not "it reached the same precondition". `open_record.rs` alone carries
/// eight distinct `debug_assert!`s. No fixture that survives in this file panics, and
/// [`assert_same_walk`] says so rather than assuming it.
struct WalkOutcome {
    /// The emitted loci, and each walker error rendered by [`render_ng_error`].
    records: Vec<Result<SampleLocusObservations, String>>,
    /// `None` when the walk panicked — there is no summary to read from a walker that
    /// unwound out from under us.
    summary: Option<SummaryCounters>,
    /// `Some(message)` when the walk panicked.
    panic_message: Option<String>,
}

/// A reference-fetch failure, said in the terms production's walker could also state — which
/// is the form the frozen error below is written in.
///
/// ng's `WalkerError::Fasta` carries a [`RefSeqError`](crate::ng::ref_seq::RefSeqError);
/// production's carried a `ChromRefFetchError`, a variant-for-variant equivalent under a
/// different name, so the two were rendered through this one shape to be comparable. **The
/// contig identifier is deliberately not part of it**: production named contigs and ng
/// numbers them, and `WalkerError::Fasta` carries `chrom_id` outside the source anyway.
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

/// ng's `WalkerError`, rendered for comparison against a frozen answer.
///
/// Only `Fasta` is rewritten, through [`FetchFailure`]. The other eight variants print no
/// module path under `{:?}`, so their rendering is the string production's walker printed
/// for the same variant. They are **listed by name rather than matched with `_`**, so a
/// variant added to the enum stops this compiling instead of silently taking the verbatim
/// path.
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

/// Render a panic payload as its message.
fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "<non-string panic payload>".to_string()
    }
}

/// Drive a walker to exhaustion, surviving a panic with whatever it emitted first.
///
/// `records` is filled inside the closure and read after the unwind: `Vec::push` is
/// exception-safe, so on a panic it holds exactly the prefix that was emitted, which is what
/// makes a comparison element-wise rather than all-or-nothing.
fn drive_ng<W>(mut walker: W, summary_of: impl FnOnce(&W) -> SummaryCounters) -> WalkOutcome
where
    W: Iterator<Item = Result<SampleLocusObservations, super::WalkerError>>,
{
    let mut records = Vec::new();
    let mut summary = None;
    let panic_message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        for item in &mut walker {
            records.push(item.map_err(|error| render_ng_error(&error)));
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

/// ng's walk over one case, every read in the one placeholder read group.
fn ng_walk(case: &Case) -> WalkOutcome {
    ng_walk_in_groups(case, 1)
}

/// ng's walk with the case's reads dealt round-robin into `groups` read groups.
///
/// The read group reaches **nothing but the observation key** — `open_record.rs` reads it only
/// to build `ObservationKey` and to order the emitted observations — so the walk is the walk of
/// `groups == 1`, and the only difference in the output is that one allele carried by reads of
/// two groups becomes two observations. A one-group fixture cannot tell an order that
/// considers the group from one that does not, which is why
/// `the_projection_orders_observations_as_the_walk_does` walks two.
fn ng_walk_in_groups(case: &Case, groups: u32) -> WalkOutcome {
    assert!(groups >= 1, "a walk has at least one read group");
    let fasta = case.fasta();
    let reads: Vec<PreparedRead> = case
        .reads
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, mut read)| {
            if groups > 1 {
                read.read_group = ReadGroupId(index as u32 % groups);
            }
            read
        })
        .collect();
    drive_ng(super::run(reads, fasta, &case.config), |walker| {
        ng_counters(walker.summary())
    })
}

/// ng's emission order for a locus's observations: by bases, then by how much of the locus the
/// read witnessed, then by read group.
///
/// The `ReadWitness` half of the comparator is **the type's own**
/// ([`ReadWitness::sort_key`]) rather than a second spelling here: `finalise` sorts
/// `KeyedObservation`s and this sorts `SequenceObservation`s, so the loop cannot be shared, but
/// the one piece that could silently drift is.
/// `the_projection_orders_observations_as_the_walk_does` covers the rest, by asserting that
/// sorting an ng locus's observations with this function leaves them where the walk emitted
/// them.
fn sort_observations(observations: &mut [SequenceObservation]) {
    observations.sort_by(|a, b| {
        a.bases
            .cmp(&b.bases)
            .then_with(|| a.read_witness.sort_key().cmp(&b.read_witness.sort_key()))
            .then_with(|| a.read_group.0.cmp(&b.read_group.0))
    });
}

/// The surface a walk is compared against a frozen answer on — **three normalisations, applied
/// to both sides**, so neither can hide a difference the other does not have.
///
/// 1. Observations sorted into ng's emission order ([`sort_observations`]). ng's already are;
///    production emitted them in the order it created its buckets.
/// 2. The two per-locus counters, `reads_without_observation` and `reads_discarded_by_cap`,
///    zeroed. Production kept neither per record.
/// 3. The chain ids of a complete observation that matches the reference cleared. Since the
///    owner's ruling of 2026-08-17 ng names every read it folds, where production withheld the
///    ids of its reference bucket.
fn comparable(locus: &SampleLocusObservations) -> ComparableLocus {
    ComparableLocus(comparable_exact_q_sum(locus))
}

/// A locus compared **field by field, exactly** — `q_sum` included, since it is an integer
/// count of steps ([`SummedLogError`]) and so no longer depends on the order it was summed in.
///
/// A newtype because `SampleLocusObservations`' own `PartialEq` is derived, and a comparison
/// against a frozen answer wants to say for itself which fields it compares and how.
struct ComparableLocus(SampleLocusObservations);

impl PartialEq for ComparableLocus {
    fn eq(&self, other: &Self) -> bool {
        let (ours, theirs) = (&self.0, &other.0);
        // Exhaustive destructure on both sides: a field added to the locus type stops this
        // compiling rather than going silently uncompared.
        let SampleLocusObservations {
            region: our_region,
            reference_bases: our_bases,
            observations: our_observations,
            reads_without_observation: our_without,
            reads_discarded_by_cap: our_capped,
            kind: our_kind,
        } = ours;
        let SampleLocusObservations {
            region: their_region,
            reference_bases: their_bases,
            observations: their_observations,
            reads_without_observation: their_without,
            reads_discarded_by_cap: their_capped,
            kind: their_kind,
        } = theirs;
        our_region == their_region
            && our_bases == their_bases
            && our_without == their_without
            && our_capped == their_capped
            && our_kind == their_kind
            && our_observations.len() == their_observations.len()
            && our_observations
                .iter()
                .zip(their_observations)
                .all(|(ours, theirs)| {
                    let SequenceObservation {
                        bases: our_observation_bases,
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
                    our_observation_bases == &theirs.bases
                        && our_witness == &theirs.read_witness
                        && our_group == &theirs.read_group
                        && our_obs == &theirs.num_obs
                        && our_fwd == &theirs.num_fwd
                        && our_mapq == &theirs.mapq_sum
                        && our_mapq_sq == &theirs.mapq_sum_sq
                        && our_placed_left == &theirs.placed_left
                        && our_ids == &theirs.chain_ids
                        && our_q_sum == &theirs.q_sum
                })
    }
}

impl std::fmt::Debug for ComparableLocus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// The normalisation [`comparable`] wraps — see there for the three rules.
fn comparable_exact_q_sum(locus: &SampleLocusObservations) -> SampleLocusObservations {
    let mut out = locus.clone();
    sort_observations(&mut out.observations);
    out.reads_without_observation = 0;
    out.reads_discarded_by_cap = 0;
    // **Complete observations only**, which is `matches_reference`'s own precondition: a
    // partial's bases stop where its read's witness stopped, so comparing them against the
    // whole locus's reference asks a question about bases the read never saw.
    let reference_bases = std::mem::take(&mut out.reference_bases);
    for observation in &mut out.observations {
        if observation.read_witness == ReadWitness::Complete
            && observation.matches_reference(&reference_bases)
        {
            observation.chain_ids.clear();
        }
    }
    out.reference_bases = reference_bases;
    out
}

/// Assert ng's walk matches a frozen answer: the same way of stopping, the same stream item for
/// item, and the same eight counters.
#[track_caller]
fn assert_same_walk(where_: &str, ours: &WalkOutcome, frozen: &WalkOutcome) {
    // **The panic comparison runs first**, and it compares the *message*. Ordering matters
    // for the diagnosis, not the verdict — with the length check first, a walk that panics
    // early reports "ng emitted 3 stream items, the frozen answer 24", which names the symptom
    // and hides the cause.
    assert_eq!(
        ours.panic_message, frozen.panic_message,
        "{where_}: the walk did not stop the way the frozen answer did",
    );

    assert_eq!(
        ours.records.len(),
        frozen.records.len(),
        "{where_}: ng emitted {} stream items, the frozen answer {}",
        ours.records.len(),
        frozen.records.len(),
    );
    for (position, (ours, frozen)) in ours.records.iter().zip(frozen.records.iter()).enumerate() {
        let normalise = |item: &Result<SampleLocusObservations, String>| {
            item.as_ref().map(comparable).map_err(String::clone)
        };
        assert_eq!(
            normalise(ours),
            normalise(frozen),
            "{where_}: stream item {position} diverged"
        );
    }

    assert_eq!(
        ours.summary, frozen.summary,
        "{where_}: the RunSummary counters diverged",
    );
}

/// Cases per seed. Small enough to stay a unit test, large enough that every behaviour
/// `the_generator_exercises_what_the_port_can_break` counts fires many times over.
const CASES_PER_SEED: usize = 400;

/// Cases per seed, overridable by `PVC_PARITY_CASES` so a soak run is one command away —
/// the convention `delimit_parity` set.
///
/// **Use `--profile soak`, not `--release`** (Cargo.toml): `[profile.release]` leaves
/// `debug-assertions` off, and most of what these walks check is asserted along the way.
/// `soak` is release-speed with the assertions and overflow checks armed:
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

/// Whether any allele at this locus is carried by more than one read group.
///
/// Two observations sharing `(bases, read_witness)` can only differ in the group, the three
/// together being the whole observation identity.
fn observations_split_by_group(locus: &SampleLocusObservations) -> bool {
    /// An observation's identity **without** its read group: the bases, and the witness as the
    /// type's own total order gives it ([`ReadWitness::sort_key`]).
    ///
    /// The witness half is a tag and the borrowed set since C2, where it was a fixed
    /// `(u8, u16, u16)` — a set has no fixed width. It still borrows from the observation, so
    /// the key lives exactly as long as the loop below.
    type ObservationIdentityWithoutGroup<'a> =
        (&'a [u8], (u8, Option<&'a WitnessedLocusPositions>));

    let mut seen: BTreeSet<ObservationIdentityWithoutGroup<'_>> = BTreeSet::new();
    for observation in &locus.observations {
        if !seen.insert((&observation.bases, observation.read_witness.sort_key())) {
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
/// Over every field of ng's locus as `{:?}` prints it, so an observation split or merged
/// differently — the axis a hash-order bug would show up on — changes the digest.
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
        let mut reads = case.reads.clone();
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
/// any number of runs in the same binary, and no other test in this file compares one ng
/// run against another on the same input. Two
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

/// **ng's walk emits a locus's observations in [`sort_observations`] order** — the order
/// `finalise` promises, checked rather than commented.
///
/// `finalise` sorts `KeyedObservation`s and [`sort_observations`] sorts `SequenceObservation`s, so
/// the two loops cannot be shared even though the `ReadWitness` comparator is
/// ([`ReadWitness::sort_key`]). What could still drift is the *rest* of the key — the bases, then
/// the group — so this walks a fixture and asserts that sorting ng's own emitted observations with
/// this function does not move them. The comparison against a frozen answer normalises through
/// the same function, so a drift here would also make that comparison order-blind.
///
/// *(The name is from when the function also laid production's records out; kept because other
/// modules' comments cite it.)*
///
/// The fixture is the general one at two read groups, because a one-group walk cannot
/// distinguish an order that considers the group from one that does not.
#[test]
fn the_projection_orders_observations_as_the_walk_does() {
    let mut loci = 0usize;
    let mut multi_observation = 0usize;
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
                    multi_observation += 1;
                }
                if observations_split_by_group(locus) {
                    grouped += 1;
                }
                let mut sorted = locus.observations.clone();
                sort_observations(&mut sorted);
                assert_eq!(
                    sorted, locus.observations,
                    "{where_}: locus {position} came out of the walk in an order this \
                     function would not have produced, so the two orders have drifted",
                );
            }
        }
    }

    assert!(
        multi_observation * 10 > loci && grouped > 0,
        "only {multi_observation} of {loci} loci carry more than one observation and {grouped} split by read \
         group — an order test needs observations to order"
    );
}

/// **Every emitted observation has a read folded into it** — which is *not* A3's eviction,
/// and D1 is where the difference was found.
///
/// This test used to claim the eviction. It ran on production's record type and asserted that
/// no emitted bucket had `num_obs == 0`, and while the walk emitted that type it was a real
/// check. **B2 made it vacuous:** ng's observations are derived from `folded_reads`, per read, so
/// a bucket no read is folded into produces no observation at all and leaves no trace in the output.
/// D1 mutated the code the test named — moving `evict_unsupported_alleles` above the
/// contributor fold loop, which strands every bucket that loop empties — and all 198 tests
/// in this module stayed green. The eviction is now enforced by a `debug_assert!` in
/// `finalise`, where the buckets still exist and where every walk in the suite checks it.
///
/// What is left here is worth keeping and is a different claim: **`num_obs` is never zero on
/// an emitted observation**, i.e. the re-derivation cannot mint an observation for nothing. That would fail if
/// `keyed_observations` ever created observations from the bucket table instead of from the reads —
/// which is precisely the B1 mistake it was written to avoid.
///
/// **Two more per-observation checks joined it at promotion step C1**, moved from the census the
/// retired differential ran them in: at least half as many chain ids as reads, and every partial
/// witness's runs inside its own locus. Both are claims about ng on its own terms.
#[test]
fn every_emitted_observation_carries_a_read() {
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
                for (index, observation) in locus.observations.iter().enumerate() {
                    assert!(
                        observation.num_obs > 0,
                        "{where_}: emitted locus {position} carries observation {index} \
                         ({:?}) with no supporting read: {locus:?}",
                        String::from_utf8_lossy(&observation.bases),
                    );
                    // Every read folded is named since 2026-08-17, and a read pair whose
                    // mates both reach this observation shares one id — no more than two
                    // reads can — so twice the ids must cover the reads.
                    assert!(
                        2 * observation.chain_ids.len() >= observation.num_obs as usize,
                        "{where_}: locus {position} observation {index} carries {} reads and \
                         only {} chain ids: {locus:?}",
                        observation.num_obs,
                        observation.chain_ids.len(),
                    );
                    // A partial witness's runs are non-empty and lie inside their own
                    // locus. Nothing on the type ties a run to the locus it ends up
                    // attached to, so this is checked on every emitted locus.
                    if let ReadWitness::Partial { positions } = &observation.read_witness {
                        for (start, end) in positions.runs() {
                            assert!(
                                start < end && u64::from(end) <= locus.region.len(),
                                "{where_}: a partial run {start}..{end} is empty or reaches \
                                 past its locus of {} positions",
                                locus.region.len(),
                            );
                        }
                    }
                }
            }
        }
    }

    // Loci where the reads disagreed, so a read has somewhere to move *to* — without them
    // the fixture never exercises a re-fold at all.
    assert!(
        multi_allele * 100 > records,
        "only {multi_allele} of {records} emitted loci carry more than two observations — this \
         fixture cannot exercise a read moving between buckets, so the property it asserts \
         is vacuous"
    );
    eprintln!(
        "{records} emitted loci, {multi_allele} of them with more than two observations; every observation \
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
/// **This test was the parity claim and the regression together.** Since promotion step C1 it
/// is the regression alone, on ng's walk: production's walker is gone. What production emitted
/// on this fixture at `d9e7b076`, after the fix, is written beside each assertion below.
#[test]
fn a_deletion_anchored_before_its_record_contributes_none_of_the_bases_it_deleted() {
    fn read(
        qname: &str,
        role: MateRole,
        start: u32,
        end: u32,
        cigar: Vec<CigarOp>,
        seq_len: usize,
    ) -> PreparedRead {
        PreparedRead {
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
            read_group: PLACEHOLDER_READ_GROUP,
        }
    }

    let case = Case {
        reference: vec!["ACGT".repeat(40), "ACGT".repeat(40)],
        reads: vec![
            // The mate that reconciliation keeps, suppressing the record at 17.
            read(
                "pair3",
                MateRole::FirstOfPair,
                4,
                22,
                vec![CigarOp::Match(6), CigarOp::Match(8), CigarOp::Match(5)],
                19,
            ),
            // The deletion anchored at 17, spanning 18–22, that opens no record.
            read(
                "pair3",
                MateRole::SecondOfPair,
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
                MateRole::FirstOfPair,
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

    let ours = ng_walk(&case);

    // Production's walker, fixed, walked this cleanly too.
    assert_eq!(
        ours.panic_message, None,
        "the input that used to trip the debug_assert! must now walk cleanly"
    );

    // The record at 19 is the one the defect corrupted: `pair9/First`'s deletion opens it
    // spanning 19..=25, and `pair3/Second` folds in carrying a deletion anchored at 17
    // whose run covers 18–22 — so of this record's seven positions it witnessed only
    // 23, 24, 25.
    let (reference, alleles) = ours
        .records
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .find(|locus| locus.region.start.get() == 19)
        .map(|locus| {
            let reference = String::from_utf8_lossy(&locus.reference_bases).to_string();
            let alleles: Vec<String> = locus
                .observations
                .iter()
                .map(|observation| String::from_utf8_lossy(&observation.bases).to_string())
                .collect();
            (reference, alleles)
        })
        .expect("pair9's deletion opens a record at 19");

    // Production's record here had the same seven reference bases.
    assert_eq!(
        reference,
        "ACGT".repeat(40)[18..25],
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
    // satisfy either alone. Production, fixed, held `AAA` and not the pre-fix spelling.
    let before_the_fix = format!("{}AAA", &reference[..1]);
    assert!(
        alleles.iter().any(|allele| allele == "AAA"),
        "the read whose deletion covers 19–22 should contribute the three bases it witnessed \
         (AAA), but the record holds {alleles:?}",
    );
    assert!(
        !alleles.contains(&before_the_fix),
        "the record still holds {before_the_fix}, the pre-fix spelling — a leading base at 19 \
         that this read had deleted. Records: {alleles:?}",
    );

    // **And this is the witnessed-extent rule, on the fixture the defect handed it.**
    // `pair3/First` matched 4–22, so of the record's seven positions it witnessed four:
    // 19, 20, 21, 22. Production folded it as `AAAAGTA` — `AAAA` plus the reference bases at
    // 23, 24 and 25, three bases it never sequenced. ng emits the four it saw.
    assert!(
        alleles.iter().any(|allele| allele == "AAAA"),
        "ng should carry only the four positions the read witnessed: {alleles:?}",
    );
    assert!(
        !alleles.iter().any(|allele| allele == "AAAAGTA"),
        "ng carries production's reference-padded spelling: {alleles:?}",
    );
}

/// **The generated cases are worth exactly what they reach, and this is what says what they
/// reach.**
///
/// A test that passes on inputs exercising none of the interesting paths is a claim dressed
/// as evidence. So the behaviours the walker's port was mutation-tested on — record widening,
/// mate-overlap reconciliation, the column depth cap, adaptor masking, chain allocation, and
/// the mate-lookup window evicting — are counted here over the generator the other tests in
/// this file use, and each is required to fire.
///
/// The counts come from ng's own `RunSummary` wherever it has one, so they measure what the
/// *walker* did rather than what the generator intended. Adaptor boundaries have no counter,
/// so they are tallied at generation. *(Until promotion step C1 they were read off
/// production's walk of the same cases.)*
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
        // Scales with `PVC_PARITY_CASES` like the other walks do: a soak that walked more
        // cases but kept the coverage floor fixed would report more confidence in the same
        // evidence.
        for _ in 0..cases_per_seed() {
            let case = generate(&mut rng);
            adaptor_boundaries += case.reads_with_live_adaptor_boundary;
            let outcome = ng_walk(&case);
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
                        // An observation is a supported allele, so "more than one" counts
                        // loci where the reads disagreed.
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

    // Each of these is a behaviour the port was mutation-tested on. A zero here means the
    // cases never reach it.
    assert!(widens > 0, "no record ever widened");
    assert!(mate_overlaps > 0, "mate-overlap reconciliation never fired");
    assert!(cap_truncations > 0, "the column depth cap never fired");
    assert!(
        adaptor_boundaries > 0,
        "no read carried a live adaptor boundary"
    );
    assert!(chain_allocations > 0, "no chain id was ever allocated");
    // The path where a pair silently degrades to two solos — a different chain id and a
    // different fold.
    assert!(
        mate_evictions > 0,
        "the mate-lookup window never evicted, so the cases never reach a pair falling back \
         to two unpaired reads"
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

/// **The `Err` half of the stream, which the main generator deliberately never produces** —
/// compared against production's walker's answers, frozen.
///
/// The main generator keeps reads in bounds and in order, so its walks never emit an error.
/// Every `WalkerError` variant is fatal and terminal for the iterator, so these cases cannot
/// live in the main generator without truncating the walks that test the fold. They get their
/// own fixtures.
///
/// # The frozen answers
///
/// Each fixture's expected outcome is what production's walker (`run` in `src/pileup/walker/`)
/// emitted on the same reads at commit `d9e7b076`, with its errors rendered into the shape
/// [`render_ng_error`] produces and its records projected into ng's locus type, then written
/// out here by hand when promotion step C1 removed the call. On all five fixtures production's
/// walk and ng's agreed at that commit. The expected error string is complete, so it is also what shows each fixture
/// **reaches** the error it is named for rather than passing on some other failure.
#[test]
fn both_walkers_report_the_same_error_on_the_same_malformed_input() {
    fn read(
        qname: &str,
        start: u32,
        end: u32,
        cigar: Vec<CigarOp>,
        seq_len: usize,
    ) -> PreparedRead {
        PreparedRead {
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
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
            read_group: PLACEHOLDER_READ_GROUP,
        }
    }

    /// Production's eight counters, in [`SummaryCounters`]' field order.
    fn counters(
        [
            reads_admitted,
            records_emitted,
            record_widen_events,
            mate_overlap_positions,
            chain_allocations,
            active_reads_high_water,
            mate_lookup_evictions,
            column_depth_truncations,
        ]: [u64; 8],
    ) -> Option<SummaryCounters> {
        Some(SummaryCounters {
            reads_admitted,
            records_emitted,
            record_widen_events,
            mate_overlap_positions,
            chain_allocations,
            active_reads_high_water,
            mate_lookup_evictions,
            column_depth_truncations,
        })
    }

    /// One of the single-base loci production emitted ahead of the fetch failure: the read's
    /// `A` over one reference base, forward strand, MAPQ 60. The error is the mapping quality's
    /// −3 nats (12,288 steps), which outweighs base quality 30's −6.9.
    /// `chain_ids` is as [`comparable`] leaves it — empty where the `A` matches the reference.
    fn frozen_locus(position: u64, reference: u8, placed_left: u32) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(position),
                end: Position(position),
            },
            reference_bases: Box::from([reference]),
            observations: vec![SequenceObservation {
                bases: Box::from(&b"A"[..]),
                read_witness: ReadWitness::Complete,
                read_group: PLACEHOLDER_READ_GROUP,
                num_obs: 1,
                num_fwd: 1,
                q_sum: SummedLogError::from_steps(-12_288),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left,
                chain_ids: if reference == b'A' {
                    Vec::new()
                } else {
                    vec![0]
                },
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    let error = |rendered: &str| WalkOutcome {
        records: vec![Err(rendered.to_string())],
        summary: counters([0; 8]),
        panic_message: None,
    };

    let reference = vec!["ACGT".repeat(40), "ACGT".repeat(40)];
    let fixtures: Vec<(&str, Vec<PreparedRead>, WalkOutcome)> = vec![
        (
            "out of order",
            vec![
                read("a", 20, 27, vec![CigarOp::Match(8)], 8),
                read("b", 4, 11, vec![CigarOp::Match(8)], 8),
            ],
            WalkOutcome {
                summary: counters([1, 0, 0, 0, 1, 1, 0, 0]),
                ..error(
                    "OutOfOrder { qname: \"b\", prev_chrom_id: 0, prev_pos: 20, chrom_id: 0, \
                     pos: 4 }",
                )
            },
        ),
        (
            // The check is `alignment_end < alignment_start`, not "the CIGAR consumes no
            // reference" — an all-insertion read whose `alignment_end` equals its start
            // sails through.
            "zero reference span",
            vec![read("i", 4, 3, vec![CigarOp::Insertion(4)], 4)],
            error("ZeroRefSpan { qname: \"i\", chrom_id: 0, pos: 4 }"),
        ),
        (
            "cigar consumes more read bases than seq provides",
            vec![read("m", 4, 11, vec![CigarOp::Match(8)], 5)],
            error(
                "MalformedRead { reason: \"CIGAR consumes 8 read bases but seq.len = 5\", \
                 qname: \"m\", chrom_id: 0, pos: 4 }",
            ),
        ),
        (
            "seq and bq of different lengths",
            vec![{
                let mut malformed = read("q", 4, 11, vec![CigarOp::Match(8)], 8);
                malformed.bq_baq.truncate(7);
                malformed
            }],
            error(
                "MalformedRead { reason: \"seq.len (8) != bq_baq.len (7)\", qname: \"q\", \
                 chrom_id: 0, pos: 4 }",
            ),
        ),
        (
            // The reference is 160 bases and this read's footprint runs past its end, so
            // `open_new` asks for bases that do not exist. `WalkerError::Fasta` is the one
            // variant whose source ng states in a different type from production, so this is
            // the fixture that exercises [`FetchFailure`]. The main generator clamps every read
            // inside its contig precisely to avoid this path.
            //
            // The failing fetch is at **161**, not at the read's start: a record is opened per
            // covered position, so the walk emits loci at 156–159 first and fails at the first
            // position past the contig.
            "fetch past the contig end",
            vec![read("far", 156, 175, vec![CigarOp::Match(20)], 20)],
            WalkOutcome {
                records: vec![
                    Ok(frozen_locus(156, b'T', 0)),
                    Ok(frozen_locus(157, b'A', 1)),
                    Ok(frozen_locus(158, b'C', 1)),
                    Ok(frozen_locus(159, b'G', 1)),
                    Err(
                        "Fasta { chrom_id: 0, start: 161, start_plus_len: 162, cause: \
                         OutOfBounds { contig_length: 160, start: 161, end: 162 } }"
                            .to_string(),
                    ),
                ],
                summary: counters([1, 4, 0, 0, 1, 1, 0, 0]),
                panic_message: None,
            },
        ),
    ];

    for (name, reads, frozen) in fixtures {
        let case = Case {
            reference: reference.clone(),
            reads,
            config: WalkerConfig::default(),
            reads_with_live_adaptor_boundary: 0,
        };
        assert_same_walk(name, &ng_walk(&case), &frozen);
    }
}

/// **The adaptor counter's real claim, checked end to end.**
///
/// `reads_with_live_adaptor_boundary` is tallied at generation and cannot see the walk, so
/// however carefully it computes liveness it remains a statement about the *input*. The
/// only honest test of "the G1 filter is on this walk" is that **removing it changes the
/// answer** — which is what the port's mutation test of adaptor masking exploited, made
/// permanent here so the property survives a generator change rather than resting on a
/// hand-run exercise. *(Read off production's walk until promotion step C1; ng's walker applies
/// the same filter.)*
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
            if ng_walk(&case).records != ng_walk(&without).records {
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
