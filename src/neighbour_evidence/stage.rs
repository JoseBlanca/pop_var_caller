//! **The pass as a stage in front of the writer** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §2, §5).
//!
//! The stage takes finished lines one at a time, in the order the writer would write them, and
//! gives them back rewritten, in the same order, once each one's window has been solved. It holds
//! the current contig's lines from the left margin of the window being filled to the last line
//! read; when the lines reaching past a window's right margin have arrived it solves the window,
//! hands back its core and forgets what no later window will read. A new contig, or the end of
//! the input, finishes the contig's last windows.
//!
//! It runs on the thread that writes, so the threads that call loci wait while a window is solved
//! and the stage never holds more than a window: the bounded buffer of §2.1, at its tightest.
//! Whether letting calling run ahead would pay is a measurement (§8.6).
//!
//! **Whether it runs at all is decided once, before the header is written** ([`decide`]): on by
//! default for a diploid cohort of at least ten samples whose median inbreeding is at least 0.5,
//! and the decision, with its reason, goes into the header either way.

use rayon::prelude::*;

use super::fit::SampleFit;
use super::line::{RewriteTally, rewrite_line, site_of};
use super::model::GenotypePairs;
use super::peak::{MemoryGuard, PeakCheck, PeakEstimate, PeakMoment, PeakTooClose};
use super::phase::{
    DISAGREEMENT_COST, EXTRA_PATHS_PER_ROUND, JoinRecord, LEAN_MIN, LastBlock, PassReadCounts,
    PathAgreement, SHARED_RECORDS_TO_JOIN, SWITCH_COST, WrittenPhase,
};
use super::records::{RecordNumbers, SiteRecord};
use super::refine::{
    FillCheck, HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW, MAX_HIDDEN_SHARE, MISMATCH, NoReadRule,
    PHASING_ROUNDS, Phasing, RefinedRecord, held_solve_bytes, record_window_bytes,
    refine_window_with_reads, round_store_bytes, solution_chunk_bytes, window_key,
    written_result_bytes,
};
use super::templates::MAX_TEMPLATES;
use super::windows::{
    MarginThinning, NEAR_MARGIN_BASES, StoreAndResults, WindowBytes, WindowSizes, plan_window,
};
use crate::read_phase::fragments::{FragmentAllele, PackedRecordFragments, RecordFragments};
use crate::types::{ContigId, Ploidy};
use crate::vcf::RecordPlace;

/// The fewest samples the pass runs on by default (`N_min`, §5) — a starting value for §8.
pub const MIN_SAMPLES: usize = 10;

/// The lowest median inbreeding the pass runs on by default (`F_min`, §5) — a starting value.
pub const MIN_MEDIAN_INBREEDING: f64 = 0.5;

/// **`--refine-genotypes`**: whether the run decides, or the operator did (§5, ruled §9 Q3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NeighbourEvidenceSwitch {
    /// Run where the cohort is large and inbred enough.
    Auto,
    /// Run wherever the pass can: two or more samples, diploid.
    On,
    /// Never run.
    Off,
}

/// **Whether the pass runs in this run, and why** — what the header states (§4.1).
#[derive(Clone, Debug, PartialEq)]
pub struct NeighbourEvidenceDecision {
    /// Whether it runs.
    pub ran: bool,
    /// Why, in words: the switch, the sample count, the median inbreeding, the ploidy.
    pub why: String,
    /// Each sample's inbreeding coefficient from the parameters, in the run's sample order — each
    /// walk's starting `F` (§3.5).
    pub inbreeding_by_sample: Vec<f64>,
    /// The window sizes the run plans with.
    pub sizes: WindowSizes,
    /// The share of drawn paths that must agree at a phase link (`--phase-min-path-agreement`,
    /// phase output spec §9); [`decide`] sets [`PathAgreement::DEFAULT`].
    pub path_agreement: PathAgreement,
    /// Whether each sample's fragments orient its haplotypes in the rounds' draws
    /// (`neighbour_evidence_reads_in_draw.md`); [`decide`] sets it. Off only where a test needs
    /// the link test's reads alone.
    pub(crate) reads_in_draw: bool,
    /// Whether the core floor was set on its own (`--refine-core-floor`). Where it was not, it is
    /// the margin floor and follows it ([`Self::with_window_floor`]), on a fresh run and on a
    /// resumed run given a new `--refine-window-floor` alike.
    pub(crate) core_floor_given: bool,
    /// **What the memory check is given** (`--refine-ignore-memory-limit`, pass spec §6, memory):
    /// the memory this run may use, and whether to go on past it. [`decide`] checks nothing;
    /// a command sets it before calling starts. It changes no output.
    pub memory_guard: MemoryGuard,
}

impl NeighbourEvidenceDecision {
    /// The same decision, a phase link kept where `agreement` of the drawn paths agree.
    #[must_use]
    pub fn with_path_agreement(self, agreement: PathAgreement) -> Self {
        Self {
            path_agreement: agreement,
            ..self
        }
    }

    /// The same decision, its windows planned within `bytes` (`--refine-window-bytes`).
    #[must_use]
    pub fn with_window_bytes(self, bytes: u64) -> Self {
        Self {
            sizes: WindowSizes {
                budget_bytes: bytes,
                ..self.sizes
            },
            ..self
        }
    }

    /// The same decision, no window's margins (nor, then, its core) shorter than `bases`
    /// (`--refine-window-floor`). A core floor not set on its own ([`Self::with_core_floor`])
    /// becomes `bases` too.
    #[must_use]
    pub fn with_window_floor(self, bases: u64) -> Self {
        Self {
            sizes: WindowSizes {
                margin_floor_bases: bases,
                core_floor_bases: if self.core_floor_given {
                    self.sizes.core_floor_bases
                } else {
                    bases
                },
                ..self.sizes
            },
            ..self
        }
    }

    /// The same decision, no window's core shorter than `bases` once its margins are at their
    /// floor (`--refine-core-floor`); at or above the margin floor, the core stops there too.
    #[must_use]
    pub fn with_core_floor(self, bases: u64) -> Self {
        Self {
            sizes: WindowSizes {
                core_floor_bases: bases,
                ..self.sizes
            },
            core_floor_given: true,
            ..self
        }
    }

    /// The same decision, no window's core longer than `bases` (the hidden `--refine-core-max`),
    /// nor its core floor.
    #[must_use]
    pub fn with_core_max(self, bases: u64) -> Self {
        Self {
            sizes: WindowSizes {
                core_bases: bases,
                core_floor_bases: self.sizes.core_floor_bases.min(bases),
                ..self.sizes
            },
            ..self
        }
    }

    /// The same decision, no window's core nor margins longer than `bases` (the hidden
    /// `--refine-window-max`), nor its floor.
    #[must_use]
    pub fn with_window_max(self, bases: u64) -> Self {
        Self {
            sizes: WindowSizes {
                core_bases: bases,
                margin_bases: bases,
                margin_floor_bases: self.sizes.margin_floor_bases.min(bases),
                ..self.sizes
            },
            ..self
        }
    }

    /// The same decision, its far margins walked at most `thinning` records a kb
    /// (`--refine-margin-density`, pass spec §6, thinned margins); `None` walks every record.
    #[must_use]
    pub fn with_margin_thinning(self, thinning: Option<MarginThinning>) -> Self {
        Self {
            sizes: WindowSizes {
                thinning,
                ..self.sizes
            },
            ..self
        }
    }

    /// The value of the `##neighbourEvidence=` header line.
    #[must_use]
    pub fn as_header_value(&self) -> String {
        let why = self.why.replace('"', "'");
        if self.ran {
            format!(
                "<ran=yes,why=\"{why}\",templates={MAX_TEMPLATES},mismatch={MISMATCH},core={},margin={},margin_floor={},window_bytes={},\
                 rounds={PHASING_ROUNDS},haplotype_templates={MAX_TEMPLATES},\
                 hidden_per_window={HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW},\
                 hidden_share={MAX_HIDDEN_SHARE},extra_paths={EXTRA_PATHS_PER_ROUND},\
                 lean_min={LEAN_MIN},path_agreement={},join_records={SHARED_RECORDS_TO_JOIN},\
                 alignment_costs={DISAGREEMENT_COST}/{SWITCH_COST},\
                 read_link_min={:.2},read_break_min={:.2},read_error_floor={},read_gq_min={},\
                 read_reach={},read_cap={}{}{}>",
                self.sizes.core_bases,
                self.sizes.margin_bases,
                self.sizes.margin_floor_bases,
                self.sizes.budget_bytes,
                self.path_agreement.get(),
                crate::read_phase::link::LINK_MIN,
                crate::read_phase::link::BREAK_MIN,
                crate::read_phase::link::ERROR_FLOOR,
                crate::read_phase::link::GQ_MIN,
                crate::read_phase::fragments::READ_REACH,
                crate::read_phase::capture::READ_ALLELES_KEPT,
                // Only where the margins are thinned, so an unthinned run's header is as it was.
                self.sizes
                    .thinning
                    .map_or(String::new(), |thinning| format!(
                        ",margin_density={thinning},margin_near={NEAR_MARGIN_BASES}"
                    )),
                // Only where a core may be shorter than the margin floor, likewise.
                if self.sizes.core_floor_bases
                    < self.sizes.margin_floor_bases.min(self.sizes.core_bases)
                {
                    format!(",core_floor={}", self.sizes.core_floor_bases)
                } else {
                    String::new()
                },
            )
        } else {
            format!("<ran=no,why=\"{why}\">")
        }
    }
}

/// **Decides whether the pass runs** (§5): `switch` from the command line, each sample's
/// inbreeding coefficient in the run's sample order, the run's ploidy.
#[must_use]
pub fn decide(
    switch: NeighbourEvidenceSwitch,
    inbreeding_by_sample: &[f64],
    ploidy: Ploidy,
) -> NeighbourEvidenceDecision {
    let samples = inbreeding_by_sample.len();
    let mut sorted = inbreeding_by_sample.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = match samples {
        0 => 0.0,
        _ if samples % 2 == 1 => sorted[samples / 2],
        _ => f64::midpoint(sorted[samples / 2 - 1], sorted[samples / 2]),
    };
    let cohort = format!("{samples} samples, median inbreeding {median:.2}");
    let (ran, why) = if switch == NeighbourEvidenceSwitch::Off {
        (false, format!("switched off; {cohort}"))
    } else if ploidy.get() != 2 {
        (
            false,
            format!("ploidy {}; the pass models two copies", ploidy.get()),
        )
    } else if samples < 2 {
        (
            false,
            "one sample: no other sample to copy from".to_string(),
        )
    } else if switch == NeighbourEvidenceSwitch::On {
        (true, format!("switched on; {cohort}"))
    } else if samples < MIN_SAMPLES {
        (false, format!("{cohort}; fewer than {MIN_SAMPLES} samples"))
    } else if median < MIN_MEDIAN_INBREEDING {
        (
            false,
            format!("{cohort}; median inbreeding below {MIN_MEDIAN_INBREEDING}"),
        )
    } else {
        (true, cohort)
    };
    NeighbourEvidenceDecision {
        ran,
        why,
        inbreeding_by_sample: inbreeding_by_sample.to_vec(),
        sizes: WindowSizes::default(),
        path_agreement: PathAgreement::DEFAULT,
        reads_in_draw: true,
        core_floor_given: false,
        memory_guard: MemoryGuard::default(),
    }
}

/// **A window the pass cannot solve within its budget** (pass spec §6): at the floor it
/// still holds more bytes than `--refine-window-bytes` allows. The run stops here rather than
/// switch refinement off or go on past the budget (owner, 2026-11-13).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowOverBudget {
    /// The window's contig, counted from 0 in the header's order, and its core's first position.
    pub contig: u32,
    pub position: u64,
    /// The records the window spans, and the samples each holds.
    pub records: usize,
    pub samples: usize,
    /// The bytes it holds, and the budget.
    pub bytes: u64,
    pub budget: u64,
}

/// **Why the stage stopped handing lines on**: a window over its budget, a window expected to need
/// more memory than the run may use, or the writer refusing a line it was handed
/// ([`NeighbourEvidenceStage::push_writing`]).
#[derive(Debug)]
pub enum StageStop<E> {
    /// A window holds more than the budget even at the smallest size.
    OverBudget(WindowOverBudget),
    /// A window's expected peak is above the share of the run's memory limit the run stops at
    /// ([`super::peak::STOP_SHARE`]), and `--refine-ignore-memory-limit` was not given.
    OverMemory(PeakTooClose),
    /// Writing a line failed.
    Write(E),
    /// **The stage had already stopped on a refused write** and was handed more. A stop drops the
    /// lines it caught part-way through a batch, so a stage that went on would write a stream with
    /// lines missing from its middle; it refuses instead. A stage stopped by its budget or the
    /// memory limit gives the same stop again.
    Stopped,
}

/// **What a stage stopped on**, kept so that it refuses every line handed in after.
#[derive(Clone, Copy, Debug)]
enum StopCause {
    OverBudget(WindowOverBudget),
    OverMemory(PeakTooClose),
    Write,
}

impl<E> From<WindowOverBudget> for StageStop<E> {
    fn from(over: WindowOverBudget) -> Self {
        Self::OverBudget(over)
    }
}

/// **What a stage that hands its lines to a vector stops on**: its budget or the memory limit,
/// never a write.
pub type VectorStageStop = StageStop<std::convert::Infallible>;

/// **How many of a core's lines are rewritten together**, on the pool, and handed on before the
/// next ones are: a core's rewritten lines were all held at once before the first was written —
/// on 2,176 tomato samples about 9 GB at the smallest window — and a chunk is a few hundred MB.
const LINES_REWRITTEN_TOGETHER: usize = 1_024;

/// **How many held lines are compressed together**, on the pool, once that many have arrived
/// since the last batch: compressing each as it arrives would put the work on the one thread every
/// line passes through.
const LINES_COMPRESSED_TOGETHER: usize = 1_024;

/// **How many arriving lines wait to be read into sites together**, on the pool, or fewer once
/// their text and read alleles reach [`BYTES_READ_TOGETHER`]. Read one at a time on the one thread
/// every line passes through, reading was the stage's largest serial step outside the solve: on 59
/// tomato samples at 2 threads, 2.7 s of the 10 s the writing pass spent outside it, with the pool
/// idle.
const LINES_READ_TOGETHER: usize = 1_024;

/// **The most bytes that wait to be read**: the lines' text and the read alleles that arrive with
/// them, 12 bytes each before they are packed. At thousands of samples, where a line is about
/// 90 kB, the queue holds about 90 lines rather than 1,024. The read alleles count too because at
/// 30 reads a position they take 8 to 10 times a line's text (on drawn cohorts: about 30 a sample
/// a record, 364 bytes, against 37 to 45 of text). Counted by its text alone, the queue held about
/// ten times this limit at 30 reads.
///
/// **What the limit does not bound**: once a batch is read, its sites — about half to two-thirds
/// of its text (24 to 26 bytes a sample a record against 37 to 51 of text on the same cohorts) —
/// and their read alleles, packed, are held until each line is taken. Lines are taken in order,
/// and a line that completes a window is taken only after that window is solved, so the rest of
/// the batch, up to 1,023 sites, is held through the solve. The budget does not count them, as it
/// does not count the walks' arrays.
const BYTES_READ_TOGETHER: usize = 8 << 20;

/// **The zstd level held lines are compressed at.** On 800 tomato samples level 1 shrank a line
/// about 5.4 times at about 440 MB/s a thread; levels −3 and 3 gave about the same.
const LINE_COMPRESSION_LEVEL: i32 = 1;

/// **A held line's text**: as it arrived, or compressed. A window holds every line it spans until
/// its core is written, and at 45 bytes a sample a line, uncompressed text was a tenth of a window's
/// memory (pass spec §6).
enum HeldText {
    Plain(Vec<u8>),
    Compressed { bytes: Vec<u8>, length: usize },
}

impl HeldText {
    /// The bytes the text holds: its length as it arrived, or compressed.
    fn held_bytes(&self) -> usize {
        match self {
            Self::Plain(line) => line.len(),
            Self::Compressed { bytes, .. } => bytes.len(),
        }
    }
}

/// **One thread's compressor and the buffer it compresses into**, made once a batch a thread
/// rather than once a line.
struct LineCompressor {
    compressor: zstd::bulk::Compressor<'static>,
    scratch: Vec<u8>,
}

impl LineCompressor {
    fn new() -> Self {
        Self {
            // PANIC-FREE: a context at a valid level fails only where allocation does.
            compressor: zstd::bulk::Compressor::new(LINE_COMPRESSION_LEVEL)
                .expect("a zstd context at level 1"),
            scratch: Vec::new(),
        }
    }
}

/// One thread's decompression context, made once a core a thread rather than once a line.
fn line_decompressor() -> zstd::bulk::Decompressor<'static> {
    // PANIC-FREE: making a context fails only where allocation does.
    zstd::bulk::Decompressor::new().expect("a zstd context")
}

impl HeldText {
    /// Compresses a plain line; a compressed line, or one already taken, stays as it is.
    ///
    /// **The compressed bytes are copied out at their own length.** zstd compresses into a buffer
    /// of its worst case, at least the line's own length, and a vector keeps its capacity: held as
    /// zstd returns it, a compressed line would hold as much memory as the plain one.
    fn compress(&mut self, with: &mut LineCompressor) {
        let Self::Plain(line) = self else {
            return;
        };
        if line.is_empty() {
            return;
        }
        with.scratch.clear();
        with.scratch
            .reserve(zstd::zstd_safe::compress_bound(line.len()));
        // PANIC-FREE: the buffer holds zstd's bound for the line, so compressing fails only where
        // allocation does.
        with.compressor
            .compress_to_buffer(line.as_slice(), &mut with.scratch)
            .expect("a line in memory compresses");
        *self = Self::Compressed {
            bytes: with.scratch.as_slice().to_vec(),
            length: line.len(),
        };
    }

    /// The line as it arrived, leaving the text empty.
    fn take(&mut self, with: &mut zstd::bulk::Decompressor<'static>) -> Vec<u8> {
        match std::mem::replace(self, Self::Plain(Vec::new())) {
            Self::Plain(line) => line,
            // PANIC-FREE: the bytes are this stage's own compression of a line of that length.
            Self::Compressed { bytes, length } => with
                .decompress(&bytes, length)
                .expect("a line this stage compressed inflates"),
        }
    }
}

/// **A line handed to the stage and not yet taken**: it waits to be read into its site with the
/// lines that arrived beside it ([`LINES_READ_TOGETHER`]).
struct ArrivedLine {
    place: RecordPlace,
    line: Vec<u8>,
    fragments: RecordFragments,
}

/// A line read into its walk record, if it has one, with its fragments packed beside it.
type ReadLine = (
    RecordPlace,
    Vec<u8>,
    Option<(SiteRecord, PackedRecordFragments)>,
);

/// One line the stage holds.
struct HeldLine {
    place: RecordPlace,
    text: HeldText,
    /// The line's walk record, as an index counted from the run's first site.
    site: Option<usize>,
}

/// **A sample's last written block** (phase output spec §2.4): its phase set, and its records
/// written within a margin of its last one, ascending.
#[derive(Clone, Debug)]
struct WrittenBlock {
    phase_set: u64,
    records: Vec<WrittenRecord>,
}

/// One record of a [`WrittenBlock`]: its run-wide site number, its position, and its phase.
#[derive(Clone, Copy, Debug)]
struct WrittenRecord {
    site: usize,
    position: u64,
    phase: WrittenPhase,
}

/// **The pass, fed one line at a time.**
pub struct NeighbourEvidenceStage {
    samples: usize,
    /// Each contig's length, by contig id, for §4.2's rule at a contig's end.
    contig_lengths: Vec<u64>,
    sizes: WindowSizes,
    path_agreement: PathAgreement,
    reads_in_draw: bool,
    fits: Vec<SampleFit>,
    contig: Option<ContigId>,
    /// Lines handed in and not yet taken, in arrival order, and the bytes of their text and read
    /// alleles ([`BYTES_READ_TOGETHER`]).
    arrived: Vec<ArrivedLine>,
    arrived_bytes: usize,
    /// How many lines wait before they are read: [`LINES_READ_TOGETHER`], or fewer in the tests
    /// that check the output does not depend on it.
    lines_read_together: usize,
    held: Vec<HeldLine>,
    /// The held lines' positions, kept beside them for the window planner.
    positions: Vec<u64>,
    /// What the held lines cost a window, as running totals from the first held line, one longer
    /// than the lines ([`WindowBytes`]): over every line a window spans, and again over its core.
    held_bytes: Vec<u64>,
    core_bytes: Vec<u64>,
    /// Of those, as running totals likewise, what `held_bytes` counts for the rounds' haplotype
    /// store and what `core_bytes` counts for the written results, which a window never holds
    /// together ([`StoreAndResults`]).
    store_bytes: Vec<u64>,
    results_bytes: Vec<u64>,
    /// **What each held line holds while no window reads it**, as a running total like `held_bytes`:
    /// what `held_bytes` counts for it less its solve state, which only a window holds, and, for a
    /// site held packed, its packed bytes in place of its numbers and reads. A window pays it for
    /// every held line it does not span ([`WindowBytes::outside`]): the stage keeps the lines
    /// within the configured margin of the next core, and those read past a window, whatever
    /// window the budget plans.
    resident_bytes: Vec<u64>,
    /// With thinned margins (pass spec §6, thinned margins): whether each held line is a record the
    /// far margins walk, and, as running totals like `held_bytes`, for each line the far margins
    /// skip, what `held_bytes` counts for its solve state and its numbers and reads
    /// (`skipped_bytes`), and the bytes it holds packed instead (`packed_bytes`), known once its
    /// line is compressed. Empty without thinning.
    far_walked: Vec<bool>,
    skipped_bytes: Vec<u64>,
    packed_bytes: Vec<u64>,
    /// With thinned margins: of `skipped_bytes`, the rounds' store, as a running total likewise.
    store_skipped_bytes: Vec<u64>,
    /// With thinned margins: the held lines some sample's last written block holds a record at
    /// ([`Self::joined_lines`]), and `skipped_bytes`, `packed_bytes` and `store_skipped_bytes`
    /// less those lines', which a
    /// window walks ([`totals_besides`]) — what a window is planned with. **Kept, not made for
    /// each line pushed**: they change only when a window is solved or held lines are dropped, and
    /// are made again then; a pushed line is in no block, so it only extends the totals. Made for
    /// every line, on 300,000 records of 10 samples at the default density they took 43 to 55 s
    /// against 8 to 12 s kept (code review 2026-11-24). Empty without thinning.
    joined: Vec<usize>,
    skipped_besides_joined: Vec<u64>,
    packed_besides_joined: Vec<u64>,
    store_skipped_besides_joined: Vec<u64>,
    /// Which records the far margins walk, decided as each arrives.
    thinner: Option<MarginThinner>,
    /// The held sites, each with its fragments; the first is site number `first_site`.
    sites: HeldSites,
    first_site: usize,
    next_site: usize,
    /// The first held line not yet written.
    core_start: usize,
    /// Held lines before this one are compressed or written.
    compressed_up_to: usize,
    tally: RewriteTally,
    /// Each sample's fill check, summed over the cores solved so far.
    fill_checks: Vec<FillCheck>,
    /// What reads did to each sample's phase, summed over the cores solved so far.
    read_counts: Vec<PassReadCounts>,
    /// Each sample's last written block on the current contig, kept across windows that write no
    /// phase for it (phase output spec §2.4).
    last_blocks: Vec<Option<WrittenBlock>>,
    /// The first position of every core solved, for the tests that score what crosses a core's
    /// edge.
    #[cfg(test)]
    core_starts: Vec<u64>,
    /// Every window planned, for the tests that read what the budget counted.
    #[cfg(test)]
    plans: Vec<super::windows::WindowPlan>,
    /// The last window planned and solved ([`Self::last_window`]).
    last_window: Option<super::windows::WindowPlan>,
    /// What the stage stopped on, if it did: it then takes no more lines ([`StageStop::Stopped`]).
    stopped: Option<StopCause>,
    /// The memory limit each window is checked against, and whether to go on past it.
    memory_guard: MemoryGuard,
    /// The highest peak expected of a window so far: the first window, and each window above
    /// every one before it, is logged.
    highest_expected_peak: u64,
}

impl NeighbourEvidenceStage {
    /// A stage for a run that decided to run the pass, over contigs of these lengths.
    #[must_use]
    pub fn new(decision: &NeighbourEvidenceDecision, contig_lengths: Vec<u64>) -> Self {
        Self {
            samples: decision.inbreeding_by_sample.len(),
            contig_lengths,
            sizes: decision.sizes,
            path_agreement: decision.path_agreement,
            reads_in_draw: decision.reads_in_draw,
            fits: decision
                .inbreeding_by_sample
                .iter()
                .map(|&inbreeding| SampleFit::new(inbreeding))
                .collect(),
            contig: None,
            arrived: Vec::new(),
            arrived_bytes: 0,
            lines_read_together: LINES_READ_TOGETHER,
            held: Vec::new(),
            positions: Vec::new(),
            held_bytes: vec![0],
            core_bytes: vec![0],
            store_bytes: vec![0],
            results_bytes: vec![0],
            resident_bytes: vec![0],
            far_walked: Vec::new(),
            skipped_bytes: vec![0],
            packed_bytes: vec![0],
            store_skipped_bytes: vec![0],
            joined: Vec::new(),
            skipped_besides_joined: vec![0],
            packed_besides_joined: vec![0],
            store_skipped_besides_joined: vec![0],
            thinner: decision.sizes.thinning.map(MarginThinner::new),
            sites: HeldSites::default(),
            first_site: 0,
            next_site: 0,
            core_start: 0,
            compressed_up_to: 0,
            tally: RewriteTally::default(),
            fill_checks: vec![FillCheck::default(); decision.inbreeding_by_sample.len()],
            read_counts: vec![PassReadCounts::default(); decision.inbreeding_by_sample.len()],
            last_blocks: vec![None; decision.inbreeding_by_sample.len()],
            #[cfg(test)]
            core_starts: Vec::new(),
            #[cfg(test)]
            plans: Vec::new(),
            last_window: None,
            stopped: None,
            memory_guard: decision.memory_guard,
            highest_expected_peak: 0,
        }
    }

    /// What the pass has done so far.
    #[must_use]
    pub fn tally(&self) -> RewriteTally {
        self.tally
    }

    /// **The last window solved**, as the budget planned it: its records and what it counted for
    /// them. For the probes that measure the budget's count against what a run holds.
    #[must_use]
    pub fn last_window(&self) -> Option<super::windows::WindowPlan> {
        self.last_window
    }

    /// Each sample's fill check so far, in the VCF's sample order.
    #[must_use]
    pub fn fill_checks(&self) -> &[FillCheck] {
        &self.fill_checks
    }

    /// What reads did to each sample's phase so far, in the VCF's sample order (read phase spec
    /// §2.7, §2.8).
    #[must_use]
    pub fn read_counts(&self) -> &[PassReadCounts] {
        &self.read_counts
    }

    /// **Takes one line**, its place already checked against the line before it by the writer,
    /// and appends to `ready` every line whose window is complete. `fragments` are the line's reads
    /// that reach another record (read phase spec §2.4). Lines are taken a batch at a time
    /// ([`Self::push_writing`]), so a window is solved, and its lines appended, when the batch
    /// holding the line that completes it is taken, or at [`Self::finish`].
    ///
    /// # Errors
    ///
    /// If a window this line completed holds more than the budget even at the smallest size, or is
    /// expected to need more memory than the run may use ([`StageStop::OverMemory`]).
    pub fn push(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        fragments: RecordFragments,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) -> Result<(), VectorStageStop> {
        self.push_writing(place, line, fragments, &mut |place, line| {
            ready.push((place, line));
            Ok(())
        })
    }

    /// **Takes one line, as [`Self::push`], and hands each line it completes to `write`** in the
    /// order they arrived, a chunk of a core at a time ([`LINES_REWRITTEN_TOGETHER`]) — so a core's
    /// rewritten lines are never all held at once.
    ///
    /// **A line is taken a little after it arrives.** Lines wait until [`LINES_READ_TOGETHER`]
    /// have arrived, or their text and read alleles reach [`BYTES_READ_TOGETHER`], then are read
    /// into sites together on the pool and taken one at a time in arrival order, exactly as each
    /// would have been taken on arrival: the running totals, the margin thinner and the window
    /// plans see the same lines in the same order, so the output does not depend on how many wait.
    /// [`Self::finish_writing`] takes those still waiting.
    ///
    /// # Errors
    ///
    /// If a window this line, or one that arrived before it and was still waiting, completed holds
    /// more than the budget even at the smallest size, or `write` refuses a line. **A stopped
    /// stage takes no more lines**: handed one, it gives the same budget stop again, or
    /// [`StageStop::Stopped`] after a refused write.
    pub fn push_writing<E>(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        fragments: RecordFragments,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        self.refuse_if_stopped()?;
        self.arrived_bytes += line.len() + fragments.len() * std::mem::size_of::<FragmentAllele>();
        self.arrived.push(ArrivedLine {
            place,
            line,
            fragments,
        });
        if self.arrived.len() >= self.lines_read_together
            || self.arrived_bytes >= BYTES_READ_TOGETHER
        {
            let taken = self.take_arrived_lines(write);
            self.note_stop(taken)?;
        }
        Ok(())
    }

    /// The stop the stage already gave, if it gave one.
    fn refuse_if_stopped<E>(&self) -> Result<(), StageStop<E>> {
        match self.stopped {
            None => Ok(()),
            Some(StopCause::OverBudget(over)) => Err(StageStop::OverBudget(over)),
            Some(StopCause::OverMemory(too_close)) => Err(StageStop::OverMemory(too_close)),
            Some(StopCause::Write) => Err(StageStop::Stopped),
        }
    }

    /// Keeps what `outcome` stopped on, if it stopped, and passes it on.
    fn note_stop<E>(&mut self, outcome: Result<(), StageStop<E>>) -> Result<(), StageStop<E>> {
        if let Err(stop) = &outcome {
            self.stopped = Some(match stop {
                StageStop::OverBudget(over) => StopCause::OverBudget(*over),
                StageStop::OverMemory(too_close) => StopCause::OverMemory(*too_close),
                StageStop::Write(_) | StageStop::Stopped => StopCause::Write,
            });
        }
        outcome
    }

    /// **Reads every waiting line into its site on the pool, then takes them one at a time**, in
    /// the order they arrived.
    fn take_arrived_lines<E>(
        &mut self,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        let samples = self.samples;
        let read: Vec<ReadLine> = std::mem::take(&mut self.arrived)
            .into_par_iter()
            .map(|arrived| {
                let site = site_of(&arrived.line, samples)
                    .map(|record| (record, PackedRecordFragments::from(arrived.fragments)));
                (arrived.place, arrived.line, site)
            })
            .collect();
        self.arrived_bytes = 0;
        for (place, line, site) in read {
            self.take_line(place, line, site, write)?;
        }
        Ok(())
    }

    /// **Takes one line already read into its site**, as [`Self::push_writing`] describes.
    fn take_line<E>(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        site: Option<(SiteRecord, PackedRecordFragments)>,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        if self.contig.is_some_and(|contig| contig != place.at.contig) {
            self.finish_contig(write)?;
        }
        self.contig = Some(place.at.contig);
        // **A line is counted at its compressed length** over every window that spans it. Until it
        // is compressed it counts at the length it arrived at, which only overstates, and no window
        // is planned on a line not yet compressed (`release`), so what a record costs does not
        // depend on when its batch was compressed. A core's lines are rewritten and handed on a
        // chunk at a time, which the budget does not count, as it does not count the walks' arrays.
        let mut held = line.len() as u64;
        let mut core = 0;
        let mut skipped = 0;
        let (mut store, mut results, mut store_skipped) = (0, 0, 0);
        let mut far_walked = false;
        let site = site.map(|(record, fragments)| {
            let (record_held, record_core) = record_window_bytes(&record);
            held += record_held + fragments.heap_bytes() as u64;
            core = record_core;
            let samples = record.sample_count() as u64;
            store = samples * round_store_bytes(record.allele_count);
            results = samples * written_result_bytes(record.allele_count)
                + solution_chunk_bytes(record.sample_count(), record.allele_count);
            if let Some(thinner) = &mut self.thinner {
                far_walked = record.is_site() && thinner.walks(record.position);
                if !far_walked {
                    store_skipped = store;
                    // Its solve state, and its numbers and reads, which it holds packed instead
                    // while no walk reads it ([`HeldSites::pack`]).
                    skipped = record.sample_count() as u64 * held_solve_bytes(record.allele_count)
                        + (record.heap_bytes() - record.sample_count()) as u64
                        + fragments.heap_bytes() as u64;
                }
            }
            self.sites.push(record, fragments);
            self.next_site += 1;
            self.next_site - 1
        });
        let last = |totals: &[u64]| totals.last().copied().unwrap_or(0);
        self.held_bytes.push(last(&self.held_bytes) + held);
        self.core_bytes.push(last(&self.core_bytes) + core);
        self.store_bytes.push(last(&self.store_bytes) + store);
        self.results_bytes.push(last(&self.results_bytes) + results);
        if self.thinner.is_some() {
            self.far_walked.push(far_walked);
            self.skipped_bytes.push(last(&self.skipped_bytes) + skipped);
            // A new line is in no last block: blocks hold only records of cores already written.
            self.skipped_besides_joined
                .push(last(&self.skipped_besides_joined) + skipped);
            self.store_skipped_bytes
                .push(last(&self.store_skipped_bytes) + store_skipped);
            self.store_skipped_besides_joined
                .push(last(&self.store_skipped_besides_joined) + store_skipped);
            // Its packed bytes are counted once it is packed, when its line is compressed.
            self.packed_bytes.push(last(&self.packed_bytes));
            self.packed_besides_joined
                .push(last(&self.packed_besides_joined));
        }
        self.positions.push(place.at.position.0);
        self.held.push(HeldLine {
            place,
            text: HeldText::Plain(line),
            site,
        });
        self.count_resident_bytes_from(self.held.len() - 1);
        if self.held.len() - self.compressed_up_to.max(self.core_start) >= LINES_COMPRESSED_TOGETHER
        {
            self.compress_held_lines();
        }
        self.release(false, write)
    }

    /// **Compresses every held line not yet compressed**, on the pool, and counts each at its
    /// compressed length from then on.
    fn compress_held_lines(&mut self) {
        let from = self.compressed_up_to.max(self.core_start);
        self.held[from..]
            .par_iter_mut()
            .for_each_init(LineCompressor::new, |with, held| held.text.compress(with));
        // With thinned margins, the sites the far margins skip are held compressed too, until a
        // window walks them (pass spec §6, thinned margins). Lines with no site among them have
        // nothing to pack, and their empty span, `0..0`, is no range of held sites past the first
        // contig.
        let sites = self.site_span(from..self.held.len());
        if self.thinner.is_some() && !sites.is_empty() {
            let first_site = self.first_site;
            let held = &self.held;
            let far_walked = &self.far_walked;
            let line_of_site: Vec<usize> = (from..held.len())
                .filter(|&line| held[line].site.is_some())
                .collect();
            let start = sites.start;
            self.sites
                .pack(sites.start - first_site..sites.end - first_site, |site| {
                    !far_walked[line_of_site[site + first_site - start]]
                });
        }
        // The running totals from the first compressed line on, less what each line saved. **The
        // saving is signed**: zstd makes a line under about 60 bytes longer, by its frame header.
        let mut saved: i64 = 0;
        for (index, held) in self.held.iter().enumerate().skip(from) {
            if let HeldText::Compressed { bytes, length } = &held.text {
                saved += *length as i64 - bytes.len() as i64;
            }
            let total = &mut self.held_bytes[index + 1];
            // PANIC-FREE: the total less the saving is the same sum with each line at its
            // compressed length, which is not negative.
            *total = total.checked_add_signed(-saved).expect(
                "a running total stays at least the bytes its lines hold, compressed or not",
            );
        }
        // What the sites packed hold so, from the same line on. The lines from `from` were pushed
        // after every block's records, so none is joined.
        if self.thinner.is_some() {
            for totals in [&mut self.packed_bytes, &mut self.packed_besides_joined] {
                totals.truncate(from + 1);
            }
            for line in from..self.held.len() {
                let own = self.held[line]
                    .site
                    .map_or(0, |site| self.sites.packed_length[site - self.first_site])
                    as u64;
                self.packed_bytes.push(self.packed_bytes[line] + own);
                self.packed_besides_joined
                    .push(self.packed_besides_joined[line] + own);
            }
        }
        self.compressed_up_to = self.held.len();
        self.count_resident_bytes_from(from);
    }

    /// **Finishes the current contig**: solves its remaining windows and appends their lines to
    /// `ready`. Called at a contig change and at the end of the input.
    ///
    /// # Errors
    ///
    /// As [`Self::push`].
    pub fn finish(
        &mut self,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) -> Result<(), VectorStageStop> {
        self.finish_writing(&mut |place, line| {
            ready.push((place, line));
            Ok(())
        })
    }

    /// **Finishes the current contig, as [`Self::finish`], handing its lines to `write`.**
    ///
    /// # Errors
    ///
    /// As [`Self::push_writing`].
    pub fn finish_writing<E>(
        &mut self,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        self.refuse_if_stopped()?;
        let finished = self
            .take_arrived_lines(write)
            .and_then(|()| self.finish_contig(write));
        self.note_stop(finished)
    }

    /// **Finishes the current contig** once every line of it handed in has been taken.
    fn finish_contig<E>(
        &mut self,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        self.release(true, write)?;
        self.held.clear();
        self.held_bytes = vec![0];
        self.core_bytes = vec![0];
        self.store_bytes = vec![0];
        self.results_bytes = vec![0];
        self.resident_bytes = vec![0];
        self.far_walked.clear();
        self.skipped_bytes = vec![0];
        self.packed_bytes = vec![0];
        self.store_skipped_bytes = vec![0];
        self.joined.clear();
        self.skipped_besides_joined = vec![0];
        self.packed_besides_joined = vec![0];
        self.store_skipped_besides_joined = vec![0];
        if let Some(thinner) = &mut self.thinner {
            *thinner = MarginThinner::new(thinner.thinning);
        }
        self.compressed_up_to = 0;
        self.positions.clear();
        self.sites.clear();
        self.first_site = self.next_site;
        self.core_start = 0;
        self.contig = None;
        // A block never crosses a contig.
        self.last_blocks.fill(None);
        Ok(())
    }

    /// Solves every window whose lines have all arrived (or all of them, when the contig is
    /// complete), hands back their cores, and forgets what no later window reads.
    fn release<E>(
        &mut self,
        contig_complete: bool,
        write: &mut impl FnMut(RecordPlace, Vec<u8>) -> Result<(), E>,
    ) -> Result<(), StageStop<E>> {
        while self.core_start < self.held.len() {
            // With thinned margins, a far-margin record some sample's last block joins at is
            // walked, so its solve state is counted (pass spec §6, thinned margins).
            #[cfg(test)]
            self.check_the_kept_totals();
            let Some(plan) = plan_window(
                &self.positions,
                WindowBytes {
                    held: &self.held_bytes,
                    core: &self.core_bytes,
                    outside: Some(&self.resident_bytes),
                    skipped: self
                        .thinner
                        .is_some()
                        .then_some(self.skipped_besides_joined.as_slice()),
                    packed: self
                        .thinner
                        .is_some()
                        .then_some(self.packed_besides_joined.as_slice()),
                    store_and_results: Some(StoreAndResults {
                        store: &self.store_bytes,
                        store_skipped: self
                            .thinner
                            .is_some()
                            .then_some(self.store_skipped_besides_joined.as_slice()),
                        results: &self.results_bytes,
                    }),
                },
                self.core_start,
                contig_complete,
                &self.sizes,
            ) else {
                return Ok(());
            };
            if self.compressed_up_to < self.held.len() {
                // Lines not yet compressed count at the length they arrived at, more than they will
                // hold; counted compressed, the window may be larger, or need lines not yet read.
                self.compress_held_lines();
                continue;
            }
            if plan.over_budget > 0 {
                return Err(StageStop::OverBudget(WindowOverBudget {
                    contig: self.contig.map_or(0, |contig| contig.0),
                    position: self.positions[plan.core_start],
                    records: plan.records(),
                    samples: self.samples,
                    bytes: plan.bytes,
                    budget: self.sizes.budget_bytes,
                }));
            }
            #[cfg(test)]
            self.core_starts.push(self.positions[plan.core_start]);
            #[cfg(test)]
            self.plans.push(plan);
            self.last_window = Some(plan);
            let began = std::time::Instant::now();
            let mut walked_records = None;
            let window_sites = self.site_span(plan.context_start..plan.context_end);
            let core_sites = self.site_span(plan.core_start..plan.core_end);
            let refined = if core_sites.is_empty() {
                Vec::new()
            } else {
                let contig_length = self
                    .contig
                    .and_then(|contig| self.contig_lengths.get(contig.0 as usize).copied());
                // The window's contig and its first core site's position key its draws: they depend
                // on the window alone, not on what else the run called.
                let first_core = &self.sites.records[core_sites.start - self.first_site];
                let key = window_key(
                    self.contig.map_or(0, |contig| contig.0),
                    first_core.position,
                );
                // Each sample's last written block, at its records in this window's left margin.
                let joins: Vec<Option<LastBlock>> = self
                    .last_blocks
                    .iter()
                    .map(|last| {
                        last.as_ref().map(|last| LastBlock {
                            phase_set: last.phase_set,
                            records: last
                                .records
                                .iter()
                                .filter(|record| {
                                    (window_sites.start..core_sites.start).contains(&record.site)
                                })
                                .map(|record| JoinRecord {
                                    record: record.site - window_sites.start,
                                    copy_0_allele: record.phase.copy_0_allele,
                                    copy_1_allele: record.phase.copy_1_allele,
                                })
                                .collect(),
                        })
                    })
                    .collect();
                let in_window =
                    window_sites.start - self.first_site..window_sites.end - self.first_site;
                let walked = self.thinner.is_some().then(|| {
                    let mut walked = vec![false; in_window.len()];
                    for (line, held) in self
                        .held
                        .iter()
                        .enumerate()
                        .take(plan.context_end)
                        .skip(plan.context_start)
                    {
                        if let Some(site) = held.site {
                            walked[site - window_sites.start] = (plan.near_start..plan.near_end)
                                .contains(&line)
                                || self.far_walked[line]
                                || self.joined.binary_search(&line).is_ok();
                        }
                    }
                    walked
                });
                // The window's walks read the numbers and reads of the records they walk, and only
                // those: the others are held compressed.
                if let Some(walked) = &walked {
                    self.sites.hold_whole(in_window.clone(), walked);
                }
                let (records, fragments) = self.sites.window(in_window);
                walked_records = walked.as_ref().map(|walked| {
                    walked
                        .iter()
                        .zip(records)
                        .filter(|(walked, record)| **walked && record.is_site())
                        .count()
                });
                // **Checked before it is solved**: what the budget counted and the records the
                // walks read are known now, and the walks' arrays are not yet made.
                let walked_count = walked_records
                    .unwrap_or_else(|| records.iter().filter(|record| record.is_site()).count());
                check_the_windows_memory(
                    self.memory_guard,
                    &mut self.highest_expected_peak,
                    PeakEstimate::for_window(
                        plan.bytes,
                        walked_count as u64,
                        self.samples,
                        rayon::current_num_threads(),
                        self.sizes.budget_bytes,
                    ),
                    PeakMoment::Window {
                        contig: self.contig.map_or(0, |contig| contig.0),
                        position: self.positions[plan.core_start],
                    },
                )
                .map_err(StageStop::OverMemory)?;
                let solved = refine_window_with_reads(
                    records,
                    Some(fragments),
                    walked.as_deref(),
                    core_sites.start - window_sites.start..core_sites.end - window_sites.start,
                    contig_length,
                    &mut self.fits,
                    MAX_TEMPLATES,
                    NoReadRule::default(),
                    Phasing::rounds_for(key)
                        .with_path_agreement(self.path_agreement)
                        .with_reads_in_draw(self.reads_in_draw),
                    &joins,
                );
                for (sum, check) in self.fill_checks.iter_mut().zip(solved.checks) {
                    *sum += check;
                }
                for (sum, counts) in self.read_counts.iter_mut().zip(solved.read_counts) {
                    *sum += counts;
                }
                self.remember_last_blocks(core_sites.start, &solved.records);
                solved.records
            };
            // **The core's lines are rewritten on the pool and handed on in their own order**, a
            // chunk at a time. Each rewrite reads only its line and its refined record, so the lines
            // are independent; done on this thread it was a quarter of the writing pass on 59
            // tomato accessions over 32 Mb, with every other thread waiting for the next window.
            // The tallies are counts, so summing them in any order gives the serial total.
            // What the core's lines are counted at, which they no longer hold once written.
            let written_text: Vec<u64> = self.held[plan.core_start..plan.core_end]
                .iter()
                .map(|held| held.text.held_bytes() as u64)
                .collect();
            let mut chunk_start = plan.core_start;
            while chunk_start < plan.core_end {
                let chunk_end = (chunk_start + LINES_REWRITTEN_TOGETHER).min(plan.core_end);
                let rewritten: Vec<(RecordPlace, Option<Vec<u8>>, RewriteTally)> = self.held
                    [chunk_start..chunk_end]
                    .par_iter_mut()
                    .map_init(
                        || (GenotypePairs::default(), line_decompressor()),
                        |(pairs, decompressor), held| {
                            let line = held.text.take(decompressor);
                            let mut tally = RewriteTally::default();
                            let written = match held.site {
                                Some(site) => rewrite_line(
                                    &line,
                                    &refined[site - core_sites.start],
                                    pairs,
                                    &mut tally,
                                ),
                                None => Some(line),
                            };
                            (held.place, written, tally)
                        },
                    )
                    .collect();
                for (place, written, tally) in rewritten {
                    self.tally.add(tally);
                    if let Some(line) = written {
                        write(place, line).map_err(StageStop::Write)?;
                    }
                }
                chunk_start = chunk_end;
            }
            // **A written line is counted at nothing**: its text is let go as it is written, but
            // the line is held on, in the next windows' left margin, for its record.
            let mut freed = 0;
            for line in plan.core_start..self.held.len() {
                if line < plan.core_end {
                    freed += written_text[line - plan.core_start];
                }
                self.held_bytes[line + 1] -= freed;
            }
            // **One line a window in the run's log**: where it was, what it held by the budget's
            // count and how long it took — what a run that slows or runs short of memory is read
            // back from.
            tracing::debug!(
                "refinement window: contig {} core {}..{}, {} records ({} in the core{}), {} \
                 samples, {:.2} GB of a {:.2} GB budget, {:.1} s",
                self.contig.map_or(0, |contig| contig.0),
                self.positions[plan.core_start],
                self.positions[plan.core_end - 1],
                plan.records(),
                plan.core_end - plan.core_start,
                walked_records.map_or(String::new(), |walked| format!(", {walked} walked")),
                self.samples,
                plan.bytes as f64 / 1e9,
                self.sizes.budget_bytes as f64 / 1e9,
                began.elapsed().as_secs_f64()
            );
            self.core_start = plan.core_end;
            self.forget_what_no_window_reads();
            self.make_the_joins_again();
            // The window unpacked the skipped records it walked, and packed again those it did not.
            self.count_resident_bytes_from(0);
        }
        Ok(())
    }

    /// **Updates each sample's last written block** from a core's records, the first being site
    /// number `first_site`: the block of its last phased record, extended where it continues the
    /// one before, and cut to the records within a margin of that record — further back, no later
    /// window's margin reaches.
    fn remember_last_blocks(&mut self, first_site: usize, core: &[RefinedRecord]) {
        for (sample, last) in self.last_blocks.iter_mut().enumerate() {
            let phased: Vec<WrittenRecord> = core
                .iter()
                .enumerate()
                .filter_map(|(offset, record)| {
                    record.phase(sample).map(|phase| {
                        let site = first_site + offset;
                        WrittenRecord {
                            site,
                            position: self.sites.records[site - self.first_site].position,
                            phase,
                        }
                    })
                })
                .collect();
            let Some(&final_record) = phased.last() else {
                continue;
            };
            let phase_set = final_record.phase.phase_set;
            let mut records = match last.take() {
                Some(previous) if previous.phase_set == phase_set => previous.records,
                _ => Vec::new(),
            };
            records.extend(
                phased
                    .iter()
                    .filter(|record| record.phase.phase_set == phase_set),
            );
            let reach = final_record
                .position
                .saturating_sub(self.sizes.margin_bases);
            records.retain(|record| record.position >= reach);
            *last = Some(WrittenBlock { phase_set, records });
        }
    }

    /// **The held lines some sample's last written block holds a record at**, ascending: the
    /// records a window's walks join at (phase output spec §2.4). Empty without thinning, where
    /// every record is walked anyway.
    fn joined_lines(&self) -> Vec<usize> {
        if self.thinner.is_none() || self.core_start == 0 {
            return Vec::new();
        }
        // Held sites are consecutive numbers from `first_site`, each on one line.
        let mut line_of_site = vec![0; self.next_site - self.first_site];
        for (line, held) in self.held.iter().enumerate() {
            if let Some(site) = held.site {
                line_of_site[site - self.first_site] = line;
            }
        }
        let mut joined: Vec<usize> = self
            .last_blocks
            .iter()
            .flatten()
            .flat_map(|block| &block.records)
            .filter(|record| record.site >= self.first_site)
            .map(|record| line_of_site[record.site - self.first_site])
            .filter(|&line| !self.far_walked[line])
            .collect();
        joined.sort_unstable();
        joined.dedup();
        joined
    }

    /// **What held line `line` holds while no window reads it**: what `held_bytes` counts for it,
    /// less the solve state a window would add for a site, and, for a site whose numbers and reads
    /// are held packed, its packed bytes in place of those — which leaves its line, its rows and
    /// its packed bytes.
    fn resident_bytes_of(&self, line: usize) -> u64 {
        let held = self.held_bytes[line + 1] - self.held_bytes[line];
        let Some(site) = self.held[line].site else {
            return held;
        };
        let index = site - self.first_site;
        if self.sites.whole[index] {
            let record = &self.sites.records[index];
            held - record.sample_count() as u64 * held_solve_bytes(record.allele_count)
        } else {
            // Packed only where the far margins skip it, so `skipped_bytes` counts it: its solve
            // state, numbers and reads.
            held - (self.skipped_bytes[line + 1] - self.skipped_bytes[line])
                + (self.packed_bytes[line + 1] - self.packed_bytes[line])
        }
    }

    /// **Counts `resident_bytes` again from held line `from` on**, the totals before it kept.
    fn count_resident_bytes_from(&mut self, from: usize) {
        self.resident_bytes.truncate(from + 1);
        for line in from..self.held.len() {
            let total = self.resident_bytes[line] + self.resident_bytes_of(line);
            self.resident_bytes.push(total);
        }
    }

    /// **Makes the kept joins and their totals again** (`joined`, `skipped_besides_joined`), once
    /// a window is solved: its core may have changed the last written blocks, and lines may
    /// have been dropped.
    fn make_the_joins_again(&mut self) {
        if self.thinner.is_none() {
            return;
        }
        self.joined = self.joined_lines();
        self.skipped_besides_joined = totals_besides(&self.skipped_bytes, &self.joined);
        self.packed_besides_joined = totals_besides(&self.packed_bytes, &self.joined);
        self.store_skipped_besides_joined = totals_besides(&self.store_skipped_bytes, &self.joined);
    }

    /// **The kept joins and resident totals are what making them again gives**: what lets them be
    /// kept rather than made for every line, checked before every window is planned in every test that runs a
    /// stage.
    #[cfg(test)]
    fn check_the_kept_totals(&self) {
        let mut resident = vec![0];
        for line in 0..self.held.len() {
            resident.push(resident[line] + self.resident_bytes_of(line));
        }
        assert_eq!(
            self.resident_bytes, resident,
            "the kept resident totals drifted"
        );
        if self.thinner.is_none() {
            // Without thinning every site is held whole, and each line is counted at the text it
            // holds — none once written — and its record's numbers, reads and solve state.
            for (line, held) in self.held.iter().enumerate() {
                let record = held.site.map_or(0, |site| {
                    let index = site - self.first_site;
                    record_window_bytes(&self.sites.records[index]).0
                        + self.sites.fragments[index].heap_bytes() as u64
                });
                assert_eq!(
                    self.held_bytes[line + 1] - self.held_bytes[line],
                    held.text.held_bytes() as u64 + record,
                    "line {line}, the core starting at {}",
                    self.core_start
                );
            }
            return;
        }
        assert_eq!(self.joined, self.joined_lines(), "the kept joins drifted");
        assert_eq!(
            self.skipped_besides_joined,
            totals_besides(&self.skipped_bytes, &self.joined),
            "the kept skipped totals drifted"
        );
        assert_eq!(
            self.packed_besides_joined,
            totals_besides(&self.packed_bytes, &self.joined),
            "the kept packed totals drifted"
        );
        assert_eq!(
            self.store_skipped_besides_joined,
            totals_besides(&self.store_skipped_bytes, &self.joined),
            "the kept skipped store totals drifted"
        );
        // A site is held in one form at a time, and a site held packed holds the bytes it is
        // counted at: packing a record again gives the bytes its first packing gave.
        for index in 0..self.sites.records.len() {
            let packed = &self.sites.packed[index];
            assert_eq!(packed.is_some(), !self.sites.whole[index], "site {index}");
            if let Some(packed) = packed {
                assert_eq!(
                    packed.len(),
                    self.sites.packed_length[index],
                    "site {index}"
                );
            }
        }
    }

    /// The sites among held lines `lines`, as a range of site numbers.
    fn site_span(&self, lines: std::ops::Range<usize>) -> std::ops::Range<usize> {
        let first = self.held[lines.clone()].iter().find_map(|held| held.site);
        let last = self.held[lines].iter().rev().find_map(|held| held.site);
        match (first, last) {
            (Some(first), Some(last)) => first..last + 1,
            _ => 0..0,
        }
    }

    /// Drops held lines further back than a margin before the next core's first line.
    fn forget_what_no_window_reads(&mut self) {
        let next = self
            .positions
            .get(self.core_start)
            .or(self.positions.last())
            .copied()
            .unwrap_or(0);
        let keep_from = self.positions[..self.core_start]
            .partition_point(|&position| position < next.saturating_sub(self.sizes.margin_bases));
        if keep_from == 0 {
            return;
        }
        let sites_dropped = self.held[..keep_from]
            .iter()
            .filter(|held| held.site.is_some())
            .count();
        self.held.drain(..keep_from);
        if self.thinner.is_some() {
            self.far_walked.drain(..keep_from);
        }
        for totals in [
            &mut self.held_bytes,
            &mut self.core_bytes,
            &mut self.skipped_bytes,
            &mut self.packed_bytes,
            &mut self.store_bytes,
            &mut self.results_bytes,
            &mut self.store_skipped_bytes,
        ] {
            if totals.len() == 1 {
                continue;
            }
            let dropped = totals[keep_from];
            totals.drain(..keep_from);
            for total in totals.iter_mut() {
                *total -= dropped;
            }
        }
        self.compressed_up_to = self.compressed_up_to.saturating_sub(keep_from);
        self.positions.drain(..keep_from);
        self.sites.drop_first(sites_dropped);
        self.first_site += sites_dropped;
        self.core_start -= keep_from;
    }
}

/// **Checks a planned window's expected peak against the run's memory limit** (pass spec §6,
/// memory), before it is solved. The first window's estimate is logged, and each one higher than
/// every window's before it in the debug log; a window above the limit's stop share stops the
/// stage, unless `--refine-ignore-memory-limit` was given, when each such new highest is a
/// warning. Without a limit nothing is stopped, and the estimates are logged the same way.
fn check_the_windows_memory(
    guard: MemoryGuard,
    highest_expected_peak: &mut u64,
    estimate: PeakEstimate,
    moment: PeakMoment,
) -> Result<(), PeakTooClose> {
    let peak = estimate.peak_bytes();
    let scale = |budget: u64| {
        // A window counted at another budget holds, and walks, as much more or less of it.
        let at = |bytes: u64| {
            u64::try_from(
                u128::from(bytes) * u128::from(budget) / u128::from(estimate.budget_bytes.max(1)),
            )
            .unwrap_or(u64::MAX)
        };
        PeakEstimate {
            counted_bytes: at(estimate.counted_bytes),
            walked_records: at(estimate.walked_records),
            budget_bytes: budget,
            ..estimate
        }
    };
    let found = super::peak::check(estimate, guard, moment, scale)?;
    if peak > *highest_expected_peak {
        let first = *highest_expected_peak == 0;
        *highest_expected_peak = peak;
        let words = super::peak::fits_in_words(&estimate, &found, guard.limit);
        let PeakMoment::Window { contig, position } = moment else {
            return Ok(());
        };
        let which = if first {
            "the first window"
        } else {
            "a window higher than any before it"
        };
        let line =
            format!("refinement memory, {which} (contig {contig} position {position}): {words}");
        match found {
            PeakCheck::Ignored(_) => tracing::warn!("{line}"),
            PeakCheck::Fits { .. } | PeakCheck::NoLimit if first => tracing::info!("{line}"),
            PeakCheck::Fits { .. } | PeakCheck::NoLimit => tracing::debug!("{line}"),
        }
    }
    Ok(())
}

/// **A running total over the held lines, less the lines in `joined`** (ascending): what a window
/// does not pay for the far margins' skipped records, or holds packed for them, besides the records
/// some sample's last written block joins at, which a window walks.
fn totals_besides(totals: &[u64], joined: &[usize]) -> Vec<u64> {
    if joined.is_empty() {
        return totals.to_vec();
    }
    let mut besides = Vec::with_capacity(totals.len());
    besides.push(0);
    let mut less = 0;
    let mut next_joined = joined.iter().peekable();
    for line in 0..totals.len() - 1 {
        if next_joined.next_if_eq(&&line).is_some() {
            less += totals[line + 1] - totals[line];
        }
        besides.push(totals[line + 1] - less);
    }
    besides
}

/// **Which records a window's far margins walk** (pass spec §6, thinned margins), decided once for each record as
/// it arrives, so it does not depend on which window reads the record: a token bucket that fills at
/// the thinning's density, a token a record, holding a kb's worth or one token if that is more,
/// and spends one on each record it walks. A stretch no denser than that is walked in full; a
/// denser one is walked at that density, evenly spaced, after a burst of at most a kb's worth.
/// Counted in millionths of a token, in integers, so the choice is the same on every machine.
#[derive(Clone, Copy, Debug)]
struct MarginThinner {
    thinning: MarginThinning,
    tokens: u64,
    last_position: Option<u64>,
}

/// A token, in [`MarginThinner`]'s units: at a density of `t` thousandths of a record a kb, a base
/// adds `t` units, so a kb adds `t / 1,000` tokens.
const THINNER_TOKEN: u64 = 1_000_000;

impl MarginThinner {
    /// A bucket at a contig's start, full.
    fn new(thinning: MarginThinning) -> Self {
        Self {
            thinning,
            tokens: Self::capacity(thinning),
            last_position: None,
        }
    }

    /// A kb's worth of tokens, or one token at a density under one a kb.
    fn capacity(thinning: MarginThinning) -> u64 {
        (u64::from(thinning.thousandths_per_kb.get()) * 1_000).max(THINNER_TOKEN)
    }

    /// Whether the record at `position`, the next site of the contig, is walked in far margins.
    fn walks(&mut self, position: u64) -> bool {
        let since = self
            .last_position
            .map_or(0, |last| position.saturating_sub(last));
        self.last_position = Some(position);
        self.tokens = self
            .tokens
            .saturating_add(since.saturating_mul(u64::from(self.thinning.thousandths_per_kb.get())))
            .min(Self::capacity(self.thinning));
        if self.tokens >= THINNER_TOKEN {
            self.tokens -= THINNER_TOKEN;
            true
        } else {
            false
        }
    }
}

/// **The held sites' walk records and, beside each, its fragments** (read phase spec §2.7): the
/// reads the link test reads. Two vectors, so a window hands the solve two slices without a copy,
/// but changed only together, here, so one cannot drift a site from the other. A line that is no
/// site keeps none, since the pass phases only sites.
///
/// **With thinned margins, a site the far margins skip is held compressed** (pass spec §6, thinned
/// margins): its numbers and reads, packed once its line is compressed, and taken out of the
/// record, which keeps its position, allele count and rows. A window puts them back for the sites
/// it walks ([`Self::hold_whole`]), dropping the packed bytes, and packs again those it does not
/// walk. **A site is held in one form at a time**: kept beside the whole numbers, the packed copy
/// was a twelfth of a 2,000-sample window at the 500-kb floor (memory review 2026-11-24). Packing
/// depends on the site's numbers and reads alone, so packing again gives the same bytes.
#[derive(Debug, Default)]
struct HeldSites {
    records: Vec<SiteRecord>,
    fragments: Vec<PackedRecordFragments>,
    /// Each site's numbers and reads, compressed — its uncompressed length, then zstd's frame —
    /// while no window walks it; `None` while it is held whole.
    packed: Vec<Option<Box<[u8]>>>,
    /// The bytes each site holds while packed, from its first packing: zero for a site the far
    /// margins walk, which is never packed.
    packed_length: Vec<usize>,
    /// Whether each site's numbers and reads are in `records` and `fragments`.
    whole: Vec<bool>,
}

impl HeldSites {
    fn push(&mut self, record: SiteRecord, fragments: PackedRecordFragments) {
        self.records.push(record);
        self.fragments.push(fragments);
        self.packed.push(None);
        self.packed_length.push(0);
        self.whole.push(true);
    }

    fn clear(&mut self) {
        self.records.clear();
        self.fragments.clear();
        self.packed.clear();
        self.packed_length.clear();
        self.whole.clear();
    }

    /// Forget the first `count` sites.
    fn drop_first(&mut self, count: usize) {
        self.records.drain(..count);
        self.fragments.drain(..count);
        self.packed.drain(..count);
        self.packed_length.drain(..count);
        self.whole.drain(..count);
    }

    /// **Packs the numbers and reads of the sites in `range` that `skipped` admits**, on the
    /// pool, the first time they are, and holds them so.
    fn pack(&mut self, range: std::ops::Range<usize>, skipped: impl Fn(usize) -> bool + Sync) {
        let start = range.start;
        self.records[range.clone()]
            .par_iter_mut()
            .zip(&mut self.fragments[range.clone()])
            .zip(&mut self.packed[range.clone()])
            .zip(&mut self.packed_length[range.clone()])
            .zip(&mut self.whole[range])
            .enumerate()
            .for_each_init(
                LineCompressor::new,
                |with, (at, ((((record, fragments), packed), packed_length), whole))| {
                    if !skipped(start + at) || *packed_length > 0 {
                        return;
                    }
                    let bytes = pack_site(with, record, fragments);
                    *packed_length = bytes.len();
                    *packed = Some(bytes);
                    *whole = false;
                },
            );
    }

    /// **Holds whole exactly the packed sites in `range` that `walked` admits** (one entry a site
    /// in `range`), on the pool: unpacks those it admits, letting their packed bytes go, and packs
    /// again those it does not admit that an earlier window unpacked.
    fn hold_whole(&mut self, range: std::ops::Range<usize>, walked: &[bool]) {
        self.records[range.clone()]
            .par_iter_mut()
            .zip(&mut self.fragments[range.clone()])
            .zip(&mut self.packed[range.clone()])
            .zip(&self.packed_length[range.clone()])
            .zip(&mut self.whole[range])
            .zip(walked)
            .for_each_init(
                || (line_decompressor(), LineCompressor::new()),
                |(decompressor, compressor),
                 (((((record, fragments), packed), &packed_length), whole), &walked)| {
                    if packed_length == 0 {
                        // Never packed: a site the far margins walk, held whole throughout.
                        return;
                    }
                    if walked && !*whole {
                        // PANIC-FREE: a packed site not held whole holds its packed bytes.
                        let packed = packed.take().expect("a site not held whole is packed");
                        let (length, frame) = packed.split_at(4);
                        let length =
                            u32::from_le_bytes(length.try_into().expect("four bytes")) as usize;
                        // PANIC-FREE: the stage reads back only frames it compressed, at their length.
                        let plain = decompressor
                            .decompress(frame, length)
                            .expect("a site's numbers compressed here decompress");
                        let mut bytes = plain.as_slice();
                        record.put_numbers(RecordNumbers::read_from(&mut bytes));
                        *fragments = PackedRecordFragments::read_from(&mut bytes);
                        *whole = true;
                    } else if !walked && *whole {
                        *packed = Some(pack_site(compressor, record, fragments));
                        *whole = false;
                    }
                },
            );
    }

    /// The sites in `range`, as the window solve takes them.
    fn window(&self, range: std::ops::Range<usize>) -> (&[SiteRecord], &[PackedRecordFragments]) {
        (&self.records[range.clone()], &self.fragments[range])
    }
}

/// **Packs a site's numbers and reads**, taking them out of it: their uncompressed length, then
/// zstd's frame, in an array of its own length.
fn pack_site(
    with: &mut LineCompressor,
    record: &mut SiteRecord,
    fragments: &mut PackedRecordFragments,
) -> Box<[u8]> {
    let numbers = record.take_numbers();
    let mut plain = Vec::new();
    numbers.write_to(&mut plain);
    std::mem::take(fragments).write_to(&mut plain);
    with.scratch.clear();
    with.scratch
        .reserve(zstd::zstd_safe::compress_bound(plain.len()));
    // PANIC-FREE: the buffer holds zstd's bound for the bytes, so compressing fails only where
    // allocation does.
    with.compressor
        .compress_to_buffer(plain.as_slice(), &mut with.scratch)
        .expect("a site's numbers in memory compress");
    let mut bytes = Vec::with_capacity(4 + with.scratch.len());
    bytes.extend_from_slice(
        &u32::try_from(plain.len())
            .expect("a site's numbers fit u32 bytes")
            .to_le_bytes(),
    );
    bytes.extend_from_slice(&with.scratch);
    bytes.into_boxed_slice()
}

#[cfg(test)]
mod tests;
