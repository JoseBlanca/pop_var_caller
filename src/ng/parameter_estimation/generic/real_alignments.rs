//! F3 — the identities that hold by construction, checked on real alignments.
//!
//! Every other test of this module builds its own loci, or fills a cell table directly.
//! These four walk a real alignment file through the real locus generator and ask the
//! accumulator questions whose answers do not depend on a truth set:
//!
//! 1. **The two tables are the same table on a single-library sample.** Fold the windowed
//!    histogram over its windows and it must equal the read-group one cell for cell
//!    (`spec/parameter_prepass_generic.md` §12.6). That is the property §1 leans on to
//!    avoid keeping a third object, and it is the reason the accumulator may build both
//!    without the multi-library machinery costing a single-library sample anything.
//! 2. **Sharded accumulation is exact.** The same territory walked as one stream, and
//!    again as sixteen shards of **whole typed regions** merged together, must give
//!    identical tables — which the integer merge makes an equality rather than a tolerance
//!    (`arch/parameter_prepass_generic.md` §9). Whole regions and never part of one: see
//!    [`deal_into_shards`] for the rule and the measurement behind it.
//! 3. **`loci_overlapping_previous` is zero.** The generic loci must partition the
//!    positions they cover, or a site enters the windowed table twice. A non-zero count
//!    here is a bug report against locus generation, not something this unit absorbs.
//! 4. **Not an identity, and here because it runs on the same walk** — the whole path end
//!    to end, the only one of the four that runs a fit, and the first time any fit in this
//!    module sees a read. Its console lines say `end to end:` where the others say
//!    `identity n:`.
//!
//! # What these cannot say
//!
//! **Nothing here is evidence about the four parameters this step fits** — the per-read
//! error rate, heterozygosity, the homozygous-non-reference rate and `F`. Identities 1–3
//! compare objects the walk already built against each other; the end-to-end test asserts
//! only that every rate was fitted from the sample's own sites and that none sits on an end
//! rung of the ladder. The values themselves are anchored in Milestone G, against the GIAB
//! truth set and against a coverage sweep.
//!
//! **Every site is scored at ploidy 2**, because `ConstantPloidy` is the only `PloidyMap`
//! this plan builds and nothing here overrides it. Both cohorts' BEDs are autosomal, so the
//! assumption is met — but a BED reaching a sex chromosome, or a polyploid organism, would
//! have its cells labelled with genotype classes they do not have, and identities 1–3 would
//! not notice, because they compare the same labels against each other.
//!
//! **The attributed arm of the cell key is not exercised.** Every sample in both cohorts
//! carries one read group — checked, not assumed, by
//! [`the_two_tables_agree_cell_for_cell`] — so `multi_library` is false throughout and no
//! site is ever keyed by the library its alternative reads came from. A real multi-library
//! alignment is what would test that, and neither cohort holds one.
//!
//! **`F` is not fitted.** The runs model refuses unless
//! [`MIN_WINDOWS_TO_FIT_INBREEDING`](super::runs::MIN_WINDOWS_TO_FIT_INBREEDING) = 3,000
//! *separate* 100 kb windows each hold at least one site — `fit_inbreeding` tests
//! `windows_holding_sites`, not a base count. What disqualifies both cohorts is therefore
//! scatter rather than size: the tomato BED is 8.0 Mb in 80 spans and the GIAB one 572 kb in
//! 100, so between them they touch a couple of hundred windows. The end-to-end test
//! therefore runs with `F` supplied.
//!
//! # Running them
//!
//! `#[ignore]`d and driven by environment, so one file serves both organisms — the
//! convention `locus_generation/pileup/parity.rs` set.
//!
//! **`--release` here is only speed, unlike in `parity.rs`.** That file needs release
//! because real paired-end data trips a reachable `debug_assert!` in *production's* walker,
//! which it runs alongside ng's; nothing here runs production, and the debug build
//! completes — measured, because copying the neighbouring justification would have stated a
//! failure mode this path does not have. What it costs is time:
//! `no_locus_overlaps_the_one_before_it` on the HG002 30x arm takes **12.4 s in release and
//! 128.7 s in debug**, 10.4 times longer, and the tomato arms walk fourteen times as much
//! sequence again (8.0 Mb of BED against 572 kb).
//!
//! The alignments live outside this worktree — `crams/` is git-ignored, so a worktree does
//! not carry them — so the run needs `DEV_EXTRA_MOUNT` to reach them and `env` to carry the
//! three variables inside, since `scripts/dev.sh` forwards `HOME` and the target directory
//! and nothing else.
//!
//! ```text
//! BENCH=/path/to/pop_var_caller/benchmarks
//!
//! # HG002 — 100 spans, 572 kb of confident regions, one library, 30x
//! DEV_EXTRA_MOUNT=$BENCH ./scripts/dev.sh env \
//!   PVC_PREPASS_FASTA=$HOME/genomes/h_sapiens/gca_grch38/GCA_….fna \
//!   PVC_PREPASS_READS=$BENCH/giab/per_sample/bam/30x/HG002.30x.seed42.bam \
//!   PVC_PREPASS_BED=$BENCH/giab/per_sample/bed/HG002_bench_azar_merged_100.bed \
//!   cargo test --release --lib parameter_estimation::generic::real_alignments \
//!     -- --ignored --nocapture --test-threads=1
//!
//! # the same 100 spans undownsampled — the one run that reaches the depth cap
//! …  PVC_PREPASS_READS=$BENCH/giab/per_sample/bam/300x/HG002_reads_selected_100_rg.cram
//!
//! # tomato — 80 spans of 100 kb, one library
//! DEV_EXTRA_MOUNT=$BENCH ./scripts/dev.sh env \
//!   PVC_PREPASS_FASTA=$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
//!   PVC_PREPASS_READS=$BENCH/tomato1/crams/SRR7279481.p1.bench.cram \
//!   PVC_PREPASS_BED=$BENCH/tomato1/regions.bed \
//!   cargo test --release --lib parameter_estimation::generic::real_alignments \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! **The 300x arm is the one worth running when only one is.** It is where the subsampling
//! of Milestone C2 fires on real reads at all: 545,863 of its 550,049 sites arrive above the
//! ladder's cap of 124 and are drawn down to it, against zero on the 30x arm and 9,273 of
//! 7,213,401 on tomato SRR7279483 — about one site in 780. On any shallower run every
//! identity below is asserted over a table no site was ever subsampled into, so the deep
//! bins of the ladder, and `CountedSite::capped`'s draw, are exercised by nothing.
//!
//! **It does not test that sharding preserves the draw, and that is now true by
//! construction rather than by test.** The draw is seeded from the locus's own contig and
//! start (`seed_at`), and shards are whole typed regions, so no locus's coordinates differ
//! between the two arms and no draw can.
//!
//! **Identity 1, by contrast, is the same tautology at every depth**, and an earlier
//! version of this comment said the opposite. It claimed the two tables are filled above
//! the cap by *different draws* — one per read group, one shared for the whole site — which
//! describes a multi-library sample. The shared-draw function, `count_whole_site_by_library`,
//! is called only when the accumulator's `multi_library` is set, and identity 1 asserts
//! there is exactly one read group. At one read group `count_whole_site` and
//! `count_by_read_group` both end in `CountedSite::capped(depth, alt, cap, seed_at(region))`
//! with the same four arguments, and the draw is a pure function of them: one draw,
//! evaluated twice. `accumulators.rs`'s
//! `the_two_tables_agree_cell_for_cell_above_the_cap_too` already pins that on synthetic
//! depths of 300 to 963, and its doc had the mechanism right all along.
//!
//! `--test-threads=1` because the four tests each walk the same alignment file and the
//! useful output is the `--nocapture` lines, which interleave otherwise.
//!
//! A BED rather than a region string, and it is required: both cohorts ship one, the reads
//! exist only inside it, and typing a whole tomato chromosome to reach 100 kb of alignment
//! would spend the run's time in region typing.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use crate::fasta::ContigList;
use crate::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use crate::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusKind, SampleLocusObservationsIterator, UnhandledReason,
};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::generic::accumulators::{
    AccumulationCounts, ConstantPloidy, GenericAccumulators, InbreedingMode,
};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::generic::estimate::GenericEstimationConfig;
use crate::ng::parameter_estimation::generic::histogram::{Cell, DepthAltHistogram};
use crate::ng::parameter_estimation::generic::{GenericSampleParameters, error_rate_ladder};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::SampleReads;
use crate::ng::read::input::read_groups::{ReadGroups, SampleReadGroups, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::left_align::LeftAlignPreparer;
use crate::ng::ref_seq::WindowedRefSeq;
use crate::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, VerificationHandle,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::{
    GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig, TypedRegionError,
    TypedRegionIterator,
};
use crate::ng::types::{GenomeRegion, InbreedingF, Ploidy, Position};

use noodles_fasta::fai;

/// How many shards the typed regions are dealt to — **a maximum**, not a promise.
/// [`deal_into_shards`] gives fewer when there are fewer regions than this.
///
/// Sixteen rather than four: a shard boundary is the only thing the sharded arm has that
/// the unsharded one does not, so the arm is worth as many of them as it can have cheaply.
/// The work is the same either way — the same regions walked once — plus fifteen more file
/// opens.
const SHARD_COUNT: usize = 16;

/// What a walk produced, counted **outside** the accumulator.
///
/// It duplicates three numbers the accumulator also keeps, and the duplication is the
/// point: a table can only be checked against the loci that filled it by counting those
/// loci somewhere else. `loci` and `positions` catch a table that dropped or double-entered
/// sites, and `overlapping` puts a second implementation behind the partition claim.
///
/// **What it does not do is watch the accumulator's detector.** With no overlap anywhere in
/// the data, a working comparison and a deleted one both answer zero — see
/// [`no_locus_overlaps_the_one_before_it`].
#[derive(Copy, Clone, Default, PartialEq, Eq, Debug)]
struct WalkTally {
    /// Generic loci the stream yielded. The others are ignored, as `add_locus` ignores
    /// them.
    loci: u64,
    /// Reference positions those loci covered — `Σ region.len()`, which exceeds `loci`
    /// wherever the pileup widened a record to an indel's reference span.
    positions: u64,
    /// Loci beginning at or before the furthest position already covered on their contig.
    overlapping: u64,
}

/// Everything the environment named, opened once.
///
/// Held together because a walk needs all of it and the sharding identity makes five:
/// one parsed `.fai`, one contig table and one `OpenReference` shared by every walk, so
/// that what the arms differ in is the region list and nothing else.
struct WalkInputs {
    fasta: PathBuf,
    contigs: Arc<ContigList>,
    index: Arc<fai::Index>,
    reference: OpenReference,
    /// The BED, resolved against this reference's contig table.
    spans: GenomeRegions,
    read_groups: ReadGroups,
    sample: SampleReadGroups,
    /// The background check that the `.fai` describes this FASTA. **Joined at the end of
    /// every test rather than dropped**: a stale index would mean every arm agreed about
    /// the wrong bases, which is a green run that proves nothing.
    verification: Option<VerificationHandle>,
    /// The alignment file and the BED, as one string — every assertion message opens with
    /// it, because a failure here is read by whoever ran the command and not by whoever
    /// wrote it.
    run_label: String,
}

impl WalkInputs {
    fn from_env() -> Self {
        let fasta = PathBuf::from(required_env_var("PVC_PREPASS_FASTA"));
        let reads = PathBuf::from(required_env_var("PVC_PREPASS_READS"));
        let bed = PathBuf::from(required_env_var("PVC_PREPASS_BED"));

        // The convenience path, for parity.rs's two reasons: it sets
        // `ReferenceInfo.fasta_path`, without which a CRAM cannot be opened at all, and
        // with a `.fai` already present it verifies in the *background* rather than
        // making a whole-genome pass before the first read is decoded.
        let cache = Arc::new(ReferenceInfoCache::new());
        let (info, verification) = read_reference_verifying_or_creating_fai(
            &cache,
            fasta.clone(),
            ReferenceCheck::VerifyAgainstIndex,
        )
        .expect("the reference is readable and has (or can derive) a .fai");
        let contigs = Arc::new(info.contig_list());
        let index =
            WindowedRefSeq::read_index(&fasta).expect("the .fai beside the reference reads");

        let bounds: Vec<_> = contigs
            .entries
            .iter()
            .map(|entry| crate::regions::ContigBounds {
                name: &entry.name,
                length: u32::try_from(entry.length).expect("a contig shorter than 4 Gb"),
            })
            .collect();
        let spans = GenomeRegions::from_bed_path(&bed, &bounds)
            .expect("the BED resolves against this reference's contigs");

        let read_groups = build_read_groups(std::slice::from_ref(&reads))
            .expect("the alignment file's header declares its read groups");
        let sample = match read_groups.read_groups_per_sample() {
            [only] => only.clone(),
            other => panic!(
                "{} holds {} samples; these identities are per sample",
                reads.display(),
                other.len()
            ),
        };

        let run_label = format!("{} over {}", reads.display(), bed.display());
        Self {
            fasta,
            contigs,
            index,
            reference: OpenReference::new(info),
            spans,
            read_groups,
            sample,
            verification,
            run_label,
        }
    }

    /// The config every accumulator in one test comes from.
    ///
    /// **One config per test, and every shard's accumulator built from it**, because
    /// `merge` proves the shards share their binning rule and their ploidy map by pointer
    /// identity rather than by value.
    fn config(&self, inbreeding: InbreedingMode) -> GenericEstimationConfig {
        GenericEstimationConfig {
            sample_name: self.sample.sample.to_string(),
            read_groups: self.sample.read_groups.clone(),
            ploidy: Arc::new(ConstantPloidy(
                Ploidy::try_new(2).expect("a positive copy number"),
            )),
            inbreeding,
            fallback_error_rates: BTreeMap::new(),
            edges: Arc::new(DepthBinEdges::new()),
            read_admission: ReadFilterConfig::default(),
        }
    }

    /// The typed-region catalog over the BED, walked **once** and kept.
    ///
    /// Both arms of the sharding identity are fed from this one list, so the only thing
    /// they differ in is which accumulator each region's loci reach.
    /// Re-typing per arm would put region typing inside the comparison, and a difference
    /// there would read as a difference in the accumulator.
    fn typed_regions(&self) -> Vec<TypedRegion> {
        let mut typed = Vec::new();
        let walk = TypedRegionIterator::over_regions(
            WindowedRefSeq::with_shared_index(
                self.fasta.clone(),
                self.contigs.clone(),
                self.index.clone(),
            ),
            self.spans.clone(),
            TypedRegionConfig::default(),
        )
        .expect("the BED's spans name contigs this reference has");
        for item in walk {
            typed.push(item.expect("the reference reads through the whole typed-region walk"));
        }
        let generic = typed
            .iter()
            .filter(|region| region.kind == RegionKind::Generic)
            .count();
        assert!(
            generic > 0,
            "{}: the catalog typed no generic region, so a generic-path identity would \
             be asserted over an empty table",
            self.run_label
        );
        typed
    }

    /// Walk `regions`, add every locus to `into`, and tally what went past **outside** the
    /// accumulator.
    ///
    /// A fresh `SampleReads` and a fresh generator per call, which is what a region shard
    /// of a real run has: the accumulators are what merge, not the readers.
    fn accumulate(&self, regions: Vec<TypedRegion>, into: &mut GenericAccumulators) -> WalkTally {
        let sample = SampleReads::open(
            &self.sample,
            &self.read_groups,
            &self.reference,
            ReadFilterConfig::default(),
            true,
        )
        .expect("the alignment file opens against this reference");

        // ng's real read preparation — left-alignment over the canonical reference. The
        // preparer holds its own accessor, separate from the ones the generator mints.
        let preparer = LeftAlignPreparer::with_default_normalizer(self.accessor());
        let generator = {
            let fasta = self.fasta.clone();
            let contigs = self.contigs.clone();
            let index = self.index.clone();
            #[allow(
                clippy::arc_with_non_send_sync,
                reason = "PileupGenerator::new is generic over the accessor and takes Arc; this one is file-backed and single-threaded, as in examples/ng_generic_loci_dump.rs"
            )]
            let reference = Arc::new(self.accessor());
            PileupGenerator::new(
                reference,
                move || {
                    WindowedRefSeq::with_shared_index(fasta.clone(), contigs.clone(), index.clone())
                },
                preparer,
                PileupGeneratorConfig::default(),
            )
            .expect("the generic generator builds against this reference")
        };
        let generators = GeneratorSet::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Generator(Box::new(generator)),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );

        let regions: Vec<Result<TypedRegion, TypedRegionError>> =
            regions.into_iter().map(Ok).collect();
        let mut stream =
            SampleLocusObservationsIterator::new(regions.into_iter(), sample, generators);
        let mut tally = WalkTally::default();
        let mut previous_end: BTreeMap<_, u64> = BTreeMap::new();
        for locus in &mut stream {
            let locus = locus.expect("the walk runs to completion on a well-formed alignment");
            into.add_locus(&locus);
            if locus.kind != LocusKind::Generic {
                continue;
            }
            tally.loci += 1;
            tally.positions += locus.region.len();
            // `max`, matching the accumulator's own rule: a locus wholly inside the span
            // already covered must still count as overlapping, so the watermark may not
            // move backwards.
            let end = locus.region.end.get();
            match previous_end.entry(locus.region.contig) {
                std::collections::btree_map::Entry::Occupied(mut seen) => {
                    if locus.region.start.get() <= *seen.get() {
                        tally.overlapping += 1;
                    }
                    let watermark = seen.get_mut();
                    *watermark = (*watermark).max(end);
                }
                std::collections::btree_map::Entry::Vacant(empty) => {
                    empty.insert(end);
                }
            }
        }
        tally
    }

    fn accessor(&self) -> WindowedRefSeq {
        WindowedRefSeq::with_shared_index(
            self.fasta.clone(),
            self.contigs.clone(),
            self.index.clone(),
        )
    }

    /// Join the background check that the `.fai` describes this FASTA.
    ///
    /// **Called before a test's assertions, not after.** `VerificationHandle`'s `Drop` is
    /// deliberately silent while unwinding, so joining at the end means a failing test
    /// abandons the check with no output — and *the index does not match the FASTA, so
    /// every arm read the wrong bases* is the first hypothesis a reader of a red identity
    /// needs ruled out. The verification runs on a background thread and a walk over half a
    /// million loci has long outlasted it, so waiting costs nothing.
    fn confirm_reference(&mut self) {
        if let Some(handle) = self.verification.take() {
            handle
                .join()
                .expect("the .fai beside the reference describes it");
        }
    }
}

fn required_env_var(variable: &str) -> String {
    std::env::var(variable)
        .unwrap_or_else(|_| panic!("set {variable} — see this module's doc comment"))
}

/// Deal `regions` into at most [`SHARD_COUNT`] groups of **contiguous blocks of whole typed
/// regions**, the way a region-sharded run divides a genome.
///
/// # Whole segments to each worker, never a segment cut
///
/// **Owner's rule, 2026-08-09** — *"Once we start parallelizing we will send whole segments
/// to each worker, never a segment shall be cut"* — **and it was settled by a measurement
/// this file produced.** An earlier version cut each generic region into thirds to
/// manufacture extra boundaries, and that lost seventeen positions on
/// `tomato1/crams/SRR7279481.p1.bench.cram`: a read's 91-base deletion spanned
/// `SL4.0ch01:32,931,518–32,931,608`, a cut fell 74 bases inside it, and the part of the
/// deletion's reference span past the cut was emitted by no region at all.
///
/// **The walk is not what needs repairing.** The genome is divided from the reference
/// alone, without ever seeing a read, so a read's deletion can cross any boundary the
/// division makes — that much is unavoidable, and where a deletion crosses into a repeat
/// tract, the generic path dropping those bases is *correct*, because they are the STR
/// path's. What is avoidable is a boundary we invent inside territory that is entirely the
/// generic path's own. Splitting between whole segments costs nothing and removes the case
/// by construction, so this test divides that way and so must the parallel driver.
///
/// Contiguous and not round-robin: a shard records only the first start and last end of
/// each contig it saw, so interleaved shards would each span nearly the whole contig and
/// `adjustments().shard_spans_overlapping` — a counter that must read zero — would count
/// every adjacent pair in the sorted span list as an overlap.
fn deal_into_shards(regions: Vec<TypedRegion>) -> Vec<Vec<TypedRegion>> {
    let per_shard = regions.len().div_ceil(SHARD_COUNT).max(1);
    regions
        .chunks(per_shard)
        .map(<[TypedRegion]>::to_vec)
        .collect()
}

/// The first cell on which two tables disagree, rendered — or, when they differ only in
/// length, that.
///
/// **A failing `assert_eq!` on two `Vec<Cell>` prints both in full**, and the measured runs
/// hold 155 to 565 occupied cells apiece. The reader is in a terminal, per this module's
/// invocation, and what they need is the one cell that moved.
///
/// Exact equality on `Cell` is sound despite its `f64` `mean_depth`, and only because that
/// field is a ratio of two integers: `CellTally::mean_depth` is `depth_sum / sites`, IEEE-754
/// division is correctly rounded, and both operands are far inside the integers an `f64`
/// represents exactly (the largest reachable depth sum is 3.1 × 10¹¹). Equal integers
/// therefore give bit-identical quotients on both sides. **Were a weighted or smoothed depth
/// ever introduced, every `assert_eq!` on cells in this file would silently become the wrong
/// comparison.**
fn first_difference(one: &[Cell], two: &[Cell]) -> String {
    one.iter()
        .zip(two)
        .find(|(left, right)| left != right)
        .map_or_else(
            || {
                format!(
                    "no cell differs, but the tables hold {} and {}",
                    one.len(),
                    two.len()
                )
            },
            |(left, right)| format!("first difference: {left:?} against {right:?}"),
        )
}

/// A one-line summary of a table, for an assertion message that has to be read in a
/// terminal rather than in a debugger.
fn one_line_summary(table: &DepthAltHistogram<u64>, ploidy: Ploidy) -> String {
    format!(
        "{} loci over {} positions, {} reads, {} occupied cells",
        table.total_loci(),
        table.total_covered_positions(),
        table.total_reads(),
        table.cells(ploidy).len()
    )
}

/// **Identity 1 — the read-group histogram is the windowed one folded over its windows**,
/// cell for cell, on a sample with one library
/// (`spec/parameter_prepass_generic.md` §12.6, `arch/…` §9).
///
/// The premise is asserted rather than assumed: if the alignment file turns out to carry
/// more than one read group the test fails saying so, because at two libraries the two
/// tables are *supposed* to differ — a site splits into one read-group entry per library
/// while entering the windowed table once, whole.
///
/// **What it cannot say.** Both tables reduce the same locus through the same counting
/// functions and bin it with the same shared ladder, so the equality is close to
/// guaranteed by the types and fails only if a field is transposed or an arm is filed
/// twice. It is a plumbing check and not evidence about any of the four numbers — spec
/// §1's own words.
#[test]
#[ignore = "needs a real BAM/CRAM, reference and BED; see the module doc comment"]
fn the_two_tables_agree_cell_for_cell() {
    let mut inputs = WalkInputs::from_env();
    inputs.confirm_reference();
    assert_eq!(
        inputs.sample.read_groups.len(),
        1,
        "{}: sample {} carries {} read groups, and this identity holds only at one — at \
         two the read-group table splits a site the windowed table enters whole",
        inputs.run_label,
        inputs.sample.sample,
        inputs.sample.read_groups.len()
    );
    let group = inputs.sample.read_groups[0];

    let config = inputs.config(InbreedingMode::Fitted);
    let mut accumulators = config.accumulators();
    let walked = inputs.accumulate(inputs.typed_regions(), &mut accumulators);

    let ploidies = accumulators.ploidies();
    assert!(
        !ploidies.is_empty(),
        "{}: the walk entered no site, so this equality would compare two empty tables",
        inputs.run_label
    );

    for ploidy in ploidies {
        let folded = accumulators.whole_sample_histogram(ploidy);
        let by_group = accumulators
            .read_group_histograms()
            .get(&(group, ploidy))
            .unwrap_or_else(|| {
                panic!(
                    "{}: the windowed table holds ploidy {ploidy} and the read-group \
                     table holds nothing for {group:?} there",
                    inputs.run_label
                )
            });

        let (folded_cells, group_cells) = (folded.cells(ploidy), by_group.cells(ploidy));
        assert!(
            folded_cells == group_cells,
            "{}: at ploidy {ploidy} the folded windowed table ({}) and the read-group \
             table ({}) differ — {}",
            inputs.run_label,
            one_line_summary(&folded, ploidy),
            one_line_summary(by_group, ploidy),
            first_difference(&folded_cells, &group_cells)
        );
        assert_eq!(
            folded.total_covered_positions(),
            by_group.total_covered_positions(),
            "{}: at ploidy {ploidy} the two tables disagree about how much reference \
             they covered",
            inputs.run_label
        );
        eprintln!(
            "identity 1: {} ploidy {ploidy} — {}",
            inputs.run_label,
            one_line_summary(&folded, ploidy)
        );
    }

    // **Both tables against the loci themselves**, so that two tables agreeing because one
    // of them was filled from the other — or because both dropped the same sites — is not
    // the same green as two tables agreeing about what the walk produced. It also gives the
    // widening a witness: `positions` exceeds `loci` only where the pileup widened a record
    // to an indel's reference span, and a table charging one position per locus would pass
    // the cell-for-cell comparison above and fail here.
    let entered: u64 = accumulators
        .read_group_histograms()
        .values()
        .map(DepthAltHistogram::total_loci)
        .sum();
    assert_eq!(
        (entered, accumulators.covered_positions()),
        (walked.loci, walked.positions),
        "{}: the walk yielded {walked:?} and the tables hold {entered} sites over {} \
         positions",
        inputs.run_label,
        accumulators.covered_positions()
    );
    assert!(
        walked.positions > walked.loci,
        "{}: every one of the {} loci covered exactly one position, so nothing here \
         exercises a locus widened to an indel's reference span and the covered-position \
         comparison above is the site count again",
        inputs.run_label,
        walked.loci
    );
}

/// **Identity 2 — the same sample walked in one set of regions and in many gives identical
/// tables** (`arch/parameter_prepass_generic.md` §9).
///
/// The unsharded arm walks the catalog's regions in one stream. The sharded arm deals the
/// same regions — **whole, never split** — into sixteen contiguous blocks, walks each in its
/// own stream with its own reader and its own generator, and merges the sixteen
/// accumulators. So it differs in both ways a real parallel run does: the generator is
/// restarted at fifteen seams, and the counts arrive in sixteen pieces that integer addition
/// has to put back together.
///
/// **Whole regions is the contract, not a convenience** (owner, 2026-08-09). An earlier
/// version cut each generic region into thirds and lost seventeen positions to a read's
/// 91-base deletion straddling one of the cuts; [`deal_into_shards`] carries the
/// measurement and the reasoning.
///
/// **What it cannot say.** It compares this walk against this walk, so a locus the generator
/// drops at *every* boundary alike is invisible to it — including a deletion that crosses a
/// boundary the *catalog* made, which no arrangement of shards can move. What it can see is
/// a locus generated differently when its region is reached by a fresh stream, and any count
/// that `merge` fails to carry.
#[test]
#[ignore = "needs a real BAM/CRAM, reference and BED; see the module doc comment"]
fn one_walk_and_sixteen_shards_of_whole_regions_agree() {
    let mut inputs = WalkInputs::from_env();
    inputs.confirm_reference();
    let config = inputs.config(InbreedingMode::Fitted);
    let regions = inputs.typed_regions();

    let mut whole = config.accumulators();
    let walked_once = inputs.accumulate(regions.clone(), &mut whole);
    // **The premise the other two identities have and this one lacked.** Every comparison
    // below is satisfied by two *empty* accumulators — two empty ploidy sets, two vacuous
    // per-ploidy loops, two empty key lists, two default counter sets. And `walked_once`
    // does not rescue it, because the tally is built from `locus.region` outside the
    // accumulator and stays non-zero when nothing is entered. Measured: gutting
    // `add_locus` to return immediately fails identities 1 and 3 and panics the end-to-end
    // test, and left this one green.
    assert!(
        walked_once.loci > 0 && whole.covered_positions() > 0,
        "{}: the walk entered no site, so the two arms would compare empty tables and this          identity would hold whatever the merge did",
        inputs.run_label
    );

    let shards = deal_into_shards(regions.clone());
    assert!(
        shards.len() > 1 && shards.iter().all(|shard| !shard.is_empty()),
        "{}: {} typed regions fell into {} shards, so this arm merged nothing",
        inputs.run_label,
        regions.len(),
        shards.len()
    );
    let shard_count = shards.len();
    let mut merged = config.accumulators();
    let mut walked_sharded = WalkTally::default();
    for shard in shards {
        let mut accumulators = config.accumulators();
        let shard_tally = inputs.accumulate(shard, &mut accumulators);
        walked_sharded.loci += shard_tally.loci;
        walked_sharded.positions += shard_tally.positions;
        walked_sharded.overlapping += shard_tally.overlapping;
        merged.merge(accumulators);
    }

    // **The loci themselves, before the accumulator sees them.** Two arms whose tables
    // agree because both dropped the same loci would pass every comparison below; this is
    // the one assertion that fails when a seam loses a locus rather than mis-files it.
    assert_eq!(
        walked_once, walked_sharded,
        "{}: one walk yielded {walked_once:?} and {shard_count} shards yielded \
         {walked_sharded:?} — the region boundaries changed which loci exist",
        inputs.run_label
    );

    assert_eq!(
        whole.ploidies(),
        merged.ploidies(),
        "{}: the two arms saw different ploidies",
        inputs.run_label
    );
    for ploidy in whole.ploidies() {
        let one = whole.whole_sample_histogram(ploidy);
        let many = merged.whole_sample_histogram(ploidy);
        let (one_cells, many_cells) = (one.cells(ploidy), many.cells(ploidy));
        assert!(
            one_cells == many_cells,
            "{}: at ploidy {ploidy} one walk gave {} and {shard_count} shards gave {} — {}",
            inputs.run_label,
            one_line_summary(&one, ploidy),
            one_line_summary(&many, ploidy),
            first_difference(&one_cells, &many_cells)
        );
        assert_eq!(
            one.total_covered_positions(),
            many.total_covered_positions(),
            "{}: at ploidy {ploidy} the two arms covered different amounts of reference",
            inputs.run_label
        );
    }

    assert_eq!(
        whole.read_group_histograms().keys().collect::<Vec<_>>(),
        merged.read_group_histograms().keys().collect::<Vec<_>>(),
        "{}: the two arms disagree about which (read group, ploidy) pairs exist",
        inputs.run_label
    );
    for (&(group, ploidy), one) in whole.read_group_histograms() {
        let many = &merged.read_group_histograms()[&(group, ploidy)];
        let (one_cells, many_cells) = (one.cells(ploidy), many.cells(ploidy));
        assert!(
            one_cells == many_cells,
            "{}: {group:?} at ploidy {ploidy} — one walk gave {} and {shard_count} \
             shards gave {} — {}",
            inputs.run_label,
            one_line_summary(one, ploidy),
            one_line_summary(many, ploidy),
            first_difference(&one_cells, &many_cells)
        );
    }

    assert_eq!(
        whole.adjustments(),
        merged.adjustments(),
        "{}: the two arms adjusted the loci differently",
        inputs.run_label
    );
    eprintln!(
        "identity 2: {} — {} typed regions as one walk and as {shard_count} shards of \
         whole regions, identical",
        inputs.run_label,
        regions.len()
    );
}

/// **Identity 3 — no generic locus begins before the one before it on its contig ended**
/// (`arch/parameter_prepass_generic.md` §9, §3).
///
/// The generic loci must partition the positions they cover: a locus may span several
/// reference bases, because the pileup widens a record to an indel's reference span, but
/// two may not cover one position or that site's evidence enters the windowed table twice.
/// **A non-zero count here is a bug report against locus generation**, not something this
/// unit repairs — which is why the accumulator counts rather than de-duplicates.
///
/// The shard-seam counter is asserted with it, since this walk is unsharded and its own
/// spans therefore cannot overlap; a non-zero value would mean the span bookkeeping is
/// reporting overlaps that are not there, and a counter that must be zero may not have
/// false positives.
///
/// **What it cannot say, and the second count does not rescue it.** The overlap is counted
/// twice — once by the accumulator and once here, over the loci themselves — so the claim
/// *these loci partition their positions* rests on two implementations rather than one.
/// **That is all the second count buys.** It does not catch a detector that never fires:
/// on data with no overlaps a deleted comparison and a working one both answer zero, and
/// both counts stay green together. Whether the detector fires when it should is C3's
/// overlapping-loci fixture's job, and that is a fixture precisely because no real cohort
/// supplies the condition.
#[test]
#[ignore = "needs a real BAM/CRAM, reference and BED; see the module doc comment"]
fn no_locus_overlaps_the_one_before_it() {
    let mut inputs = WalkInputs::from_env();
    inputs.confirm_reference();
    let config = inputs.config(InbreedingMode::Fitted);
    let mut accumulators = config.accumulators();
    let walked = inputs.accumulate(inputs.typed_regions(), &mut accumulators);

    let AccumulationCounts {
        loci_with_upstream_subsample,
        reads_without_observation,
        sites_subsampled_to_cap,
        loci_overlapping_previous,
        shard_spans_overlapping,
    } = accumulators.adjustments();

    assert!(
        accumulators.covered_positions() > 0,
        "{}: the walk covered no position, so a zero overlap count says nothing",
        inputs.run_label
    );
    // **The accumulator's answer beside one computed outside it** — two implementations of
    // the same predicate over the same loci, so the partition claim does not rest on the
    // accumulator agreeing with itself. Neither can see a detector that never fires; that
    // is C3's overlapping-loci fixture's job.
    assert_eq!(
        (loci_overlapping_previous, walked.overlapping),
        (0, 0),
        "{}: the accumulator counted {loci_overlapping_previous} generic loci beginning \
         before the previous locus on their contig ended, and an independent pass over the \
         same {} loci counted {}. A non-zero first number is a defect in locus generation \
         rather than in this accumulator; a zero first and non-zero second means the \
         accumulator's detector is silent",
        inputs.run_label,
        walked.loci,
        walked.overlapping
    );
    assert_eq!(
        shard_spans_overlapping, 0,
        "{}: one unsharded walk reported {shard_spans_overlapping} overlapping shard \
         spans, so the span bookkeeping has false positives",
        inputs.run_label
    );
    // Printed and not asserted: these three are properties of the data, and a run where
    // they are large is one whose fitted rate describes a different population of reads
    // than the caller expects (`spec/parameter_prepass.md` §2).
    eprintln!(
        "identity 3: {} — {} positions covered, 0 overlapping loci; adjustments: \
         {loci_with_upstream_subsample} loci already subsampled upstream, \
         {reads_without_observation} reads witnessed nothing, \
         {sites_subsampled_to_cap} sites subsampled to the ladder's cap",
        inputs.run_label,
        accumulators.covered_positions()
    );
}

/// **The whole generic path, end to end on real reads** — the first time any fit in this
/// module sees an alignment file.
///
/// `F` is supplied rather than fitted, because the runs model needs 3,000 windows and both
/// cohorts' BEDs are a few megabytes of scattered spans; what this exercises is the coupled
/// error-rate/frequency fit and the fallback ladder, over loci a real walk produced.
///
/// **What it asserts, and it is deliberately little**: that the fit returns, that it
/// converged, and that no read group's rate landed on either end of the error-rate ladder.
/// A railed rate is the one shape in which this estimator reports a confident wrong number
/// instead of failing, and Phred 10 to 50 is wide enough that a real library is nowhere
/// near either end. **The values themselves are Milestone G's** — bounded there against the
/// GIAB truth set and against a coverage sweep — because nothing here knows what this
/// sample's heterozygosity ought to be.
#[test]
#[ignore = "needs a real BAM/CRAM, reference and BED; see the module doc comment"]
fn the_generic_path_fits_a_real_sample_without_railing() {
    let mut inputs = WalkInputs::from_env();
    inputs.confirm_reference();
    let supplied = InbreedingF::try_new(0.0).expect("zero is a fraction");
    let config = inputs.config(InbreedingMode::Supplied(supplied));
    let mut accumulators = config.accumulators();
    inputs.accumulate(inputs.typed_regions(), &mut accumulators);

    // **§12.6 on the collapsed table, which is the arm a real run over either cohort's BED
    // must use.** A supplied `F` drops the window key and pours the sample into
    // `whole_sample`, the `u64` table `6112803b` had to add; identities 1–3 all run in
    // `Fitted` mode and never touch it. The equality is covered synthetically — two
    // accumulator tests chain the collapsed fold to the fitted fold to the read-group
    // table — but not on a table a real walk filled.
    if inputs.sample.read_groups.len() == 1 {
        let group = inputs.sample.read_groups[0];
        for ploidy in accumulators.ploidies() {
            let folded = accumulators.whole_sample_histogram(ploidy);
            let by_group = &accumulators.read_group_histograms()[&(group, ploidy)];
            let (folded_cells, group_cells) = (folded.cells(ploidy), by_group.cells(ploidy));
            assert!(
                folded_cells == group_cells,
                "{}: at ploidy {ploidy} §12.6 fails on the collapsed table a supplied F \
                 builds — {}",
                inputs.run_label,
                first_difference(&folded_cells, &group_cells)
            );
        }
    }

    // **Arch §9's other half of "no fit rails", on a table a real walk built.** A cell
    // scored below its own alternative count is what the per-*bin* mean produces, and it
    // costs 5.2 rungs and 29% of the homozygous-non-reference rate with nothing visible on
    // the outside (spec §12.10). It is asserted on synthetic tables in `histogram.rs` and,
    // until now, on nothing a walk produced.
    for ploidy in accumulators.ploidies() {
        for cell in accumulators.whole_sample_histogram(ploidy).cells(ploidy) {
            assert!(
                f64::from(cell.key.alt_reads()) <= cell.mean_depth,
                "{}: the cell {:?} is scored at depth {}, below its own {} alternative \
                 reads",
                inputs.run_label,
                cell.key,
                cell.mean_depth,
                cell.key.alt_reads()
            );
        }
    }

    let parameters = accumulators
        .estimate(&config)
        .unwrap_or_else(|error| panic!("{}: {error}", inputs.run_label));

    // **The ladder ascends in Phred and so *descends* in error rate**: rung 0 is Phred 10,
    // an error rate of 0.1, and the last rung is Phred 50, or 10⁻⁵. Naming the ends by
    // their rung rather than by their magnitude is how the first version of this check got
    // its comparison the wrong way round and called a perfectly ordinary Phred-26 fit
    // railed.
    let ladder = error_rate_ladder();
    let coarsest = ladder.first().expect("the ladder has rungs").get();
    let finest = ladder.last().expect("the ladder has rungs").get();
    assert!(
        !parameters.error_rate.is_empty(),
        "{}: the fit returned no error rate at all, so the loop below checks nothing",
        inputs.run_label
    );
    for (group, rate) in &parameters.error_rate {
        // **Rail-checking a number no scan chose says nothing.** `resolve_error_rates`
        // hands back `DEFAULT_ERROR_RATE` = 0.001 with `Provenance::Defaulted` for a read
        // group below `MIN_SITES_TO_FIT` with no sibling to lend, and 0.001 sits
        // comfortably inside the ladder — so a sample whose every group defaulted passes
        // the check below while nothing was fitted. Unreachable on a cohort of
        // single-library samples with half a million sites, and this file takes an
        // alignment path from the environment.
        assert_eq!(
            rate.provenance,
            Provenance::FittedHere,
            "{}: {group:?}'s rate came back {:?} rather than fitted from its own sites, so              the rail check below would be a statement about a scan that never ran",
            inputs.run_label,
            rate.provenance
        );
        let value = rate.value.get();
        assert!(
            value < coarsest && value > finest,
            "{}: {group:?}'s error rate came back at {value}, an end of the ladder \
             [{finest}, {coarsest}] — a railed fit, which is the one way this estimator \
             returns a confident wrong number rather than failing",
            inputs.run_label
        );
        eprintln!(
            "end to end: {} {group:?} — error rate {value:.3e} ({:?}, {} reads)",
            inputs.run_label, rate.provenance, rate.observations
        );
    }
    assert!(
        parameters.coupled_fit.converged,
        "{}: the coupled fit ran out after {} iterations",
        inputs.run_label, parameters.coupled_fit.iterations
    );
    for (ploidy, rates) in &parameters.rates {
        eprintln!(
            "end to end: {} ploidy {ploidy} — heterozygosity {}, hom-non-ref {:.6} \
             ({} sites)",
            inputs.run_label,
            rates
                .value
                .observed_heterozygosity()
                .map_or_else(|| "n/a".to_string(), |het| format!("{:.6}", het.get())),
            rates.value.homozygous_non_reference_rate().get(),
            rates.observations
        );
    }

    report_the_second_class_of_site(&parameters, &inputs.run_label, coarsest, finest);
}

/// **The second class of site, printed and checked for the three ways it can come back
/// meaningless** — the half of the end-to-end test that the noise-model milestone added.
///
/// Neither cohort has a truth set for it, so nothing here can say the share or the rate is
/// *right*; what it can say is that the fit did not answer with a shape that carries no
/// information. Three such shapes, and each is a way this fit is known to fail:
///
/// 1. **A share at 0 or 1.** Every site noisy, or none, is one class of site wearing two
///    names, and the fit is supposed to have declined instead — `fit_site_noise` returns
///    nothing when the share collapses or the likelihood gain misses its floor.
/// 2. **A noisy class no noisier than the clean one.** Handed a clean rate well above the
///    truth, `fit_site_noise` rails at the ladder's *finest* rung, because a class at 10⁻⁵ is
///    the cheapest way to absorb the all-reference sites a too-high rate cannot explain. That
///    returns a perfectly well-formed pair describing nothing, and only the ordering of the
///    two rates sees it.
/// 3. **A noisy rate on either end of the ladder**, which is the same railed-fit failure the
///    caller above checks for the clean rate — an argmax at an endpoint is where this
///    estimator returns a confident wrong number instead of failing
///    (`arch/parameter_prepass_generic.md` §9).
///
/// The clean rate is not carried out of the fit — a sample emits the share-weighted marginal
/// and the pair — so it is recovered here by inverting the marginal,
/// `ε_clean = (ε_marginal − w·ε_noisy) / (1 − w)`, which is exact and needs the share to be
/// below one, checked first.
fn report_the_second_class_of_site(
    parameters: &GenericSampleParameters,
    run_label: &str,
    coarsest: f64,
    finest: f64,
) {
    let Some(noise) = parameters.site_noise else {
        if parameters.site_noise_off_the_ladder {
            // **An outlier, and reported as one rather than failed on** (owner, 2026-08-10).
            // This sample asked for a class of site noisier than one base in ten, which is
            // outside the range the error-rate ladder covers, so the fit refused it and fell
            // back to one rate. Two of the five alignments this test is run on do it.
            eprintln!(
                "end to end: {run_label} — OUTLIER: the sample asked for a second class of \
                 site noisier than the ladder's coarsest rung, so it was refused and one rate \
                 fitted. That is a population of positions this model does not describe — a \
                 duplication the reference does not carry shows about half its reads \
                 disagreeing, five times that rung"
            );
        } else {
            eprintln!(
                "end to end: {run_label} — no second class of site: the fit declined one, so \
                 every emitted rate is a single population of sites"
            );
        }
        return;
    };

    // **Printed before anything is asserted, and the order is deliberate.** A run that fails
    // one of the three checks below is a run whose numbers someone is about to want, and an
    // assertion that fires first takes them with it: the two tomato samples whose class rails
    // reported no rate at all until this was turned round.
    let share = noise.noisy_fraction();
    let noisy = noise.noisy_error_rate().get();
    let mut clean_rates = Vec::new();
    for (group, rate) in &parameters.error_rate {
        let marginal = rate.value.get();
        let clean = share.mul_add(-noisy, marginal) / (1.0 - share);
        eprintln!(
            "end to end: {run_label} {group:?} — {:.4}% of sites noisy at {noisy:.3e}, the rest \
             clean at {clean:.3e}, emitted as {marginal:.3e}",
            100.0 * share
        );
        clean_rates.push((*group, clean));
    }

    assert!(
        share > 0.0 && share < 1.0,
        "{run_label}: the fit returned a second class of site holding {share} of them, which \
         is one class of site under two names — it was supposed to decline instead"
    );
    assert!(
        noisy < coarsest && noisy > finest,
        "{run_label}: the noisy class's rate came back at {noisy}, an end of the ladder \
         [{finest}, {coarsest}] — a railed fit, which is how this estimator returns a \
         confident wrong number rather than failing"
    );
    for (group, clean) in clean_rates {
        assert!(
            noisy > clean,
            "{run_label}: {group:?}'s noisy class disagrees with the reference at {noisy}, no \
             more often than its clean class at {clean} — the shape `fit_site_noise` returns \
             when it is handed a clean rate above the truth and absorbs the all-reference \
             sites instead"
        );
    }
}

/// [`deal_into_shards`], which needs no alignment and until now ran on none.
///
/// **Every caller above is `#[ignore]`d**, so `cargo test` compiled this function and
/// executed it never. What it decides is which regions each worker gets, and its two ways
/// of being wrong are both silent: a shard that drops regions makes the sharded arm compare
/// a smaller genome, and shards that interleave make
/// `adjustments().shard_spans_overlapping` — a counter that must read zero — report
/// overlaps that are not there. These run in the ordinary suite, in microseconds.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::ContigId;

    fn generic(start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(end),
            },
            kind: RegionKind::Generic,
        }
    }

    /// The shards hold every region exactly once, in order, in contiguous blocks — and
    /// there are never more than [`SHARD_COUNT`] of them.
    ///
    /// The counts straddle the three cases: fewer regions than shards, an exact multiple,
    /// and a remainder.
    #[test]
    fn shards_partition_the_regions_into_contiguous_blocks() {
        for count in [1usize, 2, 15, 16, 17, 32, 100] {
            let regions: Vec<_> = (0..count as u64)
                .map(|at| generic(1 + at * 10, 5 + at * 10))
                .collect();
            let shards = deal_into_shards(regions.clone());
            assert!(
                !shards.is_empty() && shards.len() <= SHARD_COUNT,
                "{count} regions gave {} shards",
                shards.len()
            );
            assert!(
                shards.iter().all(|shard| !shard.is_empty()),
                "{count} regions gave an empty shard"
            );
            assert_eq!(
                shards.concat(),
                regions,
                "{count} regions were not partitioned in order"
            );
        }
    }

    /// **No shard boundary falls inside a region** — the owner's rule of 2026-08-09, which
    /// this function exists to keep and which no assertion above could see, since every
    /// comparison in the sharded arm holds whether or not a region was split.
    #[test]
    fn no_region_is_ever_split_between_shards() {
        let regions: Vec<_> = (0..100u64)
            .map(|at| generic(1 + at * 1_000, 900 + at * 1_000))
            .collect();
        let shards = deal_into_shards(regions.clone());
        for shard in &shards {
            for region in shard {
                assert!(
                    regions.contains(region),
                    "a shard holds {region:?}, which is not one of the regions handed in — \
                     a region was split"
                );
            }
        }
        assert_eq!(
            shards.iter().map(Vec::len).sum::<usize>(),
            regions.len(),
            "the shards hold a different number of regions than were handed in"
        );
    }
}
