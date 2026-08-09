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
//! 2. **Sharded accumulation is exact.** The same territory walked as the catalog's own
//!    regions, and again with every generic region cut into pieces and the pieces dealt
//!    to four accumulators, must give identical tables — which the integer merge makes an
//!    equality rather than a tolerance (`arch/parameter_prepass_generic.md` §9).
//! 3. **`loci_overlapping_previous` is zero.** The generic loci must partition the
//!    positions they cover, or a site enters the windowed table twice. A non-zero count
//!    here is a bug report against locus generation, not something this unit absorbs.
//! 4. **The whole path runs end to end and lands nowhere near its ladder's ends** — the
//!    only one of the four that runs a fit, and the first time any fit in this module sees
//!    a read.
//!
//! # What these cannot say
//!
//! **Nothing here is evidence about any of the four numbers.** Identities 1–3 compare
//! objects the walk already built against each other; identity 4 asserts only that the fit
//! returned and did not rail. The values themselves are anchored in Milestone G, against
//! the GIAB truth set and against a coverage sweep.
//!
//! **The attributed arm of the cell key is not exercised.** Every sample in both cohorts
//! carries one read group — checked, not assumed, by
//! [`the_two_tables_agree_cell_for_cell`] — so `multi_library` is false throughout and no
//! site is ever keyed by the library its alternative reads came from. A real multi-library
//! alignment is what would test that, and neither cohort holds one.
//!
//! **`F` is not fitted.** The runs model needs
//! [`MIN_WINDOWS_TO_FIT_INBREEDING`](super::runs::MIN_WINDOWS_TO_FIT_INBREEDING) = 3,000
//! windows, which is 300 Mb that must hold sites, and both cohorts' BEDs are a few
//! megabytes of scattered spans. Identity 4 therefore runs with `F` supplied.
//!
//! # Running them
//!
//! `#[ignore]`d and driven by environment, so one file serves both organisms — the
//! convention `locus_generation/pileup/parity.rs` set. **`--release` is not an
//! optimisation**: real paired-end data hits production's reachable `debug_assert!`
//! constantly, and a debug build is also some twenty times slower over a walk this size.
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
//! **The 300x arm is the one worth running when only one is.** It is the only alignment in
//! either cohort deep enough to reach the ladder's cap of 124 reads: 545,863 of its 550,049
//! sites are subsampled there, against zero on the 30x arm and zero on two of the three
//! tomato samples. Above the cap the two tables are filled by **different draws** — one per
//! read group against each group's own depth, one shared draw for the site as a whole
//! (C2, C3) — and identity 1's claim that they nonetheless coincide at a single read group
//! is untested by any shallower run.
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
use crate::ng::parameter_estimation::generic::accumulators::{
    AccumulationCounts, ConstantPloidy, GenericAccumulators, InbreedingMode,
};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::generic::error_rate_ladder;
use crate::ng::parameter_estimation::generic::estimate::GenericEstimationConfig;
use crate::ng::parameter_estimation::generic::histogram::{Cell, DepthAltHistogram};
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

/// The widest piece the sharded arm cuts a generic region into.
///
/// **A ceiling and not the width** — see [`in_pieces`], which never cuts a region into
/// fewer than two pieces. A fixed width alone is a knob that silently does nothing: at
/// 10,000 bases it left the HG002 walk's 3,142 typed regions as 3,142 pieces, because
/// region typing had already fragmented the BED's spans to a mean of 176 bases and not one
/// of them reached the width. The arm then tested the merge and nothing about region
/// boundaries, and said so nowhere.
const PIECE_BP: u64 = 1_000;

/// How many accumulators the pieces are dealt to.
const SHARDS: usize = 4;

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
struct Target {
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
    /// How the run names itself in an assertion message — the alignment file and the BED,
    /// since a failure is read by whoever ran the command and not by whoever wrote it.
    where_: String,
}

impl Target {
    fn from_env() -> Self {
        let fasta = PathBuf::from(required("PVC_PREPASS_FASTA"));
        let reads = PathBuf::from(required("PVC_PREPASS_READS"));
        let bed = PathBuf::from(required("PVC_PREPASS_BED"));

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
        // An unknown contig or a malformed line is `from_bed_path`'s to reject, so what is
        // left for this to catch is an empty file — which would walk nothing and pass every
        // identity below.
        assert!(
            !spans.is_empty(),
            "{} resolved to no span, so the walk would cover nothing",
            bed.display()
        );

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

        let where_ = format!("{} over {}", reads.display(), bed.display());
        Self {
            fasta,
            contigs,
            index,
            reference: OpenReference::new(info),
            spans,
            read_groups,
            sample,
            verification,
            where_,
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
    /// they differ in is how the regions are cut and which accumulator each locus reaches.
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
            self.where_
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

    /// Join the background `.fai` check. Called at the end of every test that walked.
    fn confirm_reference(self) {
        if let Some(handle) = self.verification {
            handle
                .join()
                .expect("the .fai beside the reference describes it");
        }
    }
}

fn required(variable: &str) -> String {
    std::env::var(variable)
        .unwrap_or_else(|_| panic!("set {variable} — see this module's doc comment"))
}

/// Cut every generic region **in at least two**, and into pieces no wider than
/// [`PIECE_BP`], leaving every other kind whole.
///
/// **At least two, because a fixed width is a knob that can silently do nothing.** Region
/// typing fragments a BED span into runs of generic sequence between repeat tracts, and on
/// HG002 those average 176 bases — so any width worth setting for a 100 kb tomato span
/// leaves every HG002 region whole and the arm compares a walk against itself. Halving
/// puts an interior seam in every region longer than one base whatever the organism, and
/// the width then caps how long a piece may be on top of that.
///
/// **The generic ones only**, because those are the regions this step's accumulator reads
/// and the only ones whose splitting the identity is about; an STR tract cut in half is a
/// different tract, which is a question for the STR path's own plan.
fn in_pieces(regions: &[TypedRegion]) -> Vec<TypedRegion> {
    let mut pieces = Vec::new();
    for region in regions {
        if region.kind != RegionKind::Generic {
            pieces.push(region.clone());
            continue;
        }
        let (contig, start, end) = (
            region.region.contig,
            region.region.start.get(),
            region.region.end.get(),
        );
        // `.max(1)` so a one-base region yields one piece rather than a zero-width one
        // that would not advance the loop.
        let width = PIECE_BP.min((end - start + 1).div_ceil(2)).max(1);
        let mut at = start;
        while at <= end {
            // Saturating, so a width wider than what is left gives the remainder rather
            // than a `piece_end` below `at`.
            let piece_end = at.saturating_add(width - 1).min(end);
            pieces.push(TypedRegion {
                region: GenomeRegion {
                    contig,
                    start: Position(at),
                    end: Position(piece_end),
                },
                kind: RegionKind::Generic,
            });
            at = piece_end + 1;
        }
    }
    pieces
}

/// Deal `regions` to [`SHARDS`] groups in **contiguous blocks**, the way a region-sharded
/// run divides a genome.
///
/// Contiguous and not round-robin: a shard records the stretch of each contig its loci
/// covered, and `adjustments().shard_spans_overlapping` — a counter that must read zero —
/// would report every interleaved shard as overlapping every other.
fn in_shards(regions: Vec<TypedRegion>) -> Vec<Vec<TypedRegion>> {
    let per_shard = regions.len().div_ceil(SHARDS).max(1);
    regions
        .chunks(per_shard)
        .map(<[TypedRegion]>::to_vec)
        .collect()
}

/// Every cell of a table, with the ploidy each is to be scored at.
fn cells_of(table: &DepthAltHistogram<u64>, ploidy: Ploidy) -> Vec<Cell> {
    table.cells(ploidy)
}

/// A one-line summary of a table, for an assertion message that has to be read in a
/// terminal rather than in a debugger.
fn describe(table: &DepthAltHistogram<u64>, ploidy: Ploidy) -> String {
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
    let target = Target::from_env();
    assert_eq!(
        target.sample.read_groups.len(),
        1,
        "{}: sample {} carries {} read groups, and this identity holds only at one — at \
         two the read-group table splits a site the windowed table enters whole",
        target.where_,
        target.sample.sample,
        target.sample.read_groups.len()
    );
    let group = target.sample.read_groups[0];

    let config = target.config(InbreedingMode::Fitted);
    let mut accumulators = config.accumulators();
    let walked = target.accumulate(target.typed_regions(), &mut accumulators);

    let ploidies = accumulators.ploidies();
    assert!(
        !ploidies.is_empty(),
        "{}: the walk entered no site, so this equality would compare two empty tables",
        target.where_
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
                    target.where_
                )
            });

        assert_eq!(
            cells_of(&folded, ploidy),
            cells_of(by_group, ploidy),
            "{}: at ploidy {ploidy} the folded windowed table ({}) and the read-group \
             table ({}) differ",
            target.where_,
            describe(&folded, ploidy),
            describe(by_group, ploidy)
        );
        assert_eq!(
            folded.total_covered_positions(),
            by_group.total_covered_positions(),
            "{}: at ploidy {ploidy} the two tables disagree about how much reference \
             they covered",
            target.where_
        );
        eprintln!(
            "identity 1: {} ploidy {ploidy} — {}",
            target.where_,
            describe(&folded, ploidy)
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
        target.where_,
        accumulators.covered_positions()
    );
    assert!(
        walked.positions > walked.loci,
        "{}: every one of the {} loci covered exactly one position, so nothing here \
         exercises a locus widened to an indel's reference span and the covered-position \
         comparison above is the site count again",
        target.where_,
        walked.loci
    );
    target.confirm_reference();
}

/// **Identity 2 — the same sample walked in one set of regions and in many gives identical
/// tables** (`arch/parameter_prepass_generic.md` §9).
///
/// The unsharded arm walks the catalog's regions as they came. The sharded arm cuts every
/// generic one into [`PIECE_BP`]-base pieces, deals the pieces to [`SHARDS`] accumulators
/// in contiguous blocks, and merges them — so it differs in **both** ways a real run can:
/// the generator sees different region boundaries, and the counts arrive in four pieces
/// that integer addition has to put back together.
///
/// **What it cannot say.** It compares this walk against this walk, so a locus the
/// generator drops at *every* boundary alike is invisible to it. What it can see is a
/// locus generated differently at a seam, and any count that `merge` fails to carry.
#[test]
#[ignore = "needs a real BAM/CRAM, reference and BED; see the module doc comment"]
fn one_walk_and_four_shards_give_identical_tables() {
    let target = Target::from_env();
    let config = target.config(InbreedingMode::Fitted);
    let regions = target.typed_regions();

    let mut whole = config.accumulators();
    let walked_once = target.accumulate(regions.clone(), &mut whole);

    let pieces = in_pieces(&regions);
    assert!(
        pieces.len() > regions.len(),
        "{}: cutting {} typed regions gave {} pieces — nothing was cut, so this arm \
         would compare a walk against the same region boundaries and test only the merge",
        target.where_,
        regions.len(),
        pieces.len()
    );
    let piece_count = pieces.len();
    let shards = in_shards(pieces);
    assert!(
        shards.len() > 1,
        "{}: the pieces fell into one shard, so this arm merged nothing",
        target.where_
    );
    let shard_count = shards.len();
    let mut merged = config.accumulators();
    let mut walked_sharded = WalkTally::default();
    for shard in shards {
        let mut accumulators = config.accumulators();
        let shard_tally = target.accumulate(shard, &mut accumulators);
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
        target.where_
    );

    assert_eq!(
        whole.ploidies(),
        merged.ploidies(),
        "{}: the two arms saw different ploidies",
        target.where_
    );
    for ploidy in whole.ploidies() {
        let one = whole.whole_sample_histogram(ploidy);
        let many = merged.whole_sample_histogram(ploidy);
        assert_eq!(
            cells_of(&one, ploidy),
            cells_of(&many, ploidy),
            "{}: at ploidy {ploidy} one walk gave {} and {shard_count} shards gave {}",
            target.where_,
            describe(&one, ploidy),
            describe(&many, ploidy)
        );
        assert_eq!(
            one.total_covered_positions(),
            many.total_covered_positions(),
            "{}: at ploidy {ploidy} the two arms covered different amounts of reference",
            target.where_
        );
    }

    assert_eq!(
        whole.read_group_histograms().keys().collect::<Vec<_>>(),
        merged.read_group_histograms().keys().collect::<Vec<_>>(),
        "{}: the two arms disagree about which (read group, ploidy) pairs exist",
        target.where_
    );
    for (&(group, ploidy), one) in whole.read_group_histograms() {
        let many = &merged.read_group_histograms()[&(group, ploidy)];
        assert_eq!(
            cells_of(one, ploidy),
            cells_of(many, ploidy),
            "{}: {group:?} at ploidy {ploidy} — one walk gave {} and {shard_count} \
             shards gave {}",
            target.where_,
            describe(one, ploidy),
            describe(many, ploidy)
        );
    }

    assert_eq!(
        whole.adjustments(),
        merged.adjustments(),
        "{}: the two arms adjusted the loci differently",
        target.where_
    );
    eprintln!(
        "identity 2: {} — {} regions as one walk, {piece_count} pieces over \
         {shard_count} shards, identical",
        target.where_,
        regions.len()
    );
    target.confirm_reference();
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
    let target = Target::from_env();
    let config = target.config(InbreedingMode::Fitted);
    let mut accumulators = config.accumulators();
    let walked = target.accumulate(target.typed_regions(), &mut accumulators);

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
        target.where_
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
        target.where_,
        walked.loci,
        walked.overlapping
    );
    assert_eq!(
        shard_spans_overlapping, 0,
        "{}: one unsharded walk reported {shard_spans_overlapping} overlapping shard \
         spans, so the span bookkeeping has false positives",
        target.where_
    );
    // Printed and not asserted: these three are properties of the data, and a run where
    // they are large is one whose fitted rate describes a different population of reads
    // than the caller expects (`spec/parameter_prepass.md` §2).
    eprintln!(
        "identity 3: {} — {} positions covered, 0 overlapping loci; adjustments: \
         {loci_with_upstream_subsample} loci already subsampled upstream, \
         {reads_without_observation} reads witnessed nothing, \
         {sites_subsampled_to_cap} sites subsampled to the ladder's cap",
        target.where_,
        accumulators.covered_positions()
    );
    target.confirm_reference();
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
    let target = Target::from_env();
    let supplied = InbreedingF::try_new(0.0).expect("zero is a fraction");
    let config = target.config(InbreedingMode::Supplied(supplied));
    let mut accumulators = config.accumulators();
    target.accumulate(target.typed_regions(), &mut accumulators);

    let parameters = accumulators
        .estimate(&config)
        .unwrap_or_else(|error| panic!("{}: {error}", target.where_));

    // **The ladder ascends in Phred and so *descends* in error rate**: rung 0 is Phred 10,
    // an error rate of 0.1, and the last rung is Phred 50, or 10⁻⁵. Naming the ends by
    // their rung rather than by their magnitude is how the first version of this check got
    // its comparison the wrong way round and called a perfectly ordinary Phred-26 fit
    // railed.
    let ladder = error_rate_ladder();
    let coarsest = ladder.first().expect("the ladder has rungs").get();
    let finest = ladder.last().expect("the ladder has rungs").get();
    for (group, rate) in &parameters.error_rate {
        let value = rate.value.get();
        assert!(
            value < coarsest && value > finest,
            "{}: {group:?}'s error rate came back at {value}, an end of the ladder \
             [{finest}, {coarsest}] — a railed fit, which is the one way this estimator \
             returns a confident wrong number rather than failing",
            target.where_
        );
        eprintln!(
            "end to end: {} {group:?} — error rate {value:.3e} ({:?}, {} reads)",
            target.where_, rate.provenance, rate.observations
        );
    }
    assert!(
        parameters.coupled_fit.converged,
        "{}: the coupled fit ran out after {} iterations",
        target.where_, parameters.coupled_fit.iterations
    );
    for (ploidy, rates) in &parameters.rates {
        eprintln!(
            "end to end: {} ploidy {ploidy} — heterozygosity {}, hom-non-ref {:.6} \
             ({} sites)",
            target.where_,
            rates
                .value
                .observed_heterozygosity()
                .map_or_else(|| "n/a".to_string(), |het| format!("{:.6}", het.get())),
            rates.value.homozygous_non_reference_rate().get(),
            rates.observations
        );
    }
    target.confirm_reference();
}
