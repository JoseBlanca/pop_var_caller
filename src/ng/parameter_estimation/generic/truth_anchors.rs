//! G1 — the model-free anchors: what the GIAB truth set says, counted rather than fitted.
//!
//! **Every recovery test in this module generates its data from the model it then fits**, so a
//! shared misspecification cancels and the test passes. Those tests catch gross bugs and cannot
//! catch bias, which is the failure this step exists to remove
//! (`arch/parameter_prepass_generic.md` §9). This file is the first evidence that does not come
//! from the model at all: three numbers counted off the benchmark and the reads, with no fit
//! anywhere in them.
//!
//! - **The error rate** is disagreeing *reads* over total reads at loci the benchmark calls
//!   homozygous reference. No model: a read either matched the reference or it did not. Reads
//!   and not bases — a read carrying a twelve-base deletion disagrees once, not twelve times.
//! - **Heterozygosity** is the loci the benchmark calls heterozygous over the loci walked.
//! - **The homozygous-non-reference rate** is the same over the loci where no copy is the
//!   reference's.
//!
//! # What these can and cannot say
//!
//! **They bound the fitted values rather than pinning them, and the direction is not
//! symmetric.** The confident regions are the *easy* regions — the truth set excludes what it
//! could not call — so the model-free error rate is measured where reads behave and is a
//! **floor**: a fitted rate coming out *below* it is an unambiguous bug. Heterozygosity from
//! the same regions is depleted of the hard sequence where real variation concentrates, so it
//! is a floor too, and a fitted value above it is expected rather than wrong. **Nothing here is
//! a tolerance**; what is asserted is the one comparison whose failure has only one
//! explanation.
//!
//! # Per locus, not per base — owner's call, 2026-08-10
//!
//! `arch` §9 words the denominator as *the confident regions' length*, which is bases. **The
//! anchor counts loci instead**, and a locus classifies as heterozygous when **any** of its
//! reference positions carries a heterozygous benchmark record. The reason is that the fitted
//! heterozygosity this bounds is a per-locus rate — an ng locus can span several reference
//! positions — so a per-base truth and a per-locus fit would differ in their unit before they
//! differed in anything interesting. It is also the counting rule
//! `research/noise_model_overdispersion_2026-08-10.md` used to reach 9.9666 × 10⁻⁴ over 551,843
//! loci, the figure this milestone's ratios already rest on, so the anchor and those ratios are
//! one measurement rather than two that happen to agree.
//!
//! The cost is that the anchor covers the territory the alignment covers — the 100 selected
//! spans — and not the whole confident BED. For a per-read and a per-site rate that is a
//! restriction and not a hole (`arch` §9); for `F` it would void the anchor entirely, which is
//! why `F` is G3's problem and not this file's.
//!
//! # The reads are ng's own
//!
//! The research note reached its model-free error rate through `samtools mpileup`, which then
//! had to be talked into applying the same read filters as the walk. It does not have to: the
//! walk reports, per locus, how many reads there were and how many disagreed with the
//! reference, so the count here comes from **exactly the reads the estimator is fitted on**,
//! under `ReadFilterConfig::default()`, with no second tool to keep in step.
//!
//! **Above the depth cap that is fewer reads than the alignment holds, and deliberately so.**
//! `count_whole_site` subsamples a locus down to the ladder's cap of 124 by a draw seeded from
//! the locus's own coordinates. On the 300x arm that is 67.9 M reads over 549,180
//! homozygous-reference loci — 123.6 apiece, so essentially every locus is drawn down, where
//! the alignment holds something like 165 M. The count is still the right one to compare
//! against: the fit sees the identical draw from the identical seed, so both numbers describe
//! the same reads. It is simply not the alignment's read count.
//!
//! # Running it
//!
//! `#[ignore]`d and driven by environment, as `real_alignments.rs` is, plus the truth set:
//!
//! ```text
//! BENCH=/path/to/pop_var_caller/benchmarks
//!
//! DEV_EXTRA_MOUNT=$BENCH ./scripts/dev.sh env \
//!   PVC_PREPASS_FASTA=$HOME/genomes/h_sapiens/gca_grch38/GCA_….fna \
//!   PVC_PREPASS_READS=$BENCH/giab/per_sample/bam/30x/HG002.30x.seed42.bam \
//!   PVC_PREPASS_BED=$BENCH/giab/per_sample/bed/HG002_bench_azar_merged_100.bed \
//!   PVC_TRUTH_VCF=$BENCH/giab/all_bench_regions/vcfs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
//!   cargo test --release --lib parameter_estimation::generic::truth_anchors \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! **The whole-genome truth set**, not `benchmarks/ssr_hg002/`: that one is the tandem-repeat
//! benchmark, every record inside a repeat tract, and region typing routes those tracts to the
//! STR path (`arch` §9).

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use noodles_core::{Position as CorePosition, Region};
use noodles_vcf::variant::record::samples::series::Value;
use noodles_vcf::{self as vcf, io::indexed_reader::Builder as IndexedVcfBuilder};

use crate::ng::locus_generation::LocusKind;
use crate::ng::parameter_estimation::generic::accumulators::InbreedingMode;
use crate::ng::parameter_estimation::generic::depth_and_alt_reads::count_whole_site;
use crate::ng::parameter_estimation::generic::real_alignments::{WalkInputs, required_env_var};
use crate::ng::region_typing::TypedRegion;
use crate::ng::types::InbreedingF;

/// What the benchmark says about one reference position.
///
/// **Two classes and not the VCF's genotypes**, because what the fit reports is a count of
/// *non-reference copies*: zero, one, or all of them. A `1/2` record therefore lands in
/// [`EveryCopyNonReference`](TruthClass::EveryCopyNonReference) beside `1/1` — it names two
/// different alternative alleles, and neither of them is the reference.
///
/// **There is no third variant, and the absence is deliberate**: a position the benchmark makes
/// no statement about is `None` rather than a class of its own, because inside the confident
/// regions *no record here* means *no variant here* and that is the denominator, not a case.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum TruthClass {
    /// A record with one reference copy and one alternative — `0/1`, or `0|1`.
    Heterozygous,
    /// A record with no reference copy — `1/1`, or `1/2`.
    EveryCopyNonReference,
}

/// The benchmark's calls over the walked spans, by reference position.
///
/// **A record occupies its whole reference allele**, not just its first base: a deletion whose
/// `REF` is twelve bases makes a statement about twelve positions. **Measured, this changes
/// nothing on either arm** — collapsing every record to its first base gives byte-identical
/// counts — because the pileup already widens a generic locus to an indel's reference span, so
/// the locus that overlaps the record's tail is the same locus that overlaps its head. It is
/// kept because that is a property of the locus generator rather than of the truth set, and a
/// truth set is not the place to depend on it.
///
/// Keyed by `(contig, position)` because the walk meets contigs in the BED's order and there is
/// no dense position space to index into.
struct TruthCalls {
    by_position: BTreeMap<(String, u64), TruthClass>,
    /// How many **distinct** benchmark records were classified — distinct because
    /// `typed_regions()` splits the BED's 100 spans into some 3,142 typed regions and a record
    /// straddling a boundary comes back from both queries. Counting the raw returns gave 947
    /// where `bcftools` over the same BED gives 909, and the printed line is what a reader of
    /// Milestone G quotes.
    ///
    /// It counts records the classifier *placed*: a `0/0` record, or one with no `GT`, is not
    /// a statement about a variant and does not appear here.
    ///
    /// Reported because a truth set that resolved to nothing would otherwise look exactly like
    /// a sample with no variants — the failure this whole file exists to make impossible.
    records: u64,
    /// Positions where two records disagree about the class. Kept because the resolution
    /// below is a choice rather than a fact, and a reader deserves to know how often it fired.
    contested: u64,
}

impl TruthCalls {
    /// Read every benchmark record inside `spans`.
    ///
    /// **A contested position resolves to heterozygous.** Two records can touch one position —
    /// a SNP inside another record's deleted span, say — and if one says heterozygous and the
    /// other says every copy is non-reference, the position has a reference copy according to
    /// one of them. Taking the heterozygous reading keeps the every-copy-non-reference count a
    /// floor, which is the direction this anchor bounds in.
    fn read(
        vcf_path: &std::path::Path,
        spans: &[TypedRegion],
        contig_name: &dyn Fn(u32) -> String,
    ) -> Self {
        let mut reader = IndexedVcfBuilder::default()
            .build_from_path(vcf_path)
            .expect("the truth VCF opens and its .tbi is beside it");
        let header = reader.read_header().expect("the truth VCF has a header");

        let mut by_position: BTreeMap<(String, u64), TruthClass> = BTreeMap::new();
        // Scratch, not state: a record straddling two typed regions comes back from both
        // queries, and the printed count is a count of records rather than of returns.
        let mut seen_records: BTreeSet<(String, u64)> = BTreeSet::new();
        let (mut records, mut contested) = (0u64, 0u64);

        // One query per typed region. The BED's spans are the alignment's own, so a contig
        // absent from the truth set would be a mismatched reference rather than an empty
        // answer — the query says so loudly, which is why the error is not swallowed.
        for span in spans {
            let contig = contig_name(span.region.contig.0);
            let start = span.region.start.get();
            let end = span.region.end.get();
            let region = Region::new(
                contig.clone(),
                CorePosition::try_from(start as usize).expect("a 1-based start")
                    ..=CorePosition::try_from(end as usize).expect("a 1-based end"),
            );
            let query = reader
                .query(&header, &region)
                .unwrap_or_else(|error| panic!("querying {contig}:{start}-{end} failed: {error}"));
            for record in query.records() {
                let record = record.expect("a well-formed truth record");
                let Some(class) = class_of(&record, &header) else {
                    continue;
                };
                let Some(Ok(position)) = record.variant_start() else {
                    continue;
                };
                let first = position.get() as u64;
                if seen_records.insert((contig.clone(), first)) {
                    records += 1;
                }
                let last = first + record.reference_bases().len() as u64 - 1;
                for base in first..=last {
                    match by_position.entry((contig.clone(), base)) {
                        std::collections::btree_map::Entry::Occupied(mut seen) => {
                            if *seen.get() != class {
                                contested += 1;
                                seen.insert(TruthClass::Heterozygous);
                            }
                        }
                        std::collections::btree_map::Entry::Vacant(empty) => {
                            empty.insert(class);
                        }
                    }
                }
            }
        }

        Self {
            by_position,
            records,
            contested,
        }
    }

    /// What the benchmark says about a locus spanning `start..=end` — the strongest statement
    /// any of its positions carries, with heterozygous outranking every-copy-non-reference for
    /// the reason [`TruthCalls::read`] gives.
    fn classify(&self, contig: &str, start: u64, end: u64) -> Option<TruthClass> {
        let mut found = None;
        for base in start..=end {
            match self.by_position.get(&(contig.to_string(), base)) {
                Some(TruthClass::Heterozygous) => return Some(TruthClass::Heterozygous),
                Some(TruthClass::EveryCopyNonReference) => {
                    found = Some(TruthClass::EveryCopyNonReference);
                }
                None => {}
            }
        }
        found
    }
}

/// One sample's genotype, reduced to the three classes the fit reports.
///
/// `None` for anything else — a missing call, a haploid record, a record with no `GT` column.
/// Those are positions the benchmark declines to classify, and a position it declines to
/// classify is not evidence for or against anything.
fn class_of(record: &vcf::Record, header: &vcf::Header) -> Option<TruthClass> {
    let samples = record.samples();
    let genotype_column = samples.keys().iter().position(|key| key == "GT")?;
    let sample = samples.iter().next()?;
    let Some(Some(Ok(Value::Genotype(genotype)))) = sample.get_index(header, genotype_column)
    else {
        return None;
    };
    // Phasing is discarded: `0|1` and `0/1` are the same statement about copies.
    let alleles: Vec<Option<usize>> = genotype
        .iter()
        .map(|allele| allele.map(|(position, _phasing)| position))
        .collect::<std::io::Result<_>>()
        .ok()?;
    match alleles.as_slice() {
        [Some(a), Some(b)] if *a == 0 && *b == 0 => None,
        [Some(a), Some(b)] if *a == 0 || *b == 0 => Some(TruthClass::Heterozygous),
        [Some(_), Some(_)] => Some(TruthClass::EveryCopyNonReference),
        _ => None,
    }
}

/// What one walk counted, split by what the benchmark says about each locus.
#[derive(Default, Debug)]
struct ModelFreeCounts {
    loci: u64,
    heterozygous_loci: u64,
    every_copy_non_reference_loci: u64,
    /// Reads over loci the benchmark calls homozygous reference, and how many of them
    /// disagreed with it. The error rate is the second over the first, and nothing else.
    reads_at_homozygous_reference: u64,
    disagreeing_reads_at_homozygous_reference: u64,
}

impl ModelFreeCounts {
    fn error_rate(&self) -> f64 {
        self.disagreeing_reads_at_homozygous_reference as f64
            / self.reads_at_homozygous_reference as f64
    }

    fn heterozygosity(&self) -> f64 {
        self.heterozygous_loci as f64 / self.loci as f64
    }

    fn every_copy_non_reference_rate(&self) -> f64 {
        self.every_copy_non_reference_loci as f64 / self.loci as f64
    }
}

/// **The anchors, counted and then compared against the fit.**
///
/// The one assertion is the error rate's floor. Heterozygosity and the
/// every-copy-non-reference rate are printed rather than asserted, because the confident
/// regions are depleted of the sequence where real variation concentrates and a fitted value
/// above the truth is expected — there is no failure to define, only a number a reader of
/// Milestone G needs in front of them.
#[test]
#[ignore = "needs a real alignment, its reference and BED, and the GIAB truth VCF; see the module doc"]
fn the_fit_is_bounded_by_what_the_benchmark_counts() {
    let mut inputs = WalkInputs::from_env();
    inputs.confirm_reference();
    let truth_vcf = PathBuf::from(required_env_var("PVC_TRUTH_VCF"));

    let regions = inputs.typed_regions();
    let truth = TruthCalls::read(&truth_vcf, &regions, &|contig| inputs.contig_name(contig));
    assert!(
        truth.records > 0,
        "the truth set resolved to no records over these spans, which would make every locus \
         look homozygous reference and every number below meaningless"
    );

    let supplied = InbreedingF::try_new(0.0).expect("zero is a fraction");
    let config = inputs.config(InbreedingMode::Supplied(supplied));
    let mut accumulators = config.accumulators();
    let edges = config.edges.clone();

    let mut counts = ModelFreeCounts::default();
    inputs.for_each_locus(regions, &mut |locus| {
        accumulators.add_locus(locus);
        if locus.kind != LocusKind::Generic {
            return;
        }
        counts.loci += 1;
        let contig = inputs.contig_name(locus.region.contig.0);
        match truth.classify(&contig, locus.region.start.get(), locus.region.end.get()) {
            Some(TruthClass::Heterozygous) => counts.heterozygous_loci += 1,
            Some(TruthClass::EveryCopyNonReference) => counts.every_copy_non_reference_loci += 1,
            None => {
                let counted = count_whole_site(locus, &edges).counts();
                counts.reads_at_homozygous_reference += u64::from(counted.depth());
                counts.disagreeing_reads_at_homozygous_reference += u64::from(counted.alt_reads());
            }
        }
    });

    assert!(
        counts.loci > 0,
        "the walk produced no generic loci, so nothing below is a measurement"
    );
    eprintln!(
        "G1: {} — truth: {} records over {} loci; {} heterozygous, {} every-copy-non-reference, \
         {} contested positions",
        inputs.run_label,
        truth.records,
        counts.loci,
        counts.heterozygous_loci,
        counts.every_copy_non_reference_loci,
        truth.contested
    );
    eprintln!(
        "G1: {} — model-free: error rate {:.4e} ({} disagreeing of {} reads at \
         homozygous-reference loci), heterozygosity {:.4e}, every-copy-non-reference {:.4e}",
        inputs.run_label,
        counts.error_rate(),
        counts.disagreeing_reads_at_homozygous_reference,
        counts.reads_at_homozygous_reference,
        counts.heterozygosity(),
        counts.every_copy_non_reference_rate()
    );

    // **The truth's two classes must be the right way round, and nothing else here checks
    // it.** Swapping them — reading `0/1` as every-copy-non-reference and `1/1` as
    // heterozygous — leaves the error rate byte-identical, because both classes leave the
    // homozygous-reference denominator alike; only the two printed ratios move, and a printed
    // number asserts nothing. What separates them is biology: HG002 is an outbred human, so
    // heterozygous loci outnumber those with no reference copy. Measured here, 550 against
    // 317. A sample where that failed would be inbred enough that this anchor is the wrong
    // tool, not a sample this assertion should wave through.
    assert!(
        counts.heterozygous_loci > counts.every_copy_non_reference_loci,
        "{}: the benchmark calls {} loci heterozygous and {} with no reference copy. For an \
         outbred sample that is the wrong way round, and the likeliest cause is the two truth \
         classes being exchanged — which every other number here survives unchanged",
        inputs.run_label,
        counts.heterozygous_loci,
        counts.every_copy_non_reference_loci
    );

    // **One library, and this file inherits no guard from its neighbour.** The model-free
    // count is whole-site, pooled over every library; the comparison below is against each
    // read group's own rate. On a two-library sample the cleaner library could legitimately
    // sit below a pooled count and fail for no defect at all. Both cohorts carry one read
    // group, so this is checked rather than assumed.
    assert_eq!(
        inputs.read_group_count(),
        1,
        "{}: the model-free count is pooled over the whole site, so comparing it against each \
         of {} read groups' rates is not a like-for-like test",
        inputs.run_label,
        inputs.read_group_count()
    );

    let parameters = accumulators
        .estimate(&config)
        .unwrap_or_else(|error| panic!("{}: {error}", inputs.run_label));

    // **The one comparison is made without a loop, and that is the point.** It used to sit
    // inside `for … in &parameters.error_rate`, and a loop is a comparison that can be made
    // zero times: an empty map — or anything that skips the body — left the whole anchor green
    // with three console lines out of four and no assertion run at all. Asserting the map is
    // non-empty does not fix that; not looping does. There is exactly one read group here,
    // checked above, so the rate is taken directly.
    let (group, rate) = parameters
        .error_rate
        .iter()
        .next()
        .expect("a sample with one read group has one rate");
    {
        eprintln!(
            "G1: {} {group:?} — fitted error rate {:.4e} against a model-free {:.4e}, {:+.1}%",
            inputs.run_label,
            rate.value.get(),
            counts.error_rate(),
            100.0 * (rate.value.get() / counts.error_rate() - 1.0)
        );
        // **The one assertion, and its failure has a single explanation.** The confident
        // regions are the easy regions, so a rate counted there is a floor for a rate fitted
        // over the same loci: the fit sees every read this count saw and cannot see fewer
        // disagreements than were there. A fitted rate below it is not a tolerance question.
        assert!(
            rate.value.get() >= counts.error_rate(),
            "{}: {group:?}'s fitted error rate {:e} is BELOW the model-free count {:e} over \
             the same loci — the confident regions are the easy ones, so this direction has \
             no benign reading (arch §9)",
            inputs.run_label,
            rate.value.get(),
            counts.error_rate()
        );
    }

    assert!(
        !parameters.rates.is_empty(),
        "{}: the fit returned no genotype frequencies, so both ratios below would be skipped",
        inputs.run_label
    );
    for (ploidy, rates) in &parameters.rates {
        let fitted_het = rates
            .value
            .observed_heterozygosity()
            .map_or(f64::NAN, |het| het.get());
        eprintln!(
            "G1: {} ploidy {ploidy} — fitted heterozygosity {fitted_het:.4e} against a truth \
             {:.4e} ({:.3}x); fitted every-copy-non-reference {:.4e} against {:.4e} ({:.3}x)",
            inputs.run_label,
            counts.heterozygosity(),
            fitted_het / counts.heterozygosity(),
            rates.value.homozygous_non_reference_rate().get(),
            counts.every_copy_non_reference_rate(),
            rates.value.homozygous_non_reference_rate().get()
                / counts.every_copy_non_reference_rate()
        );
    }
}
