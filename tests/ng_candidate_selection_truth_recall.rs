//! **Does the support bar keep the alleles that are really there?** — candidate selection
//! scored against the GIAB trio's v4.2.1 truth set, on loci cut from real reads.
//!
//! The module's own unit tests build loci by hand and prove the rule does what it says. This one
//! has no hand-built locus in it: it runs [`select_generic`] over the merge's real tables at two
//! depths and asks what the bar cost in true alternative alleles. It is the step that would catch
//! a rule that is self-consistent and wrong on data.
//!
//! # The fixture, and what a "true allele" is in it
//!
//! `tests/data/candidate_truth/` holds one file per depth, cut from the 100 benchmark intervals
//! of `HG002_bench_azar_merged_100.bed` (572 kb) over HG002/3/4. One row per
//! `(locus, sample, allele)`: that sample's pooled reads on the allele and its compared reads at
//! the locus, which are the numerator and the denominator the rule is asked with. **Only loci
//! carrying at least one true alternative are kept**, which is what makes it 78 kB rather than a
//! few megabytes.
//!
//! **A true allele is a haplotype the truth genotypes admit, projected onto the locus's span.**
//! That is the whole of the fixture's correctness and it is not the obvious construction. A
//! homozygous record sits on both copies and a heterozygous one on a single copy, so a span with
//! no het gives one sequence and a span with one het gives two — the homozygous records alone,
//! and those plus the het. **Projecting a record on its own when another sits beside it asks for a
//! sequence the sample does not carry**, and doing exactly that is what made two correctly-called
//! alleles look lost to the bar for a day (`doc/devel/ng/spec/candidate_alleles.md` §3.3): at
//! `chr1:90667287-90667293` HG002 carries a homozygous 2-base deletion and a heterozygous
//! substitution, the caller keeps both of its haplotypes at 162 and 127 reads, and the
//! substitution-without-the-deletion sequence that a naive projection looks for has one read
//! because no read carries it. **Spans with two or more unphased heterozygous records are dropped
//! rather than guessed** — 12 of them at 300×, 11 at 30× — because enumerating both readings
//! manufactures sequences nobody carries, which is the same defect wearing a different hat.
//!
//! # What each test is for
//!
//! The recall identity alone **cannot fail on the share term**, and that is worth stating because
//! it is what the plan's own description of this step asked for. Both it and the oracle ask the
//! same predicate, so deleting the share from the caller changes both sides together. **The
//! ladders are what make this file bite:** each asserts that a stricter bar loses true alleles a
//! looser one keeps, so a term that stopped working shows up as a ladder that stopped separating.
//!
//! Both ladders are exercised at both depths, because which half of `max(floor, ceil(share ×
//! compared))` decides is a property of the locus rather than of the run. **The label "300×" is a
//! run average and not a per-locus fact**: across this fixture a sample's compared reads at a
//! locus run from 8 to 428, median 273. So the floor is *not* inert at 300× — raising it from 2
//! to 3 costs two true alleles there at every share setting — and an argument that it must be,
//! because a 5-in-100 share of 300 reads asks for 15, holds only at the median locus.

use std::collections::BTreeMap;
use std::path::Path;

use pop_var_caller::ng::calling::allele_candidates::generic::select_generic;
use pop_var_caller::ng::calling::allele_candidates::{
    CandidateSelectionConfig, MaxCandidateAlleles, SelectionScratch,
};
use pop_var_caller::ng::run::cohort_merge::build::{
    AlleleSupport, CohortObservation, SampleSupport, SupportedAllele,
};
use pop_var_caller::ng::run::cohort_merge::{MinAltObs, MinAltReadShare, MinAltReads};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};

/// One sample's evidence at one fixture locus: its reads on each allele of the merge's table,
/// and its compared reads at the locus.
struct FixtureSample {
    sample: usize,
    reads: Vec<u32>,
    compared_reads: u32,
}

/// One locus of the fixture — the merge's table, which of its alternatives are true, and the
/// covering samples.
struct FixtureLocus {
    region: GenomeRegion,
    alleles: Vec<Vec<u8>>,
    /// Parallel to `alleles`; never true at index 0, which is the reference.
    is_truth: Vec<bool>,
    samples: Vec<FixtureSample>,
}

/// A cap wide enough that it cannot bind, so every one of these tests asks about the bar alone.
fn no_cap() -> MaxCandidateAlleles {
    MaxCandidateAlleles::new(u16::MAX).expect("the widest cap the type holds")
}

fn bar_of(floor: u32, share: f64) -> MinAltReads {
    MinAltReads {
        floor: MinAltObs(std::num::NonZeroU32::new(floor).expect("a non-zero floor")),
        share: MinAltReadShare::new(share).expect("a fraction of one"),
    }
}

fn read_fixture(path: &Path) -> Vec<FixtureLocus> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("reading {}: {error}", path.display()));
    // locus id -> (start, end, allele -> (bases, is_truth), sample -> (allele -> reads, compared))
    #[allow(clippy::type_complexity)]
    let mut built: BTreeMap<
        usize,
        (
            u64,
            u64,
            BTreeMap<usize, (Vec<u8>, bool)>,
            BTreeMap<usize, (BTreeMap<usize, u32>, u32)>,
        ),
    > = BTreeMap::new();
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let field: Vec<&str> = line.split(',').collect();
        assert_eq!(
            field.len(),
            9,
            "a fixture row must hold nine fields, and this one holds {}: {line}",
            field.len()
        );
        let parse = |at: usize| -> u64 { field[at].parse().expect("a number") };
        let locus = parse(0) as usize;
        let allele = parse(3) as usize;
        let entry = built
            .entry(locus)
            .or_insert_with(|| (parse(1), parse(2), BTreeMap::new(), BTreeMap::new()));
        entry
            .2
            .insert(allele, (field[4].as_bytes().to_vec(), field[5] == "1"));
        let sample = entry
            .3
            .entry(parse(6) as usize)
            .or_insert_with(|| (BTreeMap::new(), parse(8) as u32));
        sample.0.insert(allele, parse(7) as u32);
    }

    built
        .into_values()
        .map(|(start, end, alleles, samples)| {
            let table_len = alleles.len();
            assert!(
                alleles.keys().copied().eq(0..table_len),
                "a locus's allele indices must be dense and start at the reference"
            );
            let (bases, is_truth): (Vec<Vec<u8>>, Vec<bool>) = alleles.into_values().unzip();
            FixtureLocus {
                region: GenomeRegion {
                    contig: ContigId(0),
                    start: Position(start),
                    end: Position(end),
                },
                alleles: bases,
                is_truth,
                samples: samples
                    .into_iter()
                    .map(|(sample, (reads, compared_reads))| FixtureSample {
                        sample,
                        reads: (0..table_len)
                            .map(|a| *reads.get(&a).unwrap_or(&0))
                            .collect(),
                        compared_reads,
                    })
                    .collect(),
            }
        })
        .collect()
}

/// The fixture locus as the merge would have handed it over.
///
/// **One read group per sample**, because the fixture pools them — the read-group axis is the
/// module's own unit tests' (`one_run_per_allele`), and a fixture that flattened it cannot say
/// anything about it either way.
fn observation_of(locus: &FixtureLocus) -> CohortObservation {
    CohortObservation {
        region: locus.region,
        alleles: locus
            .alleles
            .iter()
            .map(|bases| bases.clone().into_boxed_slice())
            .collect(),
        per_sample: locus
            .samples
            .iter()
            .map(|sample| SampleSupport {
                sample: sample.sample,
                supported: sample
                    .reads
                    .iter()
                    .enumerate()
                    .filter(|(_, reads)| **reads > 0)
                    .map(|(allele, reads)| SupportedAllele {
                        allele,
                        read_group: ReadGroupId(0),
                        support: AlleleSupport {
                            num_reads: *reads,
                            // Finite and negative, as a log probability is. Nothing here reads
                            // the mass: the leftover is C3's own tests', and this file is about
                            // which alleles survive.
                            q_sum: -f64::from(*reads),
                            ..AlleleSupport::default()
                        },
                    })
                    .collect(),
                // The three counts the bar's denominator deliberately excludes (spec §5.1), and
                // the partials it also excludes. Zero here because the fixture does not carry
                // them — so this file cannot notice a denominator that let one of them in, and
                // the module's `the_denominator_is_the_samples_compared_reads_and_nothing_else`
                // is what does.
                partials: Vec::new(),
                reads_without_observation: 0,
                reads_removed_as_evidence: 0,
                reads_composed_across_records: 0,
            })
            .collect(),
    }
}

/// **The admission rule, re-derived here rather than asked of the module.**
///
/// Written as the arithmetic `max(floor, ceil(share × compared))` and not through
/// `MinAltReads::reached_by`, deliberately: a test that asks the module its own question can only
/// ever agree with it. Spelled this way, a term that stopped working on one side shows up as a
/// disagreement rather than as two matching answers.
fn some_sample_reached_the_bar(
    locus: &FixtureLocus,
    allele: usize,
    floor: u32,
    share: f64,
) -> bool {
    locus.samples.iter().any(|sample| {
        let required = floor.max((share * f64::from(sample.compared_reads)).ceil() as u32);
        sample.reads[allele] >= required
    })
}

/// The alternatives selection kept, as indices into the fixture's table.
fn admitted_by_selection(locus: &FixtureLocus, floor: u32, share: f64) -> Vec<usize> {
    let config = CandidateSelectionConfig {
        min_allele_support: bar_of(floor, share),
        max_candidate_alleles: no_cap(),
    };
    let selection = select_generic(
        &observation_of(locus),
        &config,
        &mut SelectionScratch::new(),
    );
    (1..locus.alleles.len())
        .filter(|&allele| selection.remap().candidate_for(allele).is_some())
        .collect()
}

/// True alternatives selection dropped at this bar.
fn true_alleles_lost(loci: &[FixtureLocus], floor: u32, share: f64) -> usize {
    loci.iter()
        .map(|locus| {
            let kept = admitted_by_selection(locus, floor, share);
            (1..locus.alleles.len())
                .filter(|&allele| locus.is_truth[allele] && !kept.contains(&allele))
                .count()
        })
        .sum()
}

fn fixture(depth: &str) -> Vec<FixtureLocus> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/candidate_truth")
        .join(format!("giab_trio_{depth}.csv"));
    let loci = read_fixture(&path);
    assert!(
        loci.len() > 500,
        "the {depth} fixture should hold hundreds of loci and holds {}",
        loci.len()
    );
    loci
}

/// **What the bar admits is exactly what some single sample's reads reached it with** — spec §3.2
/// in one assertion, over every locus of both fixtures and at every bar the ladders below use.
///
/// **This cannot fail on the share term alone** and the file's header says why; it is here to pin
/// the *shape* of the rule — one sample suffices, the denominator is that sample's own compared
/// reads — which is what a cohort-total bar or a two-sample bar would break, and both were live
/// alternatives (spec §3.2).
#[test]
fn the_bar_admits_exactly_the_alleles_some_single_sample_reached_it_with() {
    for depth in ["30x", "300x"] {
        let loci = fixture(depth);
        for (floor, share) in [(2, 0.0), (2, 0.02), (2, 0.05), (2, 0.10), (3, 0.02)] {
            for locus in &loci {
                let kept = admitted_by_selection(locus, floor, share);
                let expected: Vec<usize> = (1..locus.alleles.len())
                    .filter(|&allele| some_sample_reached_the_bar(locus, allele, floor, share))
                    .collect();
                assert_eq!(
                    kept, expected,
                    "at {depth}, locus {}, bar {floor} reads or {share}: selection and the \
                     re-derived rule disagree about which alternatives survive",
                    locus.region
                );
            }
        }
    }
}

/// **The floor is the expensive knob, and it is expensive at both depths.**
///
/// Raising it from 2 to 3 loses true alleles that 2 keeps — 4 at 30× and 2 at 300×. **The 300×
/// half is the one worth having**, because the argument that the floor cannot bind at depth is
/// wrong: a "300×" run's loci carry between 8 and 428 compared reads in this fixture, and it is
/// the thin ones the floor decides.
#[test]
fn raising_the_floor_costs_true_alleles_at_both_depths() {
    for (depth, share) in [("30x", 0.02), ("300x", 0.05)] {
        let loci = fixture(depth);
        let at_two = true_alleles_lost(&loci, 2, share);
        let at_three = true_alleles_lost(&loci, 3, share);
        assert!(
            at_three > at_two,
            "at {depth} with a share of {share}, a floor of 3 must lose true alleles a floor of \
             2 keeps — it lost {at_three} against {at_two}"
        );
    }
}

/// **The share is the other half, and doubling it costs true alleles at both depths.**
///
/// 5 in 100 against the shipped 10 in 100: the loss goes 0 → 1 at 30× and 0 → 2 at 300×. Without
/// this the share term could be deleted and the recall identity above would still hold, because
/// it and its oracle ask the same predicate.
#[test]
fn doubling_the_share_costs_true_alleles_at_both_depths() {
    for depth in ["30x", "300x"] {
        let loci = fixture(depth);
        let at_five = true_alleles_lost(&loci, 2, 0.05);
        let at_ten = true_alleles_lost(&loci, 2, 0.10);
        assert!(
            at_ten > at_five,
            "at {depth}, a share of 10 in 100 must lose true alleles that 5 in 100 keeps — it \
             lost {at_ten} against {at_five}"
        );
    }
}

/// **The shipped bar keeps almost every true allele, and loses some** — the two guards that stop
/// every assertion above holding vacuously.
///
/// Without the second, a caller that admitted everything would satisfy the identity and both
/// ladders would be comparing zero with zero. Without the first, so would a caller that admitted
/// nothing. The counts are asserted as a proportion and a floor rather than as golden numbers, so
/// a fixture recut at another depth still means something.
#[test]
fn the_shipped_bar_keeps_nearly_every_true_allele_and_loses_at_least_one() {
    for depth in ["30x", "300x"] {
        let loci = fixture(depth);
        let bar = pop_var_caller::ng::calling::allele_candidates::DEFAULT_MIN_ALLELE_SUPPORT;
        let (floor, share) = (bar.floor.get(), bar.share.get());
        let truth: usize = loci
            .iter()
            .map(|locus| locus.is_truth.iter().filter(|is| **is).count())
            .sum();
        let lost = true_alleles_lost(&loci, floor, share);
        assert!(
            lost >= 1,
            "at {depth} the shipped bar must lose at least one true allele, or this fixture \
             cannot tell a bar from no bar at all"
        );
        assert!(
            truth >= 500 && lost * 100 < truth,
            "at {depth} the shipped bar must keep more than 99 in 100 of the {truth} true \
             alternatives, and it lost {lost}"
        );
    }
}
