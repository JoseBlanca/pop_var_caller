//! **Two questions about the window-coverage measurement, asked of a real store.** Does a
//! one-base record's head really say what its evidence says (plan step B2, spec
//! `window_coverage.md` §3.1)? And does a calling run's per-locus window equal the one a
//! straight walk of the same store computes (plan step C3, spec §10)?
//!
//! ```text
//! # the first question alone — plan step B2's command, unchanged
//! cargo run --release --example ng_window_coverage_probe -- <a store.psp> [more stores...]
//!
//! # both — the run is asked to write down what it read, then the walk checks it
//! NG_WINDOW_COVERAGE_FILE=windows.tsv pop_var_caller_exp call-from-psps --psp <a store.psp> ...
//! cargo run --release --example ng_window_coverage_probe -- \
//!     --reference <reference.fa> --windows-from-the-run windows.tsv \
//!     <a store.psp> [more stores...]
//! ```
//!
//! Without `--reference` only the first question is asked, which is what plan step B2 ran. With
//! it, the walk also recomputes **every position's window** — the whole-store recomputation of
//! spec §10 — and, given `--windows-from-the-run`, checks the run's windows against its own. The
//! run writes that file when
//! [`recorded_windows::PATH_VARIABLE`](pop_var_caller::ng::run::cohort_merge::recorded_windows::PATH_VARIABLE)
//! names one, and nothing else about the run changes.
//!
//! **The stores must be listed in the run's own sample order**, because that is what the file
//! keys its rows by; the mode-equivalence oracle writes the psps under `<out>/psps/<sample>.psp`
//! in the order the VCF's sample columns are in.
//!
//! # The first question
//!
//! The window-coverage rule takes a sample's depth at a position from one of two places. Where a
//! record covers **one base**, it takes the count the record's head already carries — the reads
//! whose whole sequence over that base was compared against the reference — and never decodes the
//! evidence behind it. Where a record covers more, it builds the evidence and reads the depth off
//! it position by position.
//!
//! The first rule is what makes a run over stored files cheap, and it is only correct if the two
//! numbers agree. The argument for them agreeing is that **at a one-base locus no read's evidence
//! can stop inside the locus**, so every observation covers the whole of it and the head's count
//! is the sum of the observations' own. That is an argument about how the walk emits records, not
//! a theorem — spec §6's seventh trap says so — so this walks a real store and checks it.
//!
//! **A counter-example is a partial witness at a one-base locus**, and finding one changes the
//! rule to "build every body" before the merge learns to call it. So the probe reports both: how
//! many one-base records disagreed, and how many carry a witness that stops short at all, which
//! is the only way the disagreement can arise.
//!
//! # The second question
//!
//! A calling run reads each sample's records through the merge's cache, which draws them in
//! *covers* — one stretch of ground at a time, ahead of the builders — and hands a builder the
//! windows that sample's accumulator has finalised so far. A window centred at `p` is finalised
//! once the sample's stream reaches `p + 250`, so **a cover that stops at its region's last base
//! leaves that region's last centres unfinalised**: the builder finds nothing there
//! and the sample reads as absent at exactly the loci nearest every region boundary. Plan step
//! C3 makes each cover draw half a window further, and nothing about that failure is visible in
//! a run's output — the VCF is unchanged either way — so it needs an oracle outside the run.
//!
//! This is that oracle. It walks the store **end to end, in one pass, with no covers and no
//! eviction**, feeding the same accumulator the same positions and depths, and compares what it
//! finalises against the windows the run recorded. The comparison is by **bit pattern**, because
//! an absent window is a pair of `NaN`s and `NaN != NaN`.
//!
//! **One difference is expected and is not a defect**: the walk calls `finish` and so finalises
//! the tail of a sample's last records, while a run does not — nothing in the merge calls
//! `finish` until plan step C5, so a centre in the last half-window of a sample's stream has no
//! later position to close it. Those are counted and reported separately from a genuine
//! disagreement.
//!
//! # What it prints
//!
//! One tab-separated block per store, then one for the total. Besides the equality itself, the
//! counts say **how much of a real store the cheap rule covers**, in records and in positions —
//! and the two are different questions, because a wide *generic* record reports depth at one
//! position while a wide *tract* reports it at every position of its span.
//!
//! # One assumption, and the day it expires
//!
//! This walk builds each body through the psp reader, which applies a record's chain-id changes
//! before decoding it; the calling run builds a kept body through its own source, **against an
//! empty live set** (`src/ng/run/psp_source.rs`'s `build`). The two agree only while the encoder
//! writes no chain ids, which is true today and stops being true at the psp path's Milestone E.
//! When it does, this measurement has to be retaken rather than cited.
//!
//! # One walk, several measurements
//!
//! [`walk`] decodes every record once and hands it to whatever [`RecordMeasurement`]s the caller
//! passed; it knows nothing about what they count. [`SingleBaseEquality`] is the first question's
//! and [`WindowRecomputation`] the second's; neither knows about the other, and a third would be
//! a third implementor and one more entry in `main`'s slice.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use pop_var_caller::fasta::ContigList;
use pop_var_caller::ng::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use pop_var_caller::ng::psp::{PspReader, RecordHead};
use pop_var_caller::ng::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::ng::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::ng::run::cohort_merge::observation_cache::LocusSummary;
use pop_var_caller::ng::run::cohort_merge::recorded_windows::read_a_row;
use pop_var_caller::ng::types::{ContigId, GenomePosition, GenomeRegion};
use pop_var_caller::ng::window_coverage::depth::{EvidenceForOneRecord, for_each_reported_depth};
use pop_var_caller::ng::window_coverage::{
    self, WindowCoverage, WindowCoverageAccumulator, WindowCoverageConfig,
};

#[cfg(test)]
use pop_var_caller::ng::locus_generation::{LocusLen, SequenceObservation, SsrDetail};
#[cfg(test)]
use pop_var_caller::ng::types::{Motif, Position, ReadGroupId, SummedLogError};

/// How many disagreeing records to print in full before printing only the count. Enough to see
/// whether counter-examples cluster in one region or are scattered.
const COUNTER_EXAMPLES_SHOWN: usize = 20;

/// Something one pass over a store's records counts.
///
/// The walk owns the decoding and the measurement owns the arithmetic and its own output, so a
/// second measurement over the same pass costs a second implementor and nothing else.
trait RecordMeasurement {
    /// One record, with its body built. `contigs` is the store's contig list, for a measurement
    /// that has to name a position a reader will go looking at.
    fn observe(&mut self, head: &RecordHead, record: &SampleLocusObservations, contigs: &[String]);

    /// Everything this measurement has to say, as `label`-prefixed tab-separated rows.
    fn print(&self, label: &str);
}

/// What one store's walk found.
#[derive(Default)]
struct WhatTheWalkFound {
    records: u64,
    /// Records covering exactly one reference base — the ones the cheap rule applies to.
    one_base_records: u64,
    /// Of those, the ones where the head's count differs from the evidence's depth.
    one_base_records_disagreeing: u64,
    /// Of those, the ones carrying at least one observation whose witness stops inside the
    /// locus — the only mechanism by which the two numbers can differ.
    one_base_records_with_a_partial_witness: u64,
    /// Of those, the ones that are a repeat tract rather than a SNP/indel site. Spec §3.1 says
    /// a single-base record is generic by construction; this is what checks it.
    one_base_tract_records: u64,
    /// Summed head counts over the one-base records — the numerator of this store's mean
    /// observation depth. **A figure quoted from this probe is a fact about a store at a given
    /// depth**, and the share of records that span one base moves with indel density and depth,
    /// so the depth has to travel beside it.
    summed_head_depth_over_one_base_records: u64,
    /// Records covering more than one base, by kind: these are the ones the rule builds.
    wide_generic_records: u64,
    wide_tract_records: u64,
    wide_bundle_records: u64,
    /// Summed reference bases over the records covering more than one base.
    bases_in_wide_records: u64,
    /// Of those bases, the ones on a **tract** record — **kept apart because the two kinds
    /// report differently.** A wide generic record reports depth at its first base alone, a
    /// tract at every position of its span, so the positions the build branch decides is
    /// `wide_generic_records + bases_in_wide_tract_records` and pooling the bases loses it.
    bases_in_wide_tract_records: u64,
    /// A wide record of a kind added to `LocusKind` since this probe was written. Zero today,
    /// and printed so that it stops being zero visibly rather than by being folded into one of
    /// the counts above.
    wide_records_of_a_kind_this_probe_does_not_know: u64,
}

impl WhatTheWalkFound {
    /// **Destructured, so that a field added here is a compile error rather than a wrong
    /// total.** The counts a step adds are exactly the ones its report then quotes.
    fn add(&mut self, other: &Self) {
        let Self {
            records,
            one_base_records,
            one_base_records_disagreeing,
            one_base_records_with_a_partial_witness,
            one_base_tract_records,
            summed_head_depth_over_one_base_records,
            wide_generic_records,
            wide_tract_records,
            wide_bundle_records,
            bases_in_wide_records,
            bases_in_wide_tract_records,
            wide_records_of_a_kind_this_probe_does_not_know,
        } = other;
        self.records += records;
        self.one_base_records += one_base_records;
        self.one_base_records_disagreeing += one_base_records_disagreeing;
        self.one_base_records_with_a_partial_witness += one_base_records_with_a_partial_witness;
        self.one_base_tract_records += one_base_tract_records;
        self.summed_head_depth_over_one_base_records += summed_head_depth_over_one_base_records;
        self.wide_generic_records += wide_generic_records;
        self.wide_tract_records += wide_tract_records;
        self.wide_bundle_records += wide_bundle_records;
        self.bases_in_wide_records += bases_in_wide_records;
        self.bases_in_wide_tract_records += bases_in_wide_tract_records;
        self.wide_records_of_a_kind_this_probe_does_not_know +=
            wide_records_of_a_kind_this_probe_does_not_know;
    }

    /// Destructured for the same reason [`add`](Self::add) is: a count nothing prints is a count
    /// nobody reads.
    fn print(&self, label: &str) {
        let Self {
            records,
            one_base_records,
            one_base_records_disagreeing,
            one_base_records_with_a_partial_witness,
            one_base_tract_records,
            summed_head_depth_over_one_base_records,
            wide_generic_records,
            wide_tract_records,
            wide_bundle_records,
            bases_in_wide_records,
            bases_in_wide_tract_records,
            wide_records_of_a_kind_this_probe_does_not_know,
        } = self;
        println!("{label}\trecords\t{records}");
        println!("{label}\tone-base-records\t{one_base_records}");
        println!(
            "{label}\tone-base-share\t{:.4}",
            *one_base_records as f64 / (*records).max(1) as f64
        );
        println!("{label}\tone-base-records-disagreeing\t{one_base_records_disagreeing}");
        println!(
            "{label}\tone-base-records-with-a-partial-witness\t\
             {one_base_records_with_a_partial_witness}"
        );
        println!("{label}\tone-base-tract-records\t{one_base_tract_records}");
        println!(
            "{label}\tmean-head-depth-at-a-one-base-record\t{:.2}",
            *summed_head_depth_over_one_base_records as f64 / (*one_base_records).max(1) as f64
        );
        println!("{label}\twide-generic-records\t{wide_generic_records}");
        println!("{label}\twide-tract-records\t{wide_tract_records}");
        println!("{label}\twide-bundle-records\t{wide_bundle_records}");
        println!("{label}\tbases-in-wide-records\t{bases_in_wide_records}");
        println!("{label}\tbases-in-wide-tract-records\t{bases_in_wide_tract_records}");
        // **Positions, not records** — the question the rule's own documentation asks. A wide
        // generic record reports at one position, a tract at every position of its span, and a
        // one-base record at its own; these three are what the two branches decide between.
        let positions_from_the_build_branch = wide_generic_records + bases_in_wide_tract_records;
        let positions_reported = one_base_records + positions_from_the_build_branch;
        println!("{label}\tpositions-reported\t{positions_reported}");
        println!("{label}\tpositions-from-the-build-branch\t{positions_from_the_build_branch}");
        println!(
            "{label}\tpositions-from-the-head-share\t{:.4}",
            *one_base_records as f64 / (positions_reported).max(1) as f64
        );
        println!(
            "{label}\twide-records-of-a-kind-this-probe-does-not-know\t\
             {wide_records_of_a_kind_this_probe_does_not_know}"
        );
    }
}

/// One record whose head and evidence disagree, in the terms a reader would go looking with.
struct Disagreement {
    contig: String,
    position: u64,
    head_says: u32,
    evidence_says: u32,
    carries_a_partial_witness: bool,
}

/// Spec §3.1's cheap rule, checked: at every one-base record, the head's count against the depth
/// the built evidence gives at the record's only position.
#[derive(Default)]
struct SingleBaseEquality {
    found: WhatTheWalkFound,
    /// **Capped at [`COUNTER_EXAMPLES_SHOWN`], because only that many are ever printed.** Were
    /// the rule to fail wholesale, keeping one of these per record would exhaust memory before
    /// the report the probe exists to print had been printed.
    disagreements: Vec<Disagreement>,
}

impl SingleBaseEquality {
    fn add(&mut self, other: &Self) {
        self.found.add(&other.found);
    }
}

impl RecordMeasurement for SingleBaseEquality {
    fn observe(&mut self, head: &RecordHead, record: &SampleLocusObservations, contigs: &[String]) {
        self.found.records += 1;
        if head.region.len() == 1 {
            self.found.one_base_records += 1;
            if !matches!(record.kind, LocusKind::Generic) {
                self.found.one_base_tract_records += 1;
            }
            let carries_a_partial_witness = carries_a_partial_witness(record);
            if carries_a_partial_witness {
                self.found.one_base_records_with_a_partial_witness += 1;
            }
            // The number the rule would have used had it built the body: the depth
            // `num_obs_along_locus` gives at the record's only position. **Destructured rather
            // than indexed**, so that the one entry a one-base locus has stops being one entry
            // loudly: `first().unwrap_or(0)` would compare a fabricated zero against a head that
            // also read zero, and count the record as checked and agreeing.
            let depth_along_locus = record.num_obs_along_locus();
            let [evidence_says] = depth_along_locus[..] else {
                panic!(
                    "a record covering one base has one depth, and this one has {}",
                    depth_along_locus.len()
                );
            };
            let head_says = head.reads_compared_with_reference;
            self.found.summed_head_depth_over_one_base_records += u64::from(head_says);
            if evidence_says != head_says {
                self.found.one_base_records_disagreeing += 1;
                if self.disagreements.len() < COUNTER_EXAMPLES_SHOWN {
                    self.disagreements.push(disagreement_at(
                        head,
                        contigs,
                        head_says,
                        evidence_says,
                        carries_a_partial_witness,
                    ));
                }
            }
        } else {
            self.found.bases_in_wide_records += head.region.len();
            match record.kind {
                LocusKind::Generic => self.found.wide_generic_records += 1,
                LocusKind::Ssr(_) => {
                    self.found.wide_tract_records += 1;
                    self.found.bases_in_wide_tract_records += head.region.len();
                }
                LocusKind::SsrBundle => {
                    self.found.wide_bundle_records += 1;
                    self.found.bases_in_wide_tract_records += head.region.len();
                }
                // REVIEW ON UPGRADE: `LocusKind` is `#[non_exhaustive]`, so a match outside the
                // crate needs this arm even though the three above are every variant there is
                // today. The rule itself is in-crate and matches without a wildcard, so a fourth
                // kind is a compile error where it decides which positions a record reports.
                _ => self.found.wide_records_of_a_kind_this_probe_does_not_know += 1,
            }
        }
    }

    fn print(&self, label: &str) {
        self.found.print(label);
        for disagreement in &self.disagreements {
            println!(
                "{label}\tcounter-example\t{}:{}\thead={}\tevidence={}\tpartial-witness={}",
                disagreement.contig,
                disagreement.position,
                disagreement.head_says,
                disagreement.evidence_says,
                disagreement.carries_a_partial_witness,
            );
        }
    }
}

/// Spec §3.3's window, recomputed for every position of one store by a straight walk — the
/// whole-store recomputation of spec §10, and the oracle plan step C3's look-ahead is proved
/// against.
///
/// **It differs from the run in how it reads the store and in nothing else.** The same depth
/// rule ([`for_each_reported_depth`]), the same accumulator with the same configuration, the
/// same reference. What it does not have is covers: it walks from the first record to the last
/// in one pass, so no centre is ever asked for before the stream has passed it, and `finish`
/// closes the tail. A run's cover can do neither, and the difference between the two is exactly
/// what this exists to measure.
struct WindowRecomputation {
    /// `None` once [`close`](Self::close) has finished the pass.
    accumulator: Option<WindowCoverageAccumulator>,
    /// Every window this walk finalised, in ascending centre order — the whole store's, since
    /// nothing here evicts.
    windows: Vec<(GenomePosition, WindowCoverage)>,
    /// The reference this sample's GC is read from — **the run's own source**, so that a
    /// disagreement is about the walk and never about which bases were counted.
    reference: WindowedRefSeq,
    /// The contigs the reference names, checked against the store's on the first record of each
    /// contig: a psp carries the reference's own contig list, so the ids index both alike, and
    /// a store written against another reference must say so rather than silently reading the
    /// wrong bases.
    reference_contigs: ContigList,
    /// The contig the last record sat on, tested before the list below — a walk crosses a contig
    /// once and stays there for millions of records.
    contig_last_checked: Option<ContigId>,
    /// Every contig id checked against `reference_contigs` so far.
    contigs_checked: Vec<ContigId>,
    /// One record's reference bases, refilled per record rather than reallocated.
    bases: Vec<u8>,
    /// One record's reported positions and their depths, cleared per record for the same reason.
    reported: Vec<(GenomePosition, u32)>,
    /// Covered positions fed to the accumulator.
    positions_observed: u64,
    /// What the comparison against a run's recorded windows found, when one was given.
    comparison: Option<ComparisonWithTheRun>,
}

/// The run's windows for one sample, and what comparing them against this walk's found.
#[derive(Default)]
struct ComparisonWithTheRun {
    /// The run's own answer at each locus it built, by the locus's first base. `None` where the
    /// run had no window there.
    run_said: HashMap<GenomePosition, Option<WindowCoverage>>,
    /// Loci where both sides have a window and the two are bit-identical.
    agreed: u64,
    /// Loci where both sides have a window and the two differ — **the failure this exists to
    /// catch**, and it must be zero.
    disagreed: u64,
    /// Loci where the run had no window and this walk has one. Not zero: a centre in the last
    /// half-window of a sample's records is closed by `finish` here and by nothing in a run.
    the_run_had_no_window: u64,
    /// Loci where neither side has a window — ground no sample covered near enough to close a
    /// centre on, and the same answer from both.
    neither_side_had_a_window: u64,
    /// The first few disagreements, in the terms a reader would go looking with.
    disagreements: Vec<String>,
    /// The first few loci the run had no window for. **Printed rather than only counted**,
    /// because where they sit is the whole of what distinguishes an expected tail from a
    /// look-ahead that is short: a tail is a handful at the end of a sample's records, and a
    /// short look-ahead is a scatter at every building region's boundary.
    where_the_run_had_no_window: Vec<String>,
}

impl WindowRecomputation {
    fn over(
        reference: &ReferenceFasta,
        run_windows: Option<HashMap<GenomePosition, Option<WindowCoverage>>>,
    ) -> Self {
        Self {
            accumulator: Some(WindowCoverageAccumulator::new(WindowCoverageConfig {
                window_bp: window_coverage::WINDOW_BP,
                gc_bins: window_coverage::GC_BINS,
                depth_bins: window_coverage::DEPTH_BINS,
                depth_scale_windows: window_coverage::DEPTH_SCALE_WINDOWS,
                depth_range_in_medians: window_coverage::DEPTH_RANGE_IN_MEDIANS,
                min_window_positions: window_coverage::MIN_WINDOW_POSITIONS,
            })),
            windows: Vec::new(),
            reference: reference.accessor(),
            reference_contigs: reference.contigs.clone(),
            contig_last_checked: None,
            contigs_checked: Vec::new(),
            bases: Vec::new(),
            reported: Vec::new(),
            positions_observed: 0,
            comparison: run_windows.map(|run_said| ComparisonWithTheRun {
                run_said,
                ..ComparisonWithTheRun::default()
            }),
        }
    }

    /// Finish the pass: close the tail, take the windows, and compare against the run's.
    ///
    /// **Separate from [`print`](RecordMeasurement::print) because `finish` consumes the
    /// accumulator**, and a measurement's `print` takes `&self`.
    fn close(&mut self) {
        let accumulator = self
            .accumulator
            .take()
            .expect("the pass is closed exactly once");
        let (tail, _histogram) = accumulator.finish();
        self.windows.extend(tail);
        // **Ascending, and asserted rather than assumed**, because the comparison below and
        // plan step D1's distribution both binary-search it.
        assert!(
            self.windows.windows(2).all(|pair| pair[0].0 < pair[1].0),
            "the walk's windows are not in ascending centre order",
        );
        self.compare_against_the_run();
    }

    /// What this walk says at `centre`, or `None` where it finalised no window there.
    fn window_at(&self, centre: GenomePosition) -> Option<WindowCoverage> {
        self.windows
            .binary_search_by_key(&centre, |(at, _)| *at)
            .ok()
            .map(|found| self.windows[found].1)
    }

    fn compare_against_the_run(&mut self) {
        let Some(mut comparison) = self.comparison.take() else {
            return;
        };
        // Sorted, so that the examples printed are the first few in genome order rather than the
        // first few in a hash map's order, which changes run to run. `GenomePosition`'s own
        // ordering *is* genome order.
        let mut loci: Vec<(GenomePosition, Option<WindowCoverage>)> = comparison
            .run_said
            .iter()
            .map(|(at, said)| (*at, *said))
            .collect();
        loci.sort_unstable_by_key(|(at, _)| *at);
        for (at, run_said) in loci {
            let walk_says = self.window_at(at);
            match (run_said, walk_says) {
                (None, None) => comparison.neither_side_had_a_window += 1,
                (None, Some(_)) => {
                    comparison.the_run_had_no_window += 1;
                    if comparison.where_the_run_had_no_window.len() < COUNTER_EXAMPLES_SHOWN {
                        comparison
                            .where_the_run_had_no_window
                            .push(one_locus_as_text(at));
                    }
                }
                // The run finalised a centre this walk did not: impossible, since the walk sees
                // every position the run does and closes strictly more. Counted as a
                // disagreement so that it cannot pass unremarked.
                (Some(run), None) => {
                    comparison.remember_a_disagreement(at, Some(run), None);
                }
                (Some(run), Some(walk)) => {
                    if the_same_window(run, walk) {
                        comparison.agreed += 1;
                    } else {
                        comparison.remember_a_disagreement(at, Some(run), Some(walk));
                    }
                }
            }
        }
        self.comparison = Some(comparison);
    }

    /// Fill [`bases`](Self::bases) with the reference over `region`, checking once per contig
    /// that the store and the reference name it the same thing.
    ///
    /// **It fills the buffer rather than returning a slice of it**, because a returned slice
    /// borrows the whole of `self` and the caller needs the accumulator at the same time. The
    /// merge's own reading of the same ground is shaped the same way, and for the same reason.
    fn read_the_bases_over(&mut self, region: GenomeRegion, store_contigs: &[String]) {
        // A walk crosses a contig once and then stays on it for millions of records, so the last
        // one is tested before the list.
        if self.contig_last_checked != Some(region.contig)
            && !self.contigs_checked.contains(&region.contig)
        {
            let index = region.contig.get() as usize;
            let in_the_reference = self
                .reference_contigs
                .entries
                .get(index)
                .unwrap_or_else(|| {
                    panic!(
                        "the store names contig {index}, which the reference does not have — \
                         the store was written against another reference"
                    )
                });
            let in_the_store = store_contigs.get(index).unwrap_or_else(|| {
                panic!("the store's header has no name for the contig its own record sits on")
            });
            assert_eq!(
                &in_the_reference.name, in_the_store,
                "contig {index} is {in_the_store} in the store and {} in the reference",
                in_the_reference.name,
            );
            self.contigs_checked.push(region.contig);
        }
        self.contig_last_checked = Some(region.contig);
        self.reference
            .fetch_into(
                region.contig,
                region.start.get(),
                region.len(),
                &mut self.bases,
            )
            .unwrap_or_else(|why| panic!("reading the reference over {region}: {why}"));
    }
}

impl ComparisonWithTheRun {
    fn remember_a_disagreement(
        &mut self,
        at: GenomePosition,
        run: Option<WindowCoverage>,
        walk: Option<WindowCoverage>,
    ) {
        self.disagreed += 1;
        if self.disagreements.len() < COUNTER_EXAMPLES_SHOWN {
            self.disagreements.push(format!(
                "{}\trun={}\twalk={}",
                one_locus_as_text(at),
                one_window_as_text(run),
                one_window_as_text(walk),
            ));
        }
    }
}

/// Whether the run's window and this walk's are the same measurement.
///
/// **By bit pattern, which is what `WindowCoverage`'s own `PartialEq` does.** An absent window is
/// a pair of `NaN`s, and `NaN != NaN` — so field-wise `==` on the floats would report a
/// disagreement at every centre the position floor silenced, which at three reads a position is a
/// large share of a store, and would refuse a correct run (spec §6 trap 4). Written as a named
/// function because `run == walk` on two `f32`-carrying structs reads as float equality and is
/// not.
fn the_same_window(run: WindowCoverage, walk: WindowCoverage) -> bool {
    run == walk
}

/// One locus as a reader would go looking for it. The contig is an id and not a name because a
/// row of this output is keyed by what the store and the reference both index by.
fn one_locus_as_text(at: GenomePosition) -> String {
    format!("contig{}:{}", at.contig.get(), at.position.get())
}

/// One window as a reader would go looking for it: both fields, or the word for what is missing.
fn one_window_as_text(window: Option<WindowCoverage>) -> String {
    match window {
        None => "none".to_owned(),
        Some(window) if window.is_absent() => "absent".to_owned(),
        Some(window) => format!("gc={} depth={}", window.gc_fraction, window.mean_depth),
    }
}

impl RecordMeasurement for WindowRecomputation {
    fn observe(&mut self, head: &RecordHead, record: &SampleLocusObservations, contigs: &[String]) {
        let region = head.region;
        // **The summary is built from the head and not derived from the record**, because that
        // is what the run's cache holds in psp mode: the head's count is the number the cheap
        // rule uses, and deriving it here would test the rule against itself.
        let summary = LocusSummary {
            region,
            non_reference_reads: head.non_reference_reads,
            reads_compared_with_reference: head.reads_compared_with_reference,
        };
        self.reported.clear();
        // **Collected before it is folded**, because the borrow the base lookup needs and the
        // borrow the accumulator needs cannot be held at once; the buffer is a field so that the
        // collection costs no allocation per record.
        let mut reported = std::mem::take(&mut self.reported);
        for_each_reported_depth::<std::convert::Infallible>(
            summary,
            EvidenceForOneRecord::InHand(record),
            |_| unreachable!("this walk builds every body before the rule sees it"),
            |at, depth| reported.push((at, depth)),
        )
        .expect("the rule cannot fail on evidence already in hand");
        if reported.is_empty() {
            self.reported = reported;
            return;
        }
        self.read_the_bases_over(region, contigs);
        let first = region.start.min(region.end).get();
        // **Every field named, none absorbed by a `..`**, so that a field this measurement gains
        // has to be answered for here rather than silently going unfed.
        let Self {
            accumulator,
            bases,
            positions_observed,
            windows,
            reference: _,
            reference_contigs: _,
            contig_last_checked: _,
            contigs_checked: _,
            reported: _,
            comparison: _,
        } = self;
        let accumulator = accumulator
            .as_mut()
            .expect("records are observed before the pass is closed");
        for &(at, depth) in reported.iter() {
            let offset = (at.position.get() - first) as usize;
            let base = *bases.get(offset).unwrap_or_else(|| {
                panic!("the record at {region} reports a depth at {at:?}, which is off its ground")
            });
            accumulator.observe(at.contig, at.position, base, depth);
            *positions_observed += 1;
        }
        while let Some(window) = accumulator.pop_ready() {
            windows.push(window);
        }
        self.reported = reported;
    }

    fn print(&self, label: &str) {
        assert!(
            self.accumulator.is_none(),
            "the pass must be closed before its windows are reported",
        );
        println!("{label}\tpositions-observed\t{}", self.positions_observed);
        println!("{label}\twindows-finalised\t{}", self.windows.len());
        let absent = self
            .windows
            .iter()
            .filter(|(_, window)| window.is_absent())
            .count();
        println!("{label}\twindows-under-the-floor\t{absent}");
        let Some(comparison) = &self.comparison else {
            return;
        };
        println!("{label}\tloci-the-run-built\t{}", comparison.run_said.len());
        println!("{label}\tloci-agreeing-bit-for-bit\t{}", comparison.agreed);
        println!("{label}\tloci-disagreeing\t{}", comparison.disagreed);
        println!(
            "{label}\tloci-the-run-had-no-window-for\t{}",
            comparison.the_run_had_no_window
        );
        println!(
            "{label}\tloci-neither-side-has-a-window-for\t{}",
            comparison.neither_side_had_a_window
        );
        for disagreement in &comparison.disagreements {
            println!("{label}\twindow-counter-example\t{disagreement}");
        }
        for at in &comparison.where_the_run_had_no_window {
            println!("{label}\tthe-run-had-no-window-at\t{at}");
        }
    }
}

/// The reference the walk reads its GC from, opened once and shared by every store's accessor.
struct ReferenceFasta {
    path: PathBuf,
    contigs: ContigList,
}

impl ReferenceFasta {
    fn open(path: &Path) -> Self {
        let cache = std::sync::Arc::new(ReferenceInfoCache::new());
        let (info, _verification) = read_reference_verifying_or_creating_fai(
            &cache,
            path.to_path_buf(),
            ReferenceCheck::VerifyAgainstIndex,
        )
        .unwrap_or_else(|why| panic!("opening {}: {why}", path.display()));
        Self {
            path: path.to_path_buf(),
            contigs: info.contig_list(),
        }
    }

    /// **One accessor per store, not one for all of them.** A `WindowedRefSeq` is a forward
    /// reader with its own cursor and window; sharing one between two walks would make them
    /// seek against each other.
    fn accessor(&self) -> WindowedRefSeq {
        WindowedRefSeq::new(self.path.clone(), self.contigs.clone())
    }
}

/// What the probe was asked to do.
struct Arguments {
    stores: Vec<PathBuf>,
    /// The reference, when the whole-store window recomputation was asked for.
    reference: Option<PathBuf>,
    /// The windows a calling run recorded, to check the recomputation against.
    windows_from_the_run: Option<PathBuf>,
}

const USAGE: &str = "usage: ng_window_coverage_probe [--reference <reference.fa>] \
                     [--windows-from-the-run <windows.tsv>] <a store.psp> [more stores...]";

impl Arguments {
    fn from_command_line() -> Self {
        let mut stores = Vec::new();
        let mut reference = None;
        let mut windows_from_the_run = None;
        let mut argv = std::env::args().skip(1);
        while let Some(argument) = argv.next() {
            match argument.as_str() {
                "--reference" => {
                    reference = Some(PathBuf::from(
                        argv.next().expect("--reference takes a path"),
                    ));
                }
                "--windows-from-the-run" => {
                    windows_from_the_run = Some(PathBuf::from(
                        argv.next().expect("--windows-from-the-run takes a path"),
                    ));
                }
                other => {
                    assert!(!other.starts_with("--"), "unknown option {other}\n{USAGE}");
                    stores.push(PathBuf::from(other));
                }
            }
        }
        assert!(!stores.is_empty(), "{USAGE}");
        assert!(
            windows_from_the_run.is_none() || reference.is_some(),
            "--windows-from-the-run compares against the recomputation, which needs --reference",
        );
        Self {
            stores,
            reference,
            windows_from_the_run,
        }
    }
}

/// The windows a calling run recorded, split by sample — one map per sample index, in the run's
/// own sample order, which is the order the stores are listed in.
///
/// **The row format is the library's**, written and read by one pair of functions in
/// [`recorded_windows`](pop_var_caller::ng::run::cohort_merge::recorded_windows), because a column
/// reordered on one side alone would still parse on the other and every locus would then read as
/// one the run had no window for — which is what a *pass* looks like in this measurement. What is
/// decided here is what the file means: which sample index is in range, and what two rows for one
/// locus are.
///
/// **A locus written twice should not happen.** A builder skips a locus it does not own rather
/// than building it, so every locus is built exactly once however the genome is divided
/// (`cohort_merge/build.rs`'s ownership rule). Two identical rows are tolerated; two that differ
/// mean either that rule broke or that two runs wrote into one file, and either way the number
/// this probe is about cannot be trusted, so it stops.
fn windows_the_run_recorded(
    recorded: &Path,
    samples: usize,
) -> Vec<HashMap<GenomePosition, Option<WindowCoverage>>> {
    let text = std::fs::read_to_string(recorded)
        .unwrap_or_else(|why| panic!("reading {}: {why}", recorded.display()));
    let mut per_sample = vec![HashMap::new(); samples];
    for (number, line) in text.lines().enumerate() {
        let row = read_a_row(line).unwrap_or_else(|why| {
            panic!("{}:{}: {why}", recorded.display(), number + 1);
        });
        assert!(
            row.sample < samples,
            "{}:{}: the file names sample {} and only {samples} stores were given — the stores \
             must be listed in the run's own sample order",
            recorded.display(),
            number + 1,
            row.sample,
        );
        if let Some(already) = per_sample[row.sample].insert(row.at, row.window) {
            assert_eq!(
                already.map(|w| (w.gc_fraction.to_bits(), w.mean_depth.to_bits())),
                row.window
                    .map(|w| (w.gc_fraction.to_bits(), w.mean_depth.to_bits())),
                "{}:{}: sample {} has two different windows at {:?} — a locus is built exactly \
                 once, so either the merge's ownership rule broke or two runs wrote into this \
                 file",
                recorded.display(),
                number + 1,
                row.sample,
                row.at,
            );
        }
    }
    per_sample
}

fn main() {
    let arguments = Arguments::from_command_line();
    let stores = &arguments.stores;
    let reference = arguments.reference.as_deref().map(ReferenceFasta::open);
    let mut run_windows = arguments
        .windows_from_the_run
        .as_deref()
        .map(|recorded| windows_the_run_recorded(recorded, stores.len()));

    println!("phase\tng-window-coverage-probe");
    let mut total = SingleBaseEquality::default();
    // **`None` where the recomputation did not run**, so its rows are absent rather than zero.
    let mut window_totals = reference.as_ref().map(|_| WindowTotals::default());
    let mut stores_walked = 0usize;
    for (sample, store) in stores.iter().enumerate() {
        let mut over_this_store = SingleBaseEquality::default();
        let mut windows = reference.as_ref().map(|reference| {
            WindowRecomputation::over(
                reference,
                run_windows
                    .as_mut()
                    .map(|per_sample| std::mem::take(&mut per_sample[sample])),
            )
        });
        match &mut windows {
            Some(windows) => walk(store, &mut [&mut over_this_store, windows]),
            None => walk(store, &mut [&mut over_this_store]),
        }
        let label = store.file_name().map_or_else(
            || store.display().to_string(),
            |name| name.to_string_lossy().into_owned(),
        );
        // **The path, beside the basename that labels every row.** Two stores can share a
        // basename, and a figure quoted from this output has to say which file it came from.
        println!("{label}\tstore\t{}", store.display());
        println!("{label}\tsample\t{sample}");
        over_this_store.print(&label);
        if let Some(windows) = &mut windows {
            windows.close();
            windows.print(&label);
            window_totals
                .as_mut()
                .expect("the totals exist wherever the recomputation does")
                .add(windows);
        }
        assert!(
            over_this_store.found.one_base_records > 0,
            "{} holds no record covering one base, so the rule under test was never exercised \
             there",
            store.display(),
        );
        total.add(&over_this_store);
        stores_walked += 1;
    }
    // **Printed whatever the store count**, so that its absence means the run did not finish. A
    // walk that panics part-way through several stores leaves the per-store blocks already
    // printed, and this line is the only marker that separates that from a complete report.
    total.print("all-stores");
    if let Some(window_totals) = &window_totals {
        window_totals.print("all-stores");
    }
    println!(
        "all-stores\tstores-walked\t{stores_walked}/{}",
        stores.len()
    );

    // **After the report, not instead of it.** A counter-example is a finding this whole probe
    // exists to surface, so every count is printed first and the failure comes last.
    //
    // **The first assertion is not a formality.** Were the walk to yield nothing, or were
    // `region.len() == 1` to stop selecting one-base records, the second would pass over an empty
    // set and the step would report a rule confirmed by no records at all.
    assert!(
        total.found.one_base_records > 0,
        "no record covering one base was checked, so this run says nothing about spec \
         window_coverage.md §3.1's cheap rule",
    );
    let disagreeing = total.found.one_base_records_disagreeing;
    assert_eq!(
        disagreeing, 0,
        "the head count and the evidence disagree at {disagreeing} records covering one base, so \
         spec window_coverage.md §3.1's cheap rule does not hold on this store and becomes \
         \"build every body\"",
    );
    if let Some(window_totals) = &window_totals {
        window_totals.refuse_a_disagreement();
    }
}

/// The window recomputation summed over the stores walked.
///
/// **Only built where the recomputation ran**, so that its rows are absent rather than zero when
/// the reference was not given — a printed zero would read as "checked and found nothing".
#[derive(Default)]
struct WindowTotals {
    positions_observed: u64,
    windows_finalised: u64,
    /// `None` where no run's recorded windows were given.
    comparison: Option<ComparisonTotals>,
}

#[derive(Default)]
struct ComparisonTotals {
    loci: u64,
    agreed: u64,
    disagreed: u64,
    the_run_had_no_window: u64,
    neither_side_had_a_window: u64,
}

impl WindowTotals {
    /// **Destructured, so that a count added to either type is a compile error rather than a
    /// total that stays zero while the per-store rows show it.** The same device
    /// `WhatTheWalkFound::add` uses, and for the same reason.
    fn add(&mut self, store: &WindowRecomputation) {
        self.positions_observed += store.positions_observed;
        self.windows_finalised += store.windows.len() as u64;
        let Some(ComparisonWithTheRun {
            run_said,
            agreed,
            disagreed,
            the_run_had_no_window,
            neither_side_had_a_window,
            // Printed per store; a total of "the first twenty of each" would mean nothing.
            disagreements: _,
            where_the_run_had_no_window: _,
        }) = &store.comparison
        else {
            return;
        };
        let totals = self
            .comparison
            .get_or_insert_with(ComparisonTotals::default);
        totals.loci += run_said.len() as u64;
        totals.agreed += agreed;
        totals.disagreed += disagreed;
        totals.the_run_had_no_window += the_run_had_no_window;
        totals.neither_side_had_a_window += neither_side_had_a_window;
    }

    fn print(&self, label: &str) {
        println!("{label}\tpositions-observed\t{}", self.positions_observed);
        println!("{label}\twindows-finalised\t{}", self.windows_finalised);
        let Some(comparison) = &self.comparison else {
            return;
        };
        println!("{label}\tloci-the-run-built\t{}", comparison.loci);
        println!("{label}\tloci-agreeing-bit-for-bit\t{}", comparison.agreed);
        println!("{label}\tloci-disagreeing\t{}", comparison.disagreed);
        println!(
            "{label}\tloci-the-run-had-no-window-for\t{}",
            comparison.the_run_had_no_window
        );
        println!(
            "{label}\tloci-neither-side-has-a-window-for\t{}",
            comparison.neither_side_had_a_window
        );
    }

    /// **The comparison's own verdict, and it refuses an empty one.** A file naming no locus, or a
    /// comparison in which nothing was ever compared, would otherwise report zero disagreements
    /// and read as a pass.
    fn refuse_a_disagreement(&self) {
        let Some(comparison) = &self.comparison else {
            return;
        };
        assert!(
            comparison.loci > 0 && comparison.agreed > 0,
            "the run recorded {} loci and {} of them were compared against a window this walk \
             finalised, so this run says nothing about spec window_coverage.md §10's \
             recomputation oracle",
            comparison.loci,
            comparison.agreed,
        );
        assert_eq!(
            comparison.disagreed, 0,
            "the run's window differs from the whole-store recomputation at {} loci",
            comparison.disagreed,
        );
    }
}

/// Walk one store, building every body, and show each record to every measurement.
fn walk(store: &Path, measurements: &mut [&mut dyn RecordMeasurement]) {
    let mut reader =
        PspReader::open(store).unwrap_or_else(|why| panic!("opening {}: {why}", store.display()));
    let contigs: Vec<String> = reader
        .header()
        .contigs
        .iter()
        .map(|contig| contig.name.clone())
        .collect();

    for streamed in reader.records().expect("walking the store") {
        let streamed = streamed.expect("a record");
        let record = streamed
            .record
            .as_ref()
            .expect("this walk builds every body");
        for measurement in measurements.iter_mut() {
            measurement.observe(&streamed.head, record, &contigs);
        }
    }
}

fn carries_a_partial_witness(record: &SampleLocusObservations) -> bool {
    record
        .observations
        .iter()
        // Matched rather than compared against `Complete`: a third kind of witness would then be
        // a compile error here, where comparing would silently file it as partial.
        .any(|observation| match &observation.read_witness {
            ReadWitness::Complete => false,
            ReadWitness::Partial { .. } => true,
        })
}

/// **`head_says` is passed in rather than read off the head here**, so that the printed
/// counter-example always shows the two values the comparison actually used. Reading it again
/// from the head would let a probe whose comparison had been changed print two equal numbers
/// beside the word `counter-example`.
fn disagreement_at(
    head: &RecordHead,
    contigs: &[String],
    head_says: u32,
    evidence_says: u32,
    carries_a_partial_witness: bool,
) -> Disagreement {
    Disagreement {
        contig: contigs
            .get(head.region.contig.get() as usize)
            .cloned()
            .unwrap_or_else(|| format!("contig-{}", head.region.contig.get())),
        position: head.region.start.get(),
        head_says,
        evidence_says,
        carries_a_partial_witness,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record at `contig 0`, `start`..`start + bases - 1`, holding `observations`.
    fn a_record_over(
        start: u64,
        bases: u64,
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(start + bases - 1),
            },
            reference_bases: vec![b'A'; bases as usize].into_boxed_slice(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    fn an_observation(num_obs: u32, read_witness: ReadWitness) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(&b"A"[..]),
            read_witness,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: num_obs,
            q_sum: SummedLogError::NONE,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn a_head_over(
        record: &SampleLocusObservations,
        reads_compared_with_reference: u32,
    ) -> RecordHead {
        RecordHead {
            region: record.region,
            non_reference_reads: 0,
            reads_compared_with_reference,
            body_bytes: 0,
        }
    }

    fn what_the_measurement_makes_of(
        record: &SampleLocusObservations,
        head: &RecordHead,
    ) -> WhatTheWalkFound {
        let mut equality = SingleBaseEquality::default();
        equality.observe(head, record, &["chr1".to_string()]);
        equality.found
    }

    /// **The case the whole step is looking for, and it has to be visible when it happens.**
    /// A read whose witness is a reach — the STR path mints these — covers the single position
    /// of a one-base locus without having spanned the locus whole, so the head does not count it
    /// and `num_obs_along_locus` does. The measured answer on the tomato stores is that this
    /// never occurs; a probe that could not see it if it did would report the same thing.
    #[test]
    fn a_partial_witness_at_a_one_base_locus_shows_up_as_a_disagreement() {
        let witness = ReadWitness::from_left(1, LocusLen::from_positions(1))
            .expect("a run covering the one position of a one-base locus");
        assert_ne!(
            witness,
            ReadWitness::Complete,
            "a witness that is a reach stays a reach even where it covers the whole locus — \
             which is what makes this shape a counter-example at all",
        );
        let record = a_record_over(
            1_000,
            1,
            vec![
                an_observation(5, ReadWitness::Complete),
                an_observation(2, witness),
            ],
        );
        // **The fixture is a genuine disagreement, not a mislabelled one**: the two quantities
        // the rule chooses between really do differ on it. The head counts only the reads
        // compared whole; the evidence counts every read whose witness covers the position.
        assert_eq!(record.non_reference_and_compared_reads().1, 5);
        assert_eq!(record.num_obs_along_locus(), vec![7]);
        let head = a_head_over(&record, 5);

        let found = what_the_measurement_makes_of(&record, &head);
        assert_eq!(found.one_base_records, 1);
        assert_eq!(found.one_base_records_with_a_partial_witness, 1);
        assert_eq!(found.one_base_records_disagreeing, 1);
    }

    /// Spec §3.1 says a single-base record is generic by construction, and the rule treats a
    /// one-base tract the same way regardless. What must not happen is the classification
    /// reaching for the kind before the span: that would send such a record down the build
    /// branch and quietly drop it from the equality being measured.
    #[test]
    fn a_one_base_tract_record_is_still_checked_against_its_head() {
        let mut record = a_record_over(1_000, 1, vec![an_observation(4, ReadWitness::Complete)]);
        record.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a two-base motif"),
            left_flank: Box::from(&b"CCC"[..]),
            right_flank: Box::from(&b"GGG"[..]),
        });
        let head = a_head_over(&record, 4);

        let found = what_the_measurement_makes_of(&record, &head);
        assert_eq!(found.one_base_records, 1);
        assert_eq!(found.one_base_tract_records, 1);
        assert_eq!(found.wide_tract_records, 0);
        assert_eq!(found.one_base_records_disagreeing, 0);
    }

    /// The ordinary record: every read spanned the base, so the head's count is the evidence's.
    #[test]
    fn a_one_base_record_of_complete_witnesses_is_checked_and_agrees() {
        let record = a_record_over(
            1_000,
            1,
            vec![
                an_observation(3, ReadWitness::Complete),
                an_observation(2, ReadWitness::Complete),
            ],
        );
        let head = a_head_over(&record, 5);

        let found = what_the_measurement_makes_of(&record, &head);
        assert_eq!(found.one_base_records, 1);
        assert_eq!(found.one_base_records_disagreeing, 0);
        assert_eq!(found.one_base_records_with_a_partial_witness, 0);
    }

    /// **A record covering more than one base is not evidence about the cheap rule**, and
    /// counting it as one would let the headline count grow while the rule went unchecked.
    #[test]
    fn a_record_covering_more_than_one_base_is_counted_as_wide_and_not_checked() {
        let record = a_record_over(1_000, 4, vec![an_observation(3, ReadWitness::Complete)]);
        let head = a_head_over(&record, 3);

        let found = what_the_measurement_makes_of(&record, &head);
        assert_eq!(found.one_base_records, 0);
        assert_eq!(found.wide_generic_records, 1);
        assert_eq!(found.bases_in_wide_records, 4);
        // A wide *generic* record reports at its first base alone, so none of its four bases is
        // a position the build branch decides — which is the split the positions row needs.
        assert_eq!(found.bases_in_wide_tract_records, 0);
    }

    /// The other half of that split: a wide tract reports at every position of its span, so all
    /// of its bases are positions the build branch decides.
    #[test]
    fn a_wide_tract_records_bases_are_all_positions_the_build_branch_decides() {
        let mut record = a_record_over(1_000, 4, vec![an_observation(3, ReadWitness::Complete)]);
        record.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a two-base motif"),
            left_flank: Box::from(&b"CCC"[..]),
            right_flank: Box::from(&b"GGG"[..]),
        });
        let head = a_head_over(&record, 3);

        let found = what_the_measurement_makes_of(&record, &head);
        assert_eq!(found.wide_tract_records, 1);
        assert_eq!(found.bases_in_wide_records, 4);
        assert_eq!(found.bases_in_wide_tract_records, 4);
    }

    /// **Two absent windows are the same measurement, and two a last bit apart are not.** The
    /// first is what breaks under float equality — an absent window is a pair of `NaN`s, so a
    /// correct run would be reported as disagreeing at every centre the position floor silenced.
    /// The second is what breaks under a comparison too loose to see the last bit, which is the
    /// resolution this oracle claims.
    #[test]
    fn absent_windows_agree_and_a_last_bit_apart_is_a_disagreement() {
        assert!(
            the_same_window(WindowCoverage::absent(), WindowCoverage::absent()),
            "two absent windows must compare equal, or the floor's own windows read as \
             disagreements",
        );
        let one = WindowCoverage {
            gc_fraction: 0.5,
            mean_depth: 3.25,
        };
        assert!(the_same_window(one, one));
        let a_bit_apart = WindowCoverage {
            gc_fraction: one.gc_fraction,
            mean_depth: f32::from_bits(one.mean_depth.to_bits() + 1),
        };
        assert!(
            !the_same_window(one, a_bit_apart),
            "the comparison is bit for bit, and this pair differs in the last bit",
        );
    }
}
