//! **Four questions about the window-coverage measurement, asked of a real store.** Does a
//! one-base record's head really say what its evidence says (plan step B2, spec
//! `window_coverage.md` §3.1)? Does a calling run's per-locus window equal the one a straight
//! walk of the same store computes (plan step C3, spec §10)? And how many covered positions is
//! each window built from, which is what spec §3.3's floor is a threshold on (plan step D1)? And
//! where does the histogram's depth axis get cut, and how much of a sample falls off the end of
//! it, which is what spec §3.4's three bin constants decide (plan step D2)?
//!
//! ```text
//! # the first question alone — plan step B2's command, unchanged
//! cargo run --release --example ng_window_coverage_probe -- <a store.psp> [more stores...]
//!
//! # the first two — the run is asked to write down what it read, then the walk checks it
//! NG_WINDOW_COVERAGE_FILE=windows.tsv pop_var_caller call-from-psps --psp <a store.psp> ...
//! cargo run --release --example ng_window_coverage_probe -- \
//!     --reference <reference.fa> --windows-from-the-run windows.tsv \
//!     <a store.psp> [more stores...]
//!
//! # the third and fourth — no run needed, only the stores and the reference their GC comes from
//! cargo run --release --example ng_window_coverage_probe -- \
//!     --reference <reference.fa> --covered-positions-per-window --bin-scheme \
//!     <a store.psp> [more stores...]
//! ```
//!
//! Without `--reference` only the first question is asked, which is what plan step B2 ran. With
//! it, the walk also recomputes **every position's window** — the whole-store recomputation of
//! spec §10 — and, given `--windows-from-the-run`, checks the run's windows against its own. The
//! run writes that file when
//! [`recorded_windows::PATH_VARIABLE`](pop_var_caller::run::cohort_merge::recorded_windows::PATH_VARIABLE)
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
//! **One difference is expected and is not a defect**, and it is about the *windows*, not the
//! histogram. Both sides call `finish`, which closes the centres in the last half-window of a
//! sample's stream — the ones no cover can reach, since there is no later position to complete
//! them. This walk keeps those; a run finishes only to take the histogram out, by which time
//! every builder has run, so nothing reads them and they reach no record. They are counted and
//! reported separately from a genuine disagreement. The histogram carries them on both sides,
//! which is why the histograms compare equal where the per-locus windows do not.
//!
//! # The third question
//!
//! A window is silenced — it reports nothing, and trains no yardstick — when it holds fewer than
//! `MIN_WINDOW_POSITIONS` covered positions (spec §3.3). Setting that floor means knowing how
//! those counts are spread over a real store's windows.
//!
//! **The spread is read off the shipped accumulator, at one accumulator per candidate floor.** One
//! such accumulator is an **arm**, and [`FloorSweep`] is the set of them. An arm's silenced count
//! at floor `F` *is* the number of windows holding fewer than `F` positions, so the arms' answers
//! in floor order are the cumulative distribution, and nothing in this file recomputes a window.
//! The arm at the floor the run actually ships must agree with the walk's own count of absent
//! windows, which is what ties the sweep to the measurement it is about.
//!
//! # The fourth question
//!
//! Every window that clears the floor is folded into a per-sample histogram, binned by GC and by
//! depth. The depth axis is not configured but **fitted to the sample**: the accumulator holds its
//! first `DEPTH_SCALE_WINDOWS` windows back, takes the median of their mean depths, and cuts
//! `DEPTH_BINS` bins spanning `DEPTH_RANGE_IN_MEDIANS` times that median. Everything above lands
//! in one overflow column, and the model fit that reads the histogram rejects the sample outright
//! once more than a fifth of its windows are there. Those three constants were starting values
//! until this measurement was run; it kept all three, and re-running it is how a later change to
//! any of them is checked.
//!
//! **The same device answers this**: one shipped accumulator per candidate setting, differing from
//! the run's in one field, fed the identical positions. Each reports the axis it fitted and how
//! much of the sample fell off the end of it, and the arm configured exactly as the run is
//! checked against the walk's own histogram.
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
//! empty live set** (`src/run/psp_source.rs`'s `build`). The two agree only while the encoder
//! writes no chain ids, which is true today and stops being true at the psp path's Milestone E.
//! When it does, this measurement has to be retaken rather than cited.
//!
//! # One walk, several measurements
//!
//! [`walk`] decodes every record once and hands it to whatever [`RecordMeasurement`]s the caller
//! passed; it knows nothing about what they count. [`SingleBaseEquality`] is the first question's
//! and [`WindowRecomputation`] the second's; neither knows about the other, and a measurement that
//! needs only the record is a third implementor and one more entry in `main`'s slice.
//!
//! **The third and fourth questions are not those, and ride inside [`WindowRecomputation`]
//! instead.** What they count is not a record but a covered position with its reference base and
//! its depth, which exists only after the recomputation has looked the bases up and decoded the
//! record's reported depths. A separate implementor would need its own copy of that — a second
//! reference accessor, a second contig check, a second call of the depth rule — and its answer
//! would then be checkable only against itself. Riding along, each is fed the identical four
//! values and each has one arm the recomputation can check: [`FloorSweep`]'s arm at the shipped
//! floor against the walk's own count of absent windows, and [`BinSchemeSweep`]'s arm at the
//! run's own configuration against the walk's own histogram.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use pop_var_caller::fasta::ContigList;
use pop_var_caller::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use pop_var_caller::paralog::coverage_model::DEFAULT_MAX_OVERFLOW_FRACTION;
use pop_var_caller::psp::{PspReader, RecordHead};
use pop_var_caller::ref_seq::{RefSeq, WindowedRefSeq};
use pop_var_caller::reference_info::{
    ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
};
use pop_var_caller::run::cohort_merge::observation_cache::LocusSummary;
use pop_var_caller::run::cohort_merge::recorded_windows::{
    histograms_beside, read_a_row, write_the_histogram,
};
use pop_var_caller::types::{ContigId, GenomePosition, GenomeRegion, Position};
use pop_var_caller::window_coverage::depth::{EvidenceForOneRecord, for_each_reported_depth};
use pop_var_caller::window_coverage::{
    self, CoverageByGcHistogram, SampleHistogram, WindowCoverage, WindowCoverageAccumulator,
    WindowCoverageConfig,
};

#[cfg(test)]
use pop_var_caller::locus_generation::{LocusLen, SequenceObservation, SsrDetail};
#[cfg(test)]
use pop_var_caller::types::{Motif, ReadGroupId, SummedLogError};

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
    /// The accumulators that measure how many covered positions each window holds, fed the same
    /// positions as the one above and differing from it only in their floor. `None` unless
    /// `--covered-positions-per-window` was given, because they repeat the window arithmetic once
    /// per floor; also `None` once [`close`](Self::close) has finished the pass, which is when the
    /// field below is filled.
    sweep: Option<FloorSweep>,
    /// What that measurement said, set by [`close`](Self::close) and `None` where the sweep was
    /// not asked for.
    windows_under_each_floor: Option<WindowsUnderEachFloor>,
    /// The accumulators that measure where the depth axis is cut, one per candidate setting, fed
    /// the same positions as the one above. `None` unless `--bin-scheme` was given, and `None`
    /// again once [`close`](Self::close) has finished the pass.
    bin_schemes: Option<BinSchemeSweep>,
    /// What that measurement said, set by [`close`](Self::close). One entry an arm, in the order
    /// [`the_bin_schemes_asked_about`] gives them.
    ///
    /// **Printed per store and never summed**, unlike [`Self::windows_under_each_floor`]: an axis
    /// is fitted from one sample's own depth, so two samples' axes have no total — the row a
    /// reader wants across stores is the worst overflow, and that is read off the rows.
    each_bin_scheme: Option<Vec<OneBinSchemeAnswer>>,
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
    /// Which sample of the run this store is — the index the recorded rows and histogram line
    /// are keyed by, and the position this store was listed at.
    sample: usize,
    /// This sample's finished histogram, or the reason it has none — `None` until
    /// [`close`](Self::close) has run.
    histogram: Option<SampleHistogram>,
    /// The same, as the run wrote it down, when a run's histograms were given. **Not `run_said`**,
    /// which in this file is the map of the run's *windows*.
    histogram_the_run_wrote: Option<String>,
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

/// The configuration a calling run gives its accumulators, **spelled with every field and no
/// `..`**, so that a setting added to [`WindowCoverageConfig`] has to be answered for here rather
/// than defaulting silently to something the run does not use.
///
/// One function rather than a literal per accumulator, because this file now builds several and a
/// walk configured unlike the run measures something else. An accumulator that differs from the
/// run in one setting names that setting and takes the rest from here with `..`, so it inherits
/// the run's value for anything added later rather than a default nobody chose.
fn the_configuration_a_run_uses() -> WindowCoverageConfig {
    WindowCoverageConfig {
        window_bp: window_coverage::WINDOW_BP,
        gc_bins: window_coverage::GC_BINS,
        depth_bins: window_coverage::DEPTH_BINS,
        depth_scale_windows: window_coverage::DEPTH_SCALE_WINDOWS,
        depth_range_in_medians: window_coverage::DEPTH_RANGE_IN_MEDIANS,
        min_window_positions: window_coverage::MIN_WINDOW_POSITIONS,
    }
}

/// **How many covered positions each window was built from, as a distribution** — the quantity
/// spec `window_coverage.md` §3.3's floor is a threshold on.
///
/// A window centred at `p` holds the sample's covered positions within 250 bases either side of
/// `p` on one contig. Choosing the floor means knowing how those counts are spread; the counts
/// themselves are private to the accumulator and no method hands one back.
///
/// **So the spread is read off the shipped accumulator rather than recomputed here.** The sweep
/// holds one accumulator per candidate floor — an **arm** — fed exactly the positions
/// [`WindowRecomputation`] is fed, and identical to it in every setting but the floor. The windows
/// an arm silences *are* the windows holding fewer than that arm's floor, so the arms' answers,
/// read in floor order, are the cumulative distribution. Nothing in this file reimplements the
/// window, which is what stops the measurement from being a check of one copy of the rule against
/// another copy of it.
///
/// **The accumulators live only for the pass.** [`close`](Self::close) consumes the sweep and
/// hands back a [`WindowsUnderEachFloor`], which is what gets printed and what a total is summed
/// from — so there is no half-closed sweep to describe, and no state that a total would have to
/// carry an empty copy of.
struct FloorSweep {
    /// One arm per candidate floor, in ascending floor order.
    arms: Vec<OneFloor>,
}

/// One arm of a [`FloorSweep`] while the pass is running: the shipped accumulator at one candidate
/// floor, and what it has silenced so far.
struct OneFloor {
    accumulator: WindowCoverageAccumulator,
    counted: WindowsUnderOneFloor,
}

/// What one candidate floor costs — over one store, or over every store summed.
#[derive(Clone)]
struct WindowsUnderOneFloor {
    /// The floor this row asks about: a window holding fewer than this many covered positions
    /// comes back absent.
    floor: u32,
    /// Windows finalised at this floor. The same number at every floor, which
    /// [`FloorSweep::close`] asserts, because the floor decides whether a window speaks and never
    /// whether it exists.
    windows: u64,
    /// Windows silenced at this floor: those holding fewer than `floor` covered positions.
    silenced: u64,
}

impl WindowsUnderOneFloor {
    /// Count one window this floor's arm finalised, and whether the floor silenced it.
    ///
    /// **One counting rule, called from both drains** — the windows [`FloorSweep::observe`] pops
    /// as the stream advances, and the tail [`FloorSweep::close`] takes from `finish`. Written
    /// out twice it could be right on one path and inverted on the other, and the row would still
    /// read as a distribution: [`FloorSweep::close`]'s checks compare arms against the walk, and
    /// every arm would be wrong the same way. A stream short enough that no centre's right edge is
    /// reached exercises only the second path, which is the shape every ten-position test here
    /// has.
    fn count(&mut self, window: WindowCoverage) {
        self.windows += 1;
        self.silenced += u64::from(window.is_absent());
    }
}

/// A closed sweep's answer — one row a floor, ascending, and the cumulative distribution of
/// covered positions per window read in that order.
#[derive(Clone)]
struct WindowsUnderEachFloor {
    floors: Vec<WindowsUnderOneFloor>,
}

/// The floors the sweep asks about — spaced most closely over the range spec §3.3 is choosing in,
/// and thinning out above it.
///
/// **The two ends are anchors rather than candidates.** A floor of 1 silences nothing, because a
/// centre always lies in its own window, so a count above zero there is a defect and not a
/// reading — [`FloorSweep::close`] asserts it. A floor of 501 is the highest the configuration
/// allows, since a 500-base window spans 501 reference positions, so its arm counts the windows
/// whose ground is *not* completely covered: every window the floor could ever be raised far
/// enough to reach.
const CANDIDATE_FLOORS: [u32; 25] = [
    1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400,
    450, 501,
];

/// [`CANDIDATE_FLOORS`] with the floor the run actually ships folded in, ascending and without
/// repeats.
///
/// **The shipped floor is always among them, wherever it is set**, so that
/// [`FloorSweep::close`]'s check of the sweep against the recomputation beside it cannot quietly
/// stop applying if [`MIN_WINDOW_POSITIONS`](window_coverage::MIN_WINDOW_POSITIONS) is ever
/// changed.
fn the_floors_asked_about() -> Vec<u32> {
    let mut floors: Vec<u32> = CANDIDATE_FLOORS.to_vec();
    floors.push(window_coverage::MIN_WINDOW_POSITIONS);
    floors.sort_unstable();
    floors.dedup();
    floors
}

impl FloorSweep {
    /// One arm per floor in [`the_floors_asked_about`], each a fresh accumulator configured as a
    /// run's but for its own floor.
    ///
    /// The `..` takes the run's own configuration as its base, so a setting added to
    /// [`WindowCoverageConfig`] reaches every arm with the value the run gives it — and
    /// [`the_configuration_a_run_uses`] still has to name it before this compiles.
    fn new() -> Self {
        Self {
            arms: the_floors_asked_about()
                .into_iter()
                .map(|floor| OneFloor {
                    accumulator: WindowCoverageAccumulator::new(WindowCoverageConfig {
                        min_window_positions: floor,
                        ..the_configuration_a_run_uses()
                    }),
                    counted: WindowsUnderOneFloor {
                        floor,
                        windows: 0,
                        silenced: 0,
                    },
                })
                .collect(),
        }
    }

    /// One covered position, into every arm.
    ///
    /// **The windows are counted and dropped**, not kept: this measurement is about how many
    /// positions each window held, and an arm's windows are the same centres the recomputation
    /// beside it already stores.
    fn observe(&mut self, contig: ContigId, position: Position, reference_base: u8, depth: u32) {
        for arm in &mut self.arms {
            arm.accumulator
                .observe(contig, position, reference_base, depth);
            while let Some((_, window)) = arm.accumulator.pop_ready() {
                arm.counted.count(window);
            }
        }
    }

    /// Finish every arm and hand back the distribution, having checked the six things that would
    /// make it a fiction.
    ///
    /// `windows_the_walk_finalised` and `silenced_at_the_shipped_floor` come from the
    /// recomputation this sweep rides along with — a separate accumulator over the same positions
    /// — so an arm disagreeing with it means the two were not fed alike, and every number here is
    /// about a different stream.
    ///
    /// **What none of the six can see is the arms' own configuration.** Every arm finalises one
    /// window per covered position whatever its window width, so a bank built with a window the
    /// run does not use passes all of these and prints a distribution of a window nothing else
    /// computes. Only the unit tests pin the width, by fixing the counts a 500-base window
    /// produces over 600 consecutive positions. The practical guard is to walk at least one store
    /// whose own windows the floor silences something in, and
    /// [`WindowsUnderEachFloor::print`]'s first row says whether a store was one.
    fn close(
        self,
        windows_the_walk_finalised: u64,
        silenced_at_the_shipped_floor: u64,
    ) -> WindowsUnderEachFloor {
        let floors: Vec<WindowsUnderOneFloor> = self
            .arms
            .into_iter()
            .map(|arm| {
                let OneFloor {
                    accumulator,
                    mut counted,
                } = arm;
                let (tail, _) = accumulator.finish();
                for (_, window) in tail {
                    counted.count(window);
                }
                counted
            })
            .collect();
        assert!(
            windows_the_walk_finalised > 0,
            "the walk finalised no window, so this sweep says nothing about how many covered \
             positions a window holds",
        );
        // **Both floor-keyed checks below have to find their floor.** Each is written as "if this
        // is the row at floor `F`", which on a floor list that has stopped holding `F` is silently
        // no check at all rather than a failure — and that list is edited by hand.
        for anchor in [1, window_coverage::MIN_WINDOW_POSITIONS] {
            assert!(
                floors.iter().any(|at| at.floor == anchor),
                "no arm at a floor of {anchor}, so the check written against it never ran",
            );
        }
        for at in &floors {
            assert_eq!(
                at.windows, windows_the_walk_finalised,
                "the arm at a floor of {} finalised {} windows and the walk beside it {}; the \
                 floor decides whether a window speaks, never whether it exists, so the two were \
                 not fed the same positions",
                at.floor, at.windows, windows_the_walk_finalised,
            );
            if at.floor == 1 {
                assert_eq!(
                    at.silenced, 0,
                    "the arm at a floor of 1 silenced {} windows; a centre always lies in its own \
                     window, so no window can hold fewer than one covered position",
                    at.silenced,
                );
            }
            if at.floor == window_coverage::MIN_WINDOW_POSITIONS {
                assert_eq!(
                    at.silenced, silenced_at_the_shipped_floor,
                    "the arm at the shipped floor of {} silenced {} windows and the walk beside \
                     it {}, over the same positions",
                    at.floor, at.silenced, silenced_at_the_shipped_floor,
                );
            }
        }
        assert!(
            floors
                .windows(2)
                .all(|pair| pair[0].silenced <= pair[1].silenced),
            "a higher floor silenced fewer windows than a lower one, which no set of window \
             counts can produce",
        );
        WindowsUnderEachFloor { floors }
    }
}

impl WindowsUnderEachFloor {
    /// One row a floor, in the same three-then-tagged shape the counter-example rows use: the
    /// floor, then how many windows it silences and how many in every 10,000.
    ///
    /// **A share in every 10,000 rather than a percentage**, because the end of this distribution
    /// the floor is being chosen for is the sparse one, where a percentage rounds to zero.
    ///
    /// The denominator is not printed: [`FloorSweep::close`] has asserted it equal to the
    /// `windows-finalised` row above.
    ///
    /// **The first row says whether the tie to the walk was a comparison or two zeroes.**
    /// [`FloorSweep::close`] checks the row at the shipped floor against the walk's own count of
    /// absent windows; on a store where the walk silences nothing, that check compares 0 with 0
    /// and would pass over arms fed anything. A store printing 0 there contributes rows that
    /// nothing outside the sweep has confirmed, and the reader has to see it beside them.
    fn print(&self, label: &str) {
        let windows = self.floors.first().map_or(0, |at| at.windows);
        let tied = self
            .floors
            .iter()
            .find(|at| at.floor == window_coverage::MIN_WINDOW_POSITIONS)
            .map_or(0, |at| at.silenced);
        println!(
            "{label}\tsweep-tied-to-the-walk\tfloor={}\twindows={tied}",
            window_coverage::MIN_WINDOW_POSITIONS,
        );
        for at in &self.floors {
            let per_ten_thousand = if windows == 0 {
                0.0
            } else {
                at.silenced as f64 * 10_000.0 / windows as f64
            };
            println!(
                "{label}\twindows-under-a-floor-of\t{}\twindows={}\tper-10000={per_ten_thousand:.2}",
                at.floor, at.silenced,
            );
        }
    }

    /// Add another store's answer into this one, floor by floor.
    ///
    /// **Destructured with every field named**, so that a count added to
    /// [`WindowsUnderOneFloor`] is a compile error rather than a total that stays zero while the
    /// per-store rows show it.
    fn add(&mut self, other: &WindowsUnderEachFloor) {
        assert_eq!(
            self.floors.len(),
            other.floors.len(),
            "two sweeps over different numbers of floors cannot be added",
        );
        for (into, from) in self.floors.iter_mut().zip(&other.floors) {
            let WindowsUnderOneFloor {
                floor,
                windows,
                silenced,
            } = from;
            assert_eq!(
                into.floor, *floor,
                "two sweeps' floors are not in the same order",
            );
            into.windows += windows;
            into.silenced += silenced;
        }
    }
}

/// **Where the histogram's depth axis is cut, and how much of the sample falls off the end of
/// it** — the three constants spec `window_coverage.md` §3.4 left to measurement, and that this
/// sweep is what settled.
///
/// The axis is not configured but fitted: the accumulator holds its first `depth_scale_windows`
/// windows back, takes the median of their mean depths, and cuts `depth_bins` bins so they span
/// `depth_range_in_medians` times that median. Everything above lands in one overflow column, and
/// the model fit that reads this histogram **rejects a sample outright** once more than a fifth of
/// its windows are in that column ([`THE_OVERFLOW_SHARE_THE_FIT_ALLOWS`]), because at that point
/// the sample's own single-copy peak has left the range and the regular bins hold noise.
///
/// **So the same device the floor sweep above uses answers all three**: one shipped accumulator
/// per setting, fed exactly the positions [`WindowRecomputation`] is fed and differing from it in
/// one field. Varying `depth_scale_windows` says how far the fitted median moves when the width is
/// fitted from more of the store, up to every window of it. Varying `depth_range_in_medians` says
/// what the factor of ten buys and costs: a wider range overflows less and resolves more coarsely,
/// and the two move in opposite directions over the same bins.
///
/// **The two answers are read together, not one at a time.** A scale sample that fits a median
/// below the whole store's costs a proportional slice off the top of the axis, so what a short
/// scale sample costs is paid in overflow — which is the range's quantity.
struct BinSchemeSweep {
    /// One arm per setting, in the order [`the_bin_schemes_asked_about`] gives them.
    arms: Vec<OneBinScheme>,
}

/// One arm of a [`BinSchemeSweep`] while the pass is running: the shipped accumulator at one
/// candidate setting.
struct OneBinScheme {
    setting: BinSchemeSetting,
    accumulator: WindowCoverageAccumulator,
}

/// Which of the two settings an arm varies, and to what — the whole of what makes one arm differ
/// from another, and what a row is labelled by, so that a reader does not have to know the arm's
/// position in a list to know what it changed.
///
/// No `Eq`, because one variant carries an `f64`. The only comparison this file makes is against
/// [`AsTheRunFitsIt`](Self::AsTheRunFitsIt), which carries nothing.
#[derive(Clone, Copy, PartialEq, Debug)]
enum BinSchemeSetting {
    /// The run's own configuration, unchanged. Exactly one arm is this, and
    /// [`BinSchemeSweep::close`] checks its answer against the walk's own histogram.
    AsTheRunFitsIt,
    /// `depth_scale_windows` moved to this many windows, everything else the run's.
    ScaleSampleOf(u32),
    /// `depth_range_in_medians` moved to this multiple, everything else the run's.
    RangeOf(f64),
}

/// What one arm's histogram came out as.
struct OneBinSchemeAnswer {
    setting: BinSchemeSetting,
    fitted: WhatOneArmFitted,
}

/// An arm either cut an axis or it has no histogram to cut one from, and **which silence it was is
/// itself a reading** — the three are three different reports about the sample, not one hole.
enum WhatOneArmFitted {
    /// The axis the arm cut, and what fell off the end of it.
    Axis(FittedAxis),
    /// No histogram, and which of [`SampleHistogram`]'s three silences it was.
    NoHistogram(&'static str),
}

/// A fitted depth axis and what fell off the end of it.
#[derive(Clone, Copy)]
struct FittedAxis {
    /// The width of one regular depth bin, in mean-depth units.
    depth_bin_width: f64,
    /// The median window depth the width was fitted from — `width * depth_bins /
    /// depth_range_in_medians`, recovered rather than re-derived, since the accumulator keeps it
    /// nowhere.
    median_depth: f64,
    /// The depth at which the regular bins stop and the overflow column begins.
    top_of_the_range: f64,
    /// Windows that landed in a regular depth bin — the counts of every cell except each GC
    /// row's last, which is `coverage_model.rs`'s `regular_total`.
    ///
    /// **Read off the cells, not off the accumulator's `windows_folded` counter.** The two agree
    /// — a folded window always reaches a cell — and [`OneArmsAnswer::of`] asserts they do; but
    /// the fraction the fit judges a sample by is a ratio of cell counts, and reading it from
    /// anywhere else would make this probe's number a near relative of the fit's rather than the
    /// fit's.
    in_a_regular_bin: u64,
    /// Windows that landed in the overflow column: their mean depth was at or above
    /// `top_of_the_range`. `coverage_model.rs`'s `overflow_total`.
    overflowed: u64,
}

impl FittedAxis {
    /// Windows folded into the histogram at all — every window that reached a cell.
    fn folded(&self) -> u64 {
        self.in_a_regular_bin + self.overflowed
    }

    /// The share the model fit's own guard is a threshold on: overflowed windows over every
    /// window that reached a cell.
    ///
    /// **Computed the way `coverage_model.rs` computes it** — the overflow column over the
    /// regular bins plus the overflow column, both summed out of the cells — so the number here
    /// is the one that would be compared against the guard, not a near relative of it.
    ///
    /// A sample with no window in any cell has no histogram at all, so this never divides by
    /// zero in a run; the guard is written anyway, and answers `0.0`, because the fit's own guard
    /// also declines to reject a sample it has counted nothing for.
    fn overflow_fraction(&self) -> f64 {
        if self.folded() == 0 {
            0.0
        } else {
            self.overflowed as f64 / self.folded() as f64
        }
    }
}

/// The largest share of a sample's windows the coverage model fit tolerates in the overflow
/// column before it rejects the sample outright — 0.20, a fifth.
///
/// **Bound to production's own constant rather than copied**, so that this measurement is read
/// against the number that would actually judge the sample and cannot drift from it in silence.
/// Reading it is not a change to production: it is a `pub const` of a `pub` module, and the local
/// name is here only to say in this file's own words what the number is.
const THE_OVERFLOW_SHARE_THE_FIT_ALLOWS: f64 = DEFAULT_MAX_OVERFLOW_FRACTION;

/// The settings the bin-scheme sweep asks about.
///
/// **The first is the run's own**, which is what ties the sweep to the walk beside it. The scale
/// samples bracket the shipped 10,000 by two orders each way and end with one no store here
/// reaches, which fits the width from every window the sample has. The ranges bracket the shipped
/// factor of ten from a quarter of it to four times.
fn the_bin_schemes_asked_about() -> Vec<BinSchemeSetting> {
    let mut settings = vec![BinSchemeSetting::AsTheRunFitsIt];
    settings.extend(
        [100, 1_000, 100_000, 1_000_000, u32::MAX]
            .into_iter()
            .map(BinSchemeSetting::ScaleSampleOf),
    );
    settings.extend(
        [2.5, 5.0, 20.0, 40.0]
            .into_iter()
            .map(BinSchemeSetting::RangeOf),
    );
    settings
}

impl BinSchemeSetting {
    /// The configuration this arm runs — the run's, with at most one field moved.
    fn configuration(self) -> WindowCoverageConfig {
        let run = the_configuration_a_run_uses();
        match self {
            BinSchemeSetting::AsTheRunFitsIt => run,
            BinSchemeSetting::ScaleSampleOf(windows) => WindowCoverageConfig {
                depth_scale_windows: windows,
                ..run
            },
            BinSchemeSetting::RangeOf(medians) => WindowCoverageConfig {
                depth_range_in_medians: medians,
                ..run
            },
        }
    }

    /// How a row labels itself: what was moved and to what.
    fn as_text(self) -> String {
        match self {
            BinSchemeSetting::AsTheRunFitsIt => "as-the-run-fits-it".to_owned(),
            BinSchemeSetting::ScaleSampleOf(u32::MAX) => "scale-sample=every-window".to_owned(),
            BinSchemeSetting::ScaleSampleOf(windows) => format!("scale-sample={windows}"),
            BinSchemeSetting::RangeOf(medians) => format!("range-in-medians={medians}"),
        }
    }
}

impl BinSchemeSweep {
    /// One arm per setting in [`the_bin_schemes_asked_about`], each a fresh accumulator configured
    /// as a run's but for the one field that setting moves — see
    /// [`BinSchemeSetting::configuration`], which is where the `..` over the run's own
    /// configuration is written.
    fn new() -> Self {
        Self {
            arms: the_bin_schemes_asked_about()
                .into_iter()
                .map(|setting| OneBinScheme {
                    setting,
                    accumulator: WindowCoverageAccumulator::new(setting.configuration()),
                })
                .collect(),
        }
    }

    /// One covered position, into every arm.
    ///
    /// **The windows are dropped as they are popped.** This measurement is about the histogram
    /// each arm ends with, and the windows themselves are the recomputation's, beside it.
    fn observe(&mut self, contig: ContigId, position: Position, reference_base: u8, depth: u32) {
        for arm in &mut self.arms {
            arm.accumulator
                .observe(contig, position, reference_base, depth);
            while arm.accumulator.pop_ready().is_some() {}
        }
    }

    /// Finish every arm and read its axis off its histogram, having checked the three things that
    /// would make the answers a fiction.
    ///
    /// `the_walks_histogram` is the recomputation's own, fitted with the run's configuration, and
    /// the arm that runs that same configuration must reproduce it exactly — the one place this
    /// sweep and the walk answer the same question, and the only check here that compares against
    /// something computed outside it. **What is compared is the printed row**, so the tolerance is
    /// the precision the row is printed at, and a field of [`FittedAxis`] enters the comparison as
    /// soon as it enters the row.
    fn close(self, the_walks_histogram: &SampleHistogram) -> Vec<OneBinSchemeAnswer> {
        let answers: Vec<OneBinSchemeAnswer> = self
            .arms
            .into_iter()
            .map(|arm| {
                let OneBinScheme {
                    setting,
                    accumulator,
                } = arm;
                let (_, histogram) = accumulator.finish();
                OneBinSchemeAnswer::of(setting, &histogram)
            })
            .collect();
        // Found rather than filtered for, so that a sweep that has lost the arm the check is
        // written against fails here instead of running no check at all.
        let as_the_run_ran_it = answers
            .iter()
            .find(|answer| answer.setting == BinSchemeSetting::AsTheRunFitsIt)
            .expect(
                "no arm runs the configuration the run uses, so the check against the walk \
                 never ran",
            );
        let the_walk =
            OneBinSchemeAnswer::of(BinSchemeSetting::AsTheRunFitsIt, the_walks_histogram);
        assert_eq!(
            as_the_run_ran_it.as_a_row(),
            the_walk.as_a_row(),
            "the arm configured as the run is fitted came out differently from the walk beside \
             it, over the same positions",
        );
        assert!(
            answers
                .iter()
                .any(|answer| answer.axis().is_some_and(|axis| axis.folded() > 0)),
            "no arm folded a window, so this sweep says nothing about where the depth axis \
             should be cut",
        );
        answers
    }
}

impl OneBinSchemeAnswer {
    /// Read one arm's axis off the histogram it finished with.
    ///
    /// **The three silences are named here rather than rendered with `Debug`**, so that a variant
    /// added to [`SampleHistogram`] is a compile error rather than a fourth silence reported under
    /// whatever name its declaration happens to carry.
    fn of(setting: BinSchemeSetting, histogram: &SampleHistogram) -> Self {
        let fitted = match histogram {
            SampleHistogram::NoWindowFinalised => {
                WhatOneArmFitted::NoHistogram("no-window-finalised")
            }
            SampleHistogram::EveryWindowUnderTheFloor => {
                WhatOneArmFitted::NoHistogram("every-window-under-the-floor")
            }
            SampleHistogram::MedianDepthNotPositive => {
                WhatOneArmFitted::NoHistogram("median-depth-not-positive")
            }
            SampleHistogram::Fitted(fitted) => {
                let CoverageByGcHistogram {
                    depth_bin_width,
                    depth_bins,
                    gc_bins,
                    windows_folded,
                    counts,
                    // Not part of the axis: how wide a window is, and how many the floor
                    // silenced, belong to the floor sweep above and the walk prints both already.
                    window_bp: _,
                    windows_under_the_floor: _,
                } = fitted;
                // Row-major `[gc_bin][depth_bin]` with one overflow column after each row's
                // regular bins, which is how `coverage_model.rs`'s own layout reads it. Both
                // totals are summed out of the cells, because the fraction the fit judges a
                // sample by is a ratio of exactly these two.
                let columns = *depth_bins as usize + 1;
                let overflowed: u64 = (0..*gc_bins as usize)
                    .map(|gc_bin| u64::from(counts[gc_bin * columns + *depth_bins as usize]))
                    .sum();
                let in_a_regular_bin: u64 = (0..*gc_bins as usize)
                    .flat_map(|gc_bin| {
                        (0..*depth_bins as usize).map(move |depth_bin| gc_bin * columns + depth_bin)
                    })
                    .map(|cell| u64::from(counts[cell]))
                    .sum();
                // The accumulator counts one folded window per cell it wrote, so the cells and
                // the counter are two readings of one number. They can part only if a cell
                // saturated at `u32::MAX`, which needs one GC-and-depth cell to hold 4.3 billion
                // windows — a reference above about 4.3 Gbp, larger than either benchmark. If it
                // ever happens the overflow fraction reads low, which is the direction that makes
                // a sample look acceptable, so it fails here rather than being reported.
                assert_eq!(
                    in_a_regular_bin + overflowed,
                    *windows_folded,
                    "the histogram's cells hold {} windows and the accumulator counted {} folded, \
                     so a cell has saturated and the overflow fraction would read low",
                    in_a_regular_bin + overflowed,
                    windows_folded,
                );
                let top_of_the_range = depth_bin_width * f64::from(*depth_bins);
                WhatOneArmFitted::Axis(FittedAxis {
                    depth_bin_width: *depth_bin_width,
                    // The width is `median * range / bins`, so the median it was fitted from
                    // is the width times bins over range — recovered, because nothing keeps
                    // it. The range is the arm's own, not the run's.
                    median_depth: depth_bin_width * f64::from(*depth_bins)
                        / setting.configuration().depth_range_in_medians,
                    top_of_the_range,
                    in_a_regular_bin,
                    overflowed,
                })
            }
        };
        Self { setting, fitted }
    }

    /// The axis this arm cut, or `None` where the sample got no histogram to cut one from.
    fn axis(&self) -> Option<FittedAxis> {
        match self.fitted {
            WhatOneArmFitted::Axis(axis) => Some(axis),
            WhatOneArmFitted::NoHistogram(_) => None,
        }
    }

    /// One arm's answer as the text a row carries — also what the check against the walk
    /// compares, so that a field added to [`FittedAxis`] is compared as soon as it is printed.
    ///
    /// **The axis is destructured with every field named**, so adding one is a compile error here
    /// rather than a field that quietly stays out of both the row and the check.
    fn as_a_row(&self) -> String {
        match self.fitted {
            WhatOneArmFitted::Axis(axis) => {
                let FittedAxis {
                    depth_bin_width,
                    median_depth,
                    top_of_the_range,
                    in_a_regular_bin,
                    overflowed,
                } = axis;
                format!(
                    "median-depth={median_depth:.4}\tdepth-bin-width={depth_bin_width:.6}\t\
                     top-of-the-range={top_of_the_range:.2}\tin-a-regular-bin={in_a_regular_bin}\t\
                     overflowed={overflowed}\tfolded={}\toverflow-fraction={:.5}",
                    axis.folded(),
                    axis.overflow_fraction(),
                )
            }
            WhatOneArmFitted::NoHistogram(silence) => format!("no-histogram={silence}"),
        }
    }

    /// Whether the model fit would reject this sample for having left the depth range.
    fn the_fit_would_reject_it(&self) -> bool {
        self.axis()
            .is_some_and(|axis| axis.overflow_fraction() > THE_OVERFLOW_SHARE_THE_FIT_ALLOWS)
    }

    fn print(&self, label: &str) {
        println!(
            "{label}\tbin-scheme\t{}\t{}\t{}",
            self.setting.as_text(),
            self.as_a_row(),
            if self.the_fit_would_reject_it() {
                "the-fit-would-reject-this-sample"
            } else {
                "the-fit-would-take-it"
            },
        );
    }
}

impl WindowRecomputation {
    fn over(
        reference: &ReferenceFasta,
        run_windows: Option<HashMap<GenomePosition, Option<WindowCoverage>>>,
        run_histogram: Option<String>,
        sample: usize,
        covered_positions_per_window: bool,
        bin_scheme: bool,
    ) -> Self {
        Self {
            accumulator: Some(WindowCoverageAccumulator::new(
                the_configuration_a_run_uses(),
            )),
            sweep: covered_positions_per_window.then(FloorSweep::new),
            windows_under_each_floor: None,
            bin_schemes: bin_scheme.then(BinSchemeSweep::new),
            each_bin_scheme: None,
            windows: Vec::new(),
            reference: reference.accessor(),
            reference_contigs: reference.contigs.clone(),
            contig_last_checked: None,
            contigs_checked: Vec::new(),
            bases: Vec::new(),
            reported: Vec::new(),
            positions_observed: 0,
            sample,
            histogram: None,
            histogram_the_run_wrote: run_histogram,
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
        let (tail, histogram) = accumulator.finish();
        self.histogram = Some(histogram);
        self.windows.extend(tail);
        // **Ascending, and asserted rather than assumed**, because the comparison below
        // binary-searches it.
        assert!(
            self.windows.windows(2).all(|pair| pair[0].0 < pair[1].0),
            "the walk's windows are not in ascending centre order",
        );
        if let Some(sweep) = self.sweep.take() {
            let absent = self
                .windows
                .iter()
                .filter(|(_, window)| window.is_absent())
                .count() as u64;
            self.windows_under_each_floor = Some(sweep.close(self.windows.len() as u64, absent));
        }
        if let Some(bin_schemes) = self.bin_schemes.take() {
            // **Closed after `self.histogram` is set**, because the arm configured as the run is
            // checked against exactly that histogram.
            let the_walks = self
                .histogram
                .as_ref()
                .expect("the walk's own histogram is taken before the bin schemes are closed");
            self.each_bin_scheme = Some(bin_schemes.close(the_walks));
        }
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
            // **"Nothing to report" is spelled two ways and means one thing.** The run writes an
            // *absent* pair — two `NaN`s — for a sample whose coverage does not begin at this
            // position and for one whose window the position floor silenced, because a cohort
            // locus collapses the two on its way to the record. This walk answers `None` where it
            // finalised no centre, and an absent pair where the floor silenced one. Comparing the
            // spellings rather than the fact reports a disagreement at every locus a sample does
            // not cover: measured, 277 of 13,866 on the tomato slice, every one of them the run
            // saying "absent" against the walk saying "none".
            let run_said = a_usable_window(run_said);
            let walk_says = a_usable_window(self.window_at(at));
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
                // The run reports a window this walk does not: impossible, since the walk sees
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

impl WindowRecomputation {
    /// This walk's histogram, and — when a run wrote its own down — whether the two are the same
    /// text.
    ///
    /// **Compared as the run's own rendering of it, not field by field.** The rendering is the
    /// library's (`recorded_windows::write_the_histogram`), so the two sides cannot disagree
    /// about what a histogram *is*; and it carries the fitted bin width as a bit pattern and every
    /// cell including the empty ones, so two histograms of different widths cannot agree by both
    /// being sparse.
    fn report_the_histogram(&self, label: &str) {
        let histogram = self
            .histogram
            .as_ref()
            .expect("the pass is closed before its histogram is reported");
        let mut walk_says = String::new();
        write_the_histogram(self.sample, histogram, &mut walk_says);
        match histogram {
            SampleHistogram::Fitted(fitted) => {
                println!(
                    "{label}\thistogram-windows-folded\t{}",
                    fitted.windows_folded
                );
                println!(
                    "{label}\thistogram-windows-under-the-floor\t{}",
                    fitted.windows_under_the_floor
                );
                println!(
                    "{label}\thistogram-depth-bin-width\t{}",
                    fitted.depth_bin_width
                );
            }
            silence => println!("{label}\thistogram\tnone: {silence:?}"),
        }
        let Some(run_said) = &self.histogram_the_run_wrote else {
            return;
        };
        if run_said.trim_end() == walk_says.trim_end() {
            println!("{label}\thistogram-agrees-with-the-run\tyes");
        } else {
            println!("{label}\thistogram-agrees-with-the-run\tno");
            println!(
                "{label}\thistogram-differs-at\t{}",
                where_the_two_lines_first_differ(run_said, &walk_says),
            );
            panic!(
                "the run's histogram for {label} is not the one a straight walk of the same \
                 store makes",
            );
        }
    }
}

/// What of two histogram lines first disagrees, **named** — a header field by its name, a cell by
/// the GC bin and depth bin it holds.
///
/// **A whole line is not a failure message.** At the shipped bin counts a fitted line is eight
/// header fields and 20,050 cells, so printing its head reports a cell difference as two
/// identical-looking truncations, and the header fields are the ones most likely to agree.
fn where_the_two_lines_first_differ(run: &str, walk: &str) -> String {
    let run: Vec<&str> = run.trim_end().split('\t').collect();
    let walk: Vec<&str> = walk.trim_end().split('\t').collect();
    // The depth bins the *run* declared: what makes a cell's index into a pair of bins. Its own
    // field is compared below, so a disagreement about it is reported before any cell is.
    let depth_bins: usize = run.get(5).and_then(|field| field.parse().ok()).unwrap_or(0);
    for (field, (run, walk)) in run.iter().zip(&walk).enumerate() {
        if run != walk {
            return format!(
                "{}\trun={run}\twalk={walk}",
                one_field_named(field, depth_bins)
            );
        }
    }
    format!(
        "no field differs in the {} they share, but the run's line has {} fields and this walk's \
         {}",
        run.len().min(walk.len()),
        run.len(),
        walk.len(),
    )
}

/// What field `field` of a histogram line holds, in the terms the histogram is written in.
///
/// The eight header fields are named; from there each field is one cell, row-major over
/// `[gc_bin][depth_bin]` with one overflow bin after each row's regular ones — so `depth_bins + 1`
/// cells to a GC row.
fn one_field_named(field: usize, depth_bins: usize) -> String {
    const HEADER: [&str; 8] = [
        "the sample",
        "the word `fitted`",
        "the window width",
        "the GC bin count",
        "the fitted depth bin width, as bits",
        "the depth bin count",
        "the windows folded",
        "the windows the floor silenced",
    ];
    if let Some(named) = HEADER.get(field) {
        return (*named).to_owned();
    }
    let cell = field - HEADER.len();
    let row = depth_bins + 1;
    if row == 0 {
        return format!("cell {cell}, whose bins the run's own depth bin count does not describe");
    }
    format!("the cell at GC bin {} depth bin {}", cell / row, cell % row,)
}

/// A window that reports something, or `None` for either way of having nothing to report.
///
/// **Absent and missing are one answer to the consumer**, which is the filter: it skips a sample
/// it has no window for, and an absent pair is what it is handed for one. Keeping them apart here
/// would make this oracle fail on a difference in spelling.
///
/// **One thing this hides, and it is not reachable today.** A window the *run* silenced against
/// one the walk measured would land in "the run had no window here" rather than in the
/// disagreements — which is where a truncated cover's centres land too, so the two would be
/// counted together. Nothing on the tomato slice silences a window (`windows-under-the-floor` is
/// 0 there), so the case does not arise; plan step D1, which sets the floor from a distribution,
/// is where it could start to.
fn a_usable_window(window: Option<WindowCoverage>) -> Option<WindowCoverage> {
    window.filter(|window| !window.is_absent())
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
            sweep,
            bin_schemes,
            bases,
            positions_observed,
            windows,
            reference: _,
            reference_contigs: _,
            contig_last_checked: _,
            contigs_checked: _,
            reported: _,
            comparison: _,
            sample: _,
            histogram: _,
            histogram_the_run_wrote: _,
            // Both filled by `close` out of the two sweeps above, and neither fed a record.
            windows_under_each_floor: _,
            each_bin_scheme: _,
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
            // **The same four values, in the same order, into the sweep** — which is what makes
            // its arms answer about this walk's windows and not about some other stream.
            if let Some(sweep) = sweep.as_mut() {
                sweep.observe(at.contig, at.position, base, depth);
            }
            if let Some(bin_schemes) = bin_schemes.as_mut() {
                bin_schemes.observe(at.contig, at.position, base, depth);
            }
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
        if let Some(under_each_floor) = &self.windows_under_each_floor {
            under_each_floor.print(label);
        }
        for answer in self.each_bin_scheme.iter().flatten() {
            answer.print(label);
        }
        self.report_the_histogram(label);
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
    /// Also report how many covered positions each window was built from, as the share of windows
    /// a floor of each candidate size would silence.
    covered_positions_per_window: bool,
    /// Also report, per sample, where the histogram's depth axis is cut and how much of the
    /// sample falls off the end of it, at the run's setting and at each candidate.
    bin_scheme: bool,
}

const USAGE: &str = "usage: ng_window_coverage_probe [--reference <reference.fa>] \
                     [--windows-from-the-run <windows.tsv>] \
                     [--covered-positions-per-window] [--bin-scheme] \
                     <a store.psp> [more stores...]";

impl Arguments {
    fn from_command_line() -> Self {
        let mut stores = Vec::new();
        let mut reference = None;
        let mut windows_from_the_run = None;
        let mut covered_positions_per_window = false;
        let mut bin_scheme = false;
        let mut argv = std::env::args().skip(1);
        while let Some(argument) = argv.next() {
            match argument.as_str() {
                "--covered-positions-per-window" => covered_positions_per_window = true,
                "--bin-scheme" => bin_scheme = true,
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
        assert!(
            !covered_positions_per_window || reference.is_some(),
            "--covered-positions-per-window counts the positions in each window, which are the \
             recomputation's, and that needs --reference",
        );
        assert!(
            !bin_scheme || reference.is_some(),
            "--bin-scheme fits a depth axis from the windows the recomputation builds, and that \
             needs --reference",
        );
        Self {
            stores,
            reference,
            windows_from_the_run,
            covered_positions_per_window,
            bin_scheme,
        }
    }
}

/// The windows a calling run recorded, split by sample — one map per sample index, in the run's
/// own sample order, which is the order the stores are listed in.
///
/// **The row format is the library's**, written and read by one pair of functions in
/// [`recorded_windows`](pop_var_caller::run::cohort_merge::recorded_windows), because a column
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

/// The histogram line a run wrote for each sample, in the run's own sample order.
///
/// **One line a sample, and that is checked**: a file with a line missing would otherwise shift
/// every later sample's histogram onto its neighbour, which is the same failure the row file's
/// sample-index check exists for.
fn histograms_the_run_recorded(recorded: &Path, samples: usize) -> Vec<String> {
    let path = histograms_beside(recorded);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|why| {
        panic!(
            "reading {}: {why} — a run writes it beside its rows, so a missing file means the \
             rows came from a run that did not finish",
            path.display()
        )
    });
    let lines: Vec<String> = text.lines().map(str::to_owned).collect();
    assert_eq!(
        lines.len(),
        samples,
        "{} holds {} histogram lines and {samples} stores were given",
        path.display(),
        lines.len(),
    );
    lines
}

fn main() {
    let arguments = Arguments::from_command_line();
    let stores = &arguments.stores;
    let reference = arguments.reference.as_deref().map(ReferenceFasta::open);
    let mut run_windows = arguments
        .windows_from_the_run
        .as_deref()
        .map(|recorded| windows_the_run_recorded(recorded, stores.len()));
    // **Beside the rows, from the file the run wrote its histograms to** — the same path with
    // `.histograms` after it, which is where `recorded_windows` puts them.
    let run_histograms = arguments
        .windows_from_the_run
        .as_deref()
        .map(|recorded| histograms_the_run_recorded(recorded, stores.len()));

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
                run_histograms
                    .as_ref()
                    .map(|per_sample| per_sample[sample].clone()),
                sample,
                arguments.covered_positions_per_window,
                arguments.bin_scheme,
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
    /// The stores' floor sweeps summed, `None` where the sweep was not asked for. **Every store
    /// has one or none**: the flag is the run's, not the store's, so a total built from some of
    /// them would be a distribution over a subset nothing names.
    windows_under_each_floor: Option<WindowsUnderEachFloor>,
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
        if let Some(this_store) = &store.windows_under_each_floor {
            // **The first store's rows are cloned and later stores are added into them**, so that
            // the total's floors are the same floors in the same order whatever the store count,
            // and so that no second constructor has to reproduce that list.
            match &mut self.windows_under_each_floor {
                Some(total) => total.add(this_store),
                nothing_yet => *nothing_yet = Some(this_store.clone()),
            }
        }
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
        if let Some(under_each_floor) = &self.windows_under_each_floor {
            under_each_floor.print(label);
        }
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

    /// **The three spellings of "no usable window" collapse into one**, which is what stops this
    /// oracle reporting a disagreement at every locus a sample does not cover: the run records an
    /// absent pair there, and a walk that finalised no centre has nothing at all.
    #[test]
    fn absent_and_missing_are_one_answer_and_a_number_is_not() {
        assert_eq!(a_usable_window(None), None);
        assert_eq!(
            a_usable_window(Some(WindowCoverage::absent())),
            None,
            "an absent pair is a window with nothing to report, like no window at all",
        );
        let measured = WindowCoverage {
            gc_fraction: 0.5,
            mean_depth: 3.25,
        };
        assert_eq!(
            a_usable_window(Some(measured)),
            Some(measured),
            "a window that reports something is not collapsed with the ones that do not",
        );
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

    /// Feed `positions` consecutive covered positions from base 1 of one contig to a floor sweep
    /// and to a plain accumulator configured as a run's, exactly as the walk feeds both, and
    /// close the sweep against what the plain one said.
    ///
    /// Every base is `A` and every depth 1, so nothing here turns on GC or on depth: what the
    /// sweep counts is how many positions fell inside each window.
    fn a_sweep_over_consecutive_positions(positions: u64) -> WindowsUnderEachFloor {
        let mut sweep = FloorSweep::new();
        let mut beside_it = WindowCoverageAccumulator::new(the_configuration_a_run_uses());
        let mut windows = Vec::new();
        for position in 1..=positions {
            sweep.observe(ContigId(0), Position(position), b'A', 1);
            beside_it.observe(ContigId(0), Position(position), b'A', 1);
            while let Some(window) = beside_it.pop_ready() {
                windows.push(window);
            }
        }
        let (tail, _) = beside_it.finish();
        windows.extend(tail);
        let absent = windows
            .iter()
            .filter(|(_, window)| window.is_absent())
            .count() as u64;
        sweep.close(windows.len() as u64, absent)
    }

    /// Windows silenced at one floor of a closed sweep's answer.
    fn silenced_at(under_each_floor: &WindowsUnderEachFloor, floor: u32) -> u64 {
        under_each_floor
            .floors
            .iter()
            .find(|at| at.floor == floor)
            .unwrap_or_else(|| panic!("the sweep has no arm at a floor of {floor}"))
            .silenced
    }

    /// Ten covered positions within 250 bases of each other, so each is a centre and each of the
    /// ten windows holds all ten: a floor of ten silences none and a floor of fifteen silences
    /// every one.
    ///
    /// The claim under test is that an arm's silenced count is the number of windows holding fewer
    /// positions than its floor, checked where both numbers are countable by hand.
    ///
    /// **Ten positions never reach any centre's right edge**, so every window here comes back from
    /// `finish` and none is popped mid-stream: the drain in [`FloorSweep::observe`] is exercised by
    /// `the_sweep_traces_a_distribution_whose_windows_hold_different_counts`, whose 600 positions
    /// close 350 centres before the stream ends, and by nothing else here.
    #[test]
    fn an_arms_silenced_count_is_the_windows_holding_fewer_positions_than_its_floor() {
        let sweep = a_sweep_over_consecutive_positions(10);
        assert_eq!(sweep.floors[0].windows, 10, "ten positions, ten centres");
        for floor in [1, 2, 3, 5, 8, 10] {
            assert_eq!(
                silenced_at(&sweep, floor),
                0,
                "every window holds all ten positions, so a floor of {floor} silences none",
            );
        }
        for floor in [15, 20, 501] {
            assert_eq!(
                silenced_at(&sweep, floor),
                10,
                "every window holds ten positions, so a floor of {floor} silences all ten windows",
            );
        }
    }

    /// Six hundred consecutive covered positions, where the counts differ between windows and are
    /// arithmetic: the window centred at `p` holds `[p − 250, p + 250]` clipped to `1..=600`, so
    /// the centres in the middle hold the full 501 and those at either end hold as few as 251.
    ///
    /// The three floors checked are the three shapes the distribution has: below every count,
    /// through the middle of it, and above every count.
    #[test]
    fn the_sweep_traces_a_distribution_whose_windows_hold_different_counts() {
        let sweep = a_sweep_over_consecutive_positions(600);
        assert_eq!(sweep.floors[0].windows, 600);
        assert_eq!(
            silenced_at(&sweep, 250),
            0,
            "the thinnest window here is the one centred at base 1, holding bases 1 to 251",
        );
        // A centre at `p <= 250` holds `p + 250`, which is under 300 for `p` in 1..=49; a centre
        // at `p >= 351` holds `851 - p`, under 300 for `p` in 552..=600. Forty-nine at each end.
        assert_eq!(silenced_at(&sweep, 300), 98);
        // The hundred centres at 251..=350 are the only ones holding the full 501, and 501 is the
        // highest floor the configuration allows.
        assert_eq!(silenced_at(&sweep, 501), 500);
    }

    /// Two identical stores summed: every count doubles — 1,200 windows against 600, 196 silenced
    /// at a floor of 300 against 98, 1,000 at 501 against 500 — and the total's floors are the
    /// same floors in the same order, which is what `main` prints the total block from.
    #[test]
    fn two_sweeps_sum_floor_by_floor() {
        let one = a_sweep_over_consecutive_positions(600);
        let other = a_sweep_over_consecutive_positions(600);
        let mut total = one.clone();
        total.add(&other);
        assert_eq!(
            total.floors.iter().map(|at| at.floor).collect::<Vec<_>>(),
            the_floors_asked_about(),
        );
        assert_eq!(total.floors[0].windows, 1_200);
        assert_eq!(silenced_at(&total, 300), 196);
        assert_eq!(silenced_at(&total, 501), 1_000);
    }

    /// A sweep fed different positions from the walk beside it is caught, rather than reporting a
    /// distribution of a stream nothing else saw.
    #[test]
    #[should_panic(expected = "not fed the same positions")]
    fn a_sweep_whose_windows_differ_from_the_walks_is_refused() {
        let mut sweep = FloorSweep::new();
        for position in 1..=10u64 {
            sweep.observe(ContigId(0), Position(position), b'A', 1);
        }
        sweep.close(11, 0);
    }

    /// The arm at the shipped floor must agree with the walk's own count of absent windows —
    /// the one place the sweep and the recomputation answer the same question.
    #[test]
    #[should_panic(expected = "the arm at the shipped floor")]
    fn a_sweep_disagreeing_with_the_walk_at_the_shipped_floor_is_refused() {
        let mut sweep = FloorSweep::new();
        for position in 1..=10u64 {
            sweep.observe(ContigId(0), Position(position), b'A', 1);
        }
        // Ten positions in one window clear the shipped floor of 50 nowhere, so the walk beside
        // this sweep would have called all ten absent; claiming none were is the mismatch.
        sweep.close(10, 0);
    }

    /// A walk that finalised nothing says nothing, and is refused rather than printing a row of
    /// zeroes at every floor, which reads as "no window was ever too thin".
    #[test]
    #[should_panic(expected = "says nothing about how many covered positions")]
    fn a_sweep_over_an_empty_store_is_refused() {
        FloorSweep::new().close(0, 0);
    }

    /// A sweep whose rows are stated outright rather than counted from a stream: one arm a
    /// `(floor, windows, silenced)` triple, its accumulator fed nothing so `finish` adds an empty
    /// tail.
    ///
    /// The two checks below are about counts no stream can produce, so a stream cannot set them
    /// up.
    fn a_sweep_saying(rows: &[(u32, u64, u64)]) -> FloorSweep {
        FloorSweep {
            arms: rows
                .iter()
                .map(|&(floor, windows, silenced)| OneFloor {
                    accumulator: WindowCoverageAccumulator::new(WindowCoverageConfig {
                        min_window_positions: floor,
                        ..the_configuration_a_run_uses()
                    }),
                    counted: WindowsUnderOneFloor {
                        floor,
                        windows,
                        silenced,
                    },
                })
                .collect(),
        }
    }

    /// A centre always lies in its own window, so no window holds fewer than one covered position
    /// and a floor of 1 can silence nothing.
    ///
    /// **It is the only check in `close` that names a number rather than comparing two quantities
    /// that can both be zero**, which is what makes it the one that cannot be vacuous on a store
    /// whose windows are all thick.
    #[test]
    #[should_panic(expected = "the arm at a floor of 1 silenced 1 windows")]
    fn an_arm_at_a_floor_of_one_silencing_a_window_is_refused() {
        a_sweep_saying(&[(1, 5, 1), (window_coverage::MIN_WINDOW_POSITIONS, 5, 1)]).close(5, 1);
    }

    /// The two floor-keyed checks are written as "if this is the row at floor `F`", so a floor
    /// list that no longer holds `F` would run neither and fail nothing. `close` refuses that
    /// instead.
    #[test]
    #[should_panic(expected = "no arm at a floor of 1")]
    fn a_sweep_that_lost_the_floor_a_check_is_written_against_is_refused() {
        a_sweep_saying(&[(window_coverage::MIN_WINDOW_POSITIONS, 5, 0)]).close(5, 0);
    }

    /// Feed a bin-scheme sweep, and a plain accumulator configured as a run's, one contig per
    /// `(depth, positions)` segment — consecutive positions from base 1, every one at that
    /// segment's depth — and close the sweep against the histogram the plain one finished with.
    ///
    /// **A contig apiece, so a window never straddles two depths**: the accumulator never spans
    /// contigs, so every window of segment *i* has mean depth exactly *i*'s, and a test can say
    /// how many windows sit where without allowing for a boundary.
    fn bin_schemes_over(segments: &[(u32, u64)]) -> Vec<OneBinSchemeAnswer> {
        let mut sweep = BinSchemeSweep::new();
        let mut beside_it = WindowCoverageAccumulator::new(the_configuration_a_run_uses());
        for (contig, &(depth, positions)) in segments.iter().enumerate() {
            let contig = ContigId(contig as u32);
            for position in 1..=positions {
                sweep.observe(contig, Position(position), b'A', depth);
                beside_it.observe(contig, Position(position), b'A', depth);
                while beside_it.pop_ready().is_some() {}
            }
        }
        let (_, histogram) = beside_it.finish();
        sweep.close(&histogram)
    }

    fn the_arm(answers: &[OneBinSchemeAnswer], setting: BinSchemeSetting) -> &OneBinSchemeAnswer {
        answers
            .iter()
            .find(|answer| answer.setting == setting)
            .unwrap_or_else(|| panic!("the sweep has no arm for {setting:?}"))
    }

    fn the_axis(answers: &[OneBinSchemeAnswer], setting: BinSchemeSetting) -> FittedAxis {
        the_arm(answers, setting)
            .axis()
            .unwrap_or_else(|| panic!("the arm at {setting:?} fitted no axis"))
    }

    /// Six hundred covered positions all at 8 reads: every window's mean depth is 8, so the
    /// fitted median is 8, the bin width is `8 × 10 / 400 = 0.2`, the regular bins stop at 80,
    /// and nothing overflows.
    ///
    /// **The arithmetic is the whole of what this pins** — the width the accumulator fits and the
    /// depth this file recovers it from are two directions of one formula, and a reader can check
    /// both against a stream whose depth never varies.
    #[test]
    fn a_stream_of_one_depth_fits_an_axis_this_file_can_recover_that_depth_from() {
        let answers = bin_schemes_over(&[(8, 600)]);
        let axis = the_axis(&answers, BinSchemeSetting::AsTheRunFitsIt);
        assert!(
            (axis.median_depth - 8.0).abs() < 1e-9,
            "{}",
            axis.median_depth
        );
        assert!(
            (axis.depth_bin_width - 0.2).abs() < 1e-9,
            "{}",
            axis.depth_bin_width
        );
        assert!(
            (axis.top_of_the_range - 80.0).abs() < 1e-9,
            "{}",
            axis.top_of_the_range
        );
        assert_eq!(axis.folded(), 600);
        assert_eq!(
            axis.overflowed, 0,
            "every window sits at 8, a tenth of the way up a range that reaches 80",
        );
    }

    /// **The range decides whether the model fit takes the sample at all**, on a stream where the
    /// windows are not all at one depth.
    ///
    /// 700 windows at 1 read a position and 300 at 20, on two contigs so neither mixes. The
    /// median is 1, so the range's multiple *is* the top of the axis: at the shipped ten it
    /// stops at 10 and the 300 deep windows overflow — 3 windows in 10, past the fit's guard of
    /// 1 in 5, so the sample is rejected and its coverage goes unused. At forty it stops at 40,
    /// nothing overflows, and the fit takes the sample; at 2.5 it stops at 2.5 and the same 300
    /// overflow.
    #[test]
    fn the_range_decides_whether_the_deep_windows_overflow_and_the_fit_refuses_the_sample() {
        let answers = bin_schemes_over(&[(1, 700), (20, 300)]);
        for (setting, overflowed, rejected) in [
            (BinSchemeSetting::RangeOf(2.5), 300, true),
            (BinSchemeSetting::AsTheRunFitsIt, 300, true),
            (BinSchemeSetting::RangeOf(40.0), 0, false),
        ] {
            let axis = the_axis(&answers, setting);
            assert!(
                (axis.median_depth - 1.0).abs() < 1e-9,
                "700 of the 1,000 windows sit at 1, so the median is 1 whatever the range: {}",
                axis.median_depth,
            );
            assert_eq!(axis.folded(), 1_000);
            assert_eq!(
                axis.overflowed, overflowed,
                "at {setting:?} the axis stops at {} and 300 windows sit at 20",
                axis.top_of_the_range,
            );
            assert_eq!(
                the_arm(&answers, setting).the_fit_would_reject_it(),
                rejected,
                "at {setting:?} the overflow fraction is {:.2} against a guard of \
                 {THE_OVERFLOW_SHARE_THE_FIT_ALLOWS}",
                axis.overflow_fraction(),
            );
        }
    }

    /// The arm configured exactly as the run must reproduce the walk's own histogram — the one
    /// place this sweep and the recomputation answer the same question, and the only check here
    /// that compares against something computed outside the sweep.
    #[test]
    #[should_panic(expected = "the arm configured as the run is fitted came out")]
    fn a_bin_scheme_sweep_disagreeing_with_the_walk_is_refused() {
        let mut sweep = BinSchemeSweep::new();
        for position in 1..=600u64 {
            sweep.observe(ContigId(0), Position(position), b'A', 8);
        }
        // The walk it is closed against saw twice that depth, so no arm can match it.
        // A walk that saw twice the depth: its median is 16 and its bins twice as wide, so the
        // arm running the run's own configuration cannot match it.
        let mut a_different_walk = WindowCoverageAccumulator::new(the_configuration_a_run_uses());
        for position in 1..=600u64 {
            a_different_walk.observe(ContigId(0), Position(position), b'A', 16);
            while a_different_walk.pop_ready().is_some() {}
        }
        let (_, histogram) = a_different_walk.finish();
        sweep.close(&histogram);
    }

    /// A sweep over a store that folded nothing says nothing about where the axis should be cut,
    /// and is refused rather than printing an arm's worth of zeroes at every setting.
    #[test]
    #[should_panic(expected = "no arm folded a window")]
    fn a_bin_scheme_sweep_over_an_empty_store_is_refused() {
        let sweep = BinSchemeSweep::new();
        let (_, histogram) =
            WindowCoverageAccumulator::new(the_configuration_a_run_uses()).finish();
        sweep.close(&histogram);
    }

    /// **The fraction's denominator is the windows in the cells, and the guard is a strict `>`**
    /// — the two places this file could disagree with the fit it measures against while every
    /// other test here stayed green.
    ///
    /// 800 windows at 1 read a position, 200 at 20, and 30 more at 1 on a contig too short to
    /// build a window the floor will take. The median is 1, so the shipped range stops the axis
    /// at 10 and the 200 deep windows overflow: **200 of the 1,000 windows in the cells, exactly
    /// the fit's guard of 1 in 5** — and the fit rejects a sample only when the share is *past*
    /// its guard, so this sample is taken.
    ///
    /// **The 30 silenced windows are what pins the denominator.** They are finalised and counted
    /// under the floor, but they reach no cell, so the fit never sees them: counting them here
    /// would put the share at 194 in 1,000 and would have taken this sample for the wrong reason.
    #[test]
    fn the_overflow_share_is_over_the_windows_in_the_cells_and_the_guard_is_not_reached_at_it() {
        let answers = bin_schemes_over(&[(1, 800), (20, 200), (1, 30)]);
        let arm = the_arm(&answers, BinSchemeSetting::AsTheRunFitsIt);
        let axis = arm.axis().expect("the axis is fitted");
        assert_eq!(
            axis.folded(),
            1_000,
            "1,000 windows reached a cell; the 30 on the short contig were silenced by the floor \
             and reached none",
        );
        assert_eq!(axis.in_a_regular_bin, 800);
        assert_eq!(axis.overflowed, 200);
        assert!(
            (axis.overflow_fraction() - THE_OVERFLOW_SHARE_THE_FIT_ALLOWS).abs() < f64::EPSILON,
            "200 of 1,000 is the guard exactly: {}",
            axis.overflow_fraction(),
        );
        assert!(
            !arm.the_fit_would_reject_it(),
            "the fit rejects a sample past its guard, not at it",
        );
    }

    /// **What the scale sample buys, on a stream whose depth changes after it has closed.**
    ///
    /// 10,300 windows at 1 read a position, then 19,700 at 100, on two contigs. The run's setting
    /// fits from the first 10,000 windows, every one of which sits at 1, so its axis stops at 10
    /// and the 19,700 deep windows all overflow — 66 windows in 100, and the model fit refuses
    /// the sample. The arm that holds every window back until the pass ends sees a median of 100,
    /// cuts an axis reaching 1,000, and overflows nothing.
    ///
    /// This is the failure the constant exists to trade against: a scale sample too small to
    /// reach the sample's own depth costs the whole sample, not a little resolution.
    #[test]
    fn a_scale_sample_that_closes_before_the_depth_changes_costs_the_sample() {
        let answers = bin_schemes_over(&[(1, 10_300), (100, 19_700)]);
        let as_the_run = the_axis(&answers, BinSchemeSetting::AsTheRunFitsIt);
        assert!(
            (as_the_run.median_depth - 1.0).abs() < 1e-9,
            "the first 10,000 windows all sit at 1: {}",
            as_the_run.median_depth,
        );
        assert_eq!(as_the_run.overflowed, 19_700);
        assert!(the_arm(&answers, BinSchemeSetting::AsTheRunFitsIt).the_fit_would_reject_it());

        let every_window = the_axis(&answers, BinSchemeSetting::ScaleSampleOf(u32::MAX));
        assert!(
            (every_window.median_depth - 100.0).abs() < 1e-9,
            "19,700 of the 30,000 windows sit at 100: {}",
            every_window.median_depth,
        );
        assert_eq!(every_window.overflowed, 0);
        assert!(
            !the_arm(&answers, BinSchemeSetting::ScaleSampleOf(u32::MAX)).the_fit_would_reject_it()
        );
    }

    /// On a stream shorter than the shipped scale sample, the arm that holds every window back
    /// answers exactly as the run does — both fit at `finish` from every window they have.
    ///
    /// It is the control for the test above: without it, that one's difference could as easily be
    /// a mis-wired arm as the scale sample doing its job.
    #[test]
    fn an_arm_that_holds_every_window_back_answers_as_the_run_does_on_a_short_stream() {
        let answers = bin_schemes_over(&[(8, 600)]);
        assert_eq!(
            the_arm(&answers, BinSchemeSetting::ScaleSampleOf(u32::MAX)).as_a_row(),
            the_arm(&answers, BinSchemeSetting::AsTheRunFitsIt).as_a_row(),
            "600 windows is under the shipped scale sample of 10,000, so both fit at `finish` \
             from every window they have",
        );
    }
}
