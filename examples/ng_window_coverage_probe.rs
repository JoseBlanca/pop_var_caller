//! **Does a one-base record's head really say what its evidence says?** — plan step B2's
//! measurement, and the claim spec `window_coverage.md` §3.1 rests its cheapest rule on.
//!
//! ```text
//! cargo run --release --example ng_window_coverage_probe -- <a store.psp> [more stores...]
//! ```
//!
//! # The question
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
//! passed; it knows nothing about what they count. Plan step C3 adds a second measurement over
//! the same pass — every position's window coverage, recomputed from the whole store — and that
//! is a new type implementing the trait plus one more entry in `main`'s slice, with nothing in
//! the walk or in [`SingleBaseEquality`] disturbed.

use std::path::{Path, PathBuf};

use pop_var_caller::ng::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use pop_var_caller::ng::psp::{PspReader, RecordHead};

#[cfg(test)]
use pop_var_caller::ng::locus_generation::{LocusLen, SequenceObservation, SsrDetail};
#[cfg(test)]
use pop_var_caller::ng::types::{
    ContigId, GenomeRegion, Motif, Position, ReadGroupId, SummedLogError,
};

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

fn main() {
    let stores: Vec<PathBuf> = std::env::args().skip(1).map(PathBuf::from).collect();
    assert!(
        !stores.is_empty(),
        "usage: ng_window_coverage_probe <a store.psp> [more stores...]",
    );

    println!("phase\tng-window-coverage-probe");
    let mut total = SingleBaseEquality::default();
    let mut stores_walked = 0usize;
    for store in &stores {
        let mut over_this_store = SingleBaseEquality::default();
        walk(store, &mut [&mut over_this_store]);
        let label = store.file_name().map_or_else(
            || store.display().to_string(),
            |name| name.to_string_lossy().into_owned(),
        );
        // **The path, beside the basename that labels every row.** Two stores can share a
        // basename, and a figure quoted from this output has to say which file it came from.
        println!("{label}\tstore\t{}", store.display());
        over_this_store.print(&label);
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
}
