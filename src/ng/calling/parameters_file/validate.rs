//! **Refusing a file that parses and means nothing** — the constraints no shape can state.
//!
//! Step C1's reader answers one question: is this text TOML that spells this shape? A file can
//! pass that and still say something no run can mean — an inbreeding coefficient of 1.7, a length
//! spectrum whose shares sum to 1.4, a contamination table in which nobody was measured. **Owner's
//! decision of 2026-08-28**, because no step of the plan owned it. Spec §9's own sentence is about
//! the other half — "a **malformed** file fails at read **with a line number**, which is what using
//! an existing parser buys" — and neither clause fits here: these files are well formed, and this
//! walk has no line to give (below).
//!
//! # Where this runs, and why it is not the reader
//!
//! After parsing and before the projection back to `RunParameters`. It is deliberately **not**
//! folded into [`ParametersFile::from_toml`](super::ParametersFile::from_toml): that entry point
//! answers a question about the *text*, and a caller that wants to read a file and inspect it —
//! this module's own tests among them — should not have to satisfy the run's constraints to do so.
//!
//! **⚑ Nothing in a run calls this yet.** The projection back to `RunParameters` is the caller
//! that must run it first, and that is the other half of step C2; until it lands, reading a
//! parameters file through this module's public entry point does not validate it.
//!
//! # What it cannot give, and says so rather than pretending
//!
//! **No line number.** The refusals of C1 come from a parser holding a byte span; these come from
//! walking a value that no longer remembers where it was written. What is given instead is **the
//! key's full path in the file's own spelling** — `inbreeding.by_sample["TS-1"]` — which a reader
//! can search for, and which survives the file being reformatted.
//!
//! **A line is available and was not taken.** The `toml` crate re-exports `serde_spanned::Spanned`,
//! so a field wrapped in `Spanned<T>` carries a byte range out of the parse this module's caller
//! already ran — no second pass. It was rejected because the wrapper would have to sit on the
//! *shape*, where it reaches `Serialize`, `PartialEq` and the round-trip equality the whole design
//! rests on (spec §1.2 goal 1), to improve a message on a file that is well formed. If the paths
//! below turn out not to be enough, that is the route.
//!
//! # Three kinds of refusal, and the third is the one that costs a day
//!
//! - **A value outside its documented range.** An inbreeding coefficient is a fraction in `[0, 1)`,
//!   a curve weight a share in `[0, 1]`, a substitution rate a probability. Also, for every float
//!   this walk reaches, a value that is not finite — which no range contains and which propagates
//!   silently through every score it touches.
//! - **A shape that is internally impossible.** A length spectrum with an even number of shares has
//!   no middle entry to be the reference length; one with fewer than three cannot express a slip in
//!   either direction; one whose shares do not sum to one is not a distribution.
//! - **A shape that says the opposite of what it means.** These are spec §5's states written
//!   longhand, and they are the expensive ones because they parse, project, and produce a plausible
//!   VCF. A contamination table in which no row was measured is an *uncontaminated run* written as
//!   a mixture with every fraction zero — the read likelihood then takes its mixture path instead
//!   of its plain one. A measurement missing either of its evidence counts is the in-memory *not
//!   measured* shape, which this file expresses by having no `measurement` key at all: it projects
//!   to the same value the absent key gives, so nothing downstream computes differently — what
//!   goes wrong is that the file claims a measurement, carries a `fitted_from_reads_of` true of
//!   nothing, and the run reports a lane as uncorrected while the file says it was measured.
//!
//! # And one refusal that is about this writer rather than about runs
//!
//! **An evidence count of exactly `i64::MAX` is the writer's saturation marker, not a
//! measurement.** TOML's integers are signed and 64-bit, so a `u64` count above `2^63 − 1` has no
//! representation; the writer saturates rather than emitting a number three different readers gave
//! three different answers to (step B2). A count arriving back at exactly that value is therefore
//! either a saturated write or a number nobody could have counted — 9.2 quintillion reads — and
//! reading it as evidence would put a number the writer knows it lost into a run's report.

use super::{
    ContaminationMeasurement, EvidenceCount, LevelSmoothing, ParametersFile, ParametersFileError,
    ShareCurve, ShareSmoothing, SlippageCurve, Warrant, WarrantedValue,
};
use crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT;
use std::collections::{BTreeMap, BTreeSet};

/// How far a length spectrum's shares may sum from one and still be a distribution.
///
/// **Wide enough for the rounding a normalised vector carries, narrow enough to catch a hand
/// edit.** The miss is in the *summation*, not in the round trip: a float written by this writer
/// and read back is the same double, but a spectrum the fit normalised is a division followed by
/// an addition of many terms, and those do not land on one.
///
/// **Measured, over geometric spectra normalised to one:** 3, 5, 9 and 41 offsets summed to
/// exactly one; 21 offsets missed by 1.1e-16 and 101 offsets by 6.7e-16. So the worst seen is
/// about seven parts in ten quadrillion, and this tolerance is 1.5 million times that — while the
/// smallest edit a person would make to one share, a hundredth, is ten million times larger than
/// the tolerance. **Thirteen orders of magnitude separate the two**, and anything from about
/// 1e-14 to 1e-4 would divide them equally well; this is not a tuned number and nothing here
/// depends on its exact value.
const SHARES_MAY_MISS_ONE_BY: f64 = 1e-9;

/// The largest integer TOML can hold, which is what the writer saturates a `u64` count to.
const A_SATURATED_COUNT: u64 = i64::MAX as u64;

impl ParametersFile {
    /// **Refuse a file that parses and says something no run can mean.**
    ///
    /// Runs after [`from_toml`](Self::from_toml) and before the projection back to
    /// `RunParameters`. See this module's header for the three kinds of refusal and for why no
    /// line number is given.
    ///
    /// # Errors
    ///
    /// [`ParametersFileError::Meaningless`], naming the key's full path in the file's own
    /// spelling and what is wrong with it. The first failure stops the walk: a file with two bad
    /// values is reported one at a time, because a reader fixing them will re-run anyway and a
    /// list of every fault in a hand-edited file is mostly one fault's consequences.
    pub fn validate(&self) -> Result<(), ParametersFileError> {
        self.the_version_is_one_this_build_reads()?;
        self.the_identity_names_a_cohort()?;
        self.every_axis_covers_what_it_is_keyed_by()?;
        self.the_batching_puts_each_sample_in_one_batch()?;
        self.every_calibration_is_a_multiplier()?;
        self.contamination_is_absent_or_measured()?;
        self.every_inbreeding_coefficient_is_a_fraction()?;
        self.the_prior_seed_is_a_pair_of_concentrations()?;
        self.the_repeat_tract_numbers_are_what_they_claim()?;
        self.every_stated_constant_is_in_range()?;
        Ok(())
    }

    fn the_version_is_one_this_build_reads(&self) -> Result<(), ParametersFileError> {
        // **Both a zero and a future version are refused**, which is the shape the project's other
        // TOML artefact already uses (`SampleSummaryError::UnsupportedVersion`). Spec §11 defers
        // what a reader *does* with an older file until there is one; refusing a version this
        // build cannot know the meaning of is not that policy, it is the guard that makes the
        // policy possible later.
        if self.format_version == 0 || self.format_version > super::FORMAT_VERSION {
            return Err(refuse(
                "format_version",
                format!(
                    "is {}, and this build reads versions 1 to {}; a higher number means the file was written by a newer build of this caller, so run that build or refit",
                    self.format_version,
                    super::FORMAT_VERSION
                ),
            ));
        }
        Ok(())
    }

    fn the_identity_names_a_cohort(&self) -> Result<(), ParametersFileError> {
        if self.ploidy == 0 {
            return Err(refuse("ploidy", "is 0, and a run calls at least one copy"));
        }
        if self.fitted_from.samples.is_empty() {
            return Err(refuse(
                "fitted_from.samples",
                "is empty, and every per-sample number in the file is keyed by a name in it",
            ));
        }
        if self.fitted_from.read_groups.is_empty() {
            return Err(refuse(
                "fitted_from.read_groups",
                "is empty, and every read of a run belongs to a read group",
            ));
        }

        // **The axis has to be dense over `0..n`.** `RunParameters` indexes calibration and
        // contamination by read-group id, so a gap drops the highest read group entirely and
        // surfaces as a panic at whichever locus first carries one of that library's reads
        // (`run_parameters.rs`, module documentation). Checked here so the message is about the
        // file. **Whether these ids are the *run's* ids is a different question** — that is one of
        // spec §6's bindings and belongs to step D2; this only asks whether the file is coherent
        // with itself.
        let mut seen: Vec<u32> = self
            .fitted_from
            .read_groups
            .iter()
            .map(|row| row.read_group)
            .collect();
        seen.sort_unstable();
        for (expected, found) in seen.iter().enumerate() {
            if *found as usize != expected {
                return Err(refuse(
                    "fitted_from.read_groups",
                    format!(
                        "has no read group {expected}: it names {}, and a run's read groups are \
                         numbered from zero with none missing. A gap silently drops the highest \
                         read group",
                        a_list_of(&seen)
                    ),
                ));
            }
        }

        // **`fitted_from.samples` is the read-group table's own first-seen sample order, exactly
        // once each.** Two separate failures, and both are silent.
        //
        // A **repeat** in the list gives two per-sample rows the same index: the projection
        // resolves both to the first slot and the second sample's coefficient and batch never
        // land anywhere, so a run scores one plant's reads under another's inbreeding.
        //
        // A **different order** is worse, because nothing anywhere would notice. The projection
        // reads the per-sample tables by *name* into this list and hands calling a vector indexed
        // by its position, while the writer writes them in the run's own sample order, which is
        // the read-group table's first-seen order. A file whose two orders disagree validates,
        // projects, and gives every sample its neighbour's inbreeding coefficient and batch.
        //
        // **This is the file's agreement with itself, not with the run** — whether these are the
        // *run's* samples is spec §6's second binding and step D2's.
        let mut first_seen: Vec<&str> = Vec::with_capacity(self.fitted_from.samples.len());
        for row in &self.fitted_from.read_groups {
            if !first_seen.contains(&row.sample.as_str()) {
                first_seen.push(row.sample.as_str());
            }
        }
        if first_seen != self.fitted_from.samples {
            let odd = first_seen
                .iter()
                .zip(&self.fitted_from.samples)
                .position(|(from_the_table, listed)| from_the_table != listed);
            return Err(refuse(
                "fitted_from.samples",
                match odd {
                    Some(at) => format!(
                        "names {:?} in position {at} where the read-group table's samples, in the \
                         order they first appear, have {:?}; every per-sample row is read by name \
                         into this list and handed to the run as a position, so a list in another \
                         order gives each sample its neighbour's numbers",
                        self.fitted_from.samples[at], first_seen[at]
                    ),
                    None => format!(
                        "holds {} names and the read-group table names {} distinct samples; the \
                         list is that table's samples, in the order they first appear, once each",
                        self.fitted_from.samples.len(),
                        first_seen.len()
                    ),
                },
            ));
        }
        Ok(())
    }

    /// **Every table keyed by the read-group axis covers it, and every per-sample table names
    /// samples the file declares.**
    ///
    /// The check above says the *declaration* runs `0..n`; this says the four tables written
    /// against that declaration agree with it. **The gap between the two is a hand edit**, which is
    /// the only way such a file arises — this writer builds all four from the dense vector itself.
    ///
    /// **And its symptom is spec §5's third row with no message.** The projection builds its
    /// vectors over `0..count` by *keyed lookup*, so a deleted calibration row does not shift the
    /// others: the missing id's slot becomes `ReadGroupCalibration::defaulted` — scale one, warrant
    /// `Defaulted`. The file's claim that this library was fitted disappears, the run reports it as
    /// a library nothing could be fitted for, and nothing anywhere says a row went missing.
    ///
    /// **A duplicate is refused for the same reason**: two rows for one read group means a keyed
    /// lookup discards one of them, and which one is an accident of order.
    fn every_axis_covers_what_it_is_keyed_by(&self) -> Result<(), ParametersFileError> {
        let read_groups = self.fitted_from.read_groups.len();
        covers_the_read_groups(
            "base_quality_calibration.by_read_group",
            self.base_quality_calibration
                .by_read_group
                .iter()
                .map(|row| row.read_group),
            read_groups,
        )?;
        if let Some(table) = &self.contamination {
            covers_the_read_groups(
                "contamination.by_read_group",
                table.by_read_group.iter().map(|row| row.read_group),
                read_groups,
            )?;
        }
        covers_the_read_groups(
            "sequencing_batches.by_read_group",
            self.sequencing_batches
                .by_read_group
                .iter()
                .map(|row| row.read_group),
            read_groups,
        )?;
        // **This one names read groups rather than covering them**, and that is the writer's
        // own rule: it writes a row only for a read group the run declared a slippage group for
        // (`from_run_parameters`'s `repeat_tracts_of`), because a read group the declaration does
        // not name has no slippage group — `StratumFits::at` answers `UnknownReadGroup` for it and
        // no slippage number is ever looked up under it. **A run with no repeat tracts at all
        // declares none**, which is Milestone E's defaults run and the single-sample case
        // `CLAUDE.md` puts first; requiring density here refused a file this caller had just
        // written.
        names_only_the_read_groups(
            "repeat_tracts.slippage_group_by_read_group",
            self.repeat_tracts
                .slippage_group_by_read_group
                .iter()
                .map(|row| row.read_group),
            read_groups,
        )?;
        // **The substitution rate is keyed by a read group too, and is likewise sparse** — a row
        // exists only where a rate was fitted for that `(read group × stratum × ploidy)`. What it
        // may not do is name a read group the identity block does not list, which is the mirror
        // of the gap refused above: the projection keys the map by that id, so a stray one is a
        // rate no locus can ever read.
        names_only_the_read_groups(
            "repeat_tracts.substitution_rate_by_stratum",
            self.repeat_tracts
                .substitution_rate_by_stratum
                .iter()
                .map(|row| row.read_group),
            read_groups,
        )?;
        self.no_repeat_tract_table_names_one_thing_twice()?;
        names_the_samples(
            "inbreeding.by_sample",
            self.inbreeding
                .by_sample
                .iter()
                .map(|row| row.sample.as_str()),
            &self.fitted_from.samples,
        )?;
        names_the_samples(
            "sequencing_batches.by_sample",
            self.sequencing_batches
                .by_sample
                .iter()
                .map(|row| row.sample.as_str()),
            &self.fitted_from.samples,
        )
    }

    /// **No repeat-tract table names one thing twice.**
    ///
    /// **A duplicate here is silent in a way a gap is not.** Each of these four tables becomes a
    /// map keyed by the row's own fields, so a second row for one key overwrites the first and
    /// nothing says which of the two the run scored under — the order the rows happen to sit in
    /// decides it. `StratumFits::over`, which builds the same map from the *fit's* outcomes,
    /// carries a release-level assert against exactly this, on the stated grounds that "the two
    /// levels can differ by a factor of five".
    ///
    /// **And this is the input path where it is plausible.** The fit cannot produce a duplicate;
    /// a person copying a `slippage_by_stratum_and_group` row to edit it and forgetting to change
    /// the stratum can. The refusal names the key, which the projection could not.
    fn no_repeat_tract_table_names_one_thing_twice(&self) -> Result<(), ParametersFileError> {
        let tracts = &self.repeat_tracts;
        names_each_key_once(
            "repeat_tracts.slippage_group_by_read_group",
            tracts
                .slippage_group_by_read_group
                .iter()
                .map(|row| format!("read_group = {}", row.read_group)),
        )?;
        names_each_key_once(
            "repeat_tracts.slippage_by_stratum_and_group",
            tracts.slippage_by_stratum_and_group.iter().map(|row| {
                format!(
                    "period = {}, reference_repeats = {}, slippage_group = {}",
                    row.period, row.reference_repeats, row.slippage_group
                )
            }),
        )?;
        names_each_key_once(
            "repeat_tracts.length_spectrum_by_stratum",
            tracts.length_spectrum_by_stratum.iter().map(|row| {
                format!(
                    "period = {}, reference_repeats = {}",
                    row.period, row.reference_repeats
                )
            }),
        )?;
        names_each_key_once(
            "repeat_tracts.length_spectrum_by_period",
            tracts
                .length_spectrum_by_period
                .iter()
                .map(|row| format!("period = {}", row.period)),
        )?;
        names_each_key_once(
            "repeat_tracts.substitution_rate_by_stratum",
            tracts.substitution_rate_by_stratum.iter().map(|row| {
                format!(
                    "read_group = {}, period = {}, reference_repeats = {}, ploidy = {}",
                    row.read_group, row.period, row.reference_repeats, row.ploidy
                )
            }),
        )
    }

    /// **A sample's libraries all ran in one batch, and it is the batch its own row names.**
    ///
    /// The batch is the population a contaminating read is drawn from (spec §3.4), and the file
    /// writes it twice — once a read group and once a sample — because the mixture reads both
    /// axes. **The two can disagree in a file and cannot in memory**:
    /// `SequencingBatches::declared` refuses a declaration in which one sample's libraries ran in
    /// two batches, naming the sample, so a file saying that is one no run could have produced
    /// and one the projection could not carry.
    ///
    /// **Its symptom is not a crash but a wrong neighbour.** A sample whose row says batch 0
    /// while its second library's row says batch 1 has its contaminant genotype drawn against one
    /// batch's frequencies and that library's reads scored against another's — two different
    /// populations, both plausible, in one sample.
    fn the_batching_puts_each_sample_in_one_batch(&self) -> Result<(), ParametersFileError> {
        let batch_of_read_group: BTreeMap<u32, u32> = self
            .sequencing_batches
            .by_read_group
            .iter()
            .map(|row| (row.read_group, row.batch))
            .collect();
        let batch_of_sample: BTreeMap<&str, u32> = self
            .sequencing_batches
            .by_sample
            .iter()
            .map(|row| (row.sample.as_str(), row.batch))
            .collect();
        for row in &self.fitted_from.read_groups {
            // **Both lookups are total** by the time this runs, and each rests on a check above:
            // the read-group one on the axis check, which refuses a batch table that does not
            // cover `0..n` once each, and the sample one on the identity check, which refuses a
            // sample list that is not the read-group table's own samples in first-seen order. The
            // second was not true until 2026-08-30, and a read-group row naming a sample the list
            // did not hold — a typo — slipped past this whole comparison. So it refuses rather
            // than skipping, and says which of the two is missing.
            let (Some(&of_read_group), Some(&of_sample)) = (
                batch_of_read_group.get(&row.read_group),
                batch_of_sample.get(row.sample.as_str()),
            ) else {
                return Err(refuse(
                    "sequencing_batches",
                    format!(
                        "does not batch read group {} and its sample {:?}: one of the two has no \
                         row, so nothing says which population that library's contaminating reads \
                         are drawn from",
                        row.read_group, row.sample
                    ),
                ));
            };
            if of_read_group != of_sample {
                return Err(refuse(
                    format!(
                        "sequencing_batches.by_read_group[read_group = {}]",
                        row.read_group
                    ),
                    format!(
                        "puts read group {} in batch {of_read_group} and \
                         `sequencing_batches.by_sample` puts its sample {:?} in batch \
                         {of_sample}; a sample's libraries all ran in one batch, because the \
                         batch is the population a contaminating read is drawn from and a sample \
                         has one",
                        row.read_group, row.sample
                    ),
                ));
            }
        }
        Ok(())
    }

    fn every_calibration_is_a_multiplier(&self) -> Result<(), ParametersFileError> {
        for row in &self.base_quality_calibration.by_read_group {
            let at = format!(
                "base_quality_calibration.by_read_group[read_group = {}].error_probability_multiplier",
                row.read_group
            );
            a_warranted_value(
                &at,
                &row.error_probability_multiplier,
                EvidenceCount::Reads(0),
            )?;
            // **Strictly above zero, and no upper bound.** A zero multiplies every read's error
            // probability to nothing, which charges the whole library the error floor — maximal
            // confidence about every base, from a number that says the fit found no errors at all;
            // `ReadGroupCalibration::from_fitted_rate` refuses it on the way in for that reason.
            // Above one is legitimate and common: it says the instrument was optimistic.
            if row.error_probability_multiplier.value <= 0.0 {
                return Err(refuse(
                    at,
                    format!(
                        "is {}, and a multiplier on an error probability is above zero — a zero \
                         charges every read of the library the error floor",
                        row.error_probability_multiplier.value
                    ),
                ));
            }
        }
        Ok(())
    }

    fn contamination_is_absent_or_measured(&self) -> Result<(), ParametersFileError> {
        let Some(table) = &self.contamination else {
            // **Absence is the uncontaminated run**, and it is a real model state rather than a
            // gap (spec §3.4). Nothing to check.
            return Ok(());
        };
        if table.by_read_group.is_empty() {
            return Err(refuse(
                "contamination.by_read_group",
                "is empty; an uncontaminated run has no [contamination] section at all, and an \
                 empty table says the section was written and then emptied",
            ));
        }
        if table
            .by_read_group
            .iter()
            .all(|row| row.measurement.is_none())
        {
            return Err(refuse(
                "contamination",
                "has a row for every read group and a measurement for none of them, which is the \
                 uncontaminated run written longhand; leave the section out instead, because a \
                 table of unmeasured rows takes the read likelihood's mixture path with every \
                 fraction zero where absence takes its plain one",
            ));
        }
        for row in &table.by_read_group {
            let Some(measurement) = &row.measurement else {
                continue;
            };
            let at = format!(
                "contamination.by_read_group[read_group = {}].measurement",
                row.read_group
            );
            // **Half-open at one, which is where the consumer's own bound is.**
            // `FrozenContamination::new` asserts `(0.0..1.0)` and says why: "a whole library of
            // another individual's DNA is not a sample of this one". A share of exactly one
            // accepted here becomes a panic several frames later, naming a read group rather than
            // a file — which is the failure this whole module exists to move earlier.
            let share_at = format!("{at}.share_of_reads_from_another_sample");
            finite(
                share_at.clone(),
                measurement.share_of_reads_from_another_sample,
            )?;
            if !(0.0..1.0).contains(&measurement.share_of_reads_from_another_sample) {
                return Err(refuse(
                    share_at,
                    format!(
                        "is {}, and a share of a lane's reads that came from somebody else is at \
                         or above zero and below one — a whole library of another individual is \
                         not a sample of this one",
                        measurement.share_of_reads_from_another_sample
                    ),
                ));
            }
            no_count_is_a_saturation(&at, measurement)?;
            // **Either count zero is the in-memory *not measured* shape**, and it is *either*
            // rather than *both* because that is what the predicate says:
            // `ContaminationView::was_measured` is `markers_with_reads > 0 && reads_on_markers >
            // 0` (`likelihood/mod.rs`), so its negation is a disjunction. A row with markers zero
            // and 90,233 reads is the worse of the two — it says *measured, 3.1%* in the file and
            // reads back in memory as never measured at all.
            //
            // **What such a row projects to is `UNMEASURED_READ_GROUP`**, the same value the
            // absent key gives, so nothing downstream computes differently; what goes wrong is
            // upstream of that. The file asserts a measurement, carries a
            // `fitted_from_reads_of` that is true of nothing, and the run reports a library as
            // uncorrected while the file says it was measured.
            if measurement.markers_with_reads == 0 || measurement.reads_on_markers == 0 {
                return Err(refuse(
                    at,
                    format!(
                        "says it was measured and carries the evidence of not being measured — \
                         markers_with_reads {} and reads_on_markers {}, where a measurement needs \
                         both above zero. Delete the whole `measurement` key instead: that is how \
                         this file says a lane could not be measured, and it is what this row \
                         already means to whatever reads it",
                        measurement.markers_with_reads, measurement.reads_on_markers
                    ),
                ));
            }
        }
        Ok(())
    }

    fn every_inbreeding_coefficient_is_a_fraction(&self) -> Result<(), ParametersFileError> {
        if self.inbreeding.by_sample.is_empty() {
            return Err(refuse(
                "inbreeding.by_sample",
                "is empty, and every sample of a run carries a coefficient",
            ));
        }
        for row in &self.inbreeding.by_sample {
            let at = format!(
                "inbreeding.by_sample[sample = {:?}].inbreeding_coefficient",
                row.sample
            );
            a_warranted_value(
                &at,
                &row.inbreeding_coefficient,
                EvidenceCount::CoveredPositions(0),
            )?;
            // **Half-open at one**, which is the range the type upstream states: a coefficient of
            // one is a sample with no heterozygous site anywhere, and the genotype prior it
            // produces has a zero where a likelihood is multiplied.
            let value = row.inbreeding_coefficient.value;
            if !(0.0..1.0).contains(&value) {
                return Err(refuse(
                    at,
                    format!(
                        "is {value}; an inbreeding coefficient is at or above zero and strictly \
                         below one — at one every heterozygote becomes impossible, and a fitted \
                         estimate should sit well below it"
                    ),
                ));
            }
        }
        Ok(())
    }

    fn the_prior_seed_is_a_pair_of_concentrations(&self) -> Result<(), ParametersFileError> {
        let seed = &self.ordinary_site_prior;
        finite(
            "ordinary_site_prior.reference_concentration",
            seed.reference_concentration,
        )?;
        finite(
            "ordinary_site_prior.alternative_concentration_total",
            seed.alternative_concentration_total,
        )?;
        if seed.reference_concentration <= 0.0 {
            return Err(refuse(
                "ordinary_site_prior.reference_concentration",
                format!(
                    "is {}, and the reference allele carries some belief at every position",
                    seed.reference_concentration
                ),
            ));
        }
        // **Zero is a real answer here and is not floored** — a cohort with no variation at all
        // (spec §3.6) — so this refuses only a negative total.
        if seed.alternative_concentration_total < 0.0 {
            return Err(refuse(
                "ordinary_site_prior.alternative_concentration_total",
                format!(
                    "is {}, and a concentration is not negative; zero is a real answer and means \
                     a cohort with no variation",
                    seed.alternative_concentration_total
                ),
            ));
        }
        Ok(())
    }

    fn the_repeat_tract_numbers_are_what_they_claim(&self) -> Result<(), ParametersFileError> {
        let tracts = &self.repeat_tracts;
        let at = "repeat_tracts.fallback_length_spectrum_concentration";
        // **A median over strata has no count of its own**, so any unit here is wrong; the
        // `defaulted` rule above already refuses one on the flat rung, and `reads` is named only
        // because the check needs a word and this key should carry no `observations` at all.
        a_warranted_value(
            at,
            &tracts.fallback_length_spectrum_concentration,
            EvidenceCount::Reads(0),
        )?;
        // **Its warrant is not free: the file's own strata decide it.** The bottom rung states
        // the median of the concentrations this run's strata fitted wherever any was fitted, and
        // the compiled-in flat constant only where none was — so `defaulted` beside a non-empty
        // `length_spectrum_by_stratum`, or `fitted_here` beside an empty one, is a claim the
        // file's own rows refute. **It is refused here because nothing downstream can**: the
        // projection carries only the number, and the writer re-derives the warrant from the same
        // strata, so a contradiction is rewritten on the way out rather than reported.
        let fitted_any = !tracts.length_spectrum_by_stratum.is_empty();
        let warrant = tracts.fallback_length_spectrum_concentration.warrant;
        if fitted_any != (warrant == Warrant::FittedHere) {
            return Err(refuse(
                at,
                if fitted_any {
                    format!(
                        "is `{}`, and {} of this file's strata were fitted on their own tracts; \
                         the bottom rung states the median of those, so its warrant is \
                         `fitted_here` — it is `defaulted` only where no stratum was fitted at all",
                        the_word_for(warrant),
                        tracts.length_spectrum_by_stratum.len()
                    )
                } else {
                    format!(
                        "is `{}`, and no stratum in this file was fitted on its own tracts, so \
                         there is no median to take; a run with nothing fitted states the \
                         compiled-in flat concentration and marks it `defaulted`",
                        the_word_for(warrant)
                    )
                },
            ));
        }
        if tracts.fallback_length_spectrum_concentration.value <= 0.0 {
            return Err(refuse(
                at,
                format!(
                    "is {}, and a concentration is above zero",
                    tracts.fallback_length_spectrum_concentration.value
                ),
            ));
        }

        for row in &tracts.slippage_by_stratum_and_group {
            let at = format!(
                "repeat_tracts.slippage_by_stratum_and_group[period = {}, reference_repeats = {}, \
                 slippage_group = {}]",
                row.period, row.reference_repeats, row.slippage_group
            );
            a_share(
                &format!("{at}.share_of_reads_that_slip"),
                row.share_of_reads_that_slip,
            )?;
            a_share(&format!("{at}.shorter_share"), row.shorter_share)?;
            // **`fall_off` is checked for being a number and not for an upper bound**, because
            // neither this file's shape nor the fit that produces it documents one. It is how fast
            // two-repeat slips fall off against one-repeat slips, so a value above one would mean
            // the larger slip is the likelier — implausible, but nothing here has established that
            // it is impossible, and refusing a fit nobody has bounded is worse than passing one.
            finite(format!("{at}.fall_off"), row.fall_off)?;
            if row.fall_off < 0.0 {
                return Err(refuse(
                    format!("{at}.fall_off"),
                    format!("is {}, and a fall-off is not negative", row.fall_off),
                ));
            }
            a_level_smoothing(&at, &row.share_of_reads_that_slip_origin.smoothing)?;
            if let Some(reads) = row.share_of_reads_that_slip_origin.expected_slipped_reads {
                finite(
                    format!("{at}.share_of_reads_that_slip_origin.expected_slipped_reads"),
                    reads,
                )?;
            }
            if let Some(shares) = &row.shorter_share_and_fall_off_origin {
                if let Some(reads) = shares.expected_slipped_reads {
                    finite(
                        format!("{at}.shorter_share_and_fall_off_origin.expected_slipped_reads"),
                        reads,
                    )?;
                }
                a_share_smoothing(
                    &format!("{at}.shorter_share_and_fall_off_origin.shorter_share_smoothing"),
                    &shares.shorter_share_smoothing,
                )?;
                a_share_smoothing(
                    &format!("{at}.shorter_share_and_fall_off_origin.fall_off_smoothing"),
                    &shares.fall_off_smoothing,
                )?;
            }
        }

        for row in &tracts.length_spectrum_by_stratum {
            let at = format!(
                "repeat_tracts.length_spectrum_by_stratum[period = {}, reference_repeats = {}]",
                row.period, row.reference_repeats
            );
            a_length_spectrum(&at, row.concentration, &row.shares_by_repeat_offset)?;
        }
        for row in &tracts.length_spectrum_by_period {
            let at = format!(
                "repeat_tracts.length_spectrum_by_period[period = {}]",
                row.period
            );
            a_length_spectrum(&at, row.concentration, &row.shares_by_repeat_offset)?;
        }

        for row in &tracts.substitution_rate_by_stratum {
            let at = format!(
                "repeat_tracts.substitution_rate_by_stratum[read_group = {}, period = {}, \
                 reference_repeats = {}, ploidy = {}].rate",
                row.read_group, row.period, row.reference_repeats, row.ploidy
            );
            a_warranted_value(&at, &row.rate, EvidenceCount::BasesCompared(0))?;
            finite(at.clone(), row.rate.value)?;
            if !(0.0..=1.0).contains(&row.rate.value) {
                return Err(refuse(
                    at,
                    format!(
                        "is {}, and a substitution rate is a probability — the chance one base \
                         inside a tract reads wrong — so it lies between zero and one",
                        row.rate.value
                    ),
                ));
            }
        }
        Ok(())
    }

    fn every_stated_constant_is_in_range(&self) -> Result<(), ParametersFileError> {
        let at = "stated_constants.repeat_tract_outlier_weight";
        let weight = &self.stated_constants.repeat_tract_outlier_weight;
        a_warranted_value(at, weight, EvidenceCount::Reads(0))?;
        // **Open at both ends, where every other share here is closed.** The scoring row asserts
        // `0 < weight < 1` (`likelihood::ssr`'s `genotype_log_likelihood_row`), so a zero or a
        // one accepted here becomes a panic several frames later naming a locus rather than the
        // file the number came from — the same shape as the contamination share of exactly one
        // this module already refuses. It is the one number in this file a person is *invited*
        // to edit (spec §3.8), so the ends are the values most worth catching here.
        if !(weight.value > 0.0 && weight.value < 1.0) {
            return Err(refuse(
                at,
                format!(
                    "is {}, and the share of a repeat tract's reads the model cannot explain is \
                     strictly inside zero and one — a zero says no read at a tract can have come \
                     from anywhere but this sample's own copies of it, and a one says none of \
                     them did",
                    weight.value
                ),
            ));
        }
        // **Two of the four warrants, and a `defaulted` one has to be the constant.** Nothing
        // fits this number, so `fitted_here` and `borrowed` are claims about it that no run
        // could make, and the in-memory shape it projects onto
        // ([`RepeatTractOutlierWeight`](crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight))
        // has nowhere to put them. **The state worth catching is the third one**: a person who
        // takes spec §1.2 goal 3 at its word — copy the file your run wrote and change one line
        // — changes the number and leaves `warrant = "defaulted"` above it. That says *this run
        // inherited 0.01* beside a number that is not 0.01, and the file's own header names the
        // fix ("change its warrant to supplied and delete its observations"), so the refusal
        // says it too.
        match weight.warrant {
            Warrant::Supplied => Ok(()),
            Warrant::Defaulted if weight.value == DEFAULT_OUTLIER_WEIGHT => Ok(()),
            Warrant::Defaulted => Err(refuse(
                at,
                format!(
                    "is {}, and its warrant is `defaulted`, which says this run inherited the \
                     compiled-in {DEFAULT_OUTLIER_WEIGHT}; a number you changed is one the run \
                     was handed, so change the warrant beside it to `supplied`",
                    weight.value
                ),
            )),
            other => Err(refuse(
                at,
                format!(
                    "has the warrant `{}`, and nothing fits this number: it is either the \
                     compiled-in {DEFAULT_OUTLIER_WEIGHT}, which is `defaulted`, or one you \
                     wrote here, which is `supplied`",
                    the_word_for(other)
                ),
            )),
        }
    }
}

/// The file's own spelling of a warrant, for a refusal that has to quote one back.
///
/// **The `serde` spelling rather than `{:?}`**, which would print the Rust variant name — a
/// word that appears nowhere in the file the reader is being sent to edit.
fn the_word_for(warrant: Warrant) -> &'static str {
    match warrant {
        Warrant::FittedHere => "fitted_here",
        Warrant::Borrowed => "borrowed",
        Warrant::Supplied => "supplied",
        Warrant::Defaulted => "defaulted",
    }
}

/// One table keyed by the read-group axis: every id from `0..axis_length`, once each.
fn covers_the_read_groups(
    at: &str,
    ids: impl Iterator<Item = u32>,
    axis_length: usize,
) -> Result<(), ParametersFileError> {
    let mut seen: Vec<u32> = ids.collect();
    seen.sort_unstable();
    if seen.len() != axis_length {
        return Err(refuse(
            at.to_string(),
            format!(
                "holds {} rows where fitted_from.read_groups declares {axis_length}; every table \
                 keyed by a read group carries one row for each, because the run reads them by id \
                 and a missing row is silently a defaulted one",
                seen.len()
            ),
        ));
    }
    for (expected, found) in seen.iter().enumerate() {
        if *found as usize != expected {
            return Err(refuse(
                at.to_string(),
                format!(
                    "names read group {found} where {expected} should be, so it either repeats one \
                     or skips one; the run reads these by id",
                ),
            ));
        }
    }
    Ok(())
}

/// One table keyed by sample name: every declared sample, once each, and no others.
fn names_the_samples<'a>(
    at: &str,
    names: impl Iterator<Item = &'a str>,
    declared: &[String],
) -> Result<(), ParametersFileError> {
    let mut seen: Vec<&str> = names.collect();
    seen.sort_unstable();
    let mut want: Vec<&str> = declared.iter().map(String::as_str).collect();
    want.sort_unstable();
    if seen == want {
        return Ok(());
    }
    // **Name the one that differs rather than printing both lists.** At the top of the committed
    // range a cohort has 3,000 samples, and two 3,000-name lists in an error message is not a
    // message. The first name in one list and not the other is what a reader needs.
    let odd_one_out = seen
        .iter()
        .find(|name| !want.contains(name))
        .or_else(|| want.iter().find(|name| !seen.contains(name)));
    Err(refuse(
        at.to_string(),
        match odd_one_out {
            Some(name) => format!(
                "does not name the same samples as fitted_from.samples: {name:?} is in one and \
                 not the other",
            ),
            // Same set, different multiplicity — a duplicated row.
            None => format!(
                "holds {} rows for {} samples, so it names one of them twice",
                seen.len(),
                want.len()
            ),
        },
    ))
}

/// A short list of read-group ids for a message — `0, 2, 3`.
///
/// **Capped, because the axis is one of the file's cohort-sized ones.** At the top of the
/// committed input range a run has 3,000 read groups (`CLAUDE.md`), and a message that printed all
/// of them would bury the id that is missing under the 2,999 that are not.
fn a_list_of(ids: &[u32]) -> String {
    const AT_MOST: usize = 12;
    let shown: Vec<String> = ids.iter().take(AT_MOST).map(u32::to_string).collect();
    if ids.len() > AT_MOST {
        format!("{} and {} more", shown.join(", "), ids.len() - AT_MOST)
    } else {
        shown.join(", ")
    }
}

/// **One table keyed by the read-group axis whose rows are a subset of it**, rather than a cover.
///
/// **Two tables of the file are sparse by construction** — the slippage-group declaration and the
/// repeat-tract substitution rate — because a row exists only where the run had something to say.
/// What neither may do is name a read group the identity block does not list: the projection keys
/// its map by that id, so a stray one is a number no locus can ever reach, and it is the mirror of
/// the gap [`covers_the_read_groups`] refuses.
fn names_only_the_read_groups(
    at: &str,
    ids: impl Iterator<Item = u32>,
    read_groups: usize,
) -> Result<(), ParametersFileError> {
    for id in ids {
        if id as usize >= read_groups {
            return Err(refuse(
                at.to_string(),
                format!(
                    "names read group {id} and this run has {read_groups}, numbered from zero; a \
                     row keyed to a library the identity block does not list is one nothing can \
                     ever read"
                ),
            ));
        }
    }
    Ok(())
}

/// **One table whose rows each name a different thing.**
///
/// The keys arrive already spelled the way the file spells them, so the refusal can quote the
/// repeated one back — which is what a reader needs and what the projection, having already lost
/// one of the two rows into a map, could not give.
fn names_each_key_once(
    at: &str,
    keys: impl Iterator<Item = String>,
) -> Result<(), ParametersFileError> {
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for key in keys {
        if !seen.insert(key.clone()) {
            return Err(refuse(
                at.to_string(),
                format!(
                    "holds two rows for `{key}`; each row of this table becomes one entry of a \
                     map keyed by exactly those fields, so the second silently replaces the first \
                     and which one a run scored under is the order they happen to sit in"
                ),
            ));
        }
    }
    Ok(())
}

/// The file's own word for a count's unit, for a refusal that has to name one.
fn the_unit_of(count: EvidenceCount) -> &'static str {
    match count {
        EvidenceCount::Reads(_) => "reads",
        EvidenceCount::CoveredPositions(_) => "covered_positions",
        EvidenceCount::BasesCompared(_) => "bases_compared",
    }
}

/// The refusal, with the key's path in the file's own spelling.
///
/// **Shared with the projection** (`to_run_parameters`), which refuses in the same words where a
/// value passes every check here and a newtype's own constructor still turns it down — a motif
/// period past the longest this build indexes, say. Two spellings of one refusal is what having
/// one function stops.
pub(super) fn refuse(field: impl Into<String>, problem: impl Into<String>) -> ParametersFileError {
    ParametersFileError::Meaningless {
        field: field.into(),
        problem: problem.into(),
    }
}

/// A number that is a number — refusing `NaN` and both infinities.
///
/// **Its own check rather than a corollary of the range tests**, because `NaN` fails every
/// comparison, so `!(0.0..=1.0).contains(&f64::NAN)` is true and would report a `NaN` as *outside
/// [0, 1]* — which sends a reader looking for a value they will not find in the file.
fn finite(at: impl Into<String>, value: f64) -> Result<(), ParametersFileError> {
    if value.is_finite() {
        return Ok(());
    }
    Err(refuse(
        at,
        format!("is {value}, which is not a number a score can be computed from"),
    ))
}

/// A number that is a share of something: finite, and within `[0, 1]`.
fn a_share(at: &str, value: f64) -> Result<(), ParametersFileError> {
    finite(at, value)?;
    if (0.0..=1.0).contains(&value) {
        return Ok(());
    }
    Err(refuse(
        at.to_string(),
        format!("is {value}, and a share is in [0, 1]"),
    ))
}

/// A value and its warrant: the number is a number, its evidence count is a count, and a
/// defaulted number has no count at all.
///
/// **That last one is the writer's own rule, enforced on the way back in.** A stated constant has
/// nothing behind it, and a count beside one would say that the constant rests on those reads —
/// so the projection out writes no `observations` for a `defaulted` value, of any quantity. A
/// file that carries one is either hand-edited or from another writer, and reading it would put
/// evidence behind a number that has none. **The fix is the one the file's own header teaches**:
/// a number you changed is `supplied`, and a supplied number keeps its count precisely because it
/// was fitted on *some* cohort.
fn a_warranted_value(
    at: &str,
    value: &WarrantedValue,
    counted_in: EvidenceCount,
) -> Result<(), ParametersFileError> {
    finite(at, value.value)?;
    let Some(observations) = value.observations else {
        return Ok(());
    };
    if value.warrant == Warrant::Defaulted {
        return Err(refuse(
            format!("{at}.observations"),
            "is written beside a `defaulted` warrant, and a stated constant has nothing behind \
             it; delete the `observations` table, or — if you changed the number — set the \
             warrant to `supplied`, which keeps its count",
        ));
    }
    // **The unit has to be the one this quantity is fitted over**, and this is the only thing
    // that says so. The projection *in* drops the unit — `Estimate<T>`'s count is a bare `u64`
    // whose unit follows the quantity — and the projection out mints it back from the call site.
    // So a calibration row that says `covered_positions = 812344` reads back as 812,344 reads and
    // is written out under a key the user did not type: a change of meaning, silently, in the one
    // direction spec §1.2 goal 1 forbids. An inbreeding coefficient is fitted over covered
    // reference positions and a repeat-tract substitution rate over bases compared; neither is a
    // read.
    if std::mem::discriminant(&observations) != std::mem::discriminant(&counted_in) {
        return Err(refuse(
            format!("{at}.observations"),
            format!(
                "counts {}, and this number is fitted over {}; the unit is not decoration — the \
                 three differ by orders of magnitude on one cohort, and a run reading this would \
                 report the count under the other name",
                the_unit_of(observations),
                the_unit_of(counted_in)
            ),
        ));
    }
    let (unit, count) = match observations {
        EvidenceCount::Reads(count) => ("reads", count),
        EvidenceCount::CoveredPositions(count) => ("covered_positions", count),
        EvidenceCount::BasesCompared(count) => ("bases_compared", count),
    };
    if count == A_SATURATED_COUNT {
        return Err(refuse(
            format!("{at}.observations.{unit}"),
            format!(
                "is {A_SATURATED_COUNT}, which is the largest integer TOML holds and is what this \
                 writer saturates a count to rather than emit one no reader agrees on; it is a \
                 lost number and not a measurement"
            ),
        ));
    }
    Ok(())
}

/// Neither of a measurement's two counts is a saturation marker.
fn no_count_is_a_saturation(
    at: &str,
    measurement: &ContaminationMeasurement,
) -> Result<(), ParametersFileError> {
    for (key, count) in [
        ("markers_with_reads", measurement.markers_with_reads),
        ("reads_on_markers", measurement.reads_on_markers),
    ] {
        if count == A_SATURATED_COUNT {
            return Err(refuse(
                format!("{at}.{key}"),
                format!(
                    "is {A_SATURATED_COUNT}, which is this writer's saturation marker and not a \
                     count"
                ),
            ));
        }
    }
    Ok(())
}

/// A length spectrum: a positive concentration, and shares that are a distribution over an odd
/// number of offsets with a middle.
fn a_length_spectrum(
    at: &str,
    concentration: f64,
    shares: &[f64],
) -> Result<(), ParametersFileError> {
    finite(format!("{at}.concentration"), concentration)?;
    if concentration <= 0.0 {
        return Err(refuse(
            format!("{at}.concentration"),
            format!("is {concentration}, and a concentration is above zero"),
        ));
    }

    // **Odd, and at least three.** The array runs from `-span` to `+span` in whole repeat units
    // from the reference length, so the middle entry *is* the reference length; an even count has
    // no middle, and a count below three cannot express a slip in both directions.
    if shares.len() < 3 || shares.len().is_multiple_of(2) {
        return Err(refuse(
            format!("{at}.shares_by_repeat_offset"),
            format!(
                "holds {} shares; it runs from -span to +span in whole repeats, so the count is \
                 odd and at least three — the middle entry is the reference length itself",
                shares.len()
            ),
        ));
    }
    let mut total = 0.0;
    for (offset, share) in shares.iter().enumerate() {
        a_share(&format!("{at}.shares_by_repeat_offset[{offset}]"), *share)?;
        total += *share;
    }
    if (total - 1.0).abs() > SHARES_MAY_MISS_ONE_BY {
        return Err(refuse(
            format!("{at}.shares_by_repeat_offset"),
            format!(
                "sums to {total}, and a length spectrum is a distribution over the lengths a \
                 tract can take; if you edited one share, the others have to give up what it took"
            ),
        ));
    }
    Ok(())
}

/// A level's smoothing: the weight a curve carried is a share, and the curve's own numbers are
/// numbers.
///
/// **The curve is provenance rather than a term in a score** — it is recorded so an interpolation
/// can be told from a measurement — so a `NaN` in it changes no genotype. It is still refused,
/// because a run that reports where a number came from should not report `nan` as the slope of the
/// line it came off, and because a reader who hand-edits a curve gets told rather than ignored.
fn a_level_smoothing(at: &str, smoothing: &LevelSmoothing) -> Result<(), ParametersFileError> {
    let at = format!("{at}.share_of_reads_that_slip_origin.smoothing");
    match smoothing {
        LevelSmoothing::ThisStratum => Ok(()),
        LevelSmoothing::ThisPeriodsCurve { curve, .. } => {
            a_slippage_curve(&format!("{at}.this_periods_curve.curve"), curve)
        }
        LevelSmoothing::Blend {
            curve_weight,
            curve,
            ..
        } => {
            a_share(&format!("{at}.blend.curve_weight"), *curve_weight)?;
            a_slippage_curve(&format!("{at}.blend.curve"), curve)
        }
    }
}

/// The same, for a share's smoothing.
fn a_share_smoothing(at: &str, smoothing: &ShareSmoothing) -> Result<(), ParametersFileError> {
    match smoothing {
        ShareSmoothing::ThisStratum => Ok(()),
        ShareSmoothing::ThisPeriodsCurve { curve, .. } => {
            a_share_curve(&format!("{at}.this_periods_curve.curve"), curve)
        }
        ShareSmoothing::Blend {
            curve_weight,
            curve,
            ..
        } => {
            a_share(&format!("{at}.blend.curve_weight"), *curve_weight)?;
            a_share_curve(&format!("{at}.blend.curve"), curve)
        }
    }
}

/// Every number a level curve records is a number.
fn a_slippage_curve(at: &str, curve: &SlippageCurve) -> Result<(), ParametersFileError> {
    for (key, value) in [
        ("rise_shape", curve.rise_shape),
        ("intercept", curve.intercept),
        ("slope", curve.slope),
        ("held_out_error", curve.held_out_error),
    ] {
        finite(format!("{at}.{key}"), value)?;
    }
    Ok(())
}

/// Every number a share curve records is a number.
fn a_share_curve(at: &str, curve: &ShareCurve) -> Result<(), ParametersFileError> {
    for (key, value) in [
        ("intercept", curve.intercept),
        ("slope", curve.slope),
        ("bend", curve.bend),
        ("centre_repeats", curve.centre_repeats),
        ("held_out_error", curve.held_out_error),
    ] {
        finite(format!("{at}.{key}"), value)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::tests::{
        THE_ROW_WHOSE_SHARES_BLEND, THE_ROW_WHOSE_SLIP_SHARE_BLENDS, a_file_using_every_shape,
    };
    use super::super::{
        ContaminationFittedFrom, ContaminationMeasurement, EvidenceCount, LevelSmoothing, SeedRung,
        ShareSmoothing, Warrant, WarrantedValue,
    };
    use super::*;

    /// The fixture with one edit, refused — returning the refusal so a test can read its field.
    #[track_caller]
    fn refused(edit: impl FnOnce(&mut ParametersFile)) -> (String, String) {
        let mut file = a_file_using_every_shape();
        edit(&mut file);
        match file.validate() {
            Ok(()) => panic!("this edit was accepted and should have been refused"),
            Err(ParametersFileError::Meaningless { field, problem }) => (field, problem),
            Err(other) => panic!("refused for the wrong reason: {other}"),
        }
    }

    /// The fixture with one edit, accepted.
    #[track_caller]
    fn accepted(edit: impl FnOnce(&mut ParametersFile)) {
        let mut file = a_file_using_every_shape();
        edit(&mut file);
        file.validate()
            .expect("this edit is legitimate and must be accepted");
    }

    /// **The file this project writes passes its own reader**, which is the check every refusal
    /// below is only meaningful against.
    #[test]
    fn the_file_a_run_writes_is_accepted() {
        a_file_using_every_shape()
            .validate()
            .expect("the writer's own fixture is a file a run can use");
    }

    /// **And so is the smallest file that still says something** — one sample, one read group,
    /// no repeat tracts, no contamination.
    ///
    /// The committed input range starts at a single sample (`CLAUDE.md`), and a `validate` written
    /// against the full fixture alone would happily refuse the bottom of that range by requiring
    /// a table that a one-sample run legitimately leaves empty.
    #[test]
    fn the_smallest_file_that_says_anything_is_accepted() {
        let mut small = a_file_using_every_shape();
        small.contamination = None;
        small.fitted_from.samples.truncate(1);
        small.fitted_from.read_groups.truncate(1);
        small.base_quality_calibration.by_read_group.truncate(1);
        small.sequencing_batches.by_read_group.truncate(1);
        small.sequencing_batches.by_sample.truncate(1);
        small.inbreeding.by_sample.truncate(1);
        small.repeat_tracts.slippage_group_by_read_group.clear();
        small.repeat_tracts.slippage_by_stratum_and_group.clear();
        small.repeat_tracts.length_spectrum_by_stratum.clear();
        small.repeat_tracts.length_spectrum_by_period.clear();
        small.repeat_tracts.substitution_rate_by_stratum.clear();
        // **A run that fitted no stratum takes the compiled-in flat concentration**, which is
        // what `defaulted` says; the fixture's `fitted_here` is a claim about a median over
        // strata this file has none of.
        small.repeat_tracts.fallback_length_spectrum_concentration = WarrantedValue {
            value: 1.0,
            warrant: Warrant::Defaulted,
            observations: None,
        };
        small
            .validate()
            .expect("one sample with no repeat tracts is the bottom of the committed range");
    }

    /// **The sample list is the read-group table's own samples, in first-seen order, once each.**
    ///
    /// Two failures, and neither has a symptom the reader would recognise. A **repeat** gives two
    /// per-sample rows one index, so the second sample's coefficient never lands and the
    /// projection panics several frames from the key. A **different order** does not panic at
    /// all: the projection reads the per-sample tables by name into this list and hands calling a
    /// vector indexed by its position, while the writer writes them in the run's own order, so a
    /// file whose two orders disagree gives every sample its neighbour's inbreeding coefficient
    /// and its neighbour's sequencing batch, silently.
    #[test]
    fn a_sample_list_that_is_not_the_read_group_tables_own_is_refused() {
        let (field, problem) = refused(|file| file.fitted_from.samples.swap(0, 1));
        assert_eq!(field, "fitted_from.samples");
        assert!(problem.contains("its neighbour's numbers"), "{problem}");

        let (field, problem) = refused(|file| {
            let repeated = file.fitted_from.samples[0].clone();
            file.fitted_from.samples.push(repeated);
        });
        assert_eq!(field, "fitted_from.samples");
        assert!(problem.contains("distinct samples"), "{problem}");

        // **A read group naming a sample the list does not hold** — a typo — used to slip past
        // the batching comparison, because the sample side of that lookup simply found nothing
        // and the row was skipped.
        let (field, _) = refused(|file| {
            file.fitted_from.read_groups[1].sample = "TS-l".into();
        });
        assert_eq!(field, "fitted_from.samples");
    }

    /// **A repeat-tract table that names one thing twice is refused, naming the key.**
    ///
    /// Each of these five tables becomes a map keyed by the row's own fields, so a second row for
    /// one key silently replaces the first and which of the two a run scored under is the order
    /// they sit in. `StratumFits::over` carries a release-level assert against exactly this on the
    /// fit's side; a file is the input path where a person copying a row to edit it can produce
    /// one.
    #[test]
    fn a_repeat_tract_table_that_names_one_thing_twice_is_refused() {
        let (field, problem) = refused(|file| {
            let row = file.repeat_tracts.slippage_by_stratum_and_group[0].clone();
            file.repeat_tracts.slippage_by_stratum_and_group.push(row);
        });
        assert_eq!(field, "repeat_tracts.slippage_by_stratum_and_group");
        assert!(
            problem.contains("period = 1, reference_repeats = 30, slippage_group = 0"),
            "the refusal quotes the repeated key back: {problem}"
        );

        for (table, at) in [
            (
                "repeat_tracts.length_spectrum_by_stratum",
                refused(|file| {
                    let row = file.repeat_tracts.length_spectrum_by_stratum[0].clone();
                    file.repeat_tracts.length_spectrum_by_stratum.push(row);
                }),
            ),
            (
                "repeat_tracts.length_spectrum_by_period",
                refused(|file| {
                    let row = file.repeat_tracts.length_spectrum_by_period[0].clone();
                    file.repeat_tracts.length_spectrum_by_period.push(row);
                }),
            ),
            (
                "repeat_tracts.substitution_rate_by_stratum",
                refused(|file| {
                    let row = file.repeat_tracts.substitution_rate_by_stratum[0].clone();
                    file.repeat_tracts.substitution_rate_by_stratum.push(row);
                }),
            ),
            (
                "repeat_tracts.slippage_group_by_read_group",
                refused(|file| {
                    let row = file.repeat_tracts.slippage_group_by_read_group[0];
                    file.repeat_tracts.slippage_group_by_read_group.push(row);
                }),
            ),
        ] {
            assert_eq!(at.0, table);
        }
    }

    /// **An evidence count spelled in another quantity's unit is refused.**
    ///
    /// The three units differ by orders of magnitude on one cohort, and the projection drops the
    /// unit on the way in — `Estimate<T>`'s count is a bare number whose unit follows the
    /// quantity — so a calibration row counting `covered_positions` reads back as that many
    /// *reads* and is written out under a key the user did not type. This is the only thing that
    /// stops it.
    #[test]
    fn an_evidence_count_in_another_quantitys_unit_is_refused() {
        let (field, problem) = refused(|file| {
            file.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .observations = Some(EvidenceCount::CoveredPositions(812_344));
        });
        assert!(field.ends_with(".observations"), "{field}");
        assert!(
            problem.contains("counts covered_positions") && problem.contains("over reads"),
            "{problem}"
        );

        let (field, _) = refused(|file| {
            file.inbreeding.by_sample[0]
                .inbreeding_coefficient
                .observations = Some(EvidenceCount::Reads(180_600_412));
        });
        assert!(field.starts_with("inbreeding.by_sample"), "{field}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0]
                .rate
                .observations = Some(EvidenceCount::Reads(40_122));
        });
        assert!(field.ends_with(".observations"), "{field}");
    }

    /// **A defaulted value carrying an evidence count is refused**, which is the writer's own
    /// rule read back in: a stated constant has nothing behind it, and a count beside one says
    /// the constant rests on those reads.
    #[test]
    fn a_defaulted_value_that_claims_evidence_is_refused() {
        let (field, problem) = refused(|file| {
            file.base_quality_calibration.by_read_group[1]
                .error_probability_multiplier
                .observations = Some(EvidenceCount::Reads(4));
        });
        assert!(field.ends_with(".observations"), "{field}");
        assert!(
            problem.contains("set the warrant to `supplied`"),
            "{problem}"
        );
    }

    /// **The fallback concentration's warrant is decided by the file's own strata**, so a claim
    /// its rows refute is refused rather than silently rewritten on the way out.
    #[test]
    fn a_fallback_warrant_the_files_own_strata_refute_is_refused() {
        let (field, problem) = refused(|file| {
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .warrant = Warrant::Defaulted;
        });
        assert_eq!(
            field,
            "repeat_tracts.fallback_length_spectrum_concentration"
        );
        assert!(problem.contains("median of those"), "{problem}");

        let (field, problem) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_stratum.clear();
        });
        assert_eq!(
            field,
            "repeat_tracts.fallback_length_spectrum_concentration"
        );
        assert!(problem.contains("no median to take"), "{problem}");
    }

    #[test]
    fn a_version_this_build_does_not_read_is_refused() {
        let (field, problem) = refused(|file| file.format_version = 0);
        assert_eq!(field, "format_version");
        assert!(problem.contains("written by a newer build"), "{problem}");
        let (field, _) = refused(|file| file.format_version = 2);
        assert_eq!(field, "format_version");
    }

    #[test]
    fn a_run_that_calls_no_copies_and_a_cohort_with_no_samples_are_refused() {
        assert_eq!(refused(|file| file.ploidy = 0).0, "ploidy");
        assert_eq!(
            refused(|file| file.fitted_from.samples.clear()).0,
            "fitted_from.samples"
        );
        assert_eq!(
            refused(|file| file.fitted_from.read_groups.clear()).0,
            "fitted_from.read_groups"
        );
        assert_eq!(
            refused(|file| file.inbreeding.by_sample.clear()).0,
            "inbreeding.by_sample"
        );
    }

    /// **A gap in the read-group ids is refused here rather than at a locus.**
    ///
    /// `RunParameters` indexes the calibration and contamination axes by read-group id, so a gap
    /// drops the highest read group entirely; its symptom is a panic at whichever locus first
    /// carries one of that library's reads, which names a locus and arrives after the pre-pass is
    /// long finished.
    #[test]
    fn a_gap_in_the_read_group_ids_is_refused() {
        let (field, problem) = refused(|file| file.fitted_from.read_groups[1].read_group = 3);
        assert_eq!(field, "fitted_from.read_groups");
        assert!(problem.contains("has no read group 1"), "{problem}");
    }

    /// **A table that does not cover the read-group axis is refused, and that is a hand edit's
    /// failure rather than a writer's.**
    ///
    /// This writer builds all four of these tables from the dense vector itself, so it cannot
    /// produce a gap. A person deleting one row can, and the symptom is silent: the projection
    /// reads these by id, so a missing calibration row does not shift the others — it becomes a
    /// defaulted scale of one, and the file's claim that the library was fitted is gone with no
    /// message. That is spec §5's third row arriving through the back door.
    #[test]
    fn a_table_that_does_not_cover_the_read_groups_is_refused() {
        let (field, problem) = refused(|file| {
            file.base_quality_calibration.by_read_group.remove(1);
        });
        assert_eq!(field, "base_quality_calibration.by_read_group");
        assert!(problem.contains("silently a defaulted one"), "{problem}");

        assert_eq!(
            refused(|file| {
                file.contamination
                    .as_mut()
                    .expect("a table")
                    .by_read_group
                    .remove(2);
            })
            .0,
            "contamination.by_read_group"
        );
        assert_eq!(
            refused(|file| {
                file.sequencing_batches.by_read_group.remove(0);
            })
            .0,
            "sequencing_batches.by_read_group"
        );
        // **The slippage-group declaration is the one read-group table that is *not* dense**,
        // and this is where that is pinned. The writer names a read group there only where the
        // run declared a slippage group for it — a run with no repeat tracts declares none at all
        // — so requiring a cover here refused a file this caller had just written. What is still
        // refused is a row naming a library the identity block does not list.
        accepted(|file| {
            file.repeat_tracts.slippage_group_by_read_group.remove(0);
        });
        let (field, problem) = refused(|file| {
            file.repeat_tracts.slippage_group_by_read_group[0].read_group = 9;
        });
        assert_eq!(field, "repeat_tracts.slippage_group_by_read_group");
        assert!(problem.contains("nothing can ever read"), "{problem}");
        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0].read_group = 9;
        });
        assert_eq!(field, "repeat_tracts.substitution_rate_by_stratum");

        // **A duplicate is the same defect wearing the right row count**, so the length check
        // cannot catch it and the ordering check has to.
        let (field, problem) = refused(|file| {
            file.base_quality_calibration.by_read_group[2].read_group = 0;
        });
        assert_eq!(field, "base_quality_calibration.by_read_group");
        assert!(problem.contains("repeats one or skips one"), "{problem}");
    }

    /// **A per-sample table that names a sample the file does not declare is refused**, and the
    /// message names the one that differs rather than printing both lists.
    ///
    /// At the top of the committed input range a cohort is 3,000 samples (`CLAUDE.md`), and an
    /// error carrying two 3,000-name lists is not a message.
    #[test]
    fn a_per_sample_table_that_names_a_stranger_is_refused() {
        let (field, problem) = refused(|file| {
            file.inbreeding.by_sample[0].sample = "a plant this run never had".into();
        });
        assert_eq!(field, "inbreeding.by_sample");
        assert!(problem.contains("a plant this run never had"), "{problem}");
        assert!(
            !problem.contains("Ailsa"),
            "and it does not print the whole cohort: {problem}"
        );

        assert_eq!(
            refused(|file| file.sequencing_batches.by_sample.truncate(1)).0,
            "sequencing_batches.by_sample"
        );

        // **A duplicated name surfaces as the missing one**, which is the more useful half: a
        // reader told "TS-1 appears twice" still has to work out whose row was overwritten, where
        // one told the name that vanished can go and look at it.
        let (_, problem) = refused(|file| {
            file.inbreeding.by_sample[1].sample = file.inbreeding.by_sample[0].sample.clone();
        });
        assert!(
            problem.contains("Ailsa"),
            "the name that went missing is the one to report: {problem}"
        );
    }

    /// **A zero multiplier is refused and a multiplier above one is not.**
    ///
    /// Above one is the ordinary case — it says the instrument was optimistic about its own
    /// qualities. Zero multiplies every read's error probability to nothing, which charges the
    /// whole library the error floor from a number that says the fit found no errors at all.
    #[test]
    fn a_calibration_that_would_charge_every_read_the_error_floor_is_refused() {
        let (field, problem) = refused(|file| {
            file.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .value = 0.0;
        });
        assert!(field.ends_with("error_probability_multiplier"), "{field}");
        assert!(problem.contains("error floor"), "{problem}");
        accepted(|file| {
            file.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .value = 4.5;
        });
    }

    /// **The uncontaminated run written longhand is refused** — the first half of what the shape
    /// used to accept.
    ///
    /// Read literally, a table whose every row is unmeasured takes the read likelihood's mixture
    /// path with every fraction zero, where absence takes its plain formula. The two are different
    /// computations over the same reads (spec §5, first row).
    #[test]
    fn a_contamination_table_nobody_was_measured_in_is_refused() {
        let (field, problem) = refused(|file| {
            for row in &mut file
                .contamination
                .as_mut()
                .expect("the fixture has a table")
                .by_read_group
            {
                row.measurement = None;
            }
        });
        assert_eq!(field, "contamination");
        assert!(problem.contains("longhand"), "{problem}");

        assert_eq!(
            refused(|file| file
                .contamination
                .as_mut()
                .expect("a table")
                .by_read_group
                .clear())
            .0,
            "contamination.by_read_group"
        );

        // And absence itself is the legitimate way to say it.
        accepted(|file| file.contamination = None);
    }

    /// **A measurement with no evidence behind it is refused** — the second half.
    ///
    /// Both counts zero is how the *in-memory* view spells *not measured*; this file spells it by
    /// having no `measurement` key. Carried across, it reads back as *measured and found clean*,
    /// and only the counts can tell those two apart.
    #[test]
    fn a_measurement_with_no_evidence_behind_it_is_refused() {
        // **One count zero is refused as well as both**, because the predicate that decides this
        // in memory is a conjunction: `ContaminationView::was_measured` is `markers > 0 && reads
        // > 0`, so its negation is a disjunction. This row is the worse of the two — it says
        // *measured, 3.1 in 100* in the file and reads back as never measured at all.
        let (field, problem) = refused(|file| {
            file.contamination.as_mut().expect("a table").by_read_group[0].measurement =
                Some(ContaminationMeasurement {
                    share_of_reads_from_another_sample: 0.031,
                    markers_with_reads: 0,
                    reads_on_markers: 90_233,
                    fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                });
        });
        assert!(field.ends_with(".measurement"), "{field}");
        assert!(problem.contains("markers_with_reads 0"), "{problem}");

        // **And a share of exactly one is refused here rather than panicking later.**
        // `FrozenContamination::new` asserts a half-open `[0, 1)` — "a whole library of another
        // individual's DNA is not a sample of this one" — several frames after this point, in a
        // message about a read group rather than about a file.
        let (field, problem) = refused(|file| {
            file.contamination.as_mut().expect("a table").by_read_group[0]
                .measurement
                .as_mut()
                .expect("read group 0 was measured")
                .share_of_reads_from_another_sample = 1.0;
        });
        assert!(
            field.ends_with("share_of_reads_from_another_sample"),
            "{field}"
        );
        assert!(problem.contains("below one"), "{problem}");

        let (field, problem) = refused(|file| {
            file.contamination.as_mut().expect("a table").by_read_group[0].measurement =
                Some(ContaminationMeasurement {
                    share_of_reads_from_another_sample: 0.0,
                    markers_with_reads: 0,
                    reads_on_markers: 0,
                    fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                });
        });
        assert!(field.ends_with(".measurement"), "{field}");
        assert!(
            problem.contains("evidence of not being measured"),
            "{problem}"
        );

        // **Measured and found clean stays legitimate**, which is the distinction being kept: a
        // zero share with counts behind it is a real answer about a real library.
        accepted(|file| {
            file.contamination.as_mut().expect("a table").by_read_group[0].measurement =
                Some(ContaminationMeasurement {
                    share_of_reads_from_another_sample: 0.0,
                    markers_with_reads: 2903,
                    reads_on_markers: 64118,
                    fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                });
        });
    }

    #[test]
    fn an_inbreeding_coefficient_outside_its_range_is_refused() {
        for value in [1.7, 1.0, -0.1] {
            let (field, problem) =
                refused(|file| file.inbreeding.by_sample[0].inbreeding_coefficient.value = value);
            assert!(field.ends_with("inbreeding_coefficient"), "{field}");
            assert!(
                problem.contains("strictly below one"),
                "at {value}: {problem}"
            );
        }
        accepted(|file| file.inbreeding.by_sample[0].inbreeding_coefficient.value = 0.999);
    }

    /// **A cohort with no variation at all is a real answer and is not refused.**
    ///
    /// Spec §3.6: an alternative total of exactly zero is a fully invariant cohort, and the
    /// flooring belongs to the per-locus expansion rather than here. A `validate` that required
    /// both concentrations positive would refuse it.
    #[test]
    fn the_priors_two_concentrations_are_checked_without_flooring_the_invariant_cohort() {
        accepted(|file| {
            file.ordinary_site_prior.alternative_concentration_total = 0.0;
            file.ordinary_site_prior.rung = SeedRung::ZeroDiversity;
        });
        assert!(
            refused(|file| file.ordinary_site_prior.alternative_concentration_total = -1.0)
                .1
                .contains("not negative")
        );
        assert!(
            refused(|file| file.ordinary_site_prior.reference_concentration = 0.0)
                .0
                .ends_with("reference_concentration")
        );
    }

    #[test]
    fn a_slippage_number_outside_its_range_is_refused() {
        let (field, _) = refused(|file| {
            file.repeat_tracts.slippage_by_stratum_and_group[0].share_of_reads_that_slip = 1.4
        });
        assert!(field.ends_with("share_of_reads_that_slip"), "{field}");
        let (field, _) = refused(|file| {
            file.repeat_tracts.slippage_by_stratum_and_group[0].shorter_share = -0.2
        });
        assert!(field.ends_with("shorter_share"), "{field}");
        let (field, problem) =
            refused(|file| file.repeat_tracts.slippage_by_stratum_and_group[0].fall_off = -0.1);
        assert!(field.ends_with("fall_off"), "{field}");
        assert!(problem.contains("not negative"), "{problem}");
    }

    /// **A number that is not a number is refused wherever the walk reaches one**, including in
    /// the curves, which are provenance rather than terms in a score.
    ///
    /// A `NaN` in a curve changes no genotype — the curve records where a slippage number came
    /// from, not what it is — but the writer can emit `nan` and the reader takes it back, so
    /// without this a run would report `nan` as the slope of the line one of its numbers came off.
    /// **These were unchecked until review**, which is how a claim of "every float" gets made about
    /// a walk that reached the row's three numbers and stopped.
    #[test]
    fn a_curve_carrying_something_that_is_not_a_number_is_refused() {
        let (field, problem) = refused(|file| {
            if let LevelSmoothing::Blend { curve, .. } = &mut file
                .repeat_tracts
                .slippage_by_stratum_and_group[THE_ROW_WHOSE_SLIP_SHARE_BLENDS]
                .share_of_reads_that_slip_origin
                .smoothing
            {
                curve.slope = f64::NAN;
            } else {
                panic!("that row blends its slip share");
            }
        });
        assert!(field.ends_with("blend.curve.slope"), "{field}");
        assert!(problem.contains("not a number"), "{problem}");

        let (field, _) = refused(|file| {
            let shares = file.repeat_tracts.slippage_by_stratum_and_group
                [THE_ROW_WHOSE_SLIP_SHARE_BLENDS]
                .shorter_share_and_fall_off_origin
                .as_mut()
                .expect("that row records its shares' origin");
            if let ShareSmoothing::ThisPeriodsCurve { curve, .. } = &mut shares.fall_off_smoothing {
                curve.centre_repeats = f64::INFINITY;
            } else {
                panic!("that row takes its fall-off from a period curve");
            }
        });
        assert!(field.ends_with("curve.centre_repeats"), "{field}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.slippage_by_stratum_and_group[THE_ROW_WHOSE_SLIP_SHARE_BLENDS]
                .share_of_reads_that_slip_origin
                .expected_slipped_reads = Some(f64::NAN);
        });
        assert!(
            field.ends_with("share_of_reads_that_slip_origin.expected_slipped_reads"),
            "{field}"
        );
    }

    /// **`fall_off` is not bounded above, and that is a decision rather than an oversight.**
    ///
    /// Neither this file's shape nor the fit that produces it documents an upper bound for it.
    /// A value above one would say a two-repeat slip is likelier than a one-repeat slip, which is
    /// implausible — but nothing has established it is impossible, and refusing a fit nobody has
    /// bounded would reject real data to enforce a guess.
    #[test]
    fn a_fall_off_above_one_is_accepted_because_nothing_documents_a_bound() {
        accepted(|file| file.repeat_tracts.slippage_by_stratum_and_group[0].fall_off = 1.5);
    }

    #[test]
    fn a_curve_weight_outside_zero_to_one_is_refused() {
        let (field, _) = refused(|file| {
            if let LevelSmoothing::Blend { curve_weight, .. } = &mut file
                .repeat_tracts
                .slippage_by_stratum_and_group[THE_ROW_WHOSE_SLIP_SHARE_BLENDS]
                .share_of_reads_that_slip_origin
                .smoothing
            {
                *curve_weight = 1.9;
            } else {
                panic!("that row blends its slip share");
            }
        });
        assert!(field.ends_with("blend.curve_weight"), "{field}");

        let (field, _) = refused(|file| {
            let shares = file.repeat_tracts.slippage_by_stratum_and_group
                [THE_ROW_WHOSE_SHARES_BLEND]
                .shorter_share_and_fall_off_origin
                .as_mut()
                .expect("that row records its shares' origin");
            if let ShareSmoothing::Blend { curve_weight, .. } = &mut shares.shorter_share_smoothing
            {
                *curve_weight = -0.5;
            } else {
                panic!("the fixture's third row blends its shorter share");
            }
        });
        assert!(field.ends_with("blend.curve_weight"), "{field}");
    }

    /// **A length spectrum that is not a distribution over an odd span is refused.**
    ///
    /// The array runs from `-span` to `+span` in whole repeats, so the middle entry *is* the
    /// reference length. An even count has no middle; fewer than three cannot express a slip in
    /// both directions; and shares that do not sum to one are not a distribution — which is what
    /// a person who edits one share and does not rebalance the rest produces.
    #[test]
    fn a_length_spectrum_that_is_not_a_distribution_is_refused() {
        // **Even and long enough**, which is what tests the parity rule rather than the length
        // rule. A two-share spectrum is refused either way — it is also below three — so a
        // fixture of two would leave a reader believing the parity check is exercised when
        // deleting it changes nothing. Measured: with `[0.5, 0.5]` here, removing the parity
        // clause left all 100 tests green.
        let (field, problem) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset =
                vec![0.25, 0.25, 0.25, 0.25];
        });
        assert!(field.ends_with("shares_by_repeat_offset"), "{field}");
        assert!(problem.contains("odd and at least three"), "{problem}");

        // And the short case, which is the other half of the same rule.
        let (_, problem) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset =
                vec![0.5, 0.5];
        });
        assert!(problem.contains("odd and at least three"), "{problem}");

        let (_, problem) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset = vec![1.0];
        });
        assert!(problem.contains("odd and at least three"), "{problem}");

        let (_, problem) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset =
                vec![0.1, 0.9, 0.1];
        });
        assert!(problem.contains("sums to"), "{problem}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.length_spectrum_by_period[0].shares_by_repeat_offset =
                vec![0.2, 0.2, 0.2];
        });
        assert!(
            field.starts_with("repeat_tracts.length_spectrum_by_period"),
            "{field}"
        );

        assert!(
            refused(|file| file.repeat_tracts.length_spectrum_by_stratum[0].concentration = 0.0)
                .0
                .ends_with("concentration")
        );
    }

    /// **A normalised spectrum sits inside the tolerance and a hand edit outside it, with six
    /// orders of magnitude of daylight on each side.**
    ///
    /// The fixture's own three-offset spectra sum to **exactly** one, so they do not exercise this
    /// at all — an earlier version of this comment claimed they missed by one unit in the last
    /// place, and that was recalled rather than measured. What does exercise it is a spectrum over
    /// many offsets, where the fit's own division and addition leave a residue.
    #[test]
    fn the_tolerance_admits_a_normalised_spectrum_and_not_an_edit() {
        fn normalised(offsets: usize) -> Vec<f64> {
            let middle = offsets / 2;
            let raw: Vec<f64> = (0..offsets)
                .map(|at| 0.6_f64.powi((at as i32 - middle as i32).abs()))
                .collect();
            let total: f64 = raw.iter().sum();
            raw.iter().map(|share| share / total).collect()
        }

        let mut worst: f64 = 0.0;
        for offsets in [3, 5, 9, 21, 41, 101] {
            let total: f64 = normalised(offsets).iter().sum();
            worst = worst.max((total - 1.0).abs());
        }
        assert!(
            worst < SHARES_MAY_MISS_ONE_BY / 1e6,
            "a normalised spectrum misses one by at most {worst:e}, which should sit a million \
             times inside the tolerance of {SHARES_MAY_MISS_ONE_BY:e}"
        );

        // The smallest edit a person makes to a share — a hundredth — is as far outside.
        let edited: f64 = [0.1, 0.81, 0.1].iter().sum();
        assert!(
            (edited - 1.0).abs() > SHARES_MAY_MISS_ONE_BY * 1e6,
            "and a hundredth moves the sum by {:e}",
            (edited - 1.0).abs()
        );
    }

    #[test]
    fn a_substitution_rate_that_is_not_a_probability_is_refused() {
        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0]
                .rate
                .value = 1.2
        });
        assert!(field.ends_with(".rate"), "{field}");
        assert!(
            refused(|file| file.stated_constants.repeat_tract_outlier_weight.value = 2.0)
                .0
                .ends_with("repeat_tract_outlier_weight")
        );
        // **And both ends of its range, which no other share in this file refuses.** The scoring
        // row asserts `0 < weight < 1`, so a zero or a one that got past here would panic at
        // whichever repeat tract came first, naming a locus rather than the file. It is also the
        // one number spec §3.8 invites a person to edit, so the ends are what a typo reaches.
        for edge in [0.0, 1.0] {
            let (field, problem) =
                refused(|file| file.stated_constants.repeat_tract_outlier_weight.value = edge);
            assert!(field.ends_with("repeat_tract_outlier_weight"), "{field}");
            assert!(
                problem.contains("strictly inside zero and one"),
                "{problem}"
            );
        }
        // A weight the scorer would accept is accepted.
        accepted(|file| {
            file.stated_constants.repeat_tract_outlier_weight = WarrantedValue {
                value: 0.5,
                warrant: Warrant::Supplied,
                observations: None,
            };
        });
    }

    /// **The outlier weight is the one key held to two of the four warrants, and the state
    /// worth catching is an edited number under a `defaulted` label.**
    ///
    /// Nothing fits this number (spec §3.8), so `fitted_here` and `borrowed` are claims about
    /// it no run could make and the shape it projects onto has nowhere to put them. The third
    /// refusal is the one a person reaches: spec §1.2 goal 3 invites them to copy the file
    /// their run wrote and change one line, and the line they change is the number — leaving
    /// `warrant = "defaulted"`, which then says the run inherited a constant it did not.
    #[test]
    fn an_outlier_weight_whose_warrant_no_run_could_mean_is_refused() {
        for unfittable in [Warrant::FittedHere, Warrant::Borrowed] {
            let (field, problem) = refused(|file| {
                file.stated_constants.repeat_tract_outlier_weight.warrant = unfittable;
            });
            assert!(field.ends_with("repeat_tract_outlier_weight"), "{field}");
            assert!(problem.contains("nothing fits this number"), "{problem}");
        }

        let (field, problem) = refused(|file| {
            file.stated_constants.repeat_tract_outlier_weight.value = 0.05;
        });
        assert!(field.ends_with("repeat_tract_outlier_weight"), "{field}");
        assert!(
            problem.contains("change the warrant beside it to `supplied`"),
            "an edited number under a defaulted warrant is told what to change: {problem}"
        );

        // The two states a run can mean are both accepted: the constant it inherited, and a
        // number somebody wrote in.
        accepted(|file| {
            file.stated_constants.repeat_tract_outlier_weight = WarrantedValue {
                value: DEFAULT_OUTLIER_WEIGHT,
                warrant: Warrant::Defaulted,
                observations: None,
            };
        });
        accepted(|file| {
            file.stated_constants.repeat_tract_outlier_weight = WarrantedValue {
                value: 0.05,
                warrant: Warrant::Supplied,
                observations: None,
            };
        });
        assert!(
            refused(|file| file
                .repeat_tracts
                .fallback_length_spectrum_concentration
                .value = 0.0)
            .0
            .ends_with("fallback_length_spectrum_concentration")
        );
    }

    /// **A count at exactly the largest integer TOML holds is the writer's saturation marker.**
    ///
    /// Step B2 measured that a `u64` above `2^63 − 1` gave three different answers from three
    /// readers, so every integer this writer emits saturates rather than emitting one nobody
    /// agrees on. A count arriving back at exactly that value is a number the writer knows it
    /// lost; read as evidence it would put 9.2 quintillion reads into a run's report.
    #[test]
    fn an_evidence_count_at_the_saturation_marker_is_refused() {
        let saturated = i64::MAX as u64;
        let (field, problem) = refused(|file| {
            file.inbreeding.by_sample[0]
                .inbreeding_coefficient
                .observations = Some(EvidenceCount::CoveredPositions(saturated));
        });
        assert!(field.ends_with("observations.covered_positions"), "{field}");
        assert!(problem.contains("saturates"), "{problem}");

        let (field, _) = refused(|file| {
            file.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .observations = Some(EvidenceCount::Reads(saturated));
        });
        assert!(field.ends_with("observations.reads"), "{field}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0]
                .rate
                .observations = Some(EvidenceCount::BasesCompared(saturated));
        });
        assert!(field.ends_with("observations.bases_compared"), "{field}");

        // One below it is a legitimate, if enormous, count.
        accepted(|file| {
            file.inbreeding.by_sample[0]
                .inbreeding_coefficient
                .observations = Some(EvidenceCount::CoveredPositions(saturated - 1));
        });
    }

    /// **A contamination measurement's counts carry the same marker**, and they are not
    /// `EvidenceCount`s, so they need their own check.
    #[test]
    fn a_contamination_count_at_the_saturation_marker_is_refused() {
        let (field, problem) = refused(|file| {
            file.contamination.as_mut().expect("a table").by_read_group[0]
                .measurement
                .as_mut()
                .expect("read group 0 was measured")
                .reads_on_markers = i64::MAX as u64;
        });
        assert!(field.ends_with("measurement.reads_on_markers"), "{field}");
        assert!(problem.contains("saturation marker"), "{problem}");
    }

    /// **A value that is not a number is reported as not a number**, not as out of range.
    ///
    /// `NaN` fails every comparison, so a range test written as `!(0.0..=1.0).contains(&value)`
    /// calls it *outside [0, 1]* — which sends a reader looking through their file for a number
    /// that is out of range, and every number they find is in range.
    #[test]
    fn a_value_that_is_not_a_number_says_so_rather_than_reporting_a_range() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let (_, problem) =
                refused(|file| file.inbreeding.by_sample[0].inbreeding_coefficient.value = value);
            assert!(
                problem.contains("not a number a score can be computed from"),
                "at {value}: {problem}"
            );
            assert!(
                !problem.contains("strictly below one"),
                "and it does not claim a range instead: {problem}"
            );
        }
    }

    /// **Every refusal names a path that occurs in the file it refuses.**
    ///
    /// A path a reader cannot find is worse than no path: it sends them looking. This walks the
    /// refusals whose field is a literal key rather than an indexed row and checks each appears
    /// in the written text.
    #[test]
    fn every_refusal_names_a_key_the_file_actually_contains() {
        let text = a_file_using_every_shape().to_toml();
        for (field, _) in [
            refused(|file| file.ploidy = 0),
            refused(|file| file.fitted_from.samples.clear()),
            refused(|file| file.fitted_from.read_groups.clear()),
            refused(|file| file.inbreeding.by_sample.clear()),
            refused(|file| file.ordinary_site_prior.reference_concentration = 0.0),
            refused(|file| {
                file.repeat_tracts
                    .fallback_length_spectrum_concentration
                    .value = 0.0
            }),
            refused(|file| file.stated_constants.repeat_tract_outlier_weight.value = 2.0),
            refused(|file| file.format_version = 9),
            // **The nested paths, which are what an earlier version of this test could not
            // see.** It compared only the last segment, so a refusal that named the row's
            // `shorter_share` where the key is `shorter_share_smoothing` passed on the strength
            // of its final `curve_weight` — and `shorter_share` is a real sibling key of that
            // same row, holding a perfectly good number. A path that names a healthy key is
            // worse than one that names nothing, because the reader stops there.
            refused(|file| {
                let shares = file.repeat_tracts.slippage_by_stratum_and_group
                    [THE_ROW_WHOSE_SHARES_BLEND]
                    .shorter_share_and_fall_off_origin
                    .as_mut()
                    .expect("that row records its shares' origin");
                if let ShareSmoothing::Blend { curve_weight, .. } =
                    &mut shares.shorter_share_smoothing
                {
                    *curve_weight = -0.5;
                } else {
                    panic!("that row blends its shorter share");
                }
            }),
            refused(|file| {
                if let LevelSmoothing::Blend { curve_weight, .. } = &mut file
                    .repeat_tracts
                    .slippage_by_stratum_and_group[THE_ROW_WHOSE_SLIP_SHARE_BLENDS]
                    .share_of_reads_that_slip_origin
                    .smoothing
                {
                    *curve_weight = 1.9;
                } else {
                    panic!("that row blends its slip share");
                }
            }),
        ] {
            // **Every segment, not just the last.** A dotted path is only as good as its
            // weakest link, and the segments that carry a row index are checked as their key
            // name — `by_read_group[0]` against `by_read_group` — because the file writes its
            // rows as inline tables inside an array and has no literal index anywhere.
            for segment in field.split('.') {
                let key = segment.split(['[', '(']).next().unwrap_or(segment);
                if key.is_empty() {
                    continue;
                }
                assert!(
                    text.contains(key),
                    "the refusal names {field}, whose segment {key} is not a key in the file"
                );
            }
        }
    }

    /// **A warrant of `supplied` is not itself grounds for refusal.**
    ///
    /// Spec §1.2 goal 3 is a person editing one line, and the file's own header tells them to mark
    /// what they changed as `supplied`. A `validate` that treated a hand-typed number as suspect
    /// would refuse the file the instructions produce.
    #[test]
    fn a_hand_edited_value_marked_supplied_is_accepted() {
        accepted(|file| {
            file.base_quality_calibration.by_read_group[0].error_probability_multiplier =
                WarrantedValue {
                    value: 1.4,
                    warrant: Warrant::Supplied,
                    observations: None,
                };
        });
    }
}
