//! **The file's shape, as the text a run writes** — spec §4's layout, emitted by hand.
//!
//! # Why this is not `serde`'s serializer
//!
//! **Spec §4 chose TOML for the comments** — each defaulted value carries where its default came
//! from, and the softest number in the file is a slippage rate fitted from one human genome
//! standing in for every chemistry a run might use. No serde serializer can emit a comment, so
//! the writer was always going to be hand-written; step B3 is where the comments land, and this
//! step is the text they attach to.
//!
//! Two things follow from writing it by hand, and both are why the derived serializer stays for
//! the tests rather than for the artefact:
//!
//! - **the layout is chosen rather than inherited.** `serde`'s TOML output puts every row of
//!   every table under its own `[[header]]`, emits a struct's table-valued fields after its
//!   scalar ones whatever the declared order, and gives a section whose fields are all tables no
//!   header at all — so the largest section of the file opens unnamed. None of that is the
//!   file's design; it is one crate's rendering of it.
//! - **the key surface needs its own golden file.** `testdata/every_shape.toml` pins the names
//!   under serde's writer and cannot see this one. `testdata/every_shape_as_written.toml` is
//!   this writer's, from the same fixture, so the two can be read side by side.
//!
//! # The layout, and where it departs from spec §4
//!
//! **One row a line, as an inline table**, for every table keyed by an axis — the read groups,
//! the samples, the batches, the `(stratum × slippage group)` pairs, the length spectra and the
//! substitution rates. Each row names its own axis value, so a line can be read on its own and a
//! person editing one sample's number can find that sample's line.
//!
//! **§4 suggests the numeric rows of §3.7 be "arrays of arrays rather than arrays of tables",
//! and this writer emits inline tables there too.** Two reasons, and the plan's own preconditions
//! put the choice here — "the *spelling* is the coder's: key names, nesting, where a table
//! starts, **whether a row is inline**":
//!
//! - **a slippage row is not flat.** It carries where its level and its two shares came from, and
//!   a number off a curve carries the whole curve — seven fields for a slippage curve, ten for a
//!   share curve, **four braces inside the row's own fields**. Flattened to a bare array, the row
//!   becomes a sequence of numbers whose meaning is a column legend the reader has to hold, and
//!   the nesting has nowhere to go.
//! - **the reason §4 gave for the suggestion no longer applies.** It says arrays of arrays
//!   "neither needs a custom encoder" — and this *is* the custom encoder, written because the
//!   comments demanded one.
//!
//! **Arrays are one element a line where they are a table, and inline where they are a number.**
//! A length spectrum's shares are three to a handful of numbers describing one thing and stay on
//! their row; `samples` and every `by_*` table get a line each, which is what spec §9 prices.
//!
//! # One naming rule, stated because it is load-bearing
//!
//! **A helper whose name starts with an article returns a piece of TOML; one that does not writes
//! into the document.** `a_toml_string`, `an_inline_table` and `a_slippage_row` hand back a
//! `String`; `section`, `scalar` and `one_a_line` take `out: &mut String` and append. The rule
//! holds for every helper here, and it is what lets a reader see at the call site whether a line
//! has been emitted yet.
//!
//! *(The sibling projection, `from_run_parameters.rs`, spends its articles differently — it has no
//! append-or-return distinction to signal, and names its builders after what they produce. Neither
//! file's rule is the crate's; each states its own.)*

use std::fmt::Write as _;

use crate::ng::alignment::StutterModel;
use crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE;

use super::{
    BaseQualityCalibrationRow, CensusTerm, ContaminationMeasurement, ContaminationRow, CurveReach,
    EvidenceCount, GroupOfNumbers, InbreedingRow, LevelOrigin, LevelSmoothing, ParametersFile,
    PeriodLengthSpectrumRow, ReadGroupBatchRow, ReadGroupRow, SampleBatchRow, SeedRung, ShareCurve,
    ShareCurveRung, ShareShape, ShareSmoothing, SharesOrigin, SlippageCurve, SlippageGroupRow,
    SlippageRow, StratumLengthSpectrumRow, SubstitutionRateRow, Warrant, WarrantedValue,
};

impl ParametersFile {
    /// **The file, as text.**
    ///
    /// Parses back to an equal value through the shape's own `Deserialize`, which is what
    /// `the_written_text_reads_back_as_the_same_file` holds — and what makes this writer and
    /// step C's reader two halves of one round trip rather than two independent readings of
    /// spec §4.
    ///
    /// **An uncontaminated run writes no `[contamination]` section**, which is spec §5's first
    /// state and the only one expressed by a whole section going missing. Every other absence is
    /// a key that is not written.
    #[must_use]
    pub fn to_toml(&self) -> String {
        let mut out = String::new();

        // **The first thing a reader meets, and it is derived from the numbers below.** A
        // defaults run's file and a fit's are the same shape by design (spec §7), so nothing on
        // the page told them apart: measured on this module's fixtures, the two opened with the
        // same 39 lines of prose and the first thing in either that said which run it was is the
        // `warrant` on line 105.
        a_derived_note(&mut out, &how_much_was_fitted(self));
        // A `#` on its own, so the derived paragraph and the written one below do not run
        // together into one block of prose.
        a_derived_note(&mut out, &[String::new()]);

        note(
            &mut out,
            &[
                "Every number this run scored its reads under, and what each one rests on.",
                "",
                "A number that could be fitted carries a `warrant`: fitted_here, borrowed, supplied or defaulted. **If you edit one, change its warrant to \"supplied\" and delete its `observations`** — otherwise this file says a number you typed was measured, and the run that reads it will report it that way. A `supplied` number that still carries `observations` came that way from another run's file, and those counts are that run's.",
                "",
                "**Two keys do not take every warrant.** `repeat_tracts.fallback_length_spectrum_concentration` is `fitted_here` only where this file holds a fitted stratum spectrum for it to be the median of, and `defaulted` only at the built-in constant; `stated_constants.repeat_tract_outlier_weight` is `defaulted` only at the built-in constant. Both take `supplied` freely, which is what you write when you change one. Anything else is refused, and says so.",
                "",
                "**What that checking reaches, and what it cannot.** It reaches those two keys and nowhere else. It cannot catch a number you changed that still says `fitted_here` or `borrowed`: nothing in this file can tell your value from a fitted one, so the run will report it as measured, with the old `observations` count still beside it. And on every other key a `defaulted` warrant is checked against nothing — including the inbreeding coefficient, whose built-in number is 0.0. The base-quality multiplier is not one of them: a `defaulted` multiplier is **not** fixed at 1.0, and the note beside `base_quality_calibration` says what it is instead. A `defaulted` substitution rate is the one to avoid writing by hand, and the reason is worth stating because a note further down looks like it contradicts this: the caller **does** default that number, at the tract and for the cells that need it, but it never writes one as a row here — so a `defaulted` warrant on a row of `substitution_rate_by_stratum` is a claim no build makes.",
                "",
                "The slippage numbers, the prior's two concentrations and the length spectrum rows carry no warrant — they say where they came from another way, and there is nowhere in them to record that you changed one. Note such an edit elsewhere.",
                "",
                "An absent key is not a zero. A missing section, a missing row and a missing key each mean the thing was not measured. A zero is not automatically the opposite: a zero **under a `fitted_here` or `borrowed` warrant** was measured and found to be zero, and a zero under any other warrant is a stated number — the inbreeding coefficient a run takes when nobody said is 0.0, `defaulted`. The sections below say which is which where it matters.",
            ],
        );

        writeln!(out, "format_version = {}", self.format_version).expect("a string never fails");
        writeln!(out, "ploidy = {}", self.ploidy).expect("a string never fails");

        section(&mut out, "fitted_from");
        // **The section's own heading says *fitted*, and on a defaults run nothing was.** A
        // geneticist reading such a file read `[fitted_from]` as a claim that the numbers under
        // it had been. The key name is the format and does not move; what it means is said here,
        // derived so that the two cannot come apart.
        a_derived_note(&mut out, &what_the_bindings_are_for(self));
        note(
            &mut out,
            &[
                "",
                "The MD5 of the reference this run ran against: every contig's bases, uppercased, run together in the order the FASTA holds them. So soft-masking and line width do not change it and contig order does. `[fitted_from.census]` below has a `reference digest` line of its own — that is this same reference seen from the evidence's side, and the key that turns a run away from the wrong reference is this one.",
            ],
        );
        scalar(
            &mut out,
            "reference_digest",
            &a_toml_string(&self.fitted_from.reference_digest),
        );
        one_a_line(
            &mut out,
            "samples",
            self.fitted_from.samples.iter().map(|s| a_toml_string(s)),
        );
        one_a_line(
            &mut out,
            "read_groups",
            self.fitted_from.read_groups.iter().map(a_read_group_row),
        );

        section(&mut out, "fitted_from.census");
        if self.fitted_from.census.terms.is_empty() {
            a_derived_note(&mut out, &this_run_had_no_census(self));
        }
        note(
            &mut out,
            &[
                "Which store of evidence these numbers were fitted from — one line for each of the things two runs' evidence has to agree on before it can be pooled, each written as a digest of what that thing was rather than as the thing itself, because several of them are whole tables and settings blocks.",
                "",
                "A run whose own evidence disagrees on any one of these lines is the demotion the section above describes, and it demotes **every** number in the file rather than the ones that line touches — the numbers were fitted together out of this one store of evidence, so a disagreement about the evidence disqualifies all of them. The run says which line differed.",
                "",
                "Do not edit these to make a run match. Nothing checks a digest against what it claims to digest, so an edit buys a run that reports numbers as fitted from its own data when they were fitted from somebody else's — and there is no way to work out the right value by hand in any case. The same goes for `reference_digest`, `samples` and `read_groups` above, where a mismatch is turned away rather than demoted.",
            ],
        );
        one_a_line(
            &mut out,
            "terms",
            self.fitted_from.census.terms.iter().map(a_census_term),
        );

        section(&mut out, "base_quality_calibration");
        note(
            &mut out,
            &[
                "What each read's own reported error probability is multiplied by, per read group. Above one says the instrument was optimistic and the reads are worse than it claimed; below one says they are better; one leaves the qualities exactly as reported. It is not a multiplier on the Phred score, which moves the other way.",
                "",
                "**A `defaulted` multiplier is not 1.0 in general**, and the rows below say `not calibrated` where it applies. A read group no error rate could be fitted for is *charged the caller's stated rate of one error in a thousand* rather than taken at the quality it reported — so its multiplier is that stated rate divided by this library's own mean reported error, and rises above one wherever the instrument claimed better than one in a thousand. A library averaging Q36 gets 4.0: every base of it is scored four times likelier to be wrong than the file it came in said. **The exception is a run that fitted nothing at all**, which read nothing to take a mean over and writes 1.0.",
            ],
        );
        one_a_line_with_notes(
            &mut out,
            "by_read_group",
            self.base_quality_calibration
                .by_read_group
                .iter()
                .map(|row| {
                    (
                        a_calibration_row(row),
                        where_it_came_from(
                            &row.error_probability_multiplier,
                            origins::CALIBRATION_MULTIPLIER,
                        ),
                    )
                }),
        );

        // **Absent means uncontaminated**, and the section is where that is said. A table of
        // zeros would claim every library was measured and found clean, which nothing measured.
        if let Some(contamination) = &self.contamination {
            section(&mut out, "contamination");
            note(
                &mut out,
                &[
                    "How much of each read group's reads came from somebody else — one row a lane, because two lanes of one library can differ: index hopping happens on a flowcell, not in a tube. Three states, three different claims:",
                    "  - this whole section absent  -> nobody identified any contamination",
                    "  - a row with no `measurement` -> this lane could not be measured",
                    "  - a zero share with non-zero counts -> measured, and found clean",
                    "To stop correcting one lane, delete its `measurement = { ... }` and leave the row; a library sequenced over several lanes has a row for each. Setting that share to zero says something else: that it was measured and found clean.",
                ],
            );
            one_a_line(
                &mut out,
                "by_read_group",
                contamination.by_read_group.iter().map(a_contamination_row),
            );
        } else {
            // **A blank line first**, so the note does not read as a remark about the calibration
            // table it follows: a reader scanning section by section met it under
            // `[base_quality_calibration]` and either mis-attributed it or skipped it looking for
            // a `[contamination]` heading that is not there. TOML has no way to head an absent
            // section, so the separation is all this can do.
            out.push('\n');
            // **⚑ The explanation of what an absent section means used to be inside the section.**
            // So the one reader who most needed it — the one holding a file with no
            // `[contamination]` — was the only reader who never saw it, and the word
            // *contamination* appeared nowhere in a defaults run's whole file. Found by a
            // geneticist reading one (step E3); the fix is that absence says its own name.
            note(
                &mut out,
                &[
                    "**Contamination: no `[contamination]` section, which means nobody identified any.** That is not the same as *measured and found clean* — nothing was measured. Your reads are scored as though none of them came from another individual, which is the read likelihood's plain formula rather than a correction of size zero. A run that fitted contamination writes a section here with one row a lane.",
                ],
            );
        }

        section(&mut out, "sequencing_batches");
        note(
            &mut out,
            &[
                "Who was sequenced beside whom — the population a contaminating read is drawn from. `batching_was_declared = false` means nobody said, so everything went in one batch. A declared batching that happens to have one batch writes identical rows, and this flag is the only thing that tells those two apart.",
                "",
                "**A sample's read groups all go in one batch**, and it is the batch its own row names — the two tables below say the same thing about a sample twice, and a file in which they disagree is refused. The grain is a lane everywhere else in this file precisely because two lanes of one plant can differ; here they cannot, because a contaminating read is drawn against one set of neighbours and a sample has one genotype to draw.",
            ],
        );
        scalar(
            &mut out,
            "batching_was_declared",
            if self.sequencing_batches.batching_was_declared {
                "true"
            } else {
                "false"
            },
        );
        one_a_line(
            &mut out,
            "by_read_group",
            self.sequencing_batches
                .by_read_group
                .iter()
                .map(a_read_group_batch_row),
        );
        one_a_line(
            &mut out,
            "by_sample",
            self.sequencing_batches
                .by_sample
                .iter()
                .map(a_sample_batch_row),
        );

        section(&mut out, "inbreeding");
        note(
            &mut out,
            &[
                "How inbred each sample is, as a fraction in [0, 1) — one row a sample, counted over the reference positions that sample's reads covered.",
            ],
        );
        one_a_line_with_notes(
            &mut out,
            "by_sample",
            self.inbreeding.by_sample.iter().map(|row| {
                (
                    an_inbreeding_row(row),
                    where_it_came_from(
                        &row.inbreeding_coefficient,
                        origins::INBREEDING_COEFFICIENT,
                    ),
                )
            }),
        );

        section(&mut out, "ordinary_site_prior");
        note(
            &mut out,
            &[
                "What the SNP and indel prior starts from at an ordinary position: how much belief the reference allele carries, and how much is shared out across whatever alternative alleles a position turns out to have. `rung` says which measurement the pair came off: a fitted population curve (`fitted_curve`), the neutral shape at a fitted heterozygosity (`neutral_shape`), a cohort with no variation at all (`zero_diversity`), or a stated heterozygosity taken from human data (`stated_heterozygosity`), which is the one that rests on nothing this run measured.",
            ],
        );
        scalar(
            &mut out,
            "reference_concentration",
            &a_toml_float(self.ordinary_site_prior.reference_concentration),
        );
        scalar(
            &mut out,
            "alternative_concentration_total",
            &a_toml_float(self.ordinary_site_prior.alternative_concentration_total),
        );
        scalar(
            &mut out,
            "rung",
            &a_toml_string(the_word_for_seed_rung(self.ordinary_site_prior.rung)),
        );

        section(&mut out, "repeat_tracts");
        let tracts = &self.repeat_tracts;
        note(
            &mut out,
            &[
                "Everything about repeat tracts. A **stratum** is a class of tract, and every row keyed by one spells it as `period` — how many bases one repeat unit is — and `reference_repeats` — how many copies of it the reference carries. A **slippage group** is a set of read groups whose reads are taken to slip alike; the run declares it, and `slippage_group_by_read_group` below is that declaration.",
                "",
                "`slippage_by_stratum_and_group` is keyed by a stratum **and** a slippage group, and a triple with no row means that group put no read in that stratum — so one stratum can have a row for one group and none for another. A stratum with no row in `length_spectrum_by_stratum` was never fitted on its own tracts and falls to its period's pooled one, or to the flat shape below. Neither absence is a zero.",
                "",
                "Three numbers a stratum: `share_of_reads_that_slip` — how often a read reports a tract length other than its allele's; `shorter_share` — of the reads that slip, the share showing a shorter tract; `fall_off` — how fast two-repeat slips fall off against one-repeat slips. `expected_slipped_reads` is fractional because it is how many reads the fitted share says slipped, not a count anybody labelled.",
                "",
                "`share_of_reads_that_slip_origin` says where the first of the three came from and `shorter_share_and_fall_off_origin` where the other two did: this stratum's own fit, its period's curve, or a blend of the two, with the curve itself written down so an interpolation can be told from a measurement. A row whose two shares were not fitted here at all has no `shorter_share_and_fall_off_origin` key. Each origin carries its own `expected_slipped_reads` where this stratum fitted a slip share of its own, and neither carries one where the number was taken whole from a curve — so a row showing the same count twice is not a duplicate, and a row showing none fitted nothing of its own.",
                "",
                "The curves under `shorter_share_and_fall_off_origin` also record `curve_fitted_on`, which says what that curve itself was fitted on: this period's own strata (`this_period`), or those same strata where there were too few to score the shape (`this_period_unscored`), or the other periods pooled (`other_periods`), or a stated constant where no period had anything to fit (`built_in_default`). The curve under `share_of_reads_that_slip_origin` is a different fit and has no such key.",
                "",
                "`fallback_length_spectrum_concentration` is what a tract falls back to where neither its own stratum nor its period was fitted: this many chromosomes' worth of belief, spread flat over whatever lengths the tract offers, so a larger number makes the prior harder for the reads to move. A run that fitted any stratum states the median of the concentrations those fits produced and marks it `fitted_here`; a run that fitted none states a built-in constant and marks it `defaulted`. A run handed this file writes back the warrant it found there, because a file's warrants survive being read — demoting them on every read would have the same cohort report different warrants depending on which way it was called. Mark it `supplied` yourself if you change the number. It carries no `observations` in any of those cases — a median over strata is not a measurement with a sample size.",
            ],
        );
        scalar_with_note(
            &mut out,
            "fallback_length_spectrum_concentration",
            &a_warranted_value(&tracts.fallback_length_spectrum_concentration),
            where_it_came_from(
                &tracts.fallback_length_spectrum_concentration,
                origins::FLAT_CONCENTRATION,
            ),
        );
        one_a_line(
            &mut out,
            "slippage_group_by_read_group",
            tracts
                .slippage_group_by_read_group
                .iter()
                .map(a_slippage_group_row),
        );
        // **An empty table and a missing row are two different claims, and the section's own
        // paragraph above covers only the second.** A *missing row* means that slippage group put
        // no read in that stratum. An *empty table* means no stratum was fitted at all, and then
        // every repeat tract of the run falls through to another caller's shipped constants —
        // which is a stronger claim about the reads than "nothing happened here", and the one a
        // geneticist reading a defaults run's file took for the weaker (step E3).
        //
        // **So the note is for the empty table only, and a partially fitted run gets none.** That
        // is not an omission: where *some* strata were fitted, how much of a given tract fell back
        // is a property of that tract and rides on its own record
        // (`TractScoringFits::cells_with_no_fitted_slippage`). A note here saying "some tracts
        // fell back" would be true of almost every run and would tell a reader nothing about
        // theirs.
        if tracts.slippage_by_stratum_and_group.is_empty() {
            for (at, paragraph) in no_stratum_was_fitted().iter().enumerate() {
                if at > 0 {
                    writeln!(out, "#").expect("a string never fails");
                }
                note_lines(&mut out, &wrapped(paragraph, ROOM_AT_THE_MARGIN));
            }
        }
        one_a_line(
            &mut out,
            "slippage_by_stratum_and_group",
            tracts
                .slippage_by_stratum_and_group
                .iter()
                .map(a_slippage_row),
        );
        note(
            &mut out,
            &[
                "`shares_by_repeat_offset` runs from -span to +span in whole repeat units from the **reference** tract length, so the middle entry is the reference length itself, the count is odd, and the shares sum to one. An array read as starting at zero shifts every length this stratum expects by one repeat.",
            ],
        );
        one_a_line(
            &mut out,
            "length_spectrum_by_stratum",
            tracts
                .length_spectrum_by_stratum
                .iter()
                .map(a_stratum_length_spectrum_row),
        );
        one_a_line(
            &mut out,
            "length_spectrum_by_period",
            tracts
                .length_spectrum_by_period
                .iter()
                .map(a_period_length_spectrum_row),
        );
        note(
            &mut out,
            &[
                "How often a base reads wrong inside a tract — per read group as well as per stratum, because that is a property of the chemistry, and per `ploidy` as well, because that is the set of genotypes the fit scored these tracts against. That third key is why a row repeats the number at the top of the file: a cohort called at two ploidies carries a row for each. Counted in bases compared, not reads: a read crossing a tract contributes as many bases as it crosses.",
            ],
        );
        // The same distinction one table down, and the same reason.
        if tracts.substitution_rate_by_stratum.is_empty() {
            note_lines(
                &mut out,
                &wrapped(&no_substitution_rate_was_fitted(), ROOM_AT_THE_MARGIN),
            );
        }
        one_a_line_with_notes(
            &mut out,
            "substitution_rate_by_stratum",
            tracts.substitution_rate_by_stratum.iter().map(|row| {
                (
                    a_substitution_rate_row(row),
                    where_it_came_from(&row.rate, origins::SUBSTITUTION_RATE),
                )
            }),
        );

        section(&mut out, "stated_constants");
        note(
            &mut out,
            &[
                "The numbers no fit produces, written out so that what a run inherited is visible and editable rather than buried in the binary.",
            ],
        );
        scalar_with_note(
            &mut out,
            "repeat_tract_outlier_weight",
            &a_warranted_value(&self.stated_constants.repeat_tract_outlier_weight),
            where_it_came_from(
                &self.stated_constants.repeat_tract_outlier_weight,
                origins::OUTLIER_WEIGHT,
            ),
        );

        out
    }
}

// ---------------------------------------------------------------------
// The three shapes every line of the file is made of
// ---------------------------------------------------------------------

/// A blank line and a `[header]` — every section opens one, including the ones whose every field
/// is a table, which `serde`'s writer leaves unnamed.
fn section(out: &mut String, name: &str) {
    writeln!(out, "\n[{name}]").expect("a string never fails");
}

/// `key = value`, where the value is already TOML.
fn scalar(out: &mut String, key: &str, value: &str) {
    writeln!(out, "{key} = {value}").expect("a string never fails");
}

/// `key = []` where there is nothing, and otherwise **one element a line**, four spaces in, each
/// with a trailing comma so that adding a row by hand is a one-line edit.
fn one_a_line(out: &mut String, key: &str, elements: impl Iterator<Item = String>) {
    one_a_line_with_notes(out, key, elements.map(|element| (element, None)));
}

/// The same, where a row may carry **a comment saying where its value came from** — which is what
/// spec §4 chose this format for, and what a `defaulted` number owes its reader (§8).
///
/// **The note is per row and only on the rows that need one**, so its cost does not scale with the
/// cohort: a run of 3,000 samples whose coefficients were all fitted writes 3,000 rows and no
/// notes.
fn one_a_line_with_notes(
    out: &mut String,
    key: &str,
    elements: impl Iterator<Item = (String, Option<&'static str>)>,
) {
    let mut elements = elements.peekable();
    if elements.peek().is_none() {
        writeln!(out, "{key} = []").expect("a string never fails");
        return;
    }
    writeln!(out, "{key} = [").expect("a string never fails");
    for (element, note) in elements {
        // **Above the row rather than after it.** A trailing comment lengthens a line that is
        // already the longest thing in the file, and TOML lets a comment sit inside a multi-line
        // array, so the note goes where a reader meets it before the row it is about.
        if let Some(note) = note {
            for line in wrapped(note, ROOM_BESIDE_A_ROW) {
                writeln!(out, "    # {line}").expect("a string never fails");
            }
        }
        writeln!(out, "    {element},").expect("a string never fails");
    }
    writeln!(out, "]").expect("a string never fails");
}

/// How wide a comment line may be, all in: its `# ` and any indent are part of it.
const COMMENT_WIDTH: usize = 80;

/// What a note at the left margin has room for, after its `# `.
const ROOM_AT_THE_MARGIN: usize = COMMENT_WIDTH - 2;

/// **What an empty `slippage_by_stratum_and_group` means**, written above it because the section's
/// own paragraph describes a *missing row* and a reader takes the empty table for the same thing.
///
/// **Every number is read off [`StutterModel::hipstr_shipped`] rather than typed here**, so the
/// sentence a user is shown cannot come to disagree with the model their tracts were actually
/// scored under.
///
/// **And it is written in the three words the table above uses**, which is what a geneticist
/// reading a produced file asked for: the section teaches
/// `share_of_reads_that_slip` / `shorter_share` / `fall_off` forty lines earlier, and a note
/// quoting four direction shares instead cannot be lined up either against that table or against
/// a later fitted run's rows — which is the comparison the reader wants.
///
/// **The inversion, and the one place the shipped model is not a point the fit could produce.**
/// [`stutter_rates_for`](crate::ng::calling::likelihood::stutter_rates::stutter_rates_for) splits
/// `share_of_reads_that_slip` — the *whole-repeat* mass — by `shorter_share`, and adds a
/// part-repeat mass of
/// [`PART_REPEAT_SHARE_OF_WHOLE`](crate::ng::calling::likelihood::stutter_rates::PART_REPEAT_SHARE_OF_WHOLE)
/// times it on top; the one-step share is the *complement* of `fall_off`. So the shipped model's
/// whole-repeat pair inverts cleanly, and its part-repeat pair does not: at a slip share of 0.10 a
/// fitted row would carry 0.005 of part-repeat mass and the shipped model carries **0.02**, four
/// times as much. That is stated rather than smoothed over, because a reader comparing the
/// defaults against their own fitted rows would otherwise find a term that does not add up.
///
/// **What the note has to say and an earlier draft did not**, all from that reader: that **one**
/// pair of numbers stands in for every stratum, so a long mononucleotide run and a short
/// tetranucleotide are scored alike where real slippage differs several-fold across that range;
/// what a *part repeat* is, a word that appeared once in the whole file and was defined nowhere;
/// and that a direction share covers slips of any size, so *a whole repeat short* means one repeat
/// **or more**.
fn no_stratum_was_fitted() -> [String; 2] {
    let shipped = StutterModel::hipstr_shipped();
    // The fit's three numbers, recovered from the model's seven — see this function's doc for the
    // arithmetic and for the one term that does not invert.
    let slips = shipped.whole_repeat_shorter_share() + shipped.whole_repeat_longer_share();
    let shorter_share = shipped.whole_repeat_shorter_share() / slips;
    let fall_off = 1.0 - shipped.whole_repeat_one_step_share();
    let part_repeat = shipped.part_repeat_shorter_share() + shipped.part_repeat_longer_share();
    // Whole percents: these are the shares a reader argues with, and none of the shipped four is a
    // fraction of a percent (`every_share_the_note_quotes_is_a_whole_percent`).
    let in_a_hundred = |share: f64| (share * 100.0).round() as u32;
    // **Two decimals, and the file's own `{:?}` would be wrong here.** `fall_off` is derived as
    // `1 - 0.95` and comes out 0.050000000000000044; the writer prints values with `Debug` so that
    // a *value* round-trips, and this is prose re-expressing a model rather than a value to copy.
    let rounded = |share: f64| format!("{share:.2}");
    [
        format!(
            "This table is empty, which is not the same as a missing row: **no stratum was \
             fitted at all**, so every repeat tract in this run was scored under the stutter \
             model this caller ships. In this section's own three numbers that is \
             `share_of_reads_that_slip` = {slips}, `shorter_share` = {shorter_share}, `fall_off` \
             = {fall_off} — so {slips_in_a_hundred} reads in 100 misreport the tract length by a \
             whole number of repeats, half of them short and half long, and \
             {fall_off_in_a_hundred} in 100 of those misreports are by more than one repeat. A \
             further {part_repeat_in_a_hundred} in 100 change the tract by a part of a repeat: an \
             insertion or deletion inside it that is not a whole number of units, so a sequencing \
             indel or an interruption rather than slippage.",
            slips = rounded(slips),
            shorter_share = rounded(shorter_share),
            fall_off = rounded(fall_off),
            slips_in_a_hundred = in_a_hundred(slips),
            fall_off_in_a_hundred = in_a_hundred(fall_off),
            part_repeat_in_a_hundred = in_a_hundred(part_repeat),
        ),
        format!(
            "**One pair of numbers stands in for every stratum.** A 20-base mononucleotide run \
             and a 5-copy tetranucleotide are scored identically here, where real slippage rises \
             steeply as the period falls and as the tract lengthens — short-period long tracts \
             are where this is furthest wrong, and are where to distrust a call first. And a \
             `shorter_share` of {shorter_share} is dead even, where real slippage usually favours \
             the shorter tract; these are HipSTR's shipped starting values, which HipSTR itself \
             replaces by fitting, and they were fitted on no organism in particular. A PCR \
             preparation generally slips more than a PCR-free one, by an amount that depends on \
             how many cycles it ran. There is nothing here to change, only rows to add: the \
             paragraphs above say what a `slippage_by_stratum_and_group` row holds, and one \
             written by hand is read back like any other. Fitting the run is the better answer, \
             and until then the calls at repeat tracts rest on somebody else's chemistry.",
            shorter_share = rounded(shorter_share),
        ),
    ]
}

/// **What an empty `substitution_rate_by_stratum` means.** Same shape as
/// [`no_stratum_was_fitted`], smaller stakes, the same reader mistake, and the number likewise
/// read off the constant a cell actually takes.
fn no_substitution_rate_was_fitted() -> String {
    format!(
        "This table is empty, which is not the same as a missing row: nothing was fitted for any \
         read group at any stratum, so every tract's cells took the caller's stated \
         {DEFAULT_SSR_SUBSTITUTION_RATE} — about one base in a thousand read wrong inside a \
         tract. Base quality inside tracts is usually worse than outside them, so on real reads \
         that number is likely optimistic."
    )
}

/// What a note indented with the row it is about has room for, after four spaces and its `# `.
const ROOM_BESIDE_A_ROW: usize = COMMENT_WIDTH - 6;

/// **A note, broken at word boundaries so that no comment line runs past [`COMMENT_WIDTH`]** —
/// counting the `# ` every line carries and the indent a row's note sits at.
///
/// The rows this file writes are long and there is nothing to be done about that; a comment is
/// prose, and prose that needs scrolling is prose nobody reads.
///
/// `room` is what is left for the words after the prefix, which is why the two callers pass
/// different numbers: a section's note starts at the left margin and a row's is indented with the
/// row it is about.
fn wrapped(note: &str, room: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut line = String::new();
    for word in note.split_whitespace() {
        if !line.is_empty() && line.chars().count() + 1 + word.chars().count() > room {
            lines.push(std::mem::take(&mut line));
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        lines.push(line);
    }
    lines
}

/// **A comment, one line each**, written where it stands rather than beside a value.
///
/// TOML's comment runs to the end of its line, so a note that has to say more than a few words
/// goes above the thing it is about — which is where a section's framing belongs anyway, since it
/// is about every row rather than one.
fn note(out: &mut String, paragraphs: &[&str]) {
    for paragraph in paragraphs {
        if paragraph.is_empty() {
            writeln!(out, "#").expect("a string never fails");
        } else if paragraph.starts_with("  ") {
            // **A line that is laid out on purpose is emitted as it stands** — the bullet lists
            // where spec §5's states are set against each other are a table in all but name, and
            // reflowing one would run the three claims together.
            writeln!(out, "# {paragraph}").expect("a string never fails");
        } else {
            note_lines(out, &wrapped(paragraph, ROOM_AT_THE_MARGIN));
        }
    }
}

/// **A note whose words were computed rather than typed** — wrapped and laid out exactly like a
/// written one, so a reader cannot tell which is which and neither can go stale against the
/// other.
fn a_derived_note(out: &mut String, paragraphs: &[String]) {
    let paragraphs: Vec<&str> = paragraphs.iter().map(String::as_str).collect();
    note(out, &paragraphs);
}

/// **How much of this file rests on a measurement of the run's own reads**, written above
/// `format_version` so that it is the first thing a reader meets.
///
/// **The question a geneticist asks first of a file they did not watch being produced**, and
/// until step F1 the file did not answer it anywhere: spec §7 makes the three sources produce
/// the same shape on purpose, and the cost is that a run that guessed every number and a run
/// that fitted every number open identically. The answer *is* in the file — one `warrant` a
/// number — but as a hundred-odd separate answers a reader has to gather, and one who does not
/// know to gather them reads a defaults run's file as a fit's.
///
/// **Derived from the numbers, never recorded beside them** ([`ParametersFile::what_the_run_fitted`]),
/// which is the rule step E3 set for the missing-slippage note and holds here for a sharper
/// reason: this file invites its reader to edit a value and its warrant (spec §1.2 goal 3), and a
/// recorded count would then be a sentence at the top contradicting the numbers below it.
fn how_much_was_fitted(file: &ParametersFile) -> Vec<String> {
    let what = file.what_the_run_fitted();
    let groups = what.groups();
    let mut note = if what.nothing_was_fitted() {
        vec![
            format!(
                "**Nothing in this file was fitted from reads** — 0 of its {groups} groups of numbers. Every number here is a constant compiled into this caller or a value somebody handed it, and your reads were all scored under them."
            ),
            String::new(),
            format!(
                "The {groups} groups: {}.",
                a_list_of_groups(&GroupOfNumbers::EVERY)
            ),
        ]
    } else if what.not_fitted().is_empty() {
        vec![format!(
            "**All {groups} groups of numbers in this file were fitted from reads**: {}.",
            a_list_of_groups(what.fitted())
        )]
    } else {
        vec![format!(
            "**{} of the {groups} groups of numbers in this file were fitted from reads**, and {} were not: {}. What a group that was not fitted holds instead is said in its own section below.",
            what.fitted().len(),
            what.not_fitted().len(),
            a_list_of_groups(what.not_fitted()),
        )]
    };

    // **⚑ *Whose* reads is a second question and the file can only half answer it**, so the half
    // it cannot answer is stated rather than left for a reader to discover. Four of the seven
    // groups carry no `warrant`, so spec §2.1's demotion cannot reach them and a file fitted over
    // another cohort still shows them as fitted. Derived from
    // `GroupOfNumbers::states_whose_reads` so that a group that gains a warrant moves out of this
    // sentence by itself.
    let (with_a_warrant, without): (Vec<_>, Vec<_>) = GroupOfNumbers::EVERY
        .into_iter()
        .partition(|group: &GroupOfNumbers| group.states_whose_reads());
    note.push(String::new());
    note.push(format!(
        "**Whose reads is a second question, and only {} of the {groups} groups answer it.** {} carry a `warrant` on every number: `fitted_here` means it was estimated from that read group's or that sample's own reads; `borrowed` means there was too little of them, so the mean of the sample's other read groups was taken; `supplied` means the run was handed the number rather than fitting it; `defaulted` means nothing could be fitted and nothing was supplied, so a stated constant was used.",
        with_a_warrant.len(),
        a_capitalised_list_of_groups(&with_a_warrant),
    ));
    note.push(String::new());
    note.push(format!(
        "The other {} — {} — say how a number was arrived at and not whose reads it came from: a smoothing origin, a rung, or which reads a contamination fraction was fitted from. **So in a file whose numbers were fitted over a different cohort those {} still read as fitted**, and only the {} above are marked down to `supplied`.",
        without.len(),
        a_list_of_groups(&without),
        without.len(),
        with_a_warrant.len(),
    ));
    note.push(String::new());
    note.push(
        "A group counts as fitted where any part of it was; within one, each number's own `warrant` or origin says which."
            .to_owned(),
    );
    note
}

/// **What `[fitted_from]` is for**, which is not the same sentence on a run that fitted nothing.
///
/// The heading is the format's and does not move; a geneticist reading a defaults run's file took
/// it for a claim that the numbers under it had been fitted from the inputs it names. **What the
/// four bindings do is the same either way** — they are what stops a file being paired with the
/// wrong cohort — so only the first clause is conditional.
fn what_the_bindings_are_for(file: &ParametersFile) -> Vec<String> {
    let what_they_do = "A run whose reference, samples or read groups do not match these is refused; one whose census does not match keeps the numbers and reports every one of them as supplied rather than fitted.";
    vec![if file.what_the_run_fitted().nothing_was_fitted() {
        format!(
            "What this run ran against. **Nothing in this file was fitted from it** — this run fitted nothing — but these lines still bind the file. {what_they_do}"
        )
    } else {
        format!("What these numbers were fitted from. {what_they_do}")
    }]
}

/// **What an empty list of census terms means**, written above it because the section's own
/// paragraphs describe a census that exists and a reader takes the empty list for a fit whose
/// evidence went unrecorded.
///
/// Two runs reach it and both are ordinary (`run_streaming.md` §2): the run that fitted nothing,
/// and any direct-mode run, which reads its evidence from the alignment files and builds no psp
/// and no census.
fn this_run_had_no_census(file: &ParametersFile) -> Vec<String> {
    let mut note = vec![
        "**This run had no census** — no store of evidence was built, either because it fitted nothing or because it read its reads straight from the alignment files. So there is nothing here for another run to line its own census up against."
            .to_owned(),
        String::new(),
        // **⚑ What this used to say was false for the reader it named.** It promised that any
        // run reading this file would find a disagreement at the first line and demote it. Two
        // empty term lists agree — `census_disagreement` zips them and finds nothing on either
        // side — and a run that has no census of its own does not compare at all. Only a run
        // that *has* a census disagrees with this one, and it disagrees because this file names
        // none. Found by a geneticist reading a produced file, 2026-08-31.
        "**What another run does with that depends on whether it has a census of its own.** A run that has one finds a disagreement here — it names twelve terms and this file names none — and reports every number in this file as supplied rather than fitted. A run that has none, which is any run reading its reads straight from the alignment files, compares nothing and takes these numbers as they stand."
            .to_owned(),
    ];
    if file.what_the_run_fitted().nothing_was_fitted() {
        note.push(String::new());
        note.push(
            "On this file neither answer changes a number's warrant: no number in it says `fitted_here` to begin with."
                .to_owned(),
        );
    }
    // **A `#` on its own at the end**, so this derived paragraph does not run into the section's
    // written one below it — a reader otherwise meets two claims about the census as one block.
    note.push(String::new());
    note
}

/// The groups, in the file's own order, as a sentence's list — **capitalised**, for a list that
/// opens a sentence.
fn a_capitalised_list_of_groups(groups: &[GroupOfNumbers]) -> String {
    let list = a_list_of_groups(groups);
    let mut letters = list.chars();
    letters.next().map_or(list.clone(), |first| {
        first.to_uppercase().collect::<String>() + letters.as_str()
    })
}

/// The groups, in the file's own order, as a sentence's list — **with an *and* before the last**,
/// because these lists are read inside sentences and a bare comma before the final item reads as
/// another clause rather than as the end of the list.
fn a_list_of_groups(groups: &[GroupOfNumbers]) -> String {
    let named: Vec<&str> = groups
        .iter()
        .map(|group| group.in_the_readers_words())
        .collect();
    match named.split_last() {
        None => String::new(),
        Some((last, &[])) => (*last).to_owned(),
        Some((last, before)) => format!("{} and {last}", before.join(", ")),
    }
}

/// The same, for lines that had to be wrapped at run time.
fn note_lines(out: &mut String, lines: &[String]) {
    for line in lines {
        writeln!(out, "# {line}").expect("a string never fails");
    }
}

/// A `key = value` whose note, where it has one, goes on the lines above it.
fn scalar_with_note(out: &mut String, key: &str, value: &str, note: Option<&'static str>) {
    if let Some(note) = note {
        note_lines(out, &wrapped(note, ROOM_AT_THE_MARGIN));
    }
    writeln!(out, "{key} = {value}").expect("a string never fails");
}

/// **Where each of the file's defaultable numbers gets its default from** (spec §8).
///
/// One place, so that a number's origin and the sentence a reader is shown cannot drift apart, and
/// so that step E1 — which turns these into named constants with their origin recorded beside
/// them — has one list to reconcile against rather than five call sites.
///
/// **§8's fourth default has no entry here and cannot have one yet.** The per-(stratum × slippage
/// group) slippage numbers are to be defaulted from the GIAB HG002 alignments, and §8 asks for
/// "which alignments, at what depth, on which date" beside them — but a slippage number carries a
/// *smoothing origin* and no `Warrant`, so the file has no state in which one is `defaulted` and
/// nothing for the comment to attach to. The measurement does not exist either (§12 question 1).
/// Both halves are recorded at Checkpoint B.
mod origins {
    /// The base-quality multiplier, where no usable rate was fitted.
    ///
    /// **⚑ It does not say *taken at face value*, and it used to.** A geneticist reading a
    /// produced file met `value = 4.0, warrant = "defaulted"` under a note saying the reads were
    /// taken at the quality the instrument reported — six Phred apart from what the run actually
    /// scored them at. The owner's ruling of 2026-08-31 is that a read group nothing could be
    /// fitted for is **charged the pre-pass's stated rate** rather than believed, so the
    /// multiplier is that rate over the library's own mean reported error and is 1.0 only by
    /// coincidence. `a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts`
    /// pins 4.0 and 2.0 for two libraries reporting 2.5 × 10⁻⁴ and 5 × 10⁻⁴.
    /// **One line, because it is written once a read group.** The explanation of what such a
    /// multiplier *is* belongs to the section and is written there once; a run at the top of the
    /// committed cohort range (`CLAUDE.md`) can have thousands of read groups, and a defaults run
    /// gets this note on every one of them.
    pub const CALIBRATION_MULTIPLIER: &str = concat!(
        "not calibrated: no usable error rate could be fitted for this read group, so it is ",
        "charged the stated rate — see the note above this table"
    );

    /// The tract ladder's bottom rung — what a run falls back to where nothing was fitted.
    pub const FLAT_CONCENTRATION: &str = concat!(
        "this run fitted no stratum on its own tracts, so there was no median to take and ",
        "this is the caller's own constant"
    );

    /// The repeat-tract outlier weight, where a run inherited it rather than being handed one.
    pub const OUTLIER_WEIGHT: &str = concat!(
        "inherited from the existing caller and never measured here: the share of ",
        "repeat-tract reads that came from nowhere the model can explain"
    );

    /// A repeat-tract substitution rate that nothing could be fitted for.
    pub const SUBSTITUTION_RATE: &str =
        "nothing was fitted for this read group at this stratum, and nothing was supplied";

    /// An inbreeding coefficient nobody stated and nothing fitted.
    ///
    /// **Rewritten 2026-08-31, and the sentence it replaced was read by a geneticist as a bug
    /// report.** It said "inbreeding has no default: a run should not be able to write this
    /// line", which was true of the *fit* and became false of the *run* when the owner ruled that
    /// a coefficient nobody states is zero. A reader whose cohort selfs met it on the one row
    /// they most needed to change and stopped to report a defect instead of changing it.
    pub const INBREEDING_COEFFICIENT: &str = concat!(
        "nobody said how inbred this sample is, so it is scored as an outcrosser: the ",
        "genotype prior multiplies its heterozygote branch by 1 - F, and at zero that ",
        "branch is left alone. If your cohort selfs, say so — a landrace near F = 0.9 ",
        "scored at zero has every homozygous stretch of its genome treated as a ",
        "surprise. Set the value here and change the warrant beside it to \"supplied\", ",
        "or declare it when you run"
    );
}

/// The note a warranted number owes its reader, which is one only where it was defaulted.
fn where_it_came_from(value: &WarrantedValue, origin: &'static str) -> Option<&'static str> {
    (value.warrant == Warrant::Defaulted).then_some(origin)
}

/// One inline table: `{ a = 1, b = 2 }`, or `{}` for one with no keys.
fn an_inline_table(fields: &[(&str, String)]) -> String {
    if fields.is_empty() {
        return "{}".to_owned();
    }
    let mut out = String::from("{ ");
    for (position, (key, value)) in fields.iter().enumerate() {
        if position > 0 {
            out.push_str(", ");
        }
        write!(out, "{key} = {value}").expect("a string never fails");
    }
    out.push_str(" }");
    out
}

// ---------------------------------------------------------------------
// The two scalar kinds that need care
// ---------------------------------------------------------------------

/// **A TOML basic string.** Backslash and quote are escaped, the five named control characters
/// take their short escapes, and anything else below a space takes `\uXXXX`.
///
/// **Everything at or above a space passes through as itself**, including non-ASCII: a TOML file
/// is UTF-8 and a sample name is whatever the sequencing centre typed. Escaping it would make
/// the one field a person is most likely to search for unfindable.
fn a_toml_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{8}' => out.push_str("\\b"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\u{c}' => out.push_str("\\f"),
            '\r' => out.push_str("\\r"),
            other if (other as u32) < 0x20 || other == '\u{7f}' => {
                write!(out, "\\u{:04X}", other as u32).expect("a string never fails");
            }
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

/// **A TOML float that reads back as the same `f64`.**
///
/// Rust's `Debug` for a float emits the shortest decimal that round-trips and always carries a
/// decimal point or an exponent, which is what TOML asks of a float: `1` is an *integer* in TOML's
/// grammar and `1.0` is a float.
///
/// **Not because `1` would fail to read back here** — measured, and it does not: the `toml`
/// crate's derived deserialiser accepts an integer where the shape asks for an `f64`, so this
/// crate's own round trip cannot see the difference. What it changes is the file's type for every
/// *other* reader: `toml::Value` types `1` as an integer, and so does every parser in every other
/// language. A file whose concentration is a float or an integer depending on who opened it is
/// not the artefact spec §1.2 goal 3 describes.
///
/// **The three special values are named rather than spelled.** TOML writes them lower-case,
/// where Rust's own formatting gives `NaN` and `inf`; a run should not produce one, and a file
/// that carries one should say so in a form a reader can parse rather than failing at it.
///
/// **Whether this round-trips every value is step C3's to establish**, on a table of adversarial
/// ones. This step owes only that what it writes is a TOML float.
fn a_toml_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return if value.is_sign_positive() {
            "inf"
        } else {
            "-inf"
        }
        .to_owned();
    }
    format!("{value:?}")
}

/// **A TOML integer, which is signed and 64 bits wide.**
///
/// The specification says "arbitrary 64-bit signed integers (from −2^63 to 2^63−1) should be
/// accepted and handled losslessly. If an integer cannot be represented losslessly, an error must
/// be thrown" — so a `u64` above `i64::MAX` has no TOML spelling at all, and writing its digits
/// gives a document whose meaning depends on who reads it. **Measured on `u64::MAX` in an
/// evidence count:** the shape's own derived reader accepts it and the whole file compares equal,
/// `toml::Value` refuses it with *u64 value was too large*, and Python's `tomllib` accepts it
/// because Python's integers are unbounded. Three readers, three answers.
///
/// **No real count reaches it; a subtraction that wrapped does.** The value saturates rather than
/// refusing, on two grounds: a panic would discard a finished run's output over a count that was
/// already wrong before this writer saw it, and `to_toml` is otherwise infallible.
///
/// **It saturates in every build, not only in release, so that the behaviour a run gets is the
/// behaviour a test can see.** A debug assertion here would make the one path this exists for
/// unreachable from the test suite, which is the defect this whole function was added to close.
/// **What a reader sees is 9,223,372,036,854,775,807**, which no evidence count can be — there
/// are no nine quintillion reads — so the saturation announces itself, and refusing a count
/// sitting at exactly `i64::MAX` is one of the things step C2's `validate` owns on the way back
/// in.
fn a_toml_integer(count: u64) -> String {
    i64::try_from(count).unwrap_or(i64::MAX).to_string()
}

// ---------------------------------------------------------------------
// One function a row, in the order the file writes them
// ---------------------------------------------------------------------

fn a_read_group_row(row: &ReadGroupRow) -> String {
    an_inline_table(&[
        ("read_group", row.read_group.to_string()),
        ("declared_id", a_toml_string(&row.declared_id)),
        ("library", a_toml_string(&row.library)),
        ("sample", a_toml_string(&row.sample)),
    ])
}

fn a_census_term(term: &CensusTerm) -> String {
    an_inline_table(&[
        ("term", a_toml_string(&term.term)),
        ("digest", a_toml_string(&term.digest)),
    ])
}

fn a_calibration_row(row: &BaseQualityCalibrationRow) -> String {
    an_inline_table(&[
        ("read_group", row.read_group.to_string()),
        (
            "error_probability_multiplier",
            a_warranted_value(&row.error_probability_multiplier),
        ),
    ])
}

fn a_contamination_row(row: &ContaminationRow) -> String {
    let mut fields = vec![
        ("read_group", row.read_group.to_string()),
        ("library", a_toml_string(&row.library)),
    ];
    // **A read group that identified nothing writes no `measurement` key at all**, where a
    // fraction of zero beside two zero counts would claim it was measured and found clean.
    if let Some(measurement) = &row.measurement {
        fields.push(("measurement", a_contamination_measurement(measurement)));
    }
    an_inline_table(&fields)
}

fn a_contamination_measurement(measurement: &ContaminationMeasurement) -> String {
    an_inline_table(&[
        (
            "share_of_reads_from_another_sample",
            a_toml_float(measurement.share_of_reads_from_another_sample),
        ),
        (
            "markers_with_reads",
            a_toml_integer(measurement.markers_with_reads),
        ),
        (
            "reads_on_markers",
            a_toml_integer(measurement.reads_on_markers),
        ),
        (
            "fitted_from_reads_of",
            a_toml_string(match measurement.fitted_from_reads_of {
                super::ContaminationFittedFrom::ThisReadGroupsOwnReads => {
                    "this_read_groups_own_reads"
                }
                super::ContaminationFittedFrom::EveryReadOfThisSample => {
                    "every_read_of_this_sample"
                }
            }),
        ),
    ])
}

fn a_read_group_batch_row(row: &ReadGroupBatchRow) -> String {
    an_inline_table(&[
        ("read_group", row.read_group.to_string()),
        ("batch", row.batch.to_string()),
    ])
}

fn a_sample_batch_row(row: &SampleBatchRow) -> String {
    an_inline_table(&[
        ("sample", a_toml_string(&row.sample)),
        ("batch", row.batch.to_string()),
    ])
}

fn an_inbreeding_row(row: &InbreedingRow) -> String {
    an_inline_table(&[
        ("sample", a_toml_string(&row.sample)),
        (
            "inbreeding_coefficient",
            a_warranted_value(&row.inbreeding_coefficient),
        ),
    ])
}

fn a_slippage_group_row(row: &SlippageGroupRow) -> String {
    an_inline_table(&[
        ("read_group", row.read_group.to_string()),
        ("slippage_group", row.slippage_group.to_string()),
    ])
}

fn a_slippage_row(row: &SlippageRow) -> String {
    let mut fields = vec![
        ("period", row.period.to_string()),
        ("reference_repeats", a_toml_integer(row.reference_repeats)),
        ("slippage_group", row.slippage_group.to_string()),
        (
            "share_of_reads_that_slip",
            a_toml_float(row.share_of_reads_that_slip),
        ),
        ("shorter_share", a_toml_float(row.shorter_share)),
        ("fall_off", a_toml_float(row.fall_off)),
        (
            "share_of_reads_that_slip_origin",
            a_level_origin(&row.share_of_reads_that_slip_origin),
        ),
    ];
    // **No shares origin is no key**, which is what the fit says about a pair whose shares were
    // never recorded — not a shares origin whose fields are empty.
    if let Some(shares) = &row.shorter_share_and_fall_off_origin {
        fields.push(("shorter_share_and_fall_off_origin", a_shares_origin(shares)));
    }
    an_inline_table(&fields)
}

fn a_level_origin(origin: &LevelOrigin) -> String {
    let mut fields = vec![("smoothing", a_level_smoothing(&origin.smoothing))];
    if let Some(expected_slipped_reads) = origin.expected_slipped_reads {
        fields.push((
            "expected_slipped_reads",
            a_toml_float(expected_slipped_reads),
        ));
    }
    an_inline_table(&fields)
}

fn a_shares_origin(origin: &SharesOrigin) -> String {
    let mut fields = Vec::new();
    if let Some(expected_slipped_reads) = origin.expected_slipped_reads {
        fields.push((
            "expected_slipped_reads",
            a_toml_float(expected_slipped_reads),
        ));
    }
    fields.push((
        "shorter_share_smoothing",
        a_share_smoothing(&origin.shorter_share_smoothing),
    ));
    fields.push((
        "fall_off_smoothing",
        a_share_smoothing(&origin.fall_off_smoothing),
    ));
    an_inline_table(&fields)
}

/// **A smoothing that used no curve is a bare word; one that used a curve is a one-key table** —
/// serde's own spelling of an enum, and the device that stops a reach being written without the
/// curve it is a reach into.
fn a_level_smoothing(smoothing: &LevelSmoothing) -> String {
    match smoothing {
        LevelSmoothing::ThisStratum => a_toml_string("this_stratum"),
        LevelSmoothing::ThisPeriodsCurve { curve, reach } => an_inline_table(&[(
            "this_periods_curve",
            an_inline_table(&[
                ("curve", a_slippage_curve(curve)),
                ("reach", a_toml_string(the_word_for_reach(*reach))),
            ]),
        )]),
        LevelSmoothing::Blend {
            curve_weight,
            curve,
            reach,
        } => an_inline_table(&[(
            "blend",
            an_inline_table(&[
                ("curve_weight", a_toml_float(*curve_weight)),
                ("curve", a_slippage_curve(curve)),
                ("reach", a_toml_string(the_word_for_reach(*reach))),
            ]),
        )]),
    }
}

fn a_share_smoothing(smoothing: &ShareSmoothing) -> String {
    match smoothing {
        ShareSmoothing::ThisStratum => a_toml_string("this_stratum"),
        ShareSmoothing::ThisPeriodsCurve { curve, reach } => an_inline_table(&[(
            "this_periods_curve",
            an_inline_table(&[
                ("curve", a_share_curve(curve)),
                ("reach", a_toml_string(the_word_for_reach(*reach))),
            ]),
        )]),
        ShareSmoothing::Blend {
            curve_weight,
            curve,
            reach,
        } => an_inline_table(&[(
            "blend",
            an_inline_table(&[
                ("curve_weight", a_toml_float(*curve_weight)),
                ("curve", a_share_curve(curve)),
                ("reach", a_toml_string(the_word_for_reach(*reach))),
            ]),
        )]),
    }
}

fn a_slippage_curve(curve: &SlippageCurve) -> String {
    an_inline_table(&[
        ("rise_shape", a_toml_float(curve.rise_shape)),
        ("intercept", a_toml_float(curve.intercept)),
        ("slope", a_toml_float(curve.slope)),
        (
            "fitted_from_repeats",
            a_toml_integer(curve.fitted_from_repeats),
        ),
        ("fitted_to_repeats", a_toml_integer(curve.fitted_to_repeats)),
        ("held_out_error", a_toml_float(curve.held_out_error)),
        ("cells", a_toml_integer(curve.cells)),
    ])
}

fn a_share_curve(curve: &ShareCurve) -> String {
    an_inline_table(&[
        (
            "shape",
            a_toml_string(the_word_for_share_shape(curve.shape)),
        ),
        ("intercept", a_toml_float(curve.intercept)),
        ("slope", a_toml_float(curve.slope)),
        ("bend", a_toml_float(curve.bend)),
        ("centre_repeats", a_toml_float(curve.centre_repeats)),
        (
            "fitted_from_repeats",
            a_toml_integer(curve.fitted_from_repeats),
        ),
        ("fitted_to_repeats", a_toml_integer(curve.fitted_to_repeats)),
        ("held_out_error", a_toml_float(curve.held_out_error)),
        ("strata", a_toml_integer(curve.strata)),
        (
            "curve_fitted_on",
            a_toml_string(the_word_for_share_curve_rung(curve.curve_fitted_on)),
        ),
    ])
}

fn a_stratum_length_spectrum_row(row: &StratumLengthSpectrumRow) -> String {
    an_inline_table(&[
        ("period", row.period.to_string()),
        ("reference_repeats", a_toml_integer(row.reference_repeats)),
        ("concentration", a_toml_float(row.concentration)),
        (
            "shares_by_repeat_offset",
            a_float_array(&row.shares_by_repeat_offset),
        ),
    ])
}

fn a_period_length_spectrum_row(row: &PeriodLengthSpectrumRow) -> String {
    an_inline_table(&[
        ("period", row.period.to_string()),
        ("concentration", a_toml_float(row.concentration)),
        (
            "shares_by_repeat_offset",
            a_float_array(&row.shares_by_repeat_offset),
        ),
    ])
}

fn a_substitution_rate_row(row: &SubstitutionRateRow) -> String {
    an_inline_table(&[
        ("read_group", row.read_group.to_string()),
        ("period", row.period.to_string()),
        ("reference_repeats", a_toml_integer(row.reference_repeats)),
        ("ploidy", row.ploidy.to_string()),
        ("rate", a_warranted_value(&row.rate)),
    ])
}

/// A length spectrum's shares — **on one line**, because they are a handful of numbers describing
/// one thing rather than a table of rows.
fn a_float_array(values: &[f64]) -> String {
    let mut out = String::from("[");
    for (position, value) in values.iter().enumerate() {
        if position > 0 {
            out.push_str(", ");
        }
        out.push_str(&a_toml_float(*value));
    }
    out.push(']');
    out
}

/// **The one shape every four-state-warranted number in the file is written in** — and its count
/// is a key that is simply not there where no fit produced one.
fn a_warranted_value(value: &WarrantedValue) -> String {
    let mut fields = vec![
        ("value", a_toml_float(value.value)),
        (
            "warrant",
            a_toml_string(the_word_for_warrant(value.warrant)),
        ),
    ];
    if let Some(observations) = value.observations {
        fields.push(("observations", an_evidence_count(observations)));
    }
    an_inline_table(&fields)
}

/// **A count that names its own unit**, because the three units differ by orders of magnitude on
/// one cohort and a reader comparing two counts without knowing which is which is not comparing
/// anything.
fn an_evidence_count(count: EvidenceCount) -> String {
    let (unit, how_many) = match count {
        EvidenceCount::Reads(reads) => ("reads", reads),
        EvidenceCount::CoveredPositions(positions) => ("covered_positions", positions),
        EvidenceCount::BasesCompared(bases) => ("bases_compared", bases),
    };
    an_inline_table(&[(unit, a_toml_integer(how_many))])
}

// ---------------------------------------------------------------------
// The file's word for each unit-variant enum
// ---------------------------------------------------------------------
//
// **Written out here rather than derived from the `Serialize` impl**, so that this writer and
// serde's produce the same spelling only if both are right. Deriving it would make the golden
// file a tautology: a renamed variant would move the writer, the golden file and the reader
// together, which is the failure `every_enum_variant_spells_as_the_file_says` exists to prevent
// one level up. `the_hand_written_words_are_serdes_words` compares the two lists.

fn the_word_for_warrant(warrant: Warrant) -> &'static str {
    match warrant {
        Warrant::FittedHere => "fitted_here",
        Warrant::Borrowed => "borrowed",
        Warrant::Supplied => "supplied",
        Warrant::Defaulted => "defaulted",
    }
}

fn the_word_for_seed_rung(rung: SeedRung) -> &'static str {
    match rung {
        SeedRung::FittedCurve => "fitted_curve",
        SeedRung::NeutralShape => "neutral_shape",
        SeedRung::ZeroDiversity => "zero_diversity",
        SeedRung::StatedHeterozygosity => "stated_heterozygosity",
    }
}

fn the_word_for_reach(reach: CurveReach) -> &'static str {
    match reach {
        CurveReach::InsideTheFittedRange => "inside_the_fitted_range",
        CurveReach::BelowTheFittedRange => "below_the_fitted_range",
        CurveReach::AboveTheFittedRange => "above_the_fitted_range",
    }
}

fn the_word_for_share_curve_rung(rung: ShareCurveRung) -> &'static str {
    match rung {
        ShareCurveRung::ThisPeriod => "this_period",
        ShareCurveRung::ThisPeriodUnscored => "this_period_unscored",
        ShareCurveRung::OtherPeriods => "other_periods",
        ShareCurveRung::BuiltInDefault => "built_in_default",
    }
}

fn the_word_for_share_shape(shape: ShareShape) -> &'static str {
    match shape {
        ShareShape::Flat => "flat",
        ShareShape::Sloping => "sloping",
        ShareShape::Turning => "turning",
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::a_file_using_every_shape;
    use super::*;
    use serde::Serialize;

    /// **What the writer emits, read back through the shape's own `Deserialize`.**
    ///
    /// This is the half of spec §1.2's goal 1 that step B can hold: what the writer puts on disk
    /// is the value it was handed. The other half — a `RunParameters` written and read back to an
    /// equal `RunParameters` — needs step C's reader and is C4's.
    #[test]
    fn the_written_text_reads_back_as_the_same_file() {
        let written = a_file_using_every_shape();
        let text = written.to_toml();
        let read: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(read, written);
    }

    /// **Every key and every spelling of the file this writer produces, pinned against a
    /// checked-in copy.**
    ///
    /// **A second golden file rather than a shared one.** `testdata/every_shape.toml` is what
    /// `serde`'s derived serializer emits, and it is the artefact's key surface written in
    /// another crate's layout — array-of-table headers, table-valued fields last, the largest
    /// section unnamed. This one is the file a run writes. The two are built from the same
    /// fixture, so a key that differs between them is a defect in one of the two writers rather
    /// than a difference of fixture.
    #[test]
    fn the_whole_shape_writes_the_documented_toml() {
        assert_eq!(
            a_file_using_every_shape().to_toml(),
            include_str!("testdata/every_shape_as_written.toml"),
            "the written file no longer matches testdata/every_shape_as_written.toml; if the \
             change is intended, regenerate that file from this fixture"
        );
    }

    /// **Rewrite `testdata/every_shape_as_written.toml` from the fixture.**
    ///
    /// Ignored, so it never runs in an ordinary suite: it makes
    /// [`the_whole_shape_writes_the_documented_toml`] pass by definition, which is the one thing
    /// that test must not do on its own. Run it deliberately, after an intended change to the
    /// layout, and read the resulting diff:
    ///
    /// ```text
    /// cargo test --lib ng::calling::parameters_file -- --ignored regenerate_the_written
    /// ```
    #[test]
    #[ignore = "rewrites the golden file; run deliberately after an intended layout change"]
    fn regenerate_the_written_golden_file() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src/ng/calling/parameters_file/testdata/every_shape_as_written.toml");
        std::fs::write(&path, a_file_using_every_shape().to_toml())
            .expect("the golden file is writable");
    }

    /// **Every row of every table is one line, and every table's rows are its own.**
    ///
    /// The layout claim spec §4 asks for, checked structurally rather than by eye: between a
    /// `key = [` and its closing bracket there are exactly as many lines as the table has rows,
    /// each indented and each ending in a comma — so a person editing one sample's coefficient
    /// edits one line, and a row can be read without holding a column legend.
    #[test]
    fn every_row_of_every_table_is_one_line() {
        let file = a_file_using_every_shape();
        let text = file.to_toml();
        let lines: Vec<&str> = text.lines().collect();

        // **Each lookup names the section it is in**, because two sections hold a `by_sample`
        // and a search from the top of the file finds the first — which happens to have the same
        // row count, so doubling the inbreeding rows would have left this test green.
        for (section, key, rows) in [
            (
                "[fitted_from]",
                "read_groups",
                file.fitted_from.read_groups.len(),
            ),
            ("[fitted_from]", "samples", file.fitted_from.samples.len()),
            // **This one already carries a note**, on the read group whose multiplier is
            // defaulted — so it is the table that proves a note is a line and not a row.
            (
                "[base_quality_calibration]",
                "by_read_group",
                file.base_quality_calibration.by_read_group.len(),
            ),
            (
                "[sequencing_batches]",
                "by_sample",
                file.sequencing_batches.by_sample.len(),
            ),
            ("[inbreeding]", "by_sample", file.inbreeding.by_sample.len()),
            (
                "[repeat_tracts]",
                "slippage_by_stratum_and_group",
                file.repeat_tracts.slippage_by_stratum_and_group.len(),
            ),
            (
                "[repeat_tracts]",
                "substitution_rate_by_stratum",
                file.repeat_tracts.substitution_rate_by_stratum.len(),
            ),
        ] {
            let header = lines
                .iter()
                .position(|line| *line == section)
                .unwrap_or_else(|| panic!("the file opens `{section}`:\n{text}"));
            let opens = lines[header..]
                .iter()
                .position(|line| *line == format!("{key} = ["))
                .unwrap_or_else(|| panic!("`{section}` opens a `{key}` table:\n{text}"))
                + header;
            let closes = lines[opens..]
                .iter()
                .position(|line| *line == "]")
                .expect("and closes it")
                + opens;
            // **A note is a line and not a row.** A defaulted value writes an indented comment
            // above its row, so the count is over the rows rather than over the lines — and the
            // stated invariant would be false for any run that defaults one.
            let written: Vec<&&str> = lines[opens + 1..closes]
                .iter()
                .filter(|line| !line.trim_start().starts_with('#'))
                .collect();
            assert_eq!(
                written.len(),
                rows,
                "`{section}`'s `{key}` holds {rows} rows and its table writes {} of them",
                written.len()
            );
            for row in written {
                assert!(
                    row.starts_with("    ") && row.ends_with(','),
                    "every row of `{key}` is one indented line ending in a comma: {row}"
                );
            }
        }

        // **And the same file with a note inside two of those tables**, which is the shape the
        // invariant above was originally stated wrongly for: a note is a line and not a row.
        let mut with_notes = a_file_using_every_shape();
        with_notes.inbreeding.by_sample[0]
            .inbreeding_coefficient
            .warrant = Warrant::Defaulted;
        with_notes.repeat_tracts.substitution_rate_by_stratum[0]
            .rate
            .warrant = Warrant::Defaulted;
        let noted = with_notes.to_toml();
        let noted_lines: Vec<&str> = noted.lines().collect();
        // Anchored to its section: two sections hold a `by_sample`, and the first from the top
        // of the file is the batching's.
        let header = noted_lines
            .iter()
            .position(|line| *line == "[inbreeding]")
            .expect("the file opens `[inbreeding]`");
        let opens = noted_lines[header..]
            .iter()
            .position(|line| *line == "by_sample = [")
            .expect("the inbreeding table opens")
            + header;
        let closes = noted_lines[opens..]
            .iter()
            .position(|line| *line == "]")
            .expect("and closes")
            + opens;
        assert!(
            noted_lines[opens + 1..closes]
                .iter()
                .any(|line| line.trim_start().starts_with('#')),
            "the defaulted coefficient's note is inside the table:\n{noted}"
        );
        assert_eq!(
            noted_lines[opens + 1..closes]
                .iter()
                .filter(|line| !line.trim_start().starts_with('#'))
                .count(),
            with_notes.inbreeding.by_sample.len(),
            "and the rows are still one a line, with the note not counted as one:\n{noted}"
        );

        assert!(
            lines.contains(&"[repeat_tracts]"),
            "the largest section opens with its own header, which serde's writer omits because \
             every field of it is a table:\n{text}"
        );
    }

    /// **A whole number is written as a float, because in TOML `1` is an integer.**
    ///
    /// A multiplier of exactly one is a legitimate fitted answer (spec §3.3), and the flat
    /// concentration is exactly one; writing either as `1` gives a document that parses and then
    /// refuses to deserialise into the shape that wrote it.
    #[test]
    fn a_whole_number_is_written_as_a_float() {
        assert_eq!(a_toml_float(1.0), "1.0");
        assert_eq!(a_toml_float(0.0), "0.0");
        assert_eq!(a_toml_float(-0.0), "-0.0");
        assert_eq!(a_toml_float(23.0), "23.0");

        let one: f64 = toml::from_str::<toml::Value>(&format!("x = {}", a_toml_float(1.0)))
            .expect("parses")
            .get("x")
            .and_then(toml::Value::as_float)
            .expect("and is a float rather than an integer");
        assert_eq!(one, 1.0);
    }

    /// **The three values a float can be that are not numbers are written in TOML's own words.**
    ///
    /// A run should produce none of them. A file that carries one should say so in a form a
    /// reader can parse, rather than in Rust's spelling — `NaN` and `inf`, neither of which TOML
    /// knows.
    #[test]
    fn the_three_values_that_are_not_numbers_are_written_in_tomls_words() {
        assert_eq!(a_toml_float(f64::NAN), "nan");
        assert_eq!(a_toml_float(f64::INFINITY), "inf");
        assert_eq!(a_toml_float(f64::NEG_INFINITY), "-inf");

        let parsed: toml::Value = toml::from_str("x = nan\ny = inf\nz = -inf").expect("TOML's own");
        assert!(
            parsed
                .get("x")
                .and_then(toml::Value::as_float)
                .expect("a float")
                .is_nan()
        );
        assert_eq!(
            parsed.get("y").and_then(toml::Value::as_float),
            Some(f64::INFINITY)
        );
    }

    /// **A sample name is whatever the sequencing centre typed, and it comes back unchanged.**
    ///
    /// A quote and a backslash are escaped; everything at or above a space passes through as
    /// itself, including non-ASCII — escaping it would make the one field a person is most
    /// likely to search for unfindable.
    #[test]
    fn a_name_that_needs_escaping_is_escaped_and_reads_back() {
        assert_eq!(a_toml_string("plain"), "\"plain\"");
        assert_eq!(
            a_toml_string("a \"quoted\" name"),
            "\"a \\\"quoted\\\" name\""
        );
        assert_eq!(a_toml_string("back\\slash"), "\"back\\\\slash\"");
        assert_eq!(a_toml_string("a\tb\nc"), "\"a\\tb\\nc\"");
        assert_eq!(a_toml_string("\u{1}"), "\"\\u0001\"");
        assert_eq!(a_toml_string("Ailsa ‘Craig’"), "\"Ailsa ‘Craig’\"");

        for name in [
            "a \"quoted\" name",
            "back\\slash",
            "a\tb\nc",
            "\u{1}",
            "Ailsa ‘Craig’ \"×2\"",
        ] {
            let text = format!("x = {}", a_toml_string(name));
            let back: toml::Value =
                toml::from_str(&text).unwrap_or_else(|e| panic!("{e} in {text}"));
            assert_eq!(
                back.get("x").and_then(toml::Value::as_str),
                Some(name),
                "{name:?} did not survive the writer"
            );
        }
    }

    /// **The writer's word for each variant is serde's word**, and neither is derived from the
    /// other.
    ///
    /// This writer spells the enums itself, so that it and the derived serializer agree only if
    /// both are right: deriving the spelling would move the writer, the golden file and the
    /// reader together on a rename, which is the failure the shape's own spelling test exists to
    /// prevent one level up.
    #[test]
    fn the_hand_written_words_are_serdes_words() {
        fn serdes_word<T: Serialize>(value: T) -> String {
            toml::Value::try_from(value)
                .expect("a unit variant is a TOML value")
                .as_str()
                .expect("a unit variant spells as a bare string")
                .to_owned()
        }

        for warrant in [
            Warrant::FittedHere,
            Warrant::Borrowed,
            Warrant::Supplied,
            Warrant::Defaulted,
        ] {
            assert_eq!(
                the_word_for_warrant(warrant),
                serdes_word(warrant),
                "{warrant:?}"
            );
        }
        for rung in [
            SeedRung::FittedCurve,
            SeedRung::NeutralShape,
            SeedRung::ZeroDiversity,
            SeedRung::StatedHeterozygosity,
        ] {
            assert_eq!(the_word_for_seed_rung(rung), serdes_word(rung), "{rung:?}");
        }
        for reach in [
            CurveReach::InsideTheFittedRange,
            CurveReach::BelowTheFittedRange,
            CurveReach::AboveTheFittedRange,
        ] {
            assert_eq!(the_word_for_reach(reach), serdes_word(reach), "{reach:?}");
        }
        for rung in [
            ShareCurveRung::ThisPeriod,
            ShareCurveRung::ThisPeriodUnscored,
            ShareCurveRung::OtherPeriods,
            ShareCurveRung::BuiltInDefault,
        ] {
            assert_eq!(
                the_word_for_share_curve_rung(rung),
                serdes_word(rung),
                "{rung:?}"
            );
        }
        for shape in [ShareShape::Flat, ShareShape::Sloping, ShareShape::Turning] {
            assert_eq!(
                the_word_for_share_shape(shape),
                serdes_word(shape),
                "{shape:?}"
            );
        }
        for source in [
            super::super::ContaminationFittedFrom::ThisReadGroupsOwnReads,
            super::super::ContaminationFittedFrom::EveryReadOfThisSample,
        ] {
            let written = a_contamination_measurement(&ContaminationMeasurement {
                share_of_reads_from_another_sample: 0.5,
                markers_with_reads: 1,
                reads_on_markers: 1,
                fitted_from_reads_of: source,
            });
            assert!(
                written.contains(&format!("\"{}\"", serdes_word(source))),
                "{source:?} is written as {written}"
            );
        }
    }

    /// **A run that declared no batching writes `false`** — and that flag is the only thing in
    /// the file that says so, because the rows a defaulted batching writes are the rows a
    /// declaration of one batch holding every library would write.
    ///
    /// Neither fixture in this module carries the state, so hard-coding the flag `true` left the
    /// whole suite green until this test.
    #[test]
    fn a_run_that_declared_no_batching_writes_the_flag_as_false() {
        let mut file = a_file_using_every_shape();
        file.sequencing_batches.batching_was_declared = false;
        let text = file.to_toml();
        assert!(
            text.contains("batching_was_declared = false"),
            "an undeclared batching says so:\n{text}"
        );
        let read: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(read, file);
    }

    /// **A shares origin whose stratum fitted nothing writes no `expected_slipped_reads` key** — spec §5's
    /// "absence, never a sentinel", on the one `Option<f64>` both fixtures always fill.
    ///
    /// A stratum that borrowed has reads of its own and no level of its own to say how many of
    /// them slipped, so a zero there is a claim nothing measured. The derived writer has its own
    /// test for this key; replacing that writer left the state untested until this one.
    #[test]
    fn a_shares_origin_that_fitted_nothing_writes_no_slipped_reads_key() {
        let mut file = a_file_using_every_shape();
        file.repeat_tracts.slippage_by_stratum_and_group[0]
            .shorter_share_and_fall_off_origin
            .as_mut()
            .expect("the first row carries a shares origin")
            .expected_slipped_reads = None;
        let text = file.to_toml();
        let row = text
            .lines()
            .find(|line| {
                line.contains("shorter_share_and_fall_off_origin")
                    && !line.trim_start().starts_with('#')
            })
            .expect("the row is written");
        assert_eq!(
            row.matches("expected_slipped_reads").count(),
            1,
            "the level origin's count stays and the shares origin's is absent, not zero: {row}"
        );
        let read: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(read, file);
    }

    /// **A file whose every table is empty still parses** — the bottom of the committed range,
    /// where a cohort has one sample and no repeat tract was fitted at all.
    ///
    /// An empty table is written `key = []` and not omitted: none of the shape's `Vec` fields
    /// carries a serde default, so a missing key is a hard parse error rather than an empty list.
    /// Both empty paths — a table's and a length spectrum's float array — were reachable by no
    /// test until this one, and deleting either left the suite green.
    #[test]
    fn a_file_with_every_table_empty_writes_and_reads_back() {
        let mut file = a_file_using_every_shape();
        file.fitted_from.read_groups.clear();
        file.fitted_from.samples.clear();
        file.fitted_from.census.terms.clear();
        file.base_quality_calibration.by_read_group.clear();
        file.contamination = None;
        file.sequencing_batches.by_read_group.clear();
        file.sequencing_batches.by_sample.clear();
        file.inbreeding.by_sample.clear();
        file.repeat_tracts.slippage_group_by_read_group.clear();
        file.repeat_tracts.slippage_by_stratum_and_group.clear();
        file.repeat_tracts.length_spectrum_by_period.clear();
        file.repeat_tracts.substitution_rate_by_stratum.clear();
        for row in &mut file.repeat_tracts.length_spectrum_by_stratum {
            row.shares_by_repeat_offset.clear();
        }

        let text = file.to_toml();
        assert!(
            text.contains("read_groups = []") && text.contains("shares_by_repeat_offset = []"),
            "an empty table is written as an empty table:\n{text}"
        );
        let read: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(read, file);
    }

    /// **Every integer in the file is one TOML defines**, which is a signed 64-bit one.
    ///
    /// A count above `i64::MAX` cannot be spelled: the shape's own derived reader takes the digits
    /// and compares equal, while `toml::Value` refuses the same document — so a round trip through
    /// this crate alone cannot see the difference. The assertion that can is a parse through
    /// `toml::Value`, which is what every other reader models the file as.
    #[test]
    fn a_count_no_toml_integer_can_hold_still_writes_a_document_every_reader_accepts() {
        let mut file = a_file_using_every_shape();
        file.base_quality_calibration.by_read_group[0]
            .error_probability_multiplier
            .observations = Some(EvidenceCount::Reads(u64::MAX));
        let text = file.to_toml();

        toml::from_str::<toml::Value>(&text)
            .unwrap_or_else(|error| panic!("every reader accepts it: {error}\n\n{text}"));
        assert!(
            text.contains(&format!("reads = {}", i64::MAX)),
            "a count above TOML's largest integer saturates at it rather than writing digits no \
             file can carry:\n{text}"
        );
    }

    /// **Every defaulted number says where its default came from, and no fitted one does.**
    ///
    /// This is what spec §4 chose TOML for: "what a person needs beside a number is where it came
    /// from and what moving it costs". The comment is on the rows that carry a default and only
    /// those, so its cost does not scale with the cohort — a run of 3,000 samples whose
    /// coefficients were all fitted writes 3,000 rows and no notes.
    #[test]
    fn every_defaulted_number_says_where_its_default_came_from() {
        let text = a_file_using_every_shape().to_toml();
        let lines: Vec<&str> = text.lines().collect();

        // The calibration's defaulted row: the note sits on the lines above it.
        let defaulted = lines
            .iter()
            .position(|line| {
                line.contains("warrant = \"defaulted\"") && line.contains("read_group = 1")
            })
            .expect("read group 1's multiplier is defaulted");
        // **However many lines the note wrapped to.** A correct writer whose origin happens to
        // wrap to one line, or to three, must not fail this.
        let note_above: String = lines[..defaulted]
            .iter()
            .rev()
            .take_while(|line| line.trim_start().starts_with('#'))
            .fold(String::new(), |all, line| format!("{line} {all}"));
        assert!(
            note_above.contains("not calibrated"),
            "the note is on the lines above the row it is about, and they say: {note_above}"
        );

        // **Each note above its own key, not merely somewhere in the file.** The two defaulted
        // scalars sit six lines apart, and swapping the two origins at their call sites leaves
        // both texts present and both wrong.
        //
        // **The fallback concentration is defaulted in a file of its own here**, because it
        // cannot be in this fixture: that one fits a stratum, and the ladder's bottom rung then
        // states the median of the run's own concentrations rather than the compiled-in flat
        // one. A run that fitted nothing is the run that writes it defaulted.
        let mut nothing_to_fit = a_file_using_every_shape();
        nothing_to_fit
            .repeat_tracts
            .length_spectrum_by_stratum
            .clear();
        nothing_to_fit
            .repeat_tracts
            .fallback_length_spectrum_concentration = WarrantedValue {
            value: 1.0,
            warrant: Warrant::Defaulted,
            observations: None,
        };
        let flat = nothing_to_fit.to_toml();
        for (key, note, in_this_text) in [
            (
                "fallback_length_spectrum_concentration",
                "this run fitted no stratum on its own tracts",
                &flat,
            ),
            (
                "repeat_tract_outlier_weight",
                "inherited from the existing caller and never measured here",
                &text,
            ),
        ] {
            let lines: Vec<&str> = in_this_text.lines().collect();
            let at = lines
                .iter()
                .position(|line| line.starts_with(key))
                .unwrap_or_else(|| panic!("`{key}` is written:\n{in_this_text}"));
            let above: String = lines[..at]
                .iter()
                .rev()
                .take_while(|line| line.starts_with('#'))
                .fold(String::new(), |all, line| format!("{line} {all}"));
            assert!(
                above.contains(note),
                "`{key}`'s note is the one about it, and the note above it is: {above}"
            );
        }

        // **The two origins no fixture defaults**, reached by defaulting them here: five origins
        // are five interchangeable strings, and a text put beside the wrong quantity is invisible
        // wherever the quantity is never defaulted.
        let mut nothing_fitted = a_file_using_every_shape();
        nothing_fitted.repeat_tracts.substitution_rate_by_stratum[0]
            .rate
            .warrant = Warrant::Defaulted;
        nothing_fitted.inbreeding.by_sample[0]
            .inbreeding_coefficient
            .warrant = Warrant::Defaulted;
        let text = nothing_fitted.to_toml();
        assert!(
            text.contains("nothing was fitted for this read group at this stratum"),
            "a defaulted substitution rate says so:\n{text}"
        );
        assert!(
            text.contains("nobody said how inbred this sample is"),
            "a defaulted inbreeding coefficient says nobody stated one:\n{text}"
        );
        assert!(
            text.contains("scored as an outcrosser"),
            "and says what taking the default costs, which is what a reader of a selfing \
             cohort's file has to act on:\n{text}"
        );

        // And the two fitted rows above and below it carry no note of their own.
        let fitted = lines
            .iter()
            .position(|line| {
                line.contains("warrant = \"fitted_here\"") && line.contains("read_group = 0")
            })
            .expect("read group 0's multiplier is fitted");
        assert!(
            !lines[fitted - 1].trim_start().starts_with('#'),
            "a fitted number has no default to explain: {}",
            lines[fitted - 1]
        );
    }

    /// **A run whose every number was fitted writes no per-row notes at all** — so the comments
    /// cost a fixed number of lines rather than one a row.
    #[test]
    fn a_run_that_defaulted_nothing_writes_no_per_row_notes() {
        let mut file = a_file_using_every_shape();
        for row in &mut file.base_quality_calibration.by_read_group {
            row.error_probability_multiplier.warrant = Warrant::FittedHere;
        }
        for row in &mut file.repeat_tracts.substitution_rate_by_stratum {
            row.rate.warrant = Warrant::FittedHere;
        }
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .warrant = Warrant::FittedHere;
        file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Supplied;

        let text = file.to_toml();
        for note in [
            "no calibration:",
            "inherited from the existing",
            "this run fitted no stratum on its own tracts",
            "nothing was fitted for this read group",
            "inbreeding has no default",
        ] {
            assert!(
                !text.contains(note),
                "nothing was defaulted, so nothing explains a default — and this one is still \
                 here: {note}\n{text}"
            );
        }
        assert!(
            text.lines()
                .all(|line| !line.trim_start().starts_with("# ") || line.starts_with("# ")),
            "the notes that remain are the sections' own, at the left margin:\n{text}"
        );
    }

    /// **The comments do not change what the file means.**
    ///
    /// A comment runs to the end of its line in TOML, so a note that landed inside a row rather
    /// than above it would silently truncate the document. Stripping every comment must leave a
    /// file that reads back as the same value.
    #[test]
    fn the_comments_change_what_a_reader_learns_and_not_what_it_reads() {
        let file = a_file_using_every_shape();
        let text = file.to_toml();

        let commented: usize = text
            .lines()
            .filter(|line| line.trim_start().starts_with('#'))
            .count();
        assert!(
            commented > 20,
            "the file carries its explanations: {commented} comment lines"
        );

        let without_comments: String = text
            .lines()
            .filter(|line| !line.trim_start().starts_with('#'))
            .map(|line| format!("{line}\n"))
            .collect();
        let read: ParametersFile = toml::from_str(&without_comments)
            .unwrap_or_else(|error| panic!("{error}\n\n{without_comments}"));
        assert_eq!(read, file, "the comments carry no meaning the values need");

        let with_comments: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(with_comments, file);
    }

    /// **No comment line runs off the side of the page**, because a comment is prose and prose
    /// that needs scrolling is prose nobody reads.
    ///
    /// The rows themselves are long and this step cannot help that — the longest is a slippage
    /// row, and why it is long is the shape's business rather than the writer's.
    #[test]
    fn no_comment_line_is_longer_than_the_prose_it_carries() {
        let text = a_file_using_every_shape().to_toml();
        for line in text
            .lines()
            .filter(|line| line.trim_start().starts_with('#'))
        {
            assert!(
                line.chars().count() <= COMMENT_WIDTH,
                "a comment line is {} characters against a width of {COMMENT_WIDTH}: {line}",
                line.chars().count()
            );
        }
    }

    /// **A note is wrapped by characters and at the width it says**, which two mutations of the
    /// wrapper survived the rest of the suite because no note in the file happens to sit on the
    /// boundary or to carry a byte that is not a character.
    #[test]
    fn a_note_wraps_by_characters_at_the_width_it_states() {
        // Exactly 78 characters fits; one more does not. The words are sized so the boundary is
        // the only thing that decides.
        let exactly = format!("{} {}", "a".repeat(39), "b".repeat(38));
        assert_eq!(exactly.chars().count(), ROOM_AT_THE_MARGIN);
        assert_eq!(
            wrapped(&exactly, ROOM_AT_THE_MARGIN).len(),
            1,
            "{ROOM_AT_THE_MARGIN} characters is one line"
        );

        let one_over = format!("{} {}", "a".repeat(40), "b".repeat(38));
        assert_eq!(one_over.chars().count(), ROOM_AT_THE_MARGIN + 1);
        assert_eq!(
            wrapped(&one_over, ROOM_AT_THE_MARGIN).len(),
            2,
            "one more character is two lines"
        );

        // **A row's note has less room**, because it carries the row's indent as well as its `# `.
        assert_eq!(ROOM_BESIDE_A_ROW + 6, COMMENT_WIDTH);
        assert_eq!(
            wrapped(&exactly, ROOM_BESIDE_A_ROW).len(),
            2,
            "what fits at the margin does not fit beside a row"
        );

        // **Characters and not bytes.** A sample name or a plant's name in a note can be
        // non-ASCII, and counting its bytes would wrap a line that fits.
        let accented = format!("{} {}", "é".repeat(39), "b".repeat(38));
        assert_eq!(accented.chars().count(), ROOM_AT_THE_MARGIN);
        assert_eq!(
            wrapped(&accented, ROOM_AT_THE_MARGIN).len(),
            1,
            "78 accented characters are 117 bytes and still one line"
        );

        // A single word longer than the width is not cut in half; it goes on its own line.
        let one_long_word = "x".repeat(200);
        assert_eq!(
            wrapped(&one_long_word, ROOM_AT_THE_MARGIN),
            vec![one_long_word.clone()]
        );
        assert!(
            wrapped("", ROOM_AT_THE_MARGIN).is_empty(),
            "an empty note writes no lines"
        );
        assert!(
            wrapped("   ", ROOM_AT_THE_MARGIN).is_empty(),
            "and neither does one of only spaces"
        );
    }

    /// **Absence is a section that is not there, and a key that is not there** — spec §5's rule,
    /// on the way onto disk.
    #[test]
    fn every_absence_is_a_missing_section_or_a_missing_key() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;
        let text = file.to_toml();
        // **The section header, not the word.** Since step E3 an uncontaminated run writes a
        // *comment* in place of the section — the explanation of what the absence means used to
        // live inside the section it explains, so the reader holding such a file never saw it —
        // and that comment names the key it is about. What must not appear is the header itself,
        // which is a line of TOML rather than of prose.
        assert!(
            text.lines()
                .filter(|line| !line.trim_start().starts_with('#'))
                .all(|line| !line.contains("[contamination]")),
            "an uncontaminated run writes no contamination section at all:\n{text}"
        );
        // Joined back into sentences first: the writer wraps a note to the file's comment width,
        // so a phrase this long is split across `# ` lines and a reader reads the sentence.
        let prose: String = text
            .lines()
            .filter_map(|line| line.trim_start().strip_prefix('#'))
            .map(str::trim)
            .collect::<Vec<_>>()
            .join(" ");
        assert!(
            prose.contains("no `[contamination]` section, which means nobody identified any"),
            "and says so where the section would have been, because the note explaining the \
             absence used to be inside the section:\n{prose}"
        );
        assert!(
            text.lines()
                .filter(|line| !line.trim_start().starts_with('#'))
                .all(|line| !line.contains("measurement")),
            "and no row of one either:\n{text}"
        );

        let with_rows = a_file_using_every_shape().to_toml();
        let unmeasured = with_rows
            .lines()
            .find(|line| line.contains("library = \"lib4\"") && !line.trim_start().starts_with('#'))
            .expect("read group 1 has a row");
        assert!(
            !unmeasured.contains("measurement"),
            "a read group that identified nothing writes no measurement key: {unmeasured}"
        );

        // **The outlier weight rather than the fallback concentration**, which this fixture
        // cannot carry as `defaulted`: it fits a stratum, and the bottom rung then states the
        // median of the run's own concentrations. The one number in the file that is defaulted
        // in every run this build can assemble is this one.
        let defaulted = with_rows
            .lines()
            .find(|line| line.starts_with("repeat_tract_outlier_weight"))
            .expect("the stated constant is written");
        assert!(
            defaulted.contains("warrant = \"defaulted\"") && !defaulted.contains("observations"),
            "a defaulted number carries its warrant and no count: {defaulted}"
        );

        // **And a run that fitted no stratum at all does write a defaulted fallback**, with the
        // note saying where the constant came from — the state this fixture cannot show beside
        // its own fitted stratum.
        let mut nothing_fitted = a_file_using_every_shape();
        nothing_fitted
            .repeat_tracts
            .length_spectrum_by_stratum
            .clear();
        nothing_fitted
            .repeat_tracts
            .fallback_length_spectrum_concentration = WarrantedValue {
            value: 1.0,
            warrant: Warrant::Defaulted,
            observations: None,
        };
        let text = nothing_fitted.to_toml();
        assert!(
            text.contains("this run fitted no stratum on its own tracts"),
            "a defaulted fallback says where the constant came from:\n{text}"
        );
    }

    /// **Every share the empty-slippage note quotes is a whole percent**, which is the assumption
    /// `no_stratum_was_fitted` renders under: it prints `(share * 100).round()` and drops whatever
    /// is left, so a shipped share of, say, 0.035 would be shown to a user as "4 in 100".
    ///
    /// **Held rather than assumed, because the note is derived and the model is not this file's.**
    /// If `StutterModel::hipstr_shipped` ever moves to a value that is not a whole percent, this
    /// fails and the renderer needs a decimal — which is better than a file quietly rounding a
    /// number a reader is being invited to argue with.
    #[test]
    fn every_share_the_note_quotes_is_a_whole_percent() {
        let shipped = StutterModel::hipstr_shipped();
        for (what, share) in [
            ("whole repeat shorter", shipped.whole_repeat_shorter_share()),
            ("whole repeat longer", shipped.whole_repeat_longer_share()),
            ("part repeat shorter", shipped.part_repeat_shorter_share()),
            ("part repeat longer", shipped.part_repeat_longer_share()),
        ] {
            let in_a_hundred = share * 100.0;
            assert!(
                (in_a_hundred - in_a_hundred.round()).abs() < 1e-9,
                "the shipped {what} share is {share}, which is {in_a_hundred} in 100 and not a \
                 whole percent; the note rounds it and would show a reader a number the model \
                 does not hold"
            );
        }
    }
}

#[cfg(test)]
mod every_float_comes_back_bit_identical {
    //! **Step C3: whether a float this caller writes reads back as the same double.**
    //!
    //! Spec §4 says plainly that this "has not been checked here", and spec §1.2 goal 1 is what
    //! rests on it: the two-mode oracle compares VCFs, so a parameter that survives a write to
    //! five decimal places shows up there as a changed call, at some locus, with nothing to say
    //! why. **The fix if it failed would have been a serialiser that formats floats for
    //! round-trip, not a different file format** — so what is measured here is the formatting, in
    //! both writers, against the one reader.
    //!
    //! # What "the same" means here, and why it is not `==`
    //!
    //! Every assertion below compares `to_bits()`. **`==` cannot see two of the three things that
    //! can go wrong**: `-0.0 == 0.0` is true, so a lost sign passes; and `NaN != NaN`, so a value
    //! that became a `NaN` fails for the wrong reason. Only the bit pattern says the double that
    //! came back is the double that went out.
    //!
    //! # The two writers, and the one that is not this module's
    //!
    //! The artefact a run writes goes through [`a_toml_float`], which is this file's. The golden
    //! `testdata/every_shape.toml` goes through `serde`'s own serialiser, which is the `toml`
    //! crate's. **Both are checked**, because the module's round-trip tests read files written by
    //! each, and a defect in either would show up as a changed number rather than as a failure to
    //! parse.

    use super::super::ParametersFile;
    use super::super::tests::a_file_using_every_shape;
    use super::a_toml_float;

    /// **The doubles most likely to break a formatter, each with the reason it is here.**
    ///
    /// A table rather than a sweep, because a sweep says *nothing broke* and a table says *this
    /// did not break*. The sweep is below, and it is the wider net.
    /// **The two subnormals come from their bit patterns rather than from decimal literals**,
    /// which is why this is a function and not a `const`: a literal for the largest subnormal is
    /// a seventeen-digit number that has to be right, and `from_bits` says what it means.
    fn the_hard_ones() -> Vec<(f64, &'static str)> {
        let mut hard = vec![
            (
                f64::from_bits(1),
                "the smallest subnormal there is, one bit set",
            ),
            (
                f64::from_bits(0x000f_ffff_ffff_ffff),
                "the largest subnormal, one bit below the smallest normal",
            ),
            (
                -f64::from_bits(1),
                "and a negative subnormal — one class the sweep is a poor net for, since a \
                 subnormal is one draw in 2,048",
            ),
        ];
        hard.extend_from_slice(THE_HARD_LITERALS);
        hard
    }

    const THE_HARD_LITERALS: &[(f64, &str)] = &[
        (0.0, "zero"),
        (-0.0, "negative zero, whose sign no comparison sees"),
        (
            1.0,
            "one, which a formatter may write without its decimal point",
        ),
        (-1.0, "minus one"),
        (0.1, "a tenth, which is not a tenth"),
        (0.1 + 0.2, "0.1 + 0.2, the classic"),
        (
            1.0 / 3.0,
            "a third, which needs all seventeen significant digits",
        ),
        (
            f64::MIN_POSITIVE,
            "the smallest normal double, 2.2250738585072014e-308",
        ),
        (f64::MAX, "the largest finite double"),
        (-f64::MAX, "and the most negative"),
        (
            9_007_199_254_740_992.0,
            "2^53, above which not every whole number is a double",
        ),
        (
            1e16,
            "a whole number a formatter may write in exponent form",
        ),
        (1e-5, "and a small one"),
        (std::f64::consts::PI, "pi"),
        (
            1.000_000_000_000_000_2,
            "the smallest double above one, one unit in the last place away",
        ),
        // **At the magnitudes the fit works in**, so the table is not only torture: a slip share,
        // a Dirichlet concentration, a held-out error and a per-base substitution rate. **These
        // are the shape of a fitted number, not a fitted number** — the plan's third adversarial
        // category is C4's, and the module header says so.
        (
            0.042_1,
            "a slip share, at the size a stratum's fit gives one",
        ),
        (3.5, "a concentration, in chromosomes"),
        (
            0.333_333_333_333_333_3,
            "a held-out error that needs every digit",
        ),
        (0.001_2, "a per-base substitution rate inside a tract"),
    ];

    /// **Every hard value, written by this module's formatter and read back by the crate.**
    #[test]
    fn this_writers_floats_read_back_as_the_same_double() {
        for (value, why) in &the_hard_ones() {
            let written = a_toml_float(*value);
            let read: toml::Value = toml::from_str(&format!("value = {written}"))
                .unwrap_or_else(|error| panic!("{why}: {written} is not TOML: {error}"));
            let read = read["value"]
                .as_float()
                .unwrap_or_else(|| panic!("{why}: {written} did not read back as a float"));
            assert_eq!(
                read.to_bits(),
                value.to_bits(),
                "{why}: {value:e} was written as {written} and came back as {read:e}"
            );
        }
    }

    /// **And a sweep over pseudo-random bit patterns**, which is what says the table above is
    /// examples rather than the whole of it.
    ///
    /// **Ten thousand doubles drawn from the whole 64-bit space**, so most are enormous, tiny or
    /// subnormal — the region a decimal formatter is least often exercised over. Non-finite
    /// patterns are skipped: `validate` refuses those, and what this measures is the formatter
    /// rather than the range check.
    ///
    /// The generator is a fixed-seed xorshift rather than a crate: a failure has to be
    /// reproducible from this file alone, and the value is printed with its bits when it fails.
    #[test]
    fn ten_thousand_arbitrary_doubles_read_back_as_themselves() {
        let mut bits: u64 = 0x2545_f491_4f6c_dd1d;
        let mut checked = 0_u32;
        for _ in 0..40_000 {
            bits ^= bits << 13;
            bits ^= bits >> 7;
            bits ^= bits << 17;
            let value = f64::from_bits(bits);
            if !value.is_finite() {
                continue;
            }
            let written = a_toml_float(value);
            let read: toml::Value = toml::from_str(&format!("value = {written}"))
                .unwrap_or_else(|error| panic!("{written} (bits {bits:#x}) is not TOML: {error}"));
            assert_eq!(
                read["value"].as_float().map(f64::to_bits),
                Some(value.to_bits()),
                "{value:e} (bits {bits:#x}) was written as {written}"
            );
            checked += 1;
            if checked == 10_000 {
                return;
            }
        }
        panic!("the sweep ran out of draws before it had checked ten thousand finite doubles");
    }

    /// **The same values inside a whole file, through both writers and the one reader.**
    ///
    /// The test above formats one number in isolation; this one puts it where a run would — a
    /// bare scalar at the top of a section, a field of a one-line inline table, an entry of an
    /// array, and a field of a table nested three deep — and reads the whole file back.
    ///
    /// **It goes through `from_toml` and not `validate`**, deliberately: most of these values are
    /// not shares, concentrations or probabilities, and what is being measured is the text rather
    /// than the meaning.
    #[test]
    fn a_hard_float_survives_a_whole_file_through_either_writer() {
        for (value, why) in &the_hard_ones() {
            let mut file = a_file_using_every_shape();
            // A bare scalar under a section header.
            file.ordinary_site_prior.reference_concentration = *value;
            // A field of a one-line inline table.
            file.repeat_tracts.length_spectrum_by_stratum[0].concentration = *value;
            // An entry of an array of floats.
            file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset[1] = *value;
            // A field of a table three deep inside a row.
            file.repeat_tracts.slippage_by_stratum_and_group[1]
                .share_of_reads_that_slip_origin
                .expected_slipped_reads = Some(*value);

            for (writer, text) in [
                ("this module's writer", file.to_toml()),
                (
                    "serde's own serialiser",
                    toml::to_string(&file).expect("the shape serialises"),
                ),
            ] {
                let read = ParametersFile::from_toml(&text)
                    .unwrap_or_else(|error| panic!("{why}, {writer}: {error}\n{text}"));
                for (at, back) in [
                    (
                        "ordinary_site_prior.reference_concentration",
                        read.ordinary_site_prior.reference_concentration,
                    ),
                    (
                        "length_spectrum_by_stratum[0].concentration",
                        read.repeat_tracts.length_spectrum_by_stratum[0].concentration,
                    ),
                    (
                        "length_spectrum_by_stratum[0].shares_by_repeat_offset[1]",
                        read.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset[1],
                    ),
                    (
                        "slippage_by_stratum_and_group[1]…expected_slipped_reads",
                        read.repeat_tracts.slippage_by_stratum_and_group[1]
                            .share_of_reads_that_slip_origin
                            .expected_slipped_reads
                            .expect("it was written"),
                    ),
                ] {
                    assert_eq!(
                        back.to_bits(),
                        value.to_bits(),
                        "{why}, through {writer}, at {at}: {value:e} came back as {back:e}"
                    );
                }
            }
        }
    }

    /// **A whole file of hard floats round-trips as one object**, which is the claim spec §1.2
    /// goal 1 makes rather than the per-field one above.
    #[test]
    fn a_file_whose_every_float_is_hard_is_equal_to_itself_after_a_trip() {
        let mut file = a_file_using_every_shape();
        // Enough of the file's floats to reach each shape the writer has: the seed's two scalars,
        // a warranted value, a spectrum's entries, a curve's coefficients and a slippage number.
        file.ordinary_site_prior.reference_concentration = 1.0 / 3.0;
        file.ordinary_site_prior.alternative_concentration_total = 5e-324;
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .value = f64::MIN_POSITIVE;
        file.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset =
            vec![0.1 + 0.2, -0.0, 1.000_000_000_000_000_2];
        file.repeat_tracts.slippage_by_stratum_and_group[1].share_of_reads_that_slip = f64::MAX;
        file.stated_constants.repeat_tract_outlier_weight.value = 1e-5;

        for (writer, text) in [
            ("this module's writer", file.to_toml()),
            (
                "serde's own serialiser",
                toml::to_string(&file).expect("the shape serialises"),
            ),
        ] {
            let read = ParametersFile::from_toml(&text)
                .unwrap_or_else(|error| panic!("{writer}: {error}\n{text}"));
            assert_eq!(read, file, "through {writer}:\n{text}");
            // **And bit-identically**, which `PartialEq` on the shape cannot say: it compares
            // `f64`s with `==`, so a sign lost off a zero would pass the line above.
            assert_eq!(
                read.repeat_tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset[1]
                    .to_bits(),
                (-0.0_f64).to_bits(),
                "a negative zero keeps its sign through {writer}, which `==` cannot see"
            );
        }
    }
}

/// **The two things a file says about itself in prose, and both are derived** — how much of it
/// was fitted from the reader's data, and what `[fitted_from]` is for on a run that fitted
/// nothing.
///
/// # Why this module exists rather than a line in the golden file
///
/// The golden file pins the file a *fitted* run writes, and both notes take their other arm on a
/// run that fitted nothing — the arm no golden file in this module covers, because the defaults
/// run's file is not checked in. These tests read the text a defaults-shaped file produces and
/// assert what a geneticist would look for in it.
#[cfg(test)]
mod what_the_file_says_about_itself {
    use super::super::tests::{a_file_using_every_shape, unwrapped_comments};
    use super::super::{GroupOfNumbers, ParametersFile, SeedRung, Warrant};

    /// Everything the file says before its first key.
    fn the_opening_of(file: &ParametersFile) -> String {
        let text = file.to_toml();
        unwrapped_comments(&text[..text.find("format_version").expect("the file has one")])
    }

    /// The every-shape fixture, stripped back to a run that fitted nothing — every warrant
    /// weakened, every fitted table emptied, the census gone.
    fn a_file_a_defaults_run_would_write() -> ParametersFile {
        let mut file = a_file_using_every_shape();
        for row in &mut file.base_quality_calibration.by_read_group {
            row.error_probability_multiplier.warrant = Warrant::Defaulted;
            row.error_probability_multiplier.observations = None;
        }
        file.contamination = None;
        for row in &mut file.inbreeding.by_sample {
            row.inbreeding_coefficient.warrant = Warrant::Defaulted;
            row.inbreeding_coefficient.observations = None;
        }
        file.ordinary_site_prior.rung = SeedRung::StatedHeterozygosity;
        file.repeat_tracts.slippage_by_stratum_and_group.clear();
        file.repeat_tracts.length_spectrum_by_stratum.clear();
        file.repeat_tracts.length_spectrum_by_period.clear();
        file.repeat_tracts.substitution_rate_by_stratum.clear();
        file.fitted_from.census.terms.clear();
        file
    }

    /// **A run that fitted nothing says so before `format_version`**, which is where a reader
    /// meets it — and it names every group, so the reader does not have to know what the seven
    /// are.
    #[test]
    fn a_defaults_runs_file_opens_by_saying_nothing_was_fitted() {
        let opening = the_opening_of(&a_file_a_defaults_run_would_write());

        assert!(
            opening.contains("Nothing in this file was fitted from reads"),
            "{opening}"
        );
        assert!(opening.contains("0 of its 7 groups"), "{opening}");
        for group in GroupOfNumbers::EVERY {
            assert!(
                opening.contains(group.in_the_readers_words()),
                "the opening does not name {}:\n{opening}",
                group.key()
            );
        }
    }

    /// **A fitted run's file opens by saying so**, and does not claim more than it can: the
    /// headline counts groups, and the sentence after it tells the reader that a number's own
    /// warrant is what says whether *that* number was measured.
    #[test]
    fn a_fitted_runs_file_opens_by_saying_every_group_was_fitted() {
        let opening = the_opening_of(&a_file_using_every_shape());

        assert!(
            opening.contains("All 7 groups of numbers in this file were fitted from reads"),
            "{opening}"
        );
        // **The maximal claim is the one that most needs checking**, so this arm names the seven
        // too — a reader counting `[section]` headings gets nine and cannot reconstruct them.
        for group in GroupOfNumbers::EVERY {
            assert!(
                opening.contains(group.in_the_readers_words()),
                "the opening does not name {}:\n{opening}",
                group.key()
            );
        }
    }

    /// **The file says which of its groups can answer *whose* reads and which cannot**, on every
    /// file rather than only where it looks needed.
    ///
    /// Four of the seven carry no `warrant`, so spec §2.1's demotion cannot reach them
    /// ([`GroupOfNumbers::states_whose_reads`]). A reader who is not told that reads a demoted
    /// file's slippage as this run's own fit — which is what three reviewers found on 2026-08-31
    /// and what this sentence exists to stop.
    #[test]
    fn the_file_says_which_groups_can_say_whose_reads() {
        for file in [
            a_file_using_every_shape(),
            a_file_a_defaults_run_would_write(),
        ] {
            let opening = the_opening_of(&file);
            assert!(
                opening.contains("only 3 of the 7 groups answer it"),
                "{opening}"
            );
            assert!(
                opening.contains("still read as fitted"),
                "the file must say what a demotion cannot reach:\n{opening}"
            );
            // Every warrant word the file uses is defined where it is introduced — `borrowed`
            // above all, which the headline's count depends on and which appeared as a value
            // twice and was explained nowhere.
            for defined in [
                "`fitted_here` means it was estimated from",
                "`borrowed` means there was too little of them",
                "`supplied` means the run was handed the number",
                "`defaulted` means nothing could be fitted",
            ] {
                assert!(
                    opening.contains(defined),
                    "{defined} is not defined:\n{opening}"
                );
            }
        }
    }

    /// **A demoted file's opening line does not claim this run fitted it.**
    ///
    /// The state spec §2.1 creates: every warrant `supplied`, and four groups whose rows carry no
    /// warrant to demote. The headline counts them as fitted — it must, since the file says a fit
    /// produced them — and the sentence below it is what stops a reader taking that for *fitted
    /// from your data*.
    #[test]
    fn a_demoted_file_says_a_fit_produced_its_numbers_and_not_that_this_run_did() {
        let opening =
            the_opening_of(&a_file_using_every_shape().demoted_to_no_better_than_supplied());
        assert!(
            !opening.contains("from your data"),
            "the claim the file cannot support:\n{opening}"
        );
        assert!(
            opening.contains("were fitted from reads"),
            "and the claim it can:\n{opening}"
        );
        assert!(opening.contains("still read as fitted"), "{opening}");
    }

    /// **A partly fitted run is the commoner case and gets the count and the list.**
    ///
    /// Two groups dropped here, so the headline is a fraction rather than either extreme and the
    /// list names exactly the two.
    #[test]
    fn a_partly_fitted_runs_file_names_the_groups_it_did_not_fit() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;
        file.repeat_tracts.substitution_rate_by_stratum.clear();
        let opening = the_opening_of(&file);

        assert!(
            opening.contains("5 of the 7 groups of numbers in this file were fitted from reads"),
            "{opening}"
        );
        assert!(
            opening.contains("2 were not: contamination and repeat-tract substitution rates"),
            "{opening}"
        );
    }

    /// **`[fitted_from]` heads a file where nothing was fitted from anything**, and the section's
    /// own note is what stops it reading as a claim. The bindings still do what they do, so only
    /// the first clause moves.
    #[test]
    fn fitted_from_does_not_claim_a_fit_on_a_run_that_had_none() {
        let text = a_file_a_defaults_run_would_write().to_toml();
        let section = unwrapped_comments(
            &text[text.find("[fitted_from]").expect("the section")
                ..text.find("[fitted_from.census]").expect("and the next")],
        );

        assert!(
            section.contains("Nothing in this file was fitted from it"),
            "{section}"
        );
        assert!(
            section.contains("these lines still bind the file"),
            "{section}"
        );
        assert!(
            !section.contains("What these numbers were fitted from"),
            "the fitted run's opening clause is still there:\n{section}"
        );
    }

    /// And a fitted run's `[fitted_from]` keeps the sentence it always had.
    #[test]
    fn fitted_from_still_says_what_the_numbers_were_fitted_from() {
        let prose = unwrapped_comments(&a_file_using_every_shape().to_toml());
        assert!(
            prose.contains("What these numbers were fitted from."),
            "{prose}"
        );
    }

    /// **An empty census is explained where it stands.** A reader holding a defaults run's file
    /// found `terms = []` under a section whose paragraphs describe a census that exists.
    #[test]
    fn an_empty_census_says_the_run_had_none() {
        let prose = unwrapped_comments(&a_file_a_defaults_run_would_write().to_toml());
        assert!(prose.contains("**This run had no census**"), "{prose}");
        // **⚑ The note used to promise a demotion that does not happen.** Two empty term lists
        // agree, and a run with no census of its own compares nothing — so *a run that reads this
        // file will find a disagreement* was false for the very reader the sentence above it
        // names. Found by a geneticist reading a produced file, 2026-08-31.
        assert!(
            prose.contains("A run that has one finds a disagreement here"),
            "the demotion is the census-holding run's, and the note must say so:\n{prose}"
        );
        assert!(
            prose.contains("A run that has none, which is any run reading its reads straight"),
            "and the other run compares nothing:\n{prose}"
        );
        assert!(
            prose.contains("no number in it says `fitted_here` to begin with"),
            "on a file that fitted nothing neither answer changes a warrant, and the note says so"
        );
    }

    /// **A run with a census gets no such note**, so the explanation appears exactly where it is
    /// needed rather than in every file.
    #[test]
    fn a_file_with_a_census_carries_no_note_about_not_having_one() {
        let prose = unwrapped_comments(&a_file_using_every_shape().to_toml());
        assert!(!prose.contains("**This run had no census**"), "{prose}");
    }

    /// **Everything the notes claim still parses.** A derived paragraph is emitted as TOML
    /// comments, and a note that forgot its `#` on one line would be read as a key.
    #[test]
    fn a_defaults_runs_file_still_reads_back_as_itself() {
        let file = a_file_a_defaults_run_would_write();
        let read = ParametersFile::from_toml(&file.to_toml()).expect("its own text parses");
        assert_eq!(read, file);
    }
}
