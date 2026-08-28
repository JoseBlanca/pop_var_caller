//! **The run's parameters as a file** — the shape a `RunParameters` takes on disk, beside what
//! each number's warrant is and which inputs it was fitted from
//! (`doc/devel/ng/spec/parameters_file.md`).
//!
//! Today those numbers exist only in memory:
//! [`RunParameters`](crate::ng::calling::run_parameters::RunParameters) assembles them from the
//! parameter pre-pass and calling reads them through a borrowed view. Nothing writes them down
//! and nothing can read them in. This module is the type that does both — the file's **shape**,
//! which is deliberately not `RunParameters` itself:
//!
//! - **`RunParameters` is dense and indexed; the file is named.** The calibration and
//!   contamination axes are `Vec`s indexed by read-group id, and the inbreeding coefficients are
//!   a `Vec` in the run's sample order. A file carrying only an order would be silently wrong
//!   against a re-ordered sample list, so every row here names its sample or its read group as
//!   well as numbering it (spec §3.5).
//! - **`RunParameters` cannot say *absent*; the file must.** The five states of spec §5 are
//!   distinctions a reader that collapses them will get wrong, and every one of them is a
//!   **missing key** here rather than a value standing in for one:
//!
//!   | what is true | how the file says it |
//!   |---|---|
//!   | no read group identified any contamination | no `[contamination]` section |
//!   | a read group was measured and found clean | a row whose `measurement` has `fraction = 0` and non-zero evidence counts |
//!   | a read group's error rate could not be fitted | a multiplier of 1.0 whose `warrant` is `defaulted` |
//!   | a stratum was furnished from its period's curves | no row for it in `length_spectrum_by_stratum` |
//!   | a slippage group put no read in a stratum | no row for that pair in `slippage_by_stratum_and_group` |
//!
//!   The rule underneath all five: **`Option<T>` is absence and never a sentinel**, and a warrant
//!   is carried rather than inferred from the value.
//!
//! # Why it lives under `calling/`
//!
//! It sits beside [`run_parameters`](crate::ng::calling::run_parameters) because it is that
//! type's serialised form, and every projection this module will grow runs `RunParameters ↔
//! file`. The dependency direction settles it: `calling` already imports from
//! `parameter_estimation` in thirteen files, and `parameter_estimation` imports nothing from
//! `calling`. So reading the pre-pass's results from here adds no edge to the module graph,
//! where a top-level peer would acquire one into each. The in-house precedent is the census
//! file, which lives inside the stage that produces it
//! (`parameter_estimation/joint/census_file.rs`).
//!
//! # The TOML tree, and that these names are provisional
//!
//! **The key names and the shape of the tree are the coder's proposal, by the owner's decision
//! of 2026-08-28**, and are revised once there is a file a fitted run actually produced to argue
//! with, rather than on paper. What the file must *contain* is spec §3 and which states must
//! stay distinguishable is spec §5; neither is a naming choice. The proposal, and what was
//! rejected, is written down in
//! `doc/devel/reports/implementations/ng_parameters_file_a1_2026-08-28.md`.
//!
//! The tree, in the one-row-a-line inline-table form **the hand-written writer of step B2 is
//! meant to emit** — which is not what `serde`'s own serializer produces today; that writes
//! `[[array of tables]]` headers instead, and `tests/testdata/every_shape.toml` is what it
//! actually emits:
//!
//! ```toml
//! format_version = 1
//! ploidy = 2
//!
//! [fitted_from]                    # §3.1 — and what a mismatch refuses
//! reference_digest = "…"
//! samples = ["TS-1", "TS-2"]
//! read_groups = [ { read_group = 0, declared_id = "HWI.3", library = "lib3", sample = "TS-1" } ]
//! [fitted_from.census]
//! terms = [ { term = "the loci actually kept", digest = "…" } ]
//!
//! [base_quality_calibration]       # §3.3
//! by_read_group = [ { read_group = 0,
//!                     error_probability_multiplier = { value = 1.0324,
//!                                                      warrant = "fitted_here",
//!                                                      observations = { reads = 812344 } } } ]
//!
//! [contamination]                  # §3.4 — the whole section absent means uncontaminated,
//!                                  #        and a row with no `measurement` was not measured
//! by_read_group = [ { read_group = 0, library = "lib3",
//!                     measurement = { fraction = 0.031, markers_with_reads = 4211,
//!                                     reads_on_markers = 90233,
//!                                     fitted_from_reads_of = "this_read_groups_own_reads" } } ]
//!
//! [sequencing_batches]             # §3.4 — declared by the run, not fitted
//! was_declared_by_the_run = false
//! by_read_group = [ { read_group = 0, batch = 0 } ]
//! by_sample = [ { sample = "TS-1", batch = 0 } ]
//!
//! [inbreeding]                     # §3.5 — the file's only cohort-sized axis
//! by_sample = [ { sample = "TS-1", inbreeding_coefficient = { value = 0.42, warrant = "fitted_here", observations = { covered_positions = 180600412 } } } ]
//!
//! [ordinary_site_prior]            # §3.6 — the seed itself, never the moments behind it
//! reference_concentration = 1.0
//! alternative_concentration_total = 0.0006
//! rung = "fitted_curve"
//!
//! [repeat_tracts]                  # §3.7
//! stated_length_spectrum_concentration = { value = 1.0, warrant = "defaulted" }
//! slippage_group_by_read_group = [ { read_group = 0, slippage_group = 0 } ]
//! slippage_by_stratum_and_group = [ … one row a (stratum × slippage group) … ]
//! length_spectrum_by_stratum = [ … only where the stratum was fitted on its own tracts … ]
//! length_spectrum_by_period = [ … only where the run asked for the middle rung … ]
//! substitution_rate_by_stratum = [ … ]
//!
//! [stated_constants]               # §3.8 — the numbers no fit produces
//! repeat_tract_outlier_weight = { value = 0.01, warrant = "defaulted" }
//! ```
//!
//! # Three conventions the whole tree keeps
//!
//! **`read_group` is always the run's own dense index**, `0..n`, and never the `@RG ID` string —
//! that name is `declared_id`, and it appears exactly once, in the read-group table of
//! `[fitted_from]`. The dense index is what
//! [`ReadGroupParameters::calibration_of`](crate::ng::calling::likelihood::ReadGroupParameters)
//! indexes by, so it is the join every other section is written against.
//!
//! **`reference_repeats` is the repeat count a stratum is the bin for**, everywhere, and never
//! a candidate allele's repeat count. Two types in the tree spell that field differently —
//! [`census::Stratum`](crate::ng::parameter_estimation::joint::census::Stratum) calls it
//! `reference_repeats` and [`ssr::Stratum`](crate::ng::parameter_estimation::ssr::Stratum) calls
//! it `repeats` — and they are the same quantity, so the file uses one word **and one width**
//! for it.
//!
//! **Every table keyed by an axis is `by_<axis>`, and every row in it names its axis value.** A
//! positional array whose meaning is its index would be silently wrong against a re-ordered
//! sample list, which is the failure the previous convention exists to prevent; the sequencing
//! batching is written as named rows for that reason rather than as two bare integer arrays.
//!
//! # What this module is not, yet
//!
//! No reading and no writing: this step is the shape alone. Three things about the eventual
//! writer and reader are already visible from here and are recorded so the next step does not
//! rediscover them:
//!
//! - **The writer cannot be `serde`'s.** Spec §4 makes the format TOML *for the comments* — each
//!   defaulted value carries where its default came from — and no serde serializer can emit a
//!   comment. The `Serialize` derives here are for tests and for a round-trip cross-check, not
//!   for the artefact a run writes.
//! - **The reader can be**, in principle: `toml::from_str` reads an inline table into a struct
//!   whatever shape wrote it. `the_documented_inline_form_parses` is what holds that, because it
//!   is the claim steps B and C meet on and the derived serializer never produces the form.
//! - **A section whose fields are all tables gets no header line of its own from `serde`.**
//!   Every field of `[repeat_tracts]` and of `[stated_constants]` is now a table or an array, so
//!   the golden file has neither header line and the largest section of the file opens unnamed.
//!   That is
//!   `serde`'s rendering and not the artefact's design: the writer of step B2 emits the section
//!   headers itself. **Read `testdata/every_shape.toml` as a record of the key names, not as the
//!   file a run will produce.**
//! - **`serde` emits a struct's table-valued fields after its scalar ones**, whatever the
//!   declared order. So `format_version` and `ploidy` open the file because they are the only
//!   scalars at the top level, not because they are declared first; and in a `Blend` row
//!   `smoothing` is emitted *last*, after `slipped_reads`, because it became a table. A
//!   hand-written writer that means to match serde's bytes has to do the same — or, better,
//!   choose its own order and pin it with its own golden file.

use serde::{Deserialize, Serialize};

/// **Which version of this format this build writes.**
///
/// Carried from the start so that a reader has something to branch on if there is ever an older
/// file; what it *does* with one is deferred until there is one (spec §11).
pub const FORMAT_VERSION: u32 = 1;

/// What a number entitles a score to claim, as the file spells it.
///
/// **The file's own spelling of [`Provenance`](crate::ng::parameter_estimation::Provenance)**,
/// which already has these four states and already ranks `Supplied` below `Borrowed` — a number
/// the run was handed says nothing about *this* data, where a borrowed one is at least a
/// measurement of a neighbouring grain (spec §2). A separate enum rather than serde attributes
/// on that one, because the file's spelling is a compatibility surface and the in-memory enum's
/// variant names are not: renaming a variant should not silently re-interpret every file on
/// disk. **`every_enum_variant_spells_as_the_file_says` is what makes that guarantee real** —
/// the round-trip test cannot, because it moves both sides of a rename at once.
///
/// **This word is reserved for these four states.** Where a repeat tract's numbers record how
/// much of a period's curve was mixed into the stratum's own answer, that is *smoothing*, not a
/// warrant, and it is spelled so ([`LevelSmoothing`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Warrant {
    /// Estimated from this grain's own data.
    FittedHere,
    /// Too little data at this grain, so the mean of the sample's other read groups was taken.
    Borrowed,
    /// The run was handed this value rather than fitting it.
    Supplied,
    /// Nothing could be fitted and nothing was supplied, so a stated constant was used.
    Defaulted,
}

/// **A number, what entitles a score to claim it, and how much data stood behind it** — the
/// file's spelling of [`Estimate<T>`](crate::ng::parameter_estimation::Estimate), and the one
/// shape every four-state-warranted number in the file is written in.
///
/// **One shape rather than a value and a warrant side by side**, so that a reader who has
/// understood one warranted number has understood all of them. The flat alternative —
/// `error_probability_multiplier = 1.0324` beside `warrant = "fitted_here"` — reads better in a
/// row that holds exactly one such number and worse everywhere else, because it needs a
/// different suffixed key each time; and two spellings of one idea is what this shape exists to
/// remove.
///
/// **The evidence count names its own unit, in the file** ([`EvidenceCount`]), because the unit
/// differs by quantity and a reader cannot be sent to the source to find it. That is not a
/// hypothetical: this comment said "reads for a per-read rate, windows for an inbreeding
/// coefficient" and **both halves were wrong** — the inbreeding coefficient counts covered
/// positions and the repeat-tract substitution rate counts bases compared, neither of which is a
/// read, and a window is 100,000 bases. A count whose unit lives only here is a count nobody can
/// compare.
///
/// **It is absent where no fit produced one**, never zero: the tract ladder's stated
/// concentration is a median over strata rather than an estimate with a sample size, and the
/// repeat-tract outlier weight is a stated constant.
///
/// # Five numbers in this file say where they came from some other way
///
/// The four-state ladder is not the only thing that can stand behind a number, and forcing it
/// onto a quantity with a better answer would be inventing one. **None of these five carries a
/// [`Warrant`]**, and each has its own settled word:
///
/// - **a contamination fraction** stands on its two evidence counts, because *measured and found
///   clean* and *not measured* are both a fraction near zero and only the counts tell them apart
///   (spec §5);
/// - **the ordinary-site prior's two concentrations** carry [`OrdinarySitePrior::rung`], which
///   says which rung of the ladder produced them;
/// - **a slippage number** carries its *origin* — [`LevelOrigin`] and [`SharesOrigin`], how much
///   of its period's curve went into it — because a level fitted from 8,000 slipped reads and a
///   level interpolated off a curve through four cells are the same `f64`;
/// - **a length spectrum** is placed rather than annotated: *which table it is in* says whether
///   it is a stratum's own fit, its period's pool, or neither;
/// - **a read group's slippage group** is the run's own declaration and is not estimated at all.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WarrantedValue {
    /// The number itself.
    pub value: f64,
    /// What entitles a score resting on it to claim what it claims.
    pub warrant: Warrant,
    /// How much data stood behind it, **naming the unit it counts in** — absent where no fit
    /// produced a count, never zero.
    pub observations: Option<EvidenceCount>,
}

/// **How much data stood behind a number, in the unit that number counts in.**
///
/// **A variant per unit rather than a bare integer**, so the unit reaches the person reading the
/// file rather than living in a doc comment they will not open: `observations = { reads = 812344 }`
/// against `observations = { covered_positions = 1806 }`. The three units differ by orders of
/// magnitude on the same cohort, so a reader comparing two numbers' evidence without knowing
/// which is which is not comparing anything.
///
/// **A new quantity with a new unit adds a variant here**, which is a deliberate act, rather than
/// silently reusing a word that means something else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceCount {
    /// Reads — what a per-read rate is fitted over, as
    /// `parameter_estimation::generic::fallback` counts them.
    Reads(u64),
    /// Reference positions the sample's reads covered — what an inbreeding coefficient is fitted
    /// over (`parameter_estimation::generic::runs`, summed across windows of 100,000 bases).
    CoveredPositions(u64),
    /// Bases inside repeat tracts that were compared against the reference — what a repeat-tract
    /// substitution rate is fitted over (`parameter_estimation::ssr`).
    BasesCompared(u64),
}

/// **Every number calling runs on, as one file.**
///
/// **Every struct in this module refuses a key it does not know**
/// (`#[serde(deny_unknown_fields)]`), and that is a decision rather than a default. The file is
/// meant to be hand-edited (spec §1.2 goal 3), and serde's ordinary behaviour is to discard an
/// unrecognised key in silence — which on an `Option` field reads a one-letter typo back as
/// *absent*, and absence is data here (spec §5). A mistyped `curv` table would say *this
/// stratum's period had no curve*, which is a fitted fact, with no diagnostic. **The cost is
/// that a file written by a later build, carrying a key this one does not know, is refused
/// rather than partially read.** For a file where a dropped key changes what a number means,
/// refusing loudly is the safer of the two; [`FORMAT_VERSION`] is where a migration policy would
/// go if one is ever needed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParametersFile {
    /// Which version of this format the file is written in — [`FORMAT_VERSION`] for a file this
    /// build wrote.
    pub format_version: u32,
    /// How many copies of the genome the run called against. **A property of the run rather
    /// than of the fit**, written so that a supplied file cannot be paired with a run at a
    /// different ploidy without saying so (spec §3.2).
    pub ploidy: u8,
    /// What these numbers were fitted from — and what a run that does not match is refused for
    /// (spec §3.1, §6).
    pub fitted_from: InputsFittedFrom,
    /// How far to trust each library's own base qualities (spec §3.3).
    pub base_quality_calibration: BaseQualityCalibration,
    /// How much of each library's DNA came from somebody else (spec §3.4).
    ///
    /// **Absent means the run is uncontaminated** — no read group identified any fraction, so the
    /// read likelihood computes its plain formula, which is the *simple* case for that model and
    /// not the weak one. **A table of zeros is a different claim** and a reader that writes one
    /// for an uncontaminated run has said every library was measured and found clean, which
    /// nothing measured (spec §5, first row). This is the first of §5's five states and the only
    /// one expressed by a whole section going missing.
    pub contamination: Option<Contamination>,
    /// Who was sequenced beside whom — the population a contaminating read is drawn from
    /// (spec §3.4).
    pub sequencing_batches: SequencingBatches,
    /// How inbred each sample is (spec §3.5).
    pub inbreeding: Inbreeding,
    /// What the SNP/indel prior is seeded from (spec §3.6).
    pub ordinary_site_prior: OrdinarySitePrior,
    /// How often a read slips a repeat, and what lengths a tract's chromosomes are spread over
    /// (spec §3.7).
    pub repeat_tracts: RepeatTracts,
    /// The numbers no fit produces, written out rather than left in the binary (spec §3.8).
    pub stated_constants: StatedConstants,
}

// ---------------------------------------------------------------------
// §3.1 — what the numbers were fitted from
// ---------------------------------------------------------------------

/// **The inputs the fit ran on, so that a file cannot be silently paired with the wrong ones.**
///
/// Four bindings and the first three refuse (spec §6): a parameters file fitted against a
/// different assembly gives a plausible VCF with wrong repeat strata; one listing samples the
/// run does not have has its inbreeding coefficients against the wrong plants; one with a gap in
/// the read-group ids drops the highest read group and surfaces as a panic at whichever locus
/// first carries one of that library's reads. The fourth — the census — **demotes rather than
/// refuses**, because a file fitted from a different census of the same cohort is still
/// interpretable, merely less warranted.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InputsFittedFrom {
    /// A digest of the reference the fit ran against. **Refuses on a mismatch.**
    pub reference_digest: String,
    /// The run's samples, by name, in the run's own order — the order every per-sample axis is
    /// written and read in. **Refuses on a mismatch**, naming the samples.
    pub samples: Vec<String>,
    /// The run's read groups, one row each, in dense-index order. **Refuses on a gap.**
    pub read_groups: Vec<ReadGroupRow>,
    /// Which census produced these numbers. **Demotes rather than refuses.**
    pub census: CensusIdentity,
}

/// One read group of the run: its dense index, and the three names the alignment file gave it.
///
/// **Identity rather than parameters** (spec §12, question 3): nothing here is fitted, and it is
/// in the file because §6 needs something to check the dense read-group axis against. It sits in
/// `[fitted_from]` rather than beside the numbers so that nobody mistakes it for something to
/// edit.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReadGroupRow {
    /// The run's own dense index, `0..n` — what every other section joins on.
    pub read_group: u32,
    /// The `@RG ID` the alignment file declares. **The only place in the file this string
    /// appears**; everywhere else a read group is its index.
    pub declared_id: String,
    /// The `@RG LB`, or the name this project synthesized when the file declared none.
    pub library: String,
    /// The `@RG SM` — which plant these reads came from.
    pub sample: String,
}

/// **Which census these numbers were fitted from, named term by term.**
///
/// The census is a store of evidence built from a psp, and several different censuses can be
/// built from one psp — so a psp has no single census to name, and the terms it was *recorded*
/// under are what identify it
/// (`doc/devel/ng/spec/parameter_prepass_joint_records.md` §6.1).
///
/// **Digested rather than held, and kept term by term rather than as one digest.** The values
/// themselves include a per-stratum locus-count table, which would be the largest thing in the
/// file for a binding whose only use is an equality. One digest over all of them would answer
/// the equality just as well; keeping them apart costs twelve short lines and lets whatever
/// reports the demotion say *which* term differed, in the words
/// [`RecordingTerms::first_disagreement`](crate::ng::parameter_estimation::joint::census::RecordingTerms::first_disagreement)
/// already uses.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CensusIdentity {
    /// One entry a recording term, in the order that type compares them in.
    pub terms: Vec<CensusTerm>,
}

/// One recording term of the census, and a digest of what it held.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CensusTerm {
    /// The term's name, in the words the census's own disagreement report uses.
    pub term: String,
    /// A digest of the term's value, lower-case hex.
    pub digest: String,
}

// ---------------------------------------------------------------------
// §3.3 — per read group, the base-quality calibration
// ---------------------------------------------------------------------

/// **One multiplier on each read's own reported error probability, per read group.**
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BaseQualityCalibration {
    /// One row a read group, in dense-index order.
    pub by_read_group: Vec<BaseQualityCalibrationRow>,
}

/// One read group's base-quality calibration.
///
/// **A multiplier of exactly 1.0 is a legitimate fitted answer as well as the default**, which
/// is why the warrant travels with the value rather than being inferred from it (spec §3.3, §5).
/// A read group with no usable rate gets a multiplier of one marked `defaulted` — never a fitted
/// zero, which would charge every read of that library the error floor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BaseQualityCalibrationRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// What each read's own reported **error probability** is multiplied by — so a value above
    /// one says the instrument was optimistic and the reads are worse than it claimed. **One**
    /// leaves the qualities exactly as reported. It is not a multiplier on the Phred score,
    /// which moves the other way.
    ///
    /// **Its evidence count is in reads, and `RunParameters` is not where to get it.**
    /// [`ReadGroupCalibration`](crate::ng::calling::likelihood::ReadGroupCalibration) is a
    /// multiplier and a provenance and keeps no count, but spec §3.3 asks for one and the fit
    /// produced one: it is on the `Estimate<ErrorRate>` that
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)
    /// reads and does not store. **So step B1's projection has to be handed that estimate**, the
    /// same motion the inbreeding coefficient needs.
    pub error_probability_multiplier: WarrantedValue,
}

// ---------------------------------------------------------------------
// §3.4 — per read group, contamination and the batching it was drawn against
// ---------------------------------------------------------------------

/// **How much of each library's DNA came from somebody else.**
///
/// **The grain is the read group, and a row names both it and its library** (spec §3.4): a
/// preparation sequenced over several lanes gives several read groups sharing one library name,
/// and two lanes of one preparation can carry genuinely different fractions, because index
/// hopping happens on a flowcell and not in a tube.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Contamination {
    /// One row a read group, in dense-index order. **Where some read group identified a
    /// fraction, every read group needs a row** — one that identified nothing gets a row with no
    /// `measurement`, which is how the file says *this library was not measured* inside a run
    /// where others were. A run where nobody identified anything has no section at all rather
    /// than a list of unmeasured rows.
    ///
    /// **Two shapes here are not states the file means to have, and both say §5's first state a
    /// second way.** An empty list says "this run is contaminated, and here are none of its read
    /// groups"; a list in which *no* row has a measurement says "contaminated, and nobody
    /// measured anything", which is the uncontaminated run written out longhand. In memory the
    /// second is collapsed on purpose —
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)
    /// turns all-absent views into an empty list and `view()` then takes the uncontaminated
    /// path — so a file spelling it as unmeasured rows and read literally would take the mixture
    /// path instead, with every fraction zero, where absence takes the plain formula.
    ///
    /// **Both are refused by step C2**, which is where the table meets the run's dense read-group
    /// axis; no shape can say a `Vec` is non-empty, and none can say *not every row is absent*.
    /// (Not step C1: that step is parsing, and both of these parse.)
    pub by_read_group: Vec<ContaminationRow>,
}

/// One read group, inside a run where **some** read group identified a fraction — and what this
/// one was found to carry, if anything.
///
/// **The second of spec §5's five states lives here**, and it is the one a reader is most likely
/// to collapse: *measured and found clean* and *not measured* both come back as a fraction near
/// zero, and only the evidence tells them apart. In memory that distinction rides on the counts —
/// *either* of them being zero means unmeasured, which
/// [`ContaminationView::was_measured`](crate::ng::calling::likelihood::ContaminationView) has to
/// be asked about; here the
/// unmeasured row simply has no measurement, so there is no fraction to misread and no evidence
/// count to compare against zero.
///
/// **It also removes a wart the in-memory type documents and cannot fix.** A read group that
/// identified nothing still has to carry a `ContaminationSource` there, and neither variant is
/// true of it — `run_parameters.rs` says so in as many words. Here it has none, because it has no
/// measurement to have one on.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContaminationRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// The library this read group's reads were prepared from — several read groups of one
    /// preparation share it, and it is written here because the whole point of the read-group
    /// grain is that two lanes of one library may differ.
    pub library: String,
    /// What the fit found for this read group. **Absent where it identified nothing**, which is
    /// not the same as a fraction of zero.
    pub measurement: Option<ContaminationMeasurement>,
}

/// What a contamination fit found for one read group: how much, and on what evidence.
///
/// **A fraction of zero here is a real answer** — this read group was measured and found clean —
/// because the evidence counts beside it are what say it was measured at all.
///
/// **A measurement whose counts are both zero is writable and is refused by step C2**, not by the
/// shape. It is exactly the in-memory `UNMEASURED_READ_GROUP`
/// ([`run_parameters.rs`](crate::ng::calling::run_parameters)), so a projection written from the
/// *view* rather than from the fit's estimate would produce one: a row saying *measured* while
/// carrying the evidence of *not measured*, with a `fitted_from_reads_of` that is true of
/// nothing. **`NonZeroU64` would make it unwritable and was not taken**, because a fit that
/// returns an estimate with no evidence behind it would then have no file to be written to at
/// all, and whether that can happen is a question about the contamination estimator rather than
/// about this shape — step C4's round trip on a real fit is what will answer it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContaminationMeasurement {
    /// The share of this read group's reads that came from another individual.
    pub fraction: f64,
    /// How many of the panel's varying positions this read group put a read on.
    pub markers_with_reads: u64,
    /// How many reads it put there.
    pub reads_on_markers: u64,
    /// Whose reads the fraction was fitted from. A fraction fitted from this library's own reads
    /// and one fitted from every read of the plant are different claims — the first can say two
    /// libraries of one sample differ, the second cannot.
    pub fitted_from_reads_of: ContaminationFittedFrom,
}

/// Whose reads a contamination fraction was fitted from — the file's spelling of
/// [`ContaminationSource`](crate::ng::parameter_estimation::joint::contamination::ContaminationSource).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ContaminationFittedFrom {
    /// This read group's own reads.
    ThisReadGroupsOwnReads,
    /// Every read of the plant, the answer then copied onto this read group.
    EveryReadOfThisSample,
}

/// **Who was sequenced beside whom, as the run was told** — the population a contaminating read
/// is drawn from.
///
/// **Written even where no contamination was fitted**, because it is a fact about the run rather
/// than about the fit (spec §3.4).
///
/// **Named rows rather than two bare integer arrays**, though the arrays would be a third the
/// size. The two axes are the same type and adjacent, so a positional pair can be exchanged
/// without any file becoming malformed — and crossing them draws a contaminating read's genotype
/// from the wrong population, with no panic and no parse error. A person moving one plant into
/// another batch also has to be able to find that plant's row.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SequencingBatches {
    /// Whether the run declared this batching, or it defaulted to one batch holding everything.
    /// **Two runs under different batchings produce fractions that are not comparable**, and
    /// this is the only thing that can tell a declared batching from an assumed one.
    pub was_declared_by_the_run: bool,
    /// One row a read group.
    pub by_read_group: Vec<ReadGroupBatchRow>,
    /// One row a sample.
    pub by_sample: Vec<SampleBatchRow>,
}

/// Which sequencing batch one read group sits in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReadGroupBatchRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// The batch it was sequenced in.
    pub batch: u32,
}

/// Which sequencing batch one sample sits in.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleBatchRow {
    /// The sample's name, as the run's read-group table spells it.
    pub sample: String,
    /// The batch it was sequenced in.
    pub batch: u32,
}

// ---------------------------------------------------------------------
// §3.5 — per sample, the inbreeding coefficient
// ---------------------------------------------------------------------

/// **How inbred each sample is** — one row a sample, and the file's only cohort-sized axis.
///
/// At the top of the committed range, 3,000 samples, this is 3,000 lines and the file is a few
/// hundred kilobytes: negligible beside the VCF it sits next to, and still openable in an editor
/// (spec §9).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Inbreeding {
    /// One row a sample, in the run's sample order. **At least one is required** — a run has at
    /// least one sample, so an empty list is a run whose sample order went missing rather than a
    /// run with none, which is what
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)
    /// asserts in release. **Refused by step C2**, alongside the contamination table's two empty
    /// shapes and for the same reason: no shape can say a `Vec` is non-empty.
    pub by_sample: Vec<InbreedingRow>,
}

/// One sample's inbreeding coefficient.
///
/// **The name is written beside the value, not just the order**, because the order is the run's
/// and a file carrying only an order would be silently wrong against a re-ordered sample list
/// (spec §3.5).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InbreedingRow {
    /// The sample's name, as the run's read-group table spells it.
    pub sample: String,
    /// The coefficient itself, a fraction in `[0, 1)`. Its evidence count is in **covered
    /// reference positions** — `chain.covered_positions_total()`, summed across windows of
    /// 100,000 bases (`parameter_estimation/generic/runs.rs`,
    /// `parameter_estimation/generic/fallback.rs`). **Not a window count**, which the runs model
    /// keeps separately and which is smaller by up to five orders of magnitude.
    ///
    /// **The warrant exists upstream and the seam into calling drops it.** The pre-pass produces
    /// an `Estimate<InbreedingF>`; the seam is
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble),
    /// which takes a bare `Vec<InbreedingF>`. So **step B1's projection has to be handed the
    /// estimates rather than that vector**, or every sample's coefficient reads as supplied when
    /// it was fitted.
    pub inbreeding_coefficient: WarrantedValue,
}

// ---------------------------------------------------------------------
// §3.6 — the ordinary-site prior's seed
// ---------------------------------------------------------------------

/// **What the SNP/indel prior is seeded from: two concentrations and which rung of the ladder
/// they came from.**
///
/// **Written as the seed, not as the moments it came from** (spec §3.6). The seed is built once
/// per run from two integrals of the fitted frequency density, and what varies per locus is only
/// how it is spread across that locus's alleles. Writing the moments instead would mean the
/// reader re-deriving the seed, and any change to that derivation would silently re-interpret
/// every existing file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OrdinarySitePrior {
    /// The reference allele's concentration.
    pub reference_concentration: f64,
    /// The concentration shared out across whatever alternative alleles a locus carries.
    /// **Exactly zero is a real answer** — a fully invariant cohort — and is not floored on the
    /// way in or out; the flooring belongs to the per-locus expansion.
    pub alternative_concentration_total: f64,
    /// Which rung of the ladder the pair came from.
    pub rung: SeedRung,
}

/// Which rung produced the seed — the file's spelling of
/// [`SeedRegime`](crate::ng::calling::genotype_prior::SeedRegime).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeedRung {
    /// Both moments came off the run's own fitted population curve — the top of the ladder.
    FittedCurve,
    /// No population curve was fitted, so the pair is the neutral shape at the heterozygosity
    /// the pre-pass **did** fit.
    NeutralShape,
    /// No heterozygosity was fitted either, so the pair rests on a stated species-range
    /// heterozygosity taken from human data. **The rung that must never be silent.**
    StatedHeterozygosity,
}

// ---------------------------------------------------------------------
// §3.7 — repeat tracts
// ---------------------------------------------------------------------

/// **The largest section: how often a read slips a repeat, what lengths a tract's chromosomes
/// are spread over, and how often a base reads wrong inside a tract.**
///
/// Its axes are what matter (spec §3.7). The slippage numbers are keyed by **slippage group**
/// and not by read group, which is what stops the section growing with the number of libraries;
/// the length spectra are keyed by stratum and by motif period, the two fitted rungs of the
/// tract ladder.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RepeatTracts {
    /// The strength the tract ladder's **bottom** rung states — a flat shape over whatever
    /// lengths the locus offers, held with this many chromosomes' worth of belief.
    ///
    /// **Its warrant is what tells the run's own median from the stated constant**, and without
    /// it the two are the same line: the median over the strata a run fitted can land on exactly
    /// 1.0, which is also `STATED_FLAT_CONCENTRATION`, and spec §8 counts this among the three
    /// numbers with an honest default that must be marked `defaulted` when it is one. It is the
    /// same reason the calibration multiplier carries a warrant beside a value of one. Its
    /// observation count is absent: a median over strata is not an estimate with a sample size.
    pub stated_length_spectrum_concentration: WarrantedValue,
    /// Which set of slippage numbers each read group's reads are drawn under. **The run's own
    /// declaration, not something inferred.**
    pub slippage_group_by_read_group: Vec<SlippageGroupRow>,
    /// One row a `(stratum × slippage group)`. **A pair with no row is a slippage group that put
    /// no read in that stratum** — never a zero slip rate (spec §5).
    pub slippage_by_stratum_and_group: Vec<SlippageRow>,
    /// One row a stratum — the tract ladder's **top** rung. **Only a stratum fitted on its own
    /// tracts has one**: a stratum furnished from its period's slippage curves carries no length
    /// spectrum at all, by construction, and that absence is what the middle rung exists to
    /// answer. So absence here is data, not a hole.
    pub length_spectrum_by_stratum: Vec<StratumLengthSpectrumRow>,
    /// One row a motif period — the tract ladder's **middle** rung, present only where the run
    /// asked for it.
    pub length_spectrum_by_period: Vec<PeriodLengthSpectrumRow>,
    /// How often a base reads wrong inside a tract, per `(read group × stratum × ploidy)`.
    pub substitution_rate_by_stratum: Vec<SubstitutionRateRow>,
}

/// Which slippage group one read group's reads are drawn under.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SlippageGroupRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// The set of read groups this one slips alike with. The default is one group holding every
    /// read group.
    pub slippage_group: u32,
}

/// **One `(stratum × slippage group)`'s three slippage numbers, and where each came from.**
///
/// A level fitted from 8,000 slipped reads and one read off a curve through four cells are the
/// same `f64`, and a consumer that weighs them alike is treating an interpolation as a
/// measurement — which is why how much curve went into each number travels beside it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SlippageRow {
    /// How many bases one repeat unit of this stratum is.
    pub period: u8,
    /// How many copies of it the reference carries — the bin this stratum is.
    pub reference_repeats: u64,
    /// Which slippage group these numbers are for.
    pub slippage_group: u32,
    /// How often a read reports a tract length other than its allele's.
    pub level: f64,
    /// Of the reads that slip, the share showing a **shorter** tract.
    pub shorter_share: f64,
    /// How fast two-repeat slips fall off against one-repeat slips.
    pub fall_off: f64,
    /// Where the level came from, and how much of this stratum's own evidence stood behind it.
    pub level_origin: LevelOrigin,
    /// Where the direction split and the fall-off came from. **Separate from the level's**,
    /// because the three numbers are smoothed on their own curves and a stratum can take its
    /// level from a curve while keeping its own shares.
    pub shares_origin: Option<SharesOrigin>,
}

/// Where a stratum's slippage **level** came from — the file's spelling of
/// [`LevelProvenance`](crate::ng::parameter_estimation::joint::ssr_fit::LevelProvenance).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LevelOrigin {
    /// This stratum's own fit, its period's curve, or a blend of the two — and, where a curve
    /// was used, which curve and whether this stratum sat inside its fitted range.
    pub smoothing: LevelSmoothing,
    /// How many of this stratum's own reads **its own fitted level** said slipped, and absent
    /// where the stratum has no level of its own because it borrowed. **Absent is not zero.**
    pub slipped_reads: Option<f64>,
}

/// Where a stratum's direction split and fall-off came from — the file's spelling of
/// [`SharesProvenance`](crate::ng::parameter_estimation::joint::ssr_fit::SharesProvenance).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SharesOrigin {
    /// How many of this stratum's own reads its own fitted level said slipped, and absent where
    /// nothing was fitted here. **Both shares are proportions over the reads that slipped**, so
    /// this one count sets how precisely the stratum holds either of them.
    ///
    /// **Written separately from [`LevelOrigin::slipped_reads`] and not yet known to be the same
    /// number.** Upstream the two fields carry near-identical documentation but state their
    /// absence conditions differently, and nothing here can settle which. Step C4's round trip
    /// on a real fit is where to compare them across every stratum: if they never differ, one of
    /// these belongs on [`SlippageRow`] instead of two here.
    pub slipped_reads: Option<f64>,
    /// Where the share of slipped reads showing a shorter tract came from.
    pub shorter_share_smoothing: ShareSmoothing,
    /// Where the fall-off came from.
    pub fall_off_smoothing: ShareSmoothing,
}

/// How much of a period's curve went into a stratum's slippage **level**, and which curve.
///
/// **The curve and the reach ride on the variant rather than beside it**, which is what stops
/// three meaningless states being writable: a reach with no curve (a claim about a fitted range
/// the file does not carry), a curve with no reach, and "its period's curve, whole" with no
/// curve recorded. In the file that is serde's own enum spelling: `smoothing = "this_stratum"`
/// for the state that carries nothing, and a single-key table for the two that do.
///
/// **Two enums rather than one shared with [`ShareSmoothing`]**, because the two carry different
/// curves. The earlier one-enum-for-both shape read better and could not hold the curve, and
/// holding the curve is what removes the illegal states.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LevelSmoothing {
    /// The stratum's own fit, with no curve in it.
    ThisStratum,
    /// Its period's curve, whole.
    ThisPeriodsCurve {
        /// The curve that supplied it.
        curve: SlippageCurve,
        /// Whether this stratum's repeat count sat inside the curve's fitted range. **A level
        /// held at a fitted end is under-stated in a known direction.**
        reach: CurveReach,
    },
    /// The two, weighed against each other.
    Blend {
        /// The share the curve carried, in `[0, 1]`.
        curve_weight: f64,
        /// The curve that was weighed in.
        curve: SlippageCurve,
        /// Whether this stratum's repeat count sat inside the curve's fitted range.
        reach: CurveReach,
    },
}

/// How much of a period's curve went into one of a stratum's two slippage **shares**, and which
/// curve. The share counterpart of [`LevelSmoothing`], and it holds the same three states.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShareSmoothing {
    /// The stratum's own fit, with no curve in it.
    ThisStratum,
    /// Its period's curve, whole.
    ThisPeriodsCurve {
        /// The curve that supplied it.
        curve: ShareCurve,
        /// Whether this stratum's repeat count sat inside the curve's fitted range. **A share
        /// held at a fitted end is the end stratum's answer**, not this stratum's.
        reach: CurveReach,
    },
    /// The two, weighed against each other.
    Blend {
        /// The share the curve carried, in `[0, 1]`.
        curve_weight: f64,
        /// The curve that was weighed in.
        curve: ShareCurve,
        /// Whether this stratum's repeat count sat inside the curve's fitted range.
        reach: CurveReach,
    },
}

/// Whether a stratum's repeat count sat inside the curve's fitted range — the file's spelling of
/// [`CurveReach`](crate::ng::parameter_estimation::joint::slippage_curve::CurveReach).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CurveReach {
    /// Inside the fitted range.
    Inside,
    /// Below the shortest stratum the curve was fitted on.
    BelowFitted,
    /// Above the longest.
    AboveFitted,
}

/// **One motif period's slippage-level curve** — the file's spelling of
/// [`SlippageCurve`](crate::ng::parameter_estimation::joint::slippage_curve::SlippageCurve).
///
/// It carries its own held-out error and how many cells stood behind it, so a consumer can tell
/// a curve through twenty-three cells from one through four.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SlippageCurve {
    /// Where between multiplying and adding the rise sits, **in `[0, 1]`**: 0 multiplies the
    /// level by a fixed factor per extra repeat, 1 adds a fixed amount. Nothing outside that
    /// range is a rise shape.
    pub rise_shape: f64,
    /// The curve's value at the low end of its fitted range.
    pub intercept: f64,
    /// How fast it rises across repeat count.
    pub slope: f64,
    /// The shortest repeat count it was fitted on.
    pub fitted_from_repeats: u64,
    /// The longest. Beyond `fitted_from_repeats ..= fitted_to_repeats` the level is held at the
    /// nearer end.
    pub fitted_to_repeats: u64,
    /// How far it missed the cells it was not fitted on.
    pub held_out_error: f64,
    /// How many cells stood behind it.
    pub cells: u64,
}

/// **One motif period's curve for one of the two slippage shares** — the file's spelling of
/// [`ShareCurve`](crate::ng::parameter_estimation::joint::share_curve::ShareCurve).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShareCurve {
    /// Flat, sloping, or turning once.
    pub shape: ShareShape,
    /// The fitted curve on the logit scale, read at `repeats`:
    /// `intercept + slope · (repeats − centre_repeats) + bend · (repeats − centre_repeats)²`.
    pub intercept: f64,
    /// See [`ShareCurve::intercept`]. Zero for a flat curve.
    pub slope: f64,
    /// See [`ShareCurve::intercept`]. Zero for anything but a turning curve.
    pub bend: f64,
    /// The repeat count the slope and the bend are measured from — the weighted mean of the
    /// strata that fed the curve.
    pub centre_repeats: f64,
    /// The lowest repeat count of a stratum that fed this curve. **A curve that does not depend
    /// on repeat count spans the whole axis**, so the two bottom rungs report every repeat count
    /// as inside.
    pub fitted_from_repeats: u64,
    /// The highest. Beyond `fitted_from_repeats ..= fitted_to_repeats` the share is held at the
    /// nearer end.
    pub fitted_to_repeats: u64,
    /// How far the curve landed from a stratum it had not seen, in logit units.
    pub held_out_error: f64,
    /// How many strata stood behind it. A curve through four strata and one through twenty-three
    /// are both curves, and a consumer must be able to tell them apart.
    pub strata: u64,
    /// Which rung of the fallback ladder produced this curve. **Not the same question as
    /// [`LevelSmoothing`]'s or [`ShareSmoothing`]'s**: those say how much curve went into a
    /// stratum's number, this says what the curve itself was fitted on.
    pub rung: ShareCurveRung,
}

/// Which rung of the fallback ladder a share curve came from — the file's spelling of
/// [`ShareCurveSource`](crate::ng::parameter_estimation::joint::share_curve::ShareCurveSource).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShareCurveRung {
    /// Fitted on this motif period's own strata, and scored.
    ThisPeriod,
    /// Fitted on this period's own strata, with too few of them to score the shape choice.
    ThisPeriodUnscored,
    /// Pooled from the other motif periods.
    OtherPeriods,
    /// A stated constant — no period had anything to fit.
    BuiltInDefault,
}

/// What shape a share curve takes — the file's spelling of
/// [`ShareShape`](crate::ng::parameter_estimation::joint::share_curve::ShareShape).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShareShape {
    /// One value for the whole motif period — the share does not move with repeat count.
    Flat,
    /// The share moves one way across repeat count and never turns back.
    Sloping,
    /// The share turns once — falls and then rises, or rises and then falls.
    Turning,
}

/// **One stratum's fitted length spectrum** — how its chromosomes are spread over tract lengths.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StratumLengthSpectrumRow {
    /// How many bases one repeat unit is.
    pub period: u8,
    /// How many copies of it the reference carries.
    pub reference_repeats: u64,
    /// How monomorphic the stratum's tracts are. Small means most tracts carry one length.
    pub concentration: f64,
    /// One share a whole-repeat offset from the **reference** tract length, running
    /// `-span ..= +span` — so the middle entry is the reference length itself, the count is odd
    /// and at least three, and the shares sum to one.
    pub shares_by_repeat_offset: Vec<f64>,
}

/// **One motif period's pooled length spectrum** — the tract ladder's middle rung.
///
/// **That the index is an offset is what makes pooling legitimate at all**: two strata of one
/// period sit at different absolute repeat counts, and it is only because every tract is
/// described relative to its own reference length that their evidence can be added up.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeriodLengthSpectrumRow {
    /// The motif length these tracts share.
    pub period: u8,
    /// How monomorphic the period's tracts are.
    pub concentration: f64,
    /// One share a whole-repeat offset from each tract's own reference length, on the same
    /// `-span ..= +span` convention as a stratum's own.
    pub shares_by_repeat_offset: Vec<f64>,
}

/// **How often a base reads wrong inside a tract**, for one read group at one stratum.
///
/// The ploidy is part of the key rather than read off the run, because it is part of the key the
/// fit stored these under — the set of genotypes each table's entries were scored against.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubstitutionRateRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// How many bases one repeat unit is.
    pub period: u8,
    /// How many copies of it the reference carries.
    pub reference_repeats: u64,
    /// How many genome copies these loci sit on.
    pub ploidy: u8,
    /// The rate itself, a probability. Its evidence count is in **bases compared** inside the
    /// tracts — `table.bases_compared()` (`parameter_estimation/ssr`) — and **not in reads**: a
    /// read spanning a stratum contributes one read and as many bases as it crosses.
    pub rate: WarrantedValue,
}

// ---------------------------------------------------------------------
// §3.8 — the constants no fit produces
// ---------------------------------------------------------------------

/// **The numbers no fit produces, written out rather than left in the binary.**
///
/// **Marking them soft is the point of writing them down** (spec §3.8). A number that appears in
/// a file the user can edit is a number the project has admitted is a guess; one that only
/// appears in the source reads as a decision.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StatedConstants {
    /// The share of repeat-tract reads that came from nowhere the model can explain. Inherited
    /// from the existing caller and never measured here.
    ///
    /// **It carries a warrant because it has two reachable states and spec §8 requires them
    /// apart.** The run's own value is `defaulted` — 0.01, from `likelihood/ssr.rs`. A value a
    /// person typed into the file is `supplied`, and spec §3.8 says a person editing it is the
    /// whole point of writing it down. Without the warrant a run reports an edited guess as the
    /// project's own constant, and spec §2.1's wholesale demotion of a mismatched file has
    /// nowhere to write itself for this one number. Its evidence count is absent: no fit
    /// produced it.
    pub repeat_tract_outlier_weight: WarrantedValue,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A plant whose name needs escaping and is not ASCII.
    ///
    /// **Sample names come from `@RG SM` and are whatever the sequencing centre typed.** Every
    /// string in an earlier version of this fixture was bare ASCII, so nothing here would have
    /// noticed a writer that mangled a quote or a non-ASCII byte — and step B replaces serde's
    /// writer, which escapes correctly, with one that has to.
    const AWKWARD_SAMPLE: &str = "Ailsa ‘Craig’ \"×2\"";

    /// A share the curve carried that needs all seventeen significant digits to recover.
    ///
    /// **One value in the fixture that a short-digit formatter cannot round-trip.** Spec §4 says
    /// whether the crate emits enough digits "has not been checked here"; step C3 is where that
    /// is settled, and this is only the fixture's part of it — with every float a short decimal,
    /// narrowing one of them from `f64` to `f32` left the emitted file byte-identical.
    const FULL_PRECISION: f64 = 1.0 / 3.0;

    /// A file with **every section non-empty and every shape used at least once**.
    ///
    /// What a shape fixture can be wrong about is leaving a variant or a row kind unexercised, so
    /// every variant of every enum appears either here or in
    /// `every_enum_variant_spells_as_the_file_says`. Three properties are deliberate:
    ///
    /// - **the two batching axes have different lengths** — three read groups over two samples —
    ///   so exchanging them cannot produce a file that still parses;
    /// - **no row has two equal numeric fields**, so a field-for-field comparison cannot pass on
    ///   two fields that were swapped. The exception is deliberate and named: read group 0 sits
    ///   in slippage group 0, and read group 2 sits in slippage group 1, so the map is not the
    ///   identity;
    /// - **one string needs escaping and one float needs full precision** (above).
    fn a_file_using_every_shape() -> ParametersFile {
        ParametersFile {
            format_version: FORMAT_VERSION,
            ploidy: 2,
            fitted_from: InputsFittedFrom {
                reference_digest: "0123456789abcdef".into(),
                samples: vec!["TS-1".into(), AWKWARD_SAMPLE.into()],
                read_groups: vec![
                    ReadGroupRow {
                        read_group: 0,
                        declared_id: "HWI.3".into(),
                        library: "lib3".into(),
                        sample: "TS-1".into(),
                    },
                    // The second library of the *same* plant — the grain the contamination
                    // fraction exists at, and the case a per-sample row would erase.
                    ReadGroupRow {
                        read_group: 1,
                        declared_id: "HWI.4".into(),
                        library: "lib4".into(),
                        sample: "TS-1".into(),
                    },
                    ReadGroupRow {
                        read_group: 2,
                        declared_id: "HWI.5".into(),
                        library: "lib5".into(),
                        sample: AWKWARD_SAMPLE.into(),
                    },
                ],
                census: CensusIdentity {
                    terms: vec![CensusTerm {
                        term: "the loci actually kept".into(),
                        digest: "fedcba9876543210".into(),
                    }],
                },
            },
            base_quality_calibration: BaseQualityCalibration {
                by_read_group: vec![
                    BaseQualityCalibrationRow {
                        read_group: 0,
                        error_probability_multiplier: WarrantedValue {
                            value: 1.0324,
                            warrant: Warrant::FittedHere,
                            observations: Some(EvidenceCount::Reads(812_344)),
                        },
                    },
                    // A multiplier of exactly one that was *not* fitted — the state the warrant
                    // exists to keep apart from a fitted one, and the read count says why: four
                    // reads is no rate to fit.
                    BaseQualityCalibrationRow {
                        read_group: 1,
                        error_probability_multiplier: WarrantedValue {
                            value: 1.0,
                            warrant: Warrant::Defaulted,
                            observations: Some(EvidenceCount::Reads(4)),
                        },
                    },
                    BaseQualityCalibrationRow {
                        read_group: 2,
                        error_probability_multiplier: WarrantedValue {
                            value: 0.87,
                            warrant: Warrant::Supplied,
                            observations: Some(EvidenceCount::Reads(640_918)),
                        },
                    },
                ],
            },
            contamination: Some(Contamination {
                by_read_group: vec![
                    ContaminationRow {
                        read_group: 0,
                        library: "lib3".into(),
                        measurement: Some(ContaminationMeasurement {
                            fraction: 0.031,
                            markers_with_reads: 4211,
                            reads_on_markers: 90233,
                            fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                        }),
                    },
                    // **Not measured**, inside a run where the others were — the same plant's
                    // second lane, which put too few reads on the panel's markers. It has no
                    // fraction to be misread as clean, and no `fitted_from_reads_of`, which
                    // would be untrue of it either way.
                    ContaminationRow {
                        read_group: 1,
                        library: "lib4".into(),
                        measurement: None,
                    },
                    // **Measured and found clean**: a fraction of zero with evidence behind it —
                    // spec §5's second row, and the state the row above is not.
                    ContaminationRow {
                        read_group: 2,
                        library: "lib5".into(),
                        measurement: Some(ContaminationMeasurement {
                            fraction: 0.0,
                            markers_with_reads: 2903,
                            reads_on_markers: 64118,
                            fitted_from_reads_of: ContaminationFittedFrom::EveryReadOfThisSample,
                        }),
                    },
                ],
            }),
            sequencing_batches: SequencingBatches {
                was_declared_by_the_run: true,
                by_read_group: vec![
                    ReadGroupBatchRow {
                        read_group: 0,
                        batch: 0,
                    },
                    ReadGroupBatchRow {
                        read_group: 1,
                        batch: 1,
                    },
                    ReadGroupBatchRow {
                        read_group: 2,
                        batch: 1,
                    },
                ],
                by_sample: vec![
                    SampleBatchRow {
                        sample: "TS-1".into(),
                        batch: 0,
                    },
                    SampleBatchRow {
                        sample: AWKWARD_SAMPLE.into(),
                        batch: 1,
                    },
                ],
            },
            inbreeding: Inbreeding {
                by_sample: vec![
                    InbreedingRow {
                        sample: "TS-1".into(),
                        inbreeding_coefficient: WarrantedValue {
                            value: 0.42,
                            warrant: Warrant::FittedHere,
                            observations: Some(EvidenceCount::CoveredPositions(180_600_412)),
                        },
                    },
                    InbreedingRow {
                        sample: AWKWARD_SAMPLE.into(),
                        inbreeding_coefficient: WarrantedValue {
                            value: 0.17,
                            warrant: Warrant::Borrowed,
                            observations: Some(EvidenceCount::CoveredPositions(9_411_027)),
                        },
                    },
                ],
            },
            ordinary_site_prior: OrdinarySitePrior {
                reference_concentration: 1.0,
                alternative_concentration_total: 0.000_6,
                rung: SeedRung::FittedCurve,
            },
            repeat_tracts: RepeatTracts {
                stated_length_spectrum_concentration: WarrantedValue {
                    value: 1.25,
                    warrant: Warrant::Defaulted,
                    observations: None,
                },
                slippage_group_by_read_group: vec![
                    SlippageGroupRow {
                        read_group: 0,
                        slippage_group: 0,
                    },
                    SlippageGroupRow {
                        read_group: 1,
                        slippage_group: 0,
                    },
                    // Not the identity map, so a row whose two fields were exchanged says
                    // something different rather than the same thing.
                    SlippageGroupRow {
                        read_group: 2,
                        slippage_group: 1,
                    },
                ],
                slippage_by_stratum_and_group: vec![
                    SlippageRow {
                        period: 2,
                        reference_repeats: 6,
                        slippage_group: 0,
                        level: 0.0421,
                        shorter_share: 0.83,
                        fall_off: 0.31,
                        level_origin: LevelOrigin {
                            smoothing: LevelSmoothing::Blend {
                                curve_weight: 0.37,
                                curve: a_slippage_curve(),
                                reach: CurveReach::Inside,
                            },
                            slipped_reads: Some(8_000.5),
                        },
                        shares_origin: Some(SharesOrigin {
                            slipped_reads: Some(8_000.5),
                            shorter_share_smoothing: ShareSmoothing::ThisStratum,
                            fall_off_smoothing: ShareSmoothing::ThisPeriodsCurve {
                                curve: a_share_curve(),
                                reach: CurveReach::AboveFitted,
                            },
                        }),
                    },
                    SlippageRow {
                        period: 2,
                        reference_repeats: 11,
                        slippage_group: 1,
                        level: 0.0913,
                        shorter_share: 0.79,
                        fall_off: 0.28,
                        level_origin: LevelOrigin {
                            smoothing: LevelSmoothing::ThisPeriodsCurve {
                                curve: a_slippage_curve(),
                                reach: CurveReach::BelowFitted,
                            },
                            // Absent, not zero: this stratum borrowed and has no level of its
                            // own to count slipped reads against.
                            slipped_reads: None,
                        },
                        shares_origin: None,
                    },
                    SlippageRow {
                        period: 1,
                        reference_repeats: 30,
                        slippage_group: 0,
                        level: 0.19,
                        shorter_share: 0.77,
                        fall_off: 0.24,
                        level_origin: LevelOrigin {
                            smoothing: LevelSmoothing::ThisStratum,
                            slipped_reads: Some(12_040.25),
                        },
                        shares_origin: Some(SharesOrigin {
                            slipped_reads: Some(12_040.25),
                            shorter_share_smoothing: ShareSmoothing::Blend {
                                curve_weight: 0.61,
                                curve: a_share_curve(),
                                reach: CurveReach::Inside,
                            },
                            fall_off_smoothing: ShareSmoothing::ThisStratum,
                        }),
                    },
                ],
                length_spectrum_by_stratum: vec![StratumLengthSpectrumRow {
                    period: 2,
                    reference_repeats: 6,
                    concentration: 3.5,
                    shares_by_repeat_offset: vec![0.1, 0.8, 0.1],
                }],
                length_spectrum_by_period: vec![PeriodLengthSpectrumRow {
                    period: 2,
                    concentration: 2.75,
                    shares_by_repeat_offset: vec![0.15, 0.7, 0.15],
                }],
                substitution_rate_by_stratum: vec![SubstitutionRateRow {
                    read_group: 0,
                    period: 2,
                    reference_repeats: 6,
                    ploidy: 2,
                    rate: WarrantedValue {
                        value: 0.0012,
                        warrant: Warrant::Borrowed,
                        observations: Some(EvidenceCount::BasesCompared(40_122)),
                    },
                }],
            },
            stated_constants: StatedConstants {
                repeat_tract_outlier_weight: WarrantedValue {
                    value: 0.01,
                    warrant: Warrant::Defaulted,
                    observations: None,
                },
            },
        }
    }

    fn a_slippage_curve() -> SlippageCurve {
        SlippageCurve {
            rise_shape: 0.55,
            intercept: 0.011,
            slope: 0.004,
            fitted_from_repeats: 5,
            fitted_to_repeats: 19,
            held_out_error: FULL_PRECISION,
            cells: 23,
        }
    }

    fn a_share_curve() -> ShareCurve {
        ShareCurve {
            shape: ShareShape::Turning,
            intercept: 1.4,
            slope: -0.09,
            bend: 0.006,
            centre_repeats: 11.5,
            fitted_from_repeats: 5,
            fitted_to_repeats: 19,
            held_out_error: 0.167,
            strata: 12,
            rung: ShareCurveRung::ThisPeriod,
        }
    }

    /// **The shape comes back the same object through `serde`.**
    ///
    /// It uses `serde`'s own serializer and **not** the file's writer, which does not exist yet.
    /// What it pins is that the tree holds no enum shape TOML cannot spell and no field that
    /// reads back as a different one.
    ///
    /// **It does not pin key names, and it cannot.** Writing and reading through the same derive
    /// moves both sides of any rename together, so renaming a field or dropping a `rename_all`
    /// leaves this green while the file on disk changes — measured, by a review that renamed
    /// `markers_with_reads` and saw every test pass.
    /// [`the_whole_shape_emits_the_documented_toml`] is what holds the key surface.
    ///
    /// **It does not pin field order either.** `toml` 1.1 emits a struct's scalar fields before
    /// its table-valued ones whatever the declared order, so no ordering of this tree can
    /// produce the "value after table" error; moving `ploidy` to the last field of
    /// [`ParametersFile`] leaves the emitted bytes identical. Whether the file *opens* with
    /// `format_version` and `ploidy` is step B's writer to decide and step B's test to hold.
    #[test]
    fn every_section_of_the_shape_survives_a_serde_round_trip() {
        let written = a_file_using_every_shape();
        let text = toml::to_string(&written).expect("the shape is expressible as TOML");
        let read: ParametersFile = toml::from_str(&text).expect("and parses back");
        assert_eq!(read, written);
    }

    /// **Every key and every spelling of the emitted file, pinned against a checked-in copy.**
    ///
    /// This is the test that makes the file's names a compatibility surface rather than an
    /// intention: a renamed field, a dropped `rename_all`, a field added or removed, or a
    /// changed nesting all fail here and nowhere else.
    ///
    /// **The golden file is `serde`'s output and not the artefact a run will write** — the
    /// hand-written writer of step B2 emits one-row-a-line inline tables and carries comments,
    /// and owes a golden file of its own. What this one holds in the meantime is the *names*,
    /// which are the same under either writer.
    #[test]
    fn the_whole_shape_emits_the_documented_toml() {
        let text = toml::to_string(&a_file_using_every_shape()).expect("serialises");
        assert_eq!(
            text,
            include_str!("testdata/every_shape.toml"),
            "the emitted file no longer matches testdata/every_shape.toml; if the change is \
             intended, regenerate that file from this fixture"
        );
    }

    /// **Rewrite `testdata/every_shape.toml` from the fixture.**
    ///
    /// Ignored, so it never runs in an ordinary suite: it makes
    /// [`the_whole_shape_emits_the_documented_toml`] pass by definition, which is the one thing
    /// that test must not do on its own. Run it deliberately, after an intended change to the
    /// shape, and read the resulting diff:
    ///
    /// ```text
    /// cargo test --lib ng::calling::parameters_file -- --ignored regenerate
    /// ```
    #[test]
    #[ignore = "rewrites the golden file; run deliberately after an intended shape change"]
    fn regenerate_the_golden_file() {
        let text = toml::to_string(&a_file_using_every_shape()).expect("serialises");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src/ng/calling/parameters_file/testdata/every_shape.toml");
        std::fs::write(&path, text).expect("the golden file is writable");
    }

    /// The file's spelling of one unit-variant enum.
    fn spelling<T: Serialize>(value: T) -> String {
        toml::Value::try_from(value)
            .expect("a unit variant is a TOML value")
            .as_str()
            .expect("a unit variant spells as a bare string")
            .to_owned()
    }

    /// **Every variant spells as the file says it does.**
    ///
    /// That spelling is the compatibility surface these enums exist to own: renaming a Rust
    /// variant must not silently re-interpret a file on disk. Neither the round trip nor the
    /// golden file above can stand in for this one — the round trip moves both sides of a rename
    /// at once, and the golden file only sees the variants the fixture happens to use, which is
    /// eleven of these twenty-one.
    #[test]
    fn every_enum_variant_spells_as_the_file_says() {
        assert_eq!(spelling(Warrant::FittedHere), "fitted_here");
        assert_eq!(spelling(Warrant::Borrowed), "borrowed");
        assert_eq!(spelling(Warrant::Supplied), "supplied");
        assert_eq!(spelling(Warrant::Defaulted), "defaulted");

        assert_eq!(
            spelling(ContaminationFittedFrom::ThisReadGroupsOwnReads),
            "this_read_groups_own_reads"
        );
        assert_eq!(
            spelling(ContaminationFittedFrom::EveryReadOfThisSample),
            "every_read_of_this_sample"
        );

        assert_eq!(spelling(SeedRung::FittedCurve), "fitted_curve");
        assert_eq!(spelling(SeedRung::NeutralShape), "neutral_shape");
        assert_eq!(
            spelling(SeedRung::StatedHeterozygosity),
            "stated_heterozygosity"
        );

        assert_eq!(spelling(CurveReach::Inside), "inside");
        assert_eq!(spelling(CurveReach::BelowFitted), "below_fitted");
        assert_eq!(spelling(CurveReach::AboveFitted), "above_fitted");

        assert_eq!(spelling(ShareCurveRung::ThisPeriod), "this_period");
        assert_eq!(
            spelling(ShareCurveRung::ThisPeriodUnscored),
            "this_period_unscored"
        );
        assert_eq!(spelling(ShareCurveRung::OtherPeriods), "other_periods");
        assert_eq!(spelling(ShareCurveRung::BuiltInDefault), "built_in_default");

        assert_eq!(spelling(ShareShape::Flat), "flat");
        assert_eq!(spelling(ShareShape::Sloping), "sloping");
        assert_eq!(spelling(ShareShape::Turning), "turning");

        assert_eq!(spelling(LevelSmoothing::ThisStratum), "this_stratum");
        assert_eq!(spelling(ShareSmoothing::ThisStratum), "this_stratum");
        // `EvidenceCount`'s three variants carry a number, so they spell as a one-key table
        // rather than a bare string; `an_evidence_count_names_its_unit_and_is_absent_where_there_is_none`
        // is what pins them, and the fixture uses all three.
    }

    /// **A smoothing that used a curve carries the curve and the reach; one that did not carries
    /// neither.**
    ///
    /// Three states that were writable when the curve, the reach and the source were three
    /// independent fields are not writable now: a reach with no curve — a claim about a fitted
    /// range the file does not carry — a curve with no reach, and "its period's curve, whole"
    /// with no curve recorded. An earlier version of this module's own fixture wrote the first
    /// of the three, and every test passed on it.
    #[test]
    fn a_smoothing_that_used_a_curve_carries_it_and_a_plain_one_carries_neither() {
        let plain = toml::Value::try_from(LevelSmoothing::ThisStratum).expect("a smoothing");
        assert_eq!(plain.as_str(), Some("this_stratum"));

        let blended = toml::Value::try_from(LevelSmoothing::Blend {
            curve_weight: 0.37,
            curve: a_slippage_curve(),
            reach: CurveReach::Inside,
        })
        .expect("a blended smoothing");
        let blend = blended.get("blend").expect("a blend writes one key");
        assert_eq!(
            blend.get("curve_weight").and_then(toml::Value::as_float),
            Some(0.37)
        );
        assert_eq!(
            blend.get("reach").and_then(toml::Value::as_str),
            Some("inside")
        );
        assert_eq!(
            blend
                .get("curve")
                .and_then(|curve| curve.get("cells"))
                .and_then(toml::Value::as_integer),
            Some(23)
        );
    }

    /// **An absent shares origin writes no key, and a present one writes it.**
    ///
    /// Spec §5's rule is that `Option<T>` is absence and never a sentinel, and the failure it
    /// guards against needs both halves: a test that only checks the absent case passes just as
    /// well for a field that is never written at all. Asserted against the parsed value's keys
    /// rather than against the rendered text, because a substring search for `curve` also
    /// matches `curve_weight` and every `_curve` in a variant's spelling.
    #[test]
    fn an_absent_shares_origin_writes_no_key_and_a_present_one_does() {
        let mut row = a_file_using_every_shape()
            .repeat_tracts
            .slippage_by_stratum_and_group[0]
            .clone();

        row.shares_origin = None;
        let value = toml::Value::try_from(&row).expect("a slippage row is a TOML value");
        assert!(
            value.get("shares_origin").is_none(),
            "an absent shares origin leaves no key, got: {value}"
        );
        assert!(
            value.get("level_origin").is_some(),
            "the level origin is not optional and is always written"
        );

        row.shares_origin = Some(SharesOrigin {
            slipped_reads: None,
            shorter_share_smoothing: ShareSmoothing::ThisStratum,
            fall_off_smoothing: ShareSmoothing::ThisStratum,
        });
        let value = toml::Value::try_from(&row).expect("a slippage row is a TOML value");
        let shares = value
            .get("shares_origin")
            .expect("a present shares origin writes its key");
        assert!(
            shares.get("slipped_reads").is_none(),
            "and its own absent count leaves no key, got: {shares}"
        );
        assert_eq!(
            shares
                .get("shorter_share_smoothing")
                .and_then(toml::Value::as_str),
            Some("this_stratum")
        );
    }

    /// Every table anywhere in the document that carries a `warrant` key, with the path it was
    /// found at.
    fn tables_carrying_a_warrant(
        value: &toml::Value,
        path: &str,
        found: &mut Vec<(String, toml::Table)>,
    ) {
        match value {
            toml::Value::Table(table) => {
                if table.contains_key("warrant") {
                    found.push((path.to_owned(), table.clone()));
                }
                for (key, child) in table {
                    tables_carrying_a_warrant(child, &format!("{path}.{key}"), found);
                }
            }
            toml::Value::Array(items) => {
                for (index, item) in items.iter().enumerate() {
                    tables_carrying_a_warrant(item, &format!("{path}[{index}]"), found);
                }
            }
            _ => {}
        }
    }

    /// The same path with every `[n]` array index removed, so an assertion about *which fields*
    /// are warranted does not move when the fixture gains a row.
    fn without_array_indices(path: &str) -> String {
        let mut stripped = path.to_owned();
        while let Some(open) = stripped.find('[') {
            let close = stripped[open..].find(']').expect("a path index is closed") + open;
            stripped.replace_range(open..=close, "");
        }
        stripped
    }

    /// **Every warranted number in the file is written the same way, and these five are all of
    /// them** — which is the whole of what this shape is for: a reader who has understood one has
    /// understood all of them.
    ///
    /// It searches the emitted document for *every* table carrying a `warrant` key rather than
    /// visiting the five by name, and **the two assertions catch different things**. A number
    /// written flat with its key spelled exactly `warrant` lands in a table that has no `value`,
    /// and the loop catches it. One written flat under a *suffixed* key —
    /// `stated_length_spectrum_warrant`, the spelling this shape replaced — is invisible to the
    /// walk, and only the list of paths catches it. A sixth warranted number, or a fifth that
    /// lost its warrant, is caught by the list too.
    #[test]
    fn every_warranted_number_is_written_the_same_way() {
        let document =
            toml::Value::try_from(a_file_using_every_shape()).expect("a file is a TOML value");
        let mut found = Vec::new();
        tables_carrying_a_warrant(&document, "", &mut found);
        assert!(!found.is_empty(), "the walk found nothing at all");

        let mut fields: Vec<String> = found
            .iter()
            .map(|(path, _)| without_array_indices(path))
            .collect();
        fields.sort();
        fields.dedup();
        assert_eq!(
            fields,
            [
                ".base_quality_calibration.by_read_group.error_probability_multiplier",
                ".inbreeding.by_sample.inbreeding_coefficient",
                ".repeat_tracts.stated_length_spectrum_concentration",
                ".repeat_tracts.substitution_rate_by_stratum.rate",
                ".stated_constants.repeat_tract_outlier_weight",
            ],
            "these five numbers carry a warrant and nothing else does"
        );

        for (path, table) in &found {
            assert!(
                table.get("value").and_then(toml::Value::as_float).is_some(),
                "{path} carries a warrant and no float value, so it is written flat rather than \
                 as a warranted number: {table:?}"
            );
            for key in table.keys() {
                assert!(
                    ["value", "warrant", "observations"].contains(&key.as_str()),
                    "{path} carries the key {key}, which is not one of the three a warranted \
                     number is written with"
                );
            }
        }
    }

    /// **An evidence count names the unit it counts in, and is absent rather than zero where no
    /// fit produced one.**
    ///
    /// The unit is what makes two counts comparable, and it is not recoverable from the number:
    /// on one cohort the same sample's covered positions and its reads differ by orders of
    /// magnitude. This module's own doc comment gave the wrong unit for two of the three
    /// quantities before a review checked them against the code that produces them, which is why
    /// the unit is in the file rather than in the comment.
    #[test]
    fn an_evidence_count_names_its_unit_and_is_absent_where_there_is_none() {
        let document =
            toml::Value::try_from(a_file_using_every_shape()).expect("a file is a TOML value");
        let mut found = Vec::new();
        tables_carrying_a_warrant(&document, "", &mut found);

        let unit_of = |needle: &str| -> Option<String> {
            let (_, table) = found
                .iter()
                .find(|(path, _)| path.contains(needle))
                .unwrap_or_else(|| panic!("a {needle} row"));
            table.get("observations").map(|count| {
                count
                    .as_table()
                    .expect("an evidence count names its unit, so it is a table")
                    .keys()
                    .next()
                    .expect("and the unit is its one key")
                    .clone()
            })
        };

        assert_eq!(
            unit_of("base_quality_calibration").as_deref(),
            Some("reads")
        );
        assert_eq!(
            unit_of("inbreeding").as_deref(),
            Some("covered_positions"),
            "an inbreeding coefficient is fitted over covered positions, not over windows"
        );
        assert_eq!(
            unit_of("substitution_rate").as_deref(),
            Some("bases_compared"),
            "a repeat-tract substitution rate is fitted over bases compared, not over reads"
        );

        // Absent, not zero: neither of these came from a fit with a sample size.
        assert_eq!(unit_of("stated_length_spectrum_concentration"), None);
        assert_eq!(unit_of("repeat_tract_outlier_weight"), None);
    }

    /// **Each of spec §5's five states is a missing key, and each is distinct from the value a
    /// reader might otherwise collapse it into.**
    ///
    /// Distinctness in the *types* is what Milestone A owes; that collapsing two of them **changes
    /// an answer** is step C5's, on fixtures built for it, and that they survive a round trip is
    /// step C4's. What this holds is the shape: for every one of the five, the absent form writes
    /// no key at all, and the form it is confusable with writes a different document.
    #[test]
    fn each_of_the_five_states_is_a_missing_key_and_not_a_value() {
        let fitted = a_file_using_every_shape();
        let document =
            |file: &ParametersFile| toml::Value::try_from(file).expect("a file is a TOML value");

        // 1. No read group identified any contamination: the section is gone, not zeroed.
        let mut uncontaminated = fitted.clone();
        uncontaminated.contamination = None;
        assert!(
            document(&uncontaminated).get("contamination").is_none(),
            "an uncontaminated run writes no contamination section"
        );
        assert!(
            document(&fitted).get("contamination").is_some(),
            "and a run that measured one writes it"
        );

        // 2. Measured and found clean against not measured. Both are near-zero fractions in
        //    memory; here one has a measurement and the other has none.
        let rows = &fitted
            .contamination
            .as_ref()
            .expect("a table")
            .by_read_group;
        let clean = rows
            .iter()
            .find(|row| {
                row.measurement
                    .as_ref()
                    .is_some_and(|found| found.fraction == 0.0)
            })
            .expect("a read group measured and found clean");
        let unmeasured = rows
            .iter()
            .find(|row| row.measurement.is_none())
            .expect("a read group that identified nothing");
        assert!(
            toml::Value::try_from(unmeasured)
                .expect("a row")
                .get("measurement")
                .is_none(),
            "an unmeasured read group writes no measurement"
        );
        assert_ne!(
            toml::Value::try_from(clean).expect("a row"),
            toml::Value::try_from(unmeasured).expect("a row")
        );

        // 3. A multiplier of exactly one that was defaulted against one that was fitted.
        let mut defaulted = fitted.clone();
        let mut as_fitted = fitted.clone();
        for (file, warrant) in [
            (&mut defaulted, Warrant::Defaulted),
            (&mut as_fitted, Warrant::FittedHere),
        ] {
            file.base_quality_calibration.by_read_group[0].error_probability_multiplier =
                WarrantedValue {
                    value: 1.0,
                    warrant,
                    observations: None,
                };
        }
        assert_ne!(
            document(&defaulted),
            document(&as_fitted),
            "a defaulted multiplier of 1.0 is not a fitted one"
        );

        // 4 and 5 are gaps the fixture already has, so both are asserted against the emitted
        // document as it stands. **Emptying a table and checking it came out empty is a weaker
        // claim** — it holds for a writer that fills in a row for every stratum it knows about,
        // which is the collapse these two states forbid.
        let table_of = |section: &str| -> Vec<toml::Value> {
            document(&fitted)
                .get("repeat_tracts")
                .and_then(|tracts| tracts.get(section))
                .and_then(toml::Value::as_array)
                .unwrap_or_else(|| panic!("the {section} table"))
                .clone()
        };
        let stratum_of = |row: &toml::Value| -> (i64, i64) {
            (
                row.get("period")
                    .and_then(toml::Value::as_integer)
                    .expect("a period"),
                row.get("reference_repeats")
                    .and_then(toml::Value::as_integer)
                    .expect("a reference repeat count"),
            )
        };

        // 4. A stratum furnished from its period's curves has no length spectrum of its own. The
        //    fixture carries slippage for three strata and a spectrum for one of them.
        let spectra = table_of("length_spectrum_by_stratum");
        let with_a_spectrum: Vec<(i64, i64)> = spectra.iter().map(stratum_of).collect();
        assert_eq!(
            with_a_spectrum,
            [(2, 6)],
            "one stratum was fitted on its own tracts and writes a spectrum"
        );
        for furnished in [(2, 11), (1, 30)] {
            assert!(
                !with_a_spectrum.contains(&furnished),
                "a stratum furnished from its period's curves has no row at all, rather than a \
                 flat one: {furnished:?}"
            );
        }

        // 5. A slippage group that put no read in a stratum has no row for that pair. Three
        //    strata crossed with two slippage groups is six pairs; three of them were fitted.
        let slippage = table_of("slippage_by_stratum_and_group");
        let written: Vec<(i64, i64, i64)> = slippage
            .iter()
            .map(|row| {
                let (period, repeats) = stratum_of(row);
                (
                    period,
                    repeats,
                    row.get("slippage_group")
                        .and_then(toml::Value::as_integer)
                        .expect("a slippage group"),
                )
            })
            .collect();
        assert_eq!(
            written.len(),
            3,
            "a row exists only where a group put reads"
        );
        assert!(
            written.contains(&(1, 30, 0)),
            "the pair that was fitted has its row"
        );
        assert!(
            !written.contains(&(1, 30, 1)),
            "the pair with no reads has no row at all, rather than a row of zeros"
        );
    }

    /// **Two shapes the file can spell that no run should mean, recorded rather than endorsed.**
    ///
    /// Both say something spec §5 requires to be said another way, and **neither can be made
    /// unspellable by a shape**: no type can say a `Vec` is non-empty, or that not every row of
    /// one is absent, or that two counts are not both zero without refusing an estimate the fit
    /// might legitimately produce. So they are step C2's to refuse — the step that takes the file
    /// to `RunParameters` and meets the run's dense read-group axis.
    ///
    /// **This test exists so the gap is visible rather than implied.** When C2 lands, these
    /// assertions invert: what round-trips today must then fail, naming the field.
    #[test]
    fn the_shape_accepts_two_things_step_c2_must_refuse() {
        // A contamination table in which nobody was measured — the uncontaminated run, written
        // longhand. Read literally it takes the mixture path with every fraction zero, where
        // absence takes the read likelihood's plain formula.
        let mut longhand = a_file_using_every_shape();
        for row in &mut longhand
            .contamination
            .as_mut()
            .expect("a table")
            .by_read_group
        {
            row.measurement = None;
        }
        let text = toml::to_string(&longhand).expect("serialises");
        assert!(
            toml::from_str::<ParametersFile>(&text).is_ok(),
            "if this now fails, step C2's refusal has landed — invert this half"
        );

        // A measurement carrying the evidence of not being measured: the in-memory
        // `UNMEASURED_READ_GROUP` shape, which a projection written from the view rather than
        // from the fit's estimate would produce.
        let mut evidenceless = a_file_using_every_shape();
        evidenceless
            .contamination
            .as_mut()
            .expect("a table")
            .by_read_group[0]
            .measurement = Some(ContaminationMeasurement {
            fraction: 0.0,
            markers_with_reads: 0,
            reads_on_markers: 0,
            fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
        });
        let text = toml::to_string(&evidenceless).expect("serialises");
        assert!(
            toml::from_str::<ParametersFile>(&text).is_ok(),
            "if this now fails, step C2's refusal has landed — invert this half"
        );
    }

    /// **A run with nothing in any table still writes a file that parses.**
    ///
    /// The empty boundary, which no other fixture reaches: a single sample with no repeat tracts
    /// is inside the committed input range (`CLAUDE.md`), and an array of tables that serialises
    /// to nothing and then fails to parse would surface only there.
    #[test]
    fn a_file_with_every_table_empty_round_trips() {
        let mut empty = a_file_using_every_shape();
        empty.fitted_from.samples.clear();
        empty.fitted_from.read_groups.clear();
        empty.fitted_from.census.terms.clear();
        empty.base_quality_calibration.by_read_group.clear();
        // Absence, not an empty table: an uncontaminated run has no section at all.
        empty.contamination = None;
        empty.sequencing_batches.by_read_group.clear();
        empty.sequencing_batches.by_sample.clear();
        empty.inbreeding.by_sample.clear();
        empty.repeat_tracts.slippage_group_by_read_group.clear();
        empty.repeat_tracts.slippage_by_stratum_and_group.clear();
        empty.repeat_tracts.length_spectrum_by_stratum.clear();
        empty.repeat_tracts.length_spectrum_by_period.clear();
        empty.repeat_tracts.substitution_rate_by_stratum.clear();

        let text = toml::to_string(&empty).expect("an empty file serialises");
        let read: ParametersFile = toml::from_str(&text).expect("and parses back");
        assert_eq!(read, empty);
    }

    /// **A stated concentration of one reads back saying whether it was fitted.**
    ///
    /// The run's own median over the strata it fitted can land on exactly 1.0, which is also the
    /// constant reached where it fitted none, and spec §8 names this among the three numbers that
    /// must be marked `defaulted` when they are. **The assertion is on what a reader recovers**,
    /// not on the two documents differing: since `warrant` is a required field of
    /// [`WarrantedValue`], two files whose warrants differ differ by construction, and a test
    /// asserting only that would pass for any shape at all.
    #[test]
    fn a_stated_concentration_of_one_reads_back_saying_whether_it_was_fitted() {
        for warrant in [Warrant::FittedHere, Warrant::Defaulted] {
            let mut file = a_file_using_every_shape();
            file.repeat_tracts.stated_length_spectrum_concentration = WarrantedValue {
                value: 1.0,
                warrant,
                observations: None,
            };
            let text = toml::to_string(&file).expect("serialises");
            let read: ParametersFile = toml::from_str(&text).expect("parses");
            assert_eq!(
                read.repeat_tracts
                    .stated_length_spectrum_concentration
                    .warrant,
                warrant,
                "a concentration of exactly 1.0 must not lose which of the two it was"
            );
        }
    }

    /// **A mistyped key is refused rather than absorbed.**
    ///
    /// The hazard is the optional fields: serde's ordinary behaviour discards an unrecognised
    /// key in silence, so `[…level_origin.smoothin]` would parse and leave `smoothing` — or, for
    /// an `Option`, an absence that is data. `deny_unknown_fields` is a type attribute with no
    /// call-site knob, so this cannot be a property the reader turns on later.
    #[test]
    fn a_mistyped_key_is_refused_rather_than_absorbed() {
        let text = toml::to_string(&a_file_using_every_shape()).expect("serialises");

        // Inserted **into** an existing table rather than under a new header: re-opening a
        // table TOML has already seen is a parse error whatever the keys are, so a test written
        // that way would pass without `deny_unknown_fields` doing anything. `replacen` and not
        // `replace`, because the fixture emits `ploidy = 2` twice — once at the top level and
        // once on a substitution-rate row — and only the first is meant to be sat beside.
        let extra = text.replacen("ploidy = 2", "ploidy = 2\nploidee = 3", 1);
        assert!(
            extra != text,
            "the fixture writes a ploidy for this test to sit beside"
        );
        let refusal = toml::from_str::<ParametersFile>(&extra)
            .expect_err("a misspelled key must not be absorbed")
            .to_string();
        assert!(
            refusal.contains("ploidee"),
            "and the refusal names the key it did not know, got: {refusal}"
        );

        let typoed = text.replace("slipped_reads", "sliped_reads");
        assert!(
            typoed != text,
            "the fixture writes a slipped-read count for this test to misspell"
        );
        assert!(
            toml::from_str::<ParametersFile>(&typoed).is_err(),
            "a misspelled optional key must not read back as absence"
        );
    }

    /// **The inline form the module documents parses**, which is the claim steps B and C meet on.
    ///
    /// The hand-written writer of step B2 emits one row a line as an inline table; `serde`'s
    /// serializer emits `[[array of tables]]` headers instead, so every other test here parses
    /// the one shape that writer will *not* produce. The nesting this reaches is four deep —
    /// `level_origin.smoothing.blend.curve` — which is the part worth knowing about before both
    /// the writer and the reader are built.
    #[test]
    fn the_documented_inline_form_parses() {
        let text = r#"
format_version = 1
ploidy = 2

[fitted_from]
reference_digest = "abc"
samples = ["TS-1"]
read_groups = [ { read_group = 0, declared_id = "HWI.3", library = "lib3", sample = "TS-1" } ]
census = { terms = [ { term = "the loci actually kept", digest = "def" } ] }

[base_quality_calibration]
by_read_group = [ { read_group = 0, error_probability_multiplier = { value = 1.0, warrant = "defaulted", observations = { reads = 4 } } } ]

# No [contamination] section at all: this run is uncontaminated, which is spec §5's first state
# and a different claim from a table of zeros.

[sequencing_batches]
was_declared_by_the_run = false
by_read_group = [ { read_group = 0, batch = 0 } ]
by_sample = [ { sample = "TS-1", batch = 0 } ]

[inbreeding]
by_sample = [ { sample = "TS-1", inbreeding_coefficient = { value = 0.42, warrant = "fitted_here", observations = { covered_positions = 180600412 } } } ]

[ordinary_site_prior]
reference_concentration = 1.0
alternative_concentration_total = 0.0006
rung = "fitted_curve"

[repeat_tracts]
stated_length_spectrum_concentration = { value = 1.0, warrant = "defaulted" }
slippage_group_by_read_group = [ { read_group = 0, slippage_group = 0 } ]
slippage_by_stratum_and_group = [ { period = 2, reference_repeats = 6, slippage_group = 0, level = 0.04, shorter_share = 0.8, fall_off = 0.3, level_origin = { smoothing = { blend = { curve_weight = 0.37, reach = "inside", curve = { rise_shape = 0.5, intercept = 0.01, slope = 0.004, fitted_from_repeats = 5, fitted_to_repeats = 19, held_out_error = 0.2, cells = 23 } } }, slipped_reads = 8000.5 } } ]
length_spectrum_by_stratum = []
length_spectrum_by_period = []
substitution_rate_by_stratum = []

[stated_constants]
repeat_tract_outlier_weight = { value = 0.01, warrant = "defaulted" }
"#;
        let file: ParametersFile = toml::from_str(text).expect("the documented inline form parses");

        let row = &file.repeat_tracts.slippage_by_stratum_and_group[0];
        assert_eq!(row.level_origin.slipped_reads, Some(8000.5));
        assert!(
            matches!(
                row.level_origin.smoothing,
                LevelSmoothing::Blend {
                    curve_weight: 0.37,
                    reach: CurveReach::Inside,
                    ..
                }
            ),
            "the four-deep inline nesting reaches the curve weight and the reach"
        );
        // The omitted key is absence, not a default: this row has no shares origin at all, and
        // the whole document has no contamination section.
        assert!(row.shares_origin.is_none());
        assert!(
            file.contamination.is_none(),
            "an omitted contamination section reads as an uncontaminated run"
        );
    }

    /// **What format version this build writes.**
    ///
    /// Pinned so that bumping it is a deliberate act with a test to change, rather than a literal
    /// somebody edits in a fixture.
    #[test]
    fn the_format_version_this_build_writes_is_one() {
        assert_eq!(FORMAT_VERSION, 1);
        assert_eq!(a_file_using_every_shape().format_version, 1);
    }
}
