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
//! - **`RunParameters` cannot say *absent*; the file must.** The five states of spec §5 — an
//!   absent contamination table, a read group measured and found clean, a defaulted calibration
//!   multiplier of exactly one, a stratum with no length spectrum, a `(stratum, slippage group)`
//!   with no row — are distinctions a reader that collapses them will get wrong, and the file's
//!   shape is where they are made unmistakable.
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
//! by_read_group = [ { read_group = 0, error_probability_multiplier = 1.0324,
//!                     warrant = "fitted_here" } ]
//!
//! [contamination]                  # §3.4 — the whole table absent means uncontaminated
//! by_read_group = [ { read_group = 0, library = "lib3", fraction = 0.031,
//!                     markers_with_reads = 4211, reads_on_markers = 90233,
//!                     fitted_from_reads_of = "this_read_groups_own_reads" } ]
//!
//! [sequencing_batches]             # §3.4 — declared by the run, not fitted
//! was_declared_by_the_run = false
//! by_read_group = [ { read_group = 0, batch = 0 } ]
//! by_sample = [ { sample = "TS-1", batch = 0 } ]
//!
//! [inbreeding]                     # §3.5 — the file's only cohort-sized axis
//! by_sample = [ { sample = "TS-1", inbreeding_coefficient = 0.42 } ]
//!
//! [ordinary_site_prior]            # §3.6 — the seed itself, never the moments behind it
//! reference_concentration = 1.0
//! alternative_concentration_total = 0.0006
//! rung = "fitted_curve"
//!
//! [repeat_tracts]                  # §3.7
//! stated_length_spectrum_concentration = 1.0
//! stated_length_spectrum_warrant = "defaulted"
//! slippage_group_by_read_group = [ { read_group = 0, slippage_group = 0 } ]
//! slippage_by_stratum_and_group = [ … one row a (stratum × slippage group) … ]
//! length_spectrum_by_stratum = [ … only where the stratum was fitted on its own tracts … ]
//! length_spectrum_by_period = [ … only where the run asked for the middle rung … ]
//! substitution_rate_by_stratum = [ … ]
//!
//! [stated_constants]               # §3.8 — the numbers no fit produces
//! repeat_tract_outlier_weight = 0.01
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
    pub contamination: Contamination,
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
    pub error_probability_multiplier: f64,
    /// Where the multiplier came from.
    pub warrant: Warrant,
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
    /// One row a read group, in dense-index order. Where **some** read group identified a
    /// fraction, **every** read group needs an entry.
    pub by_read_group: Vec<ContaminationRow>,
}

/// One read group's contamination fraction, and the evidence behind it.
///
/// **Only the counts tell *measured and found clean* from *not measured*** — a fraction near
/// zero is what both produce, and a read group that touched no marker was not measured whatever
/// its fraction says (spec §5).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContaminationRow {
    /// The run's dense read-group index.
    pub read_group: u32,
    /// The library this read group's reads were prepared from — several read groups of one
    /// preparation share it, and it is written here because the whole point of the read-group
    /// grain is that two lanes of one library may differ.
    pub library: String,
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
    /// One row a sample, in the run's sample order. **At least one is required.**
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
    /// The coefficient itself, a fraction in `[0, 1)`.
    pub inbreeding_coefficient: f64,
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
    pub stated_length_spectrum_concentration: f64,
    /// Whether that number is the run's own median over the strata it fitted, or the stated
    /// constant reached where it fitted none.
    ///
    /// **The two are the same line without it.** The run's median can land on exactly 1.0, which
    /// is also `STATED_FLAT_CONCENTRATION`, and spec §8 counts this among the three numbers with
    /// an honest default that must be marked `defaulted` when it is one. It is the same reason
    /// the calibration multiplier carries a warrant beside a value of one.
    pub stated_length_spectrum_warrant: Warrant,
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
    /// The rate itself, a probability.
    pub rate: f64,
    /// Where the rate came from.
    pub warrant: Warrant,
    /// How many reads stood behind it. **The one evidence count in the file that does not name
    /// its unit in its key**, kept for now because step A2 gives every value+warrant+count row
    /// one spelling and that is where the word is settled once rather than twice.
    pub observations: u64,
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
    pub repeat_tract_outlier_weight: f64,
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
                        error_probability_multiplier: 1.0324,
                        warrant: Warrant::FittedHere,
                    },
                    // A multiplier of exactly one that was *not* fitted — the state the warrant
                    // exists to keep apart from a fitted one.
                    BaseQualityCalibrationRow {
                        read_group: 1,
                        error_probability_multiplier: 1.0,
                        warrant: Warrant::Defaulted,
                    },
                    BaseQualityCalibrationRow {
                        read_group: 2,
                        error_probability_multiplier: 0.87,
                        warrant: Warrant::Supplied,
                    },
                ],
            },
            contamination: Contamination {
                by_read_group: vec![
                    ContaminationRow {
                        read_group: 0,
                        library: "lib3".into(),
                        fraction: 0.031,
                        markers_with_reads: 4211,
                        reads_on_markers: 90233,
                        fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                    },
                    // Measured and found clean: a fraction of zero with evidence behind it.
                    ContaminationRow {
                        read_group: 1,
                        library: "lib4".into(),
                        fraction: 0.0,
                        markers_with_reads: 3117,
                        reads_on_markers: 71004,
                        fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                    },
                    ContaminationRow {
                        read_group: 2,
                        library: "lib5".into(),
                        fraction: 0.0072,
                        markers_with_reads: 2903,
                        reads_on_markers: 64118,
                        fitted_from_reads_of: ContaminationFittedFrom::EveryReadOfThisSample,
                    },
                ],
            },
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
                        inbreeding_coefficient: 0.42,
                    },
                    InbreedingRow {
                        sample: AWKWARD_SAMPLE.into(),
                        inbreeding_coefficient: 0.17,
                    },
                ],
            },
            ordinary_site_prior: OrdinarySitePrior {
                reference_concentration: 1.0,
                alternative_concentration_total: 0.000_6,
                rung: SeedRung::FittedCurve,
            },
            repeat_tracts: RepeatTracts {
                stated_length_spectrum_concentration: 1.25,
                stated_length_spectrum_warrant: Warrant::Defaulted,
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
                    rate: 0.0012,
                    warrant: Warrant::Borrowed,
                    observations: 40_122,
                }],
            },
            stated_constants: StatedConstants {
                repeat_tract_outlier_weight: 0.01,
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
        empty.contamination.by_read_group.clear();
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

    /// **A stated concentration of one says whether it was fitted.**
    ///
    /// The run's own median over the strata it fitted can land on exactly 1.0, which is also the
    /// constant reached where it fitted none. Without the warrant beside it the two are the same
    /// line in the file — the collapse spec §5 forbids one section above, for the calibration
    /// multiplier, and spec §8 names this number among the three that must be marked
    /// `defaulted` when they are.
    #[test]
    fn a_stated_concentration_of_one_says_whether_it_was_fitted() {
        let mut fitted = a_file_using_every_shape();
        fitted.repeat_tracts.stated_length_spectrum_concentration = 1.0;
        fitted.repeat_tracts.stated_length_spectrum_warrant = Warrant::FittedHere;

        let mut defaulted = fitted.clone();
        defaulted.repeat_tracts.stated_length_spectrum_warrant = Warrant::Defaulted;

        assert_ne!(
            toml::to_string(&fitted).expect("serialises"),
            toml::to_string(&defaulted).expect("serialises"),
            "a fitted 1.0 and a stated 1.0 must not be the same file"
        );
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

        let extra = format!("{text}\n[stated_constants]\nrepeat_tract_outlier_wieght = 0.02\n");
        assert!(
            toml::from_str::<ParametersFile>(&extra).is_err(),
            "a misspelled key must not be absorbed"
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
by_read_group = [ { read_group = 0, error_probability_multiplier = 1.0, warrant = "defaulted" } ]

[contamination]
by_read_group = []

[sequencing_batches]
was_declared_by_the_run = false
by_read_group = [ { read_group = 0, batch = 0 } ]
by_sample = [ { sample = "TS-1", batch = 0 } ]

[inbreeding]
by_sample = [ { sample = "TS-1", inbreeding_coefficient = 0.42 } ]

[ordinary_site_prior]
reference_concentration = 1.0
alternative_concentration_total = 0.0006
rung = "fitted_curve"

[repeat_tracts]
stated_length_spectrum_concentration = 1.0
stated_length_spectrum_warrant = "defaulted"
slippage_group_by_read_group = [ { read_group = 0, slippage_group = 0 } ]
slippage_by_stratum_and_group = [ { period = 2, reference_repeats = 6, slippage_group = 0, level = 0.04, shorter_share = 0.8, fall_off = 0.3, level_origin = { smoothing = { blend = { curve_weight = 0.37, reach = "inside", curve = { rise_shape = 0.5, intercept = 0.01, slope = 0.004, fitted_from_repeats = 5, fitted_to_repeats = 19, held_out_error = 0.2, cells = 23 } } }, slipped_reads = 8000.5 } } ]
length_spectrum_by_stratum = []
length_spectrum_by_period = []
substitution_rate_by_stratum = []

[stated_constants]
repeat_tract_outlier_weight = 0.01
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
        // The omitted key is absence, not a default: this row has no shares origin at all.
        assert!(row.shares_origin.is_none());
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
