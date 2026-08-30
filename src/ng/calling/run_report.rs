//! **What the run scored its reads under, in a form an output can print** — the contamination
//! fraction each read group was corrected for, the batching those fractions were drawn against,
//! and the one repeat-tract constant nothing measured.
//!
//! # Why this exists at all
//!
//! A genotype computed at a contamination fraction of 3 in 100 and one computed at zero are
//! different claims about the same reads, and **nothing in the called genotype says which it
//! was** (`doc/devel/ng/spec/read_likelihoods.md` §3.6). The same holds for the one repeat-tract
//! number no fit produced; the two a *cell* can fall back to are per locus and travel on the
//! locus's own record. So the run carries what it used beside what it called, and this is the
//! type that carries it.
//!
//! Every number here is read off [`RunParameters`](super::run_parameters::RunParameters)
//! unchanged. **Nothing is computed and nothing is summarised**: a report that averaged two of a
//! sample's fractions would erase the one distinction the finer grain exists to express.
//!
//! # The grain: one row per read group, listed under its sample
//!
//! §3.6 asks for the fraction *per sample*; the fit produces one *per read group*, and a sample's
//! read groups can carry genuinely different fractions — a neighbouring library hopping its index
//! on the sequencer contaminates the run it is on and not the plant. **So a row is a read group,
//! and it names the sample it belongs to.** A per-sample line would have to pick one of the
//! fractions or average them, and both throw away the claim that they differ.
//!
//! **A read group is not always a library, and the row names both.** `@RG LB` is a grouping key: a
//! preparation sequenced over several lanes gives several read groups sharing one library name,
//! which is what [`ReadGroup::experiment`](crate::ng::read::input::read_groups::ReadGroup) exists
//! to say. Two lanes of one preparation can still differ in this number, because index hopping
//! happens on a flowcell rather than in a tube. So the row carries the read group's own `@RG ID`
//! beside its library's name, and a reader who sees one library name twice is looking at two lanes
//! rather than at a duplicated row.
//!
//! What the grain costs is nothing in the common case: a plant sequenced once, in one lane, gets
//! exactly one row — every sample of both benchmark cohorts here. What it buys is that a plant
//! whose reads came from two sequencing runs reports two fractions rather than one number that is
//! neither of them.

use super::ContaminationView;
use crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight;
use crate::ng::read::input::read_groups::NameWithOrigin;
use crate::ng::types::ReadGroupId;

/// **What one run used, as an output stage must be able to state it.**
///
/// **To be assembled** once per run, after the parameters are gathered and before any locus is
/// called — nothing here moves while the caller runs. The output stage that will print it is step
/// 10's and is not written, so today the only callers are tests.
///
/// **This is the run half. The per-locus half travels on the locus**: which rung of the tract
/// ladder a repeat tract's prior shape came from, and how many of its scoring cells fell back to
/// a stated constant, are properties of that tract and ride on
/// [`RepeatTractProvenance`](super::RepeatTractProvenance).
#[derive(Clone, PartialEq, Debug)]
pub struct RunParameterReport {
    contamination: ContaminationUsed,
    sequencing_batching: SequencingBatchingUsed,
    repeat_tract_outlier_weight: RepeatTractOutlierWeight,
}

impl RunParameterReport {
    /// Build the report from its three parts.
    ///
    /// **Plain data with no arithmetic**, so the only thing to get wrong is the pairing — three
    /// parts read off three different runs. That is why this is `pub(crate)`:
    /// [`RunParameters::report`](super::run_parameters::RunParameters::report) is the only
    /// builder that can exist, and it reads all three off one run. A consumer reads a report; it
    /// does not assemble one.
    #[must_use]
    pub(crate) fn new(
        contamination: ContaminationUsed,
        sequencing_batching: SequencingBatchingUsed,
        repeat_tract_outlier_weight: RepeatTractOutlierWeight,
    ) -> Self {
        Self {
            contamination,
            sequencing_batching,
            repeat_tract_outlier_weight,
        }
    }

    /// The contamination fraction every read group was corrected for, or that none was fitted.
    #[inline]
    #[must_use]
    pub fn contamination(&self) -> &ContaminationUsed {
        &self.contamination
    }

    /// Whether the batching those fractions were drawn against was declared or assumed.
    ///
    /// **It travels with the fractions rather than beside them**: two runs under different
    /// batchings produce contaminant frequencies that are not comparable, and the dense
    /// per-read-group view the caller holds cannot tell a declaration of one batch from a
    /// defaulted one — they are the same values
    /// (`doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
    #[inline]
    #[must_use]
    pub fn sequencing_batching(&self) -> SequencingBatchingUsed {
        self.sequencing_batching
    }

    /// **How often a read at a repeat tract is explained by none of the alleles the caller is
    /// considering** — a tract copied elsewhere in the genome, a chimera made during library
    /// preparation, a length only some of this tissue's cells carry, or a read anchored in the
    /// wrong place (`doc/devel/ng/spec/read_likelihoods.md` §4.5). **It is not contamination**,
    /// which is somebody else's DNA and has a term of its own; §4.5.1 keeps the two apart
    /// precisely because both are ways of a read not coming from this individual's two copies.
    ///
    /// **Nothing in the parameter fit measures it**, so it carries its warrant rather than a
    /// bare number: `Defaulted` where the run took the inherited 0.01, and `Supplied` where a
    /// parameters file gave it one (`doc/devel/ng/spec/parameters_file.md` §3.8, which puts it
    /// in that file so a person can change it). Those two are the only states it reaches.
    ///
    /// **Reported at the run rather than at the cell, and that is a ruling rather than a
    /// convenience.** It is one run-wide number, so folding it into a repeat tract's
    /// per-`(read group, candidate)` warrant would mark **every** tract of every run as resting
    /// on a defaulted parameter — or, where a file supplied it, on a supplied one, which the
    /// ladder ranks only a rung above — and erase the fitted-against-borrowed distinction that
    /// warrant exists to carry. Stating it once per run says the same true thing and costs
    /// nothing.
    #[inline]
    #[must_use]
    pub fn repeat_tract_outlier_weight(&self) -> RepeatTractOutlierWeight {
        self.repeat_tract_outlier_weight
    }
}

/// **What the run corrected for contamination**, per read group.
///
/// **The two arms are *nothing was estimable anywhere* and *a row for every read group*.** They
/// are not the same claim: at one sample there is no panel to compare against, so no fraction is
/// estimable at all, and that is [`Self::NoneFitted`] — not a run measured everywhere and found
/// clean. Within the second arm, [`ReadGroupContamination::was_measured`] is what tells a read
/// group that was measured from one nothing could be measured on
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
#[derive(Clone, PartialEq, Debug)]
pub enum ContaminationUsed {
    /// **The parameter fit identified no fraction anywhere**, so every locus of the run was
    /// scored on the plain formula with no mixture at all.
    ///
    /// *Absent, not a fitted zero.* At one sample it is the only possible answer, since
    /// contamination is a comparison between samples.
    NoneFitted,
    /// One row per read group of the run, in the run's sample order and, within a sample, in
    /// read-group order.
    ///
    /// **Every read group of a contaminated run has a row**, including the ones that identified
    /// nothing — those carry a zero fraction and zero evidence counts, which is what
    /// [`ReadGroupContamination::was_measured`] reads to tell them from a read group that was
    /// measured and found clean.
    PerReadGroup(Vec<ReadGroupContamination>),
}

impl ContaminationUsed {
    /// The rows, or an empty slice where nothing was fitted.
    ///
    /// **A convenience for a reader that only wants to list them**, and deliberately not the
    /// representation: an empty slice cannot say *nothing was estimable here*, which is the
    /// distinction the two arms exist for.
    #[must_use]
    pub fn rows(&self) -> &[ReadGroupContamination] {
        match self {
            Self::NoneFitted => &[],
            Self::PerReadGroup(rows) => rows,
        }
    }
}

/// **Whose reads a contamination fraction was applied to, and the fraction itself.**
///
/// The names are this type's own work; the number and the three things spec §3.6 says must
/// travel beside it are [`ContaminationView`]'s, **carried whole rather than copied out field by
/// field.** A row that spelled the fraction, the two evidence counts and their source again
/// would be a second copy of a value that already exists, and the two could then disagree about
/// what the run used.
#[derive(Clone, PartialEq, Debug)]
pub struct ReadGroupContamination {
    /// Which sample these reads belong to, as an index into **the run's sample order** — the
    /// same order [`LocusInference::per_sample`](super::LocusInference::per_sample) is in, so a
    /// reader can join a call to the fraction it was called under.
    pub sample: usize,
    /// The sample's name, as its `@RG SM` tag gives it.
    pub sample_name: Box<str>,
    /// The read group's place on the run's read-group axis — the identifier the pipeline mints,
    /// which is an index and not a name anybody wrote in a file.
    pub read_group: ReadGroupId,
    /// The read group's own `@RG ID`, verbatim — **the name that tells two lanes of one
    /// preparation apart**, since they share a library name and can carry different fractions.
    ///
    /// A label rather than an identity: the file format makes it unique within its own file and
    /// says nothing across files, so two files may each declare `ID:1`.
    pub read_group_name: Box<str>,
    /// The library's name **and whether the file declared it or this pipeline made it up**.
    ///
    /// `@RG LB` where the file carried one; otherwise a name synthesized from what the header did
    /// carry. **The origin travels with the name rather than being dropped**, because a report
    /// about the chemistry has to be able to say that a grouping is ours and not the file's —
    /// which is the reason [`NameWithOrigin`] is a pair at all. A synthesized name reported as a
    /// declared one is a claim about the run that nobody made.
    ///
    /// **Several read groups can share it**: the lanes of one preparation are one library.
    pub library: NameWithOrigin,
    /// **What the read likelihood used for these reads**, exactly as the run's parameters hold
    /// it: the fraction, how many of the panel's varying positions this read group put a read on,
    /// how many reads it put there, and whose reads the fraction was fitted from.
    ///
    /// The last of those matters because a fraction fitted from this read group's own reads and
    /// one fitted from every read of the plant and copied onto it are different claims — only the
    /// first can say that two of a plant's read groups differ.
    ///
    /// **Read [`Self::was_measured`] before reporting the source**: a read group the fit could
    /// not measure carries a source value that is true of nothing, because no variant of
    /// [`ContaminationSource`](crate::ng::parameter_estimation::joint::contamination::ContaminationSource)
    /// says *nothing was fitted here*.
    pub estimate: ContaminationView,
}

impl ReadGroupContamination {
    /// Whether anything at all stood behind this fraction.
    ///
    /// **A fraction near zero because nothing could be measured is not the claim a fraction
    /// measured and found clean is**, and the counts are the only thing that tells them apart:
    /// the search keeps zero where the likelihood barely moves with the fraction, which is the
    /// right default for a value the mixture multiplies and the wrong thing to read as *clean*
    /// (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
    ///
    /// **The predicate is [`ContaminationView::was_measured`]'s and is not restated here**, for
    /// the reason the field above gives: a rule spelled twice is a rule that can disagree with
    /// itself.
    #[must_use]
    pub fn was_measured(&self) -> bool {
        self.estimate.was_measured()
    }
}

/// **Whether the run was told who was sequenced beside whom, or nobody said.**
///
/// A contaminating read is likelier to have come from a neighbour on the same flowcell than from
/// a random member of the species, so the frequency the contaminant's allele is looked up in is
/// one batch's rather than the whole cohort's. Under the default that batch *is* the whole
/// cohort, which is the honest statement of what a run knows when nobody has said otherwise —
/// and a fraction drawn against it is the weaker kind of number.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SequencingBatchingUsed {
    /// **Nobody declared one**, so one batch held the run and every library was scored against
    /// the whole cohort's frequency. This is what both benchmark cohorts here get, and today it
    /// is what every run gets: no command-line flag carries a batching yet.
    DefaultedToOneBatch,
    /// The run declared this many batches.
    Declared {
        /// How many the declaration named. At least one, since a batching naming none is
        /// refused where it is built.
        batches: usize,
    },
}
