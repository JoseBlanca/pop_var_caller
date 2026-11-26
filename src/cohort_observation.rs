//! What the cohort merge hands the calling steps: one assembled cohort locus
//! ([`CohortObservation`]) and the evidence it holds, plus the keep rule ([`MinAltReads`]) that
//! both the merge and candidate selection ask.
//!
//! **The types live here and their constructors do not.** The merge builds a
//! [`CohortObservation`] from a closed locus in `run/cohort_merge/build.rs`
//! ([`CohortObservation::over`]), and the rows inside it are filled there too; this module holds
//! only the definitions, so that the calling library can read its input without depending on the
//! driver that produces it (`doc/devel/ng/spec/cohort_merge.md` §4.2, §4.3).

use std::num::NonZeroU32;

use crate::locus_generation::{FlankDeletion, HeldFlank, LocusKind, WitnessedLocusPositions};
use crate::types::{GenomeRegion, ReadGroupId};
use crate::window_coverage::WindowCoverage;

/// Where the reference sits in a locus's allele table: **first**, whether or not any
/// sample's reads showed it, because a cohort in which every sample is homozygous for the
/// variant still has to be genotyped against it.
pub const REFERENCE_ALLELE: usize = 0;

/// One cohort locus, assembled: the ground, the alleles the cohort showed over it, and
/// what each covering sample's reads lend each allele (arch §4).
///
/// **This is evidence, not a call.** Nothing here says which alleles are real, what they are
/// worth, or what genotype a sample has; those are the calling steps' (spec §1.2).
///
/// **Built in `run/cohort_merge/build.rs`**, by [`over`](Self::over) and
/// [`over_from_cursors`](Self::over_from_cursors), from a locus the merge closed.
#[derive(Debug)]
pub struct CohortObservation {
    /// The locus span, first position to furthest reach.
    pub region: GenomeRegion,
    /// The distinct alleles, the reference at [`REFERENCE_ALLELE`] (`doc/devel/ng/spec/cohort_merge.md` §4.2).
    pub alleles: Vec<Box<[u8]>>,
    /// The covering samples, in ascending sample order — **only** the covering ones.
    ///
    /// **A sample with no coverage over the span has no entry**, which is a different fact
    /// from an entry whose support is all reference and stays one (`doc/devel/ng/spec/cohort_merge.md` §4.2). The
    /// architecture's sketch says "indexed by the run's sample order", which reads two ways;
    /// this is the reading that keeps the distinction structural rather than resting on a
    /// zeroed row, and it is the shape the walk hands over
    /// ([`SampleMembers`](crate::run::cohort_merge::close::SampleMembers)). Each entry names its own sample.
    pub per_sample: Vec<SampleSupport>,
    /// What kind of stretch this is — a SNP/indel site ([`LocusKind::Generic`]) or a
    /// repeat tract ([`LocusKind::Ssr`]), which carries the motif the tract's read model
    /// needs and both reference flanks.
    ///
    /// **The kind was minted by the generator that produced every member of this locus**
    /// and travelled here unchanged; the walk asserts the members agree on it
    /// (`close.rs`). Without it the motif — and so the period, and so the repeat counts —
    /// was dropped at assembly, and a caller handed a cohort observation could not tell a
    /// tract from a stretch of ordinary sequence
    /// (`doc/devel/ng/spec/run_ssr_observations.md` §4).
    ///
    /// **What the clone costs, per built locus.** [`LocusKind::Generic`] owns nothing, so
    /// the SNP/indel path — the overwhelming majority of built loci — copies one
    /// discriminant. [`LocusKind::Ssr`] copies an inline motif — `MAX_MOTIF_LEN` is 6, so
    /// the bytes travel in the enum rather than on the heap — and allocates both flanks,
    /// at most 15 bytes each at the tract generator's default flank width and shorter
    /// where a contig end clamps one. Nothing is cloned for the loci the walk closes and
    /// does not build.
    pub kind: LocusKind,
    /// Each covering sample's read depth and GC over the 500-base window centred on this
    /// locus's **first base** — one entry per entry of [`per_sample`](Self::per_sample), at the
    /// same index (`doc/devel/ng/spec/window_coverage.md` §3.5).
    ///
    /// **Both fields `NaN` where the sample has no usable window there**, which is one answer to
    /// three questions and deliberately so: the sample's coverage begins a base later, or the
    /// window it would have had was built from too few covered positions and the floor silenced
    /// it, or the run took no measurement at all. The filter that reads this skips a sample it
    /// has no window for, which is what its scorer already does for an absent sample, so the
    /// three need no distinguishing here — but **compared by bit pattern, never by `==`**, which
    /// [`WindowCoverage`]'s own `PartialEq` does.
    ///
    /// **The third of those is why the drivers' agreement comparison excludes this field.** A
    /// merge run over records held whole in memory takes no measurement, so its every entry is
    /// absent; comparing it against the cache's would be comparing absence against a number.
    /// `cohort_merge::fixtures::render` says so at the one place that decides what "the same
    /// answer" means for the cohort merge.
    pub window_coverage: Vec<WindowCoverage>,
    /// **Each allele's flank deletion, at the same index as [`alleles`](Self::alleles)** —
    /// the flank bases a deletion running out of an STR segment took (owner, 2026-10-10;
    /// `str_edge_deletion_alleles.md` §5.4). **Empty means none for every allele**, which is
    /// every generic locus and every hand-built fixture; read it through
    /// [`flank_deletion_of`](Self::flank_deletion_of).
    pub allele_flank_deletions: Vec<FlankDeletion>,
}

impl CohortObservation {
    /// The flank deletion of allele `allele`: none where the table carries none.
    pub fn flank_deletion_of(&self, allele: usize) -> FlankDeletion {
        self.allele_flank_deletions
            .get(allele)
            .copied()
            .unwrap_or(FlankDeletion::NONE)
    }
}

/// One sample's evidence at one cohort locus, against that locus's allele table.
///
/// **Built in `run/cohort_merge/build.rs`**, with the allele table it indexes into, while
/// [`CohortObservation::over`] assembles the locus.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleSupport {
    /// Which sample — its index in the run's sample order, as the walk named it.
    pub sample: usize,
    /// **What this sample's reads showed, one row per `(allele, read group)`**, in ascending
    /// `(allele, read group)` order. A pair it showed no reads for has no entry, which
    /// `pooled_support_for` reads as no reads and no sums — the
    /// same answer a zeroed row would give, without a row per sample per allele per group.
    ///
    /// **The read-group axis all but spent what the sparse shape was saving, and that is worth
    /// stating rather than leaving the old sentence standing.** A row was 40 bytes and is now 48
    /// — the support is 32, the allele index 8, the group 4 and the padding 4 — where a dense
    /// record would hold a bare 32-byte support per cell. So a sample showing 2 of a locus's 3
    /// alleles from one library costs 96 bytes either way, and one showing 2 of 2 costs 96
    /// against a dense 64. **What still pays is the cohort-wide case the shape was chosen for**:
    /// a dense record is samples × alleles × *groups* cells whatever anyone showed, and it is the
    /// third factor that makes it unaffordable — a panel where one sample in eleven carries a
    /// second library (`doc/devel/ng/spec/read_groups.md` §1 — 157 of 1,707 samples in a
    /// surveyed tomato archive carry more than one) would pay for every sample the widest
    /// sample's groups.
    ///
    /// **Support is never merged across alleles**, because a genotype likelihood needs them
    /// apart (`doc/devel/ng/spec/cohort_merge.md` §4.2), **and never across read groups**, because the likelihood pools an
    /// observation's reads into one term only if every one of them would get the same number
    /// — and two reads of the same bases from two lanes have different error rates
    /// (`doc/devel/ng/spec/read_likelihoods.md` §2.3). It *is* merged within one pair: where two of the
    /// sample's own observations reached the same allele from the same group — two records'
    /// worth of one read, say — their reads and their sums are added.
    ///
    /// **A sample with one read group has exactly today's shape**, which is most samples of
    /// most runs, so the axis costs them nothing.
    pub supported: Vec<SupportedAllele>,
    /// **What this sample's reads showed of only *part* of the locus** — one entry per
    /// `(record, sequence, read group)` whose reads ran out inside it.
    ///
    /// **Ascending `(witnessed stretch, read group, bases)`, and all three parts of that key are
    /// needed.** Two entries can share a stretch — a substitution witnessed partially is exactly
    /// that — so a stretch alone does not order them, and `doc/devel/ng/spec/read_likelihoods.md`
    /// §8 makes the order a determinism requirement rather than a tidiness one: the sum over
    /// observations must run in a fixed order. The sibling field states its own key the same way
    /// and `tests::the_rows_are_ordered_by_allele_then_read_group` pins it; this one is pinned
    /// by `tests::the_partial_rows_are_ordered_by_stretch_then_read_group_then_bases`, which
    /// varies all three components so that neither dropping one nor reordering them survives.
    ///
    /// **Kept apart from [`supported`](Self::supported) rather than folded into it**, because a
    /// partial read does not say what the sample carries, it says the sample carries *at least*
    /// this, and the two claims cannot share a row: padding a partial's bases out to the locus
    /// span and interning them would put an allele in the table that no molecule carried, and it
    /// would read as a *short* allele, which is the one direction the model must not be biased
    /// in (`doc/devel/ng/spec/read_likelihoods.md` §5.1).
    ///
    /// **Filled by `build.rs`'s `partials_of_sample`, which walks the sample's records directly** rather
    /// than through [`MemberPlacement::projectable_sequences`](crate::run::cohort_merge::build::MemberPlacement::projectable_sequences) — the gate that keeps partials
    /// out of the allele derivation is exactly the one this field exists to get past. Empty
    /// where the sample's reads all spanned their records, which is the ordinary locus; the
    /// censored term of `read_likelihoods.md` §5 reads it (§5.4, corrected 2026-08-21).
    ///
    /// **It changes no locus's existence, and did not.** Whether a locus is built at
    /// all is decided from complete observations only — the filter is
    /// [`SampleLocusObservations::non_reference_and_compared_reads`](crate::locus_generation::SampleLocusObservations::non_reference_and_compared_reads), which skips every
    /// non-`Complete` observation before counting, and it lives in the locus generator rather
    /// than here.
    ///
    /// **That rule is settled for the generic path and explicitly *not* settled for repeat
    /// tracts.** §5.4.2 answers *no* under the heading "the merge's rule stays as it is **on the
    /// generic path**", and then says the opposite for tracts: a sample carrying an allele too
    /// long for a read to span shows no complete observation at all, so the filter reads it as
    /// *nothing varied here*, and "one line of the rule has to change for that to mean what it
    /// says, and only on this path". This row is on both paths — a locus of either kind that
    /// passes the verdict is assembled here — so the tract half is owed to whoever brings the
    /// STR path through the merge, and is out of this plan's scope.
    pub partials: Vec<PartialObservation>,
    /// Reads that covered one of this sample's records here and produced no observation at
    /// all — carried through from the members, not re-derived (arch §4).
    ///
    /// **Summed over the sample's records, which overstates it where a locus spans several
    /// of them**: a read silent at two of them is counted twice, since the mint records how
    /// many said nothing and not which ones. Exact where the sample has one record, which is
    /// the ordinary locus. Saturating.
    pub reads_without_observation: u32,
    /// Reads of this sample that were removed as evidence — named at some of its records
    /// inside the locus and not at all of them, or named more than once at one of them, so
    /// nothing they showed reaches `supported` (see `build.rs`'s `alleles_of_sample`). Lost depth,
    /// counted rather than inferred from an absence. Saturating.
    pub reads_removed_as_evidence: u32,
    /// Reads whose allele was composed across several of this sample's records, and whose
    /// five quality sums in `supported` are therefore **divided rather than measured** (see
    /// [`AlleleSupport`]). Zero on a one-record sample, where every sum is exact. Saturating.
    pub reads_composed_across_records: u32,
}

impl SampleSupport {
    #[cfg(test)]
    /// What this sample's reads lend `allele` **added up over its read groups** — nothing,
    /// where it showed that allele no reads, which is the answer a zeroed row would have
    /// given.
    ///
    /// **A read likelihood must not use this.** Pooling an allele's reads across groups is
    /// exactly what `doc/devel/ng/spec/read_likelihoods.md` §2.3 forbids: the reads of one term have to
    /// share an error rate, and two lanes do not. Iterate [`supported`](Self::supported)
    /// instead, which is why the pooling is in this method's name. What it is for is the
    /// questions that really are about the sample rather than about one library — how much
    /// depth the allele has here, and the fixtures that ask it.
    ///
    /// A scan, not an index, because [`supported`](Self::supported) holds only the pairs the
    /// sample showed: at an ordinary locus that is one or two entries, and it is bounded by
    /// how many alleles the whole cohort showed times how many groups this sample has.
    pub fn pooled_support_for(&self, allele: usize) -> AlleleSupport {
        self.supported
            .iter()
            .filter(|supported| supported.allele == allele)
            .fold(AlleleSupport::default(), |total, supported| {
                total.added_to(supported.support)
            })
    }
}

/// One allele of the locus **as one of the sample's read groups showed it**, and what those
/// reads lend it.
///
/// **The read group is part of the row's identity, not a label on it.** A sample whose reads
/// for one allele came from two lanes has two rows here, and nothing downstream may add them:
/// a read likelihood folds an observation's reads into one term only when every one of them
/// would get the same number, and the two lanes' error rates differ
/// (`doc/devel/ng/spec/read_likelihoods.md` §2.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportedAllele {
    /// Which allele — an index into [`CohortObservation::alleles`].
    pub allele: usize,
    /// Which of the sample's read groups showed it — one `@RG`, i.e. one lane. Carried from
    /// [`SequenceObservation::read_group`](crate::locus_generation::SequenceObservation::read_group), which keeps the groups apart at the mint for the
    /// same reason.
    pub read_group: ReadGroupId,
    /// What this sample's reads **from that group** lend it. Never zero: a pair a sample
    /// showed no reads for has no entry at all.
    pub support: AlleleSupport,
}

/// **What one sample's reads showed of part of the locus, and nothing outside that part.**
///
/// A read that entered the locus and ran off its own end shows a prefix or a suffix of what the
/// sample carries — or, since the generic fold began minting witnesses with holes, a stretch
/// with a gap in the middle. It does not say what the sample carries; it says the sample carries
/// **at least** this, over the positions it saw
/// (`doc/devel/ng/spec/read_likelihoods.md` §5.1, whose "prefix or suffix" predates the hole). Scoring it is the calling step's
/// — a censored term, which is **never less likely** than a complete observation of the same
/// bases and so can never read as evidence for a short allele (§5.2 at a repeat tract, §5.3 at
/// an ordinary locus). Carrying it is the cohort merge's, and until now the merge threw it away.
///
/// **Not "less discriminating", which this sentence said until 2026-08-25 and §5.2 no longer
/// claims.** A lower bound can separate two candidates *further* than a complete read does,
/// where one candidate is shorter than the stretch the read saw and the other is longer; §5.2's
/// correction box carries the sizes.
///
/// **Where these exist at all is decided by how wide the locus is on the reference**, and it is
/// not intuition: a single-base substitution has none, an insertion has none — its reference
/// span is its anchor base however long the inserted sequence — and a deletion has them, carrying
/// the reference allele preferentially, because a read carrying the deletion crosses every
/// deleted position without spending a read base. The class this is really for is the repeat
/// tract, where over half the overlapping reads are partial at a 60-base tract and **an allele
/// longer than a read can only ever be witnessed partially** (§5.4.1).
///
/// **One entry is one `(record, sequence, read group)`, not one allele**, and that is the point:
/// there is no allele. A partial's bases cannot be compared against a whole-span allele — that
/// comparison would report a read agreeing with the reference over everything it saw as
/// non-reference — so what a candidate is scored against is the allele's projection *restricted
/// to the positions the read witnessed* (§5.3), and both halves of that restriction are here.
///
/// **Named as the architecture names it** (`doc/devel/ng/arch/read_likelihoods.md` §2.1,
/// `GenericSampleEvidence::partials`), and owned where that sketch had a borrow: the calling
/// view can hold `&[PartialObservation]`, so there is one type rather than two.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialObservation {
    /// **Which locus positions the read witnessed, in the cohort locus's own coordinates** —
    /// zero at the locus's first position, whatever record inside it the observation came from.
    /// The axis is in the name because it is not in the type: this is the same
    /// [`WitnessedLocusPositions`] the mint puts on [`ReadWitness::Partial`](crate::locus_generation::ReadWitness::Partial), where the runs are
    /// measured against the *record*. A newtype would close that off and is the right move when
    /// a second consumer appears; with one writer and no reader it would be a type minted for
    /// nobody.
    ///
    /// The mint measures a witness against the *record* it belongs to
    /// (`SampleLocusObservations::locus_len`), and a cohort locus can hold several of one
    /// sample's records and start before any of them, so the runs are shifted onto the locus
    /// when the row is built. **A consumer indexing an allele's projection with an unshifted run
    /// would read the wrong bases and nothing would say so.**
    ///
    /// A set of runs rather than one offset and one length: the generic fold mints witnesses
    /// with holes in them, and two numbers can only describe a hole by swallowing it
    /// ([`ReadWitness::Partial`](crate::locus_generation::ReadWitness::Partial)).
    pub witnessed_in_locus: WitnessedLocusPositions,
    /// Which of the sample's read groups these reads came from — the same axis
    /// [`SupportedAllele`] keys on, and for the same reason: a censored term pools an
    /// observation's reads only if every one of them would get the same number
    /// (`doc/devel/ng/spec/read_likelihoods.md` §2.3).
    pub read_group: ReadGroupId,
    /// The bases these reads showed, exactly as the mint recorded them — **allele content, in
    /// read coordinates**, which is not the axis [`witnessed_in_locus`](Self::witnessed_in_locus) is on.
    ///
    /// **The two lengths differ by the net indel the read carried over the stretch, so equality
    /// does not license indexing one with the other.** A read carrying a two-base insertion and
    /// a two-base deletion inside the witnessed stretch comes back with as many bases as
    /// positions and is not a positional match for any of them.
    pub bases: Box<[u8]>,
    /// How many of the sample's reads showed exactly this — the same stretch, the same bases,
    /// the same read group, which is what makes them one term.
    ///
    /// **Never zero, kept so by whoever builds the row rather than by this type** — the fields
    /// are public and there is no constructor, exactly as [`SupportedAllele`]'s are, where the
    /// same rule is executed as a filter at the moment the rows are made.
    pub num_reads: u32,
    /// Σ `ln P(error)` over those reads — the mint's own sum, not divided, because a partial is
    /// not composed across records and so has nothing to apportion.
    ///
    /// **This and `num_reads` are the whole of what a partial carries, where a complete row
    /// carries six numbers**, and the plan asks for exactly these two. What the other four are
    /// for is worth knowing before anything consumes this: strand bias, the mapping-quality
    /// multi-mapper test and the read-position-bias term all read them, so at a repeat tract —
    /// where half the overlapping reads are partial at 50 reference bases and four in five at
    /// 100 (spec §5.4.1) — those filters would see the complete reads only. Adding them later
    /// means reopening whatever routes partials in.
    ///
    /// **That, and the ascending order the field on [`SampleSupport`] promises, are the builder's
    /// conventions rather than the type's**: the fields are public and there is no constructor,
    /// so `build.rs`'s `partials_of_sample` is where both hold, and where their tests point.
    pub q_sum: f64,
    /// **Which border, if any, the reads' sequence may run on past into a longer allele**, though
    /// their witness covers the whole locus. `None` for every partial the merge mints: a read that
    /// ran out inside the locus says so through its witness, and that is the only shape the walk
    /// produces.
    ///
    /// `Some` only for reads the calling run re-read after the merge
    /// ([`reread_spellings_cut_short`](crate::run::records::reread_spellings_cut_short)): reads
    /// the mapper laid straight across the reference near their end, so that the start of a longer
    /// allele came out as substitutions. Over the locus they show exactly the start of that allele
    /// followed by the reference past the locus (`Right`), or its end preceded by the reference
    /// before it (`Left`), so they are compared as a read flush to that border is — compatible
    /// with every allele whose carrier would show the same bases there
    /// (`doc/devel/ng/spec/read_likelihoods.md` §5.3).
    pub sequence_may_run_on_past: Option<LocusBorder>,
    /// **Which flank of an STR locus the reads anchored, where their witness cannot say it** —
    /// the mint's [`SequenceObservation::held_flank`](crate::locus_generation::SequenceObservation::held_flank), carried through unchanged: a partial at an
    /// STR segment or STR bundle whose witness covers the whole locus, and `None` everywhere else
    /// (`doc/devel/ng/spec/str_bundle_loci.md` §7.6).
    pub held_flank: Option<HeldFlank>,
}

/// One of a locus's two borders — see
/// [`PartialObservation::sequence_may_run_on_past`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LocusBorder {
    /// The locus's first position.
    Left,
    /// The locus's last position.
    Right,
}

/// What one sample's reads lend one allele. **Built in `run/cohort_merge/build.rs`**, which
/// also divides the sums described below.
///
/// **The read count is exact and the five sums are not, wherever a locus spans several of
/// the sample's records.** Every read is named, so it lands on exactly one allele and
/// `num_reads` counts it there. The five sums are stored by the mint already summed over the
/// reads behind one observed sequence, so when those reads take different paths across the
/// locus the sums have to be divided between the alleles — nothing recorded says how much of
/// a sum belongs to which read. [`SampleSupport::reads_composed_across_records`] says how
/// many reads at this locus were treated that way; where it is zero, every sum here is the
/// mint's own.
///
/// **How they are divided follows production's rule, and in one place not production's
/// code** — its merger faces the same question (`project_compound_scalars`,
/// `var_calling/per_group_merger.rs`). This is the owner's ruling of 2026-08-17, taken after
/// checking freebayes, which never faces it because it holds one record per read all the way
/// to the likelihood (`freebayes/src/Sample.h`, a vector of per-read observations per
/// allele):
///
/// - **the quality sum** takes, for each read, the *weakest* of the mean per-read qualities
///   among the sequences that read showed across the locus, and adds those up over the
///   reads — an allele spanning several records is no better evidenced than its weakest
///   piece, since it is wrong if any piece is. These are sums of `ln P(error)`, so the
///   weakest read is the **largest** number, and this is the place the code parts company
///   with production's: see `build.rs`'s `share_of_one_read`;
/// - **the strand and mapping-quality sums** take each read's share as the pooled mean over
///   the sequences it showed: the sums over those sequences divided by their read counts.
///   Both are properties of the read and identical at each of its sightings, so the pooling
///   estimates one number rather than mixing two;
/// - **the placement count is pooled the same way, and it is the one that mixes two
///   questions.** It is counted against each record's own position (see `placed_left`), so
///   two records of one sample ask "does the read lie mostly left of here?" about different
///   *heres*, and the pooled answer is an approximation whose disagreement is bounded by the
///   locus's width — 100 bases at the default bound.
///
/// **The strand and placement counts are left as fractions; the MAPQ sums are rounded** once
/// per row rather than per read. Rounding the strand share made the strand test depend on which
/// strand a read was *not* on: a lone composed read with a forward share of ½ rounded to one
/// forward read on either strand, and swapping every read's strand moved the strand-and-position
/// penalty at 4,260 tomato records (read-transform oracle, 2026-11-01). Both tests read the
/// counts as rates, so a fraction costs them nothing. **Per row means per `(allele, read
/// group)` since B1**: two lanes' MAPQ shares of one divided read are rounded separately, so
/// adding the rows back can differ by up to half a unit per lane from what a single row would
/// have held. `tests::a_divided_read_keeps_its_strand_share_per_read_group` pins the fractions
/// at one read divided over two lanes. (Named rather than linked: a `#[cfg(test)]` item does not
/// exist in a doc build, so a link to one is an unresolved link.)
///
/// **This row is one allele of one read group, not one allele of the sample** — see
/// [`SupportedAllele`], which carries the group. It was one row per allele until the calling
/// prerequisites split it (`doc/devel/ng/impl_plan/calling_prerequisites.md`, Milestone B step
/// B1), and the reason for splitting turned out to be general rather
/// than the STR path's alone: a read likelihood folds an observation's reads into one term only
/// when every one of them would get the same number, and reads from two lanes have different
/// error rates (`doc/devel/ng/spec/read_likelihoods.md` §2.3). Stutter, fitted per read group, is the same
/// argument on the repeat path — a locus called from a pooled row would be scored against a
/// rate belonging to no group in particular. The mint says the same thing from its own side:
/// "a per-chemistry model needs the allele × group cross **with its quality moments**"
/// ([`SequenceObservation::read_group`](crate::locus_generation::SequenceObservation)), a
/// per-group count beside one merged row giving the first and losing the second.
///
/// A sample with one read group has one row per allele, exactly as before.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AlleleSupport {
    /// How many of the sample's reads showed this allele. **Exact** — every read is named,
    /// so it is counted where it belongs rather than divided. Saturating, which needs four
    /// billion reads on one allele of one sample at one locus.
    pub num_reads: u32,
    /// Of those, how many were on the forward strand — strand bias. **A fraction, not rounded**:
    /// a read composed across records counts as the pooled forward share of the sequences it
    /// showed, and rounding that per row made a lone read with a share of ½ count forward
    /// whichever strand it was on (read-transform oracle, 2026-11-01). The strand test reads it
    /// as a rate, so nothing downstream needs it whole.
    pub num_fwd: f64,
    /// Σ per-read log-error over them — the freebayes per-read error term.
    pub q_sum: f64,
    /// Σ MAPQ over them.
    pub mapq_sum: u32,
    /// Σ MAPQ² over them; with `mapq_sum` and `num_reads` it recovers the mean and variance
    /// the MAPQ multi-mapper filter reads.
    pub mapq_sum_sq: u64,
    /// How many lay mostly left of **the record they were seen at** — freebayes'
    /// `placedLeft`, and the read-position-bias term.
    ///
    /// **Not counted against the cohort locus's own first base.** The mint counts it against
    /// each record's position — the record's anchor on the generic path
    /// (`locus_generation/pileup/`), the whole tract on the STR path
    /// (`locus_generation/ssr.rs`) — and nothing here re-anchors it, so where a locus spans
    /// several of a sample's records this mixes as many questions as it has records.
    ///
    /// A fraction for the reason [`num_fwd`](Self::num_fwd) is.
    pub placed_left: f64,
    /// How many had a low base quality inside the record they were seen at
    /// ([`SequenceObservation::num_low_quality`](crate::locus_generation::SequenceObservation::num_low_quality)) — the low-base-quality penalty's count.
    pub num_low_quality: u32,
}

impl AlleleSupport {
    #[cfg(test)]
    /// Two rows added together — **the read-group fold, and the only place it happens**.
    ///
    /// Every field is additive, so this is exact wherever the rows themselves are.
    /// [`SampleSupport::pooled_support_for`] is its one caller, and its doc says why a read
    /// likelihood must not be a second one.
    fn added_to(self, other: Self) -> Self {
        Self {
            num_reads: self.num_reads.saturating_add(other.num_reads),
            num_fwd: self.num_fwd + other.num_fwd,
            q_sum: self.q_sum + other.q_sum,
            mapq_sum: self.mapq_sum.saturating_add(other.mapq_sum),
            mapq_sum_sq: self.mapq_sum_sq.saturating_add(other.mapq_sum_sq),
            placed_left: self.placed_left + other.placed_left,
            num_low_quality: self.num_low_quality.saturating_add(other.num_low_quality),
        }
    }
}

/// The fewest non-reference **reads** one sample may show at a cohort locus and still
/// have it built — the floor half of [`MinAltReads`] (spec §4.3).
///
/// The `Obs` is production's word for a read's allele observation, and **not** the cohort
/// merge's `observation`, which is one sample's whole record over a stretch of genome
/// (spec §1.3). One observation can carry many of these reads.
///
/// Below the threshold in *every* sample the locus is **dropped**: nothing is assembled,
/// nothing is emitted, and — unlike a locus that failed [`MaxCohortLocusSpan`](crate::run::cohort_merge::MaxCohortLocusSpan) — nothing
/// is counted. A failure is ground the caller refused; this is ground it judged empty,
/// and conflating the two would stop the failed count meaning anything.
///
/// **A command-line parameter**, default [`DEFAULT_MIN_ALT_OBS`], and the reason it
/// exists is measured rather than aesthetic: in production, under production's rule
/// (below), it removes a large number of very low-quality variants and improves
/// performance substantially. Its cost is that a variant seen once, in one sample, is
/// unrecoverable — nothing downstream is emitted for that locus at all.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MinAltObs(pub NonZeroU32);

impl MinAltObs {
    /// The default threshold, [`DEFAULT_MIN_ALT_OBS`] non-reference reads.
    pub const DEFAULT: Self = Self(non_zero_default(DEFAULT_MIN_ALT_OBS));

    /// The threshold in non-reference reads.
    #[inline]
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

impl Default for MinAltObs {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 2 reads — production's number, over a rule that is now production's shape.
///
/// **ng asks the same question production asks, and asks it of a share as well.**
/// Production takes, across a group's positions, the *maximum over samples* of that
/// sample's non-reference observations (`max_nonref_obs`,
/// `cohort_integration.rs:64-78`, summed in `derive_is_kept`) — a per-sample test. ng
/// summed every sample's non-reference reads instead until 2026-08-19, and the owner
/// retired that rule when it was measured: because the bar was a fixed count and the sum
/// grows with the cohort, it stopped filtering as samples were added. Over one 100 kb
/// tomato interval it discarded 997 loci in 1,000 at one accession and only 439 in 1,000
/// at 63, where the merge then assembled a cohort observation at 546 positions per
/// kilobase.
///
/// The value is copied and the name is not. Production spelled it
/// `var_calling::DEFAULT_MIN_ALT_OBS_PER_SAMPLE`, the default of
/// `--min-alt-obs-per-sample`, and fed that one number both to the cohort keep rule ng's
/// descends from (`derive_is_kept`) and to a per-sample pre-EM filter
/// (`passes_min_alt_obs`, `variant_caller.rs`). ng's threshold is only the first, so
/// reaching the constant by name would have let a retune of production's per-sample filter
/// move ng's cohort keep. Named rather than linked: the deliberate non-reuse is the point,
/// and production is being deleted.
pub const DEFAULT_MIN_ALT_OBS: u32 = 2;

/// What share of one sample's compared reads must be non-reference — the share half of
/// [`MinAltReads`] (spec §4.3).
///
/// A fraction of 1, not a percentage: 0.02 is two reads in a hundred.
///
/// **The denominator is the sample's own reads, and that is the whole design.** A share
/// taken of the *cohort's* reads would grow its bar as samples are added while a rare
/// allele's evidence stays where it is — one carrier's handful of reads — so a large
/// study could not call the rare variation it was run to find (spec §4.3, and the
/// derivation there).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct MinAltReadShare(f64);

impl MinAltReadShare {
    /// The default share, [`DEFAULT_MIN_ALT_READ_SHARE`].
    pub const DEFAULT: Self = Self(DEFAULT_MIN_ALT_READ_SHARE);

    /// Whether `share` is a fraction of one — at least nothing, at most everything, and
    /// a real number.
    ///
    /// **One spelling for both constructors**, so that the fallible one and the `const`
    /// one cannot come to disagree about what a legal share is. That is not a
    /// hypothetical: a first draft of the pair wrote the test out twice, and a review's
    /// mutation of one copy left the other's tests green.
    ///
    /// The comparison is written out rather than asked of `(0.0..=1.0).contains(&share)`
    /// because `RangeInclusive::contains` is not a `const fn`. It refuses everything the
    /// range plus an `is_finite()` call would refuse: **a negative share and `-∞` fail
    /// the lower comparison, a share above one and `+∞` fail the upper, and `NaN` fails
    /// both.** `-0.0` is a share of nothing and passes, as it does through the range.
    const fn is_a_fraction_of_one(share: f64) -> bool {
        share >= 0.0 && share <= 1.0
    }

    /// The share, or `None` if it is not a fraction of 1 — negative, above one, or not a
    /// number. Rejecting rather than clamping: a mistyped share is a run that quietly
    /// keeps everything or nothing, and neither is visible in the output.
    pub fn new(share: f64) -> Option<Self> {
        Self::is_a_fraction_of_one(share).then_some(Self(share))
    }

    /// The same share, for a `const` that has to name one — and it **panics** where
    /// [`new`](Self::new) returns `None`.
    ///
    /// **The severity is the point, not a lapse.** It is written for `const`
    /// declarations of a rule the whole run is judged against, where a value outside
    /// `0..=1` fails the compile rather than the run and there is nobody to hand an
    /// `Option` to. Written for
    /// [`DEFAULT_MIN_ALLELE_SUPPORT`](crate::calling::allele_candidates::DEFAULT_MIN_ALLELE_SUPPORT),
    /// which needs a share of 5 in 100 and cannot reach this type's private field.
    ///
    /// **Called at runtime it aborts the process** — `[profile.release]` sets
    /// `panic = "abort"` — and nothing in the signature stops that, because a `const fn`
    /// is an ordinary function too. **A share that is not a literal in this source tree —
    /// one an operator typed, one read from a file, one computed — goes through
    /// [`new`](Self::new)**, which refuses rather than aborts.
    pub const fn new_or_panic(share: f64) -> Self {
        assert!(
            Self::is_a_fraction_of_one(share),
            "a minimum non-reference read share is a fraction of one"
        );
        Self(share)
    }

    /// The share, as a fraction of 1.
    #[inline]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl Default for MinAltReadShare {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 2 reads in a hundred — the owner's number, chosen against two measurements a
/// hundredfold apart in depth and not yet against a truth set.
///
/// **It is the half of the rule that only high depth can exercise.** At the tomato
/// panel's 11 reads a sample at a locus, two in a hundred asks for 0.22 reads and
/// [`MinAltObs`]'s floor of 2 decides everything — the share changed 12,034 kept loci to
/// 12,033. At the GIAB trio's 313 reads a sample it asks for 7, and there it does all
/// the work: 53,796 kept loci become 2,721 over 136 kb. A true heterozygote at that
/// depth shows about 156 non-reference reads, so the bar sits a factor of 20 below what
/// a real carrier produces.
pub const DEFAULT_MIN_ALT_READ_SHARE: f64 = 0.02;

/// **How much non-reference evidence one sample must show before the cohort builds a
/// locus** — the merge's keep rule, and the two numbers it is made of (spec §4.3).
///
/// A locus is built when **some single sample** shows at least
/// [`required_of`](Self::required_of) its compared reads as non-reference. Not the sum
/// across the cohort: a real variant concentrates its evidence in the samples that carry
/// it, so asking one sample is what tracks a carrier, and it is the only form of the
/// question whose answer does not drift as samples are added.
///
/// The two numbers answer the two ends of the depth axis, which is why neither alone
/// would do:
///
/// - the **floor** ([`MinAltObs`], 2 reads) is what decides at low coverage, where a
///   share of three reads rounds to nothing;
/// - the **share** ([`MinAltReadShare`], 2 in 100) is what decides at high coverage,
///   where two reads out of 300 is the error rate rather than an allele.
///
/// At one sample the rule is exactly what the cohort-sum rule it replaced was, because a
/// sum over one sample is that sample's own count.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct MinAltReads {
    /// The fewest non-reference reads that can ever suffice, whatever the depth.
    pub floor: MinAltObs,
    /// The share of a sample's compared reads asked for once depth makes it the larger.
    pub share: MinAltReadShare,
}

impl MinAltReads {
    /// The default rule: [`MinAltObs::DEFAULT`] reads or [`MinAltReadShare::DEFAULT`] of
    /// the sample's own, whichever is more.
    pub const DEFAULT: Self = Self {
        floor: MinAltObs::DEFAULT,
        share: MinAltReadShare::DEFAULT,
    };

    /// How many non-reference reads a sample that compared `reads_compared_with_reference`
    /// against the reference must show — the floor, or the share of those reads rounded
    /// **up**, whichever is larger.
    ///
    /// Rounded up because rounding down would let the share ask for less than it says: at
    /// 200 reads and a share of 2 in 100, 4 reads is the bar and 3 is not.
    ///
    /// The cast cannot go wrong at either end. The product is finite and never negative,
    /// so a share of 1 at the largest possible read count lands exactly on `u32::MAX`;
    /// anything past that saturates rather than wrapping, which is Rust's rule for
    /// float-to-integer casts since 1.45.
    #[inline]
    pub fn required_of(self, reads_compared_with_reference: u32) -> u32 {
        let asked = (self.share.get() * f64::from(reads_compared_with_reference)).ceil();
        self.floor.get().max(asked as u32)
    }

    /// Whether one sample that showed `non_reference_reads` out of
    /// `reads_compared_with_reference` reaches the rule — the question the merge asks of
    /// every contributing sample at every locus.
    ///
    /// **The numerator is the caller's, and the callers count different reads.** The
    /// merge passes a sample's non-reference reads pooled over the whole allele table,
    /// which is what the argument is named for. Candidate selection passes *one
    /// allele's* reads, so that the question "is this sequence worth calling over"
    /// is the merge's own question asked one level down
    /// (`doc/devel/ng/spec/candidate_alleles.md` §3). A discovery round passes a
    /// narrower count still — the reads the converged model is explaining as slippage
    /// (§3.4). The denominator is the same in all three: that sample's compared reads
    /// at the locus. **There is deliberately no wrapper per caller**: a second spelling
    /// of one rule is how two rules become different rules.
    ///
    /// **The floor is tested first because it is the cheap half and it cannot be
    /// skipped.** [`required_of`](Self::required_of) is never below the floor, so a
    /// sample under the floor fails whatever its depth — and answering it by integer
    /// comparison keeps the multiply off the path taken at almost every position, where
    /// the sample showed no non-reference read at all.
    #[inline]
    pub fn reached_by(self, non_reference_reads: u32, reads_compared_with_reference: u32) -> bool {
        non_reference_reads >= self.floor.get()
            && non_reference_reads >= self.required_of(reads_compared_with_reference)
    }
}

/// A default's value as the [`NonZeroU32`] the newtypes hold.
///
/// **A zero default is a build error, not a panic**, and it is the *call* that makes it
/// one: every caller is a `pub const DEFAULT: Self` item, which is evaluated when the
/// crate is compiled, so the `None` arm below never survives into a running binary.
/// Anything added later gets the same guarantee from the same line — there is no
/// separate assertion to remember to write.
///
/// Writing the defaults as `u32` rather than as `NonZeroU32` keeps them readable at
/// their declaration, which is where an operator reading the source looks for them.
pub(crate) const fn non_zero_default(default_value: u32) -> NonZeroU32 {
    match NonZeroU32::new(default_value) {
        Some(value) => value,
        None => panic!("a cohort-merge default must be non-zero"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The two halves of the keep rule, each at the depth where it decides.** The
    /// crossing point is where the share first exceeds the floor, and it is the number
    /// this rule's behaviour hinges on: below it every sample is asked for the floor, and
    /// above it the bar climbs with the sample's own depth.
    #[test]
    fn the_floor_decides_at_low_depth_and_the_share_at_high() {
        let rule = MinAltReads::DEFAULT;

        // At and below 100 reads, 2% asks for no more than the floor of 2.
        assert_eq!(
            rule.required_of(0),
            2,
            "no reads at all still asks the floor"
        );
        assert_eq!(rule.required_of(3), 2, "2% of 3 rounds up to 1");
        assert_eq!(rule.required_of(100), 2, "2% of 100 is exactly the floor");

        // Above it the share takes over, and it rounds up rather than down: 2% of 150 is
        // 3 exactly, and 2% of 151 is 3.02, which asks for 4.
        assert_eq!(
            rule.required_of(101),
            3,
            "the first depth the share decides"
        );
        assert_eq!(rule.required_of(150), 3);
        assert_eq!(rule.required_of(151), 4, "rounded up, not down");
        assert_eq!(rule.required_of(313), 7, "the GIAB trio's depth a sample");
    }

    /// `reached_by` is `required_of` with the comparison folded in — pinned separately
    /// because it carries a short-circuit on the floor, and a short-circuit that answered
    /// differently from the arithmetic it is skipping would be invisible in the merge's
    /// output.
    #[test]
    fn reaching_the_rule_agrees_with_what_it_asks_for() {
        let rule = MinAltReads::DEFAULT;

        for compared in [0u32, 1, 3, 50, 100, 101, 150, 313, 1_000] {
            for alt in 0..=12u32 {
                assert_eq!(
                    rule.reached_by(alt, compared),
                    alt >= rule.required_of(compared),
                    "{alt} non-reference of {compared} compared"
                );
            }
        }
    }

    /// A share is a fraction of one, and anything else is refused rather than clamped —
    /// a run that silently kept everything or nothing would look like a run that found
    /// everything or nothing.
    #[test]
    fn a_share_outside_a_fraction_of_one_is_refused() {
        assert_eq!(
            MinAltReadShare::new(0.0).map(MinAltReadShare::get),
            Some(0.0)
        );
        assert_eq!(
            MinAltReadShare::new(1.0).map(MinAltReadShare::get),
            Some(1.0)
        );
        assert_eq!(
            MinAltReadShare::new(0.02).map(MinAltReadShare::get),
            Some(0.02)
        );

        assert!(MinAltReadShare::new(-0.01).is_none(), "below zero");
        assert!(MinAltReadShare::new(1.01).is_none(), "above one");
        assert!(MinAltReadShare::new(f64::NAN).is_none(), "not a number");
        assert!(MinAltReadShare::new(f64::INFINITY).is_none(), "not finite");
        assert!(
            MinAltReadShare::new(f64::NEG_INFINITY).is_none(),
            "not finite, the other way"
        );
    }

    /// **The `const` constructor refuses exactly what [`MinAltReadShare::new`] refuses**,
    /// which is the whole claim `is_a_fraction_of_one` exists to make true — one test,
    /// both constructors, over the same ten values.
    ///
    /// `new_or_panic` panics where `new` returns `None`, so the two are compared through
    /// [`std::panic::catch_unwind`] rather than by reading them side by side. The values
    /// that matter here are the ones no earlier fixture reached: **a negative share**,
    /// which does not crash anything downstream but silently deletes the share half of
    /// the rule — `required_of` casts a negative product to `u32`, which saturates to
    /// zero, so `max(floor, 0)` is the floor at every depth — and **both ends of the
    /// closed range**, which must be accepted rather than refused.
    #[test]
    fn the_const_share_refuses_exactly_what_the_fallible_one_refuses() {
        for share in [
            -1.0,
            -0.05,
            -0.0,
            0.0,
            0.05,
            1.0,
            1.5,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ] {
            let accepted_by_new = MinAltReadShare::new(share).is_some();
            let accepted_by_new_or_panic =
                std::panic::catch_unwind(|| MinAltReadShare::new_or_panic(share)).is_ok();
            assert_eq!(
                accepted_by_new, accepted_by_new_or_panic,
                "the two constructors disagree about a share of {share}"
            );
        }
    }

    /// The `const` constructor's message, so the two ends of its range are pinned as
    /// *refusals* and not merely as "something went wrong".
    ///
    /// **A negative share is the one that has to be loud.** Above one and `NaN` are
    /// visibly nonsense; a negative share reads as a plausible typo and, unrefused,
    /// produces a run that keeps every allele the floor admits at every depth with
    /// nothing in the output saying so.
    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_below_zero_fails_rather_than_clamping() {
        let _ = MinAltReadShare::new_or_panic(-0.05);
    }

    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_above_one_fails_rather_than_clamping() {
        let _ = MinAltReadShare::new_or_panic(1.5);
    }

    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_that_is_not_a_number_fails() {
        let _ = MinAltReadShare::new_or_panic(f64::NAN);
    }

    /// Both ends of the closed range are legal shares — no read at all, and every read.
    #[test]
    fn a_const_share_may_sit_on_either_end_of_the_range() {
        assert_eq!(MinAltReadShare::new_or_panic(0.0).get(), 0.0);
        assert_eq!(MinAltReadShare::new_or_panic(1.0).get(), 1.0);
    }

    /// **A share of one at the largest read count lands on `u32::MAX` and does not
    /// wrap.** The cast in `required_of` is the only arithmetic in this rule that could
    /// answer a number smaller than the truth, and a keep rule that asked for *less* at
    /// the extreme would build every locus in a run that had asked for none.
    #[test]
    fn the_share_at_its_extreme_saturates_rather_than_wrapping() {
        let everything = MinAltReads {
            floor: MinAltObs::DEFAULT,
            share: MinAltReadShare::new(1.0).expect("a fraction of one"),
        };

        assert_eq!(everything.required_of(u32::MAX), u32::MAX);
        assert!(!everything.reached_by(u32::MAX - 1, u32::MAX));
        assert!(everything.reached_by(u32::MAX, u32::MAX));
    }
}
