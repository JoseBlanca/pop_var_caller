//! Assembling a cohort locus: widening every member's observed sequence to the whole
//! locus span.
//!
//! A cohort locus is one stretch of genome, and the samples in it did not all record
//! the same stretch: one sample's deletion covers five bases where another's SNP covers
//! one. Before the two can be compared they have to be written over the same ground —
//! each sequence padded with the reference bases on either side until it spans the whole
//! locus (`doc/devel/ng/spec/cohort_merge.md` §4.2, projection), and then unified into one
//! table of distinct alleles, the reference among them (§4.2, unification —
//! [`AlleleTable`]). That is what this file does; expressing each sample's *support*
//! against that table is the next step's (`doc/devel/ng/impl_plan/cohort_merge.md` B3).
//!
//! **What a sample showed is derived per read**, because a cohort locus can span several of
//! that sample's records, and a read named at some of them but not others cannot be placed
//! across it (the owner's ruling of 2026-08-17, and B0, which names every read the mint
//! folds). [`alleles_of_sample`] carries the rule, and is where the two things it does
//! **not** claim are stated: the ground outside a sample's own records is reference by
//! construction rather than by evidence, and a read is judged on its presence at those
//! records rather than on coverage of the whole locus.
//!
//! **Where the padding bases come from is a small departure from the spec's sentence.**
//! §4.2 says they "travel on the observation already", meaning
//! [`SampleLocusObservations::reference_bases`]. They do, but only over that
//! observation's own region — a SNP at one position carries one reference base and
//! knows nothing of the four the deletion beside it covers. So the locus's reference is
//! **gathered across its members** first ([`LocusReferenceBases::over`]), and each member
//! is padded from that. The observations are still the only source: nothing is fetched,
//! and no reference file is read.
//!
//! Production does the same widening from a reference it holds directly
//! (`project_local_allele`, `var_calling/per_group_merger.rs` — prefix, allele, suffix);
//! the arithmetic here is that function's, with the reference assembled rather than
//! fetched.

use ahash::AHashMap;

use super::close::{ClosedLocus, LocusCloser, SampleMembers, Verdict, span_of};
use super::{MaxCohortLocusSpan, MinAltObs};
use crate::ng::locus_generation::{ReadWitness, SampleLocusObservations, SequenceObservation};
use crate::ng::types::GenomeRegion;
use crate::pileup_record::ChainId;

/// The byte a gathered reference position starts as, so that a position no member
/// covered is a loud failure rather than a `NUL` inside every allele.
///
/// **Zero is not a base, and `N` would not do.** ng's reference fetch uppercases ACGT and
/// folds everything else to `N` (`ng/ref_seq.rs`, through `fasta::fetcher::canonicalise`),
/// so no position of a real reference can hold zero — while `N` is what every assembly gap
/// gathers, and a sentinel spelled `b'N'` would refuse loci the members do cover. A byte
/// still holding this one after every member has been copied in means the members did not
/// cover the locus they were handed with.
const NOT_COVERED: u8 = 0;

/// The reference's own bases over one whole cohort locus — what every member's
/// projection is padded from.
///
/// Gathered from the locus's members, each of which carries the reference over its own
/// region ([`SampleLocusObservations::reference_bases`]). **That they cover the locus
/// between them is what closing guarantees**, not something checked upstream: a locus
/// opens at some observation's first base and every observation that joins it starts at
/// or before the reach it had at the time, so the members' regions form one unbroken
/// stretch from the locus's first base to its last (spec §4.1). The gather asserts it
/// anyway — see [`NOT_COVERED`].
#[derive(Debug)]
pub struct LocusReferenceBases {
    /// The ground these bases cover — the closed locus's own region.
    region: GenomeRegion,
    /// One base per position of `region`, in order.
    bases: Box<[u8]>,
}

impl LocusReferenceBases {
    /// Gather the reference over `locus` from the members' own reference bases.
    ///
    /// **Only a locus the caller undertakes to build.** Closing is uncapped — how wide a
    /// locus may be is a verdict passed on it afterwards, not a limit on closing it
    /// (spec §4.1) — so a `Failed` locus can be thousands of bases wide and is never
    /// assembled (spec §3.2). The verdict is asserted here rather than left to the
    /// caller, which is also what bounds this function's allocation: a built locus is at
    /// most `max_cohort_locus_span` bases for a generic locus and the catalog's tract
    /// width for an STR one (spec §3.1).
    ///
    /// **Most of the checks here are against a caller's mistake, not against the data**,
    /// so they are assertions rather than errors: the members and the region come out of
    /// one walk ([`super::close::LocusCloser`]), which cannot produce a member outside the
    /// locus it closed. What they defend against is a later caller pairing a locus with
    /// somebody else's members — where the failure is silent, since a partly gathered
    /// reference pads every allele in the locus with `NUL` bytes and the alleles still
    /// unify, still count, and still look like sequences.
    ///
    /// **The reference-width check is the exception, and it is deliberately here anyway.**
    /// It asserts that one observation is self-consistent, which is its *producer's*
    /// guarantee rather than the pairing's — the walk never reads `reference_bases`.
    /// Today the only producer is this crate's own generator, so a mismatch is a bug here
    /// and a panic is right. **When observations are decoded from a psp file it becomes a
    /// corrupt-input failure and must become a `RunError`** beside
    /// `ObservationExceedsReachCeiling` (arch §5): it is the same class of fact about the
    /// same file.
    ///
    /// Release-level, deliberately: the release profile is the one this repo runs, and the
    /// cost is one pass over the gathered span, which the verdict above bounds.
    pub fn over(locus: &ClosedLocus<'_>) -> Self {
        assert_eq!(
            locus.verdict,
            Verdict::Build,
            "the locus {} was judged {:?}, and only a locus the caller undertakes to \
             build is assembled",
            locus.region,
            locus.verdict,
        );
        // PANIC-FREE: `usize` is 64 bits on every target this crate builds for, so the
        // conversion is total. It is not a claim about the region — nothing bounds a
        // `GenomeRegion`'s span, which is why `offset_within` below checks rather than
        // assumes.
        let span = usize::try_from(span_of(locus.region))
            .expect("a u64 span fits a usize on a 64-bit target");
        let mut bases = vec![NOT_COVERED; span];

        for sample_members in &locus.members {
            let mut previous: Option<GenomeRegion> = None;
            for member in sample_members.observations {
                assert_eq!(
                    member.region.contig, locus.region.contig,
                    "member {} of sample index {} is on another contig from the locus {}",
                    member.region, sample_members.sample, locus.region,
                );
                // **One sample's own records are disjoint and ascending, and that is
                // checked here rather than discovered while composing.** The generic mint
                // cannot produce two overlapping records — an event falling inside an open
                // record's footprint folds into that record and widens it rather than
                // opening a second one (`find_overlapping`,
                // `locus_generation/pileup/open_record.rs`) — and the walk hands each
                // sample's members over in coordinate order (`close.rs`, `SampleMembers`).
                //
                // Checked anyway, because the failure without it is the quiet one: a read's
                // allele across the locus is its sequences written side by side, so two
                // records covering the same base have no one composition, and every read
                // that showed them is dropped as evidence — a sample contributing nothing
                // at a locus where it had two changes, with nothing to say why. Composition
                // has a backstop of its own ([`MemberPlacement::compose_into`]), but it is
                // reached only when one read is named at both records, so the same defect
                // would be loud or silent depending on which reads the data carried.
                //
                // Two of a sample's records being disjoint is the **producer's** guarantee,
                // like the reference width below, so **when observations are decoded from a
                // psp file this must become a `RunError`** beside
                // `ObservationExceedsReachCeiling` (arch §5) rather than a panic.
                if let Some(previous) = previous {
                    assert!(
                        previous.end < member.region.start,
                        "the records {} and {} of sample index {} are not disjoint and in \
                         coordinate order, so no read's allele across the locus {} can be \
                         composed",
                        previous,
                        member.region,
                        sample_members.sample,
                        locus.region,
                    );
                }
                previous = Some(member.region);

                assert_eq!(
                    member.reference_bases.len() as u64,
                    span_of(member.region),
                    "member {} of sample index {} carries {} reference bases for a \
                     {}-base region",
                    member.region,
                    sample_members.sample,
                    member.reference_bases.len(),
                    span_of(member.region),
                );
                assert!(
                    member.reach() <= locus.region.end,
                    "member {} of sample index {} reaches past the locus {} it belongs to",
                    member.region,
                    sample_members.sample,
                    locus.region,
                );

                let offset = offset_within(locus.region, member.region);
                // **Written base by base so that two members overlapping cannot disagree
                // silently.** Members overlap by construction — a locus is a chain of
                // overlapping observations — and where two samples' reference bases differ
                // over shared ground the samples were called against different references.
                // Copying the second over the first would leave a reference that is a
                // mixture, and every allele in the locus would be padded from it: the same
                // plausible-looking wrong answer the other checks exist to stop. The cost
                // is one comparison per base of a span the verdict above bounds.
                for (index, (slot, base)) in bases[offset..offset + member.reference_bases.len()]
                    .iter_mut()
                    .zip(member.reference_bases.iter().copied())
                    .enumerate()
                {
                    assert!(
                        *slot == NOT_COVERED || *slot == base,
                        "two members of the locus {} disagree on the reference at \
                         position {}: {:?} and {:?} — the samples were called against \
                         different references",
                        locus.region,
                        locus
                            .region
                            .start
                            .get()
                            .saturating_add((offset + index) as u64),
                        *slot as char,
                        base as char,
                    );
                    *slot = base;
                }
            }
        }

        if let Some(gap) = bases.iter().position(|base| *base == NOT_COVERED) {
            panic!(
                "the members of the locus {} leave position {} uncovered (offset {gap} of \
                 {}), so its reference cannot be gathered",
                locus.region,
                locus.region.start.get().saturating_add(gap as u64),
                bases.len(),
            );
        }

        Self {
            region: locus.region,
            bases: bases.into_boxed_slice(),
        }
    }

    /// The ground these bases cover.
    pub fn region(&self) -> GenomeRegion {
        self.region
    }

    /// The reference's bases over the whole locus, one per position.
    ///
    /// **This is also the reference allele**: it is what a sample whose reads all matched
    /// the reference projects to, so the table of alleles the next step builds has it
    /// among them without a special case (spec §4.2). Pinned by
    /// `a_member_that_matched_the_reference_projects_to_the_locus_reference`.
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Work out where `member` sits inside the locus, once, so that every sequence that
    /// member showed is projected at the same offset over the same width.
    ///
    /// **The handle is what ties a sequence to the member it came from.** The two are only
    /// meaningful together — the member supplies the offset and the width of reference the
    /// sequence replaces, the sequence supplies the bases — and padding one member's
    /// sequence at another's offset yields a well-formed byte string that the next step
    /// would accept as an allele. Taking both as loose arguments made that expressible;
    /// going through the member first does not.
    pub fn placing<'a>(&'a self, member: &'a SampleLocusObservations) -> MemberPlacement<'a> {
        assert_eq!(
            member.region.contig, self.region.contig,
            "the member {} is on another contig from the locus {}",
            member.region, self.region,
        );
        assert!(
            member.reach() <= self.region.end,
            "the member {} reaches past the locus {} it is being projected onto",
            member.region,
            self.region,
        );

        let offset = offset_within(self.region, member.region);
        // PANIC-FREE: the assertion above puts the member's span inside the locus's, and
        // the locus's span already converted when the reference was gathered.
        let covered = usize::try_from(span_of(member.region))
            .expect("a member's span is inside the locus's, which fits a usize");

        MemberPlacement {
            reference: self,
            member,
            offset,
            covered,
        }
    }
}

/// Where one member sits inside its locus — worked out once by
/// [`LocusReferenceBases::placing`], then reused for every sequence that member showed.
#[derive(Debug)]
pub struct MemberPlacement<'a> {
    reference: &'a LocusReferenceBases,
    member: &'a SampleLocusObservations,
    /// How far into the locus the member's own region starts, in bases.
    offset: usize,
    /// How many reference bases the member's region covers — how much of the locus's
    /// reference the observed sequence replaces.
    covered: usize,
}

impl<'a> MemberPlacement<'a> {
    /// The sequences of this member that can be projected — the
    /// [`Complete`](ReadWitness::Complete) ones.
    ///
    /// **Reaching a member's sequences through here is what keeps a partial out**, and
    /// the reason is what a partial's bases are: they stop where its read's witness
    /// stopped, so padding them from the locus's reference would report that the read
    /// showed reference bases over ground it never saw — an allele no molecule carried.
    /// It is the same subset `non_reference_reads` counts over
    /// ([`SampleLocusObservations::complete_observations`]), and the cost is the same: a
    /// variant witnessed only by partial reads contributes nothing until the censored
    /// likelihood that would score one exists (`spec/locus_generation.md` §3).
    pub fn projectable_sequences(&self) -> impl Iterator<Item = &'a SequenceObservation> {
        self.member.complete_observations()
    }

    /// Widen one of this member's observed sequences to the whole locus span, into
    /// `projected`.
    ///
    /// The result is the reference before the member's own region, then the sequence the
    /// reads showed, then the reference after it — production's `project_local_allele`
    /// (`var_calling/per_group_merger.rs`). **Its length is the locus span only when the
    /// sequence matches its own region's width**: a deletion projects shorter than the
    /// span and an insertion longer, which is the point — the projected bytes are the
    /// allele, and how many bases of reference it replaces is what the member's region
    /// said.
    ///
    /// `projected` is cleared first and refilled, so one buffer serves a whole locus.
    ///
    /// A [`Partial`](ReadWitness::Partial) sequence panics rather than being padded, for
    /// the reason [`projectable_sequences`](Self::projectable_sequences) gives — that
    /// iterator is the way in that cannot reach one.
    pub fn project_into(&self, sequence: &SequenceObservation, projected: &mut Vec<u8>) {
        let bases = &self.reference.bases;
        projected.clear();
        projected.reserve(
            bases
                .len()
                .saturating_sub(self.covered)
                .saturating_add(sequence.bases.len()),
        );
        let written_to = self.compose_into(sequence, 0, projected);
        projected.extend_from_slice(&bases[written_to..]);
    }

    /// Write the reference from `written_to` up to this member's own region, then
    /// `sequence`'s bases in place of that region, and answer how far the writing has
    /// reached — an **offset into the locus's reference**, not into `composed`, since the
    /// two differ by every base an indel added or removed.
    ///
    /// **This is one substitution of a whole allele, and the caller decides how many go
    /// into one sequence.** [`project_into`](Self::project_into) makes one and closes with
    /// the reference to the end of the locus; deriving a read's allele across a locus that
    /// spans several of one sample's records makes one per record, in coordinate order,
    /// and closes the same way ([`AlleleTable`]).
    ///
    /// **A member starting before `written_to` panics rather than overwriting what is
    /// already composed**, since the read's two sequences would disagree about the bases
    /// they share and no one composition would be right; silently taking the later one
    /// would put a sequence in the table that no read showed.
    ///
    /// **It is a backstop, not the guard.** That a sample's records are disjoint and
    /// ascending is checked structurally, once per sample, before any read is consulted
    /// ([`LocusReferenceBases::over`]) — because this one is reached only when a single
    /// read is named at both of the overlapping records, and where they carry different
    /// reads the sample's evidence would vanish quietly instead. The message here says what
    /// was observed rather than diagnosing it, for the same reason: records supplied out of
    /// coordinate order trip it too.
    fn compose_into(
        &self,
        sequence: &SequenceObservation,
        written_to: usize,
        composed: &mut Vec<u8>,
    ) -> usize {
        assert_eq!(
            sequence.read_witness,
            ReadWitness::Complete,
            "only a complete observation can be projected onto the locus {}; the partial \
             at {} saw less than its own region",
            self.reference.region,
            self.member.region,
        );
        assert!(
            written_to <= self.offset,
            "the member {} of the locus {} starts at offset {} inside ground already \
             composed to {written_to}, so the sample's records are not disjoint and in \
             coordinate order",
            self.member.region,
            self.reference.region,
            self.offset,
        );

        let bases = &self.reference.bases;
        composed.extend_from_slice(&bases[written_to..self.offset]);
        composed.extend_from_slice(&sequence.bases);
        self.offset + self.covered
    }
}

/// Where the reference sits in a locus's allele table: **first**, whether or not any
/// sample's reads showed it, because a cohort in which every sample is homozygous for the
/// variant still has to be genotyped against it.
pub const REFERENCE_ALLELE: usize = 0;

/// The distinct alleles the samples showed over one cohort locus, the reference among
/// them — one table, against which every sample's support is expressed (spec §4.2).
///
/// **Two sequences that come out as the same bytes are the same allele, wherever they came
/// from.** That is the whole of unification, and it is what makes a cohort observation the
/// cohort's rather than a bundle of per-sample pieces: one sample's SNP inside another
/// sample's deletion becomes one entry once both are written over the locus's whole ground.
///
/// **What a sample showed is derived per read, not per record** — the owner's ruling of
/// 2026-08-17, and the reason B0 gives every read a chain id:
///
/// > Either we know the read covered the whole locus, and its allele is elongated with what
/// > it showed; or we know it did not cover it, and it is removed as evidence. Not being
/// > able to decide which is an error that must never happen.
///
/// [`alleles_of_sample`] is where that is turned into something decidable — presence at
/// each of *that sample's* records — and where the two branches, and their one difference,
/// are stated. Nothing here is a claim about ground the sample minted no record over.
///
/// **What this does not do is choose.** Which of these alleles are worth calling, how they
/// are written out and what they are worth are the calling steps' (spec §1.2, §13); this is
/// the evidence they read.
#[derive(Debug)]
pub struct AlleleTable {
    /// The locus's reference, gathered once and kept: every allele here was composed from
    /// it, and B3 attributes each sample's support against the same table.
    reference: LocusReferenceBases,
    /// The distinct alleles and the lookup that unified them, **the reference first**.
    alleles: AlleleLookup,
    /// How many of the cohort's reads were removed as evidence here, summed over the
    /// samples — see [`reads_removed_as_evidence`](Self::reads_removed_as_evidence).
    reads_removed_as_evidence: u32,
}

/// The alleles gathered so far, with the lookup that unifies them — one type, so that
/// pushing an allele and recording where it went cannot become two steps a caller can do
/// one of.
#[derive(Debug, Default)]
struct AlleleLookup {
    distinct: Vec<Box<[u8]>>,
    /// `distinct` keyed by its own bytes, so unification is a lookup rather than a scan.
    ///
    /// A scan would be faster at the two or three alleles most loci have, and slower where
    /// it matters: at the top of the committed cohort range a locus can hold one distinct
    /// allele per sample, and a scan makes the build quadratic in that count. Production
    /// keys the same table the same way (`per_group_merger.rs`). **Nothing iterates this
    /// map** — the table's order is `distinct`'s own — so the hasher's seed cannot reach
    /// the output.
    ///
    /// **Each distinct allele's bytes are held twice**, once here as the key and once in
    /// `distinct`. The alternative that removes it — one shared `Arc<[u8]>` — changes the
    /// type the architecture fixes for `CohortObservation::alleles` (`Vec<Box<[u8]>>`, arch
    /// §4), which would only hand the saving back at B3 if it converted. Left as it is
    /// until B3 has decided what it hands the caller, and measurable there: the cost is one
    /// extra copy of each *distinct* allele, not of each observation.
    by_bases: AHashMap<Box<[u8]>, usize>,
}

impl AlleleLookup {
    /// Add `bases` unless they are already here, and answer which allele they are.
    fn intern(&mut self, bases: &[u8]) -> usize {
        if let Some(&index) = self.by_bases.get(bases) {
            return index;
        }
        let index = self.distinct.len();
        self.distinct.push(Box::from(bases));
        self.by_bases.insert(Box::from(bases), index);
        index
    }
}

impl AlleleTable {
    /// Unify every sample's evidence over `locus` into one table.
    ///
    /// The reference goes in first, so it is allele 0 whether or not any sample's reads
    /// showed it — a cohort where every sample is homozygous for the variant still has to
    /// be genotyped against the reference. Then each sample in turn, in the run's sample
    /// order, contributes what its reads showed ([`alleles_of_sample`]), and every one of
    /// those is deduplicated against what is already there.
    ///
    /// **Unification by exact byte equality is only sound because indels were left-aligned
    /// upstream** (`LeftAlignPreparer`, `ng/read/left_align.rs`), and it is worth being
    /// exact about what that buys, because writing the same deletion at two placements
    /// *inside one locus* is not the failure. Projection already normalises that: in a
    /// homopolymer, dropping any one base gives the same sequence over the locus wherever
    /// it was anchored (pinned by
    /// `three_placements_of_one_deletion_inside_the_locus_unify`). What left-alignment
    /// buys is that the records **overlap and so chain into one locus at all**: two
    /// placements far enough apart close as two loci, each holding half the cohort's
    /// evidence, and no table ever sees them together (pinned by
    /// `two_placements_too_far_apart_to_chain_never_meet_in_one_table`). A site with two
    /// half-supported alleles reads as noisy data, not as a defect.
    pub fn over(locus: &ClosedLocus<'_>) -> Self {
        Self::assemble(locus).0
    }

    /// The table and every covering sample's support against it, in one pass over the
    /// samples — what [`CohortObservation::over`] is, with [`over`](Self::over) the same
    /// walk for a caller that wants the alleles alone.
    ///
    /// **One pass, not two.** Interning an allele answers which allele it is, so a sample's
    /// support can be accumulated as its alleles are derived; the alternative the
    /// [`index_of`](Self::index_of) doc describes — deriving again and looking the bytes up
    /// — would compose every read's allele a second time, which is the expensive half of
    /// this module. A sample's tally is grown as new alleles appear and padded to the
    /// table's width at the end, so the vectors are parallel however late an allele joins.
    fn assemble(locus: &ClosedLocus<'_>) -> (Self, Vec<SampleSupport>) {
        let reference = LocusReferenceBases::over(locus);
        let mut alleles = AlleleLookup::default();
        let reference_allele = alleles.intern(reference.bases());
        assert_eq!(
            reference_allele, REFERENCE_ALLELE,
            "the reference must be the locus's first allele",
        );

        let mut scratch = ReadAlleleScratch::default();
        let mut reads_removed_as_evidence = 0u32;
        let mut per_sample = Vec::with_capacity(locus.members.len());

        // One tally, refilled per sample rather than allocated per sample — the preference
        // the scratch beside it exists for. It cannot live *in* that scratch: the derivation
        // holds it mutably for the length of the call while the callback writes this.
        let mut tally: Vec<AlleleSupportTally> = Vec::new();

        for sample_members in &locus.members {
            let records = sample_members.observations;
            let sample = sample_members.sample;
            tally.clear();
            let mut composed_across_records = 0u32;

            let removed = alleles_of_sample(
                &reference,
                sample_members,
                &mut scratch,
                |bases, backing| {
                    let allele = alleles.intern(bases);
                    if tally.len() <= allele {
                        tally.resize(allele + 1, AlleleSupportTally::default());
                    }
                    match backing {
                        AlleleBacking::OneSequence(sequence) => tally[allele].add_whole(sequence),
                        AlleleBacking::OneRead { records, sightings } => {
                            composed_across_records = composed_across_records.saturating_add(1);
                            tally[allele]
                                .add_one_reads_share(share_of_one_read(sample, records, sightings));
                        }
                    }
                },
            );
            reads_removed_as_evidence = reads_removed_as_evidence.saturating_add(removed);

            per_sample.push(SampleSupport {
                sample,
                // **Only the alleles this sample showed.** A tally entry is reached only by
                // a sequence or by a read, and both add at least one read, so an untouched
                // entry is an allele some *other* sample introduced — and it needs no row
                // here, since asking about it answers no reads either way.
                supported: tally
                    .iter()
                    .enumerate()
                    .filter(|(_, tally)| tally.num_reads > 0)
                    .map(|(allele, tally)| SupportedAllele {
                        allele,
                        support: tally.finish(),
                    })
                    .collect(),
                reads_without_observation: records.iter().fold(0u32, |total, record| {
                    total.saturating_add(record.reads_without_observation)
                }),
                reads_removed_as_evidence: removed,
                reads_composed_across_records: composed_across_records,
            });
        }

        (
            Self {
                reference,
                alleles,
                reads_removed_as_evidence,
            },
            per_sample,
        )
    }

    /// How many reads were removed as evidence over this locus, summed across the samples.
    ///
    /// **A removal is lost depth, and it would otherwise be invisible.** The reads counted
    /// here covered some of their sample's records inside the locus and not all of them, so
    /// nothing they showed reaches the table and no later step can recover them by looking
    /// bytes up — a locus where many were removed reads as a quiet site with shallow
    /// samples. Spec §3.3 makes the same argument for the failed loci it insists are
    /// counted rather than dropped.
    ///
    /// **Where it surfaces is not decided here.** B3 owns per-sample support and C1 owns
    /// what a region reports, so this exists to carry the fact until one of them does
    /// something with it. Saturating, for the reason the walk's own total saturates: it is
    /// a diagnostic count, and reaching four billion removals at one locus is not a state
    /// worth a wider integer.
    pub fn reads_removed_as_evidence(&self) -> u32 {
        self.reads_removed_as_evidence
    }

    /// The locus's reference bases — allele 0, and what every allele here was composed
    /// from.
    pub fn reference(&self) -> &LocusReferenceBases {
        &self.reference
    }

    /// The distinct alleles, the reference first.
    pub fn alleles(&self) -> &[Box<[u8]>] {
        &self.alleles.distinct
    }

    /// Which allele `bases` is, if the table holds it.
    ///
    /// **Every allele every read showed is in the table**, so a later step can find which
    /// allele a read backs by composing the read again and looking the bytes up, rather
    /// than this build carrying an assignment for every observation. `None` means the
    /// sequence is not one the samples showed here.
    ///
    /// **It answers identity, and identity only.** How much support a sample lends an
    /// allele — the read counts and the per-read moments — is B3's, and the moments in
    /// particular cannot be recovered from the bytes: they are summed per observation, so
    /// where one observation's reads land on two alleles they have to be divided rather
    /// than looked up.
    pub fn index_of(&self, bases: &[u8]) -> Option<usize> {
        self.alleles.by_bases.get(bases).copied()
    }
}

/// What one building region delivers: the cohort observations built over it, and the ground
/// of the loci it refused (arch §4).
///
/// **Exactly one of these per region, even when every field is empty.** The organiser drains
/// regions in order and gathers what each refused, so a region that built nothing still has to
/// arrive — a missing one is a gap it cannot tell from a region with no variants
/// (spec §6.3, §3.3).
#[derive(Debug, Default)]
pub struct RegionOutcome {
    /// The survivors, in genome order. Loci are disjoint within one region, so their
    /// positions are a total order (spec §9).
    pub cohort_observations: Vec<CohortObservation>,
    /// **The ground of the loci that failed the width bound, in genome order.** A failed
    /// locus is an ordinary locus everywhere but emission (spec §3.2): nothing is built for
    /// it, and the run still has to be able to say it was refused — which is the whole
    /// meaning of the failed count (spec §3.3), the only signal that the bound is charging
    /// more than expected.
    ///
    /// **The arch also gives them a second job — displacing what overlaps them — and under
    /// the input contract stated at [`build_region`] there is nothing for that job to do.**
    /// Every builder is handed everything overlapping its ground, so a later locus that
    /// overlapped an earlier one would have chained into it and been skipped as the earlier
    /// region's; two loci owned by different regions cannot overlap at all. The spans are
    /// carried anyway, because the displacement rule is the organiser's (spec §6.1) and it
    /// is that component's business to decide whether it is a live rule or a safety net —
    /// but the reason to keep them **here** is the count, not the displacement.
    ///
    /// A locus that was **too quiet** is in neither vector: it is ground the caller examined
    /// and found empty, where a failure is ground it refused, and only one of the two is
    /// counted (spec §4.3, §1.3).
    pub failed_locus_spans: Vec<GenomeRegion>,
}

/// Build one region: close the cohort's loci over it, judge them, and assemble the survivors
/// (spec §6.2).
///
/// **A locus belongs to the region its first position falls in, and to no other.** That is
/// the whole of ownership: the same locus is closed by every builder whose observations reach
/// it, so without the rule it would be built more than once — and with it, a builder may
/// finish a locus outside its own region, which is what keeps a locus whole when a deletion
/// carries it past the end (spec §6.1). So the observations returned may reach beyond
/// `region`, and this reads past its own end to follow them.
///
/// `observations_per_sample` is one slice per sample in the run's sample order, each in
/// coordinate order, holding everything that overlaps `region` **and everything a locus
/// starting inside it can reach**. Whoever supplies them owns that guarantee: the walk cannot
/// know what it was not given, and a locus cut short by a missing observation is a wrong
/// answer rather than a failure. In the direct path the caller holds the whole stretch
/// (C2); from milestone D the observation cache draws each sample forward far enough
/// (spec §6.4).
///
/// **Deviation from the architecture's signature, and the plan's own order is why.** Arch §4
/// takes `&ObservationCache`; the cache is milestone D, and its `with_observations` hands out
/// exactly the slices this takes. Taking them directly is what lets the serial driver (C2) —
/// the oracle everything after it must reproduce — exist before the cache does.
///
/// **What it costs to hand over more than the region needs, measured.** The walk starts at
/// the beginning of the slices it is given, so every locus opening before `builder_region` is
/// closed in full and then discarded: about **3.3 µs per base of that prefix at 63 samples,
/// and 40 µs at 250** (the C1 review, release build). Handing every builder the whole stretch
/// therefore makes a run quadratic in its length — around **23 hours for a megabase at 63
/// samples**, against seconds when each builder is given a window over its own ground. That
/// is what the observation cache is for (milestone D); until it exists, the serial driver
/// hands over whole analysed regions rather than short ones, where the prefix is empty by
/// construction.
pub fn build_region(
    builder_region: GenomeRegion,
    observations_per_sample: &[&[SampleLocusObservations]],
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_obs: MinAltObs,
) -> RegionOutcome {
    let mut outcome = RegionOutcome::default();

    for locus in LocusCloser::over(observations_per_sample, max_cohort_locus_span, min_alt_obs) {
        // **Ownership, and the two ways a locus can fail to be ours.** One starting before
        // this builder's ground belongs to an earlier builder, which sees it whole; one
        // starting after its last base belongs to a later builder. Both are skipped rather
        // than built, so that every locus is built exactly once however the genome is
        // divided (spec §6.1, §9).
        //
        // **Both ends are inclusive**, and both need saying: a locus opening on the
        // builder's first base is its own, and so is one opening on its last. Widen either
        // comparison and that locus is skipped by *every* builder — lost from the run with
        // nothing to say so, and at twenty-base regions that is about one locus in twenty.
        // Pinned by `a_locus_opening_on_the_regions_first_base_belongs_to_that_region`.
        if locus.region.contig > builder_region.contig
            || (locus.region.contig == builder_region.contig
                && locus.region.start > builder_region.end)
        {
            // The walk yields loci in contig-then-position order, so nothing after this one
            // can start inside this ground either. A later builder owns the rest.
            break;
        }
        if locus.region.contig != builder_region.contig || locus.region.start < builder_region.start
        {
            continue;
        }

        match locus.verdict {
            Verdict::Build => outcome
                .cohort_observations
                .push(CohortObservation::over(&locus)),
            Verdict::Failed => outcome.failed_locus_spans.push(locus.region),
            Verdict::TooQuiet => {}
        }
    }

    outcome
}

/// One cohort locus, assembled: the ground, the alleles the cohort showed over it, and
/// what each covering sample's reads lend each allele (arch §4).
///
/// **This is evidence, not a call.** Nothing here says which alleles are real, what they are
/// worth, or what genotype a sample has; those are the calling steps' (spec §1.2).
#[derive(Debug)]
pub struct CohortObservation {
    /// The locus span, first position to furthest reach.
    pub region: GenomeRegion,
    /// The distinct alleles, the reference at [`REFERENCE_ALLELE`] (spec §4.2).
    pub alleles: Vec<Box<[u8]>>,
    /// The covering samples, in ascending sample order — **only** the covering ones.
    ///
    /// **A sample with no coverage over the span has no entry**, which is a different fact
    /// from an entry whose support is all reference and stays one (spec §4.2). The
    /// architecture's sketch says "indexed by the run's sample order", which reads two ways;
    /// this is the reading that keeps the distinction structural rather than resting on a
    /// zeroed row, and it is the shape the walk hands over
    /// ([`SampleMembers`](super::close::SampleMembers)). Each entry names its own sample.
    pub per_sample: Vec<SampleSupport>,
}

impl CohortObservation {
    /// Assemble one cohort locus: unify the alleles, then attribute every covering sample's
    /// reads to them.
    ///
    /// Only a locus the caller undertakes to build, for the reason
    /// [`LocusReferenceBases::over`] gives.
    pub fn over(locus: &ClosedLocus<'_>) -> Self {
        // Destructured rather than field-accessed, so that anything the table gains has to
        // be answered for here — dropped deliberately or carried — instead of vanishing.
        let (
            AlleleTable {
                reference: _,
                alleles,
                reads_removed_as_evidence: _,
            },
            per_sample,
        ) = AlleleTable::assemble(locus);
        Self {
            region: locus.region,
            alleles: alleles.distinct,
            per_sample,
        }
    }
}

/// One sample's evidence at one cohort locus, against that locus's allele table.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleSupport {
    /// Which sample — its index in the run's sample order, as the walk named it.
    pub sample: usize,
    /// **The alleles this sample's reads showed, and only those**, in ascending allele order
    /// (the owner's ruling at Checkpoint B). An allele it did not show has no entry, which
    /// [`support_for`](Self::support_for) reads as no reads and no sums — the same answer a
    /// zeroed row would give, without a row per sample per allele.
    ///
    /// **What that saves is proportional to how much of the table each sample missed.** A
    /// row for every allele costs samples × alleles entries of 32 bytes whatever the cohort
    /// showed; this costs one entry per allele a sample actually showed, which at an ordinary
    /// locus of two or three alleles is one or two per sample.
    ///
    /// **Support is never merged across alleles**, because a genotype likelihood needs them
    /// apart (spec §4.2). It *is* merged within one: where two of the sample's own
    /// observations reached the same allele — the same bases from two read groups, or two
    /// records' worth of one read — their reads and their sums are added.
    pub supported: Vec<SupportedAllele>,
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
    /// nothing they showed reaches `supported` (see [`alleles_of_sample`]). Lost depth,
    /// counted rather than inferred from an absence. Saturating.
    pub reads_removed_as_evidence: u32,
    /// Reads whose allele was composed across several of this sample's records, and whose
    /// five quality sums in `supported` are therefore **divided rather than measured** (see
    /// [`AlleleSupport`]). Zero on a one-record sample, where every sum is exact. Saturating.
    pub reads_composed_across_records: u32,
}

impl SampleSupport {
    /// What this sample's reads lend `allele` — **nothing, where it showed that allele no
    /// reads**, which is the answer a zeroed row would have given.
    ///
    /// A scan, not an index, because [`supported`](Self::supported) holds only the alleles
    /// the sample showed: at an ordinary locus that is one or two entries, and it is bounded
    /// by how many alleles the whole cohort showed.
    pub fn support_for(&self, allele: usize) -> AlleleSupport {
        self.supported
            .iter()
            .find(|supported| supported.allele == allele)
            .map(|supported| supported.support)
            .unwrap_or_default()
    }
}

/// One allele of the locus, and what one sample's reads lend it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportedAllele {
    /// Which allele — an index into [`CohortObservation::alleles`].
    pub allele: usize,
    /// What this sample's reads lend it. Never zero: an allele a sample showed no reads for
    /// has no entry at all.
    pub support: AlleleSupport,
}

/// What one sample's reads lend one allele.
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
///   with production's: see [`share_of_one_read`];
/// - **the strand and mapping-quality sums** take each read's share as the pooled mean over
///   the sequences it showed: the sums over those sequences divided by their read counts.
///   Both are properties of the read and identical at each of its sightings, so the pooling
///   estimates one number rather than mixing two;
/// - **the placement count is pooled the same way, and it is the one that mixes two
///   questions.** It is counted against each record's own position (see `placed_left`), so
///   two records of one sample ask "did the read start left of here?" about different
///   *heres*, and the pooled answer is an approximation whose disagreement is bounded by the
///   locus's width — 50 bases at the default bound.
///
/// Rounded once per allele rather than per read, so the counts stay as close to whole reads
/// as the division allows.
///
/// **One row per allele pools the read groups, and the STR path will want them apart again**
/// (the owner, at Checkpoint C). A sample's reads at one allele may come from several read
/// groups, and `SequenceObservation` keeps them apart deliberately — "a per-chemistry model
/// needs the allele × group cross **with its quality moments**"
/// ([`locus_generation/mod.rs`](../../locus_generation/mod.rs)). **Stutter is fitted per read
/// group**, so an STR locus called from a pooled row would be scored against a stutter rate
/// that belongs to no group in particular.
///
/// It is pooled here anyway, for now, and deliberately: the architecture's own sketch pools
/// them (arch §4, `per_allele: Vec<SequenceObservation>`, one entry per allele), nothing in
/// the generic path needs the cross, and the shape that restores it is one row per
/// `(allele, read group)`, which folds to this one — so the step that needs it can add it
/// without unpicking anything. **Whoever brings the STR path through this merge owes that
/// change.**
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AlleleSupport {
    /// How many of the sample's reads showed this allele. **Exact** — every read is named,
    /// so it is counted where it belongs rather than divided. Saturating, which needs four
    /// billion reads on one allele of one sample at one locus.
    pub num_reads: u32,
    /// Of those, how many were on the forward strand — strand bias.
    pub num_fwd: u32,
    /// Σ per-read log-error over them — the freebayes per-read error term.
    pub q_sum: f64,
    /// Σ MAPQ over them.
    pub mapq_sum: u32,
    /// Σ MAPQ² over them; with `mapq_sum` and `num_reads` it recovers the mean and variance
    /// the MAPQ multi-mapper filter reads.
    pub mapq_sum_sq: u64,
    /// How many started strictly left of **the record they were seen at** — freebayes'
    /// `placedLeft`, and the read-position-bias term.
    ///
    /// **Not counted against the cohort locus's own first base.** The mint counts it against
    /// each record's position — `alignment_start < walker_pos` on the generic path
    /// (`locus_generation/pileup/fast_column.rs`), the tract's anchor on the STR path
    /// (`locus_generation/ssr.rs`) — and nothing here re-anchors it, so where a locus spans
    /// several of a sample's records this mixes as many questions as it has records.
    pub placed_left: u32,
}

/// One sample's support as it is accumulated, before the divided sums are rounded.
///
/// The four count-like sums are gathered as `f64` because a read's share of them is a
/// fraction until every read has been added; rounding each read's share on its own would
/// lose up to half a read per read rather than per allele.
#[derive(Debug, Clone, Copy, Default)]
struct AlleleSupportTally {
    num_reads: u32,
    sums: SupportSums,
}

impl AlleleSupportTally {
    /// Everything one sequence measured, added whole — the exact case, where every read
    /// behind the sequence showed this allele and nothing has to be divided.
    fn add_whole(&mut self, sequence: &SequenceObservation) {
        self.num_reads = self.num_reads.saturating_add(sequence.num_obs);
        self.sums.add_counts_of(sequence);
        // The quality is added whole here, where a read's share takes the weakest sighting's
        // mean instead — the one place the two paths differ, stated at both of them.
        self.sums.q_sum += sequence.q_sum;
    }

    /// One read's share, divided out of the sequences it showed across the locus.
    fn add_one_reads_share(&mut self, share: SupportSums) {
        self.num_reads = self.num_reads.saturating_add(1);
        self.sums.add(share);
    }

    /// Round the divided sums back to whole numbers, once, when the allele is finished.
    ///
    /// This is the one place the five are spelled out apart, because it is the one place
    /// they differ: three are counts of reads and round, one is a sum of squares and rounds
    /// wider, and the quality is a real number that does not round at all.
    fn finish(self) -> AlleleSupport {
        AlleleSupport {
            num_reads: self.num_reads,
            num_fwd: round_to_u32(self.sums.num_fwd),
            q_sum: self.sums.q_sum,
            mapq_sum: round_to_u32(self.sums.mapq_sum),
            mapq_sum_sq: round_to_u64(self.sums.mapq_sum_sq),
            placed_left: round_to_u32(self.sums.placed_left),
        }
    }
}

/// The five quality sums, held as `f64` while they are added up.
///
/// **One type for all three jobs they do** — what a sequence measured, what one read's share
/// of several sequences comes to, and what an allele has accumulated — so the list of five
/// is written once for gathering and once for rounding, rather than at every site that
/// touches them. A sixth sum (`placed_start`, which `SequenceObservation` declines to carry
/// today) would then be added in two places, both of which fail to compile if it is
/// forgotten.
///
/// `f64` because a read's share of a sum is a fraction until every read has been added;
/// rounding each read's share on its own would lose up to half a read per read rather than
/// per allele.
#[derive(Debug, Clone, Copy, Default)]
struct SupportSums {
    num_fwd: f64,
    q_sum: f64,
    mapq_sum: f64,
    mapq_sum_sq: f64,
    placed_left: f64,
}

impl SupportSums {
    /// Add the four count-like sums one sequence measured. **The quality is left alone**,
    /// because the two callers do different things with it: everything one sequence measured
    /// is added whole, where one read's share takes the weakest sighting's mean rather than
    /// a sum. Each says so where it stands.
    fn add_counts_of(&mut self, sequence: &SequenceObservation) {
        self.num_fwd += f64::from(sequence.num_fwd);
        self.mapq_sum += f64::from(sequence.mapq_sum);
        // No `From<u64> for f64` exists because the conversion is lossy above 2^53. It
        // cannot be reached here: MAPQ is at most 60, so this would need about 10^12 reads
        // behind one sequence.
        self.mapq_sum_sq += sequence.mapq_sum_sq as f64;
        self.placed_left += f64::from(sequence.placed_left);
    }

    /// **Destructured rather than field-accessed**, so that a sixth sum added to this type
    /// is a compile error here instead of a zero nobody notices.
    fn add(&mut self, other: Self) {
        let Self {
            num_fwd,
            q_sum,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
        } = other;
        self.num_fwd += num_fwd;
        self.q_sum += q_sum;
        self.mapq_sum += mapq_sum;
        self.mapq_sum_sq += mapq_sum_sq;
        self.placed_left += placed_left;
    }

    /// Divide the four count-like sums by the reads they were measured over — **not the
    /// quality**, which is not pooled but taken from the weakest sighting
    /// ([`share_of_one_read`]). Destructured for the reason [`add`](Self::add) is.
    fn divide_counts_by(&mut self, reads: f64) {
        let Self {
            num_fwd,
            q_sum: _,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
        } = self;
        *num_fwd /= reads;
        *mapq_sum /= reads;
        *mapq_sum_sq /= reads;
        *placed_left /= reads;
    }
}

/// One read's share of the sequences it showed across the locus — the division
/// [`AlleleSupport`] documents.
///
/// `sightings` is the read's own, one per record of the sample. Every one of them names a
/// complete sequence with at least one read behind it, both of which the derivation enforces
/// before a sighting is recorded ([`alleles_of_sample`]).
fn share_of_one_read(
    sample: usize,
    records: &[SampleLocusObservations],
    sightings: &[ReadSighting],
) -> SupportSums {
    let mut weakest_mean_quality: Option<f64> = None;
    let mut reads_behind = 0f64;
    let mut share = SupportSums::default();

    for sighting in sightings {
        let sequence = &records[sighting.record as usize].observations[sighting.sequence as usize];
        // **A backstop, not the guard.** A sequence with no reads is skipped before a
        // sighting is ever recorded for it, so this cannot fire from the derivation.
        //
        // It is here because of what happens without it, which is not a failure and not an
        // infinity either: dividing by zero reads gives `-inf`, and the `max` below then
        // **discards it** in favour of the read's other sighting, so the allele comes back
        // with a plausible quality measured over reads that were not counted. Measured, by
        // removing this assertion: a read whose two sightings are a zero-read sequence and a
        // −0.5 one reports −0.5.
        //
        // That a sequence's read count agrees with the reads it names is the **producer's**
        // guarantee, like the reference width and the chain ids, so **when observations are
        // decoded from a psp file this must become a `RunError`** beside
        // `ObservationExceedsReachCeiling` (arch §5) rather than a panic.
        assert!(
            sequence.num_obs > 0,
            "the observation {:?} at {} of sample index {} names reads but counts none, so \
             no read's share of it can be taken",
            String::from_utf8_lossy(&sequence.bases),
            records[sighting.record as usize].region,
            sample,
        );
        let reads = f64::from(sequence.num_obs);
        let mean_quality = sequence.q_sum / reads;
        weakest_mean_quality = Some(match weakest_mean_quality {
            Some(weakest) => weakest.max(mean_quality),
            None => mean_quality,
        });

        reads_behind += reads;
        // The four counts only: the quality is not pooled, and adding it here would be
        // added and then overwritten below — the silent kind of dead work.
        share.add_counts_of(sequence);
    }

    // The quality is the weakest sighting's mean read, not the pooled one: a read's evidence
    // for an allele spanning several records is only as good as the least it saw well, since
    // the allele is wrong if any of its pieces is.
    //
    // **This is production's stated rule and not production's code.** Its merger plan says
    // the compound's effective quality "cannot exceed any single constituent's"
    // (`doc/devel/implementation_plans/cohort_per_group_merger.md`, step 3), and its code
    // takes `min` over the constituents' mean `q_sum`
    // (`project_compound_scalars`, `var_calling/per_group_merger.rs`). These are opposites
    // here: `q_sum` is a sum of `ln P(error)` (`pileup_record.rs`), so it is negative and a
    // *weaker* read is a *larger* number — `min` picks the constituent the read saw best,
    // which makes a compound allele look better evidenced than any single piece of it. The
    // divergence is deliberate and recorded in this step's report.
    //
    // No fallback for an empty `sightings`: a read's sightings are one group of a `chunk_by`
    // over the reads seen, which never yields an empty group, and `0.0` would not be a
    // neutral answer if it did — `ln P(error) = 0` is `P(error) = 1`, the worst quality
    // expressible, indistinguishable from a measured one.
    share.q_sum =
        weakest_mean_quality.expect("a read's sightings are one non-empty group of the reads seen");
    share.divide_counts_by(reads_behind);
    share
}

/// Round a divided count back to a whole one.
///
/// **Not a repair**: Rust's float-to-integer `as` has saturated rather than wrapped since
/// 1.45, so this answers what `value.round() as u32` answers on every input — including
/// `NaN`, both infinities, negatives and both boundaries, which is measured rather than
/// assumed (`the_rounding_of_a_divided_count_saturates_at_both_ends`). What it buys is that
/// the boundary is written where a reader can see it, and that `NaN` has a stated answer
/// rather than an inherited one.
fn round_to_u32(value: f64) -> u32 {
    let rounded = value.round();
    if rounded <= 0.0 {
        0
    } else if rounded >= f64::from(u32::MAX) {
        u32::MAX
    } else {
        rounded as u32
    }
}

/// As [`round_to_u32`], for the sum of squares.
fn round_to_u64(value: f64) -> u64 {
    let rounded = value.round();
    if rounded <= 0.0 {
        0
    } else if rounded >= u64::MAX as f64 {
        u64::MAX
    } else {
        rounded as u64
    }
}

/// Working space for deriving one sample's alleles, refilled per sample rather than
/// allocated per sample — the pattern the walk beside it uses for the same reason
/// (`close.rs`'s `cursors_at_open`). It is allocated once per locus today, and the builder
/// that owns a whole region is where it can be hoisted further (plan C1).
#[derive(Debug, Default)]
struct ReadAlleleScratch {
    /// One entry per (read, record) pair of the sample being derived, sorted.
    ///
    /// **A sorted list rather than a map of lists**, which is the shape the pileup walk
    /// settled on for the same grouping-by-chain-id question and the same reason: a map
    /// allocates a vector per read, where this is one buffer refilled for every read of
    /// every sample of the locus — and one buffer for the whole run once a builder owns the
    /// scratch (plan C1) (`resolve_mate_overlap_at_pos`,
    /// `locus_generation/pileup/genome_walk.rs`).
    by_read: Vec<ReadSighting>,
    /// The allele being composed. One buffer for the whole locus, refilled per read.
    composed: Vec<u8>,
}

/// One read, seen at one of its sample's records inside the locus, showing one sequence
/// there.
///
/// **The field order is the sort order**, and the sort is the algorithm: by read, so that a
/// read's sightings come together, then by record, so that they arrive in the order the
/// allele is composed in. Naming the fields is what keeps that a declaration rather than a
/// property of a tuple's spelling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct ReadSighting {
    /// Which read — comparable **within one sample only**, the id space being per file
    /// (`chain_id_allocator.rs`).
    read: ChainId,
    /// Which of the sample's records inside the locus, by position in
    /// [`SampleMembers::observations`].
    record: u32,
    /// Which of that record's sequences the read showed, by position in
    /// [`SampleLocusObservations::observations`].
    sequence: u32,
}

/// Every allele one sample's reads showed over the whole locus, handed to `emit` one at a
/// time and in a fixed order.
///
/// **The rule is the owner's, and it has exactly two branches** (2026-08-17): either we
/// know a read covered the locus, and its allele is what it showed elongated across the
/// locus, or we know it did not, and it is removed as evidence. What decides which is the
/// read's presence at the sample's own records inside the locus:
///
/// - **A read named at every one of them** showed something at each, and those somethings
///   are composed in coordinate order into one allele, with any ground between two records
///   filled from the locus's reference.
///
///   **That filling is unreachable on the generic path, and the reason is worth knowing.**
///   The generic mint writes a record at *every position a read covered*, not only where a
///   read departed from the reference — measured on minted records in
///   `serial.rs`'s `a_cohort_observation_is_built_from_minted_observations`, where a
///   thirty-base read yields thirty records. So a gap in a sample's records is ground its
///   reads did not cover, and a read named at two records is named at every record between
///   them: there is nothing left to fill. The code fills it anyway, because this walk does
///   not depend on which generator minted its input, and a mint that recorded only departures
///   would need it.
/// - **A read missing from any one of them is removed**, whatever the reason. It may never
///   have covered that position; a depth cap may have discarded it there, the cap acting
///   per position and leaving no identities (`reads_discarded_by_cap`); or it may have been
///   [`Partial`](ReadWitness::Partial) there, having seen only part of that record. All
///   three are the same fact — what the read showed there was never recorded — so all three
///   take the same branch.
///
/// **Where the sample has one record inside the locus, no read is consulted** and each of
/// that record's sequences is projected on its own. There is nothing to elongate across, so
/// nothing to place a read at; it is the ordinary generic locus, and it is the **only**
/// answer available on the **STR path**, which records no chain ids at all and needs none —
/// an STR locus is one record, so [`ReadWitness`] already says whether a read spanned it.
///
/// **The two branches differ on one shape, and it is not an optimisation of the other.** A
/// fragment whose mates overlap and disagree is two sequences of one record under one chain
/// id (the mint keeps both, `resolve_mate_overlap_at_pos`). With one record both are
/// emitted, each mate's sequence being complete evidence over the locus on its own; with
/// several the fragment showed no one thing at that record and is removed, because
/// composing across records needs the fragment as a unit. Pinned by
/// `a_read_showing_two_things_at_a_sole_record_still_contributes_both` beside
/// `a_read_that_showed_two_things_at_one_record_is_removed`.
///
/// **Two things this does not claim.** A read is judged on its presence at the sample's
/// records, not on coverage of the locus: the ground outside them — before the first
/// record, between two, after the last — is written from the reference because this sample
/// minted no record there, which says its reads did not depart from the reference, not that
/// this read was present. And a fragment whose two mates flank an unsequenced insert inside
/// one locus is treated as having covered it, there being nothing in the observations that
/// could say otherwise.
///
/// **Two reads that showed the same thing yield the same bytes twice** — deduplicating is
/// the table's job, not this walk's.
///
/// The order is fixed: by chain id where the reads were consulted, by the record's own
/// order of sequences where they were not.
///
/// **Answers how many of the sample's reads were removed**, so that the loss can be counted
/// rather than inferred from an absence
/// ([`AlleleTable::reads_removed_as_evidence`]). Nothing is removed on the one-record
/// branch, where no read is consulted.
fn alleles_of_sample(
    reference: &LocusReferenceBases,
    members: &SampleMembers<'_>,
    scratch: &mut ReadAlleleScratch,
    mut emit: impl FnMut(&[u8], AlleleBacking<'_>),
) -> u32 {
    let ReadAlleleScratch { by_read, composed } = scratch;
    let records = members.observations;
    by_read.clear();

    if let [only_record] = records {
        let placement = reference.placing(only_record);
        for sequence in placement.projectable_sequences() {
            // **A sequence no read is behind is not an allele the sample showed.** Both
            // branches skip it, and they have to agree: interning it here would put an
            // allele in the table that nothing supports, carrying a quality nobody measured,
            // where the branch below would refuse it — the same input answered two ways. No
            // producer emits one; both mints derive a sequence from the reads that showed
            // it.
            if sequence.num_obs == 0 {
                continue;
            }
            placement.project_into(sequence, composed);
            emit(composed, AlleleBacking::OneSequence(sequence));
        }
        return 0;
    }

    for (record_index, record) in records.iter().enumerate() {
        for (sequence_index, sequence) in record.observations.iter().enumerate() {
            if sequence.read_witness != ReadWitness::Complete || sequence.num_obs == 0 {
                continue;
            }
            // **A read the mint did not name cannot be placed, and dropping it would be
            // silent** — the observation would simply contribute no allele, and a locus
            // where every sample lost one would look like a quiet site rather than a
            // defect. B0 names every read the generic mint folds, reference-matching ones
            // included, so what this catches is a regression there.
            //
            // **It cannot fire on the STR path, whose ids are empty by design**, because a
            // locus that reaches this branch has two of one sample's records in it and an
            // STR sample cannot have two: an STR record is a segment's tract, segments are
            // the reference's own partition and no observation crosses a boundary (run spec
            // §4.3), so two of them never overlap and never chain (the same fact `close.rs`
            // rests its never-mixed-kinds check on).
            //
            // Release-level, like the walk's own checks beside it, and one comparison per
            // sequence against a derivation that already copies every sequence's bases.
            // **When observations are decoded from a psp file this becomes corrupt input
            // and must become a `RunError`** (arch §5), as the reference-width check
            // above must.
            assert!(
                !sequence.chain_ids.is_empty(),
                "the observation {:?} at {} of sample index {} carries {} reads and no \
                 chain id, so the reads that showed it cannot be placed across the locus \
                 {}",
                String::from_utf8_lossy(&sequence.bases),
                record.region,
                members.sample,
                sequence.num_obs,
                reference.region(),
            );
            for read in &sequence.chain_ids {
                // PANIC-FREE: both indices are into vectors this loop is walking, and a
                // locus with four billion records or an observation with four billion
                // sequences cannot be built — a record is at least one base wide and the
                // locus span is bounded by the verdict.
                by_read.push(ReadSighting {
                    read: *read,
                    record: record_index as u32,
                    sequence: sequence_index as u32,
                });
            }
        }
    }
    // Brings each read's sightings together and puts them in record order inside the read,
    // which is the order they are composed in. The keys are unique — a read is named at
    // most once per sequence, the ids being deduplicated as they are folded — so the sort
    // needs no stability to be deterministic.
    by_read.sort_unstable();

    // Worked out once per record rather than once per (read, record), which is what
    // [`LocusReferenceBases::placing`] says it is for. It cannot live in the scratch above:
    // a placement borrows the locus's reference and the record it places, so a reusable
    // buffer of them would put those lifetimes on the scratch and on every caller.
    let placements: Vec<MemberPlacement<'_>> = records
        .iter()
        .map(|record| reference.placing(record))
        .collect();

    let mut removed = 0u32;
    for shown in by_read.chunk_by(|left, right| left.read == right.read) {
        // **The read is evidence only if it showed something at every one of the sample's
        // records, and exactly one thing at each.** The sightings are sorted, so the first
        // half of that is `shown[i]` naming record `i` for every `i`; the second half rides
        // on it, since a read naming one record twice pushes another out of the count. That
        // second case is a fragment whose two mates overlap and disagree — the mint keeps
        // both as observations under the one chain id it gave the fragment
        // (`resolve_mate_overlap_at_pos`) — and there is no one thing it showed, so it is
        // removed like a read that was not there.
        if shown.len() != records.len()
            || shown
                .iter()
                .enumerate()
                .any(|(expected, sighting)| sighting.record as usize != expected)
        {
            removed = removed.saturating_add(1);
            continue;
        }

        composed.clear();
        let mut written_to = 0;
        for sighting in shown {
            let record = &records[sighting.record as usize];
            written_to = placements[sighting.record as usize].compose_into(
                &record.observations[sighting.sequence as usize],
                written_to,
                composed,
            );
        }
        composed.extend_from_slice(&reference.bases()[written_to..]);
        emit(
            composed,
            AlleleBacking::OneRead {
                records,
                sightings: shown,
            },
        );
    }

    removed
}

/// What backed one allele the derivation emitted — the evidence behind the bytes, so that a
/// caller attributing support need not compose them a second time to find out.
///
/// The two arms are the two branches of [`alleles_of_sample`], and the difference matters to
/// support rather than to identity: what a sequence measured is exact, where what one read
/// showed across several records has to be divided out of the sequences it was seen in
/// ([`AlleleSupport`]).
#[derive(Debug, Clone, Copy)]
enum AlleleBacking<'a> {
    /// Every read behind one sequence of the sample's sole record.
    OneSequence(&'a SequenceObservation),
    /// One read, with its sighting at each of the sample's records, in coordinate order.
    ///
    /// **The records travel with the sightings**, because a sighting is a pair of indices
    /// and means nothing without the slice they index. A caller that supplied its own slice
    /// would be right only by coincidence of expression — the same defect B1's review found
    /// in the projection, where a member and a sequence were two loose arguments.
    OneRead {
        records: &'a [SampleLocusObservations],
        sightings: &'a [ReadSighting],
    },
}

/// How far `member_region` starts into `locus_region`, in bases.
///
/// **Were the subtraction open-coded at each use**, a member starting before the locus
/// would wrap in the release profile — where overflow checks are off — and index the
/// gathered reference from somewhere near `usize::MAX`. Both are 1-based positions in a
/// type with public fields and no constructor enforcing anything, so `checked_sub` turns
/// that into a message naming both regions.
fn offset_within(locus_region: GenomeRegion, member_region: GenomeRegion) -> usize {
    let offset = member_region
        .start
        .get()
        .checked_sub(locus_region.start.get())
        .unwrap_or_else(|| {
            panic!("the region {member_region} starts before the locus {locus_region}")
        });
    // PANIC-FREE: the offset is smaller than the locus's own span, which converted when
    // the reference was gathered.
    usize::try_from(offset).expect("an offset inside a locus fits a usize")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::run::cohort_merge::close::LocusCloser;
    use crate::ng::run::cohort_merge::{MaxCohortLocusSpan, MinAltObs};
    use crate::ng::types::{ContigId, Position, ReadGroupId};

    fn region_on(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        region_on(0, start, end)
    }

    /// One distinct sequence the reads showed, complete over the region it sits in.
    ///
    /// Only `bases` and `read_witness` matter to projection, so every other field is held
    /// at a fixed value rather than varied.
    fn sequence(bases: &[u8]) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: 3,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// **One sample's observation over `region`** — the spec's *member* (§4.2): one
    /// `SampleLocusObservations`, carrying the reference over its own region and the one
    /// sequence its reads showed.
    ///
    /// Unlike the fixtures in `close.rs`, the *bases* are what these tests are about, so
    /// each is spelled out and the reference's width is checked against the region here
    /// rather than filled in — a fixture whose reference is a base short would otherwise
    /// fail inside the code under test and read as a defect in it.
    fn member(
        region: GenomeRegion,
        reference_bases: &[u8],
        observed_bases: &[u8],
    ) -> SampleLocusObservations {
        assert_eq!(
            reference_bases.len() as u64,
            span_of(region),
            "the fixture's reference must cover its own region",
        );
        SampleLocusObservations {
            region,
            reference_bases: Box::from(reference_bases),
            observations: vec![sequence(observed_bases)],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// A closed locus over `region` whose members are the given observations, one entry
    /// per sample in the order they are passed, judged `Build`.
    ///
    /// Built by hand rather than walked, so that a fixture can present ground the walk
    /// would never close — an uncovered position, a member reaching past the end — which
    /// is what the assertions in [`LocusReferenceBases::over`] exist for.
    fn closed_locus<'a>(
        region: GenomeRegion,
        per_sample: &[&'a [SampleLocusObservations]],
    ) -> ClosedLocus<'a> {
        ClosedLocus {
            region,
            members: per_sample
                .iter()
                .copied()
                .enumerate()
                .map(|(sample, observations)| SampleMembers {
                    sample,
                    observations,
                })
                .collect(),
            non_reference_reads: 0,
            verdict: Verdict::Build,
        }
    }

    /// One sequence, projected onto the locus's whole span.
    fn project(
        reference: &LocusReferenceBases,
        member: &SampleLocusObservations,
        sequence: &SequenceObservation,
    ) -> Vec<u8> {
        let mut buffer = Vec::new();
        reference
            .placing(member)
            .project_into(sequence, &mut buffer);
        buffer
    }

    /// The one sequence a single-sequence fixture carries, projected.
    fn project_only(reference: &LocusReferenceBases, member: &SampleLocusObservations) -> Vec<u8> {
        project(reference, member, &member.observations[0])
    }

    /// **A SNP inside a deletion's span projects onto the whole span** — the plan's first
    /// B1 test, and the case the merge exists for: two samples that recorded different
    /// widths of the same ground become two sequences over the same ground.
    ///
    /// The locus reads `ACGTA` at 10–14. Sample 0 deleted `CGTA` and so recorded `A` over
    /// all five bases; sample 1 has `G` → `T` at position 12 alone. Projected, the SNP is
    /// `ACTTA` — its own base with two reference bases either side — and the deletion is
    /// still `A`, since it already covered the whole locus.
    #[test]
    fn a_snp_inside_a_deletions_span_projects_onto_the_whole_span() {
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let snp = [member(region(12, 12), b"G", b"T")];
        let locus = closed_locus(region(10, 14), &[&deletion, &snp]);

        let reference = LocusReferenceBases::over(&locus);
        assert_eq!(reference.bases(), b"ACGTA");

        assert_eq!(
            project_only(&reference, &snp[0]),
            b"ACTTA",
            "two reference bases on the left, two on the right",
        );
        assert_eq!(
            project_only(&reference, &deletion[0]),
            b"A",
            "it already spanned the locus, so nothing is padded onto it",
        );
    }

    /// **An insertion's reference span stays its anchor base** — the plan's second B1
    /// test. The inserted bases lengthen the *allele*, never the ground it replaces.
    ///
    /// The same 10–14 locus. Sample 1 carries `GTTT` at position 12: the anchor base `G`
    /// plus three inserted `T`s. Projected it is `ACGTTTTA` — eight bases of allele over
    /// five bases of reference — and the locus is still five bases wide, which is why an
    /// insertion cannot push a locus past `max_cohort_locus_span` however long it is
    /// (spec §3.1).
    #[test]
    fn an_insertion_projects_at_its_anchor_base_and_leaves_the_span_alone() {
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let insertion = [member(region(12, 12), b"G", b"GTTT")];
        let locus = closed_locus(region(10, 14), &[&deletion, &insertion]);

        let reference = LocusReferenceBases::over(&locus);

        assert_eq!(project_only(&reference, &insertion[0]), b"ACGTTTTA");
        assert_eq!(span_of(reference.region()), 5, "the locus is unchanged");
    }

    /// **No single member need cover the locus** — the reference is gathered across them,
    /// which is what the spec's "the bases travel on the observation" means once a locus
    /// is wider than any one of its members.
    ///
    /// Two observations, 10–12 and 11–14, overlapping at 11–12. Between them they cover
    /// 10–14 and gather the whole `ACGTA`; neither member alone carries all five of those
    /// bases.
    #[test]
    fn the_reference_is_gathered_across_members_that_each_cover_part_of_the_locus() {
        let left_member = [member(region(10, 12), b"ACG", b"AAG")];
        let right_member = [member(region(11, 14), b"CGTA", b"CGTT")];
        let locus = closed_locus(region(10, 14), &[&left_member, &right_member]);

        let reference = LocusReferenceBases::over(&locus);

        assert_eq!(reference.bases(), b"ACGTA");
        assert_eq!(
            project_only(&reference, &left_member[0]),
            b"AAGTA",
            "the two bases on its right came from the other sample's member",
        );
    }

    /// **One sample can hold two observations inside one locus**, and both carry reference
    /// bases the gather needs — spec §4.2's "two of its own observations", which B3 sums
    /// support over. The test above spreads coverage across *samples*; this spreads it
    /// within one, which is the other loop.
    #[test]
    fn one_samples_two_observations_both_reach_the_gather() {
        let two = [
            member(region(10, 11), b"AC", b"AC"),
            member(region(12, 14), b"GTA", b"G"),
        ];
        let locus = closed_locus(region(10, 14), &[&two]);

        assert_eq!(LocusReferenceBases::over(&locus).bases(), b"ACGTA");
    }

    /// **A member that matched the reference projects to the reference allele itself** —
    /// what `bases()` claims, and what lets the next step's allele table hold the
    /// reference without a special case (spec §4.2). At a non-zero offset, so a prefix or
    /// suffix off by one would show.
    #[test]
    fn a_member_that_matched_the_reference_projects_to_the_locus_reference() {
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let matched = [member(region(12, 12), b"G", b"G")];
        let locus = closed_locus(region(10, 14), &[&deletion, &matched]);

        let reference = LocusReferenceBases::over(&locus);
        assert_eq!(project_only(&reference, &matched[0]), reference.bases());
    }

    /// A member that already covers the whole locus projects to its own bases untouched —
    /// the ordinary case at a one-position locus, and the one where padding must add
    /// nothing.
    #[test]
    fn a_member_covering_the_whole_locus_projects_to_its_own_bases() {
        let only_member = [member(region(10, 10), b"A", b"C")];
        let locus = closed_locus(region(10, 10), &[&only_member]);

        let reference = LocusReferenceBases::over(&locus);
        assert_eq!(reference.bases(), b"A");
        assert_eq!(project_only(&reference, &only_member[0]), b"C");
    }

    /// **A reference containing `N` is gathered like any other.** ng's fetch folds every
    /// non-ACGT byte to `N`, so a locus over an assembly gap carries them — and a sentinel
    /// spelled `b'N'` would refuse a locus its members do cover. Nothing else in this file
    /// distinguishes the sentinel's value from any other non-base byte.
    #[test]
    fn a_reference_containing_n_bases_is_gathered_like_any_other() {
        let over_a_gap = [member(region(10, 14), b"ACNTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&over_a_gap]);

        assert_eq!(LocusReferenceBases::over(&locus).bases(), b"ACNTA");
    }

    /// **The buffer is emptied before it is refilled.** One buffer serves a whole locus,
    /// so a projection that appended to what the last one left would grow an allele per
    /// member — and the first few members would still look right, which is what makes it
    /// worth pinning. The swing is a 13-base insertion allele followed by a 1-base
    /// deletion allele, so a rewrite that truncated around the old contents rather than
    /// clearing would show.
    #[test]
    fn a_reused_buffer_holds_only_the_latest_projection() {
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let long_insertion = [member(region(12, 12), b"G", b"GTTTTTTTT")];
        let locus = closed_locus(region(10, 14), &[&deletion, &long_insertion]);
        let reference = LocusReferenceBases::over(&locus);

        let mut buffer = Vec::new();
        reference
            .placing(&long_insertion[0])
            .project_into(&long_insertion[0].observations[0], &mut buffer);
        assert_eq!(buffer, b"ACGTTTTTTTTTA");

        reference
            .placing(&deletion[0])
            .project_into(&deletion[0].observations[0], &mut buffer);
        assert_eq!(
            buffer, b"A",
            "the previous projection is gone, not appended to or truncated around",
        );
    }

    /// **The gather and the projection both measure the span the way `close.rs` does, not
    /// through `GenomeRegion::len()`** — which computes `end + 1` before subtracting and
    /// so answers 0 at the top of the coordinate space in the release profile, where
    /// overflow checks are off (that method's own doc).
    ///
    /// A four-base locus at the ceiling. Through `len()` the gathered reference would be
    /// empty; and the member ending *at* the ceiling is the one the projection needs, since
    /// for the SNP one base short of it the two spellings agree — with `len()` its
    /// `covered` would be 0, so the whole reference would follow its bases and the allele
    /// would come back `AACGT` instead of `A`.
    #[test]
    fn a_locus_at_the_coordinate_ceiling_is_gathered_and_padded_on_its_true_width() {
        let whole = [member(region(u64::MAX - 3, u64::MAX), b"ACGT", b"A")];
        let snp = [member(region(u64::MAX - 2, u64::MAX - 2), b"C", b"T")];
        let locus = closed_locus(region(u64::MAX - 3, u64::MAX), &[&whole, &snp]);

        let reference = LocusReferenceBases::over(&locus);

        assert_eq!(reference.bases(), b"ACGT");
        assert_eq!(project_only(&reference, &snp[0]), b"ATGT");
        assert_eq!(
            project_only(&reference, &whole[0]),
            b"A",
            "it spans the locus, so nothing is padded onto it",
        );
    }

    /// A locus whose members leave one of its positions uncovered is refused: gathering
    /// the reference from them would leave a byte that is not a base, and every allele in
    /// the locus would carry it. The message names the position, which is what says which
    /// member stopped short.
    #[test]
    #[should_panic(expected = "leave position 13 uncovered")]
    fn a_locus_its_members_do_not_cover_is_refused() {
        let short_member = [member(region(10, 12), b"ACG", b"AAG")];
        let locus = closed_locus(region(10, 14), &[&short_member]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// A member starting before the locus is refused. The subtraction that places it
    /// inside the gathered span would wrap in the release profile, where overflow checks
    /// are off, and index the reference from somewhere near `usize::MAX`.
    #[test]
    #[should_panic(expected = "starts before the locus")]
    fn a_member_starting_before_the_locus_is_refused() {
        let early_member = [member(region(8, 14), b"AAACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&early_member]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// A member reaching past the locus's last base is refused by name rather than by
    /// running off the end of the gathered span — the same mistake from the other side,
    /// and the message is what says which member was wrong.
    #[test]
    #[should_panic(expected = "reaches past the locus")]
    fn a_member_reaching_past_the_locus_is_refused() {
        let overhanging_member = [member(region(12, 20), b"GTAAAAAAA", b"G")];
        let locus = closed_locus(region(10, 14), &[&overhanging_member]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// A member whose reference bases do not cover its own region is refused. Nothing
    /// upstream enforces it — `SampleLocusObservations` has public fields — and the
    /// consequence is a projection padded from the wrong bases, which is a plausible
    /// allele rather than a crash.
    #[test]
    #[should_panic(expected = "reference bases for a")]
    fn a_member_whose_reference_does_not_cover_its_region_is_refused() {
        let malformed_member = [SampleLocusObservations {
            reference_bases: Box::from(&b"AC"[..]),
            ..member(region(10, 14), b"ACGTA", b"A")
        }];
        let locus = closed_locus(region(10, 14), &[&malformed_member]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// A member on another contig is refused by the gather: positions alone would place it
    /// inside the locus, and every position exists on every contig.
    #[test]
    #[should_panic(expected = "on another contig")]
    fn a_member_on_another_contig_is_refused() {
        let member_elsewhere = [member(region_on(1, 10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region_on(0, 10, 14), &[&member_elsewhere]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// **Two members that overlap and disagree about the reference are refused.** They can
    /// only disagree if the samples were called against different references; copying the
    /// second over the first would leave a reference that is a mixture, and every allele in
    /// the locus would be padded from it — which looks like a sequence, not like a bug.
    #[test]
    #[should_panic(expected = "disagree on the reference at position 11")]
    fn two_members_disagreeing_on_the_reference_where_they_overlap_are_refused() {
        let first = [member(region(10, 12), b"ACG", b"AAG")];
        let disagreeing = [member(region(11, 14), b"TTTA", b"TTTT")];
        let locus = closed_locus(region(10, 14), &[&first, &disagreeing]);
        let _ = LocusReferenceBases::over(&locus);
    }

    /// **A locus the caller refused is never assembled** (spec §3.2). Closing is uncapped,
    /// so a failed locus can be thousands of bases wide; gathering one would allocate its
    /// whole width for a reference nothing will be built from.
    #[test]
    #[should_panic(expected = "only a locus the caller undertakes to build")]
    fn a_failed_locus_is_not_gathered() {
        let wide = [member(region(10, 14), b"ACGTA", b"A")];
        let failed = ClosedLocus {
            verdict: Verdict::Failed,
            ..closed_locus(region(10, 14), &[&wide])
        };
        let _ = LocusReferenceBases::over(&failed);
    }

    /// **The projection's own guards are not the gather's.** A caller reaches
    /// [`LocusReferenceBases::placing`] with a member it chose, and positions alone would
    /// place a member from another contig inside the locus — the result would be an allele
    /// padded from the wrong chromosome, with nothing to say so.
    #[test]
    #[should_panic(expected = "on another contig")]
    fn projecting_a_member_from_another_contig_is_refused() {
        let here = [member(region_on(0, 10, 14), b"ACGTA", b"A")];
        let there = member(region_on(1, 12, 12), b"G", b"T");
        let locus = closed_locus(region_on(0, 10, 14), &[&here]);
        let reference = LocusReferenceBases::over(&locus);

        let _ = reference.placing(&there);
    }

    /// A member reaching past the locus is refused by the projection too, and by name: the
    /// suffix slice would otherwise run off the end of the gathered reference and say
    /// nothing about which member was wrong.
    #[test]
    #[should_panic(expected = "reaches past the locus")]
    fn projecting_a_member_reaching_past_the_locus_is_refused() {
        let here = [member(region(10, 14), b"ACGTA", b"A")];
        let over_the_end = member(region(12, 20), b"GTAAAAAAA", b"G");
        let locus = closed_locus(region(10, 14), &[&here]);
        let reference = LocusReferenceBases::over(&locus);

        let _ = reference.placing(&over_the_end);
    }

    /// A member starting before the locus is refused by the projection too — the wrapping
    /// subtraction again, on the path a caller reaches directly.
    #[test]
    #[should_panic(expected = "starts before the locus")]
    fn projecting_a_member_starting_before_the_locus_is_refused() {
        let here = [member(region(10, 14), b"ACGTA", b"A")];
        let early = member(region(8, 14), b"AAACGTA", b"A");
        let locus = closed_locus(region(10, 14), &[&here]);
        let reference = LocusReferenceBases::over(&locus);

        let _ = reference.placing(&early);
    }

    /// **A partial sequence is not projectable, and the iterator is what keeps it out.**
    /// Its bases stop where the read's witness stopped, so padding them from the locus's
    /// reference would report reference bases over ground the read never saw. The member
    /// here carries one complete and one partial sequence, and only the complete one comes
    /// back.
    #[test]
    fn a_partial_sequence_is_not_offered_for_projection() {
        let mut with_a_partial = member(region(10, 14), b"ACGTA", b"A");
        let partial = SequenceObservation {
            read_witness: ReadWitness::from_left(2, with_a_partial.locus_len())
                .expect("a two-base run inside a five-base locus"),
            ..sequence(b"AC")
        };
        with_a_partial.observations.push(partial);

        let members = [with_a_partial];
        let locus = closed_locus(region(10, 14), &[&members]);
        let reference = LocusReferenceBases::over(&locus);
        let placement = reference.placing(&members[0]);

        let offered: Vec<&[u8]> = placement
            .projectable_sequences()
            .map(|sequence| &*sequence.bases)
            .collect();
        assert_eq!(offered, vec![&b"A"[..]], "the partial is not offered");
    }

    /// And projecting a partial anyway is refused rather than padded — the backstop behind
    /// the iterator above, for a caller that reaches into `observations` itself.
    #[test]
    #[should_panic(expected = "only a complete observation can be projected")]
    fn projecting_a_partial_sequence_is_refused() {
        let members = [member(region(10, 14), b"ACGTA", b"A")];
        let partial = SequenceObservation {
            read_witness: ReadWitness::from_left(2, members[0].locus_len())
                .expect("a two-base run inside a five-base locus"),
            ..sequence(b"AC")
        };
        let locus = closed_locus(region(10, 14), &[&members]);
        let reference = LocusReferenceBases::over(&locus);

        let mut buffer = Vec::new();
        reference
            .placing(&members[0])
            .project_into(&partial, &mut buffer);
    }

    // ---------------------------------------------------------------
    // Unification into one allele table (spec §4.2), and the owner's ruling of
    // 2026-08-17 on what a sample showed when a locus spans several of its records.
    // ---------------------------------------------------------------

    /// **One of a sample's records, with the reads behind each sequence named** — the
    /// fixture the per-read rule needs, where [`member`] above gives one unnamed sequence.
    ///
    /// Each entry of `shown` is one distinct sequence and the ids of the reads that showed
    /// it, and `num_obs` is set to how many those are — one id per read here, where real
    /// data collapses a read *pair* onto one id.
    fn record(
        region: GenomeRegion,
        reference_bases: &[u8],
        shown: &[(&[u8], &[ChainId])],
    ) -> SampleLocusObservations {
        let mut observations = member(region, reference_bases, b"");
        observations.observations = shown
            .iter()
            .map(|(bases, chains)| SequenceObservation {
                num_obs: chains.len() as u32,
                chain_ids: chains.to_vec(),
                ..sequence(bases)
            })
            .collect();
        observations
    }

    /// The table's alleles as byte strings, in table order.
    fn alleles_of(table: &AlleleTable) -> Vec<&[u8]> {
        table.alleles().iter().map(|allele| &**allele).collect()
    }

    /// **Two samples that showed the same change make one allele, not two.** This is the
    /// whole of unification, and the failure it guards is the quiet one: one variant
    /// written twice becomes two half-supported alleles, which reads as a noisy site
    /// rather than as a defect.
    #[test]
    fn two_samples_showing_the_same_change_make_one_allele() {
        let first = [member(region(12, 12), b"G", b"T")];
        let second = [member(region(12, 12), b"G", b"T")];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&first, &second, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"ACTTA", b"A"]);
    }

    /// **One deletion recorded at three different placements inside one locus comes out as
    /// one allele** — the plan's B2 fixture. A `CCC` run at 11–13 with one copy deleted,
    /// minted by three samples with the anchor base at 10, 11 and 12, so that the three
    /// records delete positions 11, 12 and 13 respectively.
    ///
    /// **What unifies them is the projection, not left-alignment.** Widening each record to
    /// the locus's six bases gives `ACCAG` whichever `C` was dropped, so placements *inside
    /// one locus* unify by construction. What left-alignment upstream buys is one step
    /// earlier and this test cannot show it: it is that the records overlap and therefore
    /// chain into one locus at all — see
    /// `two_placements_too_far_apart_to_chain_never_meet_in_one_table`, which is the
    /// failure that survives.
    #[test]
    fn three_placements_of_one_deletion_inside_the_locus_unify() {
        let anchored_at_10 = [member(region(10, 11), b"AC", b"A")];
        let anchored_at_11 = [member(region(11, 12), b"CC", b"C")];
        let anchored_at_12 = [member(region(12, 13), b"CC", b"C")];
        let spanning = [member(region(10, 15), b"ACCCAG", b"ACCCAG")];
        let locus = closed_locus(
            region(10, 15),
            &[&anchored_at_10, &anchored_at_11, &anchored_at_12, &spanning],
        );

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACCCAG"[..], b"ACCAG"],
            "the reference, and one deletion however it was anchored",
        );
    }

    /// **What left-alignment upstream really buys: two placements that do not overlap close
    /// as two loci, and no table ever sees them together.** The same one-base deletion in
    /// the same `CC` run, written by one sample at 11–12 and by another at 40–41 — the
    /// distance a caller that did not canonicalise placement could produce.
    ///
    /// The walk closes them as two loci because the records share no ground, so each table
    /// holds one sample's half of the evidence: exactly one alternative allele each, with
    /// the other sample's projection absent. Unification cannot repair that, which is why
    /// `LeftAlignPreparer` runs before the generator mints anything.
    #[test]
    fn two_placements_too_far_apart_to_chain_never_meet_in_one_table() {
        let here = [member(region(11, 12), b"CC", b"C")];
        let far_away = [member(region(40, 41), b"CC", b"C")];

        let loci: Vec<_> = LocusCloser::over(
            &[&here, &far_away],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs(std::num::NonZeroU32::new(1).expect("1 is non-zero")),
        )
        .collect();
        assert_eq!(loci.len(), 2, "the two placements never chain");

        let tables: Vec<Vec<Vec<u8>>> = loci
            .iter()
            .map(|locus| {
                AlleleTable::over(locus)
                    .alleles()
                    .iter()
                    .map(|allele| allele.to_vec())
                    .collect()
            })
            .collect();
        assert_eq!(
            tables,
            vec![
                vec![b"CC".to_vec(), b"C".to_vec()],
                vec![b"CC".to_vec(), b"C".to_vec()]
            ],
            "one deletion, in two tables, each backed by one sample instead of two",
        );
    }

    /// **The reference is allele 0 even where no sample's reads showed it** — a cohort in
    /// which every sample is homozygous for the variant still has to be genotyped against
    /// the reference, so it is seeded before any sample is looked at.
    #[test]
    fn the_reference_is_the_first_allele_even_when_no_sample_showed_it() {
        let one = [member(region(12, 12), b"G", b"T")];
        let another = [member(region(12, 12), b"G", b"A")];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&one, &another, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ACTTA", b"ACATA", b"A"],
        );
        assert_eq!(table.index_of(b"ACGTA"), Some(0));
    }

    /// A sample whose reads matched the reference adds no allele: its projection **is** the
    /// reference, so it lands on allele 0 (spec §4.2 — its depth is present, and that is
    /// what the next step reads it for).
    #[test]
    fn a_member_that_matched_the_reference_adds_no_allele() {
        let matched = [member(region(12, 12), b"G", b"G")];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&matched, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"A"]);
    }

    /// **A sample with two records inside one locus, and three reads that stand in
    /// different relations to them** — the fixture the owner's ruling is about, used by the
    /// four tests below.
    ///
    /// The locus is 10–14, `ACGTA`, held open by another sample's deletion. This sample
    /// minted two records: at 12, where its reads showed `A` (reads 7 and 9) or the
    /// reference `G` (read 11); and at 14, where they showed `C` (reads 7 and 11). So:
    ///
    /// - **read 7** is named at both records and showed `A` then `C`;
    /// - **read 11** is named at both and showed the reference then `C`;
    /// - **read 9** is named at 12 only — it either stopped before 14 or was capped there,
    ///   and nothing recorded which.
    fn a_sample_with_two_records() -> [SampleLocusObservations; 2] {
        [
            record(region(12, 12), b"G", &[(b"A", &[7, 9]), (b"G", &[11])]),
            record(region(14, 14), b"A", &[(b"C", &[7, 11])]),
        ]
    }

    /// **A read named at every one of the sample's records carries what it showed at each
    /// into one allele** — the ruling's first branch. Read 7 showed `A` at 12 and `C` at
    /// 14, so its allele is `ACATC`: the two changes, with the reference base at 13 between
    /// them, where this sample minted nothing because none of its reads departed there.
    #[test]
    fn a_read_named_at_every_record_composes_one_allele_across_them() {
        let two_records = a_sample_with_two_records();
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ACATC", b"ACGTC", b"A"],
            "read 7's two changes, read 11's one, and the other sample's deletion",
        );
    }

    /// **A read missing from one of the sample's records is removed as evidence** — the
    /// ruling's second branch, and the assertion that separates it from the superseded
    /// proposal.
    ///
    /// Read 9 showed `A` at 12 and is not named at 14. Crediting it with the reference
    /// there would put `ACATA` in the table — which is what projecting the record at 12 on
    /// its own gives, and what the sample's evidence would look like if presence were not
    /// consulted. It is absent, and the only allele carrying `A` at position 12 is read 7's
    /// compound.
    #[test]
    fn a_read_missing_from_one_record_is_removed_as_evidence() {
        let two_records = a_sample_with_two_records();
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            table.index_of(b"ACATA"),
            None,
            "read 9 covered 12 and nothing says it covered 14",
        );
        assert_eq!(
            table.index_of(b"ACATC"),
            Some(1),
            "and read 7, which was named at both, is still evidence",
        );
    }

    /// **The reads removed as evidence are counted, not merely dropped.** A removal is lost
    /// depth and nothing downstream can recover it: the read contributes no bytes, so a
    /// locus where many were removed would read as a quiet site with shallow samples.
    ///
    /// One read is removed in the fixture — read 9, named at 12 and not at 14 — and none in
    /// the deletion's sample, whose one record puts it on the branch that consults no read.
    #[test]
    fn the_reads_removed_as_evidence_are_counted() {
        let two_records = a_sample_with_two_records();
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        assert_eq!(AlleleTable::over(&locus).reads_removed_as_evidence(), 1);

        let nothing_removed = closed_locus(region(10, 14), &[&deletion]);
        assert_eq!(
            AlleleTable::over(&nothing_removed).reads_removed_as_evidence(),
            0,
            "a sample with one record has no read to remove",
        );
    }

    /// **A read is placed by the ids on a reference-matching observation too, and that is
    /// what B0 bought.** Read 11 agreed with the reference at 12 and showed `C` at 14; its
    /// allele `ACGTC` exists only because the record at 12 names the reads behind its
    /// reference sequence. Before B0 it carried none, so read 11 would have looked exactly
    /// like read 9 — absent — and this allele would have gone with it.
    #[test]
    fn a_read_that_agreed_with_the_reference_at_one_record_is_still_placed_by_it() {
        let two_records = a_sample_with_two_records();
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        assert_eq!(
            AlleleTable::over(&locus).index_of(b"ACGTC"),
            Some(2),
            "the reference at 12, the change at 14",
        );
    }

    /// **A read a depth cap dropped at one record is indistinguishable here from one that
    /// never covered it, and that is the point.** The cap acts per position and leaves no
    /// identities behind — `reads_discarded_by_cap` is a count, not a list — so what the
    /// read showed there was never recorded and the two arrive as the same fact: a read
    /// named at 12 and not at 14.
    ///
    /// The assertion is the sameness, so that a change which started consulting the count
    /// could not do it quietly: the table with a capped read at 14 is the table without one.
    #[test]
    fn a_read_a_depth_cap_dropped_at_one_record_is_removed_like_one_that_was_not_there() {
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let mut with_a_cap = a_sample_with_two_records();
        with_a_cap[1].reads_discarded_by_cap = 1;
        let without_a_cap = a_sample_with_two_records();
        assert_eq!(without_a_cap[1].reads_discarded_by_cap, 0, "the comparison");

        let capped = AlleleTable::over(&closed_locus(region(10, 14), &[&with_a_cap, &deletion]));
        let uncapped =
            AlleleTable::over(&closed_locus(region(10, 14), &[&without_a_cap, &deletion]));

        assert_eq!(capped.index_of(b"ACATA"), None, "read 9 is not evidence");
        assert_eq!(
            alleles_of(&capped),
            alleles_of(&uncapped),
            "nothing here consults the cap's count",
        );
    }

    /// **A read that was partial at one of the records is removed too.** Its bases stop
    /// where its witness stopped, so what it showed over the rest of that record was never
    /// recorded — the same fact as never having been there, and the same branch.
    ///
    /// Read 7 is complete at 12 and partial at 14; the compound `ACATC` it would otherwise
    /// carry is absent, and only read 11's allele survives.
    #[test]
    fn a_read_partial_at_one_record_is_removed_as_evidence() {
        let mut two_records = a_sample_with_two_records();
        let locus_len = two_records[1].locus_len();
        two_records[1].observations = vec![
            SequenceObservation {
                num_obs: 1,
                chain_ids: vec![11],
                ..sequence(b"C")
            },
            SequenceObservation {
                num_obs: 1,
                chain_ids: vec![7],
                read_witness: ReadWitness::from_left(1, locus_len)
                    .expect("a one-base run inside a one-base record"),
                ..sequence(b"C")
            },
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            table.index_of(b"ACATC"),
            None,
            "read 7 is not evidence here"
        );
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"ACGTC", b"A"]);
    }

    /// **A chain id means nothing across samples.** The identifier space is per file
    /// (`chain_id_allocator.rs`), so read 7 of one sample and read 7 of another are
    /// unrelated reads, and linking them would build a haplotype out of two different
    /// plants.
    ///
    /// Here the first sample has two records and read 7 is named at one of them only, while
    /// another sample's record at 14 also names a read 7. The first sample contributes
    /// nothing: neither the compound `ACATC` nor the single-change `ACATA` is in the table.
    #[test]
    fn the_same_chain_id_in_another_sample_does_not_complete_a_reads_evidence() {
        let two_records = [
            record(region(12, 12), b"G", &[(b"A", &[7])]),
            record(region(14, 14), b"A", &[(b"A", &[9])]),
        ];
        let elsewhere = [record(region(14, 14), b"A", &[(b"C", &[7])])];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &elsewhere, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(table.index_of(b"ACATC"), None);
        assert_eq!(table.index_of(b"ACATA"), None);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ACGTC", b"A"],
            "only the other sample's own change, and the deletion",
        );
    }

    /// Two reads that showed the same pair of changes propose the same allele, and the
    /// table holds it once — a composed allele goes through the same unification as
    /// everything else.
    #[test]
    fn two_reads_showing_the_same_pair_of_changes_make_one_allele() {
        let two_records = [
            record(region(12, 12), b"G", &[(b"A", &[7, 9])]),
            record(region(14, 14), b"A", &[(b"C", &[7, 9])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"ACATC", b"A"]);
    }

    /// A read named at three of the sample's records carries all three changes: the
    /// composition walks the records in coordinate order rather than handling a pair.
    #[test]
    fn a_read_named_at_three_records_carries_all_three_changes() {
        let three_records = [
            record(region(11, 11), b"C", &[(b"G", &[7])]),
            record(region(12, 12), b"G", &[(b"A", &[7])]),
            record(region(14, 14), b"A", &[(b"C", &[7])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&three_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"AGATC", b"A"]);
    }

    /// **A read that showed two different things at one record is removed**, the same as a
    /// read that was not there: there is no one thing it showed, so nothing can be
    /// composed. It is a fragment whose mates overlap and disagree — the mint keeps both as
    /// observations under the one chain id it gave the fragment
    /// (`resolve_mate_overlap_at_pos`).
    ///
    /// Read 7 is named at both of the record-12 sequences and at neither of record 14's, so
    /// it is named as often as the sample has records; only looking at *which* records
    /// those are refuses it. Read 11 is unaffected.
    #[test]
    fn a_read_that_showed_two_things_at_one_record_is_removed() {
        let two_records = [
            record(region(12, 12), b"G", &[(b"A", &[7]), (b"G", &[7, 11])]),
            record(region(14, 14), b"A", &[(b"C", &[11])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            table.index_of(b"ACATA"),
            None,
            "read 7 showed A and G at 12"
        );
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ACGTC", b"A"],
            "only read 11's allele, and the other sample's deletion",
        );
    }

    /// **The same read, at a sample's sole record, keeps both of the things it showed** —
    /// the one shape on which the two branches of the derivation differ, and it is a
    /// difference in kind rather than an optimisation.
    ///
    /// The fixture is the one above with the second record taken away: read 7 shows `A` and
    /// `G` at 12, its two mates having overlapped and disagreed there. With records at 12
    /// *and* 14 the fragment is removed, because composing across records needs it as a
    /// unit and it showed no one thing (`a_read_that_showed_two_things_at_one_record_is_removed`,
    /// where `ACATA` is absent). With 12 alone there is nothing to compose across, each
    /// mate's sequence already spans the locus, and `ACATA` is an allele.
    #[test]
    fn a_read_showing_two_things_at_a_sole_record_still_contributes_both() {
        let sole_record = [record(
            region(12, 12),
            b"G",
            &[(b"A", &[7]), (b"G", &[7, 11])],
        )];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&sole_record, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            table.index_of(b"ACATA"),
            Some(1),
            "with one record the disagreeing fragment is not removed",
        );
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"ACATA", b"A"]);
    }

    /// **The reference after a sample's last record closes the allele.** Every fixture
    /// above has the sample's last record at the locus's own end, where nothing is left to
    /// append; here the two records sit at 11 and 12 inside a locus reaching 14, so read
    /// 7's allele has to carry `TA` after the changes it showed.
    #[test]
    fn an_allele_composed_across_records_carries_the_reference_after_the_last_of_them() {
        let two_records = [
            record(region(11, 11), b"C", &[(b"G", &[7])]),
            record(region(12, 12), b"G", &[(b"A", &[7])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"AGATA", b"A"],
            "the anchor base, the two changes, and the reference to the locus's end",
        );
    }

    /// **An allele composed across records is measured in its own bases, and the reference
    /// it is closed with in the locus's** — the distinction every other multi-record fixture
    /// here hides, because a substitution is as wide as the reference it replaces, so the
    /// two counts stay equal at every step.
    ///
    /// The locus reads `ACGTA` at 10–14. One sample inserts a `G` after position 13 and has
    /// a SNP at 11, so its read's allele is **six** bases over a five-base locus; another
    /// deletes position 12 and has a SNP at 14, so its read's allele is **four**. Closing
    /// either from the length composed so far rather than from how much reference it has
    /// consumed loses the last base of the first and gains a spurious one on the second.
    #[test]
    fn an_allele_composed_across_records_closes_on_the_reference_it_consumed_not_its_length() {
        let with_an_insertion = [
            record(region(11, 11), b"C", &[(b"T", &[7])]),
            record(region(13, 13), b"T", &[(b"TG", &[7])]),
        ];
        let with_a_deletion = [
            record(region(11, 12), b"CG", &[(b"C", &[5])]),
            record(region(14, 14), b"A", &[(b"T", &[5])]),
        ];
        let spanning = [member(region(10, 14), b"ACGTA", b"ACGTA")];
        let locus = closed_locus(
            region(10, 14),
            &[&with_an_insertion, &with_a_deletion, &spanning],
        );

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ATGTGA", b"ACTT"],
            "six bases over five, and four over five",
        );
    }

    /// **Two samples that each hold several records are derived independently, and they
    /// carry the same read id.** That is not a contrived collision: the id space is per
    /// file, so every sample's ids start again from the same numbers, and two samples
    /// sharing id 7 is the ordinary case rather than the corner.
    ///
    /// Every earlier fixture has at most one sample with several records, so what one
    /// sample left in the working buffer could not be seen. Here it can: sample 0's read 7
    /// is sighted at two records and sample 1's read 7 at two more, so a buffer carried
    /// across samples makes sample 1's read look sighted four times over two records —
    /// present at more records than the sample has — and its allele `ACTTC` is removed as
    /// evidence instead of built.
    #[test]
    fn two_samples_that_each_hold_several_records_are_derived_independently() {
        let first = [
            record(region(11, 11), b"C", &[(b"G", &[7])]),
            record(region(12, 12), b"G", &[(b"A", &[7])]),
        ];
        let second = [
            record(region(12, 12), b"G", &[(b"T", &[7])]),
            record(region(14, 14), b"A", &[(b"C", &[7])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&first, &second, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"AGATA", b"ACTTC", b"A"],
            "both samples' compositions, neither lost to the other",
        );
        assert_eq!(
            table.reads_removed_as_evidence(),
            0,
            "neither read is short of its own sample's records",
        );
    }

    /// **The alleles come out in chain-id order, whatever order the reads were folded in.**
    /// The same two reads as the fixture above, with the ids listed the other way round in
    /// every observation: the table is identical, which is the determinism spec §9 rests on
    /// — the output cannot depend on which builder or which walk produced the ids.
    #[test]
    fn the_order_of_the_alleles_does_not_depend_on_the_order_the_ids_were_folded_in() {
        let as_folded = [
            record(region(12, 12), b"G", &[(b"A", &[7]), (b"G", &[11])]),
            record(region(14, 14), b"A", &[(b"C", &[7, 11])]),
        ];
        let folded_the_other_way = [
            record(region(12, 12), b"G", &[(b"G", &[11]), (b"A", &[7])]),
            record(region(14, 14), b"A", &[(b"C", &[11, 7])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];

        let one = AlleleTable::over(&closed_locus(region(10, 14), &[&as_folded, &deletion]));
        let other = AlleleTable::over(&closed_locus(
            region(10, 14),
            &[&folded_the_other_way, &deletion],
        ));

        assert_eq!(alleles_of(&one), alleles_of(&other));
        assert_eq!(
            alleles_of(&one),
            vec![&b"ACGTA"[..], b"ACATC", b"ACGTC", b"A"],
            "read 7 before read 11, by id and not by arrival",
        );
    }

    /// **A sample with one record contributes its sequences without any read being
    /// consulted** — which is the only answer available on the STR path, where the mint
    /// records no chain ids at all and needs none: an STR locus is one record, so the
    /// witness already says whether a read spanned it.
    ///
    /// Every sequence here carries an empty `chain_ids`, as an STR record's does, and both
    /// alleles come out.
    #[test]
    fn a_single_record_samples_alleles_need_no_chain_ids() {
        let str_record = [member(region(10, 14), b"ACGTA", b"ACGGTA")];
        let another = [member(region(10, 14), b"ACGTA", b"ACTA")];
        assert!(
            str_record[0].observations[0].chain_ids.is_empty(),
            "the fixture is the STR path's shape",
        );
        let locus = closed_locus(region(10, 14), &[&str_record, &another]);

        let table = AlleleTable::over(&locus);
        assert_eq!(alleles_of(&table), vec![&b"ACGTA"[..], b"ACGGTA", b"ACTA"]);
    }

    /// **An observation with reads and no chain id is refused where a locus spans several
    /// of a sample's records**, because there is no way to say whether those reads reached
    /// the sample's other records — the state the owner's ruling says must never be
    /// reached. Dropping them instead would be silent: the locus would simply come back
    /// with fewer alleles.
    #[test]
    #[should_panic(expected = "carries 3 reads and no chain id")]
    fn an_observation_with_reads_and_no_chain_id_is_refused_across_records() {
        let two_records = [
            record(region(12, 12), b"G", &[(b"A", &[7])]),
            member(region(14, 14), b"A", b"C"),
        ];
        assert_eq!(
            two_records[1].observations[0].num_obs, 3,
            "the fixture's unnamed reads",
        );
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let _ = AlleleTable::over(&locus);
    }

    /// **Two of one sample's records overlapping are refused whether or not a read is named
    /// at both.** The mint cannot produce them — an event inside an open record's footprint
    /// folds into that record — but the refusal has to be structural, because the
    /// composition's own backstop is only reached by a read named at both records: here the
    /// two records carry *different* reads, so every read fails the presence test, and
    /// without this check the sample would contribute nothing at all and say nothing about
    /// why.
    #[test]
    #[should_panic(expected = "are not disjoint and in coordinate order")]
    fn two_of_one_samples_records_overlapping_are_refused_even_with_no_shared_read() {
        let overlapping = [
            record(region(10, 12), b"ACG", &[(b"AAG", &[7])]),
            record(region(12, 14), b"GTA", &[(b"TTA", &[9])]),
        ];
        let locus = closed_locus(region(10, 14), &[&overlapping]);

        let _ = AlleleTable::over(&locus);
    }

    /// And with a read named at both, which is the case the composition itself would catch
    /// — the structural check is what makes the two arrive at the same refusal.
    #[test]
    #[should_panic(expected = "are not disjoint and in coordinate order")]
    fn two_of_one_samples_records_overlapping_are_refused() {
        let overlapping = [
            record(region(10, 12), b"ACG", &[(b"AAG", &[7])]),
            record(region(12, 14), b"GTA", &[(b"TTA", &[7])]),
        ];
        let locus = closed_locus(region(10, 14), &[&overlapping]);

        let _ = AlleleTable::over(&locus);
    }

    /// **A sample's records arriving out of coordinate order are refused, and named as
    /// that** rather than as an overlap. They do not overlap — 12 and 14 are a base apart —
    /// and a composition walking them in the order given would write the later one first.
    /// `SampleMembers` promises coordinate order and its fields are public, so this is the
    /// caller's-mistake class the file's other checks are written for.
    #[test]
    #[should_panic(expected = "are not disjoint and in coordinate order")]
    fn a_samples_records_out_of_coordinate_order_are_refused() {
        let out_of_order = [
            record(region(14, 14), b"A", &[(b"C", &[7])]),
            record(region(12, 12), b"G", &[(b"A", &[7])]),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&out_of_order, &deletion]);

        let _ = AlleleTable::over(&locus);
    }

    /// Every allele the samples showed is findable by its bytes, and nothing else is: the
    /// lookup is what lets B3 attribute a sample's support by composing its reads again
    /// rather than the build carrying an assignment for every one.
    #[test]
    fn index_of_finds_every_allele_and_nothing_else() {
        let snp = [member(region(12, 12), b"G", b"T")];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&snp, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(table.index_of(b"ACGTA"), Some(0));
        assert_eq!(table.index_of(b"ACTTA"), Some(1));
        assert_eq!(table.index_of(b"A"), Some(2));
        assert_eq!(
            table.index_of(b"ACGTT"),
            None,
            "a sequence no sample showed here",
        );
        assert_eq!(table.reference().bases(), b"ACGTA");
    }

    // ---------------------------------------------------------------
    // Per-sample support against the table (arch §4), and the division the owner ruled on.
    // ---------------------------------------------------------------

    /// One record of a sample, with its sequences spelled out in full — the fixture the
    /// support tests need, where [`record`] above only sets bases and read ids.
    fn record_of(
        region: GenomeRegion,
        reference_bases: &[u8],
        sequences: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        let mut observations = member(region, reference_bases, b"");
        observations.observations = sequences;
        observations
    }

    /// A sequence with its support spelled out: `reads` names the reads that showed it, and
    /// the five sums are the mint's own numbers over exactly those reads.
    fn shown_by(
        bases: &[u8],
        reads: &[ChainId],
        q_sum: f64,
        num_fwd: u32,
        mapq_sum: u32,
        mapq_sum_sq: u64,
        placed_left: u32,
    ) -> SequenceObservation {
        SequenceObservation {
            num_obs: reads.len() as u32,
            chain_ids: reads.to_vec(),
            q_sum,
            num_fwd,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
            ..sequence(bases)
        }
    }

    /// **At a sample with one record every sum is the mint's own** — nothing is divided,
    /// because each of that record's sequences is one allele and all the reads behind it
    /// showed that allele.
    ///
    /// The sample has four reads showing a SNP and two showing the reference; the SNP's
    /// allele carries exactly the four reads' numbers and the reference's exactly the two.
    #[test]
    fn a_one_record_samples_support_is_the_mints_own_numbers() {
        let one_record = [record_of(
            region(12, 12),
            b"G",
            vec![
                shown_by(b"T", &[7, 9, 11, 13], -8.0, 3, 240, 14_400, 1),
                shown_by(b"G", &[15, 17], -5.0, 1, 100, 5_000, 0),
            ],
        )];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&one_record, &deletion]);

        let observed = CohortObservation::over(&locus);
        let support = &observed.per_sample[0];

        assert_eq!(support.sample, 0);
        assert_eq!(support.reads_composed_across_records, 0, "nothing divided");
        assert_eq!(
            support.support_for(REFERENCE_ALLELE),
            AlleleSupport {
                num_reads: 2,
                num_fwd: 1,
                q_sum: -5.0,
                mapq_sum: 100,
                mapq_sum_sq: 5_000,
                placed_left: 0,
            },
            "the reads that agreed with the reference are support for allele 0",
        );
        let snp = observed
            .alleles
            .iter()
            .position(|allele| &**allele == b"ACTTA")
            .expect("the SNP is in the table");
        assert_eq!(
            support.support_for(snp),
            AlleleSupport {
                num_reads: 4,
                num_fwd: 3,
                q_sum: -8.0,
                mapq_sum: 240,
                mapq_sum_sq: 14_400,
                placed_left: 1,
            },
        );
    }

    /// **Two of a sample's own observations that reach the same allele are summed** (spec
    /// §4.2), and support is **never** merged across alleles. The same bases from two read
    /// groups are two observations and one allele.
    #[test]
    fn support_is_summed_within_an_allele_and_never_across_them() {
        let two_groups = [record_of(
            region(12, 12),
            b"G",
            vec![
                shown_by(b"T", &[7, 9], -4.0, 2, 120, 7_200, 1),
                SequenceObservation {
                    read_group: ReadGroupId(1),
                    ..shown_by(b"T", &[11], -1.0, 0, 50, 2_500, 0)
                },
                shown_by(b"A", &[13], -3.0, 1, 60, 3_600, 0),
            ],
        )];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_groups, &deletion]);

        let observed = CohortObservation::over(&locus);
        let support = &observed.per_sample[0];
        let allele_of = |bases: &[u8]| {
            observed
                .alleles
                .iter()
                .position(|allele| &**allele == bases)
                .expect("in the table")
        };

        assert_eq!(
            support.support_for(allele_of(b"ACTTA")),
            AlleleSupport {
                num_reads: 3,
                num_fwd: 2,
                q_sum: -5.0,
                mapq_sum: 170,
                mapq_sum_sq: 9_700,
                placed_left: 1,
            },
            "the two read groups' reads and sums added",
        );
        assert_eq!(
            support.support_for(allele_of(b"ACATA")),
            AlleleSupport {
                num_reads: 1,
                num_fwd: 1,
                q_sum: -3.0,
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 0,
            },
            "the other allele keeps its own, unmixed",
        );
    }

    /// **A read composed across two records has its five sums divided, and its count is
    /// still exact** — the owner's ruling of 2026-08-17.
    ///
    /// Read 7 is named at both of the sample's records. At 12 it is one of two reads behind
    /// `A`, whose sums are `q = −4.0`, 2 forward, MAPQ 120, MAPQ² 7,200, 0 placed left; at
    /// 14 it is the only read behind `C`, with `q = −0.5`, 0 forward, MAPQ 50, MAPQ² 2,500,
    /// 1 placed left. So its share is:
    ///
    /// - **quality −0.5**, the weaker of the two means (−4.0/2 = −2.0 against −0.5/1 =
    ///   −0.5). These are `ln P(error)`, so the number nearer zero is the **worse** read:
    ///   −0.5 is an error probability of about 3 in 5 and −2.0 about 1 in 7. What read 7
    ///   saw badly at 14 is what limits its evidence for an allele covering both;
    /// - **the other four pooled over the three reads behind the two sequences**: forward
    ///   2/3, MAPQ 170/3 ≈ 56.7, MAPQ² 9,700/3 ≈ 3,233.3, placed left 1/3 — each rounded
    ///   once, at the allele.
    ///
    /// Read 9, named at 12 alone, is removed, so `ACATA` is not in the table at all.
    #[test]
    fn a_read_composed_across_records_has_its_sums_divided_and_its_count_exact() {
        let two_records = [
            record_of(
                region(12, 12),
                b"G",
                vec![shown_by(b"A", &[7, 9], -4.0, 2, 120, 7_200, 0)],
            ),
            record_of(
                region(14, 14),
                b"A",
                vec![shown_by(b"C", &[7], -0.5, 0, 50, 2_500, 1)],
            ),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let observed = CohortObservation::over(&locus);
        let support = &observed.per_sample[0];
        let compound = observed
            .alleles
            .iter()
            .position(|allele| &**allele == b"ACATC")
            .expect("read 7's allele");

        assert_eq!(
            support.support_for(compound),
            AlleleSupport {
                num_reads: 1,
                num_fwd: 1,
                q_sum: -0.5,
                mapq_sum: 57,
                mapq_sum_sq: 3_233,
                placed_left: 0,
            },
        );
        assert_eq!(support.reads_composed_across_records, 1);
        assert_eq!(support.reads_removed_as_evidence, 1, "read 9");
        assert_eq!(observed.per_sample[1].reads_composed_across_records, 0);
    }

    /// **The case the division exists for: one observation's reads split across two
    /// alleles**, with a record holding two sequences and the weakest sighting first for one
    /// read and last for the other.
    ///
    /// The sample has records at 12 and 14. At 12 its reads showed `A` (reads 7 and 9, mean
    /// quality −0.5) or `T` (read 11, mean −8.0); at 14 all three showed `C` (mean −6.0). So
    /// the three reads of the `C` observation take two different paths across the locus, and
    /// its five sums have to be divided:
    ///
    /// - **reads 7 and 9** compose `ACATC`. Each takes the weaker of −0.5 and −6.0, which is
    ///   **−0.5 and comes first** — so a rule that kept the last sighting's quality would
    ///   answer −6.0 here. Their pooled counts come from 5 reads: forward 5/5, MAPQ 300/5,
    ///   MAPQ² 18,000/5, left 0/5, each ×2 reads.
    /// - **read 11** composes `ACTTC`, taking the weaker of −8.0 and −6.0, which is **−6.0
    ///   and comes last**. Its pooled counts come from 4 reads: forward 3/4 → 1, MAPQ 210/4
    ///   = 52.5 → 53, MAPQ² 11,700/4 = 2,925, left 1/4 → 0.
    ///
    /// Reading the wrong sequence of a record — the first rather than the one the read was
    /// sighted at — gives read 11 read 7's numbers instead.
    #[test]
    fn one_observations_reads_split_across_two_alleles_and_each_takes_its_own_share() {
        let two_records = [
            record_of(
                region(12, 12),
                b"G",
                vec![
                    shown_by(b"A", &[7, 9], -1.0, 2, 120, 7_200, 0),
                    shown_by(b"T", &[11], -8.0, 0, 30, 900, 1),
                ],
            ),
            record_of(
                region(14, 14),
                b"A",
                vec![shown_by(b"C", &[7, 9, 11], -18.0, 3, 180, 10_800, 0)],
            ),
        ];
        let one_read_removed = [
            record_of(
                region(11, 11),
                b"C",
                vec![shown_by(b"G", &[21, 23], -4.0, 1, 60, 3_600, 0)],
            ),
            record_of(
                region(13, 13),
                b"T",
                vec![shown_by(b"A", &[21], -2.0, 1, 30, 900, 0)],
            ),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(
            region(10, 14),
            &[&two_records, &one_read_removed, &deletion],
        );

        let observed = CohortObservation::over(&locus);
        let allele_of = |bases: &[u8]| {
            observed
                .alleles
                .iter()
                .position(|allele| &**allele == bases)
                .unwrap_or_else(|| panic!("{} is not in the table", String::from_utf8_lossy(bases)))
        };
        let split = &observed.per_sample[0];

        assert_eq!(
            split.support_for(allele_of(b"ACATC")),
            AlleleSupport {
                num_reads: 2,
                num_fwd: 2,
                q_sum: -1.0,
                mapq_sum: 120,
                mapq_sum_sq: 7_200,
                placed_left: 0,
            },
            "reads 7 and 9, each taking the weaker −0.5 of its two sightings",
        );
        assert_eq!(
            split.support_for(allele_of(b"ACTTC")),
            AlleleSupport {
                num_reads: 1,
                num_fwd: 1,
                q_sum: -6.0,
                mapq_sum: 53,
                mapq_sum_sq: 2_925,
                placed_left: 0,
            },
            "read 11, from the sequence it was sighted at rather than the record's first",
        );
        assert_eq!(
            split.reads_composed_across_records, 3,
            "all three reads were composed, not merely one",
        );
        assert_eq!(split.reads_removed_as_evidence, 0);

        let with_a_removal = &observed.per_sample[1];
        assert_eq!(
            with_a_removal.reads_removed_as_evidence, 1,
            "read 23, named at 11 and not at 13 — and this sample's count, not the cohort's",
        );
        assert_eq!(with_a_removal.reads_composed_across_records, 1);
        assert_eq!(
            with_a_removal.support_for(allele_of(b"AGGAA")).num_reads,
            1,
            "read 21 composed across both records",
        );
    }

    /// **The weakest sighting sets the quality, and this is where the code parts company
    /// with production's.** Production's merger plan says a compound's quality "cannot
    /// exceed any single constituent's" and its code takes the `min` of the constituents'
    /// mean `q_sum` — which, these being sums of `ln P(error)`, is the constituent the read
    /// saw *best*.
    ///
    /// Here the two sightings mean −6.0 and −1.0 per read, and −1.0 is the worse of the two
    /// (an error probability of about 1 in 3, against about 1 in 400). That is what the
    /// allele carries; production's `min` would give −6.0, five nats better, making an
    /// allele spanning two records look better evidenced than either piece of it.
    #[test]
    fn the_weakest_sighting_sets_a_composed_reads_quality() {
        let two_records = [
            record_of(
                region(12, 12),
                b"G",
                vec![shown_by(b"A", &[7], -6.0, 1, 60, 3_600, 0)],
            ),
            record_of(
                region(14, 14),
                b"A",
                vec![shown_by(b"C", &[7], -1.0, 1, 60, 3_600, 0)],
            ),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let observed = CohortObservation::over(&locus);
        let compound = observed
            .alleles
            .iter()
            .position(|allele| &**allele == b"ACATC")
            .expect("read 7's allele");

        assert_eq!(
            observed.per_sample[0].support_for(compound).q_sum,
            -1.0,
            "the weaker of −6.0 and −1.0, not the better",
        );
    }

    /// **The rounding of a divided count saturates at both ends and answers 0 for a
    /// number that is not one.** None of the three is reachable from this module's own
    /// arithmetic — the divided sums are non-negative and bounded by the reads behind them —
    /// so the boundaries are asserted here rather than left to a reader to derive from
    /// `as`'s rules.
    #[test]
    fn the_rounding_of_a_divided_count_saturates_at_both_ends() {
        assert_eq!(round_to_u32(0.5), 1, "half a read rounds up");
        assert_eq!(round_to_u32(0.4999), 0);
        assert_eq!(round_to_u32(-3.0), 0, "a negative count is none");
        assert_eq!(round_to_u32(f64::from(u32::MAX) * 2.0), u32::MAX);
        assert_eq!(round_to_u32(f64::INFINITY), u32::MAX);
        assert_eq!(round_to_u32(f64::NAN), 0, "not a number is not a count");
        assert_eq!(round_to_u64(2.5), 3);
        assert_eq!(round_to_u64(-1.0), 0);
        assert_eq!(round_to_u64(f64::INFINITY), u64::MAX);
        assert_eq!(round_to_u64(f64::NAN), 0);
    }

    /// **A sequence no read is behind contributes no allele, on either branch.** Nothing
    /// supports it, so interning it would put an allele in the table carrying a quality
    /// nobody measured — and the two branches must agree, or the same input is answered one
    /// way at a one-record sample and another where the locus spans several records.
    #[test]
    fn a_sequence_with_no_reads_behind_it_contributes_no_allele_on_either_branch() {
        let sole_record = [record_of(
            region(12, 12),
            b"G",
            vec![
                shown_by(b"T", &[7], -2.0, 1, 60, 3_600, 0),
                SequenceObservation {
                    num_obs: 0,
                    chain_ids: Vec::new(),
                    q_sum: -4.0,
                    ..sequence(b"A")
                },
            ],
        )];
        let two_records = [
            record_of(
                region(11, 11),
                b"C",
                vec![
                    shown_by(b"G", &[9], -2.0, 1, 60, 3_600, 0),
                    SequenceObservation {
                        num_obs: 0,
                        chain_ids: Vec::new(),
                        q_sum: -4.0,
                        ..sequence(b"T")
                    },
                ],
            ),
            record_of(
                region(13, 13),
                b"T",
                vec![shown_by(b"C", &[9], -2.0, 1, 60, 3_600, 0)],
            ),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&sole_record, &two_records, &deletion]);

        let observed = CohortObservation::over(&locus);
        assert_eq!(
            observed.alleles.len(),
            4,
            "the reference, one allele from each sample, and the deletion — got {:?}",
            observed
                .alleles
                .iter()
                .map(|allele| String::from_utf8_lossy(allele).into_owned())
                .collect::<Vec<_>>(),
        );
        assert!(
            observed.per_sample[0]
                .supported
                .iter()
                .all(|supported| supported.support.q_sum >= -2.0),
            "no allele carries the unsupported sequence's quality",
        );
    }

    /// **The divided counts are rounded once per allele, not once per read**, and the
    /// left-placed count is where the two answers differ. Reads 7, 9 and 11 are each sighted
    /// at two sequences with six reads behind them in total, of which **three** started left
    /// of the anchor and **one** was on the forward strand. So each read's share is a half a
    /// left-placed read and a sixth of a forward one.
    ///
    /// Rounded once: 1.5 left-placed becomes **2** and 0.5 forward becomes **1**. Rounded
    /// per read, each half would become a whole and the allele would claim **3** reads
    /// started left, out of three reads — twice the truth, and a number that cannot be
    /// distinguished from every read starting left.
    #[test]
    fn the_divided_counts_are_rounded_once_per_allele() {
        let two_records = [
            record_of(
                region(12, 12),
                b"G",
                vec![shown_by(b"A", &[7, 9, 11], -6.0, 1, 90, 5_400, 3)],
            ),
            record_of(
                region(14, 14),
                b"A",
                vec![shown_by(b"C", &[7, 9, 11], -6.0, 0, 90, 5_400, 0)],
            ),
        ];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let observed = CohortObservation::over(&locus);
        let compound = observed
            .alleles
            .iter()
            .position(|allele| &**allele == b"ACATC")
            .expect("the three reads' allele");

        assert_eq!(
            observed.per_sample[0].support_for(compound),
            AlleleSupport {
                num_reads: 3,
                num_fwd: 1,
                q_sum: -6.0,
                mapq_sum: 90,
                mapq_sum_sq: 5_400,
                placed_left: 2,
            },
            "1.5 left-placed reads rounded once are 2, where three halves rounded are 3",
        );
    }

    /// **A sample that did not cover the locus has no entry**, which is a different fact
    /// from an entry whose support is all reference (spec §4.2). Sample 1 here covers
    /// nothing, and the walk is what leaves it out — so the entries name their samples
    /// rather than sitting at their index.
    #[test]
    fn a_sample_with_no_coverage_has_no_support_at_all() {
        let covering = [member(region(12, 12), b"G", b"T")];
        let nothing = [];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];

        let loci: Vec<_> = LocusCloser::over(
            &[&covering, &nothing, &deletion],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        )
        .collect();
        assert_eq!(loci.len(), 1, "one locus, held open by the deletion");

        let observed = CohortObservation::over(&loci[0]);
        assert_eq!(
            observed
                .per_sample
                .iter()
                .map(|s| s.sample)
                .collect::<Vec<_>>(),
            vec![0, 2],
            "the sample that covered nothing is absent, not zeroed",
        );
        assert_eq!(observed.region, region(10, 14));
    }

    /// **A sample lists the alleles it showed and no others**, however many the cohort ended
    /// up with. Sample 0 showed one allele and never saw sample 1's `ACATA`; asking about it
    /// answers no reads and no sums, which is what a zeroed row would have said, and the
    /// three alleles it did not show cost it nothing.
    ///
    /// This is the shape the owner chose at Checkpoint B over one row per sample per allele.
    #[test]
    fn a_sample_lists_the_alleles_it_showed_and_answers_nothing_for_the_rest() {
        let first = [member(region(12, 12), b"G", b"T")];
        let second = [member(region(12, 12), b"G", b"A")];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&first, &second, &deletion]);

        let observed = CohortObservation::over(&locus);
        assert_eq!(observed.alleles.len(), 4);
        for support in &observed.per_sample {
            assert_eq!(
                support.supported.len(),
                1,
                "sample {} showed one allele and lists one",
                support.sample,
            );
            assert!(
                support.supported[0].support.num_reads > 0,
                "an entry with no reads would be a row this shape exists to avoid",
            );
        }

        let late = observed
            .alleles
            .iter()
            .position(|allele| &**allele == b"ACATA")
            .expect("sample 1's allele");
        assert_eq!(
            observed.per_sample[0].support_for(late),
            AlleleSupport::default(),
            "an allele this sample never showed answers nothing at all",
        );
        assert_eq!(
            observed.per_sample[1].support_for(late).num_reads,
            3,
            "and the sample that did show it answers its reads",
        );
    }

    // ---------------------------------------------------------------
    // One region's worth of building (spec §6.2), and the ownership rule that makes the
    // same answer come out however the genome is divided.
    // ---------------------------------------------------------------

    /// The three verdicts through `build_region`: a locus worth building comes out as an
    /// observation, a locus too wide comes out as a span and nothing else, and a locus too
    /// quiet comes out in **neither** — ground the caller examined and found empty, where a
    /// failure is ground it refused (spec §4.3).
    ///
    /// The bound is 5 bases and the threshold 2 reads. The SNP at 12 has 3 non-reference
    /// reads over one base; the chain from 20 to 27     /// at 40 is carried by a single read, under the threshold.
    #[test]
    fn a_region_yields_the_survivors_and_the_failed_spans_and_nothing_for_the_quiet() {
        let sample = [
            member(region(12, 12), b"G", b"T"),
            member(region(20, 27), b"ACGTACGT", b"A"),
            SampleLocusObservations {
                observations: vec![SequenceObservation {
                    num_obs: 1,
                    ..sequence(b"C")
                }],
                ..member(region(40, 40), b"A", b"C")
            },
        ];

        let outcome = build_region(
            region(1, 100),
            &[&sample],
            MaxCohortLocusSpan(std::num::NonZeroU32::new(5).expect("5 is non-zero")),
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            outcome
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(12, 12)],
            "the SNP is built; the wide chain and the single-read SNP are not",
        );
        assert_eq!(
            outcome.failed_locus_spans,
            vec![region(20, 27)],
            "refused, and counted"
        );
    }

    /// **A locus belongs to the region its first position falls in.** One starting before the
    /// region is an earlier builder's, even though it reaches well inside this one — the
    /// deletion at 5–30 covers most of this region and is not ours, and neither is the SNP at
    /// 12 it swallowed, because they are one locus and that locus starts at 5.
    #[test]
    fn a_locus_starting_before_the_region_belongs_to_the_earlier_builder() {
        let deletion = [member(region(5, 30), &[b'A'; 26], b"A")];
        let snp = [member(region(12, 12), b"A", b"T")];

        let outcome = build_region(
            region(10, 50),
            &[&deletion, &snp],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert!(
            outcome.cohort_observations.is_empty() && outcome.failed_locus_spans.is_empty(),
            "the locus opens at 5, which is another region's ground",
        );
    }

    /// **A locus starting inside the region is ours whole, however far past the end it
    /// reaches.** The deletion opens at 48 and runs to 60, ten bases beyond this region, and
    /// it comes out entire — that is what keeps a locus from being cut at a boundary
    /// (spec §6.1).
    #[test]
    fn a_locus_reaching_past_the_region_is_still_built_whole() {
        let deletion = [member(region(48, 60), &[b'A'; 13], b"A")];
        let snp = [member(region(55, 55), b"A", b"T")];

        let outcome = build_region(
            region(10, 50),
            &[&deletion, &snp],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            outcome
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(48, 60)],
            "one locus, reaching ten bases past the region that owns it",
        );
        assert_eq!(
            outcome.cohort_observations[0].per_sample.len(),
            2,
            "and both samples' evidence, including the one that only covers the tail",
        );
    }

    /// A locus starting after the region's last base belongs to a later builder, and the walk
    /// stops rather than reading on.
    #[test]
    fn a_locus_starting_after_the_region_is_left_to_the_later_builder() {
        let sample = [
            member(region(12, 12), b"G", b"T"),
            member(region(80, 80), b"G", b"T"),
        ];

        let outcome = build_region(
            region(1, 50),
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            outcome
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(12, 12)],
        );
    }

    /// **Every region delivers an outcome, including one with nothing in it.** The organiser
    /// drains regions in order and sums their counts, so a region that built nothing still
    /// has to arrive (spec §6.3).
    #[test]
    fn a_region_with_nothing_in_it_still_delivers_an_outcome() {
        let nothing: [SampleLocusObservations; 0] = [];

        let outcome = build_region(
            region(1, 50),
            &[&nothing],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert!(outcome.cohort_observations.is_empty());
        assert!(outcome.failed_locus_spans.is_empty());
    }

    /// **A region is one contig's ground, and a position on another contig says nothing about
    /// where the walk is.** The earlier contig's locus sits at position 80, past this
    /// region's last base of 50 — so a rule that compared positions alone would stop the walk
    /// there and lose everything on the contig the region is actually on.
    ///
    /// The walk yields loci in contig-then-position order, so contig 0's locus at 80 is seen
    /// first and contig 1's at 5 second.
    #[test]
    fn a_locus_on_another_contig_is_not_this_regions() {
        let elsewhere = [member(region_on(0, 80, 80), b"G", b"T")];
        let here = [member(region_on(1, 5, 5), b"G", b"T")];

        let outcome = build_region(
            region_on(1, 1, 50),
            &[&elsewhere, &here],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            outcome
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region_on(1, 5, 5)],
            "contig 0's locus is skipped rather than ending the walk, and contig 1's is built",
        );
    }

    /// **The keep threshold reaches the walk from here.** A builder forwards both parameters,
    /// and only one of them was ever handed a value of its own: at three non-reference reads
    /// the locus is built at the default of two and dropped at four, which is what says the
    /// argument is passed rather than a default read inside.
    #[test]
    fn the_keep_threshold_a_builder_is_given_is_the_one_the_walk_uses() {
        let sample = [member(region(12, 12), b"G", b"T")];
        let built = build_region(
            region(1, 50),
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );
        let too_quiet = build_region(
            region(1, 50),
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs(std::num::NonZeroU32::new(4).expect("4 is non-zero")),
        );

        assert_eq!(built.cohort_observations.len(), 1, "three reads reach two");
        assert!(
            too_quiet.cohort_observations.is_empty() && too_quiet.failed_locus_spans.is_empty(),
            "and fall short of four, which drops the locus without counting it",
        );
    }

    /// **A locus opening on the region's very first base is that region's.** Both ends of
    /// the ownership rule are inclusive, and this is the end that had no fixture: every
    /// other one here opens strictly inside its region, so widening the comparison to
    /// `start <= region.start` — which makes *every* builder skip such a locus, losing it
    /// from the run with nothing to say so — passed all 104 tests.
    ///
    /// It is not a rare shape. Building regions are dealt out end to end, twenty bases wide
    /// by default, so about one locus in twenty opens on a region's first base.
    #[test]
    fn a_locus_opening_on_the_regions_first_base_belongs_to_that_region() {
        let opening_on_the_boundary = [member(region(21, 21), b"G", b"T")];
        let opening_on_the_last_base = [member(region(30, 30), b"G", b"T")];

        let outcome = build_region(
            region(21, 30),
            &[&opening_on_the_boundary, &opening_on_the_last_base],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            outcome
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(21, 21), region(30, 30)],
            "the first base and the last base are both inside the region that owns them",
        );
    }

    /// **The same loci come out however the ground is divided into regions** — the property
    /// the whole parallel arrangement rests on (spec §9, §15), asserted here on the ownership
    /// rule alone, before any builder or organiser exists to test it end to end.
    ///
    /// One walk over 1–100 against ten walks over ten bases each: the same observations, in
    /// the same order, and the same failed spans. The fixture puts a locus on a boundary
    /// deliberately — the deletion at 30–34 opens in one ten-base region and ends in the next.
    #[test]
    fn dividing_the_ground_into_regions_changes_nothing() {
        let deletions = [
            member(region(12, 12), b"G", b"T"),
            member(region(30, 34), b"ACGTA", b"A"),
            member(region(77, 77), b"G", b"T"),
        ];
        let others = [
            member(region(32, 32), b"G", b"C"),
            member(region(90, 96), b"ACGTACG", b"A"),
        ];
        let samples: [&[SampleLocusObservations]; 2] = [&deletions, &others];
        let bound = MaxCohortLocusSpan(std::num::NonZeroU32::new(5).expect("5 is non-zero"));

        let whole = build_region(region(1, 100), &samples, bound, MinAltObs::DEFAULT);

        let mut in_pieces = RegionOutcome::default();
        for tenth in 0..10u64 {
            let piece = build_region(
                region(tenth * 10 + 1, tenth * 10 + 10),
                &samples,
                bound,
                MinAltObs::DEFAULT,
            );
            in_pieces
                .cohort_observations
                .extend(piece.cohort_observations);
            in_pieces
                .failed_locus_spans
                .extend(piece.failed_locus_spans);
        }

        assert_eq!(
            whole
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            in_pieces
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            "the same loci, in the same order",
        );
        assert_eq!(whole.failed_locus_spans, in_pieces.failed_locus_spans);
        for (whole, piece) in whole
            .cohort_observations
            .iter()
            .zip(&in_pieces.cohort_observations)
        {
            assert_eq!(whole.alleles, piece.alleles, "at {}", whole.region);
            assert_eq!(whole.per_sample, piece.per_sample, "at {}", whole.region);
        }
        assert_eq!(
            whole
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(12, 12), region(30, 34), region(77, 77)],
            "the two SNPs, and the deletion that swallowed the SNP at 32 into one locus",
        );
        assert_eq!(
            whole.failed_locus_spans,
            vec![region(90, 96)],
            "seven bases against five"
        );
    }

    /// The reads that said nothing are carried through from the records rather than
    /// re-derived, and summed where a sample has several of them (arch §4).
    #[test]
    fn the_reads_without_an_observation_are_carried_through_from_the_records() {
        let mut two_records = a_sample_with_two_records();
        two_records[0].reads_without_observation = 2;
        two_records[1].reads_without_observation = 3;
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&two_records, &deletion]);

        let observed = CohortObservation::over(&locus);
        assert_eq!(observed.per_sample[0].reads_without_observation, 5);
        assert_eq!(observed.per_sample[1].reads_without_observation, 0);
    }

    /// A partial sequence contributes no allele where the sample has one record either —
    /// the projection refuses it, and the table is what a reader sees the consequence in.
    #[test]
    fn a_partial_sequence_contributes_no_allele() {
        let mut with_a_partial = member(region(12, 12), b"G", b"T");
        with_a_partial.observations.push(SequenceObservation {
            read_witness: ReadWitness::from_left(1, with_a_partial.locus_len())
                .expect("a one-base run inside a one-base record"),
            ..sequence(b"A")
        });
        let members = [with_a_partial];
        let deletion = [member(region(10, 14), b"ACGTA", b"A")];
        let locus = closed_locus(region(10, 14), &[&members, &deletion]);

        let table = AlleleTable::over(&locus);
        assert_eq!(
            alleles_of(&table),
            vec![&b"ACGTA"[..], b"ACTTA", b"A"],
            "the partial's `A` would have been a fourth allele",
        );
    }
}
