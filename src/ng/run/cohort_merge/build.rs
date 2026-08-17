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

use super::close::{ClosedLocus, SampleMembers, Verdict, span_of};
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
        let reference = LocusReferenceBases::over(locus);
        let mut alleles = AlleleLookup::default();
        let reference_allele = alleles.intern(reference.bases());
        assert_eq!(
            reference_allele, REFERENCE_ALLELE,
            "the reference must be the locus's first allele",
        );

        let mut scratch = ReadAlleleScratch::default();
        let mut reads_removed_as_evidence = 0u32;
        for sample_members in &locus.members {
            let removed = alleles_of_sample(&reference, sample_members, &mut scratch, |bases| {
                alleles.intern(bases);
            });
            reads_removed_as_evidence = reads_removed_as_evidence.saturating_add(removed);
        }

        Self {
            reference,
            alleles,
            reads_removed_as_evidence,
        }
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
///   are composed in coordinate order into one allele. The ground between two records —
///   where this sample minted nothing, because none of its reads departed from the
///   reference there — is filled from the locus's reference, which is what "this sample
///   had nothing to say here" means.
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
    mut emit: impl FnMut(&[u8]),
) -> u32 {
    let ReadAlleleScratch { by_read, composed } = scratch;
    let records = members.observations;
    by_read.clear();

    if let [only_record] = records {
        let placement = reference.placing(only_record);
        for sequence in placement.projectable_sequences() {
            placement.project_into(sequence, composed);
            emit(composed);
        }
        return 0;
    }

    for (record_index, record) in records.iter().enumerate() {
        for (sequence_index, sequence) in record.observations.iter().enumerate() {
            if sequence.read_witness != ReadWitness::Complete {
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
            // The `num_obs == 0` escape is for an observation that names no reads because
            // it has none: there is nothing to place, so there is nothing to complain
            // about. No producer emits one today — both paths derive an observation from
            // the reads that showed it — so the clause guards a state rather than a caller.
            //
            // Release-level, like the walk's own checks beside it, and one comparison per
            // sequence against a derivation that already copies every sequence's bases.
            // **When observations are decoded from a psp file this becomes corrupt input
            // and must become a `RunError`** (arch §5), as the reference-width check
            // above must.
            assert!(
                sequence.num_obs == 0 || !sequence.chain_ids.is_empty(),
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
        emit(composed);
    }

    removed
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
