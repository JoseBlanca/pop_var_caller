//! Assembling a cohort locus: widening every member's observed sequence to the whole
//! locus span.
//!
//! A cohort locus is one stretch of genome, and the samples in it did not all record
//! the same stretch: one sample's deletion covers five bases where another's SNP covers
//! one. Before the two can be compared they have to be written over the same ground —
//! each sequence padded with the reference bases on either side until it spans the whole
//! locus (`doc/devel/ng/spec/cohort_merge.md` §4.2, projection). That is what this file
//! does; unifying the projections into one allele table is the next step's
//! (`doc/devel/ng/impl_plan/cohort_merge.md` B2).
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

use super::close::{ClosedLocus, Verdict, span_of};
use crate::ng::locus_generation::{ReadWitness, SampleLocusObservations, SequenceObservation};
use crate::ng::types::GenomeRegion;

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
            for member in sample_members.observations {
                assert_eq!(
                    member.region.contig, locus.region.contig,
                    "member {} of sample index {} is on another contig from the locus {}",
                    member.region, sample_members.sample, locus.region,
                );
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
        assert_eq!(
            sequence.read_witness,
            ReadWitness::Complete,
            "only a complete observation can be projected onto the locus {}; the partial \
             at {} saw less than its own region",
            self.reference.region,
            self.member.region,
        );

        let bases = &self.reference.bases;
        projected.clear();
        projected.reserve(
            bases
                .len()
                .saturating_sub(self.covered)
                .saturating_add(sequence.bases.len()),
        );
        projected.extend_from_slice(&bases[..self.offset]);
        projected.extend_from_slice(&sequence.bases);
        projected.extend_from_slice(&bases[self.offset + self.covered..]);
    }
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
    use crate::ng::run::cohort_merge::close::SampleMembers;
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
}
