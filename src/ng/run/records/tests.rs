//! What a record needs that the called locus does not carry — the padding base, the rule that
//! decides which loci reach the file, and the four numberings the per-allele counts cross.

use super::*;
use crate::ng::calling::ExpectedAlleleCopies;
use crate::ng::calling::allele_candidates::generic::select_generic;
use crate::ng::calling::allele_candidates::{CandidateSelectionConfig, SelectionScratch};
use crate::ng::calling::quality::ArtifactTestCounts;
use crate::ng::locus_generation::{LocusKind, WitnessedLocusPositions};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::ref_seq::InMemoryRefSeq;
use crate::ng::run::cohort_merge::build::{
    AlleleSupport, PartialObservation, SampleSupport, SupportedAllele,
};
use crate::ng::types::{AlleleId, ContigId, Genotype, Phred, ReadGroupId};

/// A one-contig reference spelling `bases`, so a fetch's answer is legible from the fixture.
fn reference_of(bases: &[u8]) -> InMemoryRefSeq {
    InMemoryRefSeq::from_contigs(vec![bases.to_vec()])
}

/// A region on contig 0, 1-based and inclusive at both ends.
fn region(start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(0),
        start: Position(start),
        end: Position(end),
    }
}

/// An allele table over `alleles`, the reference first.
fn table_of(alleles: &[&[u8]]) -> CandidateAlleles {
    let mut table = CandidateAlleles::new(Box::from(alleles[0]), LocusKind::Generic);
    for allele in &alleles[1..] {
        table.admit(Box::from(*allele));
    }
    table
}

// ---------------------------------------------------------------------
// The padding base
// ---------------------------------------------------------------------

/// **A record every allele of which spells bases needs no padding base**, and asking for one
/// would move its `POS` a base to the left of where the locus is.
#[test]
fn a_record_whose_alleles_all_spell_bases_takes_no_padding_base() {
    let reference = reference_of(b"ACGTACGT");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(3, 3),
        &table_of(&[b"G", b"T"]),
        &mut scratch,
    )
    .expect("a fetch that is never made cannot fail");

    assert_eq!(answer, None);
}

/// **An empty allele takes the reference base immediately to the left of the span.**
///
/// The reference is `ACGTACGT` and the span is bases 3 and 4 — `GT` — so the base to its left
/// is the `C` at position 2. A rule that read the span's own first base would answer `G`, and a
/// rule off by one the other way would answer `A`; both are legal bases, so the fixture spells
/// four different letters to tell them apart.
#[test]
fn an_empty_allele_takes_the_reference_base_to_the_left_of_the_span() {
    let reference = reference_of(b"ACGTACGT");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(3, 4),
        &table_of(&[b"GT", b""]),
        &mut scratch,
    )
    .expect("the reference holds a base at position 2");

    assert_eq!(answer, Some(PaddingBase::Left(b'C')));
}

/// **At a span starting at the contig's first base there is nothing to the left**, so the base
/// immediately to the right is appended instead and `POS` does not move (spec §5).
///
/// The span is bases 1 and 2 — `AC` — so the base to its right is the `G` at position 3.
#[test]
fn a_span_at_the_contigs_first_base_takes_the_base_to_its_right() {
    let reference = reference_of(b"ACGTACGT");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(1, 2),
        &table_of(&[b"AC", b""]),
        &mut scratch,
    )
    .expect("the reference holds a base at position 3");

    assert_eq!(answer, Some(PaddingBase::Right(b'G')));
}

/// **The shape the generic path actually produces needs no padding base**, which is why the
/// rule above is answered and never exercised on today's run.
///
/// The generic mint anchors its indels — an insertion's reference span is its anchor base alone
/// and a deletion's is the anchor plus the deleted run
/// (`ReadEvent::footprint_span`, `locus_generation/pileup/decompose.rs`) — so a five-base
/// deletion is `REF ACGTTG` against `ALT A`, and the alternative spells a base rather than
/// nothing. The empty allele spec §5 was written for is the repeat-tract path's full-tract
/// deletion, which is unbuilt.
#[test]
fn a_deletion_shaped_the_way_the_generic_mint_shapes_one_needs_no_padding_base() {
    let reference = reference_of(b"ACGTTGCA");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(1, 6),
        &table_of(&[b"ACGTTG", b"A"]),
        &mut scratch,
    )
    .expect("a fetch that is never made cannot fail");

    assert_eq!(
        answer, None,
        "the deletion's alternative is its anchor base, so no allele of the record is empty",
    );
}

/// **A span starting at the contig's *second* base still takes the base to its left.**
///
/// The contig-start branch is the one case with no base to the left, and it is exactly one
/// position wide. A rule relaxed by one — `<= 2` where it means `== 1` — appends the base to
/// the right instead and leaves `POS` where it was, which is a wrong `REF`, a wrong `ALT` and a
/// wrong position on every indel near a contig's start. Nothing but a fixture at position 2 can
/// tell the two apart.
#[test]
fn a_span_starting_at_the_second_base_of_a_contig_still_pads_from_the_left() {
    let reference = reference_of(b"ACGTACGT");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(2, 3),
        &table_of(&[b"CG", b""]),
        &mut scratch,
    )
    .expect("the reference holds a base at position 1");

    assert_eq!(
        answer,
        Some(PaddingBase::Left(b'A')),
        "the `A` at position 1, prefixed, with POS moving onto it — not the `T` at position 4",
    );
}

/// **A base the reference cannot serve is an error, never an invented `N`.**
///
/// Production's repeat-tract writer puts the letter `N` at a span it cannot read beside
/// ([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs)); spec §5 declines to port
/// that, because a base the reference does not contain, written at an unshifted position, is a
/// record that parses and lies. The state is unreachable from real reads — it needs a span
/// covering a whole contig — and the refusal is what makes that a claim rather than a hope.
#[test]
fn a_padding_base_the_reference_cannot_serve_is_refused_rather_than_invented() {
    let reference = reference_of(b"AC");
    let mut scratch = Vec::new();

    let answer = padding_base_beside(
        &reference,
        region(1, 2),
        &table_of(&[b"AC", b""]),
        &mut scratch,
    );

    assert!(
        matches!(answer, Err(RefSeqError::OutOfBounds { .. })),
        "a span covering the whole contig has no base on either side of it, and got {answer:?}",
    );
}

/// **The scratch buffer is the fetch's destination and not a place an answer can survive
/// from.** A run reuses one buffer over a genome, so a fetch that appended, or one whose answer
/// was read before the write, would hand the previous locus's base to this one.
#[test]
fn a_reused_scratch_buffer_carries_no_base_from_the_locus_before() {
    // `ACGTTGCA`: position 2 is `C` and position 6 is `G`, so the two answers differ and a
    // stale buffer would show.
    let reference = reference_of(b"ACGTTGCA");
    let mut scratch = Vec::new();

    let first = padding_base_beside(
        &reference,
        region(3, 4),
        &table_of(&[b"GT", b""]),
        &mut scratch,
    )
    .expect("a base at position 2");
    let second = padding_base_beside(
        &reference,
        region(7, 8),
        &table_of(&[b"CA", b""]),
        &mut scratch,
    )
    .expect("a base at position 6");

    assert_eq!(first, Some(PaddingBase::Left(b'C')));
    assert_eq!(
        second,
        Some(PaddingBase::Left(b'G')),
        "the base at position 6, not the one the first call left behind",
    );
}

// ---------------------------------------------------------------------
// Which loci reach the file
// ---------------------------------------------------------------------

/// A call of `genotype`, at a fixed quality, whose reads either said something or did not.
fn called(genotype: &[u16], reads_were_uninformative: bool) -> SampleGenotypeCall {
    SampleGenotypeCall::Called {
        genotype: Genotype::new(genotype.iter().map(|allele| AlleleId(*allele)).collect()),
        genotype_quality: Phred::try_new(30.0).expect("a quality"),
        reads_were_uninformative,
    }
}

/// A locus over one reference and one alternative, called `per_sample`.
fn locus_called(per_sample: Vec<SampleGenotypeCall>) -> LocusInference {
    locus_called_over(&[b"A", b"C"], per_sample, true, None)
}

/// A locus over `alleles`, called `per_sample`, with the loop's own verdict and its artifact
/// counts named.
fn locus_called_over(
    alleles: &[&[u8]],
    per_sample: Vec<SampleGenotypeCall>,
    converged: bool,
    artifact_test_counts: Option<ArtifactTestCounts>,
) -> LocusInference {
    let table = table_of(alleles);
    let copies = ExpectedAlleleCopies::new(vec![1.0; table.len()], &table);
    LocusInference::new(
        region(100, 100),
        table,
        per_sample,
        copies,
        converged,
        1,
        Provenance::Defaulted,
        None,
        Phred::try_new(60.0).expect("a quality"),
        artifact_test_counts,
    )
}

/// **A locus every sample was called homozygous-reference at establishes no variant**, so it is
/// left out of the file: its absence is what says *nothing here* (spec §9).
#[test]
fn a_locus_every_sample_is_reference_at_is_not_written() {
    let locus = locus_called(vec![called(&[0, 0], false), called(&[0, 0], false)]);

    assert!(!a_written_genotype_carries_an_alternative(&locus));
}

/// One sample carrying an alternative is what keeps the locus in the file.
#[test]
fn a_locus_one_sample_carries_an_alternative_at_is_written() {
    let locus = locus_called(vec![called(&[0, 0], false), called(&[0, 1], false)]);

    assert!(a_written_genotype_carries_an_alternative(&locus));
}

/// **The test is on the calls the file would write, not on the ones the loop made.**
///
/// A sample whose reads said nothing is written `./.` (spec §7.1), so its genotype — which the
/// prior alone produced — carries no allele into the file and cannot be what keeps a locus in
/// it. A rule reading `genotype` without asking `reads_were_uninformative` would write this
/// record with an `ALT` no written sample carries.
#[test]
fn a_heterozygote_no_sample_showed_reads_for_does_not_keep_a_locus_in_the_file() {
    let locus = locus_called(vec![called(&[0, 0], false), called(&[0, 1], true)]);

    assert!(!a_written_genotype_carries_an_alternative(&locus));
}

/// A locus where every sample was set aside by the candidate step has no genotype to write,
/// so nothing establishes a variant there either.
#[test]
fn a_locus_every_sample_was_set_aside_at_is_not_written() {
    let locus = locus_called(vec![
        SampleGenotypeCall::Missing,
        SampleGenotypeCall::Missing,
    ]);

    assert!(!a_written_genotype_carries_an_alternative(&locus));
}

// ---------------------------------------------------------------------
// What a record's per-sample and per-allele numbers are read off
// ---------------------------------------------------------------------

/// One allele's row from the sample's first read group — the shape of a sample carrying one
/// library, which is most samples of most runs.
fn row(allele: usize, num_reads: u32, mapq_sum: u32) -> SupportedAllele {
    row_from_group(allele, ReadGroupId(0), num_reads, mapq_sum)
}

/// One `(allele, read group)` row, naming the group — for the fixtures about a sample
/// sequenced from two libraries.
fn row_from_group(
    allele: usize,
    read_group: ReadGroupId,
    num_reads: u32,
    mapq_sum: u32,
) -> SupportedAllele {
    SupportedAllele {
        allele,
        read_group,
        support: AlleleSupport {
            num_reads,
            q_sum: -f64::from(num_reads),
            mapq_sum,
            ..AlleleSupport::default()
        },
    }
}

/// One covering sample, naming which sample of the run it is.
fn covering(sample: usize, rows: Vec<SupportedAllele>) -> SampleSupport {
    SampleSupport {
        sample,
        supported: rows,
        partials: Vec::new(),
        reads_without_observation: 0,
        reads_removed_as_evidence: 0,
        reads_composed_across_records: 0,
    }
}

/// A cohort locus over `alleles`, covered by `per_sample`.
fn observed(alleles: &[&[u8]], per_sample: Vec<SampleSupport>) -> CohortObservation {
    CohortObservation {
        region: region(100, 100),
        alleles: alleles.iter().map(|bases| Box::from(*bases)).collect(),
        per_sample,
    }
}

/// Narrow an observation the way a run does, and hand back what a record is built from.
fn selected(observation: &CohortObservation) -> (AlleleRemap, Vec<UnmatchedSupport>) {
    let selection = select_generic(
        observation,
        &CandidateSelectionConfig::DEFAULT,
        &mut SelectionScratch::new(),
    );
    let (_alleles, _verdict, unmatched, remap) = selection.into_parts();
    (remap, unmatched)
}

/// **A read on an allele candidate selection dropped belongs in `DP` and in no `AD` slot**, and
/// the surviving alleles' counts land on the *record's* dense ids rather than the merge's.
///
/// The merge's table here holds five sequences and the sample shows all five; the support rule
/// drops the one it showed a single read of, so the merge's indices 1, 3 and 4 become the
/// record's 1, 2 and 3, with a hole where index 2 was. **A table without a hole cannot tell a
/// correct remapping from one that merely counts up**, which is why this fixture has one: an
/// off-by-one would put allele 3's four reads under allele 4's `ALT`, and the record would look
/// ordinary.
#[test]
fn ad_is_keyed_by_the_records_allele_table_and_not_the_merges() {
    let observation = observed(
        &[b"A", b"C", b"G", b"T", b"AA"],
        vec![covering(
            0,
            vec![
                row(0, 4, 240),
                row(1, 3, 180),
                row(2, 1, 60),
                row(3, 2, 120),
                row(4, 2, 120),
            ],
        )],
    );
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(
        &[b"A", b"C", b"T", b"AA"],
        vec![called(&[0, 1], false)],
        true,
        None,
    );

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(
        evidence.samples[0].allele_reads,
        vec![4, 3, 2, 2],
        "the reference's four reads, then the three surviving alternatives in the record's own \
         order — the merge's allele 3 is the record's allele 2",
    );
    assert_eq!(
        evidence.samples[0].reads_no_written_allele_explains, 1,
        "the one read on the sequence the support rule dropped",
    );
}

/// **The cohort's mapping qualities are pooled over the same reads `AD` counts**, which is the
/// invariant `VcfRecord::new` refuses a record for breaking: `MQREF` and `MQALT` are means over
/// the reads `AD` names, so two totals that differ are two different pools.
#[test]
fn the_pooled_mapping_qualities_count_the_same_reads_as_ad() {
    let observation = observed(
        &[b"A", b"C"],
        vec![
            covering(0, vec![row(0, 4, 240), row(1, 3, 150)]),
            covering(1, vec![row(0, 2, 100), row(1, 5, 200)]),
        ],
    );
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(
        &[b"A", b"C"],
        vec![called(&[0, 1], false), called(&[1, 1], false)],
        true,
        None,
    );

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    for allele in 0..2 {
        let attributed: u64 = evidence
            .samples
            .iter()
            .map(|sample| u64::from(sample.allele_reads[allele]))
            .sum();
        assert_eq!(
            evidence.allele_mapq[allele].reads, attributed,
            "allele {allele}'s pool and the samples' AD count different reads",
        );
    }
    assert_eq!(evidence.allele_mapq[0].mapq_sum, 340);
    assert_eq!(evidence.allele_mapq[1].mapq_sum, 350);
}

/// **Each covering sample's evidence lands on the run sample it names, not on where it sits in
/// the merge's list.**
///
/// The merge holds only the samples that covered the locus, each naming its index in the run's
/// order; the record holds one entry per sample of the run. Here the run has three samples and
/// its sample 1 covered nothing, so the merge's two entries name the run's samples 0 and 2 —
/// and a record built positionally would write sample 2's seven alternative reads under sample
/// 1's name and leave sample 2 with none.
#[test]
fn a_covering_samples_reads_land_on_the_run_sample_it_names() {
    let observation = observed(
        &[b"A", b"C"],
        vec![
            covering(0, vec![row(0, 9, 540), row(1, 2, 120)]),
            covering(2, vec![row(0, 1, 60), row(1, 7, 420)]),
        ],
    );
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(
        &[b"A", b"C"],
        vec![
            called(&[0, 0], false),
            called(&[0, 0], true),
            called(&[1, 1], false),
        ],
        true,
        None,
    );

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(evidence.samples[0].allele_reads, vec![9, 2]);
    assert_eq!(
        evidence.samples[1].allele_reads,
        vec![0, 0],
        "the run's sample 1 covered nothing, and a record says so with zeros rather than with \
         somebody else's reads",
    );
    assert_eq!(evidence.samples[2].allele_reads, vec![1, 7]);
}

/// **A read that stopped inside the locus is in `DP` and in no `AD` slot.**
///
/// A partial observation says the sample carries *at least* this much and never reached the
/// allele table, so no written allele explains it — which is exactly what `DP − ΣAD` is for
/// (spec §7). It is a different fact from a dropped candidate's reads and the two are added.
#[test]
fn a_read_that_stopped_inside_the_locus_is_in_dp_and_in_no_ad_slot() {
    let mut sample = covering(0, vec![row(0, 4, 240), row(1, 3, 180)]);
    sample.partials = vec![PartialObservation {
        witnessed_in_locus: WitnessedLocusPositions::one_run_from_offset_and_length(0, 1)
            .expect("a one-position witness"),
        read_group: ReadGroupId(0),
        bases: Box::from(b"A".as_slice()),
        num_reads: 5,
        q_sum: -5.0,
    }];
    let observation = observed(&[b"A", b"C"], vec![sample]);
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(evidence.samples[0].allele_reads, vec![4, 3]);
    assert_eq!(
        evidence.samples[0].reads_no_written_allele_explains, 5,
        "the five reads that stopped inside the locus, which no allele of the record spells",
    );
}

/// **Pooled counts the two artifact tests both charge on**: 20 reads at the reference and 4 at
/// the alternative where a heterozygote expects 12, and every one of the alternative's reads on
/// one strand and one side of its record.
///
/// A fixture the tests charged nothing on could not tell the corrected quality from the
/// baseline, which is the whole of what the assertion above is for.
fn counts_the_tests_charge_on() -> ArtifactTestCounts {
    ArtifactTestCounts {
        primary_alternative: AlleleId(1),
        reference_reads: 20.0,
        reference_forward_reads: 10.0,
        reference_placed_left_reads: 10.0,
        alternative_reads: 4.0,
        alternative_forward_reads: 4.0,
        alternative_placed_left_reads: 4.0,
        total_reads: 24.0,
        genotype_expected_alternative_reads: 12.0,
    }
}

/// **The record carries the site quality after the artifact correction, and the two penalties
/// it charged** — never the loop's own field, which is the uncorrected baseline
/// (`doc/devel/ng/spec/calling_quality.md` §3.5).
#[test]
fn the_records_site_quality_is_the_corrected_one_and_the_penalties_are_beside_it() {
    let observation = observed(
        &[b"A", b"C"],
        vec![covering(0, vec![row(0, 4, 240), row(1, 3, 180)])],
    );
    let (remap, unmatched) = selected(&observation);
    let counts = counts_the_tests_charge_on();
    let locus = locus_called_over(
        &[b"A", b"C"],
        vec![called(&[0, 1], false)],
        true,
        Some(counts),
    );

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    let (expected, penalties) = correct_site_quality(locus.uncorrected_site_quality(), &counts);
    assert_eq!(evidence.corrected_site_quality, expected);
    assert_eq!(evidence.artifact_penalties, Some(penalties));
    assert!(
        evidence.corrected_site_quality < locus.uncorrected_site_quality(),
        "a fixture whose tests charge nothing could not tell the corrected number from the \
         baseline: {:?} against {:?}",
        evidence.corrected_site_quality,
        locus.uncorrected_site_quality(),
    );
}

/// **A locus that gave the two artifact tests nothing to weigh keeps its baseline**, and says
/// the tests did not run rather than writing two zeroed penalties — which a reader could not
/// tell from two tests that ran and charged nothing.
#[test]
fn a_locus_the_artifact_tests_did_not_run_on_carries_no_penalties() {
    let observation = observed(
        &[b"A", b"C"],
        vec![covering(0, vec![row(0, 4, 240), row(1, 3, 180)])],
    );
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(evidence.artifact_penalties, None);
    assert_eq!(
        evidence.corrected_site_quality,
        locus.uncorrected_site_quality()
    );
}

/// **A locus whose loop did not settle is written on the `EMNoConv` filter, not dropped** — one
/// hard locus must not cost a cohort its record, and a genotype from a loop that did not settle
/// is a weaker claim than one from a loop that did (spec §8).
#[test]
fn a_locus_whose_loop_did_not_settle_carries_the_filter_that_says_so() {
    let observation = observed(
        &[b"A", b"C"],
        vec![covering(0, vec![row(0, 4, 240), row(1, 3, 180)])],
    );
    let (remap, unmatched) = selected(&observation);

    let settled = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);
    let did_not = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], false, None);

    assert_eq!(
        evidence_for_output(&settled, &observation, &remap, &unmatched, None).filter,
        FilterVerdict::Pass,
    );
    assert_eq!(
        evidence_for_output(&did_not, &observation, &remap, &unmatched, None).filter,
        FilterVerdict::EmDidNotConverge,
    );
}

/// **One sample's reads for one allele are summed over its read groups.**
///
/// A sample sequenced from two libraries has two rows for one allele, and `AD` is a count of
/// reads however many lanes produced them — the merge keeps the rows apart because a read
/// likelihood may not pool two error rates, which is a rule about the model and not about the
/// file (`doc/devel/ng/spec/read_groups.md` §2.3). Assigning rather than summing keeps only the
/// last lane, and that does not read as a wrong number: it makes `AD` disagree with the pooled
/// mapping-quality count beside it, which `VcfRecord::new` refuses — **a run that panics at the
/// first two-library sample it meets**. Measured in a surveyed tomato archive, 157 samples in
/// 1,707 carry more than one read group, so that is about one sample in eleven.
#[test]
fn a_samples_reads_for_one_allele_are_summed_over_its_read_groups() {
    let observation = observed(
        &[b"A", b"C"],
        vec![covering(
            0,
            vec![
                row(0, 4, 240),
                row_from_group(1, ReadGroupId(0), 3, 180),
                row_from_group(1, ReadGroupId(1), 2, 100),
            ],
        )],
    );
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(
        evidence.samples[0].allele_reads,
        vec![4, 5],
        "three reads from one lane and two from the other are five reads on the alternative",
    );
    assert_eq!(evidence.allele_mapq[1].reads, 5);
    assert_eq!(evidence.allele_mapq[1].mapq_sum, 280);
}

/// **What candidate selection set aside for a sample is read off the merge's covering list, not
/// the run's sample order** — and the two differ whenever some sample of the run covered
/// nothing.
///
/// Here the run has three samples and its sample 1 covered nothing, so the merge's two entries
/// name the run's samples 0 and 2 while selection's leftover is parallel to the merge's list.
/// The covering sample that lost a sequence is the merge's entry 1, which is the run's sample 2
/// — so a leftover read by the run index would credit the reads to nobody and leave sample 2's
/// depth one short, which is a wrong number in the file and no crash anywhere.
#[test]
fn a_samples_leftover_is_read_off_the_merges_covering_list_and_not_the_runs_order() {
    let observation = observed(
        &[b"A", b"C", b"G"],
        vec![
            covering(0, vec![row(0, 5, 300), row(1, 4, 240)]),
            covering(2, vec![row(0, 5, 300), row(1, 4, 240), row(2, 1, 60)]),
        ],
    );
    let (remap, unmatched) = selected(&observation);
    assert_eq!(
        remap.candidate_for(2),
        None,
        "the fixture depends on the `G` being dropped: one read against a floor of two",
    );
    let locus = locus_called_over(
        &[b"A", b"C"],
        vec![
            called(&[0, 1], false),
            called(&[0, 0], true),
            called(&[0, 1], false),
        ],
        true,
        None,
    );

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(
        evidence.samples[2].reads_no_written_allele_explains, 1,
        "the run's sample 2 is the merge's entry 1, and it is the one that lost a sequence",
    );
    assert_eq!(evidence.samples[0].reads_no_written_allele_explains, 0);
    assert_eq!(
        evidence.samples[1].reads_no_written_allele_explains, 0,
        "the run's sample 1 covered nothing and lost nothing",
    );
}

/// **A read removed as evidence is in `DP` and in no `AD` slot** — spec §7's `DP` is every read
/// observation the sample had at the locus, and a read named at some of that sample's records
/// inside the locus and not at all of them was observed there. The merge calls it lost depth and
/// counts it for exactly this reason; leaving it out makes `DP` understate the depth at the loci
/// that span several of a sample's records.
#[test]
fn a_read_removed_as_evidence_is_in_dp_and_in_no_ad_slot() {
    let mut sample = covering(0, vec![row(0, 4, 240), row(1, 3, 180)]);
    sample.reads_removed_as_evidence = 2;
    let observation = observed(&[b"A", b"C"], vec![sample]);
    let (remap, unmatched) = selected(&observation);
    let locus = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(evidence.samples[0].allele_reads, vec![4, 3]);
    assert_eq!(
        evidence.samples[0].reads_no_written_allele_explains, 2,
        "the two reads the merge could not attribute to any of this sample's alleles",
    );
}

/// **The three kinds of unexplained read are added, not chosen between.**
///
/// A dropped candidate's reads, a read that stopped inside the locus, and a read removed as
/// evidence are three different facts about one sample, and `DP − ΣAD` is all of them. A
/// version that kept any one would understate the depth on a sample carrying the other two.
#[test]
fn the_three_kinds_of_unexplained_read_are_added_together() {
    let mut sample = covering(0, vec![row(0, 4, 240), row(1, 3, 180), row(2, 1, 60)]);
    sample.reads_removed_as_evidence = 2;
    sample.partials = vec![PartialObservation {
        witnessed_in_locus: WitnessedLocusPositions::one_run_from_offset_and_length(0, 1)
            .expect("a one-position witness"),
        read_group: ReadGroupId(0),
        bases: Box::from(b"A".as_slice()),
        num_reads: 5,
        q_sum: -5.0,
    }];
    let observation = observed(&[b"A", b"C", b"G"], vec![sample]);
    let (remap, unmatched) = selected(&observation);
    assert_eq!(
        remap.candidate_for(2),
        None,
        "the `G` is dropped at one read"
    );
    let locus = locus_called_over(&[b"A", b"C"], vec![called(&[0, 1], false)], true, None);

    let evidence = evidence_for_output(&locus, &observation, &remap, &unmatched, None);

    assert_eq!(evidence.samples[0].allele_reads, vec![4, 3]);
    assert_eq!(
        evidence.samples[0].reads_no_written_allele_explains,
        1 + 5 + 2,
        "one read on the dropped sequence, five that stopped inside the locus, two removed as \
         evidence — three counts, and the file states their sum",
    );
}
