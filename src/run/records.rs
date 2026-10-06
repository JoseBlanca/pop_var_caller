//! **What a record needs that the called locus does not carry** — the run's half of the seam
//! [`vcf::assemble`](crate::vcf::assemble) describes from the other side.
//!
//! `assemble_record` turns a [`LocusInference`] plus a [`LocusEvidenceForOutput`] into a
//! record, and its module note says the second of those is *"this module's guess at the shape
//! the stream will hand over"*. This is the stream. Everything here is read while the merge's
//! [`CohortObservation`] is still in hand, because none of it can be recovered afterwards: the
//! observation is dropped as soon as its genotypes exist, which is what calling inside the
//! builder buys (`doc/devel/ng/arch/run_streaming.md` §3.4).
//!
//! # The three numberings meet here again, and one of them is new
//!
//! [`shape_generic_locus`](crate::calling::evidence_shaping) already has to reconcile the
//! merge's covering samples, the run's sample order and the calling scratch's rows. The record
//! adds a fourth axis that is not a sample axis at all: **the merge's allele table and the
//! record's are different tables**, because candidate selection drops alleles between them. So
//! every per-allele count here is keyed by [`AlleleRemap`] and never by the merge's own index —
//! a read on an allele selection dropped belongs in `DP` and in no `AD` slot
//! (`doc/devel/ng/spec/vcf_output.md` §7), and using the merge's index would put it in the
//! wrong one instead.
//!
//! # The padding base, and why it is fetched here
//!
//! VCF cannot spell an empty allele, so a record with one is written by giving every allele a
//! flanking reference base (spec §5). **The locus does not carry that base**: the merge gathers
//! the bases of the span and no more. So it is fetched from the run's own reference — one
//! accessor for the whole run, minted beside the walkers' — at the moment the record is built,
//! which is the only moment that knows both the span and whether any allele is empty.
//!
//! **⚑ And on today's path it is never fetched, which is worth knowing before reading this as
//! live code.** Every generic record starts at an anchor base that both of its alleles keep: a
//! deletion's reference span is the anchor plus the deleted run and an insertion's is the anchor
//! plus its own inserted length (`ReadEvent::record_span`,
//! `locus_generation/pileup/decompose.rs`), so a deletion's alternative is the anchor base and
//! never nothing, an insertion's is longer still, and **no allele a generic locus is called over
//! is ever empty**. The empty allele spec §5 was written for is the repeat-tract
//! path's full-tract deletion, and that path is unbuilt. What this is, then, is the answer
//! ready for when a record needs it — and it cannot be left out: `VcfRecord::new` asserts a
//! padding base is carried **exactly** when some allele is empty, so a run that did not compute
//! one would panic at the first record that needed it rather than write a wrong one.
//!
//! **Never invented.** Production's repeat-tract writer puts the letter `N` at a span starting
//! at a contig's first base ([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs));
//! spec §5 declines to port that, so a base that cannot be read is a run that stops and says
//! which position it could not read.

use crate::calling::allele_candidates::{AlleleRemap, SelectionVerdict, UnmatchedSupport};
use crate::calling::likelihood::generic::{
    carrier_sequence_ends_with, carrier_sequence_starts_with,
};
use crate::calling::quality::artifact_correction::correct_site_quality;
use crate::calling::{CandidateAlleles, LocusInference, ReferenceBesideLocus, SampleGenotypeCall};
use crate::locus_generation::{LocusKind, LocusLen, WitnessedLocusPositions};
use crate::ref_seq::{EvictableRefSeq, RefSeq, RefSeqError};
use crate::run::cohort_merge::build::{CohortObservation, LocusBorder, PartialObservation};
use crate::types::{GenomeRegion, Position};
use crate::vcf::assemble::{LocusEvidenceForOutput, SampleEvidenceForOutput};
use crate::vcf::{FilterVerdict, MapqPool, PaddingBase, TractAnnotation};
use crate::window_coverage::WindowCoverage;

/// The first base of every contig, 1-based — the one position with no base to its left.
const FIRST_POSITION_OF_A_CONTIG: Position = Position(1);

/// **Whether this locus establishes a variant at all**, which is what decides if it reaches
/// the file (`doc/devel/ng/spec/vcf_output.md` §9).
///
/// A locus every sample was called homozygous-reference at, or no sample was called at,
/// established nothing: its absence from the file says *nothing here*, which is exactly what
/// production's generic path means by leaving it out. There is no gVCF and no reference block.
///
/// **The test is on the calls the file would write, not on the ones the loop made.** A sample
/// whose reads said nothing is written `./.` (spec §7.1) and carries no allele, so it cannot be
/// what keeps a locus in the file — which is the same rule
/// [`assemble_record`](crate::vcf::assemble::assemble_record) applies one step later.
#[must_use]
pub fn a_written_genotype_carries_an_alternative(locus: &LocusInference) -> bool {
    locus.per_sample.iter().any(|call| match call {
        SampleGenotypeCall::Missing => false,
        SampleGenotypeCall::Called {
            genotype,
            reads_were_uninformative,
            ..
        } => {
            !*reads_were_uninformative
                && genotype
                    .alleles()
                    .iter()
                    .any(|allele| !allele.is_reference())
        }
    })
}

/// The two buffers [`reference_beside_locus`] reads the flanks into, kept by a worker across
/// loci so a run allocates them once.
#[derive(Debug, Default)]
pub struct ReferenceBesideScratch {
    before: Vec<u8>,
    after: Vec<u8>,
}

/// **The reference either side of an ordinary locus, as much as its partial reads need** — what
/// a read that ran out inside the locus is compared against after the allele
/// (`doc/devel/ng/spec/read_likelihoods.md` §5.3).
///
/// A read flush to the left border shows the allele and then the reference after the locus; one
/// flush to the right shows the reference before it and then the allele. The flank is needed only
/// past the allele's end, so **each side is read as far as the longest such partial showed beyond
/// the locus's shortest allele, and not at all where none did** — the ordinary case, since most
/// loci have no partial longer than any allele. The alleles the loop is called over are always
/// among the merge's, so the shortest of the merge's bounds them all.
///
/// **A repeat tract gets nothing**: its reads are re-aligned against flanks the tract carries,
/// and this rule is the ordinary path's. A contig's end clamps the side that would run past it.
///
/// # Errors
///
/// Whatever the reference fetch refuses, other than running past the contig's end.
pub fn reference_beside_locus<'s, R>(
    reference: &R,
    observation: &CohortObservation,
    scratch: &'s mut ReferenceBesideScratch,
) -> Result<ReferenceBesideLocus<'s>, RefSeqError>
where
    R: RefSeq,
{
    scratch.before.clear();
    scratch.after.clear();
    if !matches!(observation.kind, LocusKind::Generic) {
        return Ok(ReferenceBesideLocus::NONE);
    }
    let Some(shortest_allele) = observation.alleles.iter().map(|allele| allele.len()).min() else {
        return Ok(ReferenceBesideLocus::NONE);
    };
    let locus_len = LocusLen::from_positions(observation.alleles[0].len() as u64);
    let mut needed_before = 0_usize;
    let mut needed_after = 0_usize;
    for partial in observation
        .per_sample
        .iter()
        .flat_map(|sample| &sample.partials)
    {
        let witness = &partial.witnessed_in_locus;
        let beyond_the_shortest_allele = partial.bases.len().saturating_sub(shortest_allele);
        // A row re-read after the merge is compared as a read flush to the border it is marked
        // with, whatever its witness ([`reread_spellings_cut_short`]).
        match partial.sequence_may_run_on_past {
            Some(LocusBorder::Right) => {
                needed_after = needed_after.max(beyond_the_shortest_allele);
                continue;
            }
            Some(LocusBorder::Left) => {
                needed_before = needed_before.max(beyond_the_shortest_allele);
                continue;
            }
            None => {}
        }
        // Only a single run anchored at exactly one border is compared against a flank; every
        // other shape is compared against the allele alone or against nothing.
        if witness.runs().len() != 1 {
            continue;
        }
        match (witness.is_flush_left(), witness.is_flush_right(locus_len)) {
            (true, false) => needed_after = needed_after.max(beyond_the_shortest_allele),
            (false, true) => needed_before = needed_before.max(beyond_the_shortest_allele),
            _ => {}
        }
    }

    let region = observation.region;
    if needed_before > 0 {
        let first = region.start.get();
        let start = first
            .saturating_sub(needed_before as u64)
            .max(FIRST_POSITION_OF_A_CONTIG.get());
        if start < first {
            reference.fetch_into(region.contig, start, first - start, &mut scratch.before)?;
        }
    }
    if needed_after > 0 {
        let start = region.end.get() + 1;
        match reference.fetch_into(
            region.contig,
            start,
            needed_after as u64,
            &mut scratch.after,
        ) {
            Ok(()) => {}
            Err(RefSeqError::OutOfBounds { contig_length, .. }) => {
                let available = (contig_length + 1).saturating_sub(start);
                if available > 0 {
                    reference.fetch_into(region.contig, start, available, &mut scratch.after)?;
                }
            }
            Err(other) => return Err(other),
        }
    }
    Ok(ReferenceBesideLocus {
        before: &scratch.before,
        after: &scratch.after,
    })
}

/// **The site quality below which a called locus is not written**, by default — `--min-site-quality`.
///
/// The threshold reads the quality the record would carry: the site quality after the artifact
/// correction (`calling_quality.md` §3.5, which pins that the gate and the `QUAL` column read one
/// number). A site whose two artifact penalties outweigh its evidence is floored at zero, and
/// those are most of the false calls this caller writes.
///
/// **One, not GATK's thirty, because the site quality already grows with depth.** Measured on the
/// three GIAB samples, each over its own 100 regions (false / true calls removed): below 1 removes
/// 141 / 15 at 5×, 215 / 17 at 10×, 142 / 2 at 30× and 32 / 3 at 300×; below 30 removes 212 / 117,
/// 275 / 131, 148 / 18 and 33 / 3, so the extra cut costs more than a hundred true calls at each
/// low depth for a few dozen false ones. On 50 tomato accessions at about 3×, against GATK's joint
/// calls over all 80 regions, below 1 removes 26,696 sample genotypes GATK does not call and 382
/// it does. **Zero turns the threshold off** — no quality is below zero.
pub const DEFAULT_MIN_SITE_QUALITY: f32 = 1.0;

/// **The strand bias at or above which a called locus is not written**, by default —
/// `--max-strand-bias`, in Phred, read off
/// [`strand_bias`](crate::calling::quality::artifact_correction::strand_bias).
///
/// **A cutoff and not a larger penalty, because subtraction cannot keep up with depth.** The site
/// quality grows with every variant read, and an artifact's variant reads grow with depth like a
/// real variant's; the strand penalty grows too, but more slowly. At 300× on GIAB HG004
/// chr3:107848623 all 61 variant reads are forward against 33 of 132 reference reads, the penalty
/// is 367 Phred, the baseline about 440, and 8 is left to write the site.
///
/// **Why 100.** Measured on the three GIAB samples at 5×, 10×, 30× and 300×, against GIAB v4.2.1:
/// at 300×, 100 removes 3 false calls and no true one; the largest value on a true call at any
/// depth is 75, a homozygous site, so 60 (GATK's cutoff for its own strand test) would cost it.
/// Below 300× no call reaches 40. In a cohort the reads pool, so it binds far more: on 63 tomato
/// accessions at about 3× it leaves out 3,219 of 193,893 records, and at the 1,780 of those where
/// GATK has a record of its own, GATK's own strand filters flag 1,217 (`calling_quality.md` §6.5).
/// **Zero turns the cutoff off.**
pub const DEFAULT_MAX_STRAND_BIAS: f32 = 100.0;

/// **Re-read the reads the mapper laid straight across the reference where a longer allele began**
/// — one ordinary locus, before candidate selection sees it (item 2 of the GIAB report on
/// 5d2abae4; `doc/devel/ng/spec/read_likelihoods.md` §5.3).
///
/// A read that ends a few bases into an inserted copy cannot be placed with a gap, so the mapper
/// aligns it straight across the reference, and the start of the copy comes out as substitutions.
/// It covers every position of the locus, so the merge stores it as a **complete** observation of
/// a sequence of its own — at GIAB HG003 chr15:96140584, 34 reads showing exactly the first 46
/// bases of a 24-base insertion, against 107 that carry the insertion whole. That sequence is
/// then a candidate allele no haplotype carries, and the homozygous insertion is called `1/2`.
///
/// **The rule.** An alternative allele `S` at least as long as the reference — a spelling of
/// substitutions alone, or of an insertion shorter than `A`'s — whose bases are exactly the start
/// of another allele
/// `A` followed by the reference past the locus, or exactly the end of `A` preceded by the
/// reference before it, is a spelling of `A` cut short rather than evidence of its own. Its
/// complete rows become partial rows covering the whole locus and marked to run on past that
/// border ([`PartialObservation::sequence_may_run_on_past`]), so they count for every allele a
/// carrier of which would have shown the same bases. Three guards: `A` must have at least as many
/// complete reads across the cohort as `S`, so a handful of odd reads cannot re-read a
/// well-supported allele; `A` and `S` must differ in length by two bases or more, because a
/// one-base shift inside a homopolymer is invisible and a real SNP just before one spells the
/// start of a one-base insertion; and the reference is never re-read, since it is the allele every
/// genotype is measured against and the slide in the walk already covers its ambiguous reads.
///
/// `S`'s entry in the allele table is kept with no complete reads behind it, which candidate
/// selection's admission rule then declines. A locus where nothing qualifies is left untouched,
/// which is every locus whose alleles are all one length and every locus the reference before or
/// after cannot be read for.
///
/// # Errors
///
/// Whatever the reference fetch refuses, other than running past the contig's end.
pub fn reread_spellings_cut_short<R>(
    reference: &R,
    observation: &mut CohortObservation,
    scratch: &mut ReferenceBesideScratch,
) -> Result<(), RefSeqError>
where
    R: RefSeq,
{
    if !matches!(observation.kind, LocusKind::Generic) || observation.alleles.len() < 3 {
        return Ok(());
    }
    let lengths = observation.alleles.iter().map(|allele| allele.len());
    let (shortest, longest) =
        lengths.fold((usize::MAX, 0), |(lo, hi), len| (lo.min(len), hi.max(len)));
    // A spelling and the allele it was cut from differ in length, and the flank needed to tell
    // them apart is at most that difference.
    let flank = longest - shortest;
    if flank == 0 {
        return Ok(());
    }
    let mut complete_reads = vec![0_u64; observation.alleles.len()];
    for sample in &observation.per_sample {
        for row in &sample.supported {
            complete_reads[row.allele] += u64::from(row.support.num_reads);
        }
    }
    fetch_flanks(reference, observation.region, flank, scratch)?;

    let reference_len = observation.alleles[0].len();
    let mut cut_short: Vec<Option<LocusBorder>> = vec![None; observation.alleles.len()];
    for (spelling, bases) in observation.alleles.iter().enumerate().skip(1) {
        // **Never a spelling shorter than the reference.** A deletion's own allele is a prefix of
        // the reference's whenever it removes the locus's last bases, and a read carrying it
        // crossed the locus with it — trying it here re-read 20 true deletions at 300x on GIAB.
        if complete_reads[spelling] == 0 || bases.len() < reference_len {
            continue;
        }
        for (longer, allele) in observation.alleles.iter().enumerate() {
            if longer == spelling
                || allele.len() == bases.len()
                || complete_reads[longer] < complete_reads[spelling]
            {
                continue;
            }
            // **Two bases or more apart in length.** A one-base indel shifts what follows by
            // one, which inside a homopolymer cannot be seen: a SNP just before a run of `A`s
            // spells exactly the start of a one-base insertion there. On 63 tomato accessions
            // that re-read a SNP GATK calls in 52 of them (SL4.0ch12:18047655).
            if allele.len().abs_diff(bases.len()) < 2 {
                continue;
            }
            // **A spelling that carries an insertion is read only as a shorter insertion than the
            // allele it was cut from**: a read starting inside an inserted copy, which the mapper
            // fitted with fewer inserted bases. At GIAB HG002 chr1:243535155, 37 reads showed
            // exactly the last 56 bases of a 68-base allele carrying a 19-base insertion. Against
            // a shorter allele it would be a longer insertion, which no read cut short spells.
            if bases.len() > reference_len && allele.len() < bases.len() {
                continue;
            }
            if carrier_sequence_starts_with(allele, &scratch.after, bases) {
                cut_short[spelling] = Some(LocusBorder::Right);
                break;
            }
            if carrier_sequence_ends_with(&scratch.before, allele, bases) {
                cut_short[spelling] = Some(LocusBorder::Left);
                break;
            }
        }
    }
    if cut_short.iter().all(Option::is_none) {
        return Ok(());
    }

    let positions = u16::try_from(observation.alleles[0].len())
        .expect("an ordinary locus is bounded by --max-cohort-locus-span, far inside u16");
    let whole_locus = WitnessedLocusPositions::one_run_from_offset_and_length(0, positions)
        .expect("a locus covers at least its reference base, and its length fits in u16");
    for sample in &mut observation.per_sample {
        let mut moved = false;
        sample.supported.retain(|row| {
            let Some(border) = cut_short[row.allele] else {
                return true;
            };
            sample.partials.push(PartialObservation {
                witnessed_in_locus: whole_locus.clone(),
                read_group: row.read_group,
                bases: observation.alleles[row.allele].clone(),
                num_reads: row.support.num_reads,
                q_sum: row.support.q_sum,
                sequence_may_run_on_past: Some(border),
            });
            moved = true;
            false
        });
        if moved {
            // The order the merge sorts partials in, with the mark last so the key stays total:
            // a re-read row covers the whole locus, which no row the merge minted does.
            sample.partials.sort_unstable_by(|left, right| {
                (
                    &left.witnessed_in_locus,
                    left.read_group,
                    &left.bases,
                    left.sequence_may_run_on_past,
                )
                    .cmp(&(
                        &right.witnessed_in_locus,
                        right.read_group,
                        &right.bases,
                        right.sequence_may_run_on_past,
                    ))
            });
        }
    }
    Ok(())
}

/// The reference `flank` bases either side of `region`, into `scratch`, each stopped at the
/// contig's ends.
fn fetch_flanks<R>(
    reference: &R,
    region: GenomeRegion,
    flank: usize,
    scratch: &mut ReferenceBesideScratch,
) -> Result<(), RefSeqError>
where
    R: RefSeq,
{
    scratch.before.clear();
    scratch.after.clear();
    let first = region.start.get();
    let start = first
        .saturating_sub(flank as u64)
        .max(FIRST_POSITION_OF_A_CONTIG.get());
    if start < first {
        reference.fetch_into(region.contig, start, first - start, &mut scratch.before)?;
    }
    let start = region.end.get() + 1;
    match reference.fetch_into(region.contig, start, flank as u64, &mut scratch.after) {
        Ok(()) => Ok(()),
        Err(RefSeqError::OutOfBounds { contig_length, .. }) => {
            let available = (contig_length + 1).saturating_sub(start);
            if available > 0 {
                reference.fetch_into(region.contig, start, available, &mut scratch.after)?;
            }
            Ok(())
        }
        Err(other) => Err(other),
    }
}

/// **The reference base a record with an empty allele is padded with**, or `None` where every
/// allele spells bases.
///
/// Ordinarily the base immediately to the **left** of the span, with `POS` moving one base left
/// with it. At a span starting at the contig's first base there is nothing to the left, and the
/// base immediately to the **right** is appended instead with `POS` unmoved — the VCF 4.4 rule
/// for an event at position 1 (spec §5).
///
/// `scratch` is the fetch's destination and is cleared by the fetch; it is a parameter so that
/// a run reuses one buffer over a genome rather than allocating a one-byte `Vec` per indel.
///
/// # Errors
///
/// Whatever the reference fetch refuses. The ordinary reachable cause is the FASTA becoming
/// unreadable part-way through a run; the position being out of bounds needs a span that
/// starts at a contig's first base *and* ends at its last, which no aligned read can produce —
/// a read carrying a deletion must match reference on at least one side of it.
pub fn padding_base_beside<R>(
    reference: &R,
    region: GenomeRegion,
    alleles: &CandidateAlleles,
    scratch: &mut Vec<u8>,
) -> Result<Option<PaddingBase>, RefSeqError>
where
    R: RefSeq + EvictableRefSeq,
{
    if !alleles.iter().any(<[u8]>::is_empty) {
        return Ok(None);
    }
    let at_contig_start = region.start == FIRST_POSITION_OF_A_CONTIG;
    let position = if at_contig_start {
        region.end.get() + 1
    } else {
        region.start.get() - 1
    };
    reference.fetch_into(region.contig, position, 1, scratch)?;
    let base = scratch[0];
    // **Release what the run has walked past**, so the accessor's window does not grow into a
    // resident contig over a genome. Correctness does not depend on it: an evicted position is
    // simply read again, and the loci this is driven over arrive in genome order, so nothing
    // asks for a base behind the one just fetched.
    reference.evict_before(position);
    Ok(Some(if at_contig_start {
        PaddingBase::Right(base)
    } else {
        PaddingBase::Left(base)
    }))
}

/// **Everything a record needs that the called locus does not already carry**, gathered while
/// the cohort observation is still in hand.
///
/// The four kinds of thing it collects, and where each comes from:
///
/// - **what each sample's reads showed** — `AD` per written allele, and how many of the
///   sample's reads no written allele explains, which is `DP − ΣAD` (spec §7). Read off the
///   merge's per-sample support rows through `remap`, and off what candidate selection set
///   aside;
/// - **the cohort's mapping qualities per written allele** — the same reads `AD` counts, which
///   is why they are summed from the same rows;
/// - **the site quality after the artifact correction**, and the two penalties it charged.
///   `LocusInference`'s own quality field is the **uncorrected** baseline and nothing between
///   the worker and this correction may read it as a site quality
///   (`doc/devel/ng/spec/calling_quality.md` §3.5);
/// - **the padding base**, already resolved by [`padding_base_beside`].
///
/// # Panics
///
/// On a covering sample naming a run sample the inference was not called over — the merge and
/// the calling loop are handed the same run sample count, so the two disagreeing is a wiring
/// defect in the driver that built both, which is
/// [`assemble_record`](crate::vcf::assemble::assemble_record)'s reasoning for its own
/// checks and the same choice.
///
/// And on an observation whose windows and covering samples are different lengths, for the same
/// reason: the two are one sequence held twice.
#[must_use]
pub fn evidence_for_output(
    locus: &LocusInference,
    observation: &CohortObservation,
    remap: &AlleleRemap,
    unmatched: &[UnmatchedSupport],
    selection: SelectionVerdict,
    padding_base: Option<PaddingBase>,
) -> LocusEvidenceForOutput {
    let written_alleles = locus.alleles().len();
    let run_sample_count = locus.per_sample.len();
    let mut samples: Vec<SampleEvidenceForOutput> = (0..run_sample_count)
        .map(|_| SampleEvidenceForOutput {
            allele_reads: vec![0; written_alleles],
            reads_no_written_allele_explains: 0,
            // **Absent until a covering sample's own pair is copied in below.** A sample that
            // covered nothing at this locus keeps it, which says the same thing its missing
            // entry in the merge's own list says.
            window_coverage: WindowCoverage::absent(),
        })
        .collect();
    let mut allele_mapq = vec![
        MapqPool {
            reads: 0,
            mapq_sum: 0
        };
        written_alleles
    ];

    // **The merge's two per-sample lists are one sequence held twice**, and this is the only
    // place that indexes the second by the first's index. Both fields are public and every
    // fixture builds them separately, so a locus whose lists came apart would otherwise fail as a
    // bare bounds panic below — blaming the reader for what the builder did, which is exactly
    // what the assertion inside the loop exists to prevent for the other index.
    assert_eq!(
        observation.window_coverage.len(),
        observation.per_sample.len(),
        "the locus at {} carries {} windows and {} covering samples: the two are one sequence \
         held twice, so two lengths mean they were not built together",
        observation.region,
        observation.window_coverage.len(),
        observation.per_sample.len(),
    );

    for (covering, support) in observation.per_sample.iter().enumerate() {
        assert!(
            support.sample < run_sample_count,
            "the merge's entry {covering} at {} names run sample {} and the locus was called \
             over {run_sample_count}: both are the run's sample order, so an index past its end \
             means the two were built over different cohorts",
            observation.region,
            support.sample
        );
        let sample = &mut samples[support.sample];
        // **From sparse to dense, at the same index the support came from.** The merge's windows
        // are one per covering sample and this row set is one per run sample; `covering` is the
        // index into the first and `support.sample` into the second, which is exactly the
        // translation this loop already does for the counts.
        sample.window_coverage = observation.window_coverage[covering];
        for row in &support.supported {
            // **The merge's allele index is not the record's.** A row whose allele candidate
            // selection dropped has no slot to be counted in, and its reads reach `DP` through
            // the leftover below rather than through any `AD`.
            let Some(written) = remap.candidate_for(row.allele) else {
                continue;
            };
            let written = usize::from(written.get());
            sample.allele_reads[written] =
                sample.allele_reads[written].saturating_add(row.support.num_reads);
            allele_mapq[written].reads += u64::from(row.support.num_reads);
            allele_mapq[written].mapq_sum += u64::from(row.support.mapq_sum);
        }
        // **Three kinds of read no written allele explains, and they are disjoint by
        // construction.** The leftover counts this sample's reads on sequences selection
        // dropped — a scan of the same support rows, asked by the step that dropped them. A
        // partial observation never reached the allele table at all: it says the sample carries
        // *at least* this much and is scored on its own axis. And a read removed as evidence —
        // named at some of this sample's records inside the locus and not at all of them —
        // reaches no `supported` row either, which is why the merge calls it lost depth and
        // counts it rather than leaving it to be inferred from an absence.
        //
        // **⚑ The third is a reading of spec §7, taken 2026-09-01, and the owner may overturn
        // it.** `DP` is *"every read observation the sample had at the locus, whether or not a
        // written allele explains them"*, and these reads were observed there: leaving them out
        // makes `DP` understate the depth at exactly the loci that span several of a sample's
        // records. `SampleEvidenceForOutput`'s own doc names two of the three, written before
        // this seam existed and calling its own shape provisional; `UnmatchedSupport`'s doc
        // excludes them from the *leftover* for a different reason — they carry no quality sum
        // and a read likelihood needs one — which settles the likelihood and not the file.
        // Disjointness is what makes adding them safe: the leftover is computed from
        // `supported`, and nothing these reads showed reaches it.
        let dropped = unmatched
            .get(covering)
            .map_or(0, |leftover| leftover.num_reads);
        let partial = support.partials.iter().fold(0_u32, |total, partial| {
            total.saturating_add(partial.num_reads)
        });
        let partial = partial.saturating_add(support.reads_removed_as_evidence);
        // **Added into the row rather than assigned over it**, so this count accumulates the
        // way `allele_reads` above does. The merge holds one entry per covering sample, so
        // nothing reaches a row twice; an assignment would nonetheless keep only the last of
        // two, which is a different answer from the one the alleles beside it would give.
        sample.reads_no_written_allele_explains = sample
            .reads_no_written_allele_explains
            .saturating_add(dropped.saturating_add(partial));
    }

    let (corrected_site_quality, artifact_penalties) = match locus.artifact_test_counts() {
        Some(counts) => {
            let (quality, penalties) =
                correct_site_quality(locus.uncorrected_site_quality(), &counts);
            (quality, Some(penalties))
        }
        // **A locus that gave the two tests nothing to weigh keeps its baseline**, and the
        // record says the tests did not run rather than writing two zeroed penalties, which a
        // reader could not tell from two tests that charged nothing.
        None => (locus.uncorrected_site_quality(), None),
    };

    LocusEvidenceForOutput {
        samples,
        allele_mapq,
        padding_base,
        corrected_site_quality,
        artifact_penalties,
        // **Read off the called locus's own candidate table**, which is where selection
        // stamped the kind: a repeat tract's table is `LocusKind::Ssr` and carries the motif,
        // and everything the record says about the repeat — the `STR` flag, `RU`, `PERIOD`
        // and each called allele's `REPCN` — is written from that one motif.
        repeat_tract: match locus.alleles().kind() {
            LocusKind::Ssr(detail) => Some(TractAnnotation::new(detail.motif)),
            LocusKind::Generic | LocusKind::SsrBundle => None,
        },
        filter: filter_for(locus, selection),
    }
}

/// **The one verdict this record is written on**, from the loop's answer and selection's.
///
/// # A loop that did not settle outranks everything selection could say
///
/// Not a preference: [`assemble_record`](crate::vcf::assemble::assemble_record) asserts that
/// the filter is `EMNoConv` exactly when the loop failed to converge, so the two cannot be
/// reported together and the loop's answer is the one the file has to carry. A tract that was
/// truncated *and* did not converge says `EMNoConv`, and the truncation is visible in the record
/// anyway — the alternatives it kept are the alternatives it kept.
///
/// # The tract verdicts are a tract's, and one of them ng never mints
///
/// `notPeriodic` and `tooManyAlleles` come from repeat-tract selection
/// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §6, §7). **`lowDepth` is declared and never
/// written**: production refuses a tract whose *cohort-summed* depth is under ten, and ng does
/// not port that gate — depth is asked once, upstream, per sample, by the merge's keep rule,
/// and there is no depth verdict on this path (§6). The vocabulary stays in the header because
/// a file written by an older caller can carry it.
///
/// # A truncated SNP/indel locus is not filtered, and that is unchanged
///
/// The ordinary path's cap cuts the lowest-ranked alternatives and calls the locus over the
/// rest; `tooManyAlleles` is spec §8's *tract* filter and stays one. What changed here is only
/// that a tract's truncation now reaches the column.
fn filter_for(locus: &LocusInference, selection: SelectionVerdict) -> FilterVerdict {
    if !locus.converged {
        return FilterVerdict::EmDidNotConverge;
    }
    match locus.alleles().kind() {
        LocusKind::Generic | LocusKind::SsrBundle => FilterVerdict::Pass,
        LocusKind::Ssr(_) => match selection {
            SelectionVerdict::NotPeriodic => FilterVerdict::NotPeriodic,
            SelectionVerdict::Truncated { .. } => FilterVerdict::TooManyAlleles,
            _ => FilterVerdict::Pass,
        },
    }
}

#[cfg(test)]
mod tests;
