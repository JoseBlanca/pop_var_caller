//! One locus's reads per sample, as the calling loop takes them: [`LocusEvidence`] and the
//! shapes it holds, [`GenericLocusSample`] and [`SsrLocusShape`].

use crate::calling::{
    CandidateAlleles, FrozenParameters, GenericSampleEvidence, ReferenceBesideLocus,
    SsrSampleEvidence,
};
use crate::locus_generation::{LocusKind, SsrBundleDetail, SsrDetail};
use crate::region_typing::run_parse::RunParse;
use crate::types::GenomeRegion;
use std::num::NonZeroU32;

/// One run sample's evidence at one SNP/indel locus, together with the candidate step's
/// ruling on whether the sample can be genotyped there at all.
///
/// **Two facts that must travel together, because separating them is a silent join.** The
/// evidence says what the sample's reads showed; the ruling says whether the alleles those
/// reads earned survived candidate selection's cap. A run keeps both per sample and in the
/// run's own sample order, so holding them in one entry is what stops a later step pairing
/// sample *i*'s reads with sample *j*'s ruling.
///
/// **Why the ruling is not on [`GenericSampleEvidence`] itself**, where
/// `doc/devel/ng/arch/read_likelihoods.md` §2.1 sketches it. A sample the cap ruled
/// uncallable never reaches the read likelihood: it leaves the loop before the first pass
/// and is scored against no genotype at all
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0). A field on the evidence view would
/// therefore be one the row builder is handed and never reads, and the shipped
/// `GenericSampleEvidence` does not carry it.
#[derive(Copy, Clone, Debug)]
pub struct GenericLocusSample<'a> {
    /// What this sample's reads showed, as the SNP/indel row consumes it. A sample with no
    /// coverage at the locus gets [`GenericSampleEvidence::empty`], which scores every
    /// genotype alike and leaves the prior to decide — the right answer rather than a
    /// special case (`doc/devel/ng/spec/calling_em_loop.md` §7).
    pub evidence: GenericSampleEvidence<'a>,
    /// **Whether the allele cap cut a sequence this sample's own reads had earned**, copied
    /// from `UnmatchedSupport::genotype_must_be_missing`
    /// ([`allele_candidates::UnmatchedSupport::genotype_must_be_missing`]).
    ///
    /// True means the locus is no longer called over something this sample carries, so
    /// every genotype the caller could form for it is wrong. Such a sample is set aside
    /// before the first pass: it is scored against nothing, contributes nothing to the
    /// cohort's expected copies, and its call is emitted as
    /// [`SampleGenotypeCall::Missing`] (`doc/devel/ng/spec/calling_em_loop.md` §5.0).
    ///
    /// **Nothing downstream could derive it.** The pooled error mass of the cut sequences
    /// is the same number under every genotype and cancels, so the read likelihood cannot
    /// tell a sample that lost a real allele from one that had a handful of error reads
    /// dropped — and the second is nearly every sample at nearly every locus.
    ///
    /// [`allele_candidates::UnmatchedSupport::genotype_must_be_missing`]: super::allele_candidates::UnmatchedSupport::genotype_must_be_missing
    /// [`SampleGenotypeCall::Missing`]: super::SampleGenotypeCall::Missing
    pub genotype_must_be_missing: bool,
}

impl GenericLocusSample<'_> {
    /// Whether the locus can be called on this sample at all — the negation of
    /// [`genotype_must_be_missing`](Self::genotype_must_be_missing), and **the one spelling
    /// of it**.
    ///
    /// The predicate is asked in three places — when the scratch's rows are counted, when they
    /// are claimed, and when the artifact summary decides whose reads to pool — and a fourth
    /// site that wrote it out again would be invisible to anyone changing what *callable*
    /// means.
    #[inline]
    #[must_use]
    pub fn is_callable(&self) -> bool {
        !self.genotype_must_be_missing
    }
}

/// One locus's evidence, per sample, in whichever shape its path uses.
///
/// **This enum is the only place the two calling paths meet.** Below it everything is
/// path-pure: the discriminant chooses the row builder, and it is the same discriminant
/// that chose the candidates (`doc/devel/ng/arch/calling_em_loop.md` §2).
///
/// **Every variant's `per_sample` is one entry per sample of the run, in the run's sample
/// order** — not one entry per *covering* sample. That is the order every per-sample slice
/// in [`FrozenParameters`] is indexed by and the order [`LocusInference::per_sample`] comes
/// back in, and holding one order for the whole run is what makes the loop's fixed-order
/// sum over samples reproducible (`doc/devel/ng/spec/calling_em_loop.md` §8). The merge's
/// own `CohortObservation::per_sample` holds only the covering samples and is a different
/// list; converting between them is the input edge's work, not this type's.
///
/// [`LocusInference::per_sample`]: super::LocusInference::per_sample
#[derive(Copy, Clone, Debug)]
pub enum LocusEvidence<'a> {
    /// A SNP/indel locus: what each sample showed, and which samples the candidate cap
    /// already ruled uncallable.
    Generic {
        /// Where on the reference this locus is.
        region: GenomeRegion,
        /// One entry per run sample, in run order.
        per_sample: &'a [GenericLocusSample<'a>],
        /// **The reference on either side of the locus**, as much as its partial reads need to
        /// be compared against an allele followed by what lies beyond it
        /// (`doc/devel/ng/spec/read_likelihoods.md` §5.3). Empty unless some partial showed
        /// more bases than the shortest allele has.
        reference_beside: ReferenceBesideLocus<'a>,
    },
    /// An STR segment, or an STR bundle: what each sample showed, and the repeat ground it is
    /// scored against.
    ///
    /// **No sample is set aside here**, and the absence is structural rather than an
    /// oversight. Allele discovery, a pre-pass inside candidate selection, can admit a length
    /// one sample's own top-`ploidy` cut left out, so a sample that lost one is not locked out
    /// of the locus (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1).
    Ssr {
        /// Where on the reference this locus is.
        region: GenomeRegion,
        /// One entry per run sample, in run order.
        per_sample: &'a [SsrSampleEvidence<'a>],
        /// The repeat ground — a segment's motif and flanks with each candidate's repeat count,
        /// or a bundle's runs and flanks with each candidate's run parse.
        shape: SsrLocusShape<'a>,
    },
}

/// **The repeat ground an STR locus's reads are scored against, and what each candidate carries
/// on it** — one repeat count per candidate at an STR segment, one run parse per candidate at an
/// STR bundle (`str_bundle_loci.md` §7).
#[derive(Copy, Clone, Debug)]
pub enum SsrLocusShape<'a> {
    /// An STR segment: one motif, two flanks.
    Segment {
        /// The tract's repeat unit and its two flanks — what the repeat read model scores
        /// against.
        detail: &'a SsrDetail,
        /// **How many whole repeats each candidate carries**, parallel to the locus's
        /// [`CandidateAlleles`] and in the same order.
        ///
        /// **It is not derivable from the bases, which is the whole reason it travels.** An
        /// interrupted tract — one whose repeat is broken by a substitution — holds fewer whole
        /// repeats than its length suggests, and the slippage fit is keyed by the count rather
        /// than by the length (`doc/devel/ng/spec/read_likelihoods.md` §4.4). Two candidates of
        /// equal length can therefore be drawn from different strata.
        candidate_repeat_counts: &'a [NonZeroU32],
    },
    /// An STR bundle: its reference runs and the hull's flanks.
    Bundle {
        /// The hull's flanks and its reference runs.
        detail: &'a SsrBundleDetail,
        /// **Each candidate read against the reference runs**, parallel to the locus's
        /// [`CandidateAlleles`]: its count in every run, which keys that run's stutter and prior,
        /// and each run's place in it, where a slip in that run is placed.
        candidate_parses: &'a [RunParse],
    },
}

impl SsrLocusShape<'_> {
    /// How many candidates the shape describes — one count, or one parse, per candidate.
    #[must_use]
    pub fn candidate_count(&self) -> usize {
        match self {
            Self::Segment {
                candidate_repeat_counts,
                ..
            } => candidate_repeat_counts.len(),
            Self::Bundle {
                candidate_parses, ..
            } => candidate_parses.len(),
        }
    }
}

impl<'a> LocusEvidence<'a> {
    /// A SNP/indel locus's evidence.
    ///
    /// # Panics
    ///
    /// If `per_sample` names no sample. A run has at least one sample and every sample gets
    /// an entry — an empty one is a locus whose evidence was lost on the way in, not a
    /// locus nobody covered, which is one entry per sample of
    /// [`GenericSampleEvidence::empty`].
    #[must_use]
    pub fn generic(region: GenomeRegion, per_sample: &'a [GenericLocusSample<'a>]) -> Self {
        Self::generic_beside_reference(region, per_sample, ReferenceBesideLocus::NONE)
    }

    /// A SNP/indel locus's evidence **with the reference beside it** — what the run hands the
    /// loop, so that a read that ran out inside the locus is compared against an allele
    /// followed by the reference beyond it (`doc/devel/ng/spec/read_likelihoods.md` §5.3).
    /// [`Self::generic`] is this with no flank.
    ///
    /// # Panics
    ///
    /// As [`Self::generic`].
    #[must_use]
    pub fn generic_beside_reference(
        region: GenomeRegion,
        per_sample: &'a [GenericLocusSample<'a>],
        reference_beside: ReferenceBesideLocus<'a>,
    ) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the SNP/indel evidence at {region} names no sample: a locus carries one entry \
             per sample of the run and a run has at least one sample, so an empty list is \
             evidence that went missing rather than a locus nobody covered — a sample with \
             no reads gets GenericSampleEvidence::empty()"
        );
        Self::Generic {
            region,
            per_sample,
            reference_beside,
        }
    }

    /// A repeat tract's or repeat bundle's evidence.
    ///
    /// # Panics
    ///
    /// If `per_sample` names no sample, for the reason [`Self::generic`] gives.
    #[must_use]
    pub fn ssr(
        region: GenomeRegion,
        per_sample: &'a [SsrSampleEvidence<'a>],
        detail: &'a SsrDetail,
        candidate_repeat_counts: &'a [NonZeroU32],
    ) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the repeat-tract evidence at {region} names no sample: a locus carries one \
             entry per sample of the run and a run has at least one sample, so an empty \
             list is evidence that went missing rather than a locus nobody covered"
        );
        assert!(
            !candidate_repeat_counts.is_empty(),
            "the repeat tract at {region} supplied no candidate repeat counts: a tract is \
             called over at least its reference length, and the count is what keys the \
             slippage fit. Their number is checked against the candidate table by \
             `assert_matches_locus_and_run`, which is where the two meet"
        );
        Self::Ssr {
            region,
            per_sample,
            shape: SsrLocusShape::Segment {
                detail,
                candidate_repeat_counts,
            },
        }
    }

    /// An STR bundle's evidence (`str_bundle_loci.md` §7): each candidate's run parse rides
    /// beside the samples, as a segment's repeat counts do.
    ///
    /// # Panics
    ///
    /// If `per_sample` names no sample, for the reason [`Self::generic`] gives, or no candidate
    /// was parsed.
    #[must_use]
    pub fn ssr_bundle(
        region: GenomeRegion,
        per_sample: &'a [SsrSampleEvidence<'a>],
        detail: &'a SsrBundleDetail,
        candidate_parses: &'a [RunParse],
    ) -> Self {
        assert!(
            !per_sample.is_empty(),
            "the STR bundle evidence at {region} names no sample: a locus carries one entry per \
             sample of the run"
        );
        assert!(
            !candidate_parses.is_empty(),
            "the STR bundle at {region} supplied no candidate parse: a bundle is called over at \
             least its reference allele"
        );
        Self::Ssr {
            region,
            per_sample,
            shape: SsrLocusShape::Bundle {
                detail,
                candidate_parses,
            },
        }
    }

    /// Where on the reference this locus is.
    #[inline]
    #[must_use]
    pub fn region(&self) -> GenomeRegion {
        match self {
            Self::Generic { region, .. } | Self::Ssr { region, .. } => *region,
        }
    }

    /// How many samples the run has — the length of every per-sample list at this locus.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        match self {
            Self::Generic { per_sample, .. } => per_sample.len(),
            Self::Ssr { per_sample, .. } => per_sample.len(),
        }
    }

    /// **How many of the run's samples this locus can actually be called on** — the rest were
    /// ruled uncallable by the candidate step, for having earned a sequence the allele cap cut
    /// (`doc/devel/ng/spec/candidate_alleles.md` §4.1).
    ///
    /// **It counts the run's samples and not the locus's covering ones**, which is the whole of
    /// what a caller has to get right here. A sample that covered nothing is **callable**: its
    /// evidence is empty, an empty sum is zero, every genotype scores alike and the prior
    /// decides alone. So this reaches zero only where **every sample of the run** covered the
    /// locus and every one of them lost a sequence its own reads had earned — which in a large
    /// cohort at low coverage is close to never, because some sample almost always covers
    /// nothing.
    ///
    /// **A locus can still have none.** A run must count such a locus and carry on rather than
    /// die on it (owner's ruling, 2026-09-01), so whoever drives the calling asks this before
    /// handing the locus to a genotyper — and
    /// [`LocusGenotyper::call_locus`] keeps its
    /// precondition that there is somebody to call, since its scratch cannot be prepared for
    /// no rows.
    ///
    /// **A repeat tract sets no sample aside**, so on that path this is always the run's whole
    /// sample count (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1).
    ///
    /// [`LocusGenotyper::call_locus`]: super::inference::LocusGenotyper::call_locus
    #[inline]
    #[must_use]
    pub fn callable_sample_count(&self) -> usize {
        match self {
            Self::Generic { per_sample, .. } => per_sample
                .iter()
                .filter(|sample| sample.is_callable())
                .count(),
            Self::Ssr { per_sample, .. } => per_sample.len(),
        }
    }

    /// Which calling path this evidence is on, in the words the panic messages use.
    #[inline]
    fn path_word(&self) -> &'static str {
        match self {
            Self::Generic { .. } => "SNP/indel",
            Self::Ssr { .. } => "repeat-tract",
        }
    }

    /// Check this evidence against the locus it claims to describe and the run it claims to
    /// come from — **the ordering contract, made a runtime fact rather than a promise.**
    ///
    /// Three things are checked, and each is a caller bug whose symptom is a wrong genotype
    /// rather than a crash:
    ///
    /// - **the path.** A repeat tract's evidence handed to a locus whose alleles say
    ///   SNP/indel would be scored by the wrong read model, which is a different likelihood
    ///   at every sample with nothing failing
    ///   (`doc/devel/ng/arch/calling_em_loop.md` §2). A repeat **bundle** goes to the
    ///   repeat path with a tract, which is how every other consumer of
    ///   [`LocusKind`] groups the two, and is scored over its runs
    ///   (`str_bundle_loci.md` §7; `bundle_candidates`). The one-tract candidate table,
    ///   `tract_candidates`, still refuses one by name.
    /// - **the cohort.** One run-wide sample order indexes this evidence, every per-sample
    ///   slice of `parameters`, and the calls that come back. Two lists of different lengths
    ///   are two different orders, and a positional join between them silently pairs one
    ///   sample's reads with another's inbreeding coefficient.
    /// - **the repeat counts, at a tract only.** They run parallel to the candidate table,
    ///   and a table with more candidates than counts would pair one candidate's bases with
    ///   another candidate's stratum.
    ///
    /// **The locus's genotype table is not checked here**, because this method never sees
    /// it. That check has its own home, at the one point the shape is fixed:
    /// [`CallingScratch::prepare_for_locus`].
    ///
    /// # Panics
    ///
    /// On any of the three disagreements, in release as well as debug: each is a comparison
    /// of two integers or two discriminants, against a defect that would otherwise reach the
    /// output as a genotype.
    ///
    /// [`CallingScratch::prepare_for_locus`]: super::CallingScratch::prepare_for_locus
    pub fn assert_matches_locus_and_run(
        &self,
        alleles: &CandidateAlleles,
        parameters: &FrozenParameters<'_>,
    ) {
        let paths_agree = matches!(
            (self, alleles.kind()),
            (Self::Generic { .. }, LocusKind::Generic)
                | (
                    Self::Ssr {
                        shape: SsrLocusShape::Segment { .. },
                        ..
                    },
                    LocusKind::Ssr(_)
                )
                | (
                    Self::Ssr {
                        shape: SsrLocusShape::Bundle { .. },
                        ..
                    },
                    LocusKind::SsrBundle(_)
                )
        );
        assert!(
            paths_agree,
            "the evidence at {} is on the {} path and its allele table is a {} locus, so \
             the two belong to different loci — the discriminant that chose the candidates \
             is the one that chooses the read model",
            self.region(),
            self.path_word(),
            locus_kind_word(alleles.kind())
        );
        // **The join a repeat tract adds, and it is the one no type enforces.** Its repeat
        // counts run
        // parallel to the candidate table, and a table with more candidates than counts would
        // pair one candidate's bases with another's stratum — a wrong slippage model per
        // candidate, with a well-formed genotype coming back.
        if let Self::Ssr { shape, .. } = self {
            assert_eq!(
                shape.candidate_count(),
                alleles.len(),
                "the repeat tract at {} is called over {} candidates and {} repeat counts \
                 arrived, so one of the two belongs to a different locus",
                self.region(),
                alleles.len(),
                shape.candidate_count()
            );
        }
        assert_eq!(
            self.sample_count(),
            parameters.sample_count(),
            "the evidence at {} covers {} samples and the run's frozen parameters cover {}, \
             so one of them is indexed by a different sample order",
            self.region(),
            self.sample_count(),
            parameters.sample_count()
        );
    }
}

/// What kind of locus this is, in one word, for a panic message.
///
/// **Not `{:?}` on the kind**, which for a repeat tract prints both flanks as decimal byte
/// arrays — a screenful of digits for a fact the discriminant alone carries.
pub(super) fn locus_kind_word(kind: &LocusKind) -> &'static str {
    match kind {
        LocusKind::Generic => "SNP/indel",
        LocusKind::Ssr(_) => "repeat-tract",
        LocusKind::SsrBundle(_) => "repeat-bundle",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calling::likelihood::TEST_TRACT_LEN;
    use crate::calling::tests::{
        NO_SUBSTITUTION_RATES, frozen_parameters, no_strata, one_batch, one_read_group,
        outbred_samples, region, str_two_allele_locus, two_allele_locus,
    };
    use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
    use crate::types::Motif;

    /// Two candidates' repeat counts, **different from each other**, for a tract fixture whose
    /// candidate table holds two lengths. Distinct so that a lookup reading one candidate's
    /// count for another is a different stratum rather than the same one.
    fn two_repeat_counts() -> Vec<NonZeroU32> {
        vec![
            NonZeroU32::new(6).expect("six repeats"),
            NonZeroU32::new(7).expect("seven repeats"),
        ]
    }

    fn ssr_detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }
    }

    /// Evidence and the allele table have to be on the same path: the discriminant that
    /// chose the candidates is the one that chooses the read model.
    ///
    /// **Handing a repeat tract's alleles a SNP/indel's evidence would not crash** — the
    /// generic row would score the tract's sequences as though they were substitutions,
    /// giving a different likelihood at every sample and a plausible genotype at the end of
    /// it.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn evidence_on_the_wrong_path_for_its_allele_table_is_refused() {
        let (tract_alleles, _) = str_two_allele_locus();
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&tract_alleles, &parameters);
    }

    /// The evidence and the run's per-sample parameters are indexed by one sample order, so
    /// two lists of different lengths are two different orders.
    #[test]
    #[should_panic(expected = "indexed by a different sample order")]
    fn evidence_covering_a_different_sample_count_from_the_run_is_refused() {
        let (alleles, _) = two_allele_locus();
        let calibration = one_read_group();
        let inbreeding = outbred_samples(3);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&alleles, &parameters);
    }

    /// Both paths pass their own check, and each carries the locus's region and the run's
    /// sample count.
    ///
    /// **A repeat tract's evidence carries no missing-genotype flag at all**, and the
    /// absence is the design: a tract's discovery round can put back a length the cap cut,
    /// so no sample is set aside there.
    #[test]
    fn evidence_that_matches_its_locus_and_its_run_is_accepted_on_both_paths() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(2);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let (generic_alleles, _) = two_allele_locus();
        let generic_samples = [
            GenericLocusSample {
                evidence: GenericSampleEvidence::empty(),
                genotype_must_be_missing: false,
            },
            GenericLocusSample {
                evidence: GenericSampleEvidence::empty(),
                genotype_must_be_missing: true,
            },
        ];
        let generic = LocusEvidence::generic(region(), &generic_samples);
        generic.assert_matches_locus_and_run(&generic_alleles, &parameters);
        assert_eq!(generic.region(), region());
        assert_eq!(generic.sample_count(), 2);
        match generic {
            LocusEvidence::Generic { per_sample, .. } => {
                assert!(!per_sample[0].genotype_must_be_missing);
                assert!(per_sample[1].genotype_must_be_missing);
            }
            LocusEvidence::Ssr { .. } => panic!("built on the SNP/indel path"),
        }

        let (tract_alleles, _) = str_two_allele_locus();
        let detail = ssr_detail();
        let tract_samples = [
            SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN),
            SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN),
        ];
        let counts = two_repeat_counts();
        let tract = LocusEvidence::ssr(region(), &tract_samples, &detail, &counts);
        tract.assert_matches_locus_and_run(&tract_alleles, &parameters);
        assert_eq!(tract.region(), region());
        assert_eq!(tract.sample_count(), 2);
    }

    /// A locus whose evidence list is empty is evidence that went missing, not a locus
    /// nobody covered — a sample with no reads gets an empty *entry*, and the prior decides
    /// its genotype alone.
    #[test]
    #[should_panic(expected = "evidence that went missing")]
    fn evidence_naming_no_sample_at_all_is_refused() {
        let _ = LocusEvidence::generic(region(), &[]);
    }

    /// An STR **bundle**'s allele table takes the repeat path's evidence in its bundle shape —
    /// the runs and each candidate's run parse (`str_bundle_loci.md` §7).
    ///
    /// **No fixture reached this cell before, and two mutations lived there.** Deleting the
    /// bundle arm from the accepting match left the whole suite green, and so did widening the
    /// SNP/indel arm to accept a bundle — the second of which routes a repeat bundle to the
    /// SNP/indel row: a plausible genotype at every sample, nothing failing.
    #[test]
    fn ssr_evidence_against_a_bundle_allele_table_is_accepted() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let bundle = CandidateAlleles::new(
            Box::from(b"CAGCAG".as_slice()),
            LocusKind::SsrBundle(crate::locus_generation::SsrBundleDetail::default()),
        );
        let detail = crate::locus_generation::SsrBundleDetail::default();
        let per_sample = [SsrSampleEvidence::of_bundle(&[], TEST_TRACT_LEN)];
        // One parse, because a bundle built from a reference alone is called over one allele.
        let parses = [RunParse::default()];
        let evidence = LocusEvidence::ssr_bundle(region(), &per_sample, &detail, &parses);
        evidence.assert_matches_locus_and_run(&bundle, &parameters);
    }

    /// **A segment's evidence is refused against a bundle's table**, and a bundle's against a
    /// segment's: each shape carries what its own read model scores, a repeat count per candidate
    /// or a run parse, and the other model cannot read it.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn segment_evidence_against_a_bundle_allele_table_is_refused() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );
        let bundle = CandidateAlleles::new(
            Box::from(b"CAGCAG".as_slice()),
            LocusKind::SsrBundle(crate::locus_generation::SsrBundleDetail::default()),
        );
        let detail = ssr_detail();
        let per_sample = [SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN)];
        let counts = [NonZeroU32::new(6).expect("six repeats")];
        let evidence = LocusEvidence::ssr(region(), &per_sample, &detail, &counts);
        evidence.assert_matches_locus_and_run(&bundle, &parameters);
    }

    /// The other half of the same cell: SNP/indel evidence at a repeat bundle is refused.
    ///
    /// This is the mutation that bites — an accepting arm widened to admit a bundle scores
    /// a repeat locus with the SNP/indel read model, and nothing anywhere fails.
    #[test]
    #[should_panic(expected = "belong to different loci")]
    fn generic_evidence_against_a_bundle_allele_table_is_refused() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let no_declared_batching = one_batch(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &no_declared_batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let bundle = CandidateAlleles::new(
            Box::from(b"CAGCAG".as_slice()),
            LocusKind::SsrBundle(crate::locus_generation::SsrBundleDetail::default()),
        );
        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::empty(),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(region(), &per_sample);
        evidence.assert_matches_locus_and_run(&bundle, &parameters);
    }

    /// The repeat constructor's own empty-list guard, which its SNP/indel twin's test does
    /// not reach.
    ///
    /// **Measured: this was the one assertion in the module no test noticed.** Downgrading
    /// all sixteen of the implementation's release-held checks to `debug_assert!` and
    /// running under `--release` failed a test for fifteen of them; this was the sixteenth.
    /// The two constructors are a copy-paste pair, which is the shape in which one half
    /// drifts from the other.
    #[test]
    #[should_panic(expected = "evidence that went missing")]
    fn ssr_evidence_naming_no_sample_at_all_is_refused() {
        let detail = ssr_detail();
        let _ = LocusEvidence::ssr(region(), &[], &detail, &two_repeat_counts());
    }

    /// **A tract with no candidate repeat counts at all is refused**, and it is a different
    /// failure from the empty sample list beside it: a tract is called over at least its
    /// reference length, so a run that supplied none lost the counts on the way in.
    ///
    /// **Refused here rather than left to the parameter assembly**, which would otherwise ask
    /// for a table of `read groups × 0` cells and refuse it naming the table rather than the
    /// locus.
    #[test]
    #[should_panic(expected = "supplied no candidate repeat counts")]
    fn ssr_evidence_with_no_candidate_repeat_counts_is_refused() {
        let detail = ssr_detail();
        let per_sample = [SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN)];
        let _ = LocusEvidence::ssr(region(), &per_sample, &detail, &[]);
    }

    /// **Repeat counts that are not one per candidate are refused where the two still name each
    /// other**, rather than at whichever genotype first read the wrong cell.
    ///
    /// The counts run parallel to the candidate table, so a table with more candidates than
    /// counts would pair one candidate's bases with another candidate's stratum — a wrong
    /// slippage model per candidate, with a well-formed genotype coming back and nothing
    /// failing.
    #[test]
    #[should_panic(expected = "called over 2 candidates and 1 repeat counts")]
    fn a_tract_whose_repeat_counts_are_not_one_per_candidate_is_refused() {
        let calibration = one_read_group();
        let inbreeding = outbred_samples(1);
        let strata = no_strata();
        let batching = SequencingBatches::all_together_over(1, 1);
        let parameters = frozen_parameters(
            &calibration,
            &[],
            &batching,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
        );

        let detail = ssr_detail();
        let mut alleles = CandidateAlleles::new(
            Box::from(b"ATATAT".as_slice()),
            LocusKind::Ssr(ssr_detail()),
        );
        alleles.admit(Box::from(b"ATATATAT".as_slice()));
        let per_sample = [SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN)];
        // Two candidates, one count.
        let counts = [NonZeroU32::new(6).expect("six repeats")];
        let evidence = LocusEvidence::ssr(region(), &per_sample, &detail, &counts);
        evidence.assert_matches_locus_and_run(&alleles, &parameters);
    }
}
