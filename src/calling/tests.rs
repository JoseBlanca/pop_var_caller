//! Test fixtures used by the tests of more than one type family in `calling`, and by
//! `inference`'s tests (`one_batch`).

use super::*;
use crate::calling::genotype_prior::SpectrumSeed;
use crate::locus_generation::{LocusKind, SsrDetail};
use crate::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::types::{ContigId, GenomeRegion, InbreedingF, Motif, Ploidy, Position};

pub(super) fn generic_reference() -> CandidateAlleles {
    CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic)
}

/// A run whose repeat-tract substitution rates were never fitted — the empty map, which is
/// what a fixture calling no tract needs, and what
/// `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
///
/// A `static` rather than a function, so that a call site can borrow it for as long as the
/// parameters live: `BTreeMap::new` is a `const fn`, and a temporary would be freed at the
/// end of the statement that built the view.
pub(super) static NO_SUBSTITUTION_RATES: std::collections::BTreeMap<
    crate::parameter_estimation::repeat_strata::StratumKey,
    crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
> = std::collections::BTreeMap::new();

pub(super) fn two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
    let mut alleles = generic_reference();
    alleles.admit(Box::from(b"T".as_slice()));
    let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
    (alleles, copies)
}

/// A repeat tract with two candidate lengths — the path the gene-diversity seed
/// marker belongs to, and the only one on which a locus may carry it.
pub(super) fn str_two_allele_locus() -> (CandidateAlleles, ExpectedAlleleCopies) {
    let mut alleles = CandidateAlleles::new(
        Box::from(b"ATAT".as_slice()),
        LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }),
    );
    alleles.admit(Box::from(b"ATATAT".as_slice()));
    let copies = ExpectedAlleleCopies::new(vec![2.6, 1.4], &alleles);
    (alleles, copies)
}

pub(super) fn region() -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(3),
        start: Position(940),
        end: Position(940),
    }
}

use std::collections::BTreeMap;

pub(super) fn no_strata() -> StratumFits {
    StratumFits::over(&[], BTreeMap::new())
}

pub(super) fn diploid_ploidy() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

pub(super) fn outbred_samples(samples: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("a legal coefficient"); samples]
}

pub(super) fn one_read_group() -> Vec<ReadGroupCalibration> {
    vec![ReadGroupCalibration::unscaled()]
}

pub(super) fn neutral_seed() -> SpectrumSeed {
    SpectrumSeed::new(1.0, 1e-3, genotype_prior::SeedRegime::NeutralShape)
}

/// The frozen parameters for a run of diploids, routed to whichever constructor the
/// contamination list calls for — an empty one is the run where the fit identified
/// none, which has its own named door.
///
/// **The contaminated branch takes the batching as an argument**, because it is borrowed
/// for as long as the parameters are and a helper that built one inside would hand back a
/// borrow of something already dropped. A run with no mixture has nothing to read it, so
/// the uncontaminated branch never looks at it.
pub(super) fn frozen_parameters<'a>(
    calibration: &'a [ReadGroupCalibration],
    contamination: &'a [ContaminationView],
    batching: &'a SequencingBatches,
    inbreeding: &'a [InbreedingF],
    strata: &'a StratumFits,
    substitution: &'a std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
) -> FrozenParameters<'a> {
    if contamination.is_empty() {
        FrozenParameters::uncontaminated(
            calibration,
            inbreeding,
            neutral_seed(),
            strata,
            substitution,
            diploid_ploidy(),
        )
    } else {
        FrozenParameters::new(
            calibration,
            contamination,
            batching,
            inbreeding,
            neutral_seed(),
            strata,
            substitution,
            diploid_ploidy(),
        )
    }
}

/// The default batching — one batch holding a run of `libraries` read groups over
/// `samples` samples, which is what a run that declared nothing gets.
pub(crate) fn one_batch(libraries: usize, samples: usize) -> SequencingBatches {
    let names: Vec<(String, String)> = (0..libraries)
        .map(|library| {
            (
                format!("rg{library}"),
                format!("s{}", library.min(samples.saturating_sub(1))),
            )
        })
        .collect();
    let borrowed: Vec<(&str, &str)> = names
        .iter()
        .map(|(id, sample)| (id.as_str(), sample.as_str()))
        .collect();
    let groups = crate::read::input::read_groups::ReadGroups::of_libraries(&borrowed);
    assert_eq!(
        groups.read_groups_per_sample().len(),
        samples,
        "the fixture asked for {samples} samples over {libraries} libraries"
    );
    SequencingBatches::all_together(&groups)
}
