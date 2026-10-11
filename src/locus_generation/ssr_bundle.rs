//! The STR bundle generator — one STR bundle → one locus over its hull
//! (`doc/devel/ng/spec/str_bundle_loci.md` §6.1).
//!
//! **It is the STR segment's generator run over the hull** ([`SsrGenerator`]): the same read
//! fetch and depth cap, the same per-read measurement — the spelling from the read's own
//! alignment, the flank deletion at the hull's two outer edges in its smallest form
//! (`str_edge_deletion_alleles.md` §4) — and the same per-sample tally keyed on spelling, flank
//! deletion and read group. Three things differ, and each is a field of the
//! [`MeasuredRepeat`] the hull is measured as:
//!
//! - **the delimiter** is the flat-gap one ([`SsrFlatGapAligner`]), which needs only the hull's
//!   boundaries. The segments' default prices whole-unit slips of one motif and cannot
//!   describe a hull of several. It is built **with the segments' firm-anchor test**
//!   ([`SsrFlatGapAligner::with_anchor_rule`]), so a read that ends inside the hull is not
//!   reported complete over part of it (found in review, 2026-10-11);
//! - **the edge-escape check** asks, at each edge, about the reference run reaching it, and
//!   about nothing where fixed sequence reaches it;
//! - **only reads that cross the whole hull are used**: a read that anchors one edge and ends
//!   inside is dropped and counted (`SsrGeneratorCounts::ended_inside_a_bundle`).
//!
//! What it adds is the hull's **reference runs**, found once per bundle from the hull's bases
//! alone, every motif of period 1 to 6 looked for ([`find_reference_runs`]), and carried on the
//! locus ([`SsrBundleDetail`]).

use super::ssr::{
    EdgeRunPeriods, MeasuredRepeat, SsrGenerator, SsrGeneratorConfig, SsrGeneratorConfigError,
    SsrGeneratorCounts, seed_for_repeat,
};
use super::{
    GeneratorCounts, LocusGenerationError, LocusGenerator, LocusKind, SampleLocusObservations,
    SsrBundleDetail,
};
use crate::alignment::PerQualityEmission;
use crate::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use crate::alignment::ssr_unit_robust::AnchorRule;
use crate::read::filtering::ReadFilterCounts;
use crate::read::input::SampleReads;
use crate::ref_seq::{ContigTable, EvictableRefSeq, RawRefSeq, RefSeq};
use crate::region_typing::BundleMember;
use crate::region_typing::reference_runs::{ReferenceRun, find_reference_runs};
use crate::region_typing::segment_criteria::MinCopies;
use crate::types::{Bp, GenomeRegion, Motif, ReadGroupId};

/// The STR bundle generator: turns one STR bundle into one locus over its hull.
///
/// Holds an [`SsrGenerator`] with the flat-gap delimiter, under the default [`AnchorRule`], and
/// runs it over each hull; the copy
/// floors are the typing's own (`--min-copies`), because a run is what the typing would call
/// long enough to stutter (`str_bundle_loci.md` §6.1, step 2).
pub struct SsrBundleGenerator<R: RawRefSeq + EvictableRefSeq + Send + 'static> {
    /// The segment generator this one is, run over a hull.
    hulls: SsrGenerator<R, SsrFlatGapAligner<PerQualityEmission>>,
    /// The fewest whole copies a stretch needs, per period, to be a run.
    min_copies: MinCopies,
}

impl<R> SsrBundleGenerator<R>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
{
    /// A bundle generator over `reference` (the hull's flanks and bases) and `make_reference`
    /// (the read cursor's factory), with `config` checked against `bundle_threshold` as the
    /// segment generator checks it — a hull's flanks are clean out to the same radius — and
    /// `min_copies` the typing's copy floors.
    pub fn new(
        reference: R,
        make_reference: impl FnMut() -> R + Send + 'static,
        config: SsrGeneratorConfig,
        bundle_threshold: Bp,
        min_copies: MinCopies,
    ) -> Result<Self, SsrGeneratorConfigError> {
        Ok(Self {
            hulls: SsrGenerator::new(
                reference,
                make_reference,
                SsrFlatGapAligner::with_anchor_rule(
                    PerQualityEmission::new(),
                    AnchorRule::default(),
                ),
                config,
                bundle_threshold,
            )?,
            min_copies,
        })
    }

    /// The running counts — accumulated across every bundle, readable at any point.
    pub fn counts(&self) -> &SsrGeneratorCounts {
        self.hulls.counts()
    }
}

impl<R> LocusGenerator<[BundleMember]> for SsrBundleGenerator<R>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
{
    fn begin_segment(&mut self, region: GenomeRegion) {
        LocusGenerator::<crate::region_typing::segment_criteria::SsrSegment>::begin_segment(
            &mut self.hulls,
            region,
        );
    }

    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        Some(GeneratorCounts::Ssr(self.hulls.counts()))
    }

    fn read_filter_counts(&self) -> Vec<(Option<ReadGroupId>, ReadFilterCounts)> {
        self.hulls.read_filter_counts()
    }

    fn next_locus(
        &mut self,
        members: &[BundleMember],
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        let Some(region) = self.hulls.claim_the_one_locus() else {
            return Ok(None);
        };
        // The flat-gap delimiter reads no motif; the geometry type carries one, so the first
        // member's goes there. A bundle has two members or more; an empty list is a caller's
        // mistake the locus survives.
        let delimiter_motif = members
            .first()
            .map(BundleMember::motif)
            .unwrap_or_else(|| Motif::new(b"N").expect("one base is a motif"));
        // The edges are not known until the runs are, and the runs need the hull's bases. The
        // fetch reads only the repeat's span, so the edges it is handed are a stand-in, replaced
        // from the runs found before any read is measured.
        let edges_the_fetch_ignores = EdgeRunPeriods {
            left: None,
            right: None,
        };
        let mut locus = self.hulls.fetch_locus(
            region,
            MeasuredRepeat::bundle_hull(
                region.start.get(),
                region.end.get(),
                delimiter_motif,
                edges_the_fetch_ignores,
            ),
        )?;
        let (_, hull, _) = locus.split_flanks();
        let runs = find_reference_runs(hull, &self.min_copies);
        locus.repeat = MeasuredRepeat::bundle_hull(
            region.start.get(),
            region.end.get(),
            delimiter_motif,
            edge_run_periods(&runs, hull.len() as u64),
        );

        let seed = seed_for_repeat(
            self.hulls
                .contig_name(region.contig)
                .expect("the hull's flanks were just fetched from this contig"),
            region.start.get(),
        );
        let measured = self.hulls.measure_reads(&locus, region, seed, reads)?;

        let (left_flank, reference_bases, right_flank) = locus.split_flanks();
        Ok(Some(SampleLocusObservations {
            region,
            reference_bases: reference_bases.into(),
            observations: measured.observations,
            reads_without_observation: measured.reads_without_observation,
            reads_discarded_by_cap: measured.reads_discarded_by_cap,
            kind: LocusKind::SsrBundle(SsrBundleDetail {
                left_flank: left_flank.into(),
                right_flank: right_flank.into(),
                runs: runs.into_boxed_slice(),
            }),
        }))
    }
}

/// The period of the run reaching each edge of a hull `hull_len` bases long, `None` where
/// fixed sequence reaches it.
fn edge_run_periods(runs: &[ReferenceRun], hull_len: u64) -> EdgeRunPeriods {
    EdgeRunPeriods {
        left: runs
            .first()
            .filter(|run| run.offset() == 0)
            .map(ReferenceRun::period),
        right: runs
            .last()
            .filter(|run| run.end() == hull_len)
            .map(ReferenceRun::period),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locus_generation::ReadWitness;
    use crate::read::filtering::ReadFilterConfig;
    use crate::read::input::test_fixtures::{
        fixture_reference, header, indexed_bam, matching_contigs, read_named_with_length,
    };
    use crate::ref_seq::InMemoryRefSeq;
    use crate::types::{ContigId, Position};

    /// The all-`A` fixture contigs (`chr1` = 100 bases), as the reads are opened against.
    fn fixture_ref_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_named_contigs(
            matching_contigs()
                .iter()
                .map(|(name, len, _)| ((*name).to_string(), vec![b'A'; *len]))
                .collect(),
        )
    }

    fn sample_reads_with(
        records: &[noodles_sam::alignment::RecordBuf],
    ) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    fn generator() -> SsrBundleGenerator<InMemoryRefSeq> {
        SsrBundleGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
            MinCopies::default(),
        )
        .expect("flank within the bundle threshold")
    }

    fn hull(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn members() -> Vec<BundleMember> {
        vec![
            BundleMember::new(40, 44, Motif::new(b"A").unwrap()),
            BundleMember::new(46, 52, Motif::new(b"AAT").unwrap()),
        ]
    }

    /// **One locus per bundle, over the hull, carrying the hull's flanks and reference runs**,
    /// and one only. The fixture reference is all `A`, so the hull is one 13-base `A` run, from
    /// the hull's first base to its last.
    #[test]
    fn a_bundle_yields_one_locus_over_its_hull_with_its_runs_and_flanks() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[]);
        let mut bundles = generator();
        bundles.begin_segment(hull(40, 52));
        let locus = bundles
            .next_locus(&members(), &reads)
            .expect("no error")
            .expect("one locus per bundle, even at zero coverage");
        assert_eq!(locus.region, hull(40, 52));
        assert_eq!(&*locus.reference_bases, &[b'A'; 13][..]);
        assert!(locus.observations.is_empty());
        match &locus.kind {
            LocusKind::SsrBundle(detail) => {
                assert_eq!(&*detail.left_flank, &[b'A'; 10][..]);
                assert_eq!(&*detail.right_flank, &[b'A'; 10][..]);
                assert_eq!(
                    &*detail.runs,
                    &[ReferenceRun::new(0, 13, Motif::new(b"A").unwrap()).unwrap()]
                );
            }
            other => panic!("expected a bundle locus, got {other:?}"),
        }
        assert!(
            bundles
                .next_locus(&members(), &reads)
                .expect("no error")
                .is_none(),
            "the second poll ends the bundle"
        );
    }

    /// `chr1` with a real bundle at 41..=72 (1-based): `(CA)6 TT (GAT)6`, between flanks that
    /// continue neither run. The other fixture contigs stay all `A`.
    const DESIGNED_LEFT: &[u8] = b"GTCCGATTGACCTAGGCATCGTTACGGATCAGTCCATGTG";
    const DESIGNED_HULL: &[u8] = b"CACACACACACATTGATGATGATGATGATGAT";
    const DESIGNED_RIGHT: &[u8] = b"CCTAGGACTTCGAGCTAGTCCAGTTAGC";

    fn designed_chr1() -> Vec<u8> {
        [DESIGNED_LEFT, DESIGNED_HULL, DESIGNED_RIGHT].concat()
    }

    fn designed_ref_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_named_contigs(
            matching_contigs()
                .iter()
                .map(|(name, len, _)| {
                    let bases = if *name == "chr1" {
                        designed_chr1()
                    } else {
                        vec![b'A'; *len]
                    };
                    ((*name).to_string(), bases)
                })
                .collect(),
        )
    }

    /// A read copying `chr1` from 1-based `start` for `length` bases, all aligned.
    fn designed_read(name: &str, start: usize, length: usize) -> noodles_sam::alignment::RecordBuf {
        let mut read = read_named_with_length(name, 0, start, length);
        *read.sequence_mut() = noodles_sam::alignment::record_buf::Sequence::from(
            designed_chr1()[start - 1..start - 1 + length].to_vec(),
        );
        read
    }

    /// **Reads that end inside a bundle's hull are dropped and counted by that reason; reads
    /// that cross it are complete observations of the whole hull; none is a partial.** Five reads
    /// over the designed bundle: two cross it (one copying the whole fetched window, one the
    /// delimiter's to measure), one ends inside it, one starts inside it, one lies in the left
    /// flank alone.
    #[test]
    fn a_bundle_drops_and_counts_the_reads_ending_inside_and_keeps_the_crossing_ones() {
        assert_eq!(
            designed_chr1().len(),
            100,
            "the fixture's chr1 is 100 bases"
        );
        let records = vec![
            designed_read("in_the_flank", 1, 38),
            designed_read("ends_inside", 1, 60),
            designed_read("crosses", 21, 60),
            designed_read("crosses_the_window", 25, 60),
            designed_read("starts_inside", 50, 51),
        ];
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut bundles = SsrBundleGenerator::new(
            designed_ref_bases(),
            designed_ref_bases as fn() -> InMemoryRefSeq,
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
            MinCopies::default(),
        )
        .expect("flank within the bundle threshold");
        bundles.begin_segment(hull(41, 72));
        let members = [
            BundleMember::new(41, 52, Motif::new(b"CA").unwrap()),
            BundleMember::new(55, 72, Motif::new(b"GAT").unwrap()),
        ];
        let locus = bundles
            .next_locus(&members, &reads)
            .expect("no error")
            .expect("one locus");

        let counts = bundles.counts();
        assert_eq!(counts.reads_fetched, 5);
        assert_eq!(counts.ended_inside_a_bundle, 2, "{counts:?}");
        assert_eq!(counts.outside_tract, 1, "{counts:?}");
        assert_eq!(counts.observations_partial, 0);
        assert_eq!(locus.reads_without_observation, 3);
        assert_eq!(locus.observations.len(), 1, "{:?}", locus.observations);
        let crossing = &locus.observations[0];
        assert_eq!(crossing.read_witness, ReadWitness::Complete);
        assert_eq!(&*crossing.bases, DESIGNED_HULL);
        assert_eq!(crossing.num_obs, 2);
        match &locus.kind {
            LocusKind::SsrBundle(detail) => assert_eq!(
                &*detail.runs,
                &[
                    ReferenceRun::new(0, 12, Motif::new(b"CA").unwrap()).unwrap(),
                    ReferenceRun::new(13, 19, Motif::new(b"TGA").unwrap()).unwrap(),
                ]
            ),
            other => panic!("expected a bundle locus, got {other:?}"),
        }
    }

    #[test]
    fn the_edges_take_the_period_of_the_run_reaching_them() {
        let a = Motif::new(b"A").unwrap();
        let cag = Motif::new(b"CAG").unwrap();
        let runs = [
            ReferenceRun::new(0, 9, a).unwrap(),
            ReferenceRun::new(12, 18, cag).unwrap(),
        ];
        assert_eq!(
            edge_run_periods(&runs, 30),
            EdgeRunPeriods {
                left: Some(1),
                right: Some(3)
            }
        );
        assert_eq!(
            edge_run_periods(&runs, 31),
            EdgeRunPeriods {
                left: Some(1),
                right: None
            },
            "a fixed base after the last run"
        );
        assert_eq!(
            edge_run_periods(&runs[1..], 30),
            EdgeRunPeriods {
                left: None,
                right: Some(3)
            }
        );
        assert_eq!(
            edge_run_periods(&[], 30),
            EdgeRunPeriods {
                left: None,
                right: None
            }
        );
    }
}
