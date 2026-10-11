//! **Gathering the calls: an STR record widened over a flank, and the generic records it
//! overlaps** (`doc/devel/ng/spec/str_edge_deletion_alleles.md` §6).
//!
//! A deletion running from an STR segment into a flank is one allele of the segment's record,
//! and that record is written widened over the flank bases its alleles reach
//! ([`assemble_record`](crate::vcf::assemble_record)). The generic path still calls the flank:
//! it wrote what the same reads showed there as an allele of its own record — the flank half of
//! the deletion, or, where the mapper spelled the event another way, a substitution — and it may
//! have written other variants there, another sample's SNP, say. This pass, over the finished
//! records in genome order and before anything prepares them for the filters and the VCF, makes
//! the two agree, sample by sample, because the STR record's alleles are whole events phased by
//! the reads that crossed the segment:
//!
//! 1. **every copy of a sample whose STR allele reaches a base of a generic record becomes `*`
//!    there** — "deleted by an overlapping record", as VCF writes it. Which of the generic
//!    genotype's copies it takes is the one the STR allele's own sequence spells over that record
//!    — the flank half, whatever form the generic path gave it — and only where none does, a
//!    reference copy, then any;
//! 2. **where the generic call did not say so, the STR record wins** (owner, 2026-10-10): the copy
//!    is rewritten anyway, and the sample's `GQ` there becomes the lower of the two records';
//! 3. **a generic record whose genotypes name nothing but the reference and `*` is dropped**, and
//!    an allele no genotype names any more is dropped with its reads moved to the unexplained
//!    count, as assembly does.
//!
//! **Records near an STR segment are held**, because a record widened to the left starts before
//! generic records that arrived ahead of it. Anything ending more than [`HOLD_REACH`] bases before
//! the newest record's start can no longer meet a widened record and is released, in the
//! writer's order. Records after a segment are reconciled as they arrive, against the widened
//! records whose right reach they fall in.

use crate::types::{AlleleId, GenomeRegion, Genotype, Phred};
use crate::vcf::{MapqPool, SampleCall, SampleColumn, SampleReadCounts, TractWidening, VcfRecord};
use crate::window_coverage::WindowCoverage;

/// **How far before the newest record a held record may still be met by a widened one.** A
/// record is widened by at most the flank its locus holds, and a flank deletion is at most 255
/// bases (`FlankDeletion` counts in a byte); one more for the padding base.
const HOLD_REACH: u64 = 256;

/// The `*` allele: a base this record covers that an overlapping record's allele deletes.
pub(crate) const SPANNING_DELETION: &[u8] = b"*";

/// What the pass changed, for the run's report.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeDeletionGatheringCounts {
    /// Copies of a sample's generic genotype that became `*` on the allele the STR allele spells
    /// there — the two paths agreeing.
    pub copies_starred_in_agreement: u64,
    /// Copies that became `*` against what the generic call said — the STR record winning.
    pub copies_rewritten: u64,
    /// Generic records dropped because no genotype named anything but `*` and the reference.
    pub generic_records_dropped: u64,
}

/// One widened STR record, remembered for the generic records it overlaps.
#[derive(Debug)]
struct WidenedTract {
    region: GenomeRegion,
    /// The tract's first and last base, 1-based.
    tract_first: u64,
    tract_last: u64,
    /// The record's alleles over its widened reference, the reference first.
    alleles: Vec<Box<[u8]>>,
    reach_per_allele: Vec<(u8, u8)>,
    /// Each sample's called alleles, `None` at a no-call.
    genotypes: Vec<Option<Vec<AlleleId>>>,
    genotype_qualities: Vec<Option<Phred>>,
}

impl WidenedTract {
    fn of(record: &VcfRecord, widening: &TractWidening) -> Self {
        let region = record.region();
        Self {
            region,
            tract_first: region.start.get() + u64::from(widening.left),
            tract_last: region.end.get() - u64::from(widening.right),
            alleles: record.alleles().to_vec(),
            reach_per_allele: widening.reach_per_allele.clone(),
            genotypes: record
                .sample_columns()
                .iter()
                .map(|column| column.call.genotype().map(|g| g.alleles().to_vec()))
                .collect(),
            genotype_qualities: record
                .sample_columns()
                .iter()
                .map(|column| column.call.genotype_quality())
                .collect(),
        }
    }

    /// Whether allele `allele` changes a base of `region` where the generic path looks: its
    /// reach into the left flank or the right one overlaps it.
    fn reaches(&self, allele: AlleleId, region: GenomeRegion) -> bool {
        let (left, right) = self
            .reach_per_allele
            .get(usize::from(allele.get()))
            .copied()
            .unwrap_or((0, 0));
        let (start, end) = (region.start.get(), region.end.get());
        let on_the_left =
            left > 0 && start < self.tract_first && end >= self.tract_first - u64::from(left);
        let on_the_right =
            right > 0 && end > self.tract_last && start <= self.tract_last + u64::from(right);
        on_the_left || on_the_right
    }

    /// **Does STR allele `allele` spell `generic` over `region`?** — whether a generic allele is
    /// what the reads of this STR allele showed the generic path, in whatever form the generic
    /// path wrote it. The generic record sits in one flank; the STR allele's sequence is read
    /// from that flank's end of the widened record inward, and at the generic record's own
    /// positions it must spell the generic allele. The bases between are not required to be the
    /// reference's: an STR allele can change two neighbouring flank bases the generic path called
    /// as two records (`AC`→`CT` on the stuttering benchmark's `seg00095`).
    /// `reference` is the generic record's reference, for the bases it covers outside the
    /// widened record.
    fn spells(
        &self,
        allele: AlleleId,
        region: GenomeRegion,
        reference: &[u8],
        generic: &[u8],
    ) -> bool {
        let ours = &*self.alleles[usize::from(allele.get())];
        let (start, end) = (region.start.get(), region.end.get());
        let (window_start, window_end) = (self.region.start.get(), self.region.end.get());
        if end < self.tract_first {
            // Left flank: bases of the generic record before the widened one must be unchanged.
            let outside = window_start.saturating_sub(start) as usize;
            if outside > reference.len() || !generic.starts_with(&reference[..outside]) {
                return false;
            }
            let generic = &generic[outside..];
            let between = (start.max(window_start) - window_start) as usize;
            between <= ours.len() && ours[between..].starts_with(generic)
        } else if start > self.tract_last {
            // Right flank: the mirror image.
            let outside = end.saturating_sub(window_end) as usize;
            if outside > reference.len()
                || !generic.ends_with(&reference[reference.len() - outside..])
            {
                return false;
            }
            let generic = &generic[..generic.len() - outside];
            let between = (window_end - end.min(window_end)) as usize;
            between <= ours.len() && ours[..ours.len() - between].ends_with(generic)
        } else {
            false
        }
    }

    fn overlaps(&self, region: GenomeRegion) -> bool {
        region.contig == self.region.contig
            && region.start <= self.region.end
            && region.end >= self.region.start
    }
}

/// The pass itself: feed it every finished record in the order the drivers emit them, and
/// forward what it releases.
#[derive(Debug, Default)]
pub struct EdgeDeletionGathering {
    held: Vec<(VcfRecord, Vec<WindowCoverage>)>,
    widened: Vec<WidenedTract>,
    counts: EdgeDeletionGatheringCounts,
}

impl EdgeDeletionGathering {
    /// Take one finished record, and hand back the records now safe to write, in the writer's
    /// order.
    pub fn accept(
        &mut self,
        record: VcfRecord,
        windows: Vec<WindowCoverage>,
    ) -> Vec<(VcfRecord, Vec<WindowCoverage>)> {
        let region = record.region();
        let mut released = Vec::new();
        // A new contig: everything held is behind it.
        if self
            .held
            .first()
            .is_some_and(|(held, _)| held.region().contig != region.contig)
        {
            released = self.finish_held();
            self.widened.clear();
        }
        // Widened records the stream has passed for good.
        self.widened.retain(|tract| {
            tract.region.contig == region.contig && tract.region.end >= region.start
        });

        if let Some(widening) = record.repeat_tract().and_then(|tract| tract.widening()) {
            let tract = WidenedTract::of(&record, widening);
            // The generic records already held that it overlaps — its left flank.
            let mut index = 0;
            while index < self.held.len() {
                if !self.held[index].0.is_repeat_tract()
                    && tract.overlaps(self.held[index].0.region())
                {
                    let (held, held_windows) = self.held.remove(index);
                    match reconcile(held, &tract, &mut self.counts) {
                        Some(kept) => {
                            self.held.insert(index, (kept, held_windows));
                            index += 1;
                        }
                        None => {}
                    }
                } else {
                    index += 1;
                }
            }
            self.widened.push(tract);
            self.held.push((record, windows));
        } else if record.is_repeat_tract() {
            self.held.push((record, windows));
        } else {
            // A generic record: reconciled against every widened record it falls in.
            let mut record = Some(record);
            for tract in &self.widened {
                if let Some(current) = record.take() {
                    record = if tract.overlaps(current.region()) {
                        reconcile(current, tract, &mut self.counts)
                    } else {
                        Some(current)
                    };
                }
            }
            if let Some(record) = record {
                self.held.push((record, windows));
            }
        }

        // Release what no later record can reach: ending more than the hold reach before the
        // newest record's start, in the writer's order.
        self.held.sort_by_key(|(held, _)| order_key(held));
        let horizon = region.start.get().saturating_sub(HOLD_REACH);
        let releasable = self
            .held
            .iter()
            .take_while(|(held, _)| held.region().end.get() < horizon)
            .count();
        released.extend(self.held.drain(..releasable));
        released
    }

    /// Release everything still held, at the end of the stream.
    pub fn finish(&mut self) -> Vec<(VcfRecord, Vec<WindowCoverage>)> {
        self.widened.clear();
        self.finish_held()
    }

    fn finish_held(&mut self) -> Vec<(VcfRecord, Vec<WindowCoverage>)> {
        self.held.sort_by_key(|(held, _)| order_key(held));
        std::mem::take(&mut self.held)
    }

    /// What the pass changed so far.
    pub fn counts(&self) -> EdgeDeletionGatheringCounts {
        self.counts
    }
}

/// The writer's order: written position, a generic record before a tract at a tie
/// (`vcf_output.md` §5).
fn order_key(record: &VcfRecord) -> (u64, bool) {
    (
        crate::vcf::written_position(record),
        record.is_repeat_tract(),
    )
}

/// **Make one generic record agree with one widened STR record it overlaps**, or `None` when
/// nothing is left of it to write.
fn reconcile(
    record: VcfRecord,
    tract: &WidenedTract,
    counts: &mut EdgeDeletionGatheringCounts,
) -> Option<VcfRecord> {
    let region = record.region();
    let mut alleles: Vec<Box<[u8]>> = record.alleles().to_vec();
    let mut expected_copies = record.expected_copies().to_vec();
    let mut allele_mapq = record.allele_mapq().to_vec();
    let mut columns: Vec<SampleColumn> = record.sample_columns().to_vec();
    let star_index = alleles.len();
    let mut any_star = false;

    for (sample, column) in columns.iter_mut().enumerate() {
        let Some(Some(str_genotype)) = tract.genotypes.get(sample) else {
            continue;
        };
        let SampleCall::Called {
            genotype,
            genotype_quality,
        } = &mut column.call
        else {
            continue;
        };
        let reaching: Vec<AlleleId> = str_genotype
            .iter()
            .copied()
            .filter(|&allele| tract.reaches(allele, region))
            .collect();
        if reaching.is_empty() {
            continue;
        }
        let mut copies: Vec<AlleleId> = genotype.alleles().to_vec();
        let mut taken = vec![false; copies.len()];
        let mut disagreed = false;
        for &str_allele in &reaching {
            // The copy this STR allele explains: the allele its sequence spells here, preferring
            // a variant over the reference; then the reference; then anything.
            let rank = |copy: AlleleId| -> u8 {
                let bases = &alleles[usize::from(copy.get())];
                let spelled = tract.spells(str_allele, region, &alleles[0], bases);
                match (spelled, copy.get() == 0) {
                    (true, false) => 0,
                    (true, true) => 1,
                    (false, true) => 2,
                    (false, false) => 3,
                }
            };
            let Some(chosen) = (0..copies.len())
                .filter(|&at| !taken[at])
                .min_by_key(|&at| rank(copies[at]))
            else {
                break;
            };
            disagreed |= rank(copies[chosen]) >= 2;
            taken[chosen] = true;
            copies[chosen] =
                AlleleId(u16::try_from(star_index).expect("an allele table narrower than 65,536"));
        }
        any_star = true;
        if disagreed {
            counts.copies_rewritten += 1;
            if let Some(Some(str_quality)) = tract.genotype_qualities.get(sample)
                && *str_quality < *genotype_quality
            {
                *genotype_quality = *str_quality;
            }
        } else {
            counts.copies_starred_in_agreement += 1;
        }
        *genotype = Genotype::new(copies);
    }
    if !any_star {
        return Some(record);
    }
    alleles.push(Box::from(SPANNING_DELETION));
    expected_copies.push(0.0);
    allele_mapq.push(MapqPool::default());
    for column in &mut columns {
        let mut reads = column.read_counts.allele_reads().to_vec();
        reads.push(0);
        column.read_counts = SampleReadCounts::new(reads, column.read_counts.unexplained_reads());
    }

    // A record whose genotypes name nothing but the reference and `*` says nothing of its own.
    let names_a_variant_of_its_own = columns.iter().any(|column| {
        column.call.genotype().is_some_and(|genotype| {
            genotype
                .alleles()
                .iter()
                .any(|allele| allele.get() != 0 && usize::from(allele.get()) != star_index)
        })
    });
    if !names_a_variant_of_its_own {
        counts.generic_records_dropped += 1;
        return None;
    }
    let mut no_flank_deletions = Vec::new();
    crate::vcf::drop_alternatives_no_sample_calls(
        &mut alleles,
        &mut expected_copies,
        &mut allele_mapq,
        &mut columns,
        &mut no_flank_deletions,
    );
    let padding_base = record
        .padding_base()
        .filter(|_| alleles.iter().any(|allele| allele.is_empty()));
    Some(VcfRecord::new(
        region,
        alleles,
        expected_copies,
        columns,
        allele_mapq,
        padding_base,
        record.site_quality,
        record.artifact_penalties,
        record.filter(),
        record.repeat_tract().cloned(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{ContigId, Motif, Position};
    use crate::vcf::{FilterVerdict, TractAnnotation};

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn called(alleles: &[u16], quality: f32) -> SampleCall {
        SampleCall::Called {
            genotype: Genotype::new(alleles.iter().copied().map(AlleleId).collect()),
            genotype_quality: Phred::try_new(quality).expect("a quality"),
        }
    }

    /// A record with each sample's call and one read per allele it names.
    fn record(
        at: GenomeRegion,
        alleles: &[&[u8]],
        calls: &[(&[u16], f32)],
        repeat_tract: Option<TractAnnotation>,
    ) -> VcfRecord {
        let columns: Vec<SampleColumn> = calls
            .iter()
            .map(|(genotype, quality)| {
                let mut reads = vec![0_u32; alleles.len()];
                for &allele in *genotype {
                    reads[usize::from(allele)] += 1;
                }
                SampleColumn {
                    call: called(genotype, *quality),
                    read_counts: SampleReadCounts::new(reads, 0),
                }
            })
            .collect();
        let pools = (0..alleles.len())
            .map(|allele| {
                let reads: u64 = columns
                    .iter()
                    .map(|column| u64::from(column.read_counts.allele_reads()[allele]))
                    .sum();
                MapqPool {
                    reads,
                    mapq_sum: reads * 60,
                }
            })
            .collect();
        VcfRecord::new(
            at,
            alleles.iter().map(|bases| Box::from(*bases)).collect(),
            vec![0.5; alleles.len()],
            columns,
            pools,
            None,
            Phred::try_new(50.0).expect("a quality"),
            None,
            FilterVerdict::Pass,
            repeat_tract,
        )
    }

    /// `GTT` ×3 at 101–109, widened one base right over the flank's `A` at 110, with a
    /// deletion of `TTA` taking it (allele 1, reaching one base right) and an expansion
    /// (allele 2).
    fn a_tract_widened_right(calls: &[(&[u16], f32)]) -> VcfRecord {
        record(
            region(101, 110),
            &[b"GTTGTTGTTA", b"GTTGTTG", b"GTTGTTGTTGTTA"],
            calls,
            Some(
                TractAnnotation::new(Motif::new(b"GTT").expect("a motif")).widened(TractWidening {
                    left: 0,
                    right: 1,
                    flank_bases_per_allele: vec![1, 0, 1],
                    reach_per_allele: vec![(0, 0), (0, 1), (0, 0)],
                }),
            ),
        )
    }

    fn written(gathering: &mut EdgeDeletionGathering, records: Vec<VcfRecord>) -> Vec<VcfRecord> {
        let mut out = Vec::new();
        for record in records {
            out.extend(
                gathering
                    .accept(record, Vec::new())
                    .into_iter()
                    .map(|(r, _)| r),
            );
        }
        out.extend(gathering.finish().into_iter().map(|(r, _)| r));
        out
    }

    /// **The flank half the mapper wrote as a substitution is taken by the deletion's copy,
    /// and the record goes** (the stuttering benchmark's `seg02581`, 2026-10-11). The sample
    /// carries the deletion taking the flank's `A`, and the generic path called that `A` as a
    /// `G` from the reads aligned with a mismatch: its `G` copy becomes `*`, and with nothing of
    /// its own left the generic record is dropped.
    #[test]
    fn a_flank_half_written_as_a_substitution_is_taken_and_the_record_dropped() {
        let mut gathering = EdgeDeletionGathering::default();
        let out = written(
            &mut gathering,
            vec![
                a_tract_widened_right(&[(&[1, 2], 30.0)]),
                record(region(110, 110), &[b"A", b"G"], &[(&[0, 1], 40.0)], None),
            ],
        );
        assert_eq!(out.len(), 1, "only the STR record is written");
        assert!(out[0].is_repeat_tract());
        assert_eq!(
            gathering.counts(),
            EdgeDeletionGatheringCounts {
                copies_starred_in_agreement: 1,
                copies_rewritten: 0,
                generic_records_dropped: 1,
            }
        );
    }

    /// **Another sample's own SNP stays, and the STR record wins where the generic call
    /// disagrees.** Sample 0 carries the deletion and its generic call is `0/1`: its `G` copy
    /// becomes `*`. Sample 1 carries the deletion but the generic path called it `0/0` with `GQ`
    /// 50: a reference copy becomes `*`, and `GQ` drops to the STR record's 12. Sample 2 has no
    /// deletion and a real `G`: unchanged, so the record stays, with `G` and `*`.
    #[test]
    fn a_real_snp_stays_and_a_disagreeing_call_is_rewritten() {
        let mut gathering = EdgeDeletionGathering::default();
        let out = written(
            &mut gathering,
            vec![
                a_tract_widened_right(&[(&[1, 2], 30.0), (&[0, 1], 12.0), (&[0, 0], 40.0)]),
                record(
                    region(110, 110),
                    &[b"A", b"G"],
                    &[(&[0, 1], 40.0), (&[0, 0], 50.0), (&[0, 1], 45.0)],
                    None,
                ),
            ],
        );
        assert_eq!(out.len(), 2);
        let generic = &out[1];
        assert_eq!(
            generic.alleles(),
            &[
                Box::from(&b"A"[..]),
                Box::from(&b"G"[..]),
                Box::from(SPANNING_DELETION)
            ]
        );
        let genotypes: Vec<Vec<u16>> = generic
            .sample_columns()
            .iter()
            .map(|column| {
                column
                    .call
                    .genotype()
                    .expect("called")
                    .alleles()
                    .iter()
                    .map(|a| a.get())
                    .collect()
            })
            .collect();
        assert_eq!(genotypes, vec![vec![0, 2], vec![0, 2], vec![0, 1]]);
        assert_eq!(
            generic.sample_columns()[1].call.genotype_quality(),
            Phred::try_new(12.0).ok()
        );
        assert_eq!(gathering.counts().copies_rewritten, 1);
        assert_eq!(gathering.counts().copies_starred_in_agreement, 1);
    }

    /// **A record widened to the left is written before the generic records it now starts
    /// ahead of**, and the flank half of its deletion — a generic deletion of the flank's last
    /// two bases, anchored on the base before — is taken and dropped. A SNP further left, which
    /// it does not reach, is written first, untouched.
    #[test]
    fn a_left_widened_record_is_reordered_and_its_flank_half_taken() {
        // Flank ... C(96) A(97) G(98) C(99) | tract `CA` ×3 at 100–105. The STR record is
        // widened two bases left, to 98, over a deletion of `GC` + `CA` (allele 1).
        let tract = record(
            region(98, 105),
            &[b"GCCACACA", b"CACA"],
            &[(&[0, 1], 30.0)],
            Some(
                TractAnnotation::new(Motif::new(b"CA").expect("a motif")).widened(TractWidening {
                    left: 2,
                    right: 0,
                    flank_bases_per_allele: vec![2, 0],
                    reach_per_allele: vec![(0, 0), (2, 0)],
                }),
            ),
        );
        let mut gathering = EdgeDeletionGathering::default();
        let out = written(
            &mut gathering,
            vec![
                record(region(96, 96), &[b"C", b"T"], &[(&[0, 1], 40.0)], None),
                record(region(97, 99), &[b"AGC", b"A"], &[(&[0, 1], 40.0)], None),
                tract,
            ],
        );
        let starts: Vec<u64> = out.iter().map(|r| r.region().start.get()).collect();
        assert_eq!(
            starts,
            vec![96, 98],
            "the SNP, then the widened tract; the half is gone"
        );
        assert_eq!(gathering.counts().generic_records_dropped, 1);
    }
}
