//! **The other end of the differential: production's selector, run against ng's on the same
//! tomato loci** (spec §10; arch *Test & bench shape*).
//!
//! Three of production's rules are replaced on purpose, so byte-identical output is impossible
//! by construction and a parity test with an escape clause would have no failing state. What
//! this file does instead is switch the three back **in** — a test-only re-implementation of
//! production's clear-peak nomination, its cohort-summed depth gate and its same-length sibling
//! bar, driving ng's own fold and ng's own ladder — and require the result to be production's
//! candidate set on real reads.
//!
//! **The rules are re-implemented here and not made a field of [`SsrSelectionConfig`]**: the
//! shipping binary carries one rule, and a configuration nobody should ever set does not belong
//! in it.
//!
//! # What the fixture is, and what it is not
//!
//! `testdata/tomato_tract_*.csv` are the merge's own tables and per-sample read counts at **269
//! repeat tracts of the 51-accession tomato panel**, taken from the first 400 intervals of
//! `benchmarks/ssr_tomato1/ssr_regions.bed` through `examples/ng_candidate_selection_probe`
//! (`NG_TRACT_DUMP` and `NG_TRACT_ROWS`) on 2026-09-02. Real reads, real alignment, a real
//! cohort — the panel is about three reads a position, which is the end of the depth range where
//! production's absolute-count rules actually bind.
//!
//! **Two things the fixture drops, both harmless to this comparison and neither harmless in
//! general.** Read groups are pooled into one count a sample, which is the shape production's
//! own Stage-1 evidence has; and reads that stopped inside the locus — the merge's partial
//! observations — are not carried, so a sample's *compared* read count here is the sum of its
//! rows. Production's three rules read neither: its nomination counts reads at a length, its
//! depth gate sums the cohort, and its sibling bar divides by a rung's own total.
//!
//! # The comparison is on the set, not on the order
//!
//! Production walks its nominated lengths in support order and ng's rescue returns its rungs
//! ascending, so the two push the same sequences in different orders. The order a locus's
//! alternatives come out in is the merge table's own on the shipped path and is pinned by the
//! module's unit tests; what this file asserts is that **the same sequences survive**.

use std::collections::{BTreeMap, BTreeSet};

use super::fixtures::{row, sample_showing};
use super::ssr::{SsrSelectionConfig, build_ladder, rescue_occupied_neighbours, select_ssr};
use super::{SelectionScratch, summarise_alleles};
use crate::ng::locus_generation::{LocusKind, SsrDetail};
use crate::ng::run::cohort_merge::build::CohortObservation;
use crate::ng::types::{ContigId, GenomeRegion, Motif, Ploidy, Position};

use crate::ssr::cohort::candidate_set::{CandidateCfg, assemble_candidates};
use crate::ssr::cohort::rung_ladder::{RungCfg, build_rungs};
use crate::ssr::cohort::types::{CohortLocus, LocusId, SampleEvidence, SsrQc};

/// One tract of the fixture: the merge's allele table and every covering sample's reads on it.
struct FixtureTract {
    contig: u32,
    start: u32,
    end: u32,
    motif: Vec<u8>,
    /// The merge table, allele 0 the reference tract.
    alleles: Vec<Vec<u8>>,
    /// `(sample, allele, reads)`, ascending by sample then allele.
    rows: Vec<(usize, usize, u32)>,
}

/// Read the two committed CSVs into one tract per locus, in ascending locus order.
fn fixture_tracts() -> Vec<FixtureTract> {
    let alleles_csv = include_str!("testdata/tomato_tract_alleles.csv");
    let rows_csv = include_str!("testdata/tomato_tract_rows.csv");

    let mut by_locus: BTreeMap<(u32, u32, u32), FixtureTract> = BTreeMap::new();
    for line in alleles_csv.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        assert_eq!(f.len(), 6, "the allele table has six columns: {line}");
        let key = (
            f[0].parse().expect("a contig id"),
            f[1].parse().expect("a start"),
            f[2].parse().expect("an end"),
        );
        let tract = by_locus.entry(key).or_insert_with(|| FixtureTract {
            contig: key.0,
            start: key.1,
            end: key.2,
            motif: f[3].as_bytes().to_vec(),
            alleles: Vec::new(),
            rows: Vec::new(),
        });
        let at: usize = f[4].parse().expect("an allele index");
        assert_eq!(
            at,
            tract.alleles.len(),
            "the allele table is written in index order",
        );
        tract.alleles.push(f[5].as_bytes().to_vec());
    }
    for line in rows_csv.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        assert_eq!(f.len(), 6, "an evidence row has six columns: {line}");
        let key = (
            f[0].parse().expect("a contig id"),
            f[1].parse().expect("a start"),
            f[2].parse().expect("an end"),
        );
        let tract = by_locus
            .get_mut(&key)
            .expect("every evidence row names a tract the allele table holds");
        tract.rows.push((
            f[3].parse().expect("a sample index"),
            f[4].parse().expect("an allele index"),
            f[5].parse().expect("a read count"),
        ));
    }
    by_locus.into_values().collect()
}

impl FixtureTract {
    /// The locus as ng's calling path sees it.
    fn as_ng_observation(&self) -> CohortObservation {
        let mut by_sample: BTreeMap<usize, Vec<(usize, u32)>> = BTreeMap::new();
        for &(sample, allele, reads) in &self.rows {
            by_sample.entry(sample).or_default().push((allele, reads));
        }
        let per_sample = by_sample
            .into_iter()
            .map(|(sample, mut rows)| {
                rows.sort_unstable();
                sample_showing(
                    sample,
                    rows.into_iter()
                        // The error mass is the fold's third ranking key and no rule under test
                        // reads it, so one nat a read stands in for what the fixture does not
                        // carry.
                        .map(|(allele, reads)| row(allele, reads, -f64::from(reads)))
                        .collect(),
                )
            })
            .collect();
        let mut observation = super::fixtures::locus_of(
            &self
                .alleles
                .iter()
                .map(Vec::as_slice)
                .collect::<Vec<&[u8]>>(),
            per_sample,
        );
        observation.region = GenomeRegion {
            contig: ContigId(self.contig),
            start: Position(u64::from(self.start)),
            end: Position(u64::from(self.end)),
        };
        observation.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(&self.motif).expect("the fixture's motif"),
            left_flank: Box::from(&b""[..]),
            right_flank: Box::from(&b""[..]),
        });
        observation
    }

    /// The same tract with only its lowest-numbered covering accession — a cohort of one, which
    /// is the other end of the size range every rule here has to answer at.
    fn first_sample_alone(&self) -> Self {
        let first = self
            .rows
            .iter()
            .map(|&(sample, _, _)| sample)
            .min()
            .expect("a tract the merge built has a covering sample");
        Self {
            contig: self.contig,
            start: self.start,
            end: self.end,
            motif: self.motif.clone(),
            alleles: self.alleles.clone(),
            rows: self
                .rows
                .iter()
                .copied()
                .filter(|&(sample, _, _)| sample == first)
                .collect(),
        }
    }

    /// The same locus as production's Stage-2 sees it.
    fn as_production_locus(&self) -> CohortLocus {
        let mut locus = CohortLocus::new(
            LocusId {
                chrom_id: self.contig,
                start: self.start,
                end: self.end,
            },
            crate::ssr::types::Motif::new(&self.motif).expect("the fixture's motif"),
            self.alleles[0].clone().into_boxed_slice(),
            self.alleles[0].clone().into_boxed_slice(),
        );
        let mut by_sample: BTreeMap<usize, Vec<(usize, u32)>> = BTreeMap::new();
        for &(sample, allele, reads) in &self.rows {
            by_sample.entry(sample).or_default().push((allele, reads));
        }
        for (sample, rows) in by_sample {
            // Production's contract is one entry per distinct sequence, byte-sorted.
            let mut seq_counts: Vec<(Box<[u8]>, u32)> = rows
                .into_iter()
                .map(|(allele, reads)| (self.alleles[allele].clone().into_boxed_slice(), reads))
                .collect();
            seq_counts.sort_by(|(left, _), (right, _)| left.cmp(right));
            locus.push(
                u32::try_from(sample).expect("a sample index"),
                SampleEvidence {
                    seq_counts,
                    qc: SsrQc::default(),
                },
            );
        }
        locus
    }
}

/// **Production's three replaced rules, over ng's ladder** — the candidate set production
/// assembles, expressed against the merge's table instead of against its own evidence rows.
///
/// The constants are production's own defaults, retyped here rather than read from
/// `CandidateCfg` so that a change to production's development defaults turns this red instead
/// of silently moving both sides of the comparison together.
mod productions_rules {
    use super::*;

    const PROMINENCE: u32 = 3;
    const MIN_COHORT_DEPTH: u64 = 10;
    const MAX_OUT_OF_FRAME_FRAC: f64 = 0.10;
    const MIN_SAME_LENGTH_READS: u32 = 8;
    const MIN_SAME_LENGTH_SAMPLES: u32 = 3;
    const MIN_SAME_LENGTH_FRACTION: f64 = 0.10;

    /// One sample's reads at each allele of the merge table, and the table's own read totals.
    struct Pooled {
        /// `per_sample[k][allele]` — sample `k`'s reads on that allele.
        per_sample: Vec<Vec<u32>>,
        /// Cohort reads on each allele.
        cohort: Vec<u32>,
        /// How many distinct samples showed each allele — production's recurrence term.
        samples: Vec<u32>,
    }

    fn pool(observation: &CohortObservation) -> Pooled {
        let width = observation.alleles.len();
        let mut per_sample = Vec::with_capacity(observation.per_sample.len());
        let mut cohort = vec![0_u32; width];
        let mut samples = vec![0_u32; width];
        for sample in &observation.per_sample {
            let mut here = vec![0_u32; width];
            for row in &sample.supported {
                here[row.allele] += row.support.num_reads;
            }
            for (at, reads) in here.iter().enumerate() {
                if *reads > 0 {
                    cohort[at] += reads;
                    samples[at] += 1;
                }
            }
            per_sample.push(here);
        }
        Pooled {
            per_sample,
            cohort,
            samples,
        }
    }

    /// Production's `is_periodic`: the **cohort's** off-grid read share against a grid anchored
    /// on the commonest observed length, in bases (`candidate_set.rs:114-145`).
    fn is_periodic(observation: &CohortObservation, pooled: &Pooled, period: usize) -> bool {
        if period <= 1 {
            return true;
        }
        let mut support_by_length: BTreeMap<usize, u64> = BTreeMap::new();
        for (at, reads) in pooled.cohort.iter().enumerate() {
            if *reads > 0 {
                *support_by_length
                    .entry(observation.alleles[at].len())
                    .or_insert(0) += u64::from(*reads);
            }
        }
        let total: u64 = support_by_length.values().sum();
        if total == 0 {
            return true;
        }
        let mut anchor = 0_usize;
        let mut best = 0_u64;
        for (&length, &support) in &support_by_length {
            if support > best {
                best = support;
                anchor = length;
            }
        }
        let off_grid: u64 = support_by_length
            .iter()
            .filter(|&(&length, _)| length.abs_diff(anchor) % period != 0)
            .map(|(_, &support)| support)
            .sum();
        off_grid as f64 <= MAX_OUT_OF_FRAME_FRAC * total as f64
    }

    /// Production's `is_clear_peak`: a length whose reads exceed both neighbours' by more than
    /// the prominence (`rung_ladder.rs:274-288`).
    fn clear_peaks(reads_per_rung: &[u32], repeat_counts: &[u32]) -> Vec<usize> {
        let at_count: BTreeMap<u32, u32> = repeat_counts
            .iter()
            .copied()
            .zip(reads_per_rung.iter().copied())
            .collect();
        (0..reads_per_rung.len())
            .filter(|&rung| {
                let support = reads_per_rung[rung];
                if support == 0 {
                    return false;
                }
                let count = repeat_counts[rung];
                let lower = count
                    .checked_sub(1)
                    .and_then(|c| at_count.get(&c).copied())
                    .unwrap_or(0);
                let upper = at_count.get(&(count + 1)).copied().unwrap_or(0);
                support > lower + PROMINENCE && support > upper + PROMINENCE
            })
            .collect()
    }

    /// The candidate sequences production assembles at this tract, the reference first.
    pub(super) fn candidate_set(
        observation: &CohortObservation,
        ploidy: Ploidy,
        config: &SsrSelectionConfig,
        scratch: &mut SelectionScratch,
    ) -> Vec<Vec<u8>> {
        let LocusKind::Ssr(detail) = &observation.kind else {
            panic!("the differential runs on repeat tracts only");
        };
        let period = detail.motif.period();
        let pooled = pool(observation);
        let reference = vec![observation.alleles[0].to_vec()];

        let depth: u64 = pooled.cohort.iter().map(|&reads| u64::from(reads)).sum();
        if depth < MIN_COHORT_DEPTH {
            return reference;
        }
        if !is_periodic(observation, &pooled, period) {
            return reference;
        }

        // The same fold and the same ladder the shipped path runs on.
        summarise_alleles(observation, config.shared.min_allele_support, scratch);
        build_ladder(observation, &detail.motif, scratch);
        let rungs = scratch.ladder.rung_count();
        let repeat_counts: Vec<u32> = (0..rungs)
            .map(|rung| scratch.ladder.repeat_count_at(rung))
            .collect();

        let mut admitted: Vec<Vec<u8>> = reference;
        for here_reads in &pooled.per_sample {
            let mut reads_per_rung = vec![0_u32; rungs];
            for (rung, total) in reads_per_rung.iter_mut().enumerate() {
                for &at in scratch.ladder.table_indices_at(rung) {
                    *total += here_reads[at as usize];
                }
            }
            let mut peaks = clear_peaks(&reads_per_rung, &repeat_counts);
            peaks.sort_by(|&left, &right| {
                reads_per_rung[right]
                    .cmp(&reads_per_rung[left])
                    .then_with(|| repeat_counts[left].cmp(&repeat_counts[right]))
            });
            let copies = usize::from(ploidy.get());
            let mut nominated: Vec<u32> = peaks
                .iter()
                .take(copies)
                .map(|&rung| u32::try_from(rung).expect("a rung index"))
                .collect();
            nominated.sort_unstable();
            // The one part of production's nomination ng keeps, so it is called rather than
            // written twice (spec §10's reuse map).
            rescue_occupied_neighbours(&scratch.ladder, ploidy, &mut nominated);

            for &rung in &nominated {
                let rung = rung as usize;
                let total_at_rung: u32 = scratch
                    .ladder
                    .table_indices_at(rung)
                    .iter()
                    .map(|&at| pooled.cohort[at as usize])
                    .sum();
                // Byte-sorted within the rung, as production's `seqs_at` yields.
                let mut here: Vec<usize> = scratch
                    .ladder
                    .table_indices_at(rung)
                    .iter()
                    .map(|&at| at as usize)
                    .filter(|&at| pooled.cohort[at] > 0)
                    .collect();
                here.sort_by(|&left, &right| {
                    observation.alleles[left].cmp(&observation.alleles[right])
                });
                // The representative: most reads, ties to the lexicographically smallest.
                let representative = here.iter().copied().max_by(|&left, &right| {
                    pooled.cohort[left]
                        .cmp(&pooled.cohort[right])
                        .then_with(|| observation.alleles[right].cmp(&observation.alleles[left]))
                });
                for at in here {
                    let promoted = Some(at) == representative
                        || (pooled.cohort[at] >= MIN_SAME_LENGTH_READS
                            && pooled.samples[at] >= MIN_SAME_LENGTH_SAMPLES
                            && f64::from(pooled.cohort[at])
                                >= MIN_SAME_LENGTH_FRACTION * f64::from(total_at_rung));
                    if promoted
                        && !admitted
                            .iter()
                            .any(|seq| seq == &observation.alleles[at][..])
                    {
                        admitted.push(observation.alleles[at].to_vec());
                    }
                }
            }
        }
        admitted
    }
}

/// Production's own selector's answer at this tract, through its own code.
fn productions_own_answer(tract: &FixtureTract) -> BTreeSet<Vec<u8>> {
    let locus = tract.as_production_locus();
    let rungs = build_rungs(&locus, &RungCfg::dev_default());
    let set = assemble_candidates(&locus, &rungs, 2, &CandidateCfg::dev_default());
    set.alleles.iter().map(|seq| seq.to_vec()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("two copies")
    }

    /// **The differential's production end: with the three replaced rules switched back in, the
    /// selector reproduces production's candidate set at every tract of the tomato fixture.**
    ///
    /// This is the assertion spec §10 asks for. It has a failing state at both ends: a rule
    /// re-implemented wrongly here fails against production's own code, and a change to ng's
    /// ladder or fold — which this arm drives — fails the same way.
    #[test]
    fn productions_rules_switched_in_reproduce_productions_candidate_set_on_tomato() {
        let tracts = fixture_tracts();
        assert!(
            tracts.len() > 200,
            "the fixture is the tomato panel's 269 tracts, not a handful: {}",
            tracts.len()
        );
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        let mut compared = 0_usize;
        for tract in &tracts {
            let observation = tract.as_ng_observation();
            let ours: BTreeSet<Vec<u8>> =
                productions_rules::candidate_set(&observation, diploid(), &config, &mut scratch)
                    .into_iter()
                    .collect();
            let theirs = productions_own_answer(tract);
            assert_eq!(
                ours, theirs,
                "the tract at contig {} {}-{} narrows differently under production's rules",
                tract.contig, tract.start, tract.end,
            );
            compared += 1;
        }
        assert_eq!(compared, tracts.len());
    }

    /// **The same tracts with one accession instead of fifty-one** — which is where two of
    /// production's rules actually bite, and the 51-sample comparison above cannot see either.
    ///
    /// Its cohort-summed depth gate needs ten reads *across the panel*, and fifty-one accessions
    /// at three reads a position clear that at every tract; alone, one accession does not, so
    /// production refuses most of these loci outright and offers the reference tract alone
    /// (spec §6). Its three-distinct-samples term on a same-length sibling cannot be reached at
    /// all below three samples (spec §5). **Both rules are inert in the fixture above**: raising
    /// the depth gate's constant or dropping the sample term to one leaves that test green.
    #[test]
    fn productions_rules_are_reproduced_with_one_accession_too_where_two_of_them_first_bite() {
        let tracts = fixture_tracts();
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        let mut refused_for_depth = 0_usize;
        for tract in &tracts {
            let alone = tract.first_sample_alone();
            let observation = alone.as_ng_observation();
            let ours: BTreeSet<Vec<u8>> =
                productions_rules::candidate_set(&observation, diploid(), &config, &mut scratch)
                    .into_iter()
                    .collect();
            let theirs = productions_own_answer(&alone);
            assert_eq!(
                ours, theirs,
                "the tract at contig {} {}-{} narrows differently under production's rules with \
                 one accession",
                tract.contig, tract.start, tract.end,
            );
            if theirs.len() == 1 && observation.alleles.len() > 1 {
                refused_for_depth += 1;
            }
        }
        // Pinned, because the point of this test is that the gate *fires* here: a fixture recut
        // deeper would leave it agreeing for the wrong reason. Measured 2026-09-02, and it is
        // every tract — one tomato accession never reaches ten reads at a tract.
        assert_eq!(
            refused_for_depth,
            tracts.len(),
            "one accession at three reads a position leaves production offering the reference \
             tract alone at every one of the fixture's {} tracts, where the whole panel refuses \
             none of them",
            tracts.len()
        );
    }

    /// **A same-length sibling two accessions carry and a third does not** — production's
    /// three-distinct-samples recurrence term, which neither run above can reach.
    ///
    /// At fifty-one accessions of three reads, eight cohort reads on a sequence already implies
    /// three accessions showed it, so the sample term never decides anything; at one accession
    /// the depth gate refuses the locus before the sibling bar is asked. Hand-built, therefore:
    /// the sibling carries 20 reads from two accessions — clear of the eight-read floor and of
    /// a tenth of its rung's 120 — and production refuses it for the third accession alone.
    #[test]
    fn productions_recurrence_term_is_reproduced_on_a_sibling_two_accessions_carry() {
        let tract = FixtureTract {
            contig: 0,
            start: 100,
            end: 109,
            motif: b"AT".to_vec(),
            alleles: vec![b"ATATATATAT".to_vec(), b"ATATATATAG".to_vec()],
            rows: vec![(0, 0, 50), (1, 0, 50), (2, 1, 10), (3, 1, 10)],
        };
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        let observation = tract.as_ng_observation();
        let ours: BTreeSet<Vec<u8>> =
            productions_rules::candidate_set(&observation, diploid(), &config, &mut scratch)
                .into_iter()
                .collect();
        assert_eq!(
            ours,
            productions_own_answer(&tract),
            "production refuses the sibling for recurrence and this arm has to refuse it too",
        );
        assert_eq!(
            ours,
            BTreeSet::from([b"ATATATATAT".to_vec()]),
            "and what survives is the rung's representative alone",
        );
    }

    /// **A tract whose reads do not sit on the motif grid**, hand-built because the tomato
    /// fixture holds none that production's own measure refuses — widening its off-grid share
    /// to one leaves both tests above green.
    ///
    /// Production anchors the grid on the commonest observed length in bases and refuses the
    /// locus when more than a tenth of the **cohort's** reads sit off it
    /// (`candidate_set.rs:114-145`). Eleven-base and nine-base reads against a ten-base
    /// dinucleotide mode are one base off the grid, and here they are 12 reads of 22.
    #[test]
    fn productions_periodicity_gate_is_reproduced_on_an_off_grid_tract() {
        let tract = FixtureTract {
            contig: 0,
            start: 100,
            end: 109,
            motif: b"AT".to_vec(),
            alleles: vec![
                b"ATATATATAT".to_vec(),
                b"ATATATATATA".to_vec(),
                b"ATATATATA".to_vec(),
            ],
            rows: vec![(0, 0, 10), (0, 1, 6), (0, 2, 6)],
        };
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        let observation = tract.as_ng_observation();
        let ours: BTreeSet<Vec<u8>> =
            productions_rules::candidate_set(&observation, diploid(), &config, &mut scratch)
                .into_iter()
                .collect();
        assert_eq!(
            ours,
            productions_own_answer(&tract),
            "production refuses this tract for periodicity and offers the reference alone",
        );
        assert_eq!(
            ours,
            BTreeSet::from([b"ATATATATAT".to_vec()]),
            "and what it offers is the reference tract",
        );
    }

    /// **And the differential's other end has to move**, or the comparison above would pass
    /// against a selector that had changed nothing. The shipped rules and production's disagree
    /// at a real share of the fixture's tracts, and this records how large that share is.
    #[test]
    fn the_shipped_rules_narrow_differently_from_productions_on_the_same_tracts() {
        let tracts = fixture_tracts();
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        let mut moved = 0_usize;
        let (mut ours_total, mut theirs_total) = (0_usize, 0_usize);
        for tract in &tracts {
            let observation = tract.as_ng_observation();
            let shipped: BTreeSet<Vec<u8>> = select_ssr(&observation, &config, &mut scratch)
                .selection
                .alleles()
                .iter()
                .map(<[u8]>::to_vec)
                .collect();
            let theirs = productions_own_answer(tract);
            ours_total += shipped.len();
            theirs_total += theirs.len();
            if shipped != theirs {
                moved += 1;
            }
        }
        // Pinned rather than bounded, because these three numbers are the measurement this
        // fixture exists to make and a range would let all three drift unnoticed. Measured
        // 2026-09-02 on the 269-tract fixture at ploidy 2.
        assert_eq!(
            (moved, ours_total, theirs_total),
            (184, 602, 668),
            "the shipped rules narrow 184 of the fixture's 269 tracts differently from \
             production's, to 602 candidate sequences against production's 668 — spec §4.1 and \
             §5 measure the replacement as *cheaper* in candidates as well as better in recall, \
             and on this panel it stays so",
        );
    }
}
