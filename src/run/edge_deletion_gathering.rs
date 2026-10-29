//! **Gathering the calls: an STR record widened over a flank absorbs the generic records it
//! overlaps** (`doc/devel/ng/spec/str_edge_deletion_alleles.md` §6).
//!
//! A deletion running from an STR segment into a flank is one allele of the segment's record,
//! and that record is written widened over the flank bases its alleles reach
//! ([`assemble_record`](crate::vcf::assemble_record)). The generic path still calls the flank:
//! it wrote what the same reads showed there as an allele of its own record — the flank half of
//! the deletion, or, where the mapper spelled the event another way, a substitution — and it may
//! have written other variants there, another sample's SNP, say. This pass, over the finished
//! records in genome order and before anything prepares them for the filters and the VCF,
//! **merges every generic record that overlaps a widened STR record into it**: one record over
//! the union of their spans, whose alleles are the haplotypes the samples carry there. No `*`
//! allele is written, and no two records claim one base.
//!
//! Each sample's copies are paired record by record (spec §6.2):
//!
//! 1. **an STR copy that reaches the generic record takes the generic copy its own sequence
//!    spells there** — the flank half, whatever form the generic path gave it — and only where
//!    none does, a reference copy, then any; **where the generic call did not say so, the STR
//!    record wins** (owner, 2026-10-10);
//! 2. **the STR copies left pair with the generic copies left**, which is one pairing wherever
//!    either side is homozygous;
//! 3. **where both are heterozygous, the pairing the cohort's settled samples show most** is
//!    taken, and the genotype's quality carries the chance that it is the wrong one.
//!
//! **Records near an STR segment are held**, because a record widened to the left starts before
//! generic records that arrived ahead of it. Anything ending more than [`HOLD_REACH`] bases before
//! the newest record's start can no longer meet a widened record and is released, in the
//! writer's order. A widened record is merged once the stream has passed its end, when no further
//! generic record can overlap it.

use std::collections::HashMap;

use crate::types::{AlleleId, GenomeRegion, Genotype, Phred, Position};
use crate::vcf::{
    FilterVerdict, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
    TractWidening, VcfRecord,
};
use crate::window_coverage::WindowCoverage;

/// **How far before the newest record a held record may still be met by a widened one.** A
/// record is widened by at most the flank its locus holds, and a flank deletion is at most 255
/// bases (`FlankDeletion` counts in a byte); one more for the padding base.
const HOLD_REACH: u64 = 256;

/// What the pass changed, for the run's report.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeDeletionGatheringCounts {
    /// Samples whose STR allele reaching a generic record found the generic copy it spells
    /// there — the two paths agreeing.
    pub samples_in_agreement: u64,
    /// Samples whose STR allele reaching a generic record found no generic copy it spells, and
    /// overrode one — the STR record winning.
    pub samples_overridden: u64,
    /// Generic records merged into an STR record, and so not written on their own.
    pub generic_records_merged: u64,
    /// Samples heterozygous at both records whose pairing came from the cohort (spec §6.2).
    pub samples_paired_from_the_cohort: u64,
    /// Samples a generic record had called that the merged record writes as a no-call, the STR
    /// record having none for them.
    pub generic_calls_lost_to_an_str_no_call: u64,
    /// Generic records left beside the STR record they overlap, because the merged record could
    /// not be built: an allele came out empty with no padding base kept for it, or a haplotype
    /// could not be spelled. Expected to stay zero.
    pub generic_records_left_unmerged: u64,
}

/// One widened STR record, waiting for the generic records it overlaps.
#[derive(Debug)]
struct WidenedTract {
    record: VcfRecord,
    windows: Vec<WindowCoverage>,
    /// The tract's first and last base, 1-based.
    tract_first: u64,
    tract_last: u64,
    /// The generic records that overlap it, in the order they arrived.
    generic: Vec<VcfRecord>,
}

impl WidenedTract {
    fn of(record: VcfRecord, windows: Vec<WindowCoverage>) -> Self {
        let region = record.region();
        let widening = widening_of(&record);
        Self {
            tract_first: region.start.get() + u64::from(widening.left),
            tract_last: region.end.get() - u64::from(widening.right),
            record,
            windows,
            generic: Vec::new(),
        }
    }

    fn overlaps(&self, region: GenomeRegion) -> bool {
        let own = self.record.region();
        region.contig == own.contig && region.start <= own.end && region.end >= own.start
    }

    /// Whether STR allele `allele` changes a base of `region` where the generic path looks:
    /// its reach into the left flank or the right one overlaps it.
    fn reaches(&self, allele: u16, region: GenomeRegion) -> bool {
        let (left, right) = widening_of(&self.record)
            .reach_per_allele
            .get(usize::from(allele))
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
    fn spells(&self, allele: u16, region: GenomeRegion, reference: &[u8], generic: &[u8]) -> bool {
        let ours = &*self.record.alleles()[usize::from(allele)];
        let (start, end) = (region.start.get(), region.end.get());
        let own = self.record.region();
        let (window_start, window_end) = (own.start.get(), own.end.get());
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
}

fn widening_of(record: &VcfRecord) -> &TractWidening {
    record
        .repeat_tract()
        .and_then(|tract| tract.widening())
        .expect("a widened STR record carries its widening")
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
        // A new contig: everything held and waiting is behind it.
        let contig_changed = self
            .held
            .iter()
            .map(|(held, _)| held.region().contig)
            .chain(
                self.widened
                    .iter()
                    .map(|tract| tract.record.region().contig),
            )
            .any(|contig| contig != region.contig);
        if contig_changed {
            released = self.finish();
        }
        // Widened records the stream has passed for good: nothing later can overlap them.
        let (passed, waiting): (Vec<_>, Vec<_>) = std::mem::take(&mut self.widened)
            .into_iter()
            .partition(|tract| tract.record.region().end < region.start);
        self.widened = waiting;
        for tract in passed {
            let merged = merge(tract, &mut self.counts);
            self.held.extend(merged);
        }

        let is_widened = record
            .repeat_tract()
            .is_some_and(|tract| tract.widening().is_some());
        if is_widened {
            let mut tract = WidenedTract::of(record, windows);
            // The generic records already held that it overlaps — its left flank.
            let mut index = 0;
            while index < self.held.len() {
                if !self.held[index].0.is_repeat_tract()
                    && tract.overlaps(self.held[index].0.region())
                {
                    tract.generic.push(self.held.remove(index).0);
                } else {
                    index += 1;
                }
            }
            self.widened.push(tract);
        } else if let Some(tract) = (!record.is_repeat_tract())
            .then(|| self.widened.iter_mut().find(|tract| tract.overlaps(region)))
            .flatten()
        {
            tract.generic.push(record);
        } else {
            self.held.push((record, windows));
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
        for tract in std::mem::take(&mut self.widened) {
            let merged = merge(tract, &mut self.counts);
            self.held.extend(merged);
        }
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

/// **One copy of a sample over the merged record**: the STR allele it carries, and for each
/// generic record merged so far the generic allele written over that record's bases — `0`
/// where the STR allele spells those bases itself, whether as the reference or by reaching them.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Haplotype {
    str_allele: u16,
    generic: Vec<u16>,
}

/// What a sample's column becomes while the generic records are merged one by one.
#[derive(Debug)]
enum SampleHaplotypes {
    /// The STR record has no call for this sample, so the merged record has none.
    NoCall,
    Called {
        copies: Vec<Haplotype>,
        /// The lower of the `GQ`s of the records that added something to this genotype.
        quality: f32,
        /// The probability that every pairing taken from the cohort is the right one.
        pairing_probability: f64,
    },
}

/// **Merge a widened STR record with the generic records that overlap it** (spec §6.1–§6.7),
/// returning what is to be written: the merged record, or — where it would need a padding base
/// the run did not keep — the records as they came.
fn merge(
    tract: WidenedTract,
    counts: &mut EdgeDeletionGatheringCounts,
) -> Vec<(VcfRecord, Vec<WindowCoverage>)> {
    if tract.generic.is_empty() {
        return vec![(tract.record, tract.windows)];
    }
    let mut generic = tract.generic.clone();
    generic.sort_by_key(|record| record.region().start);

    let mut samples: Vec<SampleHaplotypes> = tract
        .record
        .sample_columns()
        .iter()
        .map(|column| match &column.call {
            SampleCall::Called {
                genotype,
                genotype_quality,
            } => SampleHaplotypes::Called {
                copies: genotype
                    .alleles()
                    .iter()
                    .map(|allele| Haplotype {
                        str_allele: allele.get(),
                        generic: Vec::new(),
                    })
                    .collect(),
                quality: genotype_quality.get(),
                pairing_probability: 1.0,
            },
            SampleCall::NoCall => SampleHaplotypes::NoCall,
        })
        .collect();
    let mut tallies = EdgeDeletionGatheringCounts::default();
    for record in &generic {
        pair_with(&tract, record, &mut samples, &mut tallies);
    }

    let Some(merged) = build_merged_record(&tract, &generic, &samples) else {
        counts.generic_records_left_unmerged += generic.len() as u64;
        let mut out = vec![(tract.record, tract.windows)];
        out.extend(generic.into_iter().map(|record| (record, Vec::new())));
        return out;
    };
    counts.samples_in_agreement += tallies.samples_in_agreement;
    counts.samples_overridden += tallies.samples_overridden;
    counts.samples_paired_from_the_cohort += tallies.samples_paired_from_the_cohort;
    counts.generic_calls_lost_to_an_str_no_call += tallies.generic_calls_lost_to_an_str_no_call;
    counts.generic_records_merged += generic.len() as u64;
    vec![(merged, tract.windows)]
}

/// **Pair every sample's copies with one generic record's** (spec §6.2), extending each copy's
/// haplotype by the generic allele written over that record.
fn pair_with(
    tract: &WidenedTract,
    record: &VcfRecord,
    samples: &mut [SampleHaplotypes],
    tallies: &mut EdgeDeletionGatheringCounts,
) {
    let region = record.region();
    let generic_alleles = record.alleles();
    // Pairings settled without the cohort, counted per (haplotype so far, generic allele).
    let mut settled: HashMap<(Haplotype, u16), u32> = HashMap::new();
    // Samples heterozygous on both sides, with the copies each side has left to pair.
    let mut undecided: Vec<(usize, Vec<usize>, Vec<u16>)> = Vec::new();

    for (sample, state) in samples.iter_mut().enumerate() {
        let column = &record.sample_columns()[sample];
        let SampleHaplotypes::Called {
            copies, quality, ..
        } = state
        else {
            if column
                .call
                .genotype()
                .is_some_and(|genotype| genotype.alleles().iter().any(|a| a.get() != 0))
            {
                tallies.generic_calls_lost_to_an_str_no_call += 1;
            }
            continue;
        };
        let (Some(genotype), Some(generic_quality)) =
            (column.call.genotype(), column.call.genotype_quality())
        else {
            // No generic call: each copy is written as its STR allele spells it.
            for copy in copies.iter_mut() {
                copy.generic.push(0);
            }
            continue;
        };
        let mut generic_copies: Vec<u16> = genotype.alleles().iter().map(|a| a.get()).collect();
        let mut taken = vec![false; generic_copies.len()];
        let mut choice: Vec<Option<u16>> = vec![None; copies.len()];

        // Step 1: the copies whose STR allele reaches this record take the generic copy their
        // own sequence spells there; a variant before the reference, then the reference, then
        // anything.
        let mut reached = false;
        let mut disagreed = false;
        for (slot, copy) in copies.iter().enumerate() {
            if !tract.reaches(copy.str_allele, region) {
                continue;
            }
            reached = true;
            let rank = |generic: u16| -> u8 {
                let bases = &generic_alleles[usize::from(generic)];
                let spelled = tract.spells(copy.str_allele, region, &generic_alleles[0], bases);
                match (spelled, generic == 0) {
                    (true, false) => 0,
                    (true, true) => 1,
                    (false, true) => 2,
                    (false, false) => 3,
                }
            };
            let Some(at) = (0..generic_copies.len())
                .filter(|&at| !taken[at])
                .min_by_key(|&at| rank(generic_copies[at]))
            else {
                break;
            };
            disagreed |= rank(generic_copies[at]) >= 2;
            taken[at] = true;
            choice[slot] = Some(0);
        }
        if reached {
            if disagreed {
                tallies.samples_overridden += 1;
            } else {
                tallies.samples_in_agreement += 1;
            }
        }

        // Step 2: the copies left pair with the generic copies left.
        let open_slots: Vec<usize> = (0..copies.len()).filter(|&s| choice[s].is_none()).collect();
        let open_generic: Vec<u16> = generic_copies
            .drain(..)
            .zip(&taken)
            .filter(|(_, taken)| !**taken)
            .map(|(allele, _)| allele)
            .collect();
        if disagreed || open_generic.iter().any(|&allele| allele != 0) {
            *quality = quality.min(generic_quality.get());
        }
        let one_pairing = open_slots
            .windows(2)
            .all(|pair| copies[pair[0]] == copies[pair[1]])
            || open_generic.windows(2).all(|pair| pair[0] == pair[1])
            || open_slots.len() != open_generic.len();
        if one_pairing {
            for (index, &slot) in open_slots.iter().enumerate() {
                let allele = open_generic.get(index).copied().unwrap_or(0);
                *settled.entry((copies[slot].clone(), allele)).or_default() += 1;
                choice[slot] = Some(allele);
            }
            for (copy, chosen) in copies.iter_mut().zip(choice) {
                copy.generic.push(chosen.unwrap_or(0));
            }
        } else {
            for (copy, chosen) in copies.iter_mut().zip(&choice) {
                if let Some(chosen) = chosen {
                    copy.generic.push(*chosen);
                }
            }
            undecided.push((sample, open_slots, open_generic));
        }
    }

    // Step 3: the samples heterozygous on both sides take the pairing the settled samples show
    // most, each pair weighed by its count plus one half.
    for (sample, open_slots, open_generic) in undecided {
        let SampleHaplotypes::Called {
            copies,
            pairing_probability,
            ..
        } = &mut samples[sample]
        else {
            continue;
        };
        let weight = |haplotype: &Haplotype, allele: u16| -> f64 {
            f64::from(
                settled
                    .get(&(haplotype.clone(), allele))
                    .copied()
                    .unwrap_or(0),
            ) + 0.5
        };
        // Every distinct way of handing the generic copies to the open slots, with its weight.
        let mut pairings: Vec<Pairing> = Vec::new();
        for order in permutations(open_generic.len()) {
            let assigned: Vec<u16> = order.iter().map(|&at| open_generic[at]).collect();
            let mut pairs: Vec<(Haplotype, u16)> = open_slots
                .iter()
                .zip(&assigned)
                .map(|(&slot, &allele)| (copies[slot].clone(), allele))
                .collect();
            pairs.sort();
            if pairings.iter().any(|seen| seen.pairs == pairs) {
                continue;
            }
            let weight = pairs
                .iter()
                .map(|(haplotype, allele)| weight(haplotype, *allele))
                .product();
            pairings.push(Pairing {
                pairs,
                assigned,
                weight,
            });
        }
        let total: f64 = pairings.iter().map(|pairing| pairing.weight).sum();
        // The heaviest, the first by allele order at a tie.
        let best = pairings
            .iter()
            .fold(None::<&Pairing>, |best, candidate| match best {
                Some(best) if best.weight >= candidate.weight => Some(best),
                _ => Some(candidate),
            })
            .expect("two heterozygous sides have at least one pairing");
        let assigned = &best.assigned;
        let best = best.weight;
        *pairing_probability *= best / total;
        for (&slot, &allele) in open_slots.iter().zip(assigned) {
            copies[slot].generic.push(allele);
        }
        tallies.samples_paired_from_the_cohort += 1;
    }
}

/// One way of handing a sample's open generic copies to its open STR copies.
struct Pairing {
    /// The (haplotype so far, generic allele) pairs it makes, sorted: two orders making the same
    /// pairs are one pairing.
    pairs: Vec<(Haplotype, u16)>,
    /// The generic allele each open copy gets, in the open copies' order.
    assigned: Vec<u16>,
    /// The product of its pairs' weights.
    weight: f64,
}

/// Every order of `0..n`, the identity first.
fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn extend(prefix: &mut Vec<usize>, used: &mut [bool], out: &mut Vec<Vec<usize>>) {
        if prefix.len() == used.len() {
            out.push(prefix.clone());
            return;
        }
        for at in 0..used.len() {
            if !used[at] {
                used[at] = true;
                prefix.push(at);
                extend(prefix, used, out);
                prefix.pop();
                used[at] = false;
            }
        }
    }
    let mut out = Vec::new();
    extend(&mut Vec::with_capacity(n), &mut vec![false; n], &mut out);
    out
}

/// **The merged record** (spec §6.1, §6.3–§6.6), or `None` where an allele comes out empty and
/// the run kept no base to pad it with.
fn build_merged_record(
    tract: &WidenedTract,
    generic: &[VcfRecord],
    samples: &[SampleHaplotypes],
) -> Option<VcfRecord> {
    let str_record = &tract.record;
    let str_region = str_record.region();
    let widening = widening_of(str_record);

    // The union and its reference: at most one generic record runs past each end of the STR
    // record, as generic records are disjoint.
    let mut start = str_region.start.get();
    let mut end = str_region.end.get();
    let mut reference: Vec<u8> = str_record.alleles()[0].to_vec();
    let mut before: Option<&VcfRecord> = None;
    let mut after: Option<&VcfRecord> = None;
    for record in generic {
        let (g_start, g_end) = (record.region().start.get(), record.region().end.get());
        let bases = &record.alleles()[0];
        if g_start < start {
            let mut extended = bases[..(start - g_start) as usize].to_vec();
            extended.extend_from_slice(&reference);
            reference = extended;
            start = g_start;
            before = Some(record);
        }
        if g_end > end {
            reference.extend_from_slice(&bases[bases.len() - (g_end - end) as usize..]);
            end = g_end;
            after = Some(record);
        }
    }

    // The haplotypes the written genotypes name: the reference, the STR alleles unchanged in
    // their order, then the rest by STR allele and generic alleles.
    let unchanged = |str_allele: u16| Haplotype {
        str_allele,
        generic: vec![0; generic.len()],
    };
    let mut named: Vec<Haplotype> = samples
        .iter()
        .filter_map(|state| match state {
            SampleHaplotypes::Called { copies, .. } => Some(copies.iter().cloned()),
            SampleHaplotypes::NoCall => None,
        })
        .flatten()
        .collect();
    named.sort_by_key(|haplotype| {
        (
            haplotype.generic.iter().any(|&allele| allele != 0),
            haplotype.clone(),
        )
    });
    named.dedup();
    let mut haplotypes = vec![unchanged(0)];
    haplotypes.extend(
        named
            .into_iter()
            .filter(|haplotype| *haplotype != unchanged(0)),
    );

    // Each haplotype's sequence over the union; two haplotypes spelling one sequence are one
    // allele.
    let mut alleles: Vec<Box<[u8]>> = Vec::new();
    let mut origins: Vec<usize> = Vec::new();
    let mut index_of: HashMap<Haplotype, usize> = HashMap::new();
    for haplotype in &haplotypes {
        let bases = spell(tract, generic, &reference, start, end, haplotype)?;
        let index = match alleles.iter().position(|allele| **allele == *bases) {
            Some(index) => index,
            None => {
                alleles.push(bases.into_boxed_slice());
                origins.push(usize::from(haplotype.str_allele));
                alleles.len() - 1
            }
        };
        index_of.insert(haplotype.clone(), index);
    }

    // Bases a generic record adds beyond the STR record that every allele shares — the anchor a
    // generic indel carries — say nothing, and trimmed the record starts and ends where the STR
    // record did. A base trimmed from the front is the one an allele left empty is padded with.
    let mut trimmed_before: Option<u8> = None;
    while start < str_region.start.get()
        && alleles[0].len() > 1
        && alleles
            .iter()
            .all(|allele| allele.first() == alleles[0].first())
    {
        trimmed_before = Some(alleles[0][0]);
        for allele in &mut alleles {
            *allele = allele[1..].into();
        }
        start += 1;
    }
    while end > str_region.end.get()
        && alleles
            .iter()
            .all(|allele| allele.len() > 1 && allele.last() == alleles[0].last())
    {
        for allele in &mut alleles {
            *allele = allele[..allele.len() - 1].into();
        }
        end -= 1;
    }
    let before = before.filter(|_| start < str_region.start.get());
    let after = after.filter(|_| end > str_region.end.get());

    // Per sample: the genotype, its quality, and the STR record's reads on the new alleles.
    let built_on = |origin: usize| -> Vec<usize> {
        (0..alleles.len())
            .filter(|&index| origins[index] == origin)
            .collect()
    };
    let mut columns = Vec::with_capacity(samples.len());
    let mut called_copies = vec![0_u32; alleles.len()];
    for (state, column) in samples.iter().zip(str_record.sample_columns()) {
        let call = match state {
            SampleHaplotypes::NoCall => SampleCall::NoCall,
            SampleHaplotypes::Called {
                copies,
                quality,
                pairing_probability,
            } => {
                let ids: Vec<AlleleId> = copies
                    .iter()
                    .map(|copy| {
                        let index = index_of[copy];
                        called_copies[index] += 1;
                        AlleleId(u16::try_from(index).expect("fewer than 65,536 alleles"))
                    })
                    .collect();
                SampleCall::Called {
                    genotype: Genotype::new(ids),
                    genotype_quality: written_quality(*quality, *pairing_probability),
                }
            }
        };
        let named: Vec<usize> = call
            .genotype()
            .map(|genotype| {
                genotype
                    .alleles()
                    .iter()
                    .map(|a| usize::from(a.get()))
                    .collect()
            })
            .unwrap_or_default();
        let mut reads = vec![0_u32; alleles.len()];
        let mut unexplained = column.read_counts.unexplained_reads();
        for (origin, &count) in column.read_counts.allele_reads().iter().enumerate() {
            let candidates = built_on(origin);
            let target = if candidates.len() == 1 {
                Some(candidates[0])
            } else {
                let mut on_genotype: Vec<usize> = candidates
                    .into_iter()
                    .filter(|index| named.contains(index))
                    .collect();
                on_genotype.dedup();
                (on_genotype.len() == 1).then(|| on_genotype[0])
            };
            match target {
                Some(index) => reads[index] += count,
                None => unexplained += count,
            }
        }
        columns.push(SampleColumn {
            call,
            read_counts: SampleReadCounts::new(reads, unexplained),
        });
    }

    // The cohort-level columns, from the STR allele each new allele is built on.
    let mut expected_copies = vec![0.0; alleles.len()];
    let mut allele_mapq = vec![MapqPool::default(); alleles.len()];
    for origin in 0..str_record.alleles().len() {
        let candidates = built_on(origin);
        if candidates.is_empty() {
            continue;
        }
        let copies: u32 = candidates.iter().map(|&index| called_copies[index]).sum();
        let fitted = str_record.expected_copies()[origin];
        let pool = str_record.allele_mapq()[origin];
        for &index in &candidates {
            expected_copies[index] = if copies == 0 {
                if index == candidates[0] { fitted } else { 0.0 }
            } else {
                fitted * f64::from(called_copies[index]) / f64::from(copies)
            };
            let reads: u64 = columns
                .iter()
                .map(|column| u64::from(column.read_counts.allele_reads()[index]))
                .sum();
            allele_mapq[index] = MapqPool {
                reads,
                mapq_sum: if pool.reads == 0 {
                    0
                } else {
                    (pool.mapq_sum as f64 * reads as f64 / pool.reads as f64).round() as u64
                },
            };
        }
    }

    let padding_base = if alleles.iter().any(|allele| allele.is_empty()) {
        Some(match trimmed_before {
            Some(base) => PaddingBase::Left(base),
            None => padding_base_of(str_record, widening, before, after, start, end)?,
        })
    } else {
        None
    };

    // The site: the larger quality with its own penalties, and the STR record's filter unless
    // it passed and a generic record's did not — counting only the generic records that add an
    // allele of their own. One whose every called allele the STR alleles spell said nothing the
    // STR record does not, and leaves it as it was (spec §6.7).
    let adds_an_allele = |k: usize| haplotypes.iter().any(|haplotype| haplotype.generic[k] != 0);
    let contributing: Vec<&VcfRecord> = (0..generic.len())
        .filter(|&k| adds_an_allele(k))
        .map(|k| &generic[k])
        .collect();
    let best = std::iter::once(str_record)
        .chain(contributing.iter().copied())
        .fold(str_record, |best, record| {
            if record.site_quality.get() > best.site_quality.get() {
                record
            } else {
                best
            }
        });
    let filter = if str_record.filter() == FilterVerdict::Pass {
        contributing
            .iter()
            .map(|record| record.filter())
            .find(|filter| *filter != FilterVerdict::Pass)
            .unwrap_or(FilterVerdict::Pass)
    } else {
        str_record.filter()
    };

    let spelling_length = |origin: usize| {
        str_record.alleles()[origin].len() - usize::from(widening.flank_bases_per_allele[origin])
    };
    let rebuilt = TractWidening {
        left: u16::try_from(tract.tract_first - start).expect("a flank shorter than 65,536 bases"),
        right: u16::try_from(end - tract.tract_last).expect("a flank shorter than 65,536 bases"),
        flank_bases_per_allele: alleles
            .iter()
            .zip(&origins)
            .map(|(allele, &origin)| {
                u16::try_from(allele.len().saturating_sub(spelling_length(origin)))
                    .expect("a flank shorter than 65,536 bases")
            })
            .collect(),
        reach_per_allele: origins
            .iter()
            .map(|&origin| widening.reach_per_allele[origin])
            .collect(),
        base_before: if trimmed_before.is_some() {
            trimmed_before
        } else if start == str_region.start.get() {
            widening.base_before
        } else {
            before
                .and_then(VcfRecord::padding_base)
                .and_then(|base| match base {
                    PaddingBase::Left(base) => Some(base),
                    PaddingBase::Right(_) => None,
                })
        },
    };
    let annotation = str_record
        .repeat_tract()
        .cloned()
        .expect("a widened record is an STR record")
        .rebuilt_over(&origins, rebuilt);

    Some(VcfRecord::new(
        GenomeRegion {
            contig: str_region.contig,
            start: Position(start),
            end: Position(end),
        },
        alleles,
        expected_copies,
        columns,
        allele_mapq,
        padding_base,
        best.site_quality,
        best.artifact_penalties,
        filter,
        Some(annotation),
    ))
}

/// **A haplotype's sequence over the union** `start..=end`: its STR allele, with the reference
/// beyond the STR record and each generic allele it carries written over that record's bases.
/// `None` where the STR allele does not spell the reference over a record it is given a generic
/// allele at, which the pairing never produces.
fn spell(
    tract: &WidenedTract,
    generic: &[VcfRecord],
    reference: &[u8],
    start: u64,
    end: u64,
    haplotype: &Haplotype,
) -> Option<Vec<u8>> {
    let str_region = tract.record.region();
    let head = (str_region.start.get() - start) as usize;
    let tail = (end - str_region.end.get()) as usize;
    let mut bases = reference[..head].to_vec();
    bases.extend_from_slice(&tract.record.alleles()[usize::from(haplotype.str_allele)]);
    bases.extend_from_slice(&reference[reference.len() - tail..]);

    // Left of the tract, positions count from the start; rightmost first, so each replacement
    // leaves the positions before it untouched. Right of it, positions count from the end,
    // nearest the tract first.
    let mut order: Vec<usize> = (0..generic.len()).collect();
    order.sort_by_key(|&k| {
        let at = generic[k].region().start.get();
        if at < tract.tract_first {
            (0, u64::MAX - at)
        } else {
            (1, at)
        }
    });
    for k in order {
        let allele = haplotype.generic[k];
        if allele == 0 {
            continue;
        }
        let record = &generic[k];
        let (g_start, g_end) = (record.region().start.get(), record.region().end.get());
        let length = (g_end - g_start + 1) as usize;
        let from = if g_start < tract.tract_first {
            (g_start - start) as usize
        } else {
            bases.len().checked_sub((end - g_end) as usize + length)?
        };
        if bases.get(from..from + length)? != &*record.alleles()[0] {
            return None;
        }
        bases.splice(
            from..from + length,
            record.alleles()[usize::from(allele)].iter().copied(),
        );
    }
    Some(bases)
}

/// **The padding base of the merged record**: the base before the union, or — at a contig's
/// first base — the one after it, from whichever record holds it; `None` where none does.
fn padding_base_of(
    str_record: &VcfRecord,
    widening: &TractWidening,
    before: Option<&VcfRecord>,
    after: Option<&VcfRecord>,
    start: u64,
    end: u64,
) -> Option<PaddingBase> {
    if start == crate::vcf::FIRST_POSITION_OF_A_CONTIG.get() {
        let last = after.unwrap_or(str_record);
        return match last.padding_base() {
            Some(PaddingBase::Right(base)) if last.region().end.get() == end => {
                Some(PaddingBase::Right(base))
            }
            _ => None,
        };
    }
    match before {
        None => widening.base_before.map(PaddingBase::Left),
        Some(record) => match record.padding_base() {
            Some(PaddingBase::Left(base)) => Some(PaddingBase::Left(base)),
            _ => None,
        },
    }
}

/// **`GQ` of a merged genotype** (spec §6.3): the lower of the contributing `GQ`s, and where a
/// pairing came from the cohort, the chance that it is wrong as well.
fn written_quality(quality: f32, pairing_probability: f64) -> Phred {
    let quality = if pairing_probability < 1.0 {
        let error = crate::float::powf(10.0, -f64::from(quality) / 10.0);
        let right = (1.0 - error) * pairing_probability;
        (-10.0 * crate::float::log10(1.0 - right)) as f32
    } else {
        quality
    };
    Phred::try_new(quality.max(0.0)).unwrap_or(Phred::ZERO)
}

#[cfg(test)]
mod tests;
