//! **The phase the rounds estimate, in blocks of confident phase**
//! (`doc/devel/ng/spec/neighbour_evidence_phase_output.md`).
//!
//! Every averaged round draws its own path and [`EXTRA_PATHS_PER_ROUND`] more; each path gives a
//! sample two haplotypes. Two things are kept per sample across the rounds:
//!
//! - **the consensus** (§2.1): the round's own haplotypes, lined up with the consensus so far by
//!   the cheapest two-state path of "kept" and "swapped", added to running sums;
//! - **each path's two bits at every allele pair of every site** (§2.3): which way the path's two
//!   copies lean on the pair — copy 0 towards the pair's first allele or its second — and whether
//!   they lean by at least [`LEAN_MIN`]. A path's relative orientation of two records does not
//!   depend on which copy it calls 0, so the link test reads the bits without lining anything up.
//!
//! After the last round a sample's heterozygous calls whose consensus leans by at least
//! [`LEAN_MIN`] are linked in order, a link kept when enough paths give the two records the
//! consensus's relative orientation and lean at both; blocks of one record are written unphased.

use super::model::HaplotypePair;

/// **Paths drawn per averaged round besides the round's own** (spec §2.3): with 8, the 64-founder
/// simulated cohort's links inside blocks were wrong 0.35% of the time, against 1.43% with 4
/// (`doc/devel/ng/research/neighbour_evidence/phase_blocks_check_2026-11-02.md`).
pub const EXTRA_PATHS_PER_ROUND: usize = 8;

/// **How far a pair of copies must lean on a record's two alleles to count** (spec §2.2): the lean
/// `w − w′`, with `w = h₀(a)·h₁(b)` and `w′ = h₀(b)·h₁(a)`, at least this in size — for a
/// biallelic record the prototype's `|h₀ − h₁| ≥ 0.5`.
pub const LEAN_MIN: f64 = 0.5;

/// **The share of paths that must agree, and lean, at a link by default** (spec §9): all of them.
pub const MIN_PATH_AGREEMENT: f64 = 1.0;

/// The alignment's cost of a record whose state disagrees with it, and of a change of state
/// (spec §2.1, soft).
const DISAGREEMENT_COST: f64 = 1.0;
const SWITCH_COST: f64 = 2.0;

/// **One phased record of one sample, as written**: copy 0's allele and copy 1's, in the record's
/// alleles as the pass read them, and the block's `PS`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WrittenPhase {
    pub first: usize,
    pub second: usize,
    pub block: u64,
}

/// **Where each site's allele pairs sit** among a sample's pair entries: a site of `A` alleles has
/// `A·(A − 1)/2`, the pairs `(a, b)`, `a < b`, in order.
pub(super) struct PairLayout {
    offsets: Vec<usize>,
}

impl PairLayout {
    pub(super) fn new(allele_counts: impl Iterator<Item = usize>) -> Self {
        let mut offsets = vec![0];
        for alleles in allele_counts {
            let last = offsets.last().copied().unwrap_or(0);
            offsets.push(last + alleles * (alleles - 1) / 2);
        }
        Self { offsets }
    }

    fn entries(&self) -> usize {
        self.offsets.last().copied().unwrap_or(0)
    }

    /// The entry of pair `(a, b)`, `a < b`, at `site` of `alleles` alleles.
    fn entry(&self, site: usize, a: usize, b: usize, alleles: usize) -> usize {
        debug_assert!(a < b && b < alleles);
        // Pairs with a smaller first allele come first: (alleles − 1) + (alleles − 2) + … of them.
        self.offsets[site] + a * (2 * alleles - a - 1) / 2 + (b - a - 1)
    }
}

/// Each path's bits at one pair entry: bit `p` of `towards_first` says path `p`'s copy 0 leans to
/// the pair's first allele, bit `p` of `leans` that it leans by at least [`LEAN_MIN`].
#[derive(Clone, Copy, Debug, Default)]
struct PathBits {
    towards_first: u64,
    leans: u64,
}

/// **The lean of two copies on alleles `a` and `b`**: `h₀(a)·h₁(b) − h₀(b)·h₁(a)`, positive when
/// copy 0 carries `a` and copy 1 `b` rather than the reverse.
fn lean(first: &[f64], second: &[f64], a: usize, b: usize) -> f64 {
    first[a] * second[b] - first[b] * second[a]
}

/// **One sample's phase across a window's averaged rounds.**
pub(super) struct SamplePhase {
    /// Sums of the lined-up haplotypes, per site copy 0's probabilities of the first `A − 1`
    /// alleles then copy 1's, as `f32` — the store's layout; the last allele's is the count less the
    /// rest.
    sums: Vec<f32>,
    /// Where each site's `2·(A − 1)` sums start; one past the last at the end.
    sum_offsets: Vec<usize>,
    rounds_added: u32,
    bits: Vec<PathBits>,
    paths: u32,
}

impl SamplePhase {
    /// A sample's phase over sites of these allele counts, before any round.
    pub(super) fn new(allele_counts: &[usize], pairs: &PairLayout) -> Self {
        let mut sum_offsets = Vec::with_capacity(allele_counts.len() + 1);
        sum_offsets.push(0);
        for &alleles in allele_counts {
            let last = sum_offsets.last().copied().unwrap_or(0);
            sum_offsets.push(last + 2 * (alleles - 1));
        }
        Self {
            sums: vec![0.0; sum_offsets.last().copied().unwrap_or(0)],
            sum_offsets,
            rounds_added: 0,
            bits: vec![PathBits::default(); pairs.entries()],
            paths: 0,
        }
    }

    /// **Records one path's bits** (spec §2.3).
    ///
    /// # Panics
    ///
    /// If more than 64 paths are recorded.
    pub(super) fn record_path(&mut self, pairs: &PairLayout, haplotypes: &[HaplotypePair]) {
        assert!(self.paths < 64, "at most 64 paths a window");
        let bit = 1_u64 << self.paths;
        for (site, [first, second]) in haplotypes.iter().enumerate() {
            let alleles = first.len();
            for a in 0..alleles {
                for b in a + 1..alleles {
                    let entry = &mut self.bits[pairs.entry(site, a, b, alleles)];
                    let lean = lean(first, second, a, b);
                    if lean > 0.0 {
                        entry.towards_first |= bit;
                    }
                    if lean.abs() >= LEAN_MIN {
                        entry.leans |= bit;
                    }
                }
            }
        }
        self.paths += 1;
    }

    /// The consensus's two haplotypes at `site`, `A` probabilities a copy.
    fn consensus(&self, site: usize) -> (Vec<f64>, Vec<f64>) {
        let start = self.sum_offsets[site];
        let stored = (self.sum_offsets[site + 1] - start) / 2;
        let count = f64::from(self.rounds_added.max(1));
        let copy = |from: usize| -> Vec<f64> {
            let mut copy: Vec<f64> = self.sums[from..from + stored]
                .iter()
                .map(|&sum| f64::from(sum) / count)
                .collect();
            // f32 rounding can take the stored ones a hair past one.
            copy.push((1.0 - copy.iter().sum::<f64>()).max(0.0));
            copy
        };
        (copy(start), copy(start + stored))
    }

    /// **Lines a round's own haplotypes up with the consensus and adds them** (spec §2.1).
    pub(super) fn add_round(&mut self, haplotypes: &[HaplotypePair]) {
        let swapped = if self.rounds_added == 0 {
            vec![false; haplotypes.len()]
        } else {
            self.alignment(haplotypes)
        };
        for (site, ([first, second], &swap)) in haplotypes.iter().zip(&swapped).enumerate() {
            let (first, second) = if swap {
                (second, first)
            } else {
                (first, second)
            };
            let start = self.sum_offsets[site];
            let stored = first.len() - 1;
            for (sum, value) in self.sums[start..start + stored].iter_mut().zip(first) {
                *sum += *value as f32;
            }
            for (sum, value) in self.sums[start + stored..start + 2 * stored]
                .iter_mut()
                .zip(second)
            {
                *sum += *value as f32;
            }
        }
        self.rounds_added += 1;
    }

    /// Whether to swap the round's copies at each site: the cheapest path of two states over the
    /// sites where both the consensus and the round lean by at least [`LEAN_MIN`] on the allele the
    /// consensus's copies differ most on, ties to "kept"; each site takes the state of the last such
    /// site at or before it, the first such site's before it.
    fn alignment(&self, haplotypes: &[HaplotypePair]) -> Vec<bool> {
        let mut marked: Vec<(usize, bool)> = Vec::new();
        for (site, [first, second]) in haplotypes.iter().enumerate() {
            let (consensus_first, consensus_second) = self.consensus(site);
            let (allele, difference) = (0..first.len())
                .map(|allele| (allele, consensus_first[allele] - consensus_second[allele]))
                .fold((0, 0.0_f64), |best, next| {
                    if next.1.abs() > best.1.abs() {
                        next
                    } else {
                        best
                    }
                });
            let round_difference = first[allele] - second[allele];
            if difference.abs() >= LEAN_MIN && round_difference.abs() >= LEAN_MIN {
                marked.push((site, (difference > 0.0) == (round_difference > 0.0)));
            }
        }
        let mut swapped = vec![false; haplotypes.len()];
        let Some(&(_, first_agrees)) = marked.first() else {
            return swapped;
        };
        // Viterbi over "kept" (0) and "swapped" (1).
        let cost_at = |agrees: bool| -> [f64; 2] {
            if agrees {
                [0.0, DISAGREEMENT_COST]
            } else {
                [DISAGREEMENT_COST, 0.0]
            }
        };
        let mut cost = cost_at(first_agrees);
        let mut changed_into: Vec<[bool; 2]> = Vec::with_capacity(marked.len());
        for &(_, agrees) in &marked[1..] {
            let mut next = [0.0; 2];
            let mut change = [false; 2];
            for state in 0..2 {
                let stay = cost[state];
                let switch = cost[1 - state] + SWITCH_COST;
                change[state] = switch < stay;
                next[state] = stay.min(switch) + cost_at(agrees)[state];
            }
            changed_into.push(change);
            cost = next;
        }
        let mut state = usize::from(cost[1] < cost[0]);
        let mut states = vec![0; marked.len()];
        states[marked.len() - 1] = state;
        for (index, change) in changed_into.iter().enumerate().rev() {
            if change[state] {
                state = 1 - state;
            }
            states[index] = state;
        }
        let mut current = states[0];
        let mut next_marked = 0;
        for (site, swap) in swapped.iter_mut().enumerate() {
            while next_marked < marked.len() && marked[next_marked].0 <= site {
                current = states[next_marked];
                next_marked += 1;
            }
            *swap = current == 1;
        }
        swapped
    }

    /// **The sample's phased records**, given the candidates: per site, the written genotype's two
    /// alleles `(a, b)`, `a < b`, where it is a heterozygote the pass called, and `None` elsewhere;
    /// `positions` the sites' positions; `segment_bases` the sample's copying length `S`, beyond
    /// which no link is kept; `min_agreement` the share of paths that must agree at a link. Returns
    /// one entry per site.
    pub(super) fn phased(
        &self,
        pairs: &PairLayout,
        heterozygotes: &[Option<(usize, usize)>],
        positions: &[u64],
        segment_bases: f64,
        min_agreement: f64,
    ) -> Vec<Option<WrittenPhase>> {
        let mut out = vec![None; heterozygotes.len()];
        if self.paths == 0 {
            return out;
        }
        let all_paths = if self.paths == 64 {
            u64::MAX
        } else {
            (1_u64 << self.paths) - 1
        };
        let needed = (min_agreement * f64::from(self.paths)).ceil() as u32;
        // Candidates: the consensus leans on the genotype's pair. `towards_first` is the consensus's
        // orientation, copy 0 carrying the pair's first allele.
        let mut candidates: Vec<(usize, PathBits, bool)> = Vec::new();
        for (site, heterozygote) in heterozygotes.iter().enumerate() {
            let Some((a, b)) = *heterozygote else {
                continue;
            };
            let (first, second) = self.consensus(site);
            let lean = lean(&first, &second, a, b);
            if lean.abs() >= LEAN_MIN {
                let entry = self.bits[pairs.entry(site, a, b, first.len())];
                candidates.push((site, entry, lean > 0.0));
            }
        }
        let mut blocks: Vec<Vec<(usize, bool)>> = Vec::new();
        let mut current: Vec<(usize, bool)> = Vec::new();
        let mut previous: Option<(usize, PathBits, bool)> = None;
        for &(site, bits, towards_first) in &candidates {
            let kept = previous.is_some_and(|(before, before_bits, before_towards)| {
                let near = (positions[site] - positions[before]) as f64 <= segment_bases;
                let same_relation = towards_first == before_towards;
                // A path agrees where its two bits' relation is the consensus's, and leans at both.
                let relation = !(bits.towards_first ^ before_bits.towards_first);
                let agreeing = if same_relation { relation } else { !relation };
                let agree = agreeing & bits.leans & before_bits.leans & all_paths;
                near && agree.count_ones() >= needed
            });
            if previous.is_some()
                && !kept
                && previous.is_some_and(|(before, _, _)| positions[before] == positions[site])
            {
                // A broken link between two records at one position: the second stays unphased,
                // so two blocks never share a PS (spec §2.3).
                continue;
            }
            if !kept && !current.is_empty() {
                blocks.push(std::mem::take(&mut current));
            }
            current.push((site, towards_first));
            previous = Some((site, bits, towards_first));
        }
        if !current.is_empty() {
            blocks.push(current);
        }
        for block in blocks.iter().filter(|block| block.len() >= 2) {
            let ps = positions[block[0].0];
            for &(site, towards_first) in block {
                // PANIC-FREE: every candidate site has a heterozygote.
                let (a, b) = heterozygotes[site].expect("a candidate is a heterozygote");
                let (first, second) = if towards_first { (a, b) } else { (b, a) };
                out[site] = Some(WrittenPhase {
                    first,
                    second,
                    block: ps,
                });
            }
        }
        out
    }
}

/// **How many shared records join two windows' blocks** (spec §2.4, soft): the records of the
/// last written block nearest the core boundary that the new window also phases. With 3, no join
/// was wrong in 57 on four simulated cohorts; 2 added 13 joins, none wrong
/// (`phase_blocks_check_2026-11-02.md`).
pub const JOIN_RECORDS: usize = 3;

/// **One record of a sample's last written block, as the next window sees it**: the record's index
/// in the window, and the allele the earlier window wrote on copy 0, and on copy 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JoinRecord {
    pub record: usize,
    pub first: usize,
    pub second: usize,
}

/// **A sample's last written block** (spec §2.4): its `PS` and its records inside the window's
/// left margin, ascending.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LastBlock {
    pub block: u64,
    pub records: Vec<JoinRecord>,
}

/// **The phase a window writes in its core** (spec §2.4, §2.5), from `phased` — [`SamplePhase::phased`]
/// over the window's core and the last written block's margin records, each named by its first
/// record's position. `in_core(site)` says which sites the window writes; `last` is the last
/// written block, its `PS` and, per margin site of it, the allele written on copy 0.
///
/// The window's block that holds its first phased core record continues the last written block —
/// keeps its `PS`, its orientation swapped as a whole if need be — when the last
/// [`JOIN_RECORDS`] of the block's margin records the window phases all lie in that block and the
/// two windows' relative orientation is the same at all of them. Every other block is named by
/// its first core record's position, and a block with one core record is written unphased.
pub(super) fn written_in_core(
    phased: &[Option<WrittenPhase>],
    in_core: impl Fn(usize) -> bool,
    positions: &[u64],
    last: Option<(u64, &[(usize, usize)])>,
) -> Vec<Option<WrittenPhase>> {
    let first_core = (0..phased.len()).find(|&site| in_core(site) && phased[site].is_some());
    let joined: Option<(u64, bool)> = match (last, first_core) {
        (Some((_, written)), Some(first_core)) => {
            // PANIC-FREE: `first_core` is a phased site.
            let block = phased[first_core].expect("phased").block;
            let shared: Vec<(WrittenPhase, usize)> = written
                .iter()
                .filter_map(|&(site, first)| phased[site].map(|phase| (phase, first)))
                .collect();
            let nearest = &shared[shared.len().saturating_sub(JOIN_RECORDS)..];
            let same: Vec<bool> = nearest
                .iter()
                .map(|(phase, first)| phase.first == *first)
                .collect();
            (nearest.len() == JOIN_RECORDS
                && nearest.iter().all(|(phase, _)| phase.block == block)
                && same.iter().all(|&same_way| same_way == same[0]))
            .then(|| (block, !same[0]))
        }
        _ => None,
    };
    let mut core_members: std::collections::HashMap<u64, Vec<usize>> = Default::default();
    for site in (0..phased.len()).filter(|&site| in_core(site)) {
        if let Some(phase) = phased[site] {
            core_members.entry(phase.block).or_default().push(site);
        }
    }
    let mut out = vec![None; phased.len()];
    for (block, sites) in core_members {
        match joined {
            Some((joined_block, swap)) if joined_block == block => {
                // PANIC-FREE: a join needs a last written block.
                let (written_block, _) = last.expect("a join has a last block");
                for site in sites {
                    let phase = phased[site].expect("a member is phased");
                    let (first, second) = if swap {
                        (phase.second, phase.first)
                    } else {
                        (phase.first, phase.second)
                    };
                    out[site] = Some(WrittenPhase {
                        first,
                        second,
                        block: written_block,
                    });
                }
            }
            _ if sites.len() >= 2 => {
                let named = positions[sites[0]];
                for site in sites {
                    out[site] = phased[site].map(|phase| WrittenPhase {
                        block: named,
                        ..phase
                    });
                }
            }
            _ => {}
        }
    }
    out
}

#[cfg(test)]
mod tests;
