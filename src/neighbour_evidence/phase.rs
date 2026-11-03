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
//! A window's first block continues the sample's last written block where the two windows'
//! overlap shows how (§2.4).

use std::collections::BTreeMap;

use super::model::HaplotypePair;
use crate::read_phase::fragments::FragmentAllele;
use crate::read_phase::link::{self, Heterozygote, ReadLink, Relation, score_read_link};

/// **Paths drawn per averaged round besides the round's own** (spec §2.3): with 8, the 64-founder
/// simulated cohort's links inside blocks were wrong 0.35% of the time, against 1.43% with 4
/// (`doc/devel/ng/research/neighbour_evidence/phase_blocks_check_2026-11-02.md`).
pub const EXTRA_PATHS_PER_ROUND: usize = 8;

/// **The most paths a window can record**: one bit a path in a 64-bit mask. Rounds × paths a round
/// must stay within it; `refine::Phasing::new` checks.
pub const MAX_PATHS: usize = u64::BITS as usize;

/// **How far a pair of copies must lean on a record's two alleles to count** (spec §2.2): the lean
/// `w − w′`, with `w = h₀(a)·h₁(b)` and `w′ = h₀(b)·h₁(a)`, at least this in size — for a
/// biallelic record the prototype's `|h₀ − h₁| ≥ 0.5`.
pub const LEAN_MIN: f64 = 0.5;

/// **The share of drawn paths that must agree at a link, and lean at both its records**
/// (`--phase-min-path-agreement`, spec §9): a number in `(0, 1]`, all of them by default. The
/// count of paths it asks for is rounded up, so 0.7 of 36 paths is 26.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathAgreement(f64);

impl PathAgreement {
    /// All paths: the rule of spec §2.3.
    pub const ALL: Self = Self(1.0);

    /// The share `share`, or `None` when it is not in `(0, 1]` (NaN included).
    #[must_use]
    pub fn try_new(share: f64) -> Option<Self> {
        (share > 0.0 && share <= 1.0).then_some(Self(share))
    }

    /// The share.
    #[must_use]
    pub fn get(self) -> f64 {
        self.0
    }

    /// How many of `paths` paths must agree: the share of them, rounded up.
    fn paths_needed(self, paths: u32) -> u32 {
        (self.0 * f64::from(paths)).ceil() as u32
    }
}

impl Default for PathAgreement {
    fn default() -> Self {
        Self::ALL
    }
}

/// The alignment's cost of a marked record whose state disagrees with it (spec §2.1, soft).
pub const DISAGREEMENT_COST: f64 = 1.0;

/// The alignment's cost of a change of state between two marked records (spec §2.1, soft): with
/// 2, a stretch is swapped back only if it holds at least five marked records inside the window,
/// or three at its end.
pub const SWITCH_COST: f64 = 2.0;

/// **One phased record of one sample, as written**: the allele on copy 0 and the allele on copy 1,
/// in the record's alleles as the pass read them, and the phase set (`PS`) it belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WrittenPhase {
    pub copy_0_allele: usize,
    pub copy_1_allele: usize,
    pub phase_set: u64,
}

impl WrittenPhase {
    /// The same record with its two copies swapped.
    fn swapped(self) -> Self {
        Self {
            copy_0_allele: self.copy_1_allele,
            copy_1_allele: self.copy_0_allele,
            ..self
        }
    }
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

/// Each path's bits at one pair entry: bit `p` of `towards_lower` says path `p`'s copy 0 leans to
/// the pair's lower-numbered allele, bit `p` of `leans` that it leans by at least [`LEAN_MIN`].
#[derive(Clone, Copy, Debug, Default)]
struct PathBits {
    towards_lower: u64,
    leans: u64,
}

/// **The lean of two copies on alleles `a` and `b`**: `h₀(a)·h₁(b) − h₀(b)·h₁(a)`, positive when
/// copy 0 carries `a` and copy 1 `b` rather than the reverse.
fn lean(copy_0: &[f64], copy_1: &[f64], a: usize, b: usize) -> f64 {
    copy_0[a] * copy_1[b] - copy_0[b] * copy_1[a]
}

/// A heterozygote the consensus leans on: its site, its alleles `a < b`, the paths' bits there,
/// and whether the consensus's copy 0 carries `a`.
#[derive(Clone, Copy)]
struct Candidate {
    site: usize,
    alleles: (usize, usize),
    bits: PathBits,
    towards_lower: bool,
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
    /// If more than [`MAX_PATHS`] paths are recorded; `refine::Phasing::new` refuses round counts
    /// that would.
    pub(super) fn record_path(&mut self, pairs: &PairLayout, haplotypes: &[HaplotypePair]) {
        assert!(
            (self.paths as usize) < MAX_PATHS,
            "at most {MAX_PATHS} paths a window"
        );
        let bit = 1_u64 << self.paths;
        for (site, [copy_0, copy_1]) in haplotypes.iter().enumerate() {
            let alleles = copy_0.len();
            // A site's pairs are consecutive entries, in the order the two loops visit them.
            let mut entry = pairs.entry(site, 0, 1, alleles.max(2));
            for a in 0..alleles {
                for b in a + 1..alleles {
                    let bits = &mut self.bits[entry];
                    let lean = lean(copy_0, copy_1, a, b);
                    if lean > 0.0 {
                        bits.towards_lower |= bit;
                    }
                    if lean.abs() >= LEAN_MIN {
                        bits.leans |= bit;
                    }
                    entry += 1;
                }
            }
        }
        self.paths += 1;
    }

    /// The consensus's two haplotypes at `site`, `A` probabilities a copy, into the two buffers.
    fn consensus_into(&self, site: usize, copy_0: &mut Vec<f64>, copy_1: &mut Vec<f64>) {
        let start = self.sum_offsets[site];
        let stored = (self.sum_offsets[site + 1] - start) / 2;
        let count = f64::from(self.rounds_added.max(1));
        for (copy, from) in [(copy_0, start), (copy_1, start + stored)] {
            copy.clear();
            copy.extend(
                self.sums[from..from + stored]
                    .iter()
                    .map(|&sum| f64::from(sum) / count),
            );
            // f32 rounding can take the stored ones a hair past one.
            let rest = (1.0 - copy.iter().sum::<f64>()).max(0.0);
            copy.push(rest);
        }
    }

    /// **Lines a round's own haplotypes up with the consensus and adds them** (spec §2.1).
    pub(super) fn add_round(&mut self, haplotypes: &[HaplotypePair]) {
        let swapped = if self.rounds_added == 0 {
            vec![false; haplotypes.len()]
        } else {
            self.align(haplotypes)
        };
        for (site, ([copy_0, copy_1], &swap)) in haplotypes.iter().zip(&swapped).enumerate() {
            let (copy_0, copy_1) = if swap {
                (copy_1, copy_0)
            } else {
                (copy_0, copy_1)
            };
            let start = self.sum_offsets[site];
            let stored = copy_0.len() - 1;
            for (sum, value) in self.sums[start..start + stored].iter_mut().zip(copy_0) {
                *sum += *value as f32;
            }
            for (sum, value) in self.sums[start + stored..start + 2 * stored]
                .iter_mut()
                .zip(copy_1)
            {
                *sum += *value as f32;
            }
        }
        self.rounds_added += 1;
    }

    /// **Whether to swap the round's copies at each site**: the cheapest path of two states over
    /// the sites where both the consensus and the round lean by at least [`LEAN_MIN`] on the allele
    /// the consensus's copies differ most on, ties to "kept"; each site takes the state of the last
    /// such site at or before it, the first such site's before it.
    fn align(&self, haplotypes: &[HaplotypePair]) -> Vec<bool> {
        let (mut consensus_0, mut consensus_1) = (Vec::new(), Vec::new());
        let mut marked: Vec<(usize, bool)> = Vec::new();
        for (site, [copy_0, copy_1]) in haplotypes.iter().enumerate() {
            self.consensus_into(site, &mut consensus_0, &mut consensus_1);
            let (allele, difference) = (0..copy_0.len())
                .map(|allele| (allele, consensus_0[allele] - consensus_1[allele]))
                .fold((0, 0.0_f64), |best, next| {
                    if next.1.abs() > best.1.abs() {
                        next
                    } else {
                        best
                    }
                });
            let round_difference = copy_0[allele] - copy_1[allele];
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

    /// **Phases the sample's heterozygotes in blocks** (spec §2.2, §2.3; read phase spec §2.7).
    /// `heterozygotes` holds, per site, the written genotype's two alleles `(a, b)`, `a < b`, where
    /// it is a heterozygote the pass called, and `None` elsewhere; `positions` the sites'
    /// positions; `segment_bases` the sample's copying length `S`, beyond which no link is kept;
    /// `agreement` the share of paths that must agree at a link; `reads` the sample's fragments and
    /// written `GQ` at each site. Each block is named by its first record's position. Returns one
    /// entry per site, and what the reads did.
    ///
    /// **Reads change the link test between consecutive candidates** (read phase spec §2.7): a kept
    /// population link stands unless the reads contradict it at [`BREAK_MIN`](crate::read_phase::link::BREAK_MIN) or more, which
    /// breaks it; where the population kept no link, reads at [`LINK_MIN`](crate::read_phase::link::LINK_MIN) either way with `GQ` at
    /// least [`link::GQ_MIN`] at both ends link the two, within `S`, and when their relation is the
    /// opposite of the two consensus orientations the rest of the block is swapped as a whole.
    /// Then heterozygotes the consensus does not lean on are attached by reads, between two
    /// consecutive candidates — see [`attach_by_reads`].
    pub(super) fn phase_heterozygotes(
        &self,
        pairs: &PairLayout,
        heterozygotes: &[Option<(usize, usize)>],
        positions: &[u64],
        segment_bases: f64,
        agreement: PathAgreement,
        reads: &SiteReads<'_>,
    ) -> (Vec<Option<WrittenPhase>>, PassReadCounts) {
        let mut out = vec![None; heterozygotes.len()];
        let mut counts = PassReadCounts::default();
        if self.paths == 0 {
            return (out, counts);
        }
        let all_paths = if self.paths as usize == MAX_PATHS {
            u64::MAX
        } else {
            (1_u64 << self.paths) - 1
        };
        let needed = agreement.paths_needed(self.paths);
        let (mut consensus_0, mut consensus_1) = (Vec::new(), Vec::new());
        let candidates: Vec<Candidate> = heterozygotes
            .iter()
            .enumerate()
            .filter_map(|(site, heterozygote)| {
                let (a, b) = (*heterozygote)?;
                self.consensus_into(site, &mut consensus_0, &mut consensus_1);
                let lean = lean(&consensus_0, &consensus_1, a, b);
                (lean.abs() >= LEAN_MIN).then(|| Candidate {
                    site,
                    alleles: (a, b),
                    bits: self.bits[pairs.entry(site, a, b, consensus_0.len())],
                    towards_lower: lean > 0.0,
                })
            })
            .collect();
        let near = |before: usize, after: usize| {
            (positions[after] - positions[before]) as f64 <= segment_bases
        };
        // A link is kept within S where enough paths give the two records the consensus's relation
        // — their two bits' exclusive-or — and lean at both.
        let linked = |before: &Candidate, after: &Candidate| -> bool {
            let same_relation = !(after.bits.towards_lower ^ before.bits.towards_lower);
            let agreeing = if after.towards_lower == before.towards_lower {
                same_relation
            } else {
                !same_relation
            };
            let agree = agreeing & after.bits.leans & before.bits.leans & all_paths;
            near(before.site, after.site) && agree.count_ones() >= needed
        };
        // Each block's members, each with the orientation it is written in: whether copy 0 carries
        // the lower allele. A block's orientation is the consensus's, swapped as a whole where a
        // read-made link said the consensus had it the wrong way round. While blocks are built
        // every member is a candidate, and the open block's last one is kept beside them.
        let mut blocks: Vec<Vec<Member>> = Vec::new();
        let mut last_candidate: Option<&Candidate> = None;
        let mut swapped = false;
        for candidate in &candidates {
            let joined = blocks
                .last()
                .zip(last_candidate)
                .and_then(|(block, last_candidate)| {
                    let last = block[block.len() - 1];
                    let link = link_to_next(
                        linked(last_candidate, candidate),
                        near(last.site, candidate.site),
                        last_candidate,
                        candidate,
                        reads,
                        &mut counts,
                    )?;
                    Some(match link {
                        NextLink::Population => candidate.towards_lower ^ swapped,
                        NextLink::ByReads(relation) => relation.carry(last.towards_lower),
                    })
                });
            match (joined, blocks.last_mut()) {
                (Some(towards_lower), Some(block)) => {
                    swapped = candidate.towards_lower ^ towards_lower;
                    block.push(Member::of(candidate, towards_lower));
                    last_candidate = Some(candidate);
                }
                // A broken link between two records at one position: the second stays unphased,
                // so two blocks never share a PS (spec §2.3).
                (None, Some(block))
                    if positions[block[block.len() - 1].site] == positions[candidate.site] => {}
                _ => {
                    swapped = false;
                    blocks.push(vec![Member::of(candidate, candidate.towards_lower)]);
                    last_candidate = Some(candidate);
                }
            }
        }
        counts.records_attached_by_reads += attach_by_reads(
            &mut blocks,
            &candidates,
            heterozygotes,
            positions,
            segment_bases,
            reads,
        );
        written_phases(&mut blocks, positions, &mut out);
        (out, counts)
    }
}

/// **How a candidate is linked to the block's last candidate**, if it is.
enum NextLink {
    /// The population kept the link and the reads do not contradict it: the candidate takes the
    /// consensus's orientation, swapped as the block is.
    Population,
    /// The population kept no link and the reads made one, with this relation to the last record.
    ByReads(Relation),
}

/// **The link test between a block's last candidate and the next** (read phase spec §2.7), from
/// whether the population kept the link (`population_linked`) and whether the two lie within
/// the copying length (`near`). A kept population link stands unless the reads contradict it at
/// [`BREAK_MIN`](crate::read_phase::link::BREAK_MIN) or more; where the population kept none,
/// reads at [`LINK_MIN`](crate::read_phase::link::LINK_MIN) either way make one, within `S`, with
/// `GQ` at least [`link::GQ_MIN`] at both ends. Counts what the reads did into `counts`.
fn link_to_next(
    population_linked: bool,
    near: bool,
    last: &Candidate,
    next: &Candidate,
    reads: &SiteReads<'_>,
    counts: &mut PassReadCounts,
) -> Option<NextLink> {
    let consensus = Relation::between(last.towards_lower, next.towards_lower);
    let by_reads = reads.link(last.site, last.alleles, next.site, next.alleles);
    let confident = reads.is_confident(last.site) && reads.is_confident(next.site);
    if population_linked {
        if confident && by_reads.relation().is_some() {
            counts.population_links_reached += 1;
            counts.population_links_contradicted +=
                u32::from(by_reads.relation() != Some(consensus));
        }
        if by_reads.contradicts(consensus) {
            counts.population_links_broken += 1;
            None
        } else {
            Some(NextLink::Population)
        }
    } else if confident && near {
        let relation = by_reads.relation()?;
        counts.links_made_by_reads += 1;
        Some(NextLink::ByReads(relation))
    } else {
        None
    }
}

/// **The blocks of two or more, written**: each member's alleles in its orientation, copy 0 first,
/// and the block's phase set, its first record's position.
fn written_phases(blocks: &mut [Vec<Member>], positions: &[u64], out: &mut [Option<WrittenPhase>]) {
    for block in blocks.iter_mut().filter(|block| block.len() >= 2) {
        block.sort_by_key(|member| member.site);
        let phase_set = positions[block[0].site];
        for member in block.iter() {
            let (a, b) = member.alleles;
            let (copy_0_allele, copy_1_allele) = if member.towards_lower { (a, b) } else { (b, a) };
            out[member.site] = Some(WrittenPhase {
                copy_0_allele,
                copy_1_allele,
                phase_set,
            });
        }
    }
}

/// **One record of a block as it will be written**: the heterozygote, and whether copy 0 carries
/// its lower allele. A candidate's consensus data stay on its [`Candidate`]; a record attached by
/// reads never had any.
#[derive(Clone, Copy)]
struct Member {
    site: usize,
    alleles: (usize, usize),
    towards_lower: bool,
}

impl Member {
    fn of(candidate: &Candidate, towards_lower: bool) -> Self {
        Self {
            site: candidate.site,
            alleles: candidate.alleles,
            towards_lower,
        }
    }
}

/// **A sample's reads at a window's sites** (read phase spec §2.7): at each site its fragments, in
/// fragment order, and the `GQ` its column is written with where it is a core site the pass
/// called. Empty — no fragments and no `GQ` anywhere — for a run without reads, which phases as
/// the population alone does.
#[derive(Default)]
pub(crate) struct SiteReads<'a> {
    pub(crate) fragments: Vec<&'a [FragmentAllele]>,
    pub(crate) genotype_qualities: Vec<Option<u32>>,
}

impl SiteReads<'_> {
    /// What the fragments say about the link between two sites' heterozygotes.
    fn link(
        &self,
        one: usize,
        one_alleles: (usize, usize),
        other: usize,
        other_alleles: (usize, usize),
    ) -> ReadLink {
        let at = |site: usize| self.fragments.get(site).copied().unwrap_or(&[]);
        let heterozygote = |(a, b): (usize, usize)| {
            // PANIC-FREE: the alleles index the line's alleles, which calling numbers as
            // `AlleleId`, a `u16`; and the pass phases only a heterozygote, two different alleles.
            let allele = |allele: usize| u16::try_from(allele).expect("an allele index fits a u16");
            Heterozygote::of(allele(a), allele(b)).expect("a heterozygote's two alleles differ")
        };
        score_read_link(
            at(one),
            heterozygote(one_alleles),
            at(other),
            heterozygote(other_alleles),
        )
    }

    /// Whether a link may rest on reads at `site`: its written `GQ` is at least
    /// [`link::GQ_MIN`], the floor for a link that rests on reads.
    fn is_confident(&self, site: usize) -> bool {
        self.genotype_qualities
            .get(site)
            .copied()
            .flatten()
            .is_some_and(|quality| quality >= u32::from(link::GQ_MIN))
    }
}

/// **What reads did to one sample's phase in one window** (read phase spec §2.7, §2.8).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PassReadCounts {
    /// Population links reads reached at [`LINK_MIN`](crate::read_phase::link::LINK_MIN) or more, with `GQ` at least [`link::GQ_MIN`] at
    /// both ends: the read check's denominator (§2.8).
    pub population_links_reached: u32,
    /// Of those, the ones the reads' relation contradicted: an estimate of how often population
    /// links are wrong that does not come from the model's draws, for links under a kilobase.
    pub population_links_contradicted: u32,
    /// Population links broken because reads contradicted them at [`BREAK_MIN`](crate::read_phase::link::BREAK_MIN) or more.
    pub population_links_broken: u32,
    /// Links the population did not keep and reads made.
    pub links_made_by_reads: u32,
    /// Heterozygotes the consensus did not lean on, attached to a block by reads.
    pub records_attached_by_reads: u32,
}

impl std::ops::AddAssign for PassReadCounts {
    /// **Destructured**, so a count added to the struct fails to compile here rather than staying
    /// zero in every total.
    fn add_assign(&mut self, other: Self) {
        let Self {
            population_links_reached,
            population_links_contradicted,
            population_links_broken,
            links_made_by_reads,
            records_attached_by_reads,
        } = other;
        self.population_links_reached += population_links_reached;
        self.population_links_contradicted += population_links_contradicted;
        self.population_links_broken += population_links_broken;
        self.links_made_by_reads += links_made_by_reads;
        self.records_attached_by_reads += records_attached_by_reads;
    }
}

/// **Attach by reads the heterozygotes the consensus does not lean on** (read phase spec §2.7),
/// returning how many were attached.
///
/// Between two consecutive candidates `p` and `q`, the heterozygotes that are not candidates and
/// whose written `GQ` is at least [`link::GQ_MIN`] are taken in order as a chain from `p`, itself at
/// that `GQ`: each is attached when reads at [`LINK_MIN`](crate::read_phase::link::LINK_MIN) link it to the chain's previous record,
/// with the reads' relation, and the chain stops at the first that is not. Where `p` and `q` are in
/// one block, the chain's last record is tested against `q`, and a contradiction at [`BREAK_MIN`](crate::read_phase::link::BREAK_MIN)
/// of the relation the block implies leaves the whole chain unphased. Attaching never changes an
/// existing link or a population orientation.
fn attach_by_reads(
    blocks: &mut [Vec<Member>],
    candidates: &[Candidate],
    heterozygotes: &[Option<(usize, usize)>],
    positions: &[u64],
    segment_bases: f64,
    reads: &SiteReads<'_>,
) -> u32 {
    // Each member's block and written orientation, by site; and which sites are candidates.
    let mut block_of: Vec<Option<(usize, bool)>> = vec![None; heterozygotes.len()];
    for (index, block) in blocks.iter().enumerate() {
        for member in block {
            block_of[member.site] = Some((index, member.towards_lower));
        }
    }
    let mut is_candidate = vec![false; heterozygotes.len()];
    for candidate in candidates {
        is_candidate[candidate.site] = true;
    }
    let mut attached_count = 0;
    for pair in candidates.windows(2) {
        let (p, q) = (pair[0], pair[1]);
        let Some((block, p_towards_lower)) = block_of[p.site] else {
            continue;
        };
        if !reads.is_confident(p.site) {
            continue;
        }
        let mut chain: Vec<Member> = Vec::new();
        let mut previous = Member {
            site: p.site,
            alleles: p.alleles,
            towards_lower: p_towards_lower,
        };
        for site in p.site + 1..q.site {
            let Some(alleles) = heterozygotes[site] else {
                continue;
            };
            if is_candidate[site] || !reads.is_confident(site) {
                continue;
            }
            let near = (positions[site] - positions[previous.site]) as f64 <= segment_bases;
            let Some(relation) = reads
                .link(previous.site, previous.alleles, site, alleles)
                .relation()
                .filter(|_| near)
            else {
                break;
            };
            previous = Member {
                site,
                alleles,
                towards_lower: relation.carry(previous.towards_lower),
            };
            chain.push(previous);
        }
        if let (Some(last), Some((q_block, q_towards_lower))) = (chain.last(), block_of[q.site])
            && q_block == block
        {
            let implied = Relation::between(last.towards_lower, q_towards_lower);
            if reads
                .link(last.site, last.alleles, q.site, q.alleles)
                .contradicts(implied)
            {
                chain.clear();
            }
        }
        attached_count += chain.len() as u32;
        blocks[block].extend(chain);
    }
    attached_count
}

/// **How many shared records join two windows' blocks** (spec §2.4, soft): the records of the
/// last written block nearest the core boundary that the new window also phases. With 3, no join
/// was wrong in 53 on four simulated cohorts; 2 added 11 joins, none wrong
/// (`phase_blocks_check_2026-11-02.md`).
pub const SHARED_RECORDS_TO_JOIN: usize = 3;

/// **One record of a sample's last written block, as the next window sees it**: the record's index
/// in the window, and the alleles the earlier window wrote on copy 0 and copy 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JoinRecord {
    pub record: usize,
    pub copy_0_allele: usize,
    pub copy_1_allele: usize,
}

/// **A sample's last written block** (spec §2.4): its phase set and its records inside the
/// window's left margin, ascending.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LastBlock {
    pub phase_set: u64,
    pub records: Vec<JoinRecord>,
}

/// **Joins a window's blocks to the last written block and names them** (spec §2.4, §2.5), from
/// `phased` — [`SamplePhase::phase_heterozygotes`] over the window's core and the last written
/// block's margin records. `in_core(site)` says which sites the window writes; `site_of(record)`
/// turns a [`JoinRecord::record`] into a site.
///
/// The window's block that holds its first phased core record continues `last` — keeps its phase
/// set, its orientation swapped as a whole if need be — when the [`SHARED_RECORDS_TO_JOIN`] of
/// `last`'s records nearest the core that the window phases all lie in that block, and the two
/// windows' relative orientation is the same at all of them. Every other block is named by its
/// first core record's position, and a block with one core record is written unphased. Returns
/// one entry per site, `None` outside the core.
pub(super) fn join_and_name_core_blocks(
    phased: &[Option<WrittenPhase>],
    in_core: impl Fn(usize) -> bool,
    positions: &[u64],
    last: Option<&LastBlock>,
    site_of: impl Fn(usize) -> usize,
) -> Vec<Option<WrittenPhase>> {
    // Each block's core records, keyed by the window's name for the block.
    let mut core_blocks: BTreeMap<u64, Vec<(usize, WrittenPhase)>> = BTreeMap::new();
    for (site, phase) in phased.iter().enumerate() {
        if let (true, Some(phase)) = (in_core(site), phase) {
            core_blocks
                .entry(phase.phase_set)
                .or_default()
                .push((site, *phase));
        }
    }
    let first_core_block = phased
        .iter()
        .enumerate()
        .find_map(|(site, phase)| phase.filter(|_| in_core(site)))
        .map(|phase| phase.phase_set);
    // The join: the block it continues in this window, whether it is swapped, and the phase set.
    let join: Option<(u64, bool, u64)> = match (last, first_core_block) {
        (Some(last), Some(block)) => {
            let shared: Vec<(WrittenPhase, usize)> = last
                .records
                .iter()
                .filter_map(|record| {
                    phased[site_of(record.record)].map(|phase| (phase, record.copy_0_allele))
                })
                .collect();
            let nearest = &shared[shared.len().saturating_sub(SHARED_RECORDS_TO_JOIN)..];
            let in_block = nearest.iter().all(|(phase, _)| phase.phase_set == block);
            match nearest {
                [(phase, written), rest @ ..]
                    if nearest.len() == SHARED_RECORDS_TO_JOIN && in_block =>
                {
                    let same_way = phase.copy_0_allele == *written;
                    rest.iter()
                        .all(|(phase, written)| (phase.copy_0_allele == *written) == same_way)
                        .then_some((block, !same_way, last.phase_set))
                }
                _ => None,
            }
        }
        _ => None,
    };
    let mut out = vec![None; phased.len()];
    for (block, members) in core_blocks {
        match join {
            Some((joined, swap, phase_set)) if joined == block => {
                for (site, phase) in members {
                    let phase = if swap { phase.swapped() } else { phase };
                    out[site] = Some(WrittenPhase { phase_set, ..phase });
                }
            }
            _ if members.len() >= 2 => {
                let phase_set = positions[members[0].0];
                for (site, phase) in members {
                    out[site] = Some(WrittenPhase { phase_set, ..phase });
                }
            }
            _ => {}
        }
    }
    out
}

#[cfg(test)]
mod tests;
