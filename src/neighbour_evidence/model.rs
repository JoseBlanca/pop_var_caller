//! **One sample's walk along a window** — the copying model with an autozygous state
//! (`doc/devel/ng/spec/neighbour_evidence_pass.md` §3).
//!
//! At each record the hidden state is either *the two copies copy templates `k₁` and `k₂`* (K²
//! states) or *the two copies are one, copying template `k`* (K states). A forward and a backward
//! walk give each state's probability given the whole window; the genotype posterior is the
//! sample's own likelihood weighted by what the states say the alleles are (§3.6).
//!
//! Diploid only: the pair state is two copies (§5). The arithmetic is the prototype's
//! (`doc/devel/ng/research/neighbour_evidence/hmm.py`, `haplotypes`), checked against it to 1e-9 by
//! this module's tests.

use super::random::{Uniforms, pick, pick_at, pick_in_running};
use crate::calling::genotype_table::{GenotypeIdx, GenotypeTable};
use crate::float;
use crate::types::Ploidy;
use std::collections::HashMap;

/// **The numbers a walk runs with** (spec §3.4, §3.5).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WalkParameters {
    /// Mean length, in bases, of a stretch one copy copies from one template (`S`).
    pub segment_bases: f64,
    /// Mean length, in bases, of an autozygous run (`R`).
    pub run_bases: f64,
    /// Long-run share of the sample that is autozygous (`F`), strictly between 0 and 1.
    pub autozygous_share: f64,
    /// The chance a template's allele is not the one the copy carries (`ε`).
    pub mismatch: f64,
    /// **The chance that a copy changing template lands on the last one** — the background, which
    /// the window solve puts last (spec §3.3, §3.5) — the rest shared alike by the others. `None`:
    /// every template alike, `1/K`. Ignored when the walk has one template.
    pub background_landing_share: Option<f64>,
}

impl WalkParameters {
    /// **Where a copy that changes template lands**: one weight a template, summing to one.
    /// `template_weights`, when given, are the relative weights of the first `templates − 1` — the
    /// sample templates (spec §3.3) — on any positive scale; without them those weigh alike. The
    /// last template, the background, takes `background_landing_share`, or without one the mean
    /// of the others' relative weights.
    ///
    /// # Panics
    ///
    /// If `template_weights` is not one weight for each template but the last.
    fn landing_weights(&self, templates: usize, template_weights: Option<&[f64]>) -> Vec<f64> {
        if templates == 1 {
            return vec![1.0];
        }
        let mut weights: Vec<f64> = match template_weights {
            Some(given) => {
                assert_eq!(
                    given.len(),
                    templates - 1,
                    "one landing weight for each template but the background"
                );
                given.to_vec()
            }
            None => vec![1.0; templates - 1],
        };
        let others: f64 = weights.iter().sum();
        match self.background_landing_share {
            Some(background) => {
                for weight in &mut weights {
                    *weight *= (1.0 - background) / others;
                }
                weights.push(background);
            }
            None => {
                weights.push(others / (templates - 1) as f64);
                let total: f64 = weights.iter().sum();
                for weight in &mut weights {
                    *weight /= total;
                }
            }
        }
        weights
    }
}

/// **One record as one sample's walk sees it.**
#[derive(Clone, Copy, Debug)]
pub struct WalkSite<'a> {
    /// The record's position on the contig; sites must ascend.
    pub position: u64,
    /// How many alleles the record writes.
    pub allele_count: usize,
    /// The sample's likelihood of each genotype, in the genotype table's order, on any common
    /// scale — read off `PL` (§3.2).
    pub likelihoods: &'a [f64],
    /// Each template's allele shares, template by template: `templates × allele_count` numbers,
    /// each template's summing to one — read off its `GP` (§3.2). The genotype prediction and the
    /// drawn haplotypes read them.
    pub template_shares: &'a [f64],
    /// **The shares the states' emissions read, where they differ from `template_shares`** —
    /// definite alleles drawn from them (phasing spec §2.3), laid out the same way. `None`: the
    /// emissions read `template_shares`.
    pub emission_shares: Option<&'a [f64]>,
}

/// **A run of numbers at each site, every site's in one array**: site `i`'s at
/// `values[starts[i]..starts[i + 1]]`. A walk returns two of these; as a vector a site they were
/// two allocations a site a walk, the cost the drawn paths' one array removed (2026-11-13).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct PerSite {
    values: Vec<f64>,
    starts: Vec<usize>,
}

impl PerSite {
    /// Zeros, `lengths` numbers at each site in turn.
    fn zeroed(lengths: impl Iterator<Item = usize>) -> Self {
        let mut starts = vec![0];
        for length in lengths {
            let last = starts.last().copied().unwrap_or(0);
            starts.push(last + length);
        }
        Self {
            values: vec![0.0; starts.last().copied().unwrap_or(0)],
            starts,
        }
    }

    /// How many sites it holds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.starts.len().saturating_sub(1)
    }

    /// Whether it holds no site.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Each site's numbers, in order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &[f64]> + '_ {
        self.starts
            .windows(2)
            .map(|bounds| &self.values[bounds[0]..bounds[1]])
    }

    fn site_mut(&mut self, site: usize) -> &mut [f64] {
        &mut self.values[self.starts[site]..self.starts[site + 1]]
    }
}

impl std::ops::Index<usize> for PerSite {
    type Output = [f64];

    fn index(&self, site: usize) -> &[f64] {
        &self.values[self.starts[site]..self.starts[site + 1]]
    }
}

impl<'a> IntoIterator for &'a PerSite {
    type Item = &'a [f64];
    type IntoIter = Box<dyn ExactSizeIterator<Item = &'a [f64]> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter())
    }
}

/// **What a walk gives back**, one entry per site.
#[derive(Clone, Debug, PartialEq)]
pub struct WalkResult {
    /// The sample's genotype posterior at each site, in the genotype table's order.
    pub posteriors: PerSite,
    /// **What the sample's neighbours alone say at each site**: the genotype distribution the
    /// states give before the sample's own likelihood is multiplied in (§3.6), normalised. At a
    /// site where the sample has no reads it is the posterior; where it has reads, comparing the
    /// two measures how well the neighbours predict this sample.
    pub neighbour_predictions: PerSite,
    /// The probability the sample is autozygous at each site.
    pub autozygous: Vec<f64>,
    /// **The probability that a copy copies the last template** at each site — at least one of the
    /// two, or the one an autozygous state copies. The window solve puts the background template
    /// last (§3.3), so this is the chance the sample has no sampled haplotype to copy there.
    pub last_template_copy: Vec<f64>,
    /// The log-likelihood of the sample's evidence over the window under these parameters — what
    /// the parameter grid is chosen by (§3.5).
    pub log_likelihood: f64,
}

/// The two copies' alleles behind each genotype index, read off the shared genotype table — never
/// an index formula of our own (spec §7) — and kept per allele count.
#[derive(Default)]
pub(crate) struct GenotypePairs {
    by_allele_count: HashMap<usize, Vec<(usize, usize)>>,
}

impl GenotypePairs {
    pub(crate) fn of(&mut self, allele_count: usize) -> &[(usize, usize)] {
        self.by_allele_count.entry(allele_count).or_insert_with(|| {
            let diploid = Ploidy::try_new(2).expect("two copies is a ploidy");
            let table = GenotypeTable::build(diploid, allele_count);
            let table = table.view();
            (0..table.genotype_count())
                .map(|genotype| {
                    // PANIC-FREE: the index ranges over the table's own count.
                    let copies = table
                        .allele_counts_of(GenotypeIdx(genotype as u32))
                        .expect("an index below the count names a genotype");
                    let mut alleles = copies
                        .iter()
                        .enumerate()
                        .flat_map(|(allele, &count)| std::iter::repeat_n(allele, count as usize));
                    let first = alleles.next().expect("a diploid genotype has a first copy");
                    let second = alleles
                        .next()
                        .expect("a diploid genotype has a second copy");
                    (first, second)
                })
                .collect()
        })
    }

    /// The pairs of an allele count [`Self::of`] has already been asked for; `None` otherwise.
    pub(crate) fn known(&self, allele_count: usize) -> Option<&[(usize, usize)]> {
        self.by_allele_count.get(&allele_count).map(Vec::as_slice)
    }
}

/// **One site's blurred shares and likelihoods**, kept from the forward walk for the backward one
/// and the draw; its emissions are in the window's [`Emissions`].
struct SiteTerms {
    /// `q̃`: `templates × allele_count`, each template's shares blurred by the mismatch rate — the
    /// shares the prediction and the drawn haplotypes read, whatever the emissions read.
    blurred: Vec<f64>,
    /// The sample's likelihood as a symmetric `A × A` matrix, `L(a/b)` at `[a][b]` and `[b][a]`;
    /// every path's haplotypes read it.
    by_alleles: Vec<f64>,
}

/// Scratch the emissions are computed in: a template's shares weighed by the likelihood, the shares
/// by allele, and the site's own emission shares blurred where it gives them.
#[derive(Default)]
struct EmissionScratch {
    weighted: Vec<f64>,
    by_allele: Vec<f64>,
    own: Vec<f64>,
}

/// One site's terms, its emissions written over `autozygous` (`K`) and `pair_emissions` (`K × K`,
/// row `k₁`).
fn site_terms(
    site: &WalkSite<'_>,
    templates: usize,
    mismatch: f64,
    pairs: &[(usize, usize)],
    emissions: (&mut [f64], &mut [f64]),
    scratch: &mut EmissionScratch,
) -> SiteTerms {
    let alleles = site.allele_count;
    assert_eq!(
        site.likelihoods.len(),
        pairs.len(),
        "a site at {} has {} alleles, so {} diploid genotypes, and {} likelihoods",
        site.position,
        alleles,
        pairs.len(),
        site.likelihoods.len()
    );
    assert_eq!(
        site.template_shares.len(),
        templates * alleles,
        "a site at {} needs {templates} templates × {alleles} allele shares",
        site.position
    );
    let blurred: Vec<f64> = site
        .template_shares
        .iter()
        .map(|&share| blur(share, mismatch, alleles))
        .collect();
    // The likelihood as a symmetric allele-by-allele matrix, so the pair emission is two matrix
    // products: Σ_a Σ_b q̃_k₁(a) · L(a/b) · q̃_k₂(b).
    let mut by_alleles = vec![0.0; alleles * alleles];
    for (&(first, second), &likelihood) in pairs.iter().zip(site.likelihoods) {
        by_alleles[first * alleles + second] = likelihood;
        by_alleles[second * alleles + first] = likelihood;
    }
    let terms = SiteTerms {
        blurred,
        by_alleles,
    };
    emissions_into(site, &terms, templates, mismatch, emissions, scratch);
    terms
}

/// A template's share blurred by the mismatch rate: `(1 − ε)·q + ε/A`.
fn blur(share: f64, mismatch: f64, alleles: usize) -> f64 {
    (1.0 - mismatch) * share + mismatch / alleles as f64
}

/// **One site's emissions from its terms**, written over `autozygous` (`K`) and `pair_emissions`
/// (`K × K`, row `k₁`): what [`site_terms`] writes, computed again where a walk needs a site's
/// emissions it no longer holds.
fn emissions_into(
    site: &WalkSite<'_>,
    terms: &SiteTerms,
    templates: usize,
    mismatch: f64,
    (autozygous, pair_emissions): (&mut [f64], &mut [f64]),
    scratch: &mut EmissionScratch,
) {
    let alleles = site.allele_count;
    let EmissionScratch {
        weighted,
        by_allele,
        own,
    } = scratch;
    // The emissions read their own shares where the site gives them, the blurred shares otherwise.
    let emitting: &[f64] = match site.emission_shares {
        Some(shares) => {
            assert_eq!(
                shares.len(),
                site.template_shares.len(),
                "a site at {} gives emission shares laid out unlike its template shares",
                site.position
            );
            own.clear();
            own.extend(shares.iter().map(|&share| blur(share, mismatch, alleles)));
            own
        }
        None => &terms.blurred,
    };
    let by_alleles = &terms.by_alleles;
    weighted.clear();
    weighted.resize(templates * alleles, 0.0);
    for (out, shares) in weighted
        .chunks_exact_mut(alleles)
        .zip(emitting.chunks_exact(alleles))
    {
        for (b, value) in out.iter_mut().enumerate() {
            *value = shares
                .iter()
                .zip(by_alleles.iter().skip(b).step_by(alleles))
                .map(|(share, likelihood)| share * likelihood)
                .sum();
        }
    }
    // Each first template's row of emissions is built one allele at a time across every second
    // template, so the inner loop runs over templates and vectorises; each emission still adds its
    // alleles' terms in order, so the bits are those of a sum per pair.
    by_allele.clear();
    by_allele.resize(alleles * templates, 0.0);
    for (template, shares) in emitting.chunks_exact(alleles).enumerate() {
        for (allele, &share) in shares.iter().enumerate() {
            by_allele[allele * templates + template] = share;
        }
    }
    for (out, first) in pair_emissions
        .chunks_exact_mut(templates)
        .zip(weighted.chunks_exact(alleles))
    {
        for (allele, (&weight, second)) in first
            .iter()
            .zip(by_allele.chunks_exact(templates))
            .enumerate()
        {
            for (value, &share) in out.iter_mut().zip(second) {
                if allele == 0 {
                    *value = weight * share;
                } else {
                    *value += weight * share;
                }
            }
        }
    }
    for (k, emission) in autozygous.iter_mut().enumerate() {
        *emission = (0..alleles)
            .map(|a| emitting[k * alleles + a] * by_alleles[a * alleles + a])
            .sum();
    }
}

/// The three persistence probabilities between two sites `gap` bases apart: a copy keeps its
/// template, the autozygous state persists, the other state persists.
///
/// **The two states are one two-state Markov chain** (spec §3.5): runs of mean `R`, gaps of mean
/// `R′ = R · (1 − F) / F`, so the chance of staying is `F + (1 − F) · e^(−d/τ)` in a run and
/// `(1 − F) + F · e^(−d/τ)` outside one, with `1/τ = 1/R + 1/R′`. Across a gap far longer than
/// both the state returns to its long-run share `F` — not, as two independent `e^(−d/R)` stays
/// would make it, to the opposite state (found by review, 2026-11-01).
fn persistence(gap: f64, parameters: &WalkParameters) -> (f64, f64, f64) {
    let share = parameters.autozygous_share;
    let gap_bases = parameters.run_bases * (1.0 - share) / share;
    let decay = float::exp(-gap * (1.0 / parameters.run_bases + 1.0 / gap_bases));
    (
        float::exp(-gap / parameters.segment_bases),
        share + (1.0 - share) * decay,
        (1.0 - share) + share * decay,
    )
}

/// **One sample's walk over a window's sites** (spec §3). `templates` is `K`; every site carries
/// that many templates' shares. `template_weights` weighs where a copy that changes template lands,
/// one relative weight for each template but the background (`WalkParameters::landing_weights`);
/// `None` weighs them alike. Returns each site's genotype posterior and autozygous probability.
///
/// # Panics
///
/// If a site's likelihoods or shares do not match its allele count, the autozygous share is not
/// strictly between 0 and 1, or `template_weights` is not `templates − 1` long.
#[must_use]
pub fn walk(
    sites: &[WalkSite<'_>],
    templates: usize,
    template_weights: Option<&[f64]>,
    parameters: WalkParameters,
) -> WalkResult {
    walk_with_draw(
        sites,
        templates,
        template_weights,
        parameters,
        None,
        None,
        None,
    )
    .0
}

/// **A sample's read link between two walk sites**: `L` between `one` and `other`, each with its
/// heterozygote's two alleles, the lower first; `None` when no fragment of the sample counts at both.
pub(crate) type ReadLinkOf<'a> =
    dyn Fn(usize, (usize, usize), usize, (usize, usize)) -> Option<f64> + 'a;

/// **What the reads give the haplotype draw** (`neighbour_evidence_reads_in_draw.md` §2.2): one
/// sample's read link between two walk sites, and the stream the orientations are drawn from.
pub(crate) struct DrawReads<'a> {
    /// The sample's read links.
    pub link: &'a ReadLinkOf<'a>,
    /// Whether the sample's reads at a walk site are kept out of the draw: one of its hidden
    /// records (§5) or of its records held out for the read check (§2.6). Such a site never joins a
    /// chain.
    pub kept_out: &'a dyn Fn(usize) -> bool,
    /// `HET_MASS_MIN`: the smallest heterozygous mass of an orientable record.
    pub het_mass_min: f64,
    /// The orientations' stream, read only where a read-linked group exists.
    pub uniforms: Uniforms,
}

/// **One haplotype pair at one site**: each copy's probability of carrying each allele, `A`
/// numbers a copy.
#[cfg(test)]
pub(crate) type HaplotypePair = [Vec<f64>; 2];

/// **One drawn path's two haplotypes at every site of a walk** (phasing spec §2.1): at each site,
/// copy one's `A` allele probabilities and then copy two's, every site's in one array. One array
/// rather than two vectors a site: a round draws nine paths at every site of a window, and the
/// vectors' allocations were a tenth of the pass's time and a fifth of a walk's memory.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct PathHaplotypes {
    /// Where each site's numbers start in `values`.
    starts: Vec<usize>,
    values: Vec<f64>,
}

impl PathHaplotypes {
    fn with_capacity(sites: usize) -> Self {
        Self {
            starts: Vec::with_capacity(sites),
            values: Vec::with_capacity(4 * sites),
        }
    }

    /// Appends a site of `alleles` alleles, its two copies at zero, to be written.
    fn push_site(&mut self, alleles: usize) -> (&mut [f64], &mut [f64]) {
        let start = self.values.len();
        self.starts.push(start);
        self.values.resize(start + 2 * alleles, 0.0);
        self.values[start..].split_at_mut(alleles)
    }

    /// How many sites the path holds.
    pub(crate) fn len(&self) -> usize {
        self.starts.len()
    }

    /// The two copies' allele probabilities at `site`.
    pub(crate) fn site(&self, site: usize) -> [&[f64]; 2] {
        let start = self.starts[site];
        let end = self
            .starts
            .get(site + 1)
            .copied()
            .unwrap_or(self.values.len());
        let (first, second) = self.values[start..end].split_at((end - start) / 2);
        [first, second]
    }

    /// Every site's two copies, in order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = [&[f64]; 2]> + '_ {
        (0..self.len()).map(|site| self.site(site))
    }
}

#[cfg(test)]
impl From<Vec<HaplotypePair>> for PathHaplotypes {
    fn from(pairs: Vec<HaplotypePair>) -> Self {
        let mut path = Self::with_capacity(pairs.len());
        for [first, second] in pairs {
            let (one, two) = path.push_site(first.len());
            one.copy_from_slice(&first);
            two.copy_from_slice(&second);
        }
        path
    }
}

/// **Every site's joint of the two copies' alleles along one path**, `A × A` numbers a site in one
/// array, and whether the site has one (`None` in [`joint_along`]'s terms). Kept across a walk's
/// paths and refilled for each.
#[derive(Default)]
struct Joints {
    starts: Vec<usize>,
    values: Vec<f64>,
    present: Vec<bool>,
}

impl Joints {
    fn clear(&mut self) {
        self.starts.clear();
        self.values.clear();
        self.present.clear();
    }

    fn len(&self) -> usize {
        self.starts.len()
    }

    fn range(&self, site: usize) -> std::ops::Range<usize> {
        let start = self.starts[site];
        start
            ..self
                .starts
                .get(site + 1)
                .copied()
                .unwrap_or(self.values.len())
    }

    fn get(&self, site: usize) -> Option<&[f64]> {
        self.present[site].then(|| &self.values[self.range(site)])
    }

    fn get_mut(&mut self, site: usize) -> Option<&mut [f64]> {
        let range = self.range(site);
        self.present[site].then(move || &mut self.values[range])
    }
    #[cfg(test)]
    fn from_options(joints: &[Option<Vec<f64>>]) -> Self {
        let mut out = Self::default();
        for joint in joints {
            out.starts.push(out.values.len());
            out.present.push(joint.is_some());
            match joint {
                Some(values) => out.values.extend_from_slice(values),
                // A missing joint still takes its site's room; 2 × 2 serves the tests' sites.
                None => out.values.extend_from_slice(&[0.0; 4]),
            }
        }
        out
    }

    #[cfg(test)]
    fn to_options(&self) -> Vec<Option<Vec<f64>>> {
        (0..self.len())
            .map(|site| self.get(site).map(<[f64]>::to_vec))
            .collect()
    }
}

/// **A round's own drawn path**: the two haplotypes along it, and at each site the walk template
/// each copy follows there — copy one's, then copy two's; the same template twice in an
/// autozygous state. The last template is the background. The next round reads which templates
/// were followed to keep a sample from copying its own reads back (`neighbour_evidence_phasing.md`
/// §2.3).
#[derive(Debug, PartialEq)]
pub(crate) struct DrawnPath {
    pub haplotypes: PathHaplotypes,
    pub followed: Vec<[usize; 2]>,
}

/// **The walk, and one state path drawn from it with the two haplotypes along it**
/// (`neighbour_evidence_phasing.md` §2.1). The walk is [`walk`]'s, unchanged. The path is drawn
/// backwards from the last site, each state given the one after it and the forward probabilities,
/// one uniform a site from `uniforms`, the states in the order autozygous `0..K` then pairs row by
/// row. Along the path, copy one carries allele `a` and copy two `b` with probability
/// `∝ q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})` in a pair state and `∝ q̃_k(a) · L({a, a})`, `a = b`, in
/// an autozygous one; each copy's haplotype is that joint's margin.
///
/// # Panics
///
/// As [`walk`]; and if every state at a site has zero weight — every genotype's likelihood zero
/// there, which a likelihood row normalised to its largest value never is.
#[must_use]
pub(crate) fn walk_and_draw(
    sites: &[WalkSite<'_>],
    templates: usize,
    template_weights: Option<&[f64]>,
    parameters: WalkParameters,
    uniforms: &mut Uniforms,
    reads: Option<&mut DrawReads<'_>>,
) -> (WalkResult, DrawnPath) {
    let (result, drawn, _) = walk_with_draw(
        sites,
        templates,
        template_weights,
        parameters,
        Some(uniforms),
        None,
        reads,
    );
    (
        result,
        // PANIC-FREE: `walk_with_draw` draws whenever it is given a stream, as it is here.
        drawn.expect("a stream was given, so a path was drawn"),
    )
}

/// **[`walk_and_draw`], and `extra` more paths drawn from the same forward probabilities**
/// (`neighbour_evidence_phase_output.md` §2.3): the round's own path from `uniforms`, exactly as
/// [`walk_and_draw`] draws it, and the extra paths one after another from `extra_uniforms`, so
/// drawing them changes neither the round's path nor anything the round passes on. The extra paths
/// serve the phase's link test only.
///
/// # Panics
///
/// As [`walk_and_draw`].
#[must_use]
#[allow(clippy::too_many_arguments)] // the walk's own inputs and the two draws' streams
pub(crate) fn walk_and_draw_paths(
    sites: &[WalkSite<'_>],
    templates: usize,
    template_weights: Option<&[f64]>,
    parameters: WalkParameters,
    uniforms: &mut Uniforms,
    extra_uniforms: &mut Uniforms,
    extra: usize,
    reads: Option<&mut DrawReads<'_>>,
) -> (WalkResult, DrawnPath, Vec<PathHaplotypes>) {
    let (result, drawn, extra_paths) = walk_with_draw(
        sites,
        templates,
        template_weights,
        parameters,
        Some(uniforms),
        Some((extra_uniforms, extra)),
        reads,
    );
    (
        result,
        // PANIC-FREE: as in `walk_and_draw`.
        drawn.expect("a stream was given, so a path was drawn"),
        extra_paths,
    )
}

/// **The forward walk's probabilities, a few sites' in one array of each kind** rather than a
/// vector a site: slot `i`'s `K` autozygous states at `autozygous[i·K..]`, its `K²` pairs row by
/// row at `pairs[i·K²..]`, and those pairs' row and column sums at `rows[i·K..]` and
/// `columns[i·K..]`. **A slot is not a window's site**: the forward pass uses two slots, the site and
/// the one before, and the backward walk one slot a site of the block it recomputes
/// ([`FORWARD_BLOCK`]), slot 0 at the block's first site.
///
/// **Every number is written before it is read**, so a walk takes the arrays as the last walk on
/// its thread left them ([`WALK_ARRAYS`]) rather than allocating them afresh each walk.
#[derive(Default)]
struct Forward {
    k: usize,
    autozygous: Vec<f64>,
    pairs: Vec<f64>,
    rows: Vec<f64>,
    columns: Vec<f64>,
}

/// **Some sites' emissions, in one array of each kind** as [`Forward`] holds its probabilities:
/// slot `i`'s `K` autozygous emissions at `autozygous[i·K..]` and its `K²` pair emissions, row
/// `k₁`, at `pairs[i·K²..]`. A walk holds a block's ([`FORWARD_BLOCK`]) and one more: the window's
/// every site, at `K + K²` numbers each, was a third of the pass's memory on 400 tomato samples and
/// 48 threads (pass spec §6).
#[derive(Default)]
struct Emissions {
    k: usize,
    autozygous: Vec<f64>,
    pairs: Vec<f64>,
}

impl Emissions {
    /// Room for `sites` sites of `k` templates, holding whatever the last walk left.
    fn resize(&mut self, sites: usize, k: usize) {
        self.k = k;
        self.autozygous.resize(sites * k, 0.0);
        self.pairs.resize(sites * k * k, 0.0);
    }

    /// Site `site`'s two kinds of emission, to be written.
    fn site_mut(&mut self, site: usize) -> (&mut [f64], &mut [f64]) {
        let (k, states) = (self.k, self.k * self.k);
        (
            &mut self.autozygous[site * k..(site + 1) * k],
            &mut self.pairs[site * states..(site + 1) * states],
        )
    }

    fn autozygous(&self, site: usize) -> &[f64] {
        &self.autozygous[site * self.k..(site + 1) * self.k]
    }

    fn pairs(&self, site: usize) -> &[f64] {
        let states = self.k * self.k;
        &self.pairs[site * states..(site + 1) * states]
    }

    /// Copies slot `from`'s emissions into slot `to`.
    fn carry(&mut self, from: usize, to: usize) {
        let (k, states) = (self.k, self.k * self.k);
        self.autozygous
            .copy_within(from * k..(from + 1) * k, to * k);
        self.pairs
            .copy_within(from * states..(from + 1) * states, to * states);
    }
}

/// **How many sites a block of the backward walk holds** (see [`walk_with_draw`]): its forward
/// probabilities and emissions, about 330 KB at 17 templates, within a core's share of the L2
/// cache on the machine it was measured on (2 MB a core, two threads).
const FORWARD_BLOCK: usize = 64;

/// A walk's arrays: a block's emissions and forward probabilities, and the forward probabilities
/// at every block's first site.
#[derive(Default)]
struct WalkArrays {
    emissions: Emissions,
    forward: Forward,
    checkpoints: Checkpoints,
}

/// **The forward probabilities at the first site of every block**, `K + K²` numbers a block.
#[derive(Default)]
struct Checkpoints {
    k: usize,
    autozygous: Vec<f64>,
    pairs: Vec<f64>,
}

impl Checkpoints {
    fn resize(&mut self, blocks: usize, k: usize) {
        self.k = k;
        self.autozygous.resize(blocks * k, 0.0);
        self.pairs.resize(blocks * k * k, 0.0);
    }

    /// Keeps `forward`'s slot `slot` as block `block`'s.
    fn save(&mut self, block: usize, forward: &Forward, slot: usize) {
        let (k, states) = (self.k, self.k * self.k);
        self.autozygous[block * k..(block + 1) * k].copy_from_slice(forward.autozygous(slot));
        self.pairs[block * states..(block + 1) * states].copy_from_slice(forward.pairs(slot));
    }

    /// Writes block `block`'s into `forward`'s slot `slot`.
    fn restore(&self, block: usize, forward: &mut Forward, slot: usize) {
        let (k, states) = (self.k, self.k * self.k);
        forward.autozygous[slot * k..(slot + 1) * k]
            .copy_from_slice(&self.autozygous[block * k..(block + 1) * k]);
        forward.pairs[slot * states..(slot + 1) * states]
            .copy_from_slice(&self.pairs[block * states..(block + 1) * states]);
    }
}

/// **`values`' chunk `read` to read and chunk `write` to write**, chunks of `size`; two different
/// chunks.
fn read_and_write(
    values: &mut [f64],
    size: usize,
    read: usize,
    write: usize,
) -> (&[f64], &mut [f64]) {
    assert_ne!(read, write, "a step reads one slot and writes another");
    if read < write {
        let (before, after) = values.split_at_mut(write * size);
        (&before[read * size..(read + 1) * size], &mut after[..size])
    } else {
        let (before, after) = values.split_at_mut(read * size);
        (
            &after[..size],
            &mut before[write * size..(write + 1) * size],
        )
    }
}

/// **One forward step**, written into `forward`'s slot `at`: from the landing weights at a walk's
/// first site, or from slot `previous` across the persistence `(keep, stay autozygous, stay
/// other)` between the two sites — writing the previous slot's row and column sums on the way —
/// then times the site's `emissions` and normalised. Returns the scale it was divided by. The
/// forward walk and the backward walk's blocks both step through this one function, so a block
/// computed again holds the same bits.
fn forward_step(
    forward: &mut Forward,
    at: usize,
    previous: Option<(usize, (f64, f64, f64))>,
    landing: &[f64],
    share: f64,
    (emission_autozygous, emission_pairs): (&[f64], &[f64]),
) -> f64 {
    let k = forward.k;
    match previous {
        None => {
            let autozygous = &mut forward.autozygous[at * k..(at + 1) * k];
            for (value, &weight) in autozygous.iter_mut().zip(landing) {
                *value = share * weight;
            }
            let pairs = &mut forward.pairs[at * k * k..(at + 1) * k * k];
            for (row, &first) in pairs.chunks_exact_mut(k).zip(landing) {
                for (value, &second) in row.iter_mut().zip(landing) {
                    *value = (1.0 - share) * first * second;
                }
            }
        }
        Some((before, (keep, stay_autozygous, stay_other))) => {
            let (previous_autozygous, autozygous) =
                read_and_write(&mut forward.autozygous, k, before, at);
            let (previous_pairs, pairs) = read_and_write(&mut forward.pairs, k * k, before, at);
            let mass_autozygous: f64 = previous_autozygous.iter().sum();
            let mass_pairs = lane_sum(previous_pairs);
            // The previous site's margins, which this step needs and the draw reads again.
            let rows = &mut forward.rows[before * k..(before + 1) * k];
            let columns = &mut forward.columns[before * k..(before + 1) * k];
            margins(previous_pairs, k, rows, columns);
            for ((value, &previous), &weight) in
                autozygous.iter_mut().zip(previous_autozygous).zip(landing)
            {
                let within = keep * previous + (1.0 - keep) * mass_autozygous * weight;
                *value = stay_autozygous * within + (1.0 - stay_other) * mass_pairs * weight;
            }
            forward_pair_step(
                previous_pairs,
                rows,
                columns,
                landing,
                PairTransition::new(keep, mass_pairs, stay_other),
                (1.0 - stay_autozygous) * mass_autozygous,
                pairs,
            );
        }
    }
    let autozygous = &mut forward.autozygous[at * k..(at + 1) * k];
    for (value, emission) in autozygous.iter_mut().zip(emission_autozygous) {
        *value *= emission;
    }
    let pairs = &mut forward.pairs[at * k * k..(at + 1) * k * k];
    for (value, emission) in pairs.iter_mut().zip(emission_pairs) {
        *value *= emission;
    }
    let scale: f64 = autozygous.iter().sum::<f64>() + lane_sum(pairs);
    let inverse = 1.0 / scale;
    scale_by(autozygous, inverse);
    scale_by(pairs, inverse);
    scale
}

thread_local! {
    /// The arrays a walk on this thread fills — see [`Forward`].
    static WALK_ARRAYS: std::cell::Cell<WalkArrays> = std::cell::Cell::new(WalkArrays::default());
}

impl Forward {
    /// Room for `sites` sites of `k` templates, holding whatever the last walk left.
    fn resize(&mut self, sites: usize, k: usize) {
        self.k = k;
        self.autozygous.resize(sites * k, 0.0);
        self.pairs.resize(sites * k * k, 0.0);
        self.rows.resize(sites * k, 0.0);
        self.columns.resize(sites * k, 0.0);
    }

    fn autozygous(&self, site: usize) -> &[f64] {
        &self.autozygous[site * self.k..(site + 1) * self.k]
    }

    fn pairs(&self, site: usize) -> &[f64] {
        let states = self.k * self.k;
        &self.pairs[site * states..(site + 1) * states]
    }

    fn rows(&self, site: usize) -> &[f64] {
        &self.rows[site * self.k..(site + 1) * self.k]
    }

    fn columns(&self, site: usize) -> &[f64] {
        &self.columns[site * self.k..(site + 1) * self.k]
    }
}

/// A state of the walk: autozygous on one template, or a pair of templates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum State {
    Autozygous(usize),
    Pair(usize, usize),
}

fn walk_with_draw(
    sites: &[WalkSite<'_>],
    templates: usize,
    template_weights: Option<&[f64]>,
    parameters: WalkParameters,
    uniforms: Option<&mut Uniforms>,
    extra: Option<(&mut Uniforms, usize)>,
    mut reads: Option<&mut DrawReads<'_>>,
) -> (WalkResult, Option<DrawnPath>, Vec<PathHaplotypes>) {
    assert!(
        sites
            .windows(2)
            .all(|pair| pair[0].position <= pair[1].position),
        "a walk's sites must not descend in position"
    );
    let share = parameters.autozygous_share;
    assert!(
        share > 0.0 && share < 1.0 && templates > 0,
        "a walk needs an autozygous share strictly between 0 and 1 and at least one template: \
         got {share} and {templates}"
    );
    let k = templates;
    let landing = parameters.landing_weights(k, template_weights);
    let mut genotype_pairs = GenotypePairs::default();
    let WalkArrays {
        mut emissions,
        mut forward,
        mut checkpoints,
    } = WALK_ARRAYS.take();
    // **The forward probabilities are kept at the first site of every block of
    // [`FORWARD_BLOCK`] sites, not at every site.** The backward walk computes a block's again from
    // its first site, by the same steps from the same inputs, just before it reads them: a
    // block's probabilities and emissions fit a core's cache, where the window's did not, and the
    // walk spent most of its time waiting for memory (perf: 0.6 instructions a cycle).
    let blocks = sites.len().div_ceil(FORWARD_BLOCK);
    forward.resize(FORWARD_BLOCK, k);
    emissions.resize(FORWARD_BLOCK + 1, k);
    checkpoints.resize(blocks, k);
    let mut scratch = EmissionScratch::default();
    let mut terms: Vec<SiteTerms> = Vec::with_capacity(sites.len());
    // Each site's persistence from the site before, entry 0 unused: read by the forward step, the
    // block's recomputation, the draw and the backward step, so computed once rather than four
    // times — three `exp`s each.
    let persistences: Vec<(f64, f64, f64)> = std::iter::once((1.0, 1.0, 1.0))
        .chain(
            sites
                .windows(2)
                .map(|pair| persistence((pair[1].position - pair[0].position) as f64, &parameters)),
        )
        .collect();
    let transition = |index: usize| persistences[index];

    // Forward: each site's state distribution given the sites up to it, normalised; the scales
    // multiply to the window's likelihood. Two slots, the site's and the one before.
    let mut log_likelihood = 0.0;
    for (index, site) in sites.iter().enumerate() {
        terms.push(site_terms(
            site,
            k,
            parameters.mismatch,
            genotype_pairs.of(site.allele_count),
            emissions.site_mut(0),
            &mut scratch,
        ));
        let previous = (index > 0).then(|| ((index - 1) % 2, transition(index)));
        let scale = forward_step(
            &mut forward,
            index % 2,
            previous,
            &landing,
            share,
            (emissions.autozygous(0), emissions.pairs(0)),
        );
        log_likelihood += float::ln(scale);
        if index % FORWARD_BLOCK == 0 {
            checkpoints.save(index / FORWARD_BLOCK, &forward, index % 2);
        }
    }

    // The round's own path from `uniforms`, then each extra path from `extra_uniforms`, one
    // uniform a site from the last site back. All are drawn in the backward walk's one sweep, so
    // each site's totals are summed once for every path.
    let mut streams: Vec<Vec<f64>> = Vec::new();
    let draws_own = uniforms.is_some();
    if let Some(uniforms) = uniforms {
        streams.push((0..sites.len()).map(|_| uniforms.next_uniform()).collect());
    }
    if let Some((extra_uniforms, count)) = extra {
        for _ in 0..count {
            streams.push(
                (0..sites.len())
                    .map(|_| extra_uniforms.next_uniform())
                    .collect(),
            );
        }
    }
    let mut state_paths: Vec<Vec<State>> = vec![Vec::with_capacity(sites.len()); streams.len()];

    // Backward, right to left, and each site's posterior as soon as its backward vector exists.
    let mut posteriors = PerSite::zeroed(
        sites
            .iter()
            .map(|site| genotype_pairs.of(site.allele_count).len()),
    );
    let mut neighbour_predictions = posteriors.clone();
    let mut autozygous_probability = vec![0.0; sites.len()];
    let mut last_template_copy = vec![0.0; sites.len()];
    let mut backward_autozygous = vec![1.0; k];
    let mut backward_pairs = vec![1.0; k * k];
    // Every site's scratch, reused so the backward walk allocates only what it returns.
    let mut weighted_autozygous = vec![0.0; k];
    let mut weighted_pairs = vec![0.0; k * k];
    let (mut rows, mut columns) = (vec![0.0; k], vec![0.0; k]);
    let mut state_autozygous = vec![0.0; k];
    let mut by_second = vec![0.0; k * k];
    let (mut second_copy, mut pair_weight, mut autozygous_weight) =
        (Vec::new(), Vec::new(), Vec::new());
    let mut pairs_of_alleles = genotype_pairs;
    for block in (0..blocks).rev() {
        // The block's forward probabilities from its first site, and its emissions; slot
        // `FORWARD_BLOCK` of the emissions holds the next block's first site's.
        let start = block * FORWARD_BLOCK;
        let end = (start + FORWARD_BLOCK).min(sites.len());
        checkpoints.restore(block, &mut forward, 0);
        for index in start..end {
            let local = index - start;
            emissions_into(
                &sites[index],
                &terms[index],
                k,
                parameters.mismatch,
                emissions.site_mut(local),
                &mut scratch,
            );
            if index > start {
                forward_step(
                    &mut forward,
                    local,
                    Some((local - 1, transition(index))),
                    &landing,
                    share,
                    (emissions.autozygous(local), emissions.pairs(local)),
                );
            }
        }
        let last = end - start - 1;
        margins(
            &forward.pairs[last * k * k..(last + 1) * k * k],
            k,
            &mut forward.rows[last * k..(last + 1) * k],
            &mut forward.columns[last * k..(last + 1) * k],
        );
        for index in (start..end).rev() {
            let local = index - start;
            let next_slot = if index + 1 < end {
                local + 1
            } else {
                FORWARD_BLOCK
            };
            if !streams.is_empty() {
                let draw = SiteDraw::new(
                    forward.autozygous(local),
                    forward.pairs(local),
                    forward.rows(local),
                    forward.columns(local),
                    k,
                );
                let persistence = (index + 1 < sites.len()).then(|| transition(index + 1));
                for (path, stream) in state_paths.iter_mut().zip(&streams) {
                    let next = persistence.zip(path.last().copied());
                    path.push(draw.draw(next, stream[sites.len() - 1 - index], &landing));
                }
            }
            if index + 1 < sites.len() {
                let (keep, stay_autozygous, stay_other) = transition(index + 1);
                for ((weighted, value), emission) in weighted_autozygous
                    .iter_mut()
                    .zip(&backward_autozygous)
                    .zip(emissions.autozygous(next_slot))
                {
                    *weighted = value * emission;
                }
                for ((weighted, value), emission) in weighted_pairs
                    .iter_mut()
                    .zip(&backward_pairs)
                    .zip(emissions.pairs(next_slot))
                {
                    *weighted = value * emission;
                }
                // Read backwards, a jump lands where `landing` says, so what a state can reach is
                // weighed by it: each row and column summed with the landing weights, and the totals.
                let total_autozygous: f64 = weighted_autozygous
                    .iter()
                    .zip(&landing)
                    .map(|(value, weight)| value * weight)
                    .sum();
                landing_margins(&weighted_pairs, &landing, &mut rows, &mut columns);
                let total_pairs: f64 = rows
                    .iter()
                    .zip(&landing)
                    .map(|(row, weight)| row * weight)
                    .sum();
                for (value, weighted) in backward_autozygous.iter_mut().zip(&weighted_autozygous) {
                    let within = keep * weighted + (1.0 - keep) * total_autozygous;
                    *value = stay_autozygous * within + (1.0 - stay_autozygous) * total_pairs;
                }
                backward_pair_step(
                    &weighted_pairs,
                    &rows,
                    &columns,
                    PairTransition::new(keep, total_pairs, stay_other),
                    (1.0 - stay_other) * total_autozygous,
                    &mut backward_pairs,
                );
                let total: f64 =
                    backward_autozygous.iter().sum::<f64>() + lane_sum(&backward_pairs);
                let inverse = 1.0 / total;
                scale_by(&mut backward_autozygous, inverse);
                scale_by(&mut backward_pairs, inverse);
            }
            let site = &sites[index];
            let site_terms = &terms[index];
            let alleles = site.allele_count;

            // The state posterior with this site's own emission divided back out, so the sample's
            // likelihood enters once, below (§3.6).
            //
            // **A state whose emission is exactly zero contributes zero**, not 0/0. It happens where the
            // reads rule out every genotype the state allows — a deep heterozygote's homozygous
            // likelihoods underflow to 0 at 30× — and then each of the emission's non-negative terms,
            // q̃ · q̃ · L(g), is zero, so the state adds nothing to any genotype whatever its weight.
            let without_emission = |both: f64, emission: f64| {
                if emission > 0.0 { both / emission } else { 0.0 }
            };
            let (forward_autozygous, forward_pairs) =
                (forward.autozygous(local), forward.pairs(local));
            for (state, value) in state_autozygous.iter_mut().enumerate() {
                *value = without_emission(
                    forward_autozygous[state] * backward_autozygous[state],
                    emissions.autozygous(local)[state],
                );
            }
            // Column by column, `[second][first]`, for the spreading over alleles below.
            // Read row by row, so the three inputs are read in order, and written to their columns.
            // Two loops rather than one over the states split by `% k` and `/ k`: those were two
            // integer divisions a state, by a `k` known only at run time, and the line was an eighth
            // of the pass's time (perf, 800 tomato samples).
            let pair_emissions = emissions.pairs(local);
            for first in 0..k {
                let row = first * k..(first + 1) * k;
                for (second, ((&forward, &backward), &emission)) in forward_pairs[row.clone()]
                    .iter()
                    .zip(&backward_pairs[row.clone()])
                    .zip(&pair_emissions[row])
                    .enumerate()
                {
                    by_second[second * k + first] = without_emission(forward * backward, emission);
                }
            }
            let both_autozygous: f64 = (0..k)
                .map(|state| forward_autozygous[state] * backward_autozygous[state])
                .sum();
            let both_pairs = lane_dot(forward_pairs, &backward_pairs);
            autozygous_probability[index] = both_autozygous / (both_autozygous + both_pairs);
            let last = k - 1;
            let with_last: f64 = forward_autozygous[last] * backward_autozygous[last]
                + (0..k)
                    .map(|other| {
                        let row =
                            forward_pairs[last * k + other] * backward_pairs[last * k + other];
                        let column = if other == last {
                            0.0
                        } else {
                            forward_pairs[other * k + last] * backward_pairs[other * k + last]
                        };
                        row + column
                    })
                    .sum::<f64>();
            last_template_copy[index] = with_last / (both_autozygous + both_pairs);

            // The pair states' weight on each ordered allele pair, Q̃ᵀ · W · Q̃ (§3.6): first each
            // first-copy template's weight spread over the second copy's alleles (K² · A), then over
            // the first copy's (K · A²) — not K² · A² by summing every state over every pair.
            let blurred = &site_terms.blurred;
            // `second_copy[b][first]`: the second templates taken in order, each adding its share of
            // allele `b` to every first template at once — the same terms in the same order as a sum
            // per first template, but the inner loop runs over templates and vectorises.
            second_copy.clear();
            second_copy.resize(alleles * k, 0.0);
            for (column, shares) in by_second.chunks_exact(k).zip(blurred.chunks_exact(alleles)) {
                for (out, &share) in second_copy.chunks_exact_mut(k).zip(shares) {
                    for (value, &weight) in out.iter_mut().zip(column) {
                        *value += weight * share;
                    }
                }
            }
            pair_weight.clear();
            pair_weight.resize(alleles * alleles, 0.0);
            for (first, first_shares) in blurred.chunks_exact(alleles).enumerate() {
                for (&first_share, out) in first_shares
                    .iter()
                    .zip(pair_weight.chunks_exact_mut(alleles))
                {
                    for (value, spread) in out.iter_mut().zip(second_copy.chunks_exact(k)) {
                        *value += first_share * spread[first];
                    }
                }
            }
            autozygous_weight.clear();
            autozygous_weight.resize(alleles, 0.0);
            for (&weight, shares) in state_autozygous.iter().zip(blurred.chunks_exact(alleles)) {
                for (value, share) in autozygous_weight.iter_mut().zip(shares) {
                    *value += weight * share;
                }
            }
            let prediction = neighbour_predictions.site_mut(index);
            for (value, &(a, b)) in prediction.iter_mut().zip(pairs_of_alleles.of(alleles)) {
                let mut weight = pair_weight[a * alleles + b];
                if a == b {
                    weight += autozygous_weight[a];
                } else {
                    weight += pair_weight[b * alleles + a];
                }
                *value = weight;
            }
            let posterior = posteriors.site_mut(index);
            for ((value, weight), likelihood) in posterior
                .iter_mut()
                .zip(prediction.iter())
                .zip(site.likelihoods)
            {
                *value = weight * likelihood;
            }
            let total: f64 = posterior.iter().sum();
            for value in posterior.iter_mut() {
                *value /= total;
            }
            let total: f64 = prediction.iter().sum();
            for value in prediction.iter_mut() {
                *value /= total;
            }
        }
        emissions.carry(0, FORWARD_BLOCK);
    }
    WALK_ARRAYS.set(WalkArrays {
        emissions,
        forward,
        checkpoints,
    });

    // The paths were drawn from the last site back.
    for path in &mut state_paths {
        path.reverse();
    }
    let mut joints = Joints::default();
    let mut paths: Vec<PathHaplotypes> = state_paths
        .iter()
        .map(|states| haplotypes_along(sites, &terms, states, reads.as_deref_mut(), &mut joints))
        .collect();
    let drawn = draws_own.then(|| DrawnPath {
        haplotypes: paths.remove(0),
        followed: state_paths[0]
            .iter()
            .map(|&state| match state {
                State::Autozygous(template) => [template, template],
                State::Pair(first, second) => [first, second],
            })
            .collect(),
    });
    let extra_paths = paths;

    (
        WalkResult {
            posteriors,
            neighbour_predictions,
            autozygous: autozygous_probability,
            last_template_copy,
            log_likelihood,
        },
        drawn,
        extra_paths,
    )
}

/// **One site's forward probabilities, arranged so a state can be drawn given the next in about
/// `K` steps** rather than a scan of all `K + K²` weights.
///
/// Given the next state, a state's weight is its forward probability times the transition, and
/// the transition is a share of every state of a kind, plus — towards a pair — a share of the next
/// pair's row and of its column, plus a share of the next state itself (spec §3.5). A uniform
/// first picks one of those parts by its mass, then a state inside it: by bisection over the
/// pairs' running totals for a share of every pair, by a scan of `K` for a row, a column or the
/// autozygous states. That draws from exactly the distribution the full scan does. The margins
/// are the forward pass's own; the running totals are summed the first time a path needs them —
/// rarely, because a jump of both copies at once is the least likely transition — and serve every
/// path drawn through the site.
struct SiteDraw<'a> {
    autozygous: &'a [f64],
    pairs: &'a [f64],
    k: usize,
    rows: &'a [f64],
    columns: &'a [f64],
    autozygous_mass: f64,
    /// The pairs' total, the sum of the row totals.
    pair_mass: f64,
    /// The pairs' running totals, row by row, once a path has needed them.
    running: std::cell::OnceCell<Vec<f64>>,
}

impl<'a> SiteDraw<'a> {
    fn new(
        autozygous: &'a [f64],
        pairs: &'a [f64],
        rows: &'a [f64],
        columns: &'a [f64],
        k: usize,
    ) -> Self {
        Self {
            autozygous,
            pairs,
            k,
            rows,
            columns,
            autozygous_mass: autozygous.iter().sum(),
            pair_mass: rows.iter().sum(),
            running: std::cell::OnceCell::new(),
        }
    }

    fn pair_mass(&self) -> f64 {
        self.pair_mass
    }

    fn running(&self) -> &[f64] {
        self.running.get_or_init(|| {
            let mut total = 0.0;
            self.pairs
                .iter()
                .map(|value| {
                    total += value;
                    total
                })
                .collect()
        })
    }

    /// The state a uniform `u` draws here, given the next site's state and the persistence
    /// `(keep, stay autozygous, stay other)` between the two; `None` at the last site. A copy that
    /// jumped landed on template `j` with `landing[j]`.
    fn draw(&self, next: Option<((f64, f64, f64), State)>, u: f64, landing: &[f64]) -> State {
        let k = self.k;
        let any_template =
            |rest: f64| State::Autozygous(pick_at(self.autozygous.iter().copied(), rest));
        let any_pair = |rest: f64| {
            let position = pick_in_running(self.running(), rest);
            State::Pair(position / k, position % k)
        };
        match next {
            None => match choose_part(&[self.autozygous_mass, self.pair_mass()], u) {
                (0, rest) => any_template(rest),
                (_, rest) => any_pair(rest),
            },
            Some(((keep, stay_autozygous, stay_other), State::Pair(first, second))) => {
                let both_landed = landing[first] * landing[second];
                let to_autozygous = (1.0 - stay_autozygous) * both_landed;
                let none_kept = stay_other * (1.0 - keep) * (1.0 - keep) * both_landed;
                // Copy one kept its template and copy two landed on `second`, or the other way.
                let first_kept = stay_other * keep * (1.0 - keep) * landing[second];
                let second_kept = stay_other * keep * (1.0 - keep) * landing[first];
                let both_kept = stay_other * keep * keep;
                let masses = [
                    self.autozygous_mass * to_autozygous,
                    self.pair_mass() * none_kept,
                    self.rows[first] * first_kept,
                    self.columns[second] * second_kept,
                    self.pairs[first * k + second] * both_kept,
                ];
                match choose_part(&masses, u) {
                    (0, rest) => any_template(rest / to_autozygous),
                    (1, rest) => any_pair(rest / none_kept),
                    (2, rest) => {
                        let row = &self.pairs[first * k..][..k];
                        State::Pair(first, pick_at(row.iter().copied(), rest / first_kept))
                    }
                    (3, rest) => {
                        let column = self.pairs.iter().skip(second).step_by(k).copied();
                        State::Pair(pick_at(column, rest / second_kept), second)
                    }
                    _ => State::Pair(first, second),
                }
            }
            Some(((keep, stay_autozygous, stay_other), State::Autozygous(template))) => {
                let jumped = stay_autozygous * (1.0 - keep) * landing[template];
                let to_pairs = (1.0 - stay_other) * landing[template];
                let masses = [
                    self.autozygous_mass * jumped,
                    self.autozygous[template] * stay_autozygous * keep,
                    self.pair_mass() * to_pairs,
                ];
                match choose_part(&masses, u) {
                    (0, rest) => any_template(rest / jumped),
                    (1, _) => State::Autozygous(template),
                    (_, rest) => any_pair(rest / to_pairs),
                }
            }
        }
    }
}

/// **Which of a site's weight parts a uniform `u` falls in**, by their `masses` in order, and how
/// far into that part's mass it falls; a part of zero mass is never chosen.
///
/// # Panics
///
/// If every mass is zero.
fn choose_part(masses: &[f64], u: f64) -> (usize, f64) {
    let total: f64 = masses.iter().sum();
    assert!(total > 0.0, "a draw needs a positive weight");
    let target = u * total;
    let part = pick_at(masses.iter().copied(), target);
    let before: f64 = masses[..part].iter().sum();
    (part, (target - before).max(0.0))
}

/// **The two haplotypes along one drawn path** (phasing spec §2.1), each site's joint oriented by
/// the sample's reads where `reads` links it (reads-in-draw spec §2.2).
fn haplotypes_along(
    sites: &[WalkSite<'_>],
    terms: &[SiteTerms],
    states: &[State],
    reads: Option<&mut DrawReads<'_>>,
    joints: &mut Joints,
) -> PathHaplotypes {
    joints.clear();
    for ((site, site_terms), &state) in sites.iter().zip(terms).zip(states) {
        let alleles = site.allele_count;
        let start = joints.values.len();
        joints.starts.push(start);
        joints.values.resize(start + alleles * alleles, 0.0);
        let present = joint_into(
            state,
            &site_terms.blurred,
            &site_terms.by_alleles,
            alleles,
            &mut joints.values[start..],
        );
        joints.present.push(present);
    }
    if let Some(reads) = reads {
        orient_by_reads(sites, joints, reads);
    }
    let mut path = PathHaplotypes::with_capacity(sites.len());
    for (index, site) in sites.iter().enumerate() {
        copies_into(
            joints.get(index),
            site.allele_count,
            path.push_site(site.allele_count),
        );
    }
    path
}

/// **The heterozygote a chain reads at a site**, its two alleles the lower first, or `None` when
/// the site is not orientable: its heterozygous mass is below `het_mass_min` (reads-in-draw spec
/// §2.2). The pair is the one with the most mass over both orders, ties to the first visited.
fn orientable_pair(joint: &[f64], alleles: usize, het_mass_min: f64) -> Option<(usize, usize)> {
    let mut heterozygous = 0.0;
    let mut best = None;
    let mut most = f64::NEG_INFINITY;
    for a in 0..alleles {
        for b in a + 1..alleles {
            let mass = joint[a * alleles + b] + joint[b * alleles + a];
            heterozygous += mass;
            if mass > most {
                most = mass;
                best = Some((a, b));
            }
        }
    }
    best.filter(|_| heterozygous >= het_mass_min)
}

/// **Draws the orientation of every read-linked group along one path and moves each grouped
/// site's heterozygous mass onto it** (reads-in-draw spec §2.2). A group is a run of consecutive
/// orientable sites, each linked to the next by the sample's fragments; a site outside every group
/// of two or more keeps its joint.
fn orient_by_reads(sites: &[WalkSite<'_>], joints: &mut Joints, reads: &mut DrawReads<'_>) {
    let orientable: Vec<(usize, (usize, usize))> = (0..joints.len())
        .filter(|&site| !(reads.kept_out)(site))
        .filter_map(|site| {
            let joint = joints.get(site)?;
            orientable_pair(joint, sites[site].allele_count, reads.het_mass_min)
                .map(|pair| (site, pair))
        })
        .collect();
    let mut start = 0;
    let mut links = Vec::new();
    for index in 0..=orientable.len() {
        let link = (index > start && index < orientable.len())
            .then(|| {
                let (one, one_pair) = orientable[index - 1];
                let (other, other_pair) = orientable[index];
                (reads.link)(one, one_pair, other, other_pair)
            })
            .flatten();
        match link {
            Some(log_odds) => links.push(log_odds),
            None if index > start => {
                if index - start > 1 {
                    orient_group(
                        joints,
                        sites,
                        &orientable[start..index],
                        &links,
                        &mut reads.uniforms,
                    );
                }
                start = index;
                links.clear();
            }
            None => {}
        }
    }
}

/// **One group's chain** (reads-in-draw spec §2.2): own weights from each joint's two orders, the
/// link weight `e^{±L/2}` between neighbours, a forward pass normalised at each site, and a draw
/// backwards from the last; then each site's heterozygous mass moves onto the drawn order.
/// `links[i]` joins `group[i]` and `group[i + 1]`.
fn orient_group(
    joints: &mut Joints,
    sites: &[WalkSite<'_>],
    group: &[(usize, (usize, usize))],
    links: &[f64],
    uniforms: &mut Uniforms,
) {
    debug_assert_eq!(links.len() + 1, group.len(), "one link between neighbours");
    // Own weights: [straight, swapped], straight putting the lower allele on copy one.
    let own: Vec<[f64; 2]> = group
        .iter()
        .map(|&(site, _)| {
            let alleles = sites[site].allele_count;
            // PANIC-FREE: only sites with a joint are orientable.
            let joint = joints.get(site).expect("an orientable site has a joint");
            let mut weights = [0.0; 2];
            for a in 0..alleles {
                for b in a + 1..alleles {
                    weights[0] += joint[a * alleles + b];
                    weights[1] += joint[b * alleles + a];
                }
            }
            weights
        })
        .collect();
    let same = |log_odds: f64, one: usize, other: usize| {
        float::exp(if one == other { log_odds } else { -log_odds } / 2.0)
    };
    let mut forward: Vec<[f64; 2]> = Vec::with_capacity(group.len());
    for (index, weights) in own.iter().enumerate() {
        let mut next = *weights;
        if let Some(previous) = forward.last() {
            for (orientation, value) in next.iter_mut().enumerate() {
                *value *= (0..2)
                    .map(|before| previous[before] * same(links[index - 1], before, orientation))
                    .sum::<f64>();
            }
        }
        let total = next[0] + next[1];
        forward.push([next[0] / total, next[1] / total]);
    }
    let mut drawn = vec![0; group.len()];
    for index in (0..group.len()).rev() {
        let weights = match drawn.get(index + 1) {
            Some(&after) => [0, 1].map(|orientation| {
                forward[index][orientation] * same(links[index], orientation, after)
            }),
            None => forward[index],
        };
        drawn[index] = pick(&weights, uniforms.next_uniform());
    }
    for (&(site, _), &orientation) in group.iter().zip(&drawn) {
        let alleles = sites[site].allele_count;
        // PANIC-FREE: as above.
        let joint = joints
            .get_mut(site)
            .expect("an orientable site has a joint");
        for a in 0..alleles {
            for b in a + 1..alleles {
                let (kept, emptied) = if orientation == 0 {
                    (a * alleles + b, b * alleles + a)
                } else {
                    (b * alleles + a, a * alleles + b)
                };
                joint[kept] += joint[emptied];
                joint[emptied] = 0.0;
            }
        }
    }
}

/// **The two copies' haplotypes at one site, given the state the path is in there**: the joint of
/// copy one's allele `a` and copy two's `b`, `q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})` in a pair state
/// and `q̃_k(a) · L({a, a})` on the diagonal in an autozygous one, normalised; each copy is its
/// margin. `blurred` holds the templates' blurred shares, `by_alleles` the likelihood as a
/// symmetric `A × A` matrix. Where every pair is ruled out each copy is flat.
#[cfg(test)]
fn copies_along(
    state: State,
    blurred: &[f64],
    by_alleles: &[f64],
    alleles: usize,
) -> HaplotypePair {
    let mut joint = vec![0.0; alleles * alleles];
    let present = joint_into(state, blurred, by_alleles, alleles, &mut joint);
    let mut path = PathHaplotypes::with_capacity(1);
    copies_into(
        present.then_some(&joint[..]),
        alleles,
        path.push_site(alleles),
    );
    let [first, second] = path.site(0);
    [first.to_vec(), second.to_vec()]
}

/// **Each copy's haplotype, the margins of `joint`**, written over the two copies, which hold
/// zeros; flat where there is no joint.
fn copies_into(
    joint: Option<&[f64]>,
    alleles: usize,
    (first_copy, second_copy): (&mut [f64], &mut [f64]),
) {
    let Some(joint) = joint else {
        // Every allele pair ruled out — the state's emission was zero (§3.6's guard).
        first_copy.fill(1.0 / alleles as f64);
        second_copy.fill(1.0 / alleles as f64);
        return;
    };
    for a in 0..alleles {
        for b in 0..alleles {
            let value = joint[a * alleles + b];
            first_copy[a] += value;
            second_copy[b] += value;
        }
    }
}

/// **The joint of [`copies_along`] at one site**, normalised to sum to one, written over `joint`,
/// `A × A` zeros; `false` where every allele pair is ruled out.
fn joint_into(
    state: State,
    blurred: &[f64],
    by_alleles: &[f64],
    alleles: usize,
    joint: &mut [f64],
) -> bool {
    match state {
        State::Pair(first, second) => {
            for a in 0..alleles {
                for b in 0..alleles {
                    joint[a * alleles + b] = blurred[first * alleles + a]
                        * blurred[second * alleles + b]
                        * by_alleles[a * alleles + b];
                }
            }
        }
        State::Autozygous(template) => {
            for a in 0..alleles {
                joint[a * alleles + a] =
                    blurred[template * alleles + a] * by_alleles[a * alleles + a];
            }
        }
    }
    let total: f64 = joint.iter().sum();
    if total <= 0.0 {
        return false;
    }
    for value in joint.iter_mut() {
        *value /= total;
    }
    true
}

/// **A pair state's transition, split by how many copies keep their template**: both
/// (`keep²`), one (`keep · (1 − keep)`, the other jumping), or neither (`(1 − keep)²`, times the
/// mass it comes from), all times the chance of staying out of the autozygous state. Where a jump
/// lands is the landing weights' business, applied by the step that reads this.
#[derive(Clone, Copy)]
struct PairTransition {
    stay_other: f64,
    both_kept: f64,
    one_kept: f64,
    none_kept: f64,
}

impl PairTransition {
    fn new(keep: f64, pair_mass: f64, stay_other: f64) -> Self {
        Self {
            stay_other,
            both_kept: keep * keep,
            one_kept: keep * (1.0 - keep),
            none_kept: (1.0 - keep) * (1.0 - keep) * pair_mass,
        }
    }
}

/// **One forward step of the pair states** (spec §3.5), written over `out`: each pair's mass from
/// `values` (`k × k`, row by row) with its row and column sums, moved by `transition`, a jumping
/// copy landing on template `j` with `landing[j]`, plus `from_autozygous`, the mass leaving the
/// autozygous states, which lands on pair `(a, b)` with `landing[a] · landing[b]`.
fn forward_pair_step(
    values: &[f64],
    rows: &[f64],
    columns: &[f64],
    landing: &[f64],
    transition: PairTransition,
    from_autozygous: f64,
    out: &mut [f64],
) {
    let PairTransition {
        stay_other,
        both_kept,
        one_kept,
        none_kept,
    } = transition;
    let k = columns.len();
    let fresh = stay_other * none_kept + from_autozygous;
    for (((row, &row_sum), &first), out_row) in values
        .chunks_exact(k)
        .zip(rows)
        .zip(landing)
        .zip(out.chunks_exact_mut(k))
    {
        // Pair (a, b): its own mass kept, its row's mass with copy two landing on b, its column's
        // with copy one landing on a, and fresh mass landing on both.
        let row_term = stay_other * one_kept * row_sum + fresh * first;
        let column_factor = stay_other * one_kept * first;
        for (((out, &value), &column), &second) in
            out_row.iter_mut().zip(row).zip(columns).zip(landing)
        {
            *out = stay_other * both_kept * value + row_term * second + column_factor * column;
        }
    }
}

/// **One backward step of the pair states**, written over `out`: `values` the next site's
/// backward probabilities times its emissions, `rows` and `columns` their sums weighed by the
/// landing weights ([`landing_margins`]), `transition` built on their landing-weighed total, and
/// `from_autozygous` what the autozygous states, landing-weighed, add to every pair.
fn backward_pair_step(
    values: &[f64],
    rows: &[f64],
    columns: &[f64],
    transition: PairTransition,
    from_autozygous: f64,
    out: &mut [f64],
) {
    let PairTransition {
        stay_other,
        both_kept,
        one_kept,
        none_kept,
    } = transition;
    let k = columns.len();
    for ((row, &row_sum), out_row) in values
        .chunks_exact(k)
        .zip(rows)
        .zip(out.chunks_exact_mut(k))
    {
        let row_term = one_kept * row_sum + none_kept;
        for ((out, &value), &column) in out_row.iter_mut().zip(row).zip(columns) {
            *out =
                stay_other * (both_kept * value + row_term + one_kept * column) + from_autozygous;
        }
    }
}

/// Row and column sums of a `k × k` matrix stored row by row, each entry weighed by the landing
/// weight of the template it would land on: row `a` sums `m[a][b] · landing[b]`, column `b` sums
/// `landing[a] · m[a][b]`.
fn landing_margins(matrix: &[f64], landing: &[f64], rows: &mut [f64], columns: &mut [f64]) {
    columns.fill(0.0);
    for ((row, out), &first) in matrix
        .chunks_exact(landing.len())
        .zip(rows.iter_mut())
        .zip(landing)
    {
        let mut row_total = 0.0;
        for ((column, &value), &second) in columns.iter_mut().zip(row).zip(landing) {
            row_total += value * second;
            *column += first * value;
        }
        *out = row_total;
    }
}

/// Every value times `factor`. A loop of its own over one slice, which vectorises where a loop
/// over two slices chained together does not.
fn scale_by(values: &mut [f64], factor: f64) {
    for value in values {
        *value *= factor;
    }
}

/// How many running totals [`lane_sum`] and [`lane_dot`] keep.
const LANES: usize = 8;

/// **A sum in [`LANES`] lanes**: each lane adds every eighth value, then the lanes are added in
/// order and the values past the last whole eight after them. The order is fixed by the length
/// alone, so the bits are the same on every machine and pool width, and the lanes vectorise where
/// one running total — each addition waiting for the one before — cannot (`float_portability.md`).
fn lane_sum(values: &[f64]) -> f64 {
    let (chunks, rest) = values.as_chunks::<LANES>();
    let mut lanes = [0.0; LANES];
    for chunk in chunks {
        for (lane, value) in lanes.iter_mut().zip(chunk) {
            *lane += value;
        }
    }
    lanes.iter().sum::<f64>() + rest.iter().sum::<f64>()
}

/// **[`lane_sum`] of the products of `one` and `other`**, pair by pair.
fn lane_dot(one: &[f64], other: &[f64]) -> f64 {
    let ((one_chunks, one_rest), (other_chunks, other_rest)) =
        (one.as_chunks::<LANES>(), other.as_chunks::<LANES>());
    let mut lanes = [0.0; LANES];
    for (a, b) in one_chunks.iter().zip(other_chunks) {
        for ((lane, a), b) in lanes.iter_mut().zip(a).zip(b) {
            *lane += a * b;
        }
    }
    let rest: f64 = one_rest.iter().zip(other_rest).map(|(a, b)| a * b).sum();
    lanes.iter().sum::<f64>() + rest
}

/// Row and column sums of a `k × k` matrix stored row by row, each summed in index order, written
/// over `rows` and `columns`.
fn margins(matrix: &[f64], k: usize, rows: &mut [f64], columns: &mut [f64]) {
    columns.fill(0.0);
    for (row, row_total) in matrix.chunks_exact(k).zip(rows) {
        let mut total = 0.0;
        for (column, &value) in columns.iter_mut().zip(row) {
            total += value;
            *column += value;
        }
        *row_total = total;
    }
}

#[cfg(test)]
mod tests;
