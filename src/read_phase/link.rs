//! **A link's log-odds from reads** — how much more likely the fragments covering two heterozygous
//! records of a sample are under *cis* than under *trans* (read phase spec §2.5).
//!
//! A link is **cis** when the lower-numbered alleles of the two genotypes lie on one copy, and
//! **trans** otherwise. A fragment that showed allele `xᵢ` at one record and `xⱼ` at the other says
//! cis when `xᵢ` is the lower allele exactly when `xⱼ` is, and says trans otherwise; one that showed
//! an allele outside either genotype says nothing. Its weight is HapCUT2's edge weight: the log of
//! the chance it shows what it shows if its claim is true, over the chance if false. Both phase
//! writers use this module — the reads-only stage and the neighbour-evidence pass.

use std::sync::OnceLock;

use super::capture::LARGEST_ERROR_PHRED;
use super::keeper::FragmentAllele;

/// **`T_link`, ln 199**: reads make or keep a link at odds of at least 199 to 1, one wrong link in
/// 200 if the error model holds (§2.5). With the error floor at 0.01 one fragment is worth at most
/// 3.90 nats, so a link needs two good fragments that agree.
pub const LINK_MIN: f64 = 5.293_304_824_724_492;

/// **`T_break`, ln 9**: reads that contradict a population link at odds of 9 to 1 or more break it
/// (§2.5, §2.7). Lower than [`LINK_MIN`] because breaking a link costs block length and a wrong
/// link costs much more; one fragment can reach it.
pub const BREAK_MIN: f64 = 2.197_224_577_336_219_6;

/// **`e_floor`**: the error per read allele that no base or mapping quality sees — a chimeric
/// fragment, a PCR error, a base the geometric mean of §2.2 rates too well — and that strikes
/// reads independently. A starting value, to be fitted on deep data (§7.2). **Not** a guard
/// against a collapsed duplication, whose reads err together: the hidden-duplication filter is.
pub const ERROR_FLOOR: f64 = 0.01;

/// **`GQ_link`**: the written `GQ` both ends of a link resting on reads must have (§2.5). A
/// heterozygote that is really homozygous gives a confident link of random sign, since every
/// fragment shows the same allele there; `GQ` 30 bounds each end's chance of that at 1 in 1,000.
pub const GQ_MIN: u8 = 30;

/// **How a read allele's Phred becomes its error**: `min(max(10^(−q/10), e_floor), ½)`.
///
/// The ceiling matters as much as the floor. Above one half a read allele would vote for the
/// opposite relation, and Phred 0 to 3 is real input — the walk gives the losing mate of an
/// overlapping pair base quality 0. Clamped at one half, such a read allele weighs exactly zero.
#[must_use]
pub fn error_of(error_phred: u8) -> f64 {
    let error = crate::float::powf(10.0, -f64::from(error_phred) / 10.0);
    error.clamp(ERROR_FLOOR, 0.5)
}

/// **One fragment's weight**, `ln([(1 − eᵢ)(1 − eⱼ) + eᵢeⱼ] / [eᵢ(1 − eⱼ) + eⱼ(1 − eᵢ)])`: never
/// negative, since neither error exceeds one half, and zero when either is one half.
#[must_use]
pub fn weight(error_phred_i: u8, error_phred_j: u8) -> f64 {
    let table = WEIGHTS.get_or_init(|| {
        let mut table = [[0.0; PHREDS]; PHREDS];
        for (i, row) in table.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let (e_i, e_j) = (error_of(i as u8), error_of(j as u8));
                let agree = (1.0 - e_i) * (1.0 - e_j) + e_i * e_j;
                let disagree = e_i * (1.0 - e_j) + e_j * (1.0 - e_i);
                *cell = crate::float::ln(agree / disagree);
            }
        }
        table
    });
    let clamp = |phred: u8| usize::from(phred.min(LARGEST_ERROR_PHRED));
    table[clamp(error_phred_i)][clamp(error_phred_j)]
}

/// Every stored Phred, 0 to [`LARGEST_ERROR_PHRED`].
const PHREDS: usize = LARGEST_ERROR_PHRED as usize + 1;

/// Each pair of Phreds' weight, computed once: the link sums look a weight up per fragment.
static WEIGHTS: OnceLock<[[f64; PHREDS]; PHREDS]> = OnceLock::new();

/// A heterozygous genotype's two alleles, the lower first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Heterozygote {
    pub lower: u16,
    pub upper: u16,
}

impl Heterozygote {
    /// The two alleles in either order; `None` when they are one allele.
    #[must_use]
    pub fn of(one: u16, other: u16) -> Option<Self> {
        match one.cmp(&other) {
            std::cmp::Ordering::Less => Some(Self {
                lower: one,
                upper: other,
            }),
            std::cmp::Ordering::Greater => Some(Self {
                lower: other,
                upper: one,
            }),
            std::cmp::Ordering::Equal => None,
        }
    }

    fn holds(self, allele: u16) -> bool {
        allele == self.lower || allele == self.upper
    }
}

/// **What the reads say about one link**: the log-odds of cis against trans, and how many
/// fragments said anything.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadLink {
    /// `L`: positive favours cis, negative trans.
    pub log_odds: f64,
    /// The fragments that covered both records with an allele of each genotype.
    pub fragments: u32,
}

/// The relation a link states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Relation {
    /// The lower alleles on one copy.
    Cis,
    /// The lower allele of one with the upper allele of the other.
    Trans,
}

impl ReadLink {
    /// **The relation the reads make**: `Some` when `|L|` reaches [`LINK_MIN`].
    #[must_use]
    pub fn relation(self) -> Option<Relation> {
        if self.log_odds >= LINK_MIN {
            Some(Relation::Cis)
        } else if self.log_odds <= -LINK_MIN {
            Some(Relation::Trans)
        } else {
            None
        }
    }

    /// **Whether the reads contradict `relation`** at [`BREAK_MIN`] or more.
    #[must_use]
    pub fn contradicts(self, relation: Relation) -> bool {
        match relation {
            Relation::Cis => self.log_odds <= -BREAK_MIN,
            Relation::Trans => self.log_odds >= BREAK_MIN,
        }
    }
}

/// **`L` for one sample's link between two records**, from that sample's fragments at each, both
/// in fragment order as [`RecordFragments::of_sample`](super::keeper::RecordFragments::of_sample)
/// gives them.
///
/// The sum runs in fragment order, in `f64`, so it is the same at any thread count (§2.5).
#[must_use]
pub fn read_link(
    at_one: &[FragmentAllele],
    one: Heterozygote,
    at_other: &[FragmentAllele],
    other: Heterozygote,
) -> ReadLink {
    let mut link = ReadLink {
        log_odds: 0.0,
        fragments: 0,
    };
    let (mut i, mut j) = (0, 0);
    while i < at_one.len() && j < at_other.len() {
        let (x, y) = (at_one[i], at_other[j]);
        match x.fragment.cmp(&y.fragment) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                if one.holds(x.allele) && other.holds(y.allele) {
                    let w = weight(x.error_phred, y.error_phred);
                    let cis = (x.allele == one.lower) == (y.allele == other.lower);
                    link.log_odds += if cis { w } else { -w };
                    link.fragments += 1;
                }
                i += 1;
                j += 1;
            }
        }
    }
    link
}

#[cfg(test)]
mod tests;
