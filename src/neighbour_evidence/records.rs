//! **One record as the pass holds it** (`doc/devel/ng/spec/neighbour_evidence_pass.md` §3.2): what
//! every sample's walk reads at it — each sample's genotype likelihoods from its `PL` and allele
//! shares from its `GP` — and the background template's shares from the record's `AF` (§3.3).
//!
//! The numbers are kept as `f32`, sample after sample in one array a field, because a window holds
//! one of these for every sample at every record and that is what bounds the pass's memory (§6). The
//! inputs carry far less precision than `f32` does: `PL` is in whole Phred units, `GP` at three
//! decimals.

use crate::calling::genotype_table::GenotypeTable;
use crate::float;
use crate::types::Ploidy;

/// **What one sample's column says at a record.**
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleRow {
    /// Written `.` in both `PL` and `GP` — set aside by the candidate step. The sample is not
    /// re-called here, and as a template it carries equal shares, no information (§3.2).
    Absent,
    /// A `PL` with every genotype equal: the sample's reads said nothing here, and calling wrote
    /// `./.` with its prior as `GP` (§4.2).
    Flat,
    /// A `PL` that is not flat: the sample's own reads said something — an *informative* record of
    /// the sample (§1.3).
    Informative,
}

/// One sample's `PL` and `GP` at a record, as parsed from its column; `None` for a column written
/// `.`.
#[derive(Clone, Copy, Debug)]
pub struct SampleColumns<'a> {
    /// `PL`, in the genotype table's order.
    pub phred_likelihoods: &'a [u32],
    /// `GP`, in the genotype table's order.
    pub genotype_posteriors: &'a [f64],
}

/// **One PASS record with rows, as every walk reads it.**
#[derive(Clone, Debug)]
pub struct SiteRecord {
    /// The record's written `POS`.
    pub position: u64,
    /// How many alleles the record writes, the reference included.
    pub allele_count: usize,
    /// The background template's share of each allele: the cohort's frequencies (§3.3).
    background_shares: Vec<f32>,
    /// What each sample's column says.
    rows: Vec<SampleRow>,
    /// Each sample's genotype likelihoods, `10^(−PL/10)`, `genotype_count` numbers a sample; zeros
    /// where the sample is [`SampleRow::Absent`].
    likelihoods: Vec<f32>,
    /// Each sample's expected share of its copies carrying each allele, from its `GP`,
    /// `allele_count` numbers a sample; zeros where the sample is [`SampleRow::Absent`].
    shares: Vec<f32>,
}

impl SiteRecord {
    /// Builds a record from its `POS`, its `AF` (one frequency per alternative allele) and each
    /// sample's `PL` and `GP` — `None` for a column written `.`.
    ///
    /// # Panics
    ///
    /// If a sample's `PL` or `GP` does not hold one number per diploid genotype over
    /// `alt_frequencies.len() + 1` alleles.
    #[must_use]
    pub fn new(
        position: u64,
        alt_frequencies: &[f64],
        samples: &[Option<SampleColumns<'_>>],
    ) -> Self {
        let allele_count = alt_frequencies.len() + 1;
        let diploid = Ploidy::try_new(2).expect("two copies is a ploidy");
        let table = GenotypeTable::build(diploid, allele_count);
        let table = table.view();
        let genotype_count = table.genotype_count();
        let copies_of = table.genotype_allele_counts();

        let mut background = Vec::with_capacity(allele_count);
        background_shares(alt_frequencies, &mut background);
        let mut rows = Vec::with_capacity(samples.len());
        let mut likelihoods = vec![0.0_f32; samples.len() * genotype_count];
        let mut shares = vec![0.0_f32; samples.len() * allele_count];
        for (sample, columns) in samples.iter().enumerate() {
            let Some(columns) = columns else {
                rows.push(SampleRow::Absent);
                continue;
            };
            assert!(
                columns.phred_likelihoods.len() == genotype_count
                    && columns.genotype_posteriors.len() == genotype_count,
                "a record at {position} with {allele_count} alleles has {genotype_count} genotypes; \
                 sample {sample} has {} PL and {} GP values",
                columns.phred_likelihoods.len(),
                columns.genotype_posteriors.len()
            );
            let lowest = columns.phred_likelihoods.iter().copied().min().unwrap_or(0);
            let flat = columns
                .phred_likelihoods
                .iter()
                .all(|&phred| phred == lowest);
            rows.push(if flat {
                SampleRow::Flat
            } else {
                SampleRow::Informative
            });
            for (out, &phred) in likelihoods[sample * genotype_count..(sample + 1) * genotype_count]
                .iter_mut()
                .zip(columns.phred_likelihoods)
            {
                *out = float::powf(10.0, -f64::from(phred - lowest) / 10.0) as f32;
            }
            // q(a) = Σ_g GP_g · copies_a(g) / 2, rescaled because GP is written rounded.
            let mut sample_shares = vec![0.0_f64; allele_count];
            for (genotype, &posterior) in columns.genotype_posteriors.iter().enumerate() {
                let copies = &copies_of[genotype * allele_count..(genotype + 1) * allele_count];
                for (share, &count) in sample_shares.iter_mut().zip(copies) {
                    *share += posterior * f64::from(count) / 2.0;
                }
            }
            let total: f64 = sample_shares.iter().sum();
            for (out, share) in shares[sample * allele_count..(sample + 1) * allele_count]
                .iter_mut()
                .zip(&sample_shares)
            {
                *out = if total > 0.0 {
                    (share / total) as f32
                } else {
                    1.0 / allele_count as f32
                };
            }
        }
        Self {
            position,
            allele_count,
            background_shares: background.iter().map(|&share| share as f32).collect(),
            rows,
            likelihoods,
            shares,
        }
    }

    /// How many diploid genotypes the record has.
    #[must_use]
    pub fn genotype_count(&self) -> usize {
        self.allele_count * (self.allele_count + 1) / 2
    }

    /// **The bytes the record holds on the heap** — what a window pays for it besides its line and
    /// its reads (pass spec §6). Counted from the arrays' lengths, not their capacities, and a row
    /// at one byte: the count places the windows, so it must not depend on how the allocator or
    /// the compiler sized anything.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        (self.background_shares.len() + self.likelihoods.len() + self.shares.len())
            * std::mem::size_of::<f32>()
            + self.rows.len()
    }

    /// How many samples the record has a column for.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.rows.len()
    }

    /// What sample `sample`'s column says.
    #[must_use]
    pub fn row(&self, sample: usize) -> SampleRow {
        self.rows[sample]
    }

    /// **Whether any walk uses this record** (§5): it needs at least two samples with rows, one to
    /// re-call and one to copy from.
    #[must_use]
    pub fn is_site(&self) -> bool {
        self.rows
            .iter()
            .filter(|&&row| row != SampleRow::Absent)
            .take(2)
            .count()
            == 2
    }

    /// Sample `sample`'s genotype likelihoods; flat where the sample is absent, which in a walk is
    /// an emission of one in every state.
    pub fn likelihoods(&self, sample: usize) -> impl Iterator<Item = f64> + '_ {
        let genotypes = self.genotype_count();
        let absent = self.rows[sample] == SampleRow::Absent;
        self.likelihoods[sample * genotypes..(sample + 1) * genotypes]
            .iter()
            .map(move |&likelihood| if absent { 1.0 } else { f64::from(likelihood) })
    }

    /// Sample `sample`'s allele shares as a template; `1/A` each where the sample is absent
    /// (§3.2).
    pub fn template_shares(&self, sample: usize) -> impl Iterator<Item = f64> + '_ {
        let alleles = self.allele_count;
        let absent = self.rows[sample] == SampleRow::Absent;
        self.shares[sample * alleles..(sample + 1) * alleles]
            .iter()
            .map(move |&share| {
                if absent {
                    1.0 / alleles as f64
                } else {
                    f64::from(share)
                }
            })
    }

    /// The background template's allele shares (§3.3).
    pub fn background_shares(&self) -> impl Iterator<Item = f64> + '_ {
        self.background_shares.iter().map(|&share| f64::from(share))
    }
}

/// **The background template's allele shares** at a record (spec §3.3): the cohort's allele
/// frequencies, the reference taking what the alternatives leave. `AF` is written rounded, so the
/// alternatives can sum past one; the reference is then 0 and the shares are rescaled to sum to one.
/// Written into `shares`, which is cleared first.
pub fn background_shares(alt_frequencies: &[f64], shares: &mut Vec<f64>) {
    shares.clear();
    let alternatives: f64 = alt_frequencies.iter().sum();
    shares.push((1.0 - alternatives).max(0.0));
    shares.extend(alt_frequencies.iter().map(|frequency| frequency.max(0.0)));
    let total: f64 = shares.iter().sum();
    let count = shares.len() as f64;
    for share in shares.iter_mut() {
        *share = if total > 0.0 {
            *share / total
        } else {
            1.0 / count
        };
    }
}

#[cfg(test)]
mod tests;
