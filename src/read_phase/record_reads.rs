//! **Every sample's read alleles at one record** — what rides on a written record from calling to
//! the fragment keeper, numbered as the record numbers its alleles (read phase spec §2.3).

use super::capture::ReadAllele;

/// One read allele and the run sample whose read it was.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SampleReadAllele {
    /// The sample's index in the run's sample order.
    pub sample: u32,
    /// What the sample's read showed.
    pub read_allele: ReadAllele,
}

/// **Every sample's read alleles at one record**, in one list sorted by sample and then by read.
///
/// **One flat list rather than one per sample**, because most of a cohort's samples carry none at
/// most records — a sample without reads there, or an STR record, which names no reads — and a
/// list per run sample would cost its header for each of them, 72 kB a record at 3,000 samples.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RecordReadAlleles {
    entries: Vec<SampleReadAllele>,
}

impl RecordReadAlleles {
    /// **The merge's read alleles, in the record's numbering**: `entries` holds every covering
    /// sample's read alleles, each naming its run sample, sorted by sample and then by read, as
    /// the merge's capture gives them; `written` takes a merge allele index to the record's,
    /// `None` where candidate selection dropped the allele, and such a read allele is dropped
    /// with it (§2.3, step 1).
    ///
    /// # Panics
    ///
    /// If the entries are not sorted by sample and then by read: the merge captures them that
    /// way, so anything else is a wiring defect.
    pub fn from_merge(
        mut entries: Vec<SampleReadAllele>,
        written: impl Fn(u32) -> Option<u32>,
    ) -> Self {
        assert!(
            entries.is_sorted(),
            "read alleles out of order: the merge lists them by sample and then by read"
        );
        entries.retain_mut(|entry| match written(entry.read_allele.allele) {
            Some(allele) => {
                entry.read_allele.allele = allele;
                true
            }
            None => false,
        });
        Self { entries }
    }

    /// **Renumber the alleles as the record renumbered its own** (§2.3, step 2): `renumbered[old]`
    /// is the new index, or `None` where the record dropped that allele, and a read allele on it is
    /// dropped too.
    pub fn renumber(&mut self, renumbered: &[Option<usize>]) {
        self.entries.retain_mut(|entry| {
            match renumbered
                .get(entry.read_allele.allele as usize)
                .copied()
                .flatten()
            {
                Some(new) => {
                    entry.read_allele.allele = new as u32;
                    true
                }
                None => false,
            }
        });
    }

    /// Every entry, sorted by sample and then by read.
    #[must_use]
    pub fn entries(&self) -> &[SampleReadAllele] {
        &self.entries
    }

    /// One sample's read alleles, in read order.
    #[must_use]
    pub fn of_sample(&self, sample: u32) -> &[SampleReadAllele] {
        let start = self.entries.partition_point(|entry| entry.sample < sample);
        let end = self.entries.partition_point(|entry| entry.sample <= sample);
        &self.entries[start..end]
    }

    /// How many read alleles the record carries, over every sample.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the record carries none.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests;
