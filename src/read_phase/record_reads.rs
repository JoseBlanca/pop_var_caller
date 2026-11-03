//! **Every sample's read alleles at one record** — what rides on a written record from calling to
//! the fragment keeper, numbered as the record numbers its alleles (read phase spec §2.3).

use super::capture::ReadAllele;

/// One read allele and the run sample whose read it was.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SampleReadAllele {
    /// The sample's index in the run's sample order.
    pub sample: u32,
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
    /// **The merge's lists, in the record's numbering**: `lists` holds `(run sample, that
    /// sample's read alleles)` in ascending sample order, each list in read order, as the merge's
    /// per-sample lists give them; `written` takes a merge allele index to
    /// the record's, `None` where candidate selection dropped the allele, and such a read allele
    /// is dropped with it (§2.3, step 1).
    ///
    /// # Panics
    ///
    /// If the samples are not in strictly ascending order: the merge lists covering samples that
    /// way, so anything else is a wiring defect.
    pub fn from_merge(
        lists: impl IntoIterator<Item = (usize, Vec<ReadAllele>)>,
        written: impl Fn(u32) -> Option<u32>,
    ) -> Self {
        let mut entries = Vec::new();
        let mut previous: Option<usize> = None;
        for (sample, list) in lists {
            assert!(
                previous.is_none_or(|previous| previous < sample),
                "covering samples out of order: {previous:?} then {sample}"
            );
            previous = Some(sample);
            let sample = u32::try_from(sample).expect("a run's samples are counted in a u32");
            entries.extend(list.into_iter().filter_map(|read_allele| {
                written(read_allele.allele).map(|allele| SampleReadAllele {
                    sample,
                    read_allele: ReadAllele {
                        allele,
                        ..read_allele
                    },
                })
            }));
        }
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

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests;
