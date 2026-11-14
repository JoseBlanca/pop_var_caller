//! **Fragments** — a sample's read alleles once the fragment keeper has numbered them, the form in
//! which reads travel from the ordered thread to whichever step writes phase (read phase spec §2.4).
//!
//! These are the interchange between stages: the keeper (`run::fragment_keeper`) makes them, the
//! writer, the paralog spill and both phase writers read them.

/// **`R`, the reach**: the longest distance in bases between a read's first read allele and a later
/// one that is still kept (§2.4), a starting value adopted by the owner on 2026-11-03. It covers
/// ordinary paired-end inserts — the founder simulator draws fragments of 350 ± 50 bases, and
/// HapCUT2 joins mates up to 1,000 bases by default.
pub const READ_REACH: u64 = 1_000;

/// One sample's read allele at one record, its read renamed to a fragment number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FragmentAllele {
    /// The sample's index in the run's sample order.
    pub sample: u32,
    /// **The read, as a number counted from 0 within each contig**, in the order reads were first
    /// seen. Never written anywhere; it means nothing outside its sample and contig.
    pub fragment: u32,
    /// The allele, numbered as the record numbers its alleles.
    pub allele: u16,
    /// The read's error as a Phred integer (§2.2).
    pub error_phred: u8,
}

/// **Every sample's fragments at one record**, sorted by sample and then fragment — the order the
/// link sums run in, so they are the same at any thread count (§2.5).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RecordFragments {
    entries: Vec<FragmentAllele>,
}

impl RecordFragments {
    /// From entries in any order; they are sorted here.
    #[must_use]
    pub fn from_entries(mut entries: Vec<FragmentAllele>) -> Self {
        entries.sort_unstable();
        Self { entries }
    }

    /// Every entry, sorted by sample and then fragment.
    #[must_use]
    pub fn entries(&self) -> &[FragmentAllele] {
        &self.entries
    }

    /// One sample's fragments here, in fragment order.
    #[must_use]
    pub fn of_sample(&self, sample: u32) -> &[FragmentAllele] {
        let start = self.entries.partition_point(|entry| entry.sample < sample);
        let end = self.entries.partition_point(|entry| entry.sample <= sample);
        &self.entries[start..end]
    }

    /// How many fragment alleles the record carries, over every sample.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the record carries none: no sample's read here reached another record.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// **One read allele as a window of the neighbour-evidence pass holds it**: a [`FragmentAllele`]
/// without its sample, which [`PackedRecordFragments`] keeps once per sample rather than once per
/// read allele. Packed to seven bytes, against [`FragmentAllele`]'s twelve: a window holds every
/// read allele of every sample at every record it spans, and on tomato that was the largest single
/// share of its memory (`neighbour_evidence_pass.md` §6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C, packed)]
pub struct ReadAllele {
    fragment: u32,
    allele: u16,
    error_phred: u8,
}

impl ReadAllele {
    /// The read, numbered as [`FragmentAllele::fragment`] numbers it.
    #[must_use]
    pub fn fragment(self) -> u32 {
        self.fragment
    }

    /// The allele, numbered as the record numbers its alleles.
    #[must_use]
    pub fn allele(self) -> u16 {
        self.allele
    }

    /// The read's error as a Phred integer.
    #[must_use]
    pub fn error_phred(self) -> u8 {
        self.error_phred
    }
}

/// **Every sample's read alleles at one record, packed**: [`RecordFragments`] with each sample
/// named once, beside where its entries start, instead of in every entry. The same entries in the
/// same order, so every sum over them is the same.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PackedRecordFragments {
    /// Each sample with entries, ascending, and the index of its first entry.
    samples: Vec<(u32, u32)>,
    entries: Vec<ReadAllele>,
}

impl PackedRecordFragments {
    /// One sample's read alleles here, in fragment order.
    #[must_use]
    pub fn of_sample(&self, sample: u32) -> &[ReadAllele] {
        let Ok(at) = self
            .samples
            .binary_search_by_key(&sample, |&(sample, _)| sample)
        else {
            return &[];
        };
        let start = self.samples[at].1 as usize;
        let end = self
            .samples
            .get(at + 1)
            .map_or(self.entries.len(), |&(_, start)| start as usize);
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

    /// The bytes it holds on the heap, from the arrays' lengths: the count places refinement's
    /// windows, so it must not depend on how the allocator sized them.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        self.samples.len() * 2 * std::mem::size_of::<u32>()
            + self.entries.len() * std::mem::size_of::<ReadAllele>()
    }
}

impl From<RecordFragments> for PackedRecordFragments {
    fn from(fragments: RecordFragments) -> Self {
        let mut samples: Vec<(u32, u32)> = Vec::new();
        let mut entries = Vec::with_capacity(fragments.len());
        for entry in fragments.entries() {
            if samples
                .last()
                .is_none_or(|&(sample, _)| sample != entry.sample)
            {
                // PANIC-FREE: a record holds far fewer than 2³² read alleles.
                let start =
                    u32::try_from(entries.len()).expect("a record's read alleles fit a u32");
                samples.push((entry.sample, start));
            }
            entries.push(ReadAllele {
                fragment: entry.fragment,
                allele: entry.allele,
                error_phred: entry.error_phred,
            });
        }
        samples.shrink_to_fit();
        Self { samples, entries }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Packing keeps every sample's entries, in order, and loses nothing else**: each sample's
    /// packed read alleles are its unpacked ones less the sample, a sample with none gets none,
    /// and the seven-byte size the window's memory is planned with holds.
    #[test]
    fn packing_keeps_each_samples_read_alleles_in_order() {
        let entry = |sample: u32, fragment: u32, allele: u16, error_phred: u8| FragmentAllele {
            sample,
            fragment,
            allele,
            error_phred,
        };
        let fragments = RecordFragments::from_entries(vec![
            entry(7, 3, 1, 20),
            entry(2, 9, 0, 30),
            entry(7, 1, 0, 40),
            entry(2, 4, 300, 5),
            entry(11, 0, 2, 35),
        ]);
        let packed = PackedRecordFragments::from(fragments.clone());
        assert_eq!(packed.len(), fragments.len());
        for sample in 0..13 {
            let unpacked: Vec<(u32, u16, u8)> = fragments
                .of_sample(sample)
                .iter()
                .map(|entry| (entry.fragment, entry.allele, entry.error_phred))
                .collect();
            let repacked: Vec<(u32, u16, u8)> = packed
                .of_sample(sample)
                .iter()
                .map(|entry| (entry.fragment(), entry.allele(), entry.error_phred()))
                .collect();
            assert_eq!(repacked, unpacked, "sample {sample}");
        }
        assert_eq!(std::mem::size_of::<ReadAllele>(), 7);
        assert!(PackedRecordFragments::default().of_sample(0).is_empty());
    }
}
