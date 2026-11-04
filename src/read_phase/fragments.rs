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
