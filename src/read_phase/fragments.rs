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
///
/// **Three bytes a read allele where they fit** (pass spec §6): a sample's fragment numbers at one
/// record are close together — its reads there were first seen within about a read's reach — so
/// each is held as its distance from the sample's first, in a `u16`, beside the allele in two bits
/// and the Phred in six ([`TightReadAllele`]). A record with an allele past 3 holds them in four
/// ([`CompactReadAllele`]: a `u8` allele and a `u8` Phred); one where some sample's fragments span
/// more than a `u16`, or an allele is numbered past a `u8`, holds every read allele whole
/// ([`ReadAllele`], seven bytes). Nothing is lost either way.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PackedRecordFragments {
    /// Each sample with entries, ascending.
    samples: Vec<PackedSample>,
    entries: PackedEntries,
}

/// One sample's place among a [`PackedRecordFragments`]'s entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PackedSample {
    sample: u32,
    /// The index of its first entry.
    start: u32,
    /// The fragment its compact entries count from: its first entry's.
    first_fragment: u32,
}

/// A record's read alleles, in the smallest form every one of them fits.
#[derive(Debug, Clone, PartialEq, Eq)]
enum PackedEntries {
    Tight(Vec<TightReadAllele>),
    Compact(Vec<CompactReadAllele>),
    Whole(Vec<ReadAllele>),
}

impl Default for PackedEntries {
    fn default() -> Self {
        Self::Tight(Vec::new())
    }
}

/// **A read allele in three bytes**: its fragment's distance from its sample's first at the record,
/// then its allele in the top two bits of the last byte and its error Phred in the other six.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TightReadAllele {
    fragment_offset: [u8; 2],
    allele_and_phred: u8,
}

/// The largest allele a [`TightReadAllele`] holds: two bits.
const TIGHT_MAX_ALLELE: u8 = 3;

/// The largest error Phred a [`TightReadAllele`] holds: six bits, past
/// [`LARGEST_ERROR_PHRED`](super::capture::LARGEST_ERROR_PHRED)'s 60.
const TIGHT_MAX_PHRED: u8 = 63;

impl TightReadAllele {
    /// The entry in three bytes, if its allele and Phred fit.
    fn of(compact: CompactReadAllele) -> Option<Self> {
        (compact.allele <= TIGHT_MAX_ALLELE && compact.error_phred <= TIGHT_MAX_PHRED).then(|| {
            Self {
                fragment_offset: compact.fragment_offset.to_le_bytes(),
                allele_and_phred: (compact.allele << 6) | compact.error_phred,
            }
        })
    }
}

/// **A read allele in four bytes**: its fragment as the distance from its sample's first at the
/// record, its allele and its Phred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompactReadAllele {
    fragment_offset: u16,
    allele: u8,
    error_phred: u8,
}

/// **One sample's read alleles at one record**, in fragment order, as a
/// [`PackedRecordFragments`] holds them: read each with [`SampleReadAlleles::get`].
#[derive(Debug, Clone, Copy)]
pub enum SampleReadAlleles<'a> {
    /// Three-byte entries and the fragment their offsets count from.
    Tight {
        first_fragment: u32,
        entries: &'a [TightReadAllele],
    },
    /// Compact entries and the fragment their offsets count from.
    Compact {
        first_fragment: u32,
        entries: &'a [CompactReadAllele],
    },
    /// Whole entries.
    Whole(&'a [ReadAllele]),
}

impl<'a> SampleReadAlleles<'a> {
    /// No read alleles.
    pub const EMPTY: SampleReadAlleles<'static> = SampleReadAlleles::Whole(&[]);

    /// How many there are.
    #[must_use]
    pub fn len(self) -> usize {
        match self {
            Self::Tight { entries, .. } => entries.len(),
            Self::Compact { entries, .. } => entries.len(),
            Self::Whole(entries) => entries.len(),
        }
    }

    /// Whether there are none.
    #[must_use]
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// The `index`th, whole.
    ///
    /// # Panics
    ///
    /// If `index` is not below [`Self::len`].
    #[must_use]
    pub fn get(self, index: usize) -> ReadAllele {
        match self {
            Self::Tight {
                first_fragment,
                entries,
            } => {
                let entry = entries[index];
                ReadAllele {
                    fragment: first_fragment + u32::from(u16::from_le_bytes(entry.fragment_offset)),
                    allele: u16::from(entry.allele_and_phred >> 6),
                    error_phred: entry.allele_and_phred & TIGHT_MAX_PHRED,
                }
            }
            Self::Compact {
                first_fragment,
                entries,
            } => {
                let entry = entries[index];
                ReadAllele {
                    fragment: first_fragment + u32::from(entry.fragment_offset),
                    allele: u16::from(entry.allele),
                    error_phred: entry.error_phred,
                }
            }
            Self::Whole(entries) => entries[index],
        }
    }

    /// Every one, whole, in order.
    pub fn iter(self) -> impl Iterator<Item = ReadAllele> + 'a {
        (0..self.len()).map(move |index| self.get(index))
    }
}

impl PackedRecordFragments {
    /// One sample's read alleles here, in fragment order.
    #[must_use]
    pub fn of_sample(&self, sample: u32) -> SampleReadAlleles<'_> {
        let Ok(at) = self
            .samples
            .binary_search_by_key(&sample, |packed| packed.sample)
        else {
            return SampleReadAlleles::EMPTY;
        };
        let start = self.samples[at].start as usize;
        let end = self
            .samples
            .get(at + 1)
            .map_or(self.len(), |next| next.start as usize);
        match &self.entries {
            PackedEntries::Tight(entries) => SampleReadAlleles::Tight {
                first_fragment: self.samples[at].first_fragment,
                entries: &entries[start..end],
            },
            PackedEntries::Compact(entries) => SampleReadAlleles::Compact {
                first_fragment: self.samples[at].first_fragment,
                entries: &entries[start..end],
            },
            PackedEntries::Whole(entries) => SampleReadAlleles::Whole(&entries[start..end]),
        }
    }

    /// How many read alleles the record carries, over every sample.
    #[must_use]
    pub fn len(&self) -> usize {
        match &self.entries {
            PackedEntries::Tight(entries) => entries.len(),
            PackedEntries::Compact(entries) => entries.len(),
            PackedEntries::Whole(entries) => entries.len(),
        }
    }

    /// Whether the record carries none.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether its read alleles are held in four bytes each or fewer.
    #[must_use]
    pub fn is_compact(&self) -> bool {
        matches!(
            self.entries,
            PackedEntries::Tight(_) | PackedEntries::Compact(_)
        )
    }

    /// Whether its read alleles are held in three bytes each.
    #[must_use]
    pub fn is_tight(&self) -> bool {
        matches!(self.entries, PackedEntries::Tight(_))
    }

    /// The bytes it holds on the heap, from the arrays' lengths: the count places refinement's
    /// windows, so it must not depend on how the allocator sized them.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        let entry_bytes = match &self.entries {
            PackedEntries::Tight(_) => std::mem::size_of::<TightReadAllele>(),
            PackedEntries::Compact(_) => std::mem::size_of::<CompactReadAllele>(),
            PackedEntries::Whole(_) => std::mem::size_of::<ReadAllele>(),
        };
        self.samples.len() * 3 * std::mem::size_of::<u32>() + self.len() * entry_bytes
    }
}

impl From<RecordFragments> for PackedRecordFragments {
    fn from(fragments: RecordFragments) -> Self {
        let mut samples: Vec<PackedSample> = Vec::new();
        for (index, entry) in fragments.entries().iter().enumerate() {
            if samples
                .last()
                .is_none_or(|last| last.sample != entry.sample)
            {
                // PANIC-FREE: a record holds far fewer than 2³² read alleles.
                let start = u32::try_from(index).expect("a record's read alleles fit a u32");
                samples.push(PackedSample {
                    sample: entry.sample,
                    start,
                    first_fragment: entry.fragment,
                });
            }
        }
        samples.shrink_to_fit();
        // Entries are sorted by sample and then fragment, so each sample's first is its lowest.
        let mut sample_at = 0;
        let mut compact = Vec::with_capacity(fragments.len());
        for (index, entry) in fragments.entries().iter().enumerate() {
            while samples
                .get(sample_at + 1)
                .is_some_and(|next| next.start as usize <= index)
            {
                sample_at += 1;
            }
            let offset = u16::try_from(entry.fragment - samples[sample_at].first_fragment);
            let allele = u8::try_from(entry.allele);
            let (Ok(fragment_offset), Ok(allele)) = (offset, allele) else {
                let whole = fragments
                    .entries()
                    .iter()
                    .map(|entry| ReadAllele {
                        fragment: entry.fragment,
                        allele: entry.allele,
                        error_phred: entry.error_phred,
                    })
                    .collect();
                return Self {
                    samples,
                    entries: PackedEntries::Whole(whole),
                };
            };
            compact.push(CompactReadAllele {
                fragment_offset,
                allele,
                error_phred: entry.error_phred,
            });
        }
        // Three bytes an entry where every entry's allele and Phred fit, as most records' do.
        let tight: Option<Vec<TightReadAllele>> =
            compact.iter().copied().map(TightReadAllele::of).collect();
        let entries = match tight {
            Some(tight) => PackedEntries::Tight(tight),
            None => PackedEntries::Compact(compact),
        };
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

    /// **A record whose read alleles fit four bytes holds them so, and one that does not holds
    /// them whole; both give back what went in**: fragments up to 65,535 past a sample's first and
    /// alleles up to 255 fit; a sample whose fragments span one more, or an allele of 256, makes the
    /// record whole. The bytes counted follow the form.
    #[test]
    fn read_alleles_are_held_in_four_bytes_where_they_fit() {
        let entry = |sample: u32, fragment: u32, allele: u16| FragmentAllele {
            sample,
            fragment,
            allele,
            error_phred: 33,
        };
        let check = |entries: Vec<FragmentAllele>, compact: bool| {
            let fragments = RecordFragments::from_entries(entries);
            let packed = PackedRecordFragments::from(fragments.clone());
            assert_eq!(packed.is_compact(), compact);
            for sample in 0..4 {
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
            let samples = 3;
            let entry_bytes = if packed.is_tight() {
                3
            } else if compact {
                4
            } else {
                7
            };
            assert_eq!(
                packed.heap_bytes(),
                samples * 12 + packed.len() * entry_bytes
            );
        };
        let far = 4_000_000_000;
        check(
            vec![
                entry(0, 10, 0),
                entry(0, 10 + 65_535, 255),
                entry(2, far, 1),
                entry(3, 7, 2),
            ],
            true,
        );
        check(
            vec![
                entry(0, 10, 0),
                entry(0, 10 + 65_536, 1),
                entry(2, far, 1),
                entry(3, 7, 2),
            ],
            false,
        );
        check(
            vec![entry(0, 10, 0), entry(2, far, 256), entry(3, 7, 2)],
            false,
        );
        assert_eq!(std::mem::size_of::<CompactReadAllele>(), 4);
        assert_eq!(std::mem::size_of::<TightReadAllele>(), 3);
    }

    /// **Three bytes a read allele where every allele is at most 3 and every Phred at most 63**,
    /// four where one is not, and every entry read back as it was.
    #[test]
    fn read_alleles_are_held_in_three_bytes_where_their_alleles_and_phreds_fit() {
        let entry = |sample: u32, fragment: u32, allele: u16, error_phred: u8| FragmentAllele {
            sample,
            fragment,
            allele,
            error_phred,
        };
        let check = |entries: Vec<FragmentAllele>, tight: bool| {
            let fragments = RecordFragments::from_entries(entries);
            let packed = PackedRecordFragments::from(fragments.clone());
            assert_eq!(packed.is_tight(), tight);
            assert!(packed.is_compact());
            for sample in 0..3 {
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
            let entry_bytes = if tight { 3 } else { 4 };
            assert_eq!(packed.heap_bytes(), 2 * 12 + packed.len() * entry_bytes);
        };
        check(
            vec![
                entry(0, 5, 0, 0),
                entry(0, 5 + 65_535, 3, 63),
                entry(2, 9, 1, 60),
            ],
            true,
        );
        check(vec![entry(0, 5, 4, 30), entry(2, 9, 1, 60)], false);
        check(vec![entry(0, 5, 1, 64), entry(2, 9, 1, 60)], false);
    }
}
