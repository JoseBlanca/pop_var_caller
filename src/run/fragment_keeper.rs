//! **The fragment keeper** — on the ordered thread, keep only the read alleles of reads that reach
//! a second record, and number those reads as per-sample fragments (read phase spec §2.4).
//!
//! Most read alleles belong to reads seen at one record, and a read seen once says nothing about
//! phase. Finding the others needs the records in genome order and memory across them, so this
//! runs where every driver puts records back in order: after the edge-deletion gathering releases
//! each record, before it is prepared for the VCF or the spill. Workers produce read alleles per
//! locus and keep nothing across loci, so neither the building regions nor the thread count can
//! change what is kept here.

use std::collections::VecDeque;

use ahash::AHashMap;

use crate::read_phase::capture::{READ_ALLELES_KEPT, ReadAllele};
use crate::read_phase::fragments::{FragmentAllele, READ_REACH, RecordFragments};
use crate::types::{ChainId, ContigId};
use crate::vcf::VcfRecord;

/// What the keeper did with the read alleles it was handed, for the run's report.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct KeeperCounts {
    /// Read alleles the cap of [`READ_ALLELES_KEPT`] a sample a record removed at capture, before
    /// the keeper: what the cap did, at the deep end of the depth range.
    pub read_alleles_over_the_cap: u64,
    /// Read alleles handed in, over every record and sample.
    pub read_alleles_in: u64,
    /// Kept: their read reached another record within [`READ_REACH`] of its first.
    pub read_alleles_kept: u64,
    /// Dropped because their read had no other read allele within reach: a read seen once.
    pub read_alleles_alone: u64,
    /// Dropped because they lay more than [`READ_REACH`] bases past their read's first one.
    pub read_alleles_beyond_reach: u64,
}

/// **The keeper's counts in the run log**, once every record has left it.
pub fn log_counts(counts: KeeperCounts) {
    tracing::info!(
        "read-backed phase: {} read alleles over the cap of {} a sample a record dropped at \
         capture; {} of {} read alleles kept, their reads seen at two records or more within {} \
         bases; {} dropped as their read's only one, {} as more than {} bases past their read's \
         first",
        counts.read_alleles_over_the_cap,
        READ_ALLELES_KEPT,
        counts.read_alleles_kept,
        counts.read_alleles_in,
        READ_REACH,
        counts.read_alleles_alone,
        counts.read_alleles_beyond_reach,
        READ_REACH,
    );
}

/// What the keeper remembers of one read.
#[derive(Debug, Clone, Copy)]
struct ReadSeen {
    /// The written position of the record its first read allele was at.
    first_position: u64,
    /// Records within [`READ_REACH`] of `first_position` it has a read allele at, the first
    /// included.
    records_within_reach: u32,
    fragment: u32,
}

/// One read in the order the keeper first saw it.
#[derive(Debug, Clone, Copy)]
struct FirstSeen {
    position: u64,
    read: ChainId,
}

/// One sample's reads seen on the current contig.
#[derive(Debug, Default)]
struct SampleReads {
    by_read: AHashMap<ChainId, ReadSeen>,
    /// The same reads in the order they were first seen, which is the order of their
    /// `first_position`, so forgetting pops from the front rather than scanning the map.
    first_seen: VecDeque<FirstSeen>,
    next_fragment: u32,
}

/// One read allele of a held record, still within reach of its read's first.
#[derive(Debug, Clone, Copy)]
struct HeldReadAllele {
    sample: u32,
    read: ChainId,
    allele: u16,
    error_phred: u8,
}

/// One record waiting for the reads after it, with its read alleles taken off it.
#[derive(Debug)]
struct HeldRecord<P> {
    record: VcfRecord,
    payload: P,
    /// The record's written position, the writer's order key.
    position: u64,
    read_alleles: Vec<HeldReadAllele>,
}

/// **The fragment keeper** (§2.4). Records go in in the writer's order; each comes out once the
/// newest record is more than [`READ_REACH`] bases past it, or its contig ends, in the order it
/// went in, with its read alleles taken off it and handed back as [`RecordFragments`].
///
/// **Positions are written positions** — `POS` after the padding rule moved it — because that is
/// what the writer orders records by. A record's region start is not: a generic record padded on
/// the left starts one base after its written `POS`, so it can precede an STR record that starts
/// at that `POS` and is written there too, the one tie the writer admits.
///
/// **The rule.** A read allele is kept when its read has a read allele at another record no more
/// than `R` bases from the read's first one. One more than `R` from it is dropped and counted.
/// The read cap was applied before this, at capture; the order matters, because capping first
/// and testing reach second keeps different reads than the reverse.
///
/// **What it forgets, and why that loses evidence but never invents it.** A read is forgotten once
/// its first record lies more than `R` behind the oldest record still held: none of its read
/// alleles can be kept after that. A read met again later — a pair with an insert over a
/// kilobase, say — is numbered as a new fragment. Its read alleles there are still one molecule's,
/// so the worst this does is make one fragment two.
///
/// `P` is whatever travels beside each record — the driver's coverage windows.
///
/// **A run whose output reads no fragments uses [`Self::passing_through`]**, which keeps none:
/// see that constructor for what it does instead.
#[derive(Debug)]
pub struct FragmentKeeper<P> {
    /// Whether it keeps fragments at all; `false` for [`Self::passing_through`].
    keeps_fragments: bool,
    held: VecDeque<HeldRecord<P>>,
    samples: Vec<SampleReads>,
    contig: Option<ContigId>,
    counts: KeeperCounts,
}

impl<P> Default for FragmentKeeper<P> {
    fn default() -> Self {
        Self {
            keeps_fragments: true,
            held: VecDeque::new(),
            samples: Vec::new(),
            contig: None,
            counts: KeeperCounts::default(),
        }
    }
}

/// One record out of the keeper: the record, what travelled beside it, and its fragments.
pub type KeptRecord<P> = (VcfRecord, P, RecordFragments);

impl<P> FragmentKeeper<P> {
    /// **A keeper for a run whose output reads no fragments**
    /// ([`writer_reads_fragments`](crate::vcf::writer_reads_fragments) is `false`): each record
    /// leaves the moment it arrives, with its read alleles taken off it as the keeping keeper takes
    /// them, and with no fragments.
    ///
    /// **Why the record is not held.** Holding is what lets a later record's reads join this one's
    /// into a fragment; a run that will read no fragment gains nothing from it and pays for it in
    /// records parked on the one thread that puts them in genome order. Released at once, the
    /// records leave in the order they came, as the keeping keeper's do.
    ///
    /// **Why the read alleles are still taken off.** A run that gathers no fragments asks the merge
    /// for no read alleles either, so there are normally none to take; taking them anyway means a
    /// record leaves this keeper in one shape whichever way it was built. Its counts stay at zero,
    /// and the driver logs none.
    #[must_use]
    pub fn passing_through() -> Self {
        Self {
            keeps_fragments: false,
            ..Self::default()
        }
    }

    /// The keeper a run uses: [`Self::default`] where its output reads fragments,
    /// [`Self::passing_through`] where it does not.
    #[must_use]
    pub fn for_a_run(gathers_fragments: bool) -> Self {
        if gathers_fragments {
            Self::default()
        } else {
            Self::passing_through()
        }
    }

    /// Whether this keeper keeps fragments: `false` only for [`Self::passing_through`].
    #[must_use]
    pub fn keeps_fragments(&self) -> bool {
        self.keeps_fragments
    }

    /// **Take one record**, appending to `released` every held record it releases, in the order
    /// they came.
    ///
    /// # Panics
    ///
    /// If a record is written before the one handed in before it on the same contig: the writer's
    /// order is non-decreasing in written position, so that is a wiring defect.
    pub fn accept(&mut self, mut record: VcfRecord, payload: P, released: &mut Vec<KeptRecord<P>>) {
        if !self.keeps_fragments {
            let _ = record.replace_read_alleles(Default::default());
            released.push((record, payload, RecordFragments::default()));
            return;
        }
        let contig = record.region().contig;
        let position = crate::vcf::encode::written_position(&record);
        if self.contig != Some(contig) {
            self.release_all_into(released);
            for sample in &mut self.samples {
                sample.by_read.clear();
                sample.first_seen.clear();
                sample.next_fragment = 0;
            }
            self.contig = Some(contig);
        }
        if let Some(newest) = self.held.back() {
            assert!(
                newest.position <= position,
                "records out of order: {} after {}",
                position,
                newest.position
            );
        }
        while let Some(oldest) = self.held.front()
            && position - oldest.position > READ_REACH
        {
            // PANIC-FREE: the loop condition just saw a front.
            let oldest = self.held.pop_front().expect("a front, just seen");
            released.push(self.release(oldest));
        }
        self.forget_reads_out_of_reach(
            self.held.front().map_or(position, |oldest| oldest.position),
        );

        let read_alleles = record.replace_read_alleles(Default::default());
        self.counts.read_alleles_over_the_cap += read_alleles.dropped_by_the_cap();
        let mut within_reach = Vec::with_capacity(read_alleles.len());
        for entry in read_alleles.entries() {
            self.counts.read_alleles_in += 1;
            let sample = entry.sample as usize;
            if self.samples.len() <= sample {
                self.samples.resize_with(sample + 1, SampleReads::default);
            }
            let reads = &mut self.samples[sample];
            let ReadAllele {
                read,
                allele,
                error_phred,
            } = entry.read_allele;
            let next_fragment = &mut reads.next_fragment;
            let first_seen = &mut reads.first_seen;
            let seen = reads.by_read.entry(read).or_insert_with(|| {
                let fragment = *next_fragment;
                *next_fragment += 1;
                first_seen.push_back(FirstSeen { position, read });
                ReadSeen {
                    first_position: position,
                    records_within_reach: 0,
                    fragment,
                }
            });
            if position - seen.first_position > READ_REACH {
                self.counts.read_alleles_beyond_reach += 1;
                continue;
            }
            seen.records_within_reach += 1;
            within_reach.push(HeldReadAllele {
                sample: entry.sample,
                read,
                // PANIC-FREE: the allele indexes the record's alleles, which calling numbers as
                // `AlleleId`, a `u16`.
                allele: u16::try_from(allele).expect("a record holds fewer than 65,536 alleles"),
                error_phred,
            });
        }
        self.held.push_back(HeldRecord {
            record,
            payload,
            position,
            read_alleles: within_reach,
        });
    }

    /// **Release everything held**, at the end of the run, into `released`.
    pub fn finish(&mut self, released: &mut Vec<KeptRecord<P>>) {
        self.release_all_into(released);
    }

    /// What the keeper has done so far.
    #[must_use]
    pub fn counts(&self) -> KeeperCounts {
        self.counts
    }

    fn release_all_into(&mut self, released: &mut Vec<KeptRecord<P>>) {
        while let Some(oldest) = self.held.pop_front() {
            released.push(self.release(oldest));
        }
    }

    /// One record out: every record within reach of each of its reads' first ones has arrived,
    /// since the newest is more than `R` past this one, which is at or after the read's first.
    fn release(&mut self, held: HeldRecord<P>) -> KeptRecord<P> {
        let mut entries = Vec::new();
        for HeldReadAllele {
            sample,
            read,
            allele,
            error_phred,
        } in held.read_alleles
        {
            // A held read allele was kept within reach, `position − first_position ≤ R`, and a
            // read is forgotten only once `first_position + R` is behind the oldest held record,
            // which is this one or later; so a held record's reads are never forgotten.
            let Some(seen) = self.samples[sample as usize].by_read.get(&read) else {
                // UNREACHABLE: by the forgetting rule just stated.
                unreachable!("a held record's read {read} of sample {sample} was forgotten");
            };
            if seen.records_within_reach >= 2 {
                self.counts.read_alleles_kept += 1;
                entries.push(FragmentAllele {
                    sample,
                    fragment: seen.fragment,
                    allele,
                    error_phred,
                });
            } else {
                self.counts.read_alleles_alone += 1;
            }
        }
        (
            held.record,
            held.payload,
            RecordFragments::from_entries(entries),
        )
    }

    /// Forget the reads none of whose read alleles can still be kept: their first record lies more
    /// than `R` behind `oldest_position`, the oldest record held (or the one arriving, when none
    /// is), so every read allele they have at a held record was already dropped as beyond reach,
    /// and every later one would be too.
    fn forget_reads_out_of_reach(&mut self, oldest_position: u64) {
        for sample in &mut self.samples {
            while let Some(&FirstSeen { position, read }) = sample.first_seen.front()
                && position + READ_REACH < oldest_position
            {
                sample.first_seen.pop_front();
                sample.by_read.remove(&read);
            }
        }
    }
}

#[cfg(test)]
mod tests;
