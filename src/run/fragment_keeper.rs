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

use crate::read_phase::capture::ReadAllele;
use crate::read_phase::fragments::{FragmentAllele, READ_REACH, RecordFragments};
use crate::types::{ChainId, ContigId};
use crate::vcf::VcfRecord;

/// What the keeper did with the read alleles it was handed, for the run's report.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct KeeperCounts {
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
        "read-backed phase: {} of {} read alleles kept, their reads seen at two records or more \
         within {} bases; {} dropped as their read's only one, {} as more than {} bases past \
         their read's first",
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
    /// The position of the record its first read allele was at.
    first_at: u64,
    /// Records within [`READ_REACH`] of `first_at` it has a read allele at, the first included.
    records_within_reach: u32,
    fragment: u32,
}

/// One sample's reads seen on the current contig.
#[derive(Debug, Default)]
struct SampleReads {
    by_read: AHashMap<ChainId, ReadSeen>,
    /// The same reads in the order they were first seen, which is the order of their `first_at`,
    /// so forgetting pops from the front rather than scanning the map.
    first_seen: VecDeque<(u64, ChainId)>,
    next_fragment: u32,
}

/// One record waiting for the reads after it, with its read alleles taken off it.
#[derive(Debug)]
struct Held<P> {
    record: VcfRecord,
    payload: P,
    at: u64,
    /// The record's read alleles still within reach of their read's first, as `(sample, read,
    /// allele, Phred)`.
    read_alleles: Vec<(u32, ChainId, u16, u8)>,
}

/// **The fragment keeper** (§2.4). Records go in in the writer's order; each comes out once the
/// newest record is more than [`READ_REACH`] bases past it, or its contig ends, in the order it
/// went in, with its read alleles taken off it and handed back as [`RecordFragments`].
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
#[derive(Debug)]
pub struct FragmentKeeper<P> {
    held: VecDeque<Held<P>>,
    samples: Vec<SampleReads>,
    contig: Option<ContigId>,
    counts: KeeperCounts,
}

impl<P> Default for FragmentKeeper<P> {
    fn default() -> Self {
        Self {
            held: VecDeque::new(),
            samples: Vec::new(),
            contig: None,
            counts: KeeperCounts::default(),
        }
    }
}

impl<P> FragmentKeeper<P> {
    /// **Take one record**, handing back every held record it releases, in the order they came.
    ///
    /// # Panics
    ///
    /// If a record starts before the one handed in before it on the same contig: the writer's
    /// order is non-decreasing, so that is a wiring defect.
    pub fn accept(
        &mut self,
        mut record: VcfRecord,
        payload: P,
    ) -> Vec<(VcfRecord, P, RecordFragments)> {
        let region = record.region();
        let at = region.start.get();
        let mut released = Vec::new();
        if self.contig != Some(region.contig) {
            self.release_all_into(&mut released);
            for sample in &mut self.samples {
                sample.by_read.clear();
                sample.first_seen.clear();
                sample.next_fragment = 0;
            }
            self.contig = Some(region.contig);
        }
        if let Some(newest) = self.held.back() {
            assert!(
                newest.at <= at,
                "records out of order: {} after {}",
                at,
                newest.at
            );
        }
        while self
            .held
            .front()
            .is_some_and(|oldest| at - oldest.at > READ_REACH)
        {
            let oldest = self.held.pop_front().expect("a front, just seen");
            released.push(self.release(oldest));
        }
        self.forget_reads_out_of_reach(self.held.front().map_or(at, |oldest| oldest.at));

        let read_alleles = record.replace_read_alleles(Default::default());
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
                first_seen.push_back((at, read));
                ReadSeen {
                    first_at: at,
                    records_within_reach: 0,
                    fragment,
                }
            });
            if at - seen.first_at > READ_REACH {
                self.counts.read_alleles_beyond_reach += 1;
                continue;
            }
            seen.records_within_reach += 1;
            within_reach.push((
                entry.sample,
                read,
                u16::try_from(allele).expect("a record holds fewer than 65,536 alleles"),
                error_phred,
            ));
        }
        self.held.push_back(Held {
            record,
            payload,
            at,
            read_alleles: within_reach,
        });
        released
    }

    /// **Release everything held**, at the end of the run.
    pub fn finish(&mut self) -> Vec<(VcfRecord, P, RecordFragments)> {
        let mut released = Vec::new();
        self.release_all_into(&mut released);
        released
    }

    /// What the keeper has done so far.
    #[must_use]
    pub fn counts(&self) -> KeeperCounts {
        self.counts
    }

    fn release_all_into(&mut self, released: &mut Vec<(VcfRecord, P, RecordFragments)>) {
        while let Some(oldest) = self.held.pop_front() {
            released.push(self.release(oldest));
        }
    }

    /// One record out: every record within reach of each of its reads' first ones has arrived,
    /// since the newest is more than `R` past this one, which is at or after the read's first.
    fn release(&mut self, held: Held<P>) -> (VcfRecord, P, RecordFragments) {
        let mut entries = Vec::new();
        for (sample, read, allele, error_phred) in held.read_alleles {
            let seen = self.samples[sample as usize].by_read[&read];
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
    /// than `R` behind `oldest`, the oldest record held (or the one arriving, when none is), so
    /// every read allele they have at a held record was already dropped as beyond reach, and every
    /// later one would be too.
    fn forget_reads_out_of_reach(&mut self, oldest: u64) {
        for sample in &mut self.samples {
            while let Some(&(first_at, read)) = sample.first_seen.front()
                && first_at + READ_REACH < oldest
            {
                sample.first_seen.pop_front();
                sample.by_read.remove(&read);
            }
        }
    }
}

#[cfg(test)]
mod tests;
