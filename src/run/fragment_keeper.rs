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
use rayon::prelude::*;

use crate::read_phase::capture::{READ_ALLELES_KEPT, ReadAllele};
use crate::read_phase::fragments::{FragmentAllele, READ_REACH, RecordFragments};
use crate::read_phase::record_reads::SampleReadAllele;
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

/// **How many samples one task of the keeper handles.** Small enough that a cohort of a few
/// thousand splits into dozens of tasks for the pool, large enough that a task's work — tens of
/// microseconds at three reads a position — outweighs handing it to another thread.
pub const SAMPLES_A_TASK: usize = 64;

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
/// **Each sample's reads are its own, so the samples are shared among the threads.** Taking a
/// record in, forgetting reads out of reach and releasing a record each touch one sample's reads
/// at a time and nothing else, so each is done over groups of [`SAMPLES_A_TASK`] samples at once
/// and the groups' results joined in sample order: the result is the one-thread result, read
/// allele for read allele, whatever the thread count. It matters at the large end of the cohort
/// range: at 2,176 tomato samples and about three reads a position a record carried about 18,000
/// read alleles, and the keeper on one thread cost 1.8 ms a record, a third of the time the one
/// thread putting records in genome order spent on each. Below `SAMPLES_A_TASK` samples there is one
/// group and the work stays on the calling thread.
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
    /// they came. The same as [`Self::accept_batch`] with a batch of one.
    ///
    /// # Panics
    ///
    /// If a record is written before the one handed in before it on the same contig: the writer's
    /// order is non-decreasing in written position, so that is a wiring defect.
    pub fn accept(&mut self, record: VcfRecord, payload: P, released: &mut Vec<KeptRecord<P>>) {
        self.accept_batch(vec![(record, payload)], released);
    }

    /// **Take several records, in order**, appending to `released` every held record they
    /// release, in the order they came: the same records, read alleles and counts as taking them
    /// one at a time.
    ///
    /// **Why a batch.** Which held records a record releases, and how far back reads are forgotten,
    /// depend on record positions alone, so they are planned here on one thread for the whole
    /// batch. Each group of [`SAMPLES_A_TASK`] samples then replays the plan in order on its own
    /// reads, and the groups' results are joined in sample order. Taken one at a time, each record
    /// paid up to three rounds of handing its groups to the pool and waiting for them; at 2,176
    /// tomato samples that wait, not the work, was most of the keeper's time.
    ///
    /// # Panics
    ///
    /// As [`Self::accept`].
    pub fn accept_batch(&mut self, batch: Vec<(VcfRecord, P)>, released: &mut Vec<KeptRecord<P>>) {
        if !self.keeps_fragments {
            for (mut record, payload) in batch {
                let _ = record.replace_read_alleles(Default::default());
                released.push((record, payload, RecordFragments::default()));
            }
            return;
        }
        // **The plan, from positions alone**, in the order taking the records one at a time
        // performs it. `combined` numbers the held records as they will stand: those held
        // before the batch first, then the batch's own, so releasing is always the next number.
        let already_held = self.held.len();
        let mut queue: VecDeque<u64> = self.held.iter().map(|held| held.position).collect();
        let mut next_released = 0_usize;
        let mut newest = self.held.back().map(|held| held.position);
        let mut contig = self.contig;
        let mut steps = Vec::with_capacity(batch.len() * 3);
        let mut incoming = Vec::with_capacity(batch.len());
        for (mut record, payload) in batch {
            let record_contig = record.region().contig;
            let position = crate::vcf::encode::written_position(&record);
            if contig != Some(record_contig) {
                while queue.pop_front().is_some() {
                    steps.push(Step::Release(next_released));
                    next_released += 1;
                }
                steps.push(Step::Reset);
                contig = Some(record_contig);
                newest = None;
            }
            if let Some(newest) = newest {
                assert!(
                    newest <= position,
                    "records out of order: {position} after {newest}"
                );
            }
            while let Some(&oldest) = queue.front()
                && position - oldest > READ_REACH
            {
                queue.pop_front();
                steps.push(Step::Release(next_released));
                next_released += 1;
            }
            steps.push(Step::Forget(queue.front().copied().unwrap_or(position)));
            steps.push(Step::TakeIn(incoming.len(), position));
            queue.push_back(position);
            newest = Some(position);
            let read_alleles = record.replace_read_alleles(Default::default());
            self.counts.read_alleles_over_the_cap += read_alleles.dropped_by_the_cap();
            // The entries are sorted by sample, so the last names the highest.
            if let Some(last) = read_alleles.entries().last()
                && self.samples.len() <= last.sample as usize
            {
                self.samples
                    .resize_with(last.sample as usize + 1, SampleReads::default);
            }
            incoming.push((record, payload, position, read_alleles));
        }
        self.contig = contig;

        // **Each group of samples replays the plan.**
        // Only the read alleles go to the pool, not the records or what travels beside them.
        let held_read_alleles: Vec<&[HeldReadAllele]> = self
            .held
            .iter()
            .take(next_released.min(already_held))
            .map(|held| held.read_alleles.as_slice())
            .collect();
        let incoming_entries: Vec<&[SampleReadAllele]> = incoming
            .iter()
            .map(|(.., read_alleles)| read_alleles.entries())
            .collect();
        let held = &held_read_alleles;
        let incoming_ref = &incoming_entries;
        let steps = &steps;
        let samples = &mut self.samples;
        let per_task: Vec<TaskOutcome> = samples
            .par_chunks_mut(SAMPLES_A_TASK)
            .enumerate()
            .map(|(task, samples)| {
                let first = task * SAMPLES_A_TASK;
                let mut outcome = TaskOutcome {
                    taken: Vec::with_capacity(incoming_ref.len()),
                    released: Vec::new(),
                };
                for step in steps {
                    match *step {
                        Step::Reset => {
                            for sample in samples.iter_mut() {
                                sample.by_read.clear();
                                sample.first_seen.clear();
                                sample.next_fragment = 0;
                            }
                        }
                        Step::Forget(oldest_position) => {
                            forget_out_of_reach(samples, oldest_position);
                        }
                        Step::TakeIn(at, position) => {
                            let entries = incoming_ref[at];
                            let in_task =
                                of_samples(entries, first, samples.len(), |entry| entry.sample);
                            outcome
                                .taken
                                .push(take_in(samples, first, in_task, position));
                        }
                        Step::Release(number) => {
                            let read_alleles: &[HeldReadAllele] = if number < already_held {
                                of_samples(held[number], first, samples.len(), |held| held.sample)
                            } else {
                                &outcome.taken[number - already_held].0
                            };
                            let fragments = fragments_of(samples, first, read_alleles);
                            outcome.released.push(fragments);
                        }
                    }
                }
                outcome
            })
            .collect();

        // **Joined in sample order**: the batch's records go on the back of the held list, and
        // the releases come off its front, in the plan's order.
        let mut per_task = per_task;
        for (at, (record, payload, position, read_alleles)) in incoming.into_iter().enumerate() {
            let entries = read_alleles.entries();
            self.counts.read_alleles_in += entries.len() as u64;
            let mut within_reach = Vec::with_capacity(entries.len());
            for task in &mut per_task {
                let (kept, beyond_reach) = std::mem::take(&mut task.taken[at]);
                within_reach.extend(kept);
                self.counts.read_alleles_beyond_reach += beyond_reach;
            }
            self.held.push_back(HeldRecord {
                record,
                payload,
                position,
                read_alleles: within_reach,
            });
        }
        for release in 0..next_released {
            // PANIC-FREE: the plan released no more records than were held or taken in.
            let oldest = self
                .held
                .pop_front()
                .expect("a held record the plan released");
            let mut entries = Vec::with_capacity(oldest.read_alleles.len());
            for task in &mut per_task {
                let (kept, alone) = std::mem::take(&mut task.released[release]);
                self.counts.read_alleles_kept += kept.len() as u64;
                self.counts.read_alleles_alone += alone;
                entries.extend(kept);
            }
            released.push((
                oldest.record,
                oldest.payload,
                RecordFragments::from_entries(entries),
            ));
        }
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
        let read_alleles = &held.read_alleles;
        let released: Vec<(Vec<FragmentAllele>, u64)> = self
            .samples
            .par_chunks(SAMPLES_A_TASK)
            .enumerate()
            .map(|(task, samples)| {
                let first = task * SAMPLES_A_TASK;
                let in_task = of_samples(read_alleles, first, samples.len(), |held| held.sample);
                fragments_of(samples, first, in_task)
            })
            .collect();
        let mut entries = Vec::with_capacity(read_alleles.len());
        for (kept, alone) in released {
            self.counts.read_alleles_kept += kept.len() as u64;
            self.counts.read_alleles_alone += alone;
            entries.extend(kept);
        }
        (
            held.record,
            held.payload,
            RecordFragments::from_entries(entries),
        )
    }
}

/// One step of a batch's plan ([`FragmentKeeper::accept_batch`]), in the order taking the
/// records one at a time performs it.
#[derive(Debug, Clone, Copy)]
enum Step {
    /// Release held record number `n`: those held before the batch first, then the batch's own.
    Release(usize),
    /// A new contig: forget every read.
    Reset,
    /// Forget the reads whose first record lies more than `R` behind this position.
    Forget(u64),
    /// Take in the batch's record `at`, written at `position`.
    TakeIn(usize, u64),
}

/// What one group of samples made of a batch's plan: per record taken in, its read alleles
/// within reach and the count beyond; per record released, its fragment alleles and the count
/// alone.
struct TaskOutcome {
    taken: Vec<(Vec<HeldReadAllele>, u64)>,
    released: Vec<(Vec<FragmentAllele>, u64)>,
}

/// Forget the reads none of whose read alleles can still be kept: their first record lies more
/// than `R` behind `oldest_position`, the oldest record held (or the one arriving, when none
/// is), so every read allele they have at a held record was already dropped as beyond reach, and
/// every later one would be too.
fn forget_out_of_reach(samples: &mut [SampleReads], oldest_position: u64) {
    for sample in samples {
        while let Some(&FirstSeen { position, read }) = sample.first_seen.front()
            && position + READ_REACH < oldest_position
        {
            sample.first_seen.pop_front();
            sample.by_read.remove(&read);
        }
    }
}

/// **The part of `sorted` — sorted by sample — whose samples are `first..first + count`**, the
/// samples one task of the keeper holds.
fn of_samples<T>(sorted: &[T], first: usize, count: usize, sample: impl Fn(&T) -> u32) -> &[T] {
    let start = sorted.partition_point(|item| (sample(item) as usize) < first);
    let end =
        start + sorted[start..].partition_point(|item| (sample(item) as usize) < first + count);
    &sorted[start..end]
}

/// **One task's samples take in a record's read alleles at `position`**: `samples[i]` is run
/// sample `first + i`, and `entries` are those samples' read alleles, sorted by sample. Returns the
/// read alleles within reach of their read's first, in the order given, and how many were not.
fn take_in(
    samples: &mut [SampleReads],
    first: usize,
    entries: &[SampleReadAllele],
    position: u64,
) -> (Vec<HeldReadAllele>, u64) {
    let mut within_reach = Vec::with_capacity(entries.len());
    let mut beyond_reach = 0;
    for entry in entries {
        let reads = &mut samples[entry.sample as usize - first];
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
            beyond_reach += 1;
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
    (within_reach, beyond_reach)
}

/// **One task's samples release a held record's read alleles**: `samples[i]` is run sample
/// `first + i`, and `held` are those samples' read alleles. Returns the fragment alleles of reads
/// seen at two records or more, and how many read alleles were their read's only one.
fn fragments_of(
    samples: &[SampleReads],
    first: usize,
    held: &[HeldReadAllele],
) -> (Vec<FragmentAllele>, u64) {
    let mut kept = Vec::with_capacity(held.len());
    let mut alone = 0;
    for &HeldReadAllele {
        sample,
        read,
        allele,
        error_phred,
    } in held
    {
        // A held read allele was kept within reach, `position − first_position ≤ R`, and a read is
        // forgotten only once `first_position + R` is behind the oldest held record, which is this
        // one or later; so a held record's reads are never forgotten.
        let Some(seen) = samples[sample as usize - first].by_read.get(&read) else {
            // UNREACHABLE: by the forgetting rule just stated.
            unreachable!("a held record's read {read} of sample {sample} was forgotten");
        };
        if seen.records_within_reach >= 2 {
            kept.push(FragmentAllele {
                sample,
                fragment: seen.fragment,
                allele,
                error_phred,
            });
        } else {
            alone += 1;
        }
    }
    (kept, alone)
}

#[cfg(test)]
mod tests;
