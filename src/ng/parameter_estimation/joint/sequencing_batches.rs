//! **Who was sequenced beside whom, as the run was told** — the grouping the contaminating
//! population is drawn from.
//!
//! A **sequencing batch** is the set of libraries that ran together: a flowcell, a plate, a
//! submission. It matters for one reason. A contaminating read is far likelier to have come
//! from a neighbour on the same run than from a random member of the species, so when the read
//! likelihood asks *how often does the contaminating population show this allele here*, the
//! population it means is the samples that ran beside this one and not the whole cohort
//! (`doc/devel/ng/spec/read_likelihoods.md` §3.6;
//! `doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
//!
//! **Nothing declares one yet.** [`SequencingBatches::declared`] has no caller outside this
//! module's tests: the command-line flag that would carry a batching is unbuilt, so every run
//! today gets [`SequencingBatches::all_together`] and every contaminated read group is scored
//! against the whole cohort's frequency. That is the honest default and costs such a run
//! nothing; what it means for a reader is that the refusals below describe a shape no input can
//! currently reach.
//!
//! **It is stated, never inferred.** The grouping is absent from both benchmark cohorts'
//! alignments — the tomato archive's `@RG` lines carry no platform unit, and SRA rewrote the
//! read names — so a pipeline that guessed it from what survives would be wrong in silence. The
//! default is [`SequencingBatches::all_together`], one batch holding the whole run, which is the
//! honest statement of what a run knows when nobody has said otherwise: every read group then
//! reads the cohort frequency, and a run that declares no batching loses nothing it had.
//!
//! # Two views of one partition, and the loop needs both
//!
//! The batching is declared over **read groups**, because that is the grain a file header gives
//! and because one sample's libraries can in principle run on different flowcells. The calling
//! loop reads it two ways:
//!
//! - [`BatchOfEachReadGroup`] — which row of the contaminant-frequency table each *library's*
//!   reads are scored against, since the contamination fraction beside it is per read group;
//! - [`BatchOfEachSample`] — which batch each *sample's* expected allele copies are added into
//!   when those frequencies are built.
//!
//! They are two wrapper types over one slice shape for the reason
//! [`BatchOfEachReadGroup`](crate::ng::types::BatchOfEachReadGroup) records: at one library per
//! sample the two agree in length, so a transposition passes every shape check and comes back as
//! a wrong contaminant frequency.
//!
//! # A sample whose libraries ran in different batches is refused
//!
//! The sample-keyed view needs one batch per sample, and a sample split across two batches has
//! no single answer: the choices are to pick one, to average the two populations, or to refuse.
//! **This refuses** ([`SequencingBatchError::SampleSequencedInSeveralBatches`]), because the
//! read likelihood deliberately declined to invent a rule here and because picking a majority
//! batch would score half a sample's reads against the wrong neighbours with nothing said.
//! **Under the shipped default it cannot arise** — one batch holds everything — so what this
//! refusal costs is a run that declares a batching splitting a sample, and what it buys is that
//! nobody discovers the rule by reading the genotypes.

use std::collections::{BTreeMap, BTreeSet};

use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::types::{BatchId, BatchOfEachReadGroup, BatchOfEachSample, ReadGroupId};

/// Why a declared batching does not describe this run.
///
/// Every variant is a refusal rather than a repair, and each names the thing the run would have
/// got instead: a wrong contaminant population for some subset of the reads, which changes
/// genotypes and appears nowhere in the output.
#[non_exhaustive]
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum SequencingBatchError {
    /// Read groups the declared batching left out.
    ///
    /// **A run that names any batch must name them all.** A user who lists three plates and
    /// forgets four libraries would otherwise get whatever the unlisted ones defaulted to —
    /// the cohort frequency, or the first batch's — for exactly those libraries, and nothing
    /// would say so (`doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
    #[error(
        "the declared sequencing batches leave read group(s) {read_groups:?} out; a run that \
         names any batch must name them all, because an unlisted library would be scored \
         against a contaminating population nobody chose for it"
    )]
    ReadGroupNotBatched { read_groups: Vec<u32> },

    /// One read group named by two batches.
    ///
    /// A batching is a partition, so this is not a wider claim but a contradictory one: the
    /// library cannot have run on both plates.
    #[error(
        "read group {read_group} is named by sequencing batches {first} and {second}; a \
         library ran on one of them, so this batching is not a partition of the run"
    )]
    ReadGroupInSeveralBatches {
        read_group: u32,
        first: usize,
        second: usize,
    },

    /// A batch naming a read group the run does not have.
    ///
    /// The batching and the run's read groups were minted from different inputs — the id means
    /// nothing here, and silently dropping it would shrink a batch by an unknown amount.
    #[error(
        "sequencing batch {batch} names read group {read_group}, and the run has {read_groups} \
         read groups; the batching and the run's read groups came from different inputs"
    )]
    UnknownReadGroup {
        batch: usize,
        read_group: u32,
        read_groups: usize,
    },

    /// A declared batch holding no read groups.
    ///
    /// It would become a row of the contaminant-frequency table that no library reads and no
    /// sample fills — and a row nobody filled leaves through the frequency's no-evidence
    /// fallback, which is indistinguishable from a batch that was really sequenced and really
    /// showed nothing.
    #[error(
        "sequencing batch {batch} holds no read groups; a batch nobody was sequenced in is a \
         batching that does not describe this run"
    )]
    EmptyBatch { batch: usize },

    /// A batching that declares no batches at all.
    ///
    /// *No batching was declared* is [`SequencingBatches::all_together`] and not an empty list —
    /// one named way to say it, so that a caller reaches the decision rather than the shortest
    /// thing that compiles.
    #[error(
        "a declared batching holds no batches at all; a run that was told nothing about who ran \
         beside whom is `SequencingBatches::all_together`, which is one batch holding all of it"
    )]
    NoBatches,

    /// A sample whose libraries were declared in more than one batch.
    ///
    /// See this module's own documentation: the rule for such a sample was deliberately never
    /// invented, and refusing is what keeps it from being decided by accident.
    #[error(
        "sample {sample} was sequenced in batches {batches:?}; the contaminating population a \
         sample's reads are scored against is the batch it ran in, and no rule has been settled \
         for a sample whose libraries ran in several — so this run is refused rather than given \
         one of them"
    )]
    SampleSequencedInSeveralBatches { sample: String, batches: Vec<u32> },
}

/// Which read groups were sequenced together, as the run was told — dense over both axes the
/// calling loop reads.
///
/// **A partition, checked at construction**: every read group of the run in exactly one batch,
/// and the batch ids are `0..batch_count` with nothing missing, because both consumers index a
/// table by them. The default is one batch holding everything, so the type is never optional and
/// no consumer branches on its absence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequencingBatches {
    /// Entry *i* is read group *i*'s batch.
    of_each_read_group: Vec<BatchId>,
    /// Entry *i* is sample *i*'s batch, in the run's sample order — the order
    /// [`ReadGroups::read_groups_per_sample`] mints.
    of_each_sample: Vec<BatchId>,
    batch_count: usize,
    /// Whether this came from [`Self::all_together`] rather than from a declaration.
    defaulted: bool,
}

impl SequencingBatches {
    /// **Every read group in one batch — the default**, and the honest statement of what a run
    /// knows when nobody has said otherwise.
    ///
    /// # Panics
    ///
    /// On a run with no read groups or no samples. Every read of a run belongs to a read group
    /// and every read group names a sample, so an empty axis is a run whose read groups went
    /// missing rather than a run with none.
    #[must_use]
    pub fn all_together(groups: &ReadGroups) -> Self {
        let (read_group_count, sample_count) = checked_axes(groups);
        Self::all_together_over(read_group_count, sample_count)
    }

    /// **One batch holding a run of `libraries` read groups over `samples` samples** — the same
    /// default as [`all_together`](Self::all_together), where the run's read-group table is not
    /// to hand.
    ///
    /// The default needs nothing from that table but its two axis lengths, and there are places
    /// that have the lengths and not the table: a benchmark, and the allocation test, which
    /// links the library from outside and cannot reach a real header.
    ///
    /// **Two counts in a row is a swap waiting to happen, and here it is a harmless one.**
    /// [`crate::ng::calling::FrozenParameters::new`] checks the batching against *both* of the
    /// run's axes — one calibration per read group, one inbreeding coefficient per sample — so a
    /// transposed pair is refused there unless the two counts are equal, and where they are
    /// equal the transposed value is the same value. That is not true of the two *views*, which
    /// carry different meanings at the same length, and is why those stay separate types.
    ///
    /// # Panics
    ///
    /// On a run of no read groups or no samples, as [`all_together`](Self::all_together).
    #[must_use]
    pub fn all_together_over(libraries: usize, samples: usize) -> Self {
        assert!(
            libraries > 0,
            "every read of a run belongs to a read group and a run has at least one, so a run of \
             no read groups is one whose read groups went missing"
        );
        assert!(
            samples > 0,
            "every read group of a run names a sample, so a run with read groups and no samples \
             is one whose sample order went missing"
        );
        Self {
            of_each_read_group: vec![BatchId::ALL_TOGETHER; libraries],
            of_each_sample: vec![BatchId::ALL_TOGETHER; samples],
            batch_count: 1,
            defaulted: true,
        }
    }

    /// The batching the run was given — batch *i* is the *i*-th entry of `batches`.
    ///
    /// # Errors
    ///
    /// Every variant of [`SequencingBatchError`]; each one is a way the declaration fails to
    /// describe this run, and each is refused rather than repaired.
    ///
    /// # Panics
    ///
    /// On a run with no read groups or no samples, as [`Self::all_together`].
    pub fn declared(
        groups: &ReadGroups,
        batches: &[BTreeSet<ReadGroupId>],
    ) -> Result<Self, SequencingBatchError> {
        let (read_group_count, sample_count) = checked_axes(groups);
        if batches.is_empty() {
            return Err(SequencingBatchError::NoBatches);
        }

        let mut of_each_read_group: BTreeMap<ReadGroupId, usize> = BTreeMap::new();
        for (batch, members) in batches.iter().enumerate() {
            if members.is_empty() {
                return Err(SequencingBatchError::EmptyBatch { batch });
            }
            for &read_group in members {
                if read_group.get() as usize >= read_group_count {
                    return Err(SequencingBatchError::UnknownReadGroup {
                        batch,
                        read_group: read_group.get(),
                        read_groups: read_group_count,
                    });
                }
                if let Some(&first) = of_each_read_group.get(&read_group) {
                    return Err(SequencingBatchError::ReadGroupInSeveralBatches {
                        read_group: read_group.get(),
                        first,
                        second: batch,
                    });
                }
                of_each_read_group.insert(read_group, batch);
            }
        }

        // **Every read group of the run, not merely every read group somebody named.** The
        // check above only sees the ids the declaration mentions.
        let unbatched: Vec<u32> = (0..read_group_count)
            .map(|group| ReadGroupId(group as u32))
            .filter(|group| !of_each_read_group.contains_key(group))
            .map(ReadGroupId::get)
            .collect();
        if !unbatched.is_empty() {
            return Err(SequencingBatchError::ReadGroupNotBatched {
                read_groups: unbatched,
            });
        }

        let of_each_read_group: Vec<BatchId> = (0..read_group_count)
            .map(|group| BatchId(of_each_read_group[&ReadGroupId(group as u32)] as u32))
            .collect();

        // **The sample axis, and the one place a declaration can be well-formed and still have
        // no answer.** A sample's libraries must all have run together, because the frequency
        // its contaminant is drawn against is one batch's.
        let mut of_each_sample = Vec::with_capacity(sample_count);
        for sample in groups.read_groups_per_sample() {
            let mut sample_batches: BTreeSet<u32> = BTreeSet::new();
            for &read_group in &sample.read_groups {
                sample_batches.insert(of_each_read_group[read_group.get() as usize].get());
            }
            let mut sample_batches = sample_batches.into_iter();
            let first = sample_batches
                .next()
                .expect("every sample of a run has at least one read group");
            if sample_batches.next().is_some() {
                let mut batches: Vec<u32> = sample
                    .read_groups
                    .iter()
                    .map(|&group| of_each_read_group[group.get() as usize].get())
                    .collect();
                batches.sort_unstable();
                batches.dedup();
                return Err(SequencingBatchError::SampleSequencedInSeveralBatches {
                    sample: sample.sample.to_string(),
                    batches,
                });
            }
            of_each_sample.push(BatchId(first));
        }

        Ok(Self {
            of_each_read_group,
            of_each_sample,
            batch_count: batches.len(),
            defaulted: false,
        })
    }

    /// Which batch each read group ran in, in read-group id order.
    #[inline]
    #[must_use]
    pub fn of_each_read_group(&self) -> BatchOfEachReadGroup<'_> {
        BatchOfEachReadGroup(&self.of_each_read_group)
    }

    /// Which batch each sample ran in, in the run's sample order.
    #[inline]
    #[must_use]
    pub fn of_each_sample(&self) -> BatchOfEachSample<'_> {
        BatchOfEachSample(&self.of_each_sample)
    }

    /// How many batches the run declares — **one** under the default.
    #[inline]
    #[must_use]
    pub fn batch_count(&self) -> usize {
        self.batch_count
    }

    /// How many read groups this batching covers.
    #[inline]
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.of_each_read_group.len()
    }

    /// How many samples this batching covers.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.of_each_sample.len()
    }

    /// Whether this is the default — one batch holding everything, because nobody said
    /// otherwise.
    ///
    /// **It travels with the numbers drawn against it.** Two runs under different batchings
    /// produce contaminant frequencies that are not comparable, and a frequency taken over one
    /// batch holding the whole cohort is the weaker kind — so an output that reports a
    /// contamination fraction has to be able to say which of the two it was
    /// (`doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6). The dense
    /// [`BatchOfEachReadGroup`] cannot answer it: a defaulted batching and a declaration of one
    /// batch holding every library are the same value.
    #[inline]
    #[must_use]
    pub fn is_default(&self) -> bool {
        self.defaulted
    }
}

/// The two axis lengths, refused where either is empty.
fn checked_axes(groups: &ReadGroups) -> (usize, usize) {
    assert!(
        !groups.is_empty(),
        "every read of a run belongs to a read group and a run has at least one, so a run whose \
         read-group table is empty is one whose read groups went missing"
    );
    let sample_count = groups.read_groups_per_sample().len();
    // **Held in debug only, and that is not an oversight.** `ReadGroups` groups by sample, so a
    // non-empty read-group table has at least one sample and the check above has already refused
    // an empty one — no input reaches this, so no test can, and a release check no test can
    // reach is one the suite cannot keep honest. What it guards is `ReadGroups` changing how it
    // builds the by-sample view.
    debug_assert!(
        sample_count > 0,
        "every read group of a run names a sample, so a run with read groups and no samples is \
         one whose sample order went missing"
    );
    (groups.len(), sample_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A run whose read groups are `libraries`, each naming the sample beside it.
    ///
    /// **The sample order is first-seen**, which is what `ReadGroups` mints and what the
    /// sample-keyed view is in — so a fixture whose samples appear in a different order from
    /// its read groups is the one that catches a view built by read-group order.
    fn run_of(libraries: &[(&str, &str)]) -> ReadGroups {
        ReadGroups::of_libraries(libraries)
    }

    #[test]
    fn the_default_puts_every_read_group_and_every_sample_in_one_batch() {
        let groups = run_of(&[("a", "s1"), ("b", "s2"), ("c", "s2")]);
        let batches = SequencingBatches::all_together(&groups);
        assert_eq!(
            batches.of_each_read_group(),
            BatchOfEachReadGroup(&[BatchId::ALL_TOGETHER; 3])
        );
        assert_eq!(
            batches.of_each_sample(),
            BatchOfEachSample(&[BatchId::ALL_TOGETHER; 2]),
            "the run has three libraries and two samples, so the two views differ in length"
        );
        assert_eq!(batches.batch_count(), 1);
        assert!(batches.is_default());
    }

    /// **The two views are different lengths and are keyed by different things**, which is the
    /// whole reason they are two types. A run of three libraries over two samples returns three
    /// entries on one axis and two on the other.
    #[test]
    fn the_two_views_are_keyed_by_what_they_say_they_are() {
        // Read groups 0 and 2 are sample `late`; read group 1 is sample `early`, so the run's
        // sample order is [late, early] — **first seen, not alphabetical**, which is what
        // `ReadGroups` mints. A sample-keyed view built in alphabetical order would come back
        // [batch 1, batch 0]: the same two values, the same length, and every shape check
        // passed. The names are chosen so that the two orders disagree.
        let groups = run_of(&[("a", "late"), ("b", "early"), ("c", "late")]);
        let batches = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(0), ReadGroupId(2)]),
                BTreeSet::from([ReadGroupId(1)]),
            ],
        )
        .expect("a partition of the run");
        assert_eq!(
            batches.of_each_read_group(),
            BatchOfEachReadGroup(&[BatchId(0), BatchId(1), BatchId(0)])
        );
        assert_eq!(
            batches.of_each_sample(),
            BatchOfEachSample(&[BatchId(0), BatchId(1)]),
            "sample `late` is read groups 0 and 2, both in batch 0; sample `early` is read \
             group 1, in batch 1"
        );
        assert_eq!(batches.batch_count(), 2);
        assert!(!batches.is_default());
    }

    #[test]
    fn a_read_group_nobody_batched_is_refused_by_name() {
        let groups = run_of(&[("a", "s1"), ("b", "s2"), ("c", "s3")]);
        let refusal = SequencingBatches::declared(
            &groups,
            &[BTreeSet::from([ReadGroupId(0), ReadGroupId(2)])],
        )
        .expect_err("read group 1 was left out");
        assert_eq!(
            refusal,
            SequencingBatchError::ReadGroupNotBatched {
                read_groups: vec![1]
            }
        );
    }

    #[test]
    fn a_read_group_in_two_batches_is_refused() {
        let groups = run_of(&[("a", "s1"), ("b", "s2")]);
        let refusal = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
                BTreeSet::from([ReadGroupId(1)]),
            ],
        )
        .expect_err("read group 1 is in both");
        assert_eq!(
            refusal,
            SequencingBatchError::ReadGroupInSeveralBatches {
                read_group: 1,
                first: 0,
                second: 1
            }
        );
    }

    #[test]
    fn a_batch_naming_a_read_group_the_run_does_not_have_is_refused() {
        let groups = run_of(&[("a", "s1")]);
        let refusal = SequencingBatches::declared(
            &groups,
            &[BTreeSet::from([ReadGroupId(0), ReadGroupId(7)])],
        )
        .expect_err("read group 7 is not this run's");
        assert_eq!(
            refusal,
            SequencingBatchError::UnknownReadGroup {
                batch: 0,
                read_group: 7,
                read_groups: 1
            }
        );
    }

    #[test]
    fn an_empty_batch_is_refused() {
        let groups = run_of(&[("a", "s1")]);
        let refusal = SequencingBatches::declared(
            &groups,
            &[BTreeSet::new(), BTreeSet::from([ReadGroupId(0)])],
        )
        .expect_err("batch 0 holds nothing");
        assert_eq!(refusal, SequencingBatchError::EmptyBatch { batch: 0 });
    }

    #[test]
    fn a_declaration_of_no_batches_is_refused_and_names_the_default() {
        let groups = run_of(&[("a", "s1")]);
        let refusal =
            SequencingBatches::declared(&groups, &[]).expect_err("no batches were declared");
        assert_eq!(refusal, SequencingBatchError::NoBatches);
    }

    /// **The rule nobody settled, refused loudly.** A sample with two libraries on two plates
    /// has no single contaminating population, and the alternatives — pick the majority, or
    /// average the two — would score half its reads against the wrong neighbours with nothing
    /// said.
    #[test]
    fn a_sample_sequenced_in_two_batches_is_refused_by_name() {
        let groups = run_of(&[("a", "split"), ("b", "split"), ("c", "whole")]);
        let refusal = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(0)]),
                BTreeSet::from([ReadGroupId(1), ReadGroupId(2)]),
            ],
        )
        .expect_err("sample `split` ran on both plates");
        assert_eq!(
            refusal,
            SequencingBatchError::SampleSequencedInSeveralBatches {
                sample: "split".to_string(),
                batches: vec![0, 1]
            }
        );
    }

    /// The same run, batched so that each sample's libraries stay together, is accepted — so
    /// the refusal above is about the *split*, not about a sample having two libraries.
    #[test]
    fn a_sample_with_two_libraries_in_one_batch_is_accepted() {
        let groups = run_of(&[("a", "split"), ("b", "split"), ("c", "whole")]);
        let batches = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
                BTreeSet::from([ReadGroupId(2)]),
            ],
        )
        .expect("each sample's libraries ran together");
        assert_eq!(
            batches.of_each_sample(),
            BatchOfEachSample(&[BatchId(0), BatchId(1)])
        );
    }

    /// **A declaration of one batch holding everything is not the default**, and only
    /// `is_default` can tell them apart: the dense views are identical.
    #[test]
    fn one_declared_batch_holding_the_run_is_not_the_default() {
        let groups = run_of(&[("a", "s1"), ("b", "s2")]);
        let declared = SequencingBatches::declared(
            &groups,
            &[BTreeSet::from([ReadGroupId(0), ReadGroupId(1)])],
        )
        .expect("a partition of the run");
        let defaulted = SequencingBatches::all_together(&groups);
        assert_eq!(
            declared.of_each_read_group(),
            defaulted.of_each_read_group(),
            "the two dense views are the same value, which is why the distinction needs a \
             field of its own"
        );
        assert!(!declared.is_default());
        assert!(defaulted.is_default());
    }

    /// **The two doors to the default agree**, which is what lets the second exist: one takes
    /// the run's read-group table and the other its two axis lengths, and a run of three
    /// libraries over two samples gets the same value from either.
    #[test]
    fn the_two_doors_to_the_default_agree() {
        let groups = run_of(&[("a", "s1"), ("b", "s2"), ("c", "s2")]);
        assert_eq!(
            SequencingBatches::all_together(&groups),
            SequencingBatches::all_together_over(3, 2)
        );
    }

    /// **A run with no read groups is one whose read groups went missing**, not a run with
    /// none: every read of a run belongs to one. Both doors refuse it, because the dense views
    /// they build would otherwise come back empty and every consumer would read a zero-length
    /// axis as a legal answer.
    #[test]
    #[should_panic(expected = "read-group table is empty")]
    fn the_default_over_a_run_with_no_read_groups_is_refused() {
        let _ = SequencingBatches::all_together(&run_of(&[]));
    }

    /// And the declared door, which is the one a user reaches.
    #[test]
    #[should_panic(expected = "read-group table is empty")]
    fn a_declaration_over_a_run_with_no_read_groups_is_refused() {
        let _ = SequencingBatches::declared(&run_of(&[]), &[BTreeSet::from([ReadGroupId(0)])]);
    }

    /// **A run with read groups always has samples**, which is why the sample-axis check beside
    /// the one above is debug-only: this is what makes it unreachable rather than merely
    /// untested.
    #[test]
    fn a_run_with_read_groups_always_has_at_least_one_sample() {
        for libraries in [
            vec![("a", "s1")],
            vec![("a", "s1"), ("b", "s1")],
            vec![("a", "s1"), ("b", "s2"), ("c", "s1")],
        ] {
            let groups = run_of(&libraries);
            assert!(!groups.read_groups_per_sample().is_empty());
        }
    }

    /// Batch ids are the declaration's own order, so a run may declare its plates in any order
    /// and the frequency table's rows follow it.
    #[test]
    fn batch_ids_are_the_declarations_own_order() {
        let groups = run_of(&[("a", "s1"), ("b", "s2"), ("c", "s3")]);
        let batches = SequencingBatches::declared(
            &groups,
            &[
                BTreeSet::from([ReadGroupId(2)]),
                BTreeSet::from([ReadGroupId(0)]),
                BTreeSet::from([ReadGroupId(1)]),
            ],
        )
        .expect("a partition of the run");
        assert_eq!(
            batches.of_each_read_group(),
            BatchOfEachReadGroup(&[BatchId(1), BatchId(2), BatchId(0)])
        );
        assert_eq!(batches.batch_count(), 3);
    }
}
