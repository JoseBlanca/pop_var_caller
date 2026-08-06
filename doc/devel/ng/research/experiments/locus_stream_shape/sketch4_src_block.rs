//! **SKETCH — throwaway.** The columnar alternative to one owned locus per covered base.
//!
//! Built for `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md` sketch 1 and
//! deleted after the decision. Nothing here commits an interface.
//!
//! A [`LocusBlock`] is many loci with each field held as one array shared across the
//! block: flat byte buffers plus offset arrays for the ragged fields, and an index
//! saying where each locus's observations begin. A [`LocusSink`] is what the walk
//! emits into — either today's `VecDeque` of owned records (arm A) or a block
//! (arms B and C).

use std::collections::VecDeque;

use crate::ng::locus_generation::{ReadWitness, SampleLocusObservations};
use crate::ng::types::ReadGroupId;
use crate::pileup_record::ChainId;

/// A run of positions a partial witness covered, in locus coordinates.
pub type WitnessRun = (u16, u16);

/// Many loci, each field one array.
///
/// Every ragged field is a flat buffer plus an offsets array carrying `n + 1` entries,
/// the first of which is zero — so locus `i`'s slice is `off[i]..off[i + 1]` and no
/// length is stored twice.
#[derive(Debug, Default)]
pub struct LocusBlock {
    // ---- per locus ----
    pub contig: Vec<u32>,
    /// 1-based inclusive, as `GenomeRegion` is.
    pub start: Vec<u32>,
    pub end: Vec<u32>,
    pub reads_without_observation: Vec<u32>,
    pub reads_discarded_by_cap: Vec<u32>,
    /// `n + 1` entries into `ref_bases`.
    pub ref_bases_off: Vec<u32>,
    pub ref_bases: Vec<u8>,
    /// `n + 1` entries into the per-observation arrays.
    pub obs_off: Vec<u32>,
    /// **SKETCH 4.** The cheap per-locus summary the cohort merge folds: summed
    /// `num_obs` over the observations whose bases differ from this locus's reference
    /// bases. Empty unless the block was built with `summary` on.
    ///
    /// It exists because nothing else in the block answers "is this locus worth
    /// materialising" without touching `obs_bases` — the heaviest array here. Sketch 2's
    /// `.psp` fold read the equivalent straight out of the file, which stores its alleles
    /// reference-first; ng's observations carry raw bases and no reference marker, so the
    /// comparison has to happen somewhere, and the cheapest somewhere is at fill time,
    /// when the reference bases are already in cache.
    pub locus_nonref_obs: Vec<u32>,

    // ---- per observation ----
    /// `m + 1` entries into `obs_bases`.
    pub obs_bases_off: Vec<u32>,
    pub obs_bases: Vec<u8>,
    /// `0` complete, `1` partial. A partial's runs are `obs_runs[obs_runs_off[j]..]`.
    pub obs_witness_kind: Vec<u8>,
    /// `m + 1` entries into `obs_runs`.
    pub obs_runs_off: Vec<u32>,
    pub obs_runs: Vec<WitnessRun>,
    pub obs_read_group: Vec<u32>,
    pub obs_num_obs: Vec<u32>,
    pub obs_num_fwd: Vec<u32>,
    pub obs_q_sum: Vec<f64>,
    pub obs_mapq_sum: Vec<u32>,
    pub obs_mapq_sum_sq: Vec<u64>,
    pub obs_placed_left: Vec<u32>,
    /// `m + 1` entries into `obs_chain_ids`.
    pub obs_chain_off: Vec<u32>,
    pub obs_chain_ids: Vec<ChainId>,

    /// Running payload size, kept incrementally so the budget test is a compare.
    bytes: usize,
    /// The block is full at or above this many payload bytes. Zero means unbounded.
    budget_bytes: usize,
    /// **SKETCH 4.** Whether `locus_nonref_obs` is maintained. Constant for the life of a
    /// run, so the branch it guards predicts perfectly — the same device as
    /// `LocusSink::wants_records`.
    summary: bool,
    /// **SKETCH 4.** Where the open locus's reference bases start in `ref_bases`.
    open_ref_lo: usize,
}

/// Bytes one observation adds to `bytes`, excluding its ragged parts.
const OBSERVATION_FIXED_BYTES: usize = 4 + 1 + 4 + 4 + 4 + 4 + 8 + 4 + 8 + 4 + 4;
/// Bytes one locus adds to `bytes`, excluding its ragged parts.
const LOCUS_FIXED_BYTES: usize = 4 + 4 + 4 + 4 + 4 + 4 + 4;

impl LocusBlock {
    /// An empty block that reports itself full at `budget_bytes` of payload.
    pub fn with_budget(budget_bytes: usize) -> Self {
        Self::with_budget_and_summary(budget_bytes, false)
    }

    /// **SKETCH 4.** As [`with_budget`](Self::with_budget), maintaining the cheap
    /// per-locus summary column when `summary` is set.
    pub fn with_budget_and_summary(budget_bytes: usize, summary: bool) -> Self {
        let mut block = Self {
            budget_bytes,
            summary,
            ..Self::default()
        };
        block.clear();
        block
    }

    /// Empty the block, keeping every allocation. The offset arrays keep their leading
    /// zero, which is what makes `off[i]..off[i + 1]` total.
    pub fn clear(&mut self) {
        self.contig.clear();
        self.start.clear();
        self.end.clear();
        self.reads_without_observation.clear();
        self.reads_discarded_by_cap.clear();
        self.ref_bases.clear();
        self.obs_bases.clear();
        self.obs_runs.clear();
        self.obs_witness_kind.clear();
        self.obs_read_group.clear();
        self.obs_num_obs.clear();
        self.obs_num_fwd.clear();
        self.obs_q_sum.clear();
        self.obs_mapq_sum.clear();
        self.obs_mapq_sum_sq.clear();
        self.obs_placed_left.clear();
        self.obs_chain_ids.clear();
        self.locus_nonref_obs.clear();
        self.open_ref_lo = 0;
        for offsets in [
            &mut self.ref_bases_off,
            &mut self.obs_off,
            &mut self.obs_bases_off,
            &mut self.obs_runs_off,
            &mut self.obs_chain_off,
        ] {
            offsets.clear();
            offsets.push(0);
        }
        self.bytes = 0;
    }

    /// How many loci the block holds.
    pub fn len(&self) -> usize {
        self.contig.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contig.is_empty()
    }

    /// How many observations the block holds, over all its loci.
    pub fn observation_count(&self) -> usize {
        self.obs_num_obs.len()
    }

    /// Payload bytes currently held — what the budget is judged against.
    pub fn payload_bytes(&self) -> usize {
        self.bytes
    }

    /// Bytes the block's arrays have reserved, whether used or not.
    pub fn capacity_bytes(&self) -> usize {
        self.contig.capacity() * 4
            + self.start.capacity() * 4
            + self.end.capacity() * 4
            + self.reads_without_observation.capacity() * 4
            + self.reads_discarded_by_cap.capacity() * 4
            + self.ref_bases_off.capacity() * 4
            + self.ref_bases.capacity()
            + self.obs_off.capacity() * 4
            + self.obs_bases_off.capacity() * 4
            + self.obs_bases.capacity()
            + self.obs_witness_kind.capacity()
            + self.obs_runs_off.capacity() * 4
            + self.obs_runs.capacity() * 4
            + self.obs_read_group.capacity() * 4
            + self.obs_num_obs.capacity() * 4
            + self.obs_num_fwd.capacity() * 4
            + self.obs_q_sum.capacity() * 8
            + self.obs_mapq_sum.capacity() * 4
            + self.obs_mapq_sum_sq.capacity() * 8
            + self.obs_placed_left.capacity() * 4
            + self.obs_chain_off.capacity() * 4
            + self.obs_chain_ids.capacity() * 8
            + self.locus_nonref_obs.capacity() * 4
    }

    /// **SKETCH 4.** Locus `i`'s summed non-reference observations — the cheap column the
    /// cohort merge folds. Zero when the summary is off, which is why the fold arm asks
    /// for it explicitly.
    #[inline]
    pub fn nonref_obs_of(&self, i: usize) -> u32 {
        self.locus_nonref_obs[i]
    }

    /// Whether the block has reached its budget. Tested after each whole locus, never
    /// inside one — a locus is never split across two blocks.
    pub fn is_full(&self) -> bool {
        self.budget_bytes != 0 && self.bytes >= self.budget_bytes
    }

    /// Start a locus. Every `begin_locus` must be closed by exactly one `end_locus`,
    /// with its observations pushed between the two.
    pub fn begin_locus(
        &mut self,
        contig: u32,
        start: u32,
        end: u32,
        reference_bases: &[u8],
        reads_without_observation: u32,
        reads_discarded_by_cap: u32,
    ) {
        self.contig.push(contig);
        self.start.push(start);
        self.end.push(end);
        self.reads_without_observation
            .push(reads_without_observation);
        self.reads_discarded_by_cap.push(reads_discarded_by_cap);
        if self.summary {
            self.open_ref_lo = self.ref_bases.len();
            self.locus_nonref_obs.push(0);
            self.bytes += 4;
        }
        self.ref_bases.extend_from_slice(reference_bases);
        self.ref_bases_off.push(self.ref_bases.len() as u32);
        self.bytes += LOCUS_FIXED_BYTES + reference_bases.len();
    }

    /// Append one observation to the locus opened by `begin_locus`.
    #[allow(clippy::too_many_arguments)]
    pub fn push_observation(
        &mut self,
        bases: &[u8],
        witness_runs: Option<&[WitnessRun]>,
        read_group: ReadGroupId,
        num_obs: u32,
        num_fwd: u32,
        q_sum: f64,
        mapq_sum: u32,
        mapq_sum_sq: u64,
        placed_left: u32,
        chain_ids: &[ChainId],
    ) {
        if self.summary && bases != &self.ref_bases[self.open_ref_lo..] {
            // PANIC-FREE: `begin_locus` pushed the entry this observation belongs to.
            let slot = self
                .locus_nonref_obs
                .last_mut()
                .expect("an observation is only pushed inside an open locus");
            *slot = slot.saturating_add(num_obs);
        }
        self.obs_bases.extend_from_slice(bases);
        self.obs_bases_off.push(self.obs_bases.len() as u32);
        match witness_runs {
            None => self.obs_witness_kind.push(0),
            Some(runs) => {
                self.obs_witness_kind.push(1);
                self.obs_runs.extend_from_slice(runs);
            }
        }
        self.obs_runs_off.push(self.obs_runs.len() as u32);
        self.obs_read_group.push(read_group.0);
        self.obs_num_obs.push(num_obs);
        self.obs_num_fwd.push(num_fwd);
        self.obs_q_sum.push(q_sum);
        self.obs_mapq_sum.push(mapq_sum);
        self.obs_mapq_sum_sq.push(mapq_sum_sq);
        self.obs_placed_left.push(placed_left);
        self.obs_chain_ids.extend_from_slice(chain_ids);
        self.obs_chain_off.push(self.obs_chain_ids.len() as u32);
        self.bytes += OBSERVATION_FIXED_BYTES
            + bases.len()
            + witness_runs.map_or(0, |runs| runs.len() * 4)
            + chain_ids.len() * 8;
    }

    /// Close the locus opened by `begin_locus`.
    pub fn end_locus(&mut self) {
        self.obs_off.push(self.obs_num_obs.len() as u32);
    }

    /// Copy a whole owned record in — the path the general fold takes when it has not
    /// been converted, and the path the SSR generator would take.
    pub fn push_record(&mut self, locus: &SampleLocusObservations) {
        self.begin_locus(
            locus.region.contig.0,
            locus.region.start.get() as u32,
            locus.region.end.get() as u32,
            &locus.reference_bases,
            locus.reads_without_observation,
            locus.reads_discarded_by_cap,
        );
        for observation in &locus.observations {
            let runs: Option<&[WitnessRun]> = match &observation.read_witness {
                ReadWitness::Complete => None,
                ReadWitness::Partial { positions } => Some(positions.runs_slice()),
            };
            self.push_observation(
                &observation.bases,
                runs,
                observation.read_group,
                observation.num_obs,
                observation.num_fwd,
                observation.q_sum,
                observation.mapq_sum,
                observation.mapq_sum_sq,
                observation.placed_left,
                &observation.chain_ids,
            );
        }
        self.end_locus();
    }

    /// The observation index range of locus `i`.
    #[inline]
    pub fn observations_of(&self, i: usize) -> std::ops::Range<usize> {
        self.obs_off[i] as usize..self.obs_off[i + 1] as usize
    }

    /// Locus `i`'s reference bases.
    #[inline]
    pub fn reference_bases_of(&self, i: usize) -> &[u8] {
        &self.ref_bases[self.ref_bases_off[i] as usize..self.ref_bases_off[i + 1] as usize]
    }

    /// Observation `j`'s bases.
    #[inline]
    pub fn bases_of(&self, j: usize) -> &[u8] {
        &self.obs_bases[self.obs_bases_off[j] as usize..self.obs_bases_off[j + 1] as usize]
    }

    /// Observation `j`'s chain ids.
    #[inline]
    pub fn chain_ids_of(&self, j: usize) -> &[ChainId] {
        &self.obs_chain_ids[self.obs_chain_off[j] as usize..self.obs_chain_off[j + 1] as usize]
    }

    /// Observation `j`'s witness runs, empty for a complete witness.
    #[inline]
    pub fn runs_of(&self, j: usize) -> &[WitnessRun] {
        &self.obs_runs[self.obs_runs_off[j] as usize..self.obs_runs_off[j + 1] as usize]
    }

    /// Rebuild locus `i` as the owned record the record path would have emitted — what
    /// the random-locus sampler does with a locus it keeps, since a view cannot outlive
    /// its block. Returns the bytes copied, which the sketch reports.
    pub fn copy_out(&self, i: usize) -> (SampleLocusObservations, usize) {
        use crate::ng::locus_generation::{
            LocusKind, SequenceObservation, WitnessedLocusPositions,
        };
        use crate::ng::types::{ContigId, GenomeRegion, Position};

        let mut copied = self.reference_bases_of(i).len();
        let mut observations = Vec::with_capacity(self.observations_of(i).len());
        for j in self.observations_of(i) {
            let runs = self.runs_of(j);
            let read_witness = if self.obs_witness_kind[j] == 0 {
                ReadWitness::Complete
            } else {
                ReadWitness::Partial {
                    // The stored runs are already canonical, so this round-trips.
                    positions: WitnessedLocusPositions::from_half_open_runs(
                        runs.iter().copied(),
                    )
                    .expect("a stored partial witness is never empty"),
                }
            };
            copied += self.bases_of(j).len() + runs.len() * 4 + self.chain_ids_of(j).len() * 8;
            observations.push(SequenceObservation {
                bases: self.bases_of(j).to_vec().into_boxed_slice(),
                read_witness,
                read_group: ReadGroupId(self.obs_read_group[j]),
                num_obs: self.obs_num_obs[j],
                num_fwd: self.obs_num_fwd[j],
                q_sum: self.obs_q_sum[j],
                mapq_sum: self.obs_mapq_sum[j],
                mapq_sum_sq: self.obs_mapq_sum_sq[j],
                placed_left: self.obs_placed_left[j],
                chain_ids: self.chain_ids_of(j).to_vec(),
            });
        }
        (
            SampleLocusObservations {
                region: GenomeRegion {
                    contig: ContigId(self.contig[i]),
                    start: Position(u64::from(self.start[i])),
                    end: Position(u64::from(self.end[i])),
                },
                reference_bases: self.reference_bases_of(i).to_vec().into_boxed_slice(),
                observations,
                reads_without_observation: self.reads_without_observation[i],
                reads_discarded_by_cap: self.reads_discarded_by_cap[i],
                kind: LocusKind::Generic,
            },
            copied,
        )
    }
}

/// Where the walk puts what it produces: today's owned records, or a block.
///
/// One enum rather than a type parameter, deliberately. Making the walk generic over its
/// sink would ripple a parameter through `PileupGenerator`, `GeneratorSlot`,
/// `GeneratorSet` and `SampleLocusObservationsIterator` — including the STR generator,
/// which this sketch does not touch. The cost is one predictable branch per locus, the
/// same branch in every arm, against a per-locus cost measured in thousands of
/// instructions.
#[derive(Debug)]
pub struct LocusSink {
    records: VecDeque<SampleLocusObservations>,
    block: LocusBlock,
    columns: bool,
    /// The region the generator clamps loci to. The record path applies this in
    /// `PileupGenerator::next_locus`, after the walker has yielded; the columnar path has
    /// no such moment, so it applies it here.
    clamp: Option<(u32, u32)>,
    /// Loci the clamp rejected — the columnar path's share of
    /// `PileupGeneratorCounts::records_outside_region`.
    pub records_outside_region: u64,
}

impl LocusSink {
    /// A sink that collects owned records: the shipped behaviour.
    pub fn records() -> Self {
        Self {
            records: VecDeque::new(),
            block: LocusBlock::default(),
            columns: false,
            clamp: None,
            records_outside_region: 0,
        }
    }

    /// A sink that fills a block, full at `budget_bytes` of payload.
    pub fn columns(budget_bytes: usize) -> Self {
        Self::columns_with_summary(budget_bytes, false)
    }

    /// **SKETCH 4.** As [`columns`](Self::columns), with the cohort merge's cheap
    /// per-locus summary column maintained.
    pub fn columns_with_summary(budget_bytes: usize, summary: bool) -> Self {
        Self {
            records: VecDeque::new(),
            block: LocusBlock::with_budget_and_summary(budget_bytes, summary),
            columns: true,
            clamp: None,
            records_outside_region: 0,
        }
    }

    /// Whether the emit sites should build owned records. Constant for the life of a run,
    /// so the branch it guards predicts perfectly.
    #[inline(always)]
    pub fn wants_records(&self) -> bool {
        !self.columns
    }

    /// Point the columnar path's clamp at `region`.
    pub fn set_clamp(&mut self, start: u32, end: u32) {
        self.clamp = Some((start, end));
    }

    /// Whether a locus anchored at `start` is inside the clamp. Loci outside are counted
    /// and dropped, exactly as `PileupGenerator::next_locus` does on the record path.
    #[inline]
    pub fn accepts(&mut self, start: u32) -> bool {
        match self.clamp {
            None => true,
            Some((from, to)) => {
                let inside = start >= from && start <= to;
                if !inside {
                    self.records_outside_region += 1;
                }
                inside
            }
        }
    }

    /// Push an owned record — the record path, and the columnar path's fallback for
    /// producers this sketch did not convert.
    #[inline]
    pub fn push_record(&mut self, locus: SampleLocusObservations) {
        if self.columns {
            if self.accepts(locus.region.start.get() as u32) {
                self.block.push_record(&locus);
            }
        } else {
            self.records.push_back(locus);
        }
    }

    /// The block, for the columnar path's direct writes.
    #[inline]
    pub fn block_mut(&mut self) -> &mut LocusBlock {
        &mut self.block
    }

    pub fn block(&self) -> &LocusBlock {
        &self.block
    }

    /// The record queue, for the walk's own draining.
    #[inline]
    pub fn records_mut(&mut self) -> &mut VecDeque<SampleLocusObservations> {
        &mut self.records
    }

    /// Whether anything is waiting — a record on the record path, a locus on the
    /// columnar one.
    #[inline]
    pub fn has_output(&self) -> bool {
        if self.columns {
            !self.block.is_empty()
        } else {
            !self.records.is_empty()
        }
    }

    /// Whether the walk should stop and hand the block over.
    #[inline]
    pub fn is_full(&self) -> bool {
        self.columns && self.block.is_full()
    }

    /// Whether the walk's fill loop has produced enough to return.
    ///
    /// The record path returns as soon as one locus exists, because its consumer takes
    /// them one at a time. The columnar path keeps walking until the block reaches its
    /// budget — that is what the block *is*.
    #[inline]
    pub fn should_yield(&self) -> bool {
        if self.columns {
            self.block.is_full()
        } else {
            !self.records.is_empty()
        }
    }

    #[inline]
    pub fn pop_record(&mut self) -> Option<SampleLocusObservations> {
        self.records.pop_front()
    }

    pub fn clear_records(&mut self) {
        self.records.clear();
    }
}
