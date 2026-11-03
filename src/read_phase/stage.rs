//! **Phase from reads alone** — the writer's stage when refinement does not run, at ploidy two
//! (read phase spec §2.6).
//!
//! It sits in the writer, not in calling, because `FILTER` is final only there: with the
//! hidden-duplication filter on, the verdict is written in pass three, which feeds the spilled
//! lines back through the writer. Like the neighbour-evidence pass, it reads `GT`, `GQ` and
//! `FILTER` from each line and rewrites `GT`.
//!
//! **Which records of a sample it may phase.** A record is *eligible* when it is `PASS`, the
//! sample's `GT` is heterozygous and its `GQ` is at least [`GQ_MIN`]. An eligible record is
//! *phaseable* when one of its fragments also has a read allele at another eligible record.
//! Eligible records without such a fragment are skipped, not treated as gaps.
//!
//! **Blocks.** A sample's phaseable records are walked left to right, and each is linked to the
//! one before it when the reads' log-odds reaches [`LINK_MIN`](super::link::LINK_MIN) either way. Before it joins, it is
//! tested against every earlier record of the block within [`READ_REACH`] — a pair whose mates
//! cover records 1 and 3, with 2 in the gap between them, links 1 and 3 but is never counted for
//! 1–2 or 2–3 — and a contradiction at [`BREAK_MIN`](super::link::BREAK_MIN) ends the block
//! instead. A block of one record
//! is written unphased. `PS` is the `POS` of the block's first record, which is written `a|b` with
//! its lower allele on copy 0.
//!
//! **Holding.** A record's phaseability needs its fragments' later read alleles, up to `R` ahead,
//! and whether its block grows past it can need another `R`; so a record is *settled* — made
//! phaseable or not, and linked — once the stream is more than `R` past it, and *released* once it
//! is more than `2R` past it. Lines leave in the order they came.

use std::collections::VecDeque;
use std::fmt::Write as _;

use ahash::AHashMap;

use super::fragments::{FragmentAllele, READ_REACH, RecordFragments};
use super::link::{GQ_MIN, Heterozygote, Relation, score_read_link};
use crate::types::ContigId;
use crate::vcf::RecordPlace;
use crate::vcf::line_fields::{LineFields, parse_genotype};

/// What the stage did, for the run's report.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct ReadPhaseTally {
    /// Lines that went through the stage.
    pub lines: u64,
    /// Sample columns written phased, `|` with a `PS`.
    pub columns_phased: u64,
    /// Blocks of two records or more written.
    pub blocks: u64,
    /// Links reads made between consecutive phaseable records.
    pub links: u64,
    /// Links reads made but refused because a fragment over an earlier record of the block
    /// contradicted the relation the block implies.
    pub links_refused_by_an_earlier_record: u64,
    /// Lines the stage could not read — not text, or fewer than the nine fixed columns — written
    /// as they came, without `PS`. Every line comes from this program's encoder, so this stays
    /// zero; the report warns if it does not.
    pub lines_not_read: u64,
}

/// Where settling put one sample at one line.
#[derive(Debug, Clone, Copy)]
struct BlockPlacement {
    block: u64,
    /// Whether copy 0 carries the upper allele.
    upper_on_copy_0: bool,
    phase_set: u64,
}

#[derive(Debug)]
struct HeldLine {
    place: RecordPlace,
    line: Vec<u8>,
    fragments: RecordFragments,
    /// Per sample, its heterozygote where the column is eligible.
    eligible: Vec<Option<Heterozygote>>,
    /// Per sample, where settling placed it; filled when the line settles.
    placed: Vec<Option<BlockPlacement>>,
    is_settled: bool,
}

/// One record of a sample's open block, kept while it is within reach of what comes next.
#[derive(Debug)]
struct BlockRecord {
    position: u64,
    heterozygote: Heterozygote,
    upper_on_copy_0: bool,
    fragments: Vec<FragmentAllele>,
}

#[derive(Debug)]
struct OpenBlock {
    id: u64,
    phase_set: u64,
    size: u32,
    /// Its records within reach of the newest, oldest first; the newest is the last.
    records: VecDeque<BlockRecord>,
}

/// A block that has ended, kept while a held line may still name it.
#[derive(Debug, Clone, Copy)]
struct ClosedBlock {
    id: u64,
    size: u32,
}

/// One fragment in the order it was first seen.
#[derive(Debug, Clone, Copy)]
struct FirstSeen {
    position: u64,
    fragment: u32,
}

/// One sample's state along the contig.
#[derive(Debug, Default)]
struct SampleReadBlocks {
    /// Each fragment's reach so far: how many eligible records it has a read allele at.
    eligible_records: AHashMap<u32, u32>,
    /// The same fragments in first-seen order, to forget them from the front.
    first_seen: VecDeque<FirstSeen>,
    block: Option<OpenBlock>,
    /// Closed blocks, oldest first.
    closed: VecDeque<ClosedBlock>,
}

/// What testing a phaseable record against its sample's open block found.
enum Join {
    /// It joins the block, oriented as the reads say.
    Extends {
        block: u64,
        phase_set: u64,
        upper_on_copy_0: bool,
    },
    /// The reads link it to the block's newest record, but a fragment over an earlier record of
    /// the block contradicts the relation the block implies: the block ends.
    RefusedByAnEarlierRecord,
    /// No open block, no record of it within reach, or no link at [`LINK_MIN`]: the block ends.
    ///
    /// [`LINK_MIN`]: super::link::LINK_MIN
    NoLink,
}

/// **Test a phaseable record against the open block** (§2.6): linked to the block's newest record
/// within reach, and contradicted by no earlier one within reach. Drops the block's records more
/// than `R` behind `position` first; a pair whose mates cover records 1 and 3, with 2 in the gap,
/// links 1 and 3 but is never counted for 1–2 or 2–3.
fn join_test(
    block: &mut OpenBlock,
    position: u64,
    heterozygote: Heterozygote,
    here: &[FragmentAllele],
) -> Join {
    while block
        .records
        .front()
        .is_some_and(|record| record.position + READ_REACH < position)
    {
        block.records.pop_front();
    }
    let Some(newest) = block.records.back() else {
        return Join::NoLink;
    };
    let Some(relation) =
        score_read_link(&newest.fragments, newest.heterozygote, here, heterozygote).relation()
    else {
        return Join::NoLink;
    };
    let upper_on_copy_0 = relation.carry(newest.upper_on_copy_0);
    let contradicted = block.records.iter().any(|record| {
        let implied = Relation::between(record.upper_on_copy_0, upper_on_copy_0);
        score_read_link(&record.fragments, record.heterozygote, here, heterozygote)
            .contradicts(implied)
    });
    if contradicted {
        Join::RefusedByAnEarlierRecord
    } else {
        Join::Extends {
            block: block.id,
            phase_set: block.phase_set,
            upper_on_copy_0,
        }
    }
}

/// **The stage** — see the module documentation.
#[derive(Debug)]
pub struct ReadPhaseStage {
    sample_count: usize,
    held: VecDeque<HeldLine>,
    contig: Option<ContigId>,
    per_sample: Vec<SampleReadBlocks>,
    next_block: u64,
    tally: ReadPhaseTally,
    /// Scratch for [`Self::release`]: per sample, what the released line writes phased.
    phased: Vec<Option<(Heterozygote, BlockPlacement)>>,
}

impl ReadPhaseStage {
    /// A stage over `sample_count` sample columns.
    #[must_use]
    pub fn new(sample_count: usize) -> Self {
        Self {
            sample_count,
            held: VecDeque::new(),
            contig: None,
            per_sample: (0..sample_count)
                .map(|_| SampleReadBlocks::default())
                .collect(),
            next_block: 0,
            tally: ReadPhaseTally::default(),
            phased: Vec::with_capacity(sample_count),
        }
    }

    /// What the stage has done so far.
    #[must_use]
    pub fn tally(&self) -> ReadPhaseTally {
        self.tally
    }

    /// **Takes one line**, its place already checked by the writer, with the record's fragments,
    /// and appends to `ready` every line it releases, rewritten.
    pub fn push(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        fragments: RecordFragments,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) {
        if self.contig.is_some_and(|contig| contig != place.at.contig) {
            self.finish(ready);
        }
        self.contig = Some(place.at.contig);
        let position = place.at.position.get();

        let eligible = eligible_columns(&line, self.sample_count);
        for (sample, heterozygote) in eligible.iter().enumerate() {
            if heterozygote.is_none() {
                continue;
            }
            let state = &mut self.per_sample[sample];
            for entry in fragments.of_sample(sample as u32) {
                let records = state
                    .eligible_records
                    .entry(entry.fragment)
                    .or_insert_with(|| {
                        state.first_seen.push_back(FirstSeen {
                            position,
                            fragment: entry.fragment,
                        });
                        0
                    });
                *records += 1;
            }
        }
        self.held.push_back(HeldLine {
            place,
            line,
            fragments,
            placed: vec![None; eligible.len()],
            eligible,
            is_settled: false,
        });

        // Settled: every line more than `R` behind this one. Released: more than `2R`.
        for index in 0..self.held.len() {
            if self.held[index].place.at.position.get() + READ_REACH >= position {
                break;
            }
            if !self.held[index].is_settled {
                self.settle(index);
            }
        }
        while let Some(oldest) = self.held.front()
            && oldest.place.at.position.get() + 2 * READ_REACH < position
        {
            // PANIC-FREE: the loop condition just saw a front.
            let oldest = self.held.pop_front().expect("a front, just seen");
            ready.push(self.release(oldest));
        }
        // A fragment first seen more than `2R` back cannot be at a line still to settle, which
        // starts no earlier than `R` back, its fragments no earlier than `R` before that.
        for state in &mut self.per_sample {
            while let Some(&FirstSeen {
                position: first_position,
                fragment,
            }) = state.first_seen.front()
                && first_position + 2 * READ_REACH < position
            {
                state.first_seen.pop_front();
                state.eligible_records.remove(&fragment);
            }
        }
    }

    /// **Settles and releases everything held**: at a contig change and at the end of the run.
    /// Every sample's fragment counts, open block and closed blocks are cleared, since fragment
    /// numbers restart on each contig.
    pub fn finish(&mut self, ready: &mut Vec<(RecordPlace, Vec<u8>)>) {
        for index in 0..self.held.len() {
            if !self.held[index].is_settled {
                self.settle(index);
            }
        }
        while let Some(oldest) = self.held.pop_front() {
            ready.push(self.release(oldest));
        }
        for state in &mut self.per_sample {
            if let Some(block) = state.block.take() {
                self.tally.blocks += u64::from(block.size >= 2);
            }
            state.eligible_records.clear();
            state.first_seen.clear();
            state.closed.clear();
        }
    }

    /// **Settle one held line**: for each sample eligible there, whether it is phaseable, and if
    /// so, which block it joins or opens.
    fn settle(&mut self, index: usize) {
        let HeldLine {
            place,
            fragments,
            eligible,
            placed,
            is_settled,
            ..
        } = &mut self.held[index];
        *is_settled = true;
        let position = place.at.position.get();
        for (sample, heterozygote) in eligible.iter().enumerate() {
            let Some(heterozygote) = *heterozygote else {
                continue;
            };
            let state = &mut self.per_sample[sample];
            let here = fragments.of_sample(sample as u32);
            let phaseable = here.iter().any(|entry| {
                state
                    .eligible_records
                    .get(&entry.fragment)
                    .is_some_and(|&records| records >= 2)
            });
            if !phaseable {
                continue;
            }
            let join = match state.block.as_mut() {
                Some(block) => join_test(block, position, heterozygote, here),
                None => Join::NoLink,
            };
            let record = |upper_on_copy_0| BlockRecord {
                position,
                heterozygote,
                upper_on_copy_0,
                fragments: here.to_vec(),
            };
            match (join, state.block.as_mut()) {
                (
                    Join::Extends {
                        block,
                        phase_set,
                        upper_on_copy_0,
                    },
                    Some(open),
                ) => {
                    self.tally.links += 1;
                    open.size += 1;
                    open.records.push_back(record(upper_on_copy_0));
                    placed[sample] = Some(BlockPlacement {
                        block,
                        upper_on_copy_0,
                        phase_set,
                    });
                }
                (join, _) => {
                    if matches!(join, Join::RefusedByAnEarlierRecord) {
                        self.tally.links_refused_by_an_earlier_record += 1;
                    }
                    if let Some(closed) = state.block.take() {
                        self.tally.blocks += u64::from(closed.size >= 2);
                        state.closed.push_back(ClosedBlock {
                            id: closed.id,
                            size: closed.size,
                        });
                    }
                    let id = self.next_block;
                    self.next_block += 1;
                    state.block = Some(OpenBlock {
                        id,
                        phase_set: position,
                        size: 1,
                        records: VecDeque::from([record(false)]),
                    });
                    placed[sample] = Some(BlockPlacement {
                        block: id,
                        upper_on_copy_0: false,
                        phase_set: position,
                    });
                }
            }
        }
    }

    /// **Release one settled line, rewritten**: `PS` in every column, `.` where the column is not
    /// phased; a phased column's `GT` written `a|b`, copy 0 first.
    ///
    /// A placed column's block is the open one or a closed one still in `closed`: blocks are
    /// numbered in the order they open, and `closed` loses only blocks numbered below the one this
    /// line, the oldest held, names — which no later line can name either.
    fn release(&mut self, held: HeldLine) -> (RecordPlace, Vec<u8>) {
        self.tally.lines += 1;
        self.phased.clear();
        self.phased.resize(self.sample_count, None);
        for (sample, placed) in held.placed.iter().enumerate() {
            let Some(placed) = *placed else {
                continue;
            };
            let state = &mut self.per_sample[sample];
            while state
                .closed
                .front()
                .is_some_and(|closed| closed.id < placed.block)
            {
                state.closed.pop_front();
            }
            let size = match &state.block {
                Some(open) if open.id == placed.block => open.size,
                _ => state
                    .closed
                    .iter()
                    .find(|closed| closed.id == placed.block)
                    .map_or(1, |closed| closed.size),
            };
            if size >= 2 {
                // PANIC-FREE: settling places a sample only where its column is eligible.
                let heterozygote = held.eligible[sample].expect("placed only where eligible");
                self.phased[sample] = Some((heterozygote, placed));
                self.tally.columns_phased += 1;
            }
        }
        match with_phase(&held.line, &self.phased) {
            Some(line) => (held.place, line),
            None => {
                self.tally.lines_not_read += 1;
                (held.place, held.line)
            }
        }
    }
}

/// **Each sample's column, read for eligibility**: its heterozygote where the line is `PASS`, the
/// column's `GT` names two different alleles and its `GQ` is at least [`GQ_MIN`]. A line with no
/// sample columns, or one the stage cannot read, has none eligible.
fn eligible_columns(line: &[u8], sample_count: usize) -> Vec<Option<Heterozygote>> {
    let mut eligible = vec![None; sample_count];
    let Some(fields) = LineFields::parse(line) else {
        return eligible;
    };
    if !fields.is_pass() || fields.column_count() != sample_count {
        return eligible;
    }
    let (Some(genotype_at), Some(quality_at)) = (fields.format_key("GT"), fields.format_key("GQ"))
    else {
        return eligible;
    };
    let last_needed = genotype_at.max(quality_at);
    for (slot, column) in eligible.iter_mut().zip(fields.columns()) {
        let (mut genotype, mut genotype_quality) = (None, None);
        for (index, value) in column.split(':').enumerate().take(last_needed + 1) {
            if index == genotype_at {
                genotype = Some(value);
            }
            if index == quality_at {
                genotype_quality = Some(value);
            }
        }
        let (Some(genotype), Some(Ok(genotype_quality))) =
            (genotype, genotype_quality.map(str::parse::<u32>))
        else {
            continue;
        };
        if genotype_quality < u32::from(GQ_MIN) {
            continue;
        }
        let Some((one, other)) = parse_genotype(genotype) else {
            continue;
        };
        let (Ok(one), Ok(other)) = (u16::try_from(one), u16::try_from(other)) else {
            continue;
        };
        *slot = Heterozygote::of(one, other);
    }
    eligible
}

/// **The line with `PS` in every column**, or `None` for a line the stage cannot read.
///
/// `PS` takes its place in `FORMAT` if the line has one, and goes after the last key otherwise.
/// Each column is padded with `.` to `FORMAT`'s length, so `PS` lands in its own place; a phased
/// column's `GT`, at `GT`'s place in `FORMAT`, becomes `a|b`, copy 0's allele first, and its `PS`
/// its block's `POS`. Every other column's `PS` is `.`.
fn with_phase(line: &[u8], phased: &[Option<(Heterozygote, BlockPlacement)>]) -> Option<Vec<u8>> {
    let fields = LineFields::parse(line)?;
    if fields.column_count() == 0 {
        // A line without sample columns has nothing to carry `PS`.
        return Some(line.to_vec());
    }
    let keys = fields.format_keys().count();
    let genotype_at = fields.format_key("GT");
    let phase_set_at = fields.format_key("PS");
    let field_count = keys.max(phase_set_at.map_or(keys + 1, |at| at + 1));
    let phase_set_at = phase_set_at.unwrap_or(keys);
    let mut out = String::with_capacity(line.len() + 8 * fields.column_count() + 3);
    for (index, field) in fields.fixed.iter().enumerate() {
        if index > 0 {
            out.push('\t');
        }
        out.push_str(field);
    }
    if phase_set_at == keys {
        out.push_str(":PS");
    }
    for (sample, column) in fields.columns().enumerate() {
        out.push('\t');
        let phase = phased.get(sample).copied().flatten();
        let write_value = |out: &mut String, index: usize, value: &str| {
            match phase {
                Some((_, placed)) if index == phase_set_at => {
                    // PANIC-FREE: writing to a `String` cannot fail.
                    write!(out, "{}", placed.phase_set).expect("a String takes any text");
                }
                None if index == phase_set_at => out.push('.'),
                Some((heterozygote, placed)) if Some(index) == genotype_at => {
                    let (copy_0, copy_1) = heterozygote.on_copies(placed.upper_on_copy_0);
                    // PANIC-FREE: as above.
                    write!(out, "{copy_0}|{copy_1}").expect("a String takes any text");
                }
                _ => out.push_str(value),
            }
        };
        let mut present = 0;
        for (index, value) in column.split(':').enumerate() {
            if index > 0 {
                out.push(':');
            }
            write_value(&mut out, index, value);
            present = index + 1;
        }
        for index in present..field_count {
            out.push(':');
            write_value(&mut out, index, ".");
        }
    }
    Some(out.into_bytes())
}

#[cfg(test)]
mod tests;
