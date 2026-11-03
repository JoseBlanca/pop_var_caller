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

use ahash::AHashMap;

use super::keeper::{FragmentAllele, READ_REACH, RecordFragments};
use super::link::{GQ_MIN, Heterozygote, Relation, read_link};
use crate::types::ContigId;
use crate::vcf::RecordPlace;

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
}

/// One sample's column at one line, as the stage reads it.
#[derive(Debug, Clone, Copy)]
struct Eligible {
    het: Heterozygote,
}

/// What settling decided for one sample at one line.
#[derive(Debug, Clone, Copy)]
struct Placed {
    block: u64,
    /// Whether copy 0 carries the upper allele.
    flipped: bool,
    phase_set: u64,
}

#[derive(Debug)]
struct HeldLine {
    place: RecordPlace,
    line: Vec<u8>,
    fragments: RecordFragments,
    /// Per sample, its heterozygote when the column is eligible.
    eligible: Vec<Option<Eligible>>,
    /// Per sample, where settling placed it; filled when the line settles.
    placed: Vec<Option<Placed>>,
    settled: bool,
}

/// One record of a sample's open block, kept while it is within reach of what comes next.
#[derive(Debug)]
struct Member {
    at: u64,
    het: Heterozygote,
    flipped: bool,
    fragments: Vec<FragmentAllele>,
}

#[derive(Debug)]
struct OpenBlock {
    id: u64,
    phase_set: u64,
    size: u32,
    /// Its records within reach of the newest, oldest first; the newest is the last.
    members: VecDeque<Member>,
}

/// One sample's state along the contig.
#[derive(Debug, Default)]
struct SamplePhase {
    /// Each fragment's first position and how many eligible records it has a read allele at.
    eligible_count: AHashMap<u32, (u64, u32)>,
    /// The same fragments in first-seen order, to forget them from the front.
    first_seen: VecDeque<(u64, u32)>,
    block: Option<OpenBlock>,
    /// Closed blocks and their sizes, oldest first, while a held line may still name them.
    closed: VecDeque<(u64, u32)>,
}

/// **The stage** — see the module documentation.
#[derive(Debug)]
pub struct ReadPhaseStage {
    samples: usize,
    held: VecDeque<HeldLine>,
    contig: Option<ContigId>,
    per_sample: Vec<SamplePhase>,
    next_block: u64,
    tally: ReadPhaseTally,
}

impl ReadPhaseStage {
    /// A stage over `samples` sample columns.
    #[must_use]
    pub fn new(samples: usize) -> Self {
        Self {
            samples,
            held: VecDeque::new(),
            contig: None,
            per_sample: (0..samples).map(|_| SamplePhase::default()).collect(),
            next_block: 0,
            tally: ReadPhaseTally::default(),
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
        let at = place.at.position.get();

        let eligible = eligible_columns(&line, self.samples);
        for (sample, column) in eligible.iter().enumerate() {
            if column.is_none() {
                continue;
            }
            let state = &mut self.per_sample[sample];
            for entry in fragments.of_sample(sample as u32) {
                let counted = state
                    .eligible_count
                    .entry(entry.fragment)
                    .or_insert_with(|| {
                        state.first_seen.push_back((at, entry.fragment));
                        (at, 0)
                    });
                counted.1 += 1;
            }
        }
        self.held.push_back(HeldLine {
            place,
            line,
            fragments,
            placed: vec![None; eligible.len()],
            eligible,
            settled: false,
        });

        // Settled: every line more than `R` behind this one. Released: more than `2R`.
        for index in 0..self.held.len() {
            let held_at = self.held[index].place.at.position.get();
            if held_at + READ_REACH >= at {
                break;
            }
            if !self.held[index].settled {
                self.settle(index);
            }
        }
        while self
            .held
            .front()
            .is_some_and(|oldest| oldest.place.at.position.get() + 2 * READ_REACH < at)
        {
            let oldest = self.held.pop_front().expect("a front, just seen");
            ready.push(self.release(oldest));
        }
        // A fragment first seen more than `2R` back cannot be at a line still to settle, which
        // starts no earlier than `R` back, its fragments no earlier than `R` before that.
        for state in &mut self.per_sample {
            while let Some(&(first_at, fragment)) = state.first_seen.front()
                && first_at + 2 * READ_REACH < at
            {
                state.first_seen.pop_front();
                state.eligible_count.remove(&fragment);
            }
        }
    }

    /// **Settles and releases everything held**: at a contig change and at the end of the run.
    pub fn finish(&mut self, ready: &mut Vec<(RecordPlace, Vec<u8>)>) {
        for index in 0..self.held.len() {
            if !self.held[index].settled {
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
            state.eligible_count.clear();
            state.first_seen.clear();
            state.closed.clear();
        }
    }

    /// **Settle one held line**: for each sample eligible there, whether it is phaseable, and if
    /// so, which block it joins or opens.
    fn settle(&mut self, index: usize) {
        let held = &mut self.held[index];
        held.settled = true;
        let at = held.place.at.position.get();
        let phase_set_here = at;
        for sample in 0..self.samples {
            let Some(Eligible { het }) = held.eligible[sample] else {
                continue;
            };
            let state = &mut self.per_sample[sample];
            let here: Vec<FragmentAllele> = held.fragments.of_sample(sample as u32).to_vec();
            let phaseable = here.iter().any(|entry| {
                state
                    .eligible_count
                    .get(&entry.fragment)
                    .is_some_and(|&(_, count)| count >= 2)
            });
            if !phaseable {
                continue;
            }

            let joined = state.block.as_mut().and_then(|block| {
                while block
                    .members
                    .front()
                    .is_some_and(|member| member.at + READ_REACH < at)
                {
                    block.members.pop_front();
                }
                let last = block.members.back()?;
                let relation = read_link(&last.fragments, last.het, &here, het).relation()?;
                let flipped = last.flipped ^ (relation == Relation::Trans);
                // Every earlier member within reach, tested against the relation the block
                // implies between it and this record.
                let contradicted = block.members.iter().any(|member| {
                    let implied = if member.flipped == flipped {
                        Relation::Cis
                    } else {
                        Relation::Trans
                    };
                    read_link(&member.fragments, member.het, &here, het).contradicts(implied)
                });
                Some((block.id, block.phase_set, flipped, contradicted))
            });
            match joined {
                Some((block, phase_set, flipped, false)) => {
                    self.tally.links += 1;
                    let open = state.block.as_mut().expect("joined an open block");
                    open.size += 1;
                    open.members.push_back(Member {
                        at,
                        het,
                        flipped,
                        fragments: here,
                    });
                    held.placed[sample] = Some(Placed {
                        block,
                        flipped,
                        phase_set,
                    });
                }
                joined => {
                    if matches!(joined, Some((_, _, _, true))) {
                        self.tally.links_refused_by_an_earlier_record += 1;
                    }
                    if let Some(closed) = state.block.take() {
                        self.tally.blocks += u64::from(closed.size >= 2);
                        state.closed.push_back((closed.id, closed.size));
                    }
                    let id = self.next_block;
                    self.next_block += 1;
                    state.block = Some(OpenBlock {
                        id,
                        phase_set: phase_set_here,
                        size: 1,
                        members: VecDeque::from([Member {
                            at,
                            het,
                            flipped: false,
                            fragments: here,
                        }]),
                    });
                    held.placed[sample] = Some(Placed {
                        block: id,
                        flipped: false,
                        phase_set: phase_set_here,
                    });
                }
            }
        }
    }

    /// **Release one settled line, rewritten**: `PS` appended to `FORMAT` and to every column,
    /// `.` where the column is not phased; a phased column's `GT` written `a|b`, copy 0 first.
    fn release(&mut self, held: HeldLine) -> (RecordPlace, Vec<u8>) {
        self.tally.lines += 1;
        let mut phased: Vec<Option<(Heterozygote, Placed)>> = vec![None; self.samples];
        for (sample, placed) in held.placed.iter().enumerate() {
            let Some(placed) = *placed else {
                continue;
            };
            let state = &mut self.per_sample[sample];
            while state
                .closed
                .front()
                .is_some_and(|&(block, _)| block < placed.block)
            {
                state.closed.pop_front();
            }
            let size = match &state.block {
                Some(open) if open.id == placed.block => open.size,
                _ => state
                    .closed
                    .iter()
                    .find(|&&(block, _)| block == placed.block)
                    .map_or(1, |&(_, size)| size),
            };
            if size >= 2 {
                let het = held.eligible[sample]
                    .expect("placed only where eligible")
                    .het;
                phased[sample] = Some((het, placed));
                self.tally.columns_phased += 1;
            }
        }
        (held.place, with_phase(&held.line, &phased))
    }
}

/// **Each sample's column, read for eligibility**: its heterozygote where the line is `PASS`, the
/// column's `GT` names two different alleles and its `GQ` is at least [`GQ_MIN`]. A line with no
/// sample columns, or one the stage cannot read, has none eligible.
fn eligible_columns(line: &[u8], samples: usize) -> Vec<Option<Eligible>> {
    let mut eligible = vec![None; samples];
    let Ok(text) = std::str::from_utf8(line) else {
        return eligible;
    };
    let fields: Vec<&str> = text.split('\t').collect();
    if fields.len() != 9 + samples || fields[6] != "PASS" {
        return eligible;
    }
    let format: Vec<&str> = fields[8].split(':').collect();
    let (Some(genotype_at), Some(quality_at)) = (
        format.iter().position(|key| *key == "GT"),
        format.iter().position(|key| *key == "GQ"),
    ) else {
        return eligible;
    };
    for (sample, column) in fields[9..].iter().enumerate() {
        let values: Vec<&str> = column.split(':').collect();
        let (Some(genotype), Some(quality)) = (values.get(genotype_at), values.get(quality_at))
        else {
            continue;
        };
        let Ok(quality) = quality.parse::<u32>() else {
            continue;
        };
        if quality < u32::from(GQ_MIN) {
            continue;
        }
        let mut alleles = genotype.split(['/', '|']).map(str::parse::<u16>);
        let (Some(Ok(one)), Some(Ok(other)), None) =
            (alleles.next(), alleles.next(), alleles.next())
        else {
            continue;
        };
        eligible[sample] = Heterozygote::of(one, other).map(|het| Eligible { het });
    }
    eligible
}

/// **The line with `PS` on every column**: `FORMAT` gains `PS` after its last key; each column is
/// padded with `.` to `FORMAT`'s length, so `PS` lands in its own place, and gains its block's
/// `POS` where it is phased — its `GT` then `a|b`, copy 0's allele first — and `.` elsewhere.
fn with_phase(line: &[u8], phased: &[Option<(Heterozygote, Placed)>]) -> Vec<u8> {
    let Ok(text) = std::str::from_utf8(line) else {
        return line.to_vec();
    };
    let fields: Vec<&str> = text.split('\t').collect();
    if fields.len() < 10 {
        return line.to_vec();
    }
    let keys = fields[8].split(':').count();
    let mut out = String::with_capacity(line.len() + 3 * (fields.len() - 8));
    for (at, field) in fields[..8].iter().enumerate() {
        if at > 0 {
            out.push('\t');
        }
        out.push_str(field);
    }
    out.push('\t');
    out.push_str(fields[8]);
    out.push_str(":PS");
    for (sample, column) in fields[9..].iter().enumerate() {
        out.push('\t');
        let present = column.split(':').count();
        match phased.get(sample).copied().flatten() {
            Some((het, placed)) => {
                let (copy_0, copy_1) = if placed.flipped {
                    (het.upper, het.lower)
                } else {
                    (het.lower, het.upper)
                };
                // `GT` is the first key of every column calling writes.
                let rest = column.split_once(':').map_or("", |(_, rest)| rest);
                out.push_str(&format!("{copy_0}|{copy_1}"));
                if !rest.is_empty() {
                    out.push(':');
                    out.push_str(rest);
                }
                for _ in present..keys {
                    out.push_str(":.");
                }
                out.push(':');
                out.push_str(&placed.phase_set.to_string());
            }
            None => {
                out.push_str(column);
                for _ in present..keys {
                    out.push_str(":.");
                }
                out.push_str(":.");
            }
        }
    }
    out.into_bytes()
}

#[cfg(test)]
mod tests;
