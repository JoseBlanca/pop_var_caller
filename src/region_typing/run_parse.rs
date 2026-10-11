//! How one allele of an STR bundle reads against the hull's reference runs: how many copies it
//! carries in each run, and where in the allele each run lies (`str_bundle_loci.md` §7.1).
//!
//! A bundle's hull is a few runs — stretches of whole copies of one motif, found once in the
//! reference (`reference_runs.rs`) — with fixed sequence between and around them. An allele the
//! reads show is that hull with some runs longer or shorter, and perhaps a substituted base. The
//! **run parse** recovers that reading: it spells the allele as
//!
//! ```text
//! fixed₀  motif₁ × c₁  part₁  fixed₁  motif₂ × c₂  part₂  …  fixed_R
//! ```
//!
//! — each run as whole copies of its motif followed by the part copy the reference run ends on,
//! the fixed sequence as the reference has it — choosing the counts `c₁ … c_R` that need **the
//! fewest edited bases** (a substitution, an inserted base or a deleted one each cost one). Ties go
//! to the counts **nearest the reference's**, summed over the runs.
//!
//! The counts are the allele's **run counts**: the bundle's counterpart of a segment's repeat
//! count, which keys the stutter lookup and the prior and is written as `REPCNB`. The place of
//! each run in the allele is where a stutter slip in that run is placed (§7.3).
//!
//! An allele with a substitution or a part-unit change still parses — the edits are its
//! substitutions, and a part-unit change shows as an inserted or deleted base; it simply costs
//! more edits than a clean allele, and the count of edits is the purity a candidate tie-break
//! reads (§7.2). A bundle with no run parses every allele as fixed sequence alone: no counts.
//!
//! # The method
//!
//! The template above is a small automaton — a chain of states for each fixed piece and part copy,
//! and for each run a chain of copies with an exit after every whole copy — and the parse is an
//! edit-distance alignment of the allele against it, by dynamic programming over (template state,
//! allele position). A hull is at most 100 bases (`--max-str-bundle-len`) and an allele about as
//! long, so a parse is tens of thousands of cells.

use crate::region_typing::reference_runs::ReferenceRun;

/// Where one reference run lies in an allele, and how many whole copies of its motif the allele
/// carries there.
///
/// `start..end` is a range of the allele's own bases. It holds the run's whole copies and its part
/// copy, and any base the parse had to insert where the run meets fixed sequence — so the runs'
/// places plus the fixed pieces between them tile the allele, and a run whose span is the whole
/// hull has the whole allele as its place.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunInAllele {
    /// The run's first base in the allele.
    pub start: u32,
    /// One past the run's last base in the allele.
    pub end: u32,
    /// Whole copies of the run's motif the allele carries there — the run count.
    pub copies: u32,
}

/// One allele read against a bundle's reference runs: each run's place and count, in hull order,
/// and how many bases the reading had to edit.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RunParse {
    runs: Box<[RunInAllele]>,
    edited_bases: u32,
}

impl RunParse {
    /// Each reference run's place and count in the allele, in hull order — one per run of the
    /// bundle, empty for a bundle with no run.
    #[inline]
    #[must_use]
    pub fn runs(&self) -> &[RunInAllele] {
        &self.runs
    }

    /// How many bases the reading edited — substituted, inserted or deleted. Zero for the
    /// reference allele, and for any allele that differs from it only by whole copies.
    #[inline]
    #[must_use]
    pub fn edited_bases(&self) -> u32 {
        self.edited_bases
    }

    /// The run counts, in hull order.
    pub fn copies(&self) -> impl ExactSizeIterator<Item = u32> + '_ {
        self.runs.iter().map(|run| run.copies)
    }
}

/// A step of the template automaton: from state `from` to the state it belongs to, either
/// spelling one base or spelling nothing (the exit from a run after a whole copy).
#[derive(Debug, Clone, Copy)]
struct TemplateEdge {
    from: u32,
    /// The base this step spells, or `None` for a step that spells nothing.
    base: Option<u8>,
    /// For a run's exit: which run, and after how many whole copies. Its distance from the
    /// reference's count is the tie-break.
    exit: Option<RunExit>,
}

#[derive(Debug, Clone, Copy)]
struct RunExit {
    run: u32,
    copies: u32,
    distance_from_reference: u32,
}

/// The template as states, each with the edges that reach it. States are numbered so every edge
/// goes from a lower number to a higher one, which is the order the alignment fills them.
struct Template {
    /// Edges into each state, indexed by state.
    edges_into: Vec<Vec<TemplateEdge>>,
    /// Each run's entry state and exit state, in hull order.
    run_ends: Vec<(u32, u32)>,
    /// Whether each run comes straight after another run, with no fixed base between.
    run_follows_a_run: Vec<bool>,
}

impl Template {
    fn new() -> Self {
        Self {
            edges_into: vec![Vec::new()],
            run_ends: Vec::new(),
            run_follows_a_run: Vec::new(),
        }
    }

    fn last_state(&self) -> u32 {
        (self.edges_into.len() - 1) as u32
    }

    /// A new state reached from `from` by spelling `base`.
    fn spell(&mut self, from: u32, base: u8) -> u32 {
        self.edges_into.push(vec![TemplateEdge {
            from,
            base: Some(base),
            exit: None,
        }]);
        self.last_state()
    }
}

/// Build the template of `hull`'s fixed pieces and `runs`, allowing each run up to `most_copies`
/// whole copies.
fn build_template(hull: &[u8], runs: &[ReferenceRun], allele_len: usize) -> Template {
    let mut template = Template::new();
    let mut at = 0u32;
    let mut hull_at = 0usize;
    let mut previous_was_run = false;
    for (index, run) in runs.iter().enumerate() {
        let run_start = run.offset() as usize;
        for &base in &hull[hull_at..run_start] {
            at = template.spell(at, base);
        }
        template
            .run_follows_a_run
            .push(previous_was_run && hull_at == run_start);
        let entry = at;
        let motif = run.motif();
        let unit = motif.as_bytes();
        let period = unit.len();
        let reference_copies = run.copies() as u32;
        // More copies than the allele has bases for are never cheaper than fewer: each would be
        // `period` deleted bases.
        let most_copies = (allele_len / period) as u32 + 1;
        let mut after_copies = Vec::with_capacity(most_copies as usize + 1);
        after_copies.push(entry);
        let mut state = entry;
        for _ in 0..most_copies {
            for &base in unit {
                state = template.spell(state, base);
            }
            after_copies.push(state);
        }
        // One join state, reached from the end of every whole number of copies.
        template.edges_into.push(
            after_copies
                .iter()
                .enumerate()
                .map(|(copies, &from)| TemplateEdge {
                    from,
                    base: None,
                    exit: Some(RunExit {
                        run: index as u32,
                        copies: copies as u32,
                        distance_from_reference: (copies as u32).abs_diff(reference_copies),
                    }),
                })
                .collect(),
        );
        let mut exit = template.last_state();
        let part_copy = &hull[run_start + run.copies() as usize * period..run.end() as usize];
        for &base in part_copy {
            exit = template.spell(exit, base);
        }
        template.run_ends.push((entry, exit));
        at = exit;
        hull_at = run.end() as usize;
        previous_was_run = true;
    }
    for &base in &hull[hull_at..] {
        at = template.spell(at, base);
    }
    template
}

/// What the alignment keeps per cell: edited bases, then distance of the counts from the
/// reference's, compared in that order.
type Cost = (u32, u32);

/// How a cell was reached, for the walk back.
#[derive(Debug, Clone, Copy)]
enum Move {
    Start,
    /// An allele base with no template base: from the same state, one allele base earlier.
    Inserted,
    /// The edge at this index of the state's incoming list, spelling an allele base.
    Spelled(u32),
    /// The edge at this index, spelling a template base the allele has not got.
    Deleted(u32),
    /// The edge at this index, spelling nothing (a run's exit).
    Exited(u32),
}

/// **Read `allele` against the bundle's reference runs** — `hull` is the reference's own hull
/// bases, `runs` its reference runs in hull order (`str_bundle_loci.md` §7.1).
///
/// The parse with the fewest edited bases, ties to the run counts nearest the reference's; the
/// remaining ties go to whichever reading the alignment meets first, which is fixed, so a parse is
/// deterministic. The reference hull itself parses to its runs' own copy counts with no edit.
///
/// # Panics
///
/// If a run reaches past the hull, or the runs are out of hull order or overlap — the reference
/// runs a bundle locus carries never do (`find_reference_runs`).
#[must_use]
pub fn parse_into_runs(allele: &[u8], hull: &[u8], runs: &[ReferenceRun]) -> RunParse {
    assert!(
        runs.windows(2)
            .all(|pair| pair[0].end() <= pair[1].offset())
            && runs.last().is_none_or(|run| run.end() <= hull.len() as u64),
        "a bundle's reference runs lie inside its hull, in order and without overlap"
    );
    if runs.is_empty() {
        return RunParse {
            runs: Box::new([]),
            edited_bases: edit_distance(allele, hull),
        };
    }
    let template = build_template(hull, runs, allele.len());
    let states = template.edges_into.len();
    let columns = allele.len() + 1;
    let mut cost: Vec<Cost> = vec![(u32::MAX, u32::MAX); states * columns];
    let mut how: Vec<Move> = vec![Move::Start; states * columns];
    let cell = |state: usize, column: usize| state * columns + column;
    cost[cell(0, 0)] = (0, 0);

    for column in 0..columns {
        for state in 0..states {
            let mut best = cost[cell(state, column)];
            let mut best_move = how[cell(state, column)];
            let mut offer = |candidate: Cost, candidate_move: Move| {
                if candidate < best {
                    best = candidate;
                    best_move = candidate_move;
                }
            };
            for (index, edge) in template.edges_into[state].iter().enumerate() {
                let from = edge.from as usize;
                match edge.base {
                    Some(base) => {
                        if column > 0 {
                            let (edits, distance) = cost[cell(from, column - 1)];
                            if edits != u32::MAX {
                                let mismatch = u32::from(allele[column - 1] != base);
                                offer((edits + mismatch, distance), Move::Spelled(index as u32));
                            }
                        }
                        let (edits, distance) = cost[cell(from, column)];
                        if edits != u32::MAX {
                            offer((edits + 1, distance), Move::Deleted(index as u32));
                        }
                    }
                    None => {
                        let (edits, distance) = cost[cell(from, column)];
                        if edits != u32::MAX {
                            let added = edge.exit.map_or(0, |exit| exit.distance_from_reference);
                            offer((edits, distance + added), Move::Exited(index as u32));
                        }
                    }
                }
            }
            if column > 0 {
                let (edits, distance) = cost[cell(state, column - 1)];
                if edits != u32::MAX {
                    offer((edits + 1, distance), Move::Inserted);
                }
            }
            cost[cell(state, column)] = best;
            how[cell(state, column)] = best_move;
        }
    }

    // Walk back from the last state at the allele's end, noting each run's count and the
    // allele positions at which the walk stood on each run's entry and exit states.
    let mut placed = vec![
        RunInAllele {
            start: 0,
            end: 0,
            copies: 0,
        };
        runs.len()
    ];
    let mut first_column_at = vec![u32::MAX; states];
    let mut last_column_at = vec![0u32; states];
    let (mut state, mut column) = (states - 1, allele.len());
    let edited_bases = cost[cell(state, column)].0;
    loop {
        first_column_at[state] = first_column_at[state].min(column as u32);
        last_column_at[state] = last_column_at[state].max(column as u32);
        match how[cell(state, column)] {
            Move::Start => break,
            Move::Inserted => column -= 1,
            Move::Spelled(index) => {
                state = template.edges_into[state][index as usize].from as usize;
                column -= 1;
            }
            Move::Deleted(index) => {
                state = template.edges_into[state][index as usize].from as usize;
            }
            Move::Exited(index) => {
                let edge = template.edges_into[state][index as usize];
                if let Some(exit) = edge.exit {
                    placed[exit.run as usize].copies = exit.copies;
                }
                state = edge.from as usize;
            }
        }
    }
    for (run, &(entry, exit)) in template.run_ends.iter().enumerate() {
        // Bases inserted where a run meets fixed sequence are the run's; where it meets the run
        // before it, the earlier run's.
        placed[run].start = if template.run_follows_a_run[run] {
            last_column_at[entry as usize]
        } else {
            first_column_at[entry as usize]
        };
        placed[run].end = last_column_at[exit as usize];
    }
    RunParse {
        runs: placed.into_boxed_slice(),
        edited_bases,
    }
}

/// Edited bases between two sequences: substitutions, insertions and deletions, one each.
fn edit_distance(one: &[u8], other: &[u8]) -> u32 {
    let mut previous: Vec<u32> = (0..=other.len() as u32).collect();
    let mut current = vec![0u32; other.len() + 1];
    for (row, &base) in one.iter().enumerate() {
        current[0] = row as u32 + 1;
        for (column, &other_base) in other.iter().enumerate() {
            current[column + 1] = (previous[column] + u32::from(base != other_base))
                .min(previous[column + 1] + 1)
                .min(current[column] + 1);
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[other.len()]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Motif;

    fn run(offset: u64, length: u64, motif: &[u8]) -> ReferenceRun {
        ReferenceRun::new(offset, length, Motif::new(motif).unwrap()).unwrap()
    }

    fn places(parse: &RunParse) -> Vec<(u32, u32, u32)> {
        parse
            .runs()
            .iter()
            .map(|run| (run.start, run.end, run.copies))
            .collect()
    }

    /// `(CA)8 GT (TG)5` style: two runs with two fixed bases between them, and a flank piece of
    /// fixed sequence at each end of the hull.
    fn two_run_hull() -> (Vec<u8>, Vec<ReferenceRun>) {
        let hull = [&b"GG"[..], &b"CA".repeat(8), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let runs = vec![run(2, 16, b"CA"), run(20, 15, b"TTG")];
        (hull, runs)
    }

    #[test]
    fn the_reference_hull_parses_to_its_own_runs_with_no_edit() {
        let (hull, runs) = two_run_hull();
        let parse = parse_into_runs(&hull, &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(2, 18, 8), (20, 35, 5)]);
    }

    #[test]
    fn a_whole_copy_more_in_one_run_is_counted_there_with_no_edit() {
        let (hull, runs) = two_run_hull();
        let allele = [&b"GG"[..], &b"CA".repeat(8), b"GT", &b"TTG".repeat(7), b"C"].concat();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(2, 18, 8), (20, 41, 7)]);

        let allele = [&b"GG"[..], &b"CA".repeat(5), b"GT", &b"TTG".repeat(5), b"C"].concat();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(2, 12, 5), (14, 29, 5)]);
    }

    /// A substitution inside a run keeps its count and costs one edit.
    #[test]
    fn a_substitution_inside_a_run_costs_one_edit_and_keeps_the_count() {
        let (hull, runs) = two_run_hull();
        let mut allele = hull.clone();
        allele[7] = b'G'; // the third copy's A
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 1);
        assert_eq!(places(&parse), vec![(2, 18, 8), (20, 35, 5)]);
    }

    /// A one-base insertion inside a dinucleotide run is a part-unit change: one edit, the count
    /// nearest the reference's, and the extra base inside the run's place.
    #[test]
    fn a_part_unit_change_is_one_edit_and_stays_in_the_runs_place() {
        let (hull, runs) = two_run_hull();
        let allele = [
            &b"GG"[..],
            &b"CA".repeat(8),
            b"C",
            b"GT",
            &b"TTG".repeat(5),
            b"C",
        ]
        .concat();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 1);
        assert_eq!(parse.copies().collect::<Vec<_>>(), vec![8, 5]);
        assert_eq!(places(&parse), vec![(2, 19, 8), (21, 36, 5)]);
    }

    /// An allele that lost the hull's first two (fixed) bases to a deletion across the left edge:
    /// the bases' absence is two edits, and the runs keep their counts.
    #[test]
    fn a_deletion_of_fixed_bases_at_the_edge_costs_its_bases_and_keeps_the_counts() {
        let (hull, runs) = two_run_hull();
        let allele = hull[2..].to_vec();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 2);
        assert_eq!(places(&parse), vec![(0, 16, 8), (18, 33, 5)]);
    }

    #[test]
    fn a_bundle_with_no_run_parses_to_no_counts_and_its_edit_distance() {
        let hull = b"ACGTTGCAAGCT".to_vec();
        let parse = parse_into_runs(b"ACGTTGCAGCT", &hull, &[]);
        assert!(parse.runs().is_empty());
        assert_eq!(parse.edited_bases(), 1);
        assert_eq!(parse_into_runs(&hull, &hull, &[]).edited_bases(), 0);
    }

    /// A run that ends on a part copy: `CACACAC` is three `CA` copies and a `C`.
    #[test]
    fn a_run_with_a_part_copy_counts_whole_copies_only() {
        let hull = b"CACACAC".to_vec();
        let runs = vec![run(0, 7, b"CA")];
        let parse = parse_into_runs(&hull, &hull, &runs);
        assert_eq!(places(&parse), vec![(0, 7, 3)]);
        let parse = parse_into_runs(b"CACACACAC", &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(0, 9, 4)]);
    }

    /// One more base than a whole-copy change explains: the tie between one inserted base and
    /// one deleted base a copy further goes to the count nearest the reference's, and the run's
    /// place still covers the whole allele — what lets a one-run bundle score as a segment.
    #[test]
    fn a_tie_goes_to_the_count_nearest_the_reference_and_the_place_covers_the_allele() {
        let hull = b"CACACAC".to_vec();
        let runs = vec![run(0, 7, b"CA")];
        let parse = parse_into_runs(b"CACACACA", &hull, &runs);
        assert_eq!(parse.edited_bases(), 1);
        assert_eq!(places(&parse), vec![(0, 8, 3)]);
    }

    /// Two adjacent runs with no fixed base between them, the §6.1 example `(CT)20 (AT)8 T8`.
    #[test]
    fn adjacent_runs_split_the_allele_between_them() {
        let hull = [b"CT".repeat(20), b"AT".repeat(8), b"T".repeat(8)].concat();
        let runs = vec![run(0, 40, b"CT"), run(40, 16, b"AT"), run(56, 8, b"T")];
        let parse = parse_into_runs(&hull, &hull, &runs);
        assert_eq!(places(&parse), vec![(0, 40, 20), (40, 56, 8), (56, 64, 8)]);
        let allele = [b"CT".repeat(21), b"AT".repeat(8), b"T".repeat(7)].concat();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(0, 42, 21), (42, 58, 8), (58, 65, 7)]);
    }

    #[test]
    fn a_run_can_lose_every_copy() {
        let (hull, runs) = two_run_hull();
        let allele = [&b"GG"[..], b"GT", &b"TTG".repeat(5), b"C"].concat();
        let parse = parse_into_runs(&allele, &hull, &runs);
        assert_eq!(parse.edited_bases(), 0);
        assert_eq!(places(&parse), vec![(2, 2, 0), (4, 19, 5)]);
    }

    #[test]
    fn the_edit_distance_counts_each_kind_of_edit_once() {
        assert_eq!(edit_distance(b"", b""), 0);
        assert_eq!(edit_distance(b"ACGT", b"ACGT"), 0);
        assert_eq!(edit_distance(b"ACGT", b"AGGT"), 1);
        assert_eq!(edit_distance(b"ACGT", b"ACGGT"), 1);
        assert_eq!(edit_distance(b"ACGT", b"AGT"), 1);
        assert_eq!(edit_distance(b"", b"ACG"), 3);
    }
}
