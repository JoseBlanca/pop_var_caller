//! The reference runs of an STR bundle's hull: where in the reference sequence a slip of the
//! polymerase can happen, and in which motif (`str_bundle_loci.md` §4, §6.1).
//!
//! A bundle's members are what make it a bundle, but they overlap and nest (`str_bundle_loci.md`
//! §5.1: of HG002's 2,410 bundles, 551 hold one member wholly inside another), so they cannot say
//! where stutter goes. The runs can: each is a stretch of pure, whole copies of one motif of
//! period 1 to 6, at least the typing's copy floor for that period, and no two overlap. Every
//! hull base outside the runs is **fixed sequence**, where no stutter is modelled.
//!
//! **Every motif is looked for, not only the members'** (owner, 2026-10-11, in review). The
//! members' motifs missed the spec's own example: at GRCh38 chr1:8,007,443 the typing lists `CT`
//! and `T` and no `AT`, so the 16-base `AT` run between them came out as fixed sequence.

use crate::region_typing::segment_criteria::MinCopies;
use crate::types::{MAX_MOTIF_LEN, Motif};

/// One stretch of pure copies of one motif inside an STR bundle's hull: where it starts, how long
/// it is, and its repeat unit as the reference spells it at the run's first base.
///
/// **Whole copies, then at most one part copy at its end**, as a reference STR segment keeps the
/// part copy it ends on: `CACACAC` is a run of three `CA` copies and a `C`. Its copy count is its
/// length over its period, rounded down.
///
/// The motif is **phase-faithful** — the run's own first `period` bases — so it can be any
/// rotation of the unit: a run of `CT` copies that starts on a `T` has the motif `TC`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceRun {
    offset: u64,
    length: u64,
    motif: Motif,
}

impl ReferenceRun {
    /// A run starting `offset` bases into the hull, `length` bases long, of `motif`.
    ///
    /// `None` when the run holds no whole copy (`length` below the motif's length): such a
    /// stretch is not a run of anything, and a file claiming one is refused by its reader.
    pub fn new(offset: u64, length: u64, motif: Motif) -> Option<Self> {
        (length >= motif.period() as u64).then_some(Self {
            offset,
            length,
            motif,
        })
    }

    /// The run's first base, as bases from the hull's first base.
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// The run's length in bases, part copy included.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// One past the run's last base, as bases from the hull's first base.
    pub fn end(&self) -> u64 {
        self.offset + self.length
    }

    /// The repeat unit, phase-faithful: the run's own first `period` bases.
    pub fn motif(&self) -> Motif {
        self.motif
    }

    /// The motif's length in bases.
    pub fn period(&self) -> usize {
        self.motif.period()
    }

    /// Whole copies of the motif in the run: its length over its period, rounded down.
    pub fn copies(&self) -> u64 {
        self.length / self.motif.period() as u64
    }
}

/// The hull's reference runs, in hull order, by the three steps of `str_bundle_loci.md` §6.1.
///
/// 1. **For every period from 1 to 6**, every maximal stretch of `hull` that repeats with that
///    period. A stretch is **counted under its smallest period**: one whose unit is itself
///    repeats of a shorter unit (`ATAT` at period 4, `AA` at period 2) is the shorter period's
///    stretch and is skipped here. So a motif and its rotations are one motif, found once.
///    A unit holding any base other than `A`, `C`, `G` or `T` is no motif: no run of `N`.
/// 2. **Kept if it reaches the copy floor** for its period, `min_copies` — the typing's own
///    `--min-copies`, so a run is what the typing would call long enough to stutter.
/// 3. **Chosen longest first**, in bases, a tie going to the earlier start, then the shorter
///    period. A stretch that overlaps a run already chosen gives up the shared bases. Where it
///    gave them up at its end, it is cut back to whole copies; where at its start, it is re-read
///    from its new first base, which is still a whole copy of a rotation, and its natural end
///    keeps its part copy. It is kept if it still reaches the floor. The chosen runs never
///    overlap.
///
/// **A stretch is at least one whole copy and one more base**: one isolated copy of a motif,
/// with nothing repeating, is not found as a stretch. That only matters under a copy floor of 1,
/// which no default sets.
///
/// **One stretch can be split in two** by a run already chosen inside it — possible only when the
/// earlier choice was itself cut back and so is shorter than this stretch. Each piece is then
/// judged on its own. The spec does not name this case; the pieces are what "gives up the shared
/// bases" leaves.
///
/// An empty answer is a legitimate one: a cluster of impure repeats, none with a stretch at the
/// floor, has no run, and is called with no stutter modelled anywhere in it.
///
/// The bundle's members are not consulted (owner, 2026-10-11, in review): they overlap and nest,
/// and listing only their motifs missed runs no member names.
pub fn find_reference_runs(hull: &[u8], min_copies: &MinCopies) -> Vec<ReferenceRun> {
    let mut stretches: Vec<ReferenceRun> = Vec::new();
    for period in 1..=MAX_MOTIF_LEN {
        let floor = u64::from(min_copies.for_period(period as u8)).max(1);
        stretches.extend(
            periodic_stretches_of(hull, period).filter(|stretch| stretch.copies() >= floor),
        );
    }

    // Length, start and period decide every tie: two stretches alike in all three have the same
    // first `period` bases of the hull, so the same motif.
    stretches.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then(a.offset.cmp(&b.offset))
            .then(a.period().cmp(&b.period()))
    });

    let mut chosen: Vec<ReferenceRun> = Vec::new();
    for stretch in stretches {
        let floor = u64::from(min_copies.for_period(stretch.period() as u8)).max(1);
        for (start, end) in uncovered_pieces(stretch.offset, stretch.end(), &chosen) {
            let period = stretch.period() as u64;
            let length = if end < stretch.end() {
                // Cut at its end by a run chosen there: back to whole copies.
                (end - start) / period * period
            } else {
                end - start
            };
            if length / period < floor {
                continue;
            }
            let at = start as usize;
            // PANIC-FREE: a kept piece holds at least one whole copy (`length >= period`) and lies
            // inside its stretch, whose units are valid motifs.
            let motif = Motif::new(&hull[at..at + period as usize])
                .expect("a piece of a stretch holds one whole copy of a motif that was valid");
            chosen.push(ReferenceRun {
                offset: start,
                length,
                motif,
            });
        }
    }
    chosen.sort_by_key(|run| run.offset);
    chosen
}

/// Every maximal stretch of `hull` that repeats with `period` and whose unit is a motif of that
/// smallest period — not a power of a shorter unit, and of `A`, `C`, `G` and `T` only — at least
/// one copy and one more base long.
fn periodic_stretches_of(hull: &[u8], period: usize) -> impl Iterator<Item = ReferenceRun> + '_ {
    let last_comparable = hull.len().saturating_sub(period);
    // A base repeats the one a period further on: a maximal run of such bases from `start` to
    // `stop` (exclusive) is a stretch from `start` to `stop + period`.
    let mut at = 0;
    std::iter::from_fn(move || {
        while at < last_comparable {
            if hull[at] != hull[at + period] {
                at += 1;
                continue;
            }
            let start = at;
            while at < last_comparable && hull[at] == hull[at + period] {
                at += 1;
            }
            let unit = &hull[start..start + period];
            if is_a_motif_of_its_own_period(unit) {
                // PANIC-FREE: `unit` is 1 to 6 bases, the range `Motif::new` accepts.
                let motif = Motif::new(unit).expect("a unit of 1 to 6 bases is a motif");
                return Some(ReferenceRun {
                    offset: start as u64,
                    length: (at + period - start) as u64,
                    motif,
                });
            }
        }
        None
    })
}

/// Whether `unit` is a motif counted at its own length: every base one of `A`, `C`, `G`, `T`, and
/// not a repeat of a shorter unit — `AT` and `AAT` are, `AA` and `ATAT` are not.
fn is_a_motif_of_its_own_period(unit: &[u8]) -> bool {
    let period = unit.len();
    unit.iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
        && !(1..period).any(|shorter| {
            period.is_multiple_of(shorter)
                && (shorter..period).all(|at| unit[at] == unit[at - shorter])
        })
}

/// The parts of `[start, end)` that no run of `chosen` covers, in order.
fn uncovered_pieces(start: u64, end: u64, chosen: &[ReferenceRun]) -> Vec<(u64, u64)> {
    let mut covered: Vec<(u64, u64)> = chosen
        .iter()
        .filter(|run| run.offset < end && run.end() > start)
        .map(|run| (run.offset, run.end()))
        .collect();
    covered.sort_unstable();
    let mut pieces = Vec::new();
    let mut from = start;
    for (covered_start, covered_end) in covered {
        if covered_start > from {
            pieces.push((from, covered_start.min(end)));
        }
        from = from.max(covered_end);
    }
    if from < end {
        pieces.push((from, end));
    }
    pieces
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `(offset, length, motif)` of each run, for compact assertions.
    fn spelled(runs: &[ReferenceRun]) -> Vec<(u64, u64, String)> {
        runs.iter()
            .map(|run| {
                (
                    run.offset(),
                    run.length(),
                    String::from_utf8(run.motif().as_bytes().to_vec()).unwrap(),
                )
            })
            .collect()
    }

    fn calling_floors() -> MinCopies {
        MinCopies::default()
    }

    /// GRCh38 chr1:8,007,443, the first example of `str_bundle_loci.md` §5.1: `(CT)20`, `(AT)8`,
    /// then eight `T`s — the 64 bases as the reference has them. The typing lists only `CT` and `T`
    /// as members there, and the `AT` run is found all the same, because every motif is looked
    /// for. The `T` where the `AT` and `T` stretches meet goes to the longer `AT` stretch.
    #[test]
    fn the_ct_at_t_bundle_comes_out_as_three_runs_and_the_shared_t_goes_to_the_longer() {
        let hull = [b"CT".repeat(20), b"AT".repeat(8), b"T".repeat(8)].concat();
        assert_eq!(hull.len(), 64);
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![
                (0, 40, "CT".to_string()),
                (40, 16, "AT".to_string()),
                (56, 8, "T".to_string()),
            ]
        );
        assert_eq!(
            runs.iter().map(ReferenceRun::copies).collect::<Vec<_>>(),
            vec![20, 8, 8]
        );
    }

    /// GRCh38 chr1:28,487,365, the second example of §5.1: an `A` member over the first nine bases
    /// and an `AAAAAT` member over the rest. The `AAAAAT` stretches reach two whole copies, below
    /// the floor for period 6, so the hull is `A9` and fixed sequence.
    #[test]
    fn the_a_aaaaat_bundle_is_a9_and_fixed_sequence() {
        let hull = b"AAAAAAAAATAAAAATAAAATAAAAATAAA";
        let runs = find_reference_runs(hull, &calling_floors());
        assert_eq!(spelled(&runs), vec![(0, 9, "A".to_string())]);

        // With a floor of two copies for period 6, the 16-base `AAAAAT` stretch is chosen first;
        // the `A9` keeps four bases, under its floor of 8; the 14-base second stretch gives up the
        // four it shares and keeps ten, one copy.
        let floors = MinCopies::new([8, 6, 6, 6, 5, 2], 3);
        let runs = find_reference_runs(hull, &floors);
        assert_eq!(spelled(&runs), vec![(4, 16, "AAAAAT".to_string())]);
    }

    /// A cluster of impure repeats with no stretch at the floor has no run at all.
    #[test]
    fn a_bundle_of_impure_repeats_has_no_run() {
        let hull = b"ACACTCACACGTAGAGTAGAGA";
        let runs = find_reference_runs(hull, &calling_floors());
        assert!(runs.is_empty(), "{:?}", spelled(&runs));
    }

    /// A stretch keeps the reference's own phase as its motif: `TC`, not `CT`.
    #[test]
    fn a_run_keeps_its_own_phase_as_its_motif() {
        let hull = [b"GG".as_slice(), &b"TC".repeat(7), b"GG"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        // The stretch starts at the first T and ends on the last C: TCTC…TC, 14 bases.
        assert_eq!(spelled(&runs), vec![(2, 14, "TC".to_string())]);
    }

    /// A part copy at a stretch's end stays in the run and is not a copy.
    #[test]
    fn a_part_copy_at_the_end_stays_in_the_run() {
        // Eight CAG copies then CA, flanked by bases that continue nothing.
        let hull = [b"T".as_slice(), &b"CAG".repeat(8), b"CA", b"T"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(spelled(&runs), vec![(1, 26, "CAG".to_string())]);
        assert_eq!(runs[0].copies(), 8);
    }

    /// A stretch that shares bases with a longer one gives them up and keeps the rest, at either
    /// end. Period 1 and whole copies: what the two cut rules below have in common.
    #[test]
    fn a_stretch_gives_up_the_bases_it_shares_with_a_longer_run() {
        // (AC)7 A then A×12: the AC stretch is ACAC…ACA (15 bases), the A stretch AAAAAAAAAAAAA
        // from the AC stretch's last A (13 bases). Longest first: AC (15) is chosen; the A
        // stretch gives up its first base and keeps 12.
        let hull = [b"AC".repeat(7).as_slice(), b"A", &b"A".repeat(12), b"G"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![(0, 15, "AC".to_string()), (15, 12, "A".to_string())]
        );

        // Now the A stretch is longer than the AC one, so it is chosen first and the AC stretch
        // gives up its last A: 14 bases, (AC)7.
        let hull = [b"AC".repeat(7).as_slice(), b"A", &b"A".repeat(20), b"G"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![(0, 14, "AC".to_string()), (14, 21, "A".to_string())]
        );
    }

    /// **A stretch cut at its end goes back to whole copies.** `T (CAG)8 G×30 T`: the `G`
    /// stretch runs from the last `CAG`'s `G` (31 bases) and is chosen first; the `CAG` stretch
    /// (24 bases) gives up that `G`, leaving 23, and goes back to whole copies: `(CAG)7`, 21
    /// bases. Without the cut-back it would be 23.
    #[test]
    fn a_stretch_cut_at_its_end_goes_back_to_whole_copies() {
        let hull = [b"T".as_slice(), &b"CAG".repeat(8), &b"G".repeat(30), b"T"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![(1, 21, "CAG".to_string()), (24, 31, "G".to_string())]
        );
    }

    /// **A stretch cut only at its start keeps its natural end's part copy.**
    /// `T G×30 (CAG)8 CA T`: the period-3 stretch starts on the last `G` (`GCA…`, 27 bases); the
    /// 30-base `G` stretch is chosen first; the other is re-read from its next base as `CAG` and
    /// keeps the `CA` at its end: 26 bases. Always cutting back would give 24.
    #[test]
    fn a_stretch_cut_only_at_its_start_keeps_its_part_copy() {
        let hull = [
            b"T".as_slice(),
            &b"G".repeat(30),
            &b"CAG".repeat(8),
            b"CA",
            b"T",
        ]
        .concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![(1, 30, "G".to_string()), (31, 26, "CAG".to_string())]
        );
    }

    /// A stretch cut back below its floor is dropped.
    #[test]
    fn a_stretch_cut_below_the_floor_is_dropped() {
        // (AT)6 then twenty Ts: the AT stretch is ATAT…AT (12 bases, 6 copies, at the floor of
        // 6), the T stretch the AT run's last T and the twenty after it (21 bases). The T stretch
        // is chosen first; the AT stretch gives up its last T, is cut back to whole copies (10
        // bases, 5 copies), and falls below the floor.
        let hull = [b"G".as_slice(), &b"AT".repeat(6), &b"T".repeat(20), b"G"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(spelled(&runs), vec![(12, 21, "T".to_string())]);
    }

    /// **A stretch is counted under its smallest period**: `(AT)8` is one period-2 run, not also
    /// a period-4 or period-6 one, and twelve `A`s one period-1 run, not a period-2 `AA` run.
    #[test]
    fn a_stretch_is_counted_under_its_smallest_period() {
        let floors = MinCopies::uniform(2);
        let hull = [b"G".as_slice(), &b"AT".repeat(8), b"G"].concat();
        assert_eq!(
            spelled(&find_reference_runs(&hull, &floors)),
            vec![(1, 16, "AT".to_string())]
        );
        let hull = [b"G".as_slice(), &b"A".repeat(12), b"G"].concat();
        assert_eq!(
            spelled(&find_reference_runs(&hull, &floors)),
            vec![(1, 12, "A".to_string())]
        );
    }

    /// No run is made of `N`, or of any base but `A`, `C`, `G`, `T`, whatever the floor.
    #[test]
    fn no_run_is_found_in_n_bases() {
        let hull = [b"G".as_slice(), &b"N".repeat(20), &b"NA".repeat(8), b"G"].concat();
        assert!(find_reference_runs(&hull, &MinCopies::uniform(1)).is_empty());
    }

    /// A hull shorter than any period, or empty, has no run, and no index runs past it.
    #[test]
    fn a_hull_too_short_to_repeat_has_no_run() {
        assert!(find_reference_runs(b"", &MinCopies::uniform(1)).is_empty());
        assert!(find_reference_runs(b"C", &MinCopies::uniform(1)).is_empty());
        assert!(find_reference_runs(b"CA", &MinCopies::uniform(1)).is_empty());
    }

    /// Equal lengths: the earlier start wins the shared bases.
    #[test]
    fn a_tie_in_length_goes_to_the_earlier_start() {
        // Fifteen As, then (GA)7. The A stretch is the fifteen As; the GA stretch starts on the
        // last of them (A G A G … A) and is fifteen bases too.
        let hull = [b"C".as_slice(), &b"A".repeat(15), &b"GA".repeat(7), b"C"].concat();
        let runs = find_reference_runs(&hull, &calling_floors());
        assert_eq!(
            spelled(&runs),
            vec![(1, 15, "A".to_string()), (16, 14, "GA".to_string())],
            "the A stretch starts first and keeps the shared A; the GA stretch, re-read from its \
             next base, is (GA)7"
        );
    }

    /// A stretch can contain a run chosen before it, when that run was itself cut back; each side
    /// is then judged on its own.
    #[test]
    fn a_stretch_split_by_a_run_inside_it_keeps_each_piece_at_the_floor() {
        let pieces = uncovered_pieces(
            0,
            40,
            &[ReferenceRun::new(10, 6, Motif::new(b"A").unwrap()).unwrap()],
        );
        assert_eq!(pieces, vec![(0, 10), (16, 40)]);
        let pieces = uncovered_pieces(
            5,
            20,
            &[
                ReferenceRun::new(0, 8, Motif::new(b"A").unwrap()).unwrap(),
                ReferenceRun::new(18, 8, Motif::new(b"A").unwrap()).unwrap(),
            ],
        );
        assert_eq!(pieces, vec![(8, 18)]);
        assert!(
            uncovered_pieces(
                5,
                20,
                &[ReferenceRun::new(0, 30, Motif::new(b"A").unwrap()).unwrap()]
            )
            .is_empty()
        );
    }

    #[test]
    fn a_unit_is_a_motif_only_at_its_own_period_and_of_acgt() {
        assert!(is_a_motif_of_its_own_period(b"A"));
        assert!(is_a_motif_of_its_own_period(b"AT"));
        assert!(is_a_motif_of_its_own_period(b"AAT"));
        assert!(is_a_motif_of_its_own_period(b"AATG"));
        assert!(is_a_motif_of_its_own_period(b"ATATA"));
        assert!(is_a_motif_of_its_own_period(b"AAAAAT"));
        assert!(!is_a_motif_of_its_own_period(b"AA"));
        assert!(!is_a_motif_of_its_own_period(b"ATAT"));
        assert!(!is_a_motif_of_its_own_period(b"ATATAT"));
        assert!(!is_a_motif_of_its_own_period(b"AATAAT"));
        assert!(!is_a_motif_of_its_own_period(b"N"));
        assert!(!is_a_motif_of_its_own_period(b"ca"));
    }

    /// A run with no whole copy cannot be built.
    #[test]
    fn a_run_shorter_than_its_motif_is_refused() {
        assert!(ReferenceRun::new(0, 2, Motif::new(b"CAG").unwrap()).is_none());
        assert!(ReferenceRun::new(0, 3, Motif::new(b"CAG").unwrap()).is_some());
    }

    /// The runs never overlap and lie inside the hull, on a hull built to make them compete.
    #[test]
    fn the_chosen_runs_never_overlap() {
        let hull = b"TTTTTTTTTATATATATATATAAAAAAAAAGAGAGAGAGAGGGGGGGGGTCTCTCTCTCTCTTTTTTTT";
        let runs = find_reference_runs(hull, &MinCopies::uniform(3));
        assert!(!runs.is_empty());
        for pair in runs.windows(2) {
            assert!(pair[0].end() <= pair[1].offset(), "{:?}", spelled(&runs));
        }
        for run in &runs {
            assert!(run.end() <= hull.len() as u64);
            let bases = &hull[run.offset() as usize..run.end() as usize];
            let motif = run.motif();
            let unit = motif.as_bytes();
            assert!(
                bases
                    .iter()
                    .enumerate()
                    .all(|(at, base)| *base == unit[at % unit.len()]),
                "a run is pure copies of its own motif: {:?}",
                spelled(&runs)
            );
            assert!(run.copies() >= 3);
        }
    }
}
