//! **A read's misread rate after a long run of one base, in the direction it was sequenced**
//! (`doc/devel/ng/spec/str_context_error.md`, `str_context_error_fit.md`).
//!
//! A read sequenced through a long run of one base misreads the bases after it: most often the
//! next base, read as one more base of the run, and every base for at least 20 after it somewhat
//! more than elsewhere. A [`MisreadTable`] holds, per run length and distance after the run, the
//! **excess**: the misread rate in reads sequenced through the run minus the rate in reads reaching
//! the same base from the other side. It is added to the read group's fitted rate, because the
//! fitted rate is the background both groups of reads share and the excess is what the run adds on
//! top of it.
//!
//! **Which entry of the table a base belongs to** is one rule ([`contexts_along`]), used by the
//! parameter pre-pass on the reference to count and by calling on the candidate to charge, so the
//! two agree on what an entry means (`str_context_error_fit.md` §2.3).
//!
//! **Calling reads the context off the candidate**, with the locus's flanks around it: the question
//! is how likely this read is to have misread the candidate, and the candidate is what the read was
//! copied from.

use crate::float;
use crate::types::ErrorRate;

/// The direction a read was sequenced in, on the reference: a forward read from left to right,
/// a reverse read from right to left.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadDirection {
    Forward,
    Reverse,
}

/// How many bases after a run its excess reaches. The measurement stopped at 20, flat to the end.
pub const BASES_CHARGED_AFTER_A_RUN: usize = 20;

/// The shortest run, in bases, whose following bases are charged. Below it the measured excess is
/// under 1 in 1,000 (`str_context_error.md` §2).
pub const SHORTEST_RUN_CHARGED: usize = 8;

/// The shortest run of the open top length class at the next base: 23 bases or more.
const LONGEST_NEXT_BASE_CLASS: usize = 23;

/// How many entries a [`MisreadTable`] has: 16 run lengths at the next base (8 to 22, then 23 or
/// more) and 3 length classes further on (8 to 11, 12 to 15, 16 or more).
pub const ENTRY_COUNT: usize = 16 + 3;

/// **One entry of the table**: a range of run lengths in bases, and whether the base is the next
/// one after the run or one of the 19 after that. Indexes [`MisreadTable`]'s entries, in the order
/// [`RunEntry::all`] lists them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RunEntry(u8);

/// Where a base sits after the run its entry names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceAfterRun {
    /// The base right after the run's last base.
    NextBase,
    /// One of the 2nd to 20th bases after it.
    FurtherOn,
}

impl RunEntry {
    /// The entry for a run of `length` bases, `distance` bases before the base (1 for the next
    /// base). `None` for a run too short to be charged or a distance past
    /// [`BASES_CHARGED_AFTER_A_RUN`].
    #[must_use]
    pub fn of(length: usize, distance: usize) -> Option<Self> {
        if length < SHORTEST_RUN_CHARGED || !(1..=BASES_CHARGED_AFTER_A_RUN).contains(&distance) {
            return None;
        }
        let index = if distance == 1 {
            length.min(LONGEST_NEXT_BASE_CLASS) - SHORTEST_RUN_CHARGED
        } else {
            16 + match length {
                8..=11 => 0,
                12..=15 => 1,
                _ => 2,
            }
        };
        Some(Self(index as u8))
    }

    /// Every entry, in table order.
    pub fn all() -> impl Iterator<Item = Self> {
        (0..ENTRY_COUNT as u8).map(Self)
    }

    /// The entry's place in the table.
    #[must_use]
    pub fn index(self) -> usize {
        usize::from(self.0)
    }

    /// The entry at a place in the table; `None` past its end.
    #[must_use]
    pub fn at(index: usize) -> Option<Self> {
        (index < ENTRY_COUNT).then_some(Self(index as u8))
    }

    /// Whether the entry is the next base after the run or one further on.
    #[must_use]
    pub fn distance(self) -> DistanceAfterRun {
        if self.0 < 16 {
            DistanceAfterRun::NextBase
        } else {
            DistanceAfterRun::FurtherOn
        }
    }

    /// The run lengths, in bases, the entry covers: the shortest, and the longest or `None` for an
    /// open top class.
    #[must_use]
    pub fn lengths(self) -> (usize, Option<usize>) {
        let i = usize::from(self.0);
        if i < 16 {
            let length = i + SHORTEST_RUN_CHARGED;
            return (length, (length < LONGEST_NEXT_BASE_CLASS).then_some(length));
        }
        [(8, Some(11)), (12, Some(15)), (16, None)][i - 16]
    }
}

/// **What reads that crossed a run misread more than reads that did not**, at one entry: the
/// chance of reading the base as one more base of the run, and the chance of reading any other
/// wrong base. Further-on entries keep their whole excess in `other_wrong`, since the measurement
/// found no favourite wrong base there.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MisreadExcess {
    pub as_the_run: f64,
    pub other_wrong: f64,
}

/// **One read group's misread excess after runs**, one [`MisreadExcess`] per [`RunEntry`].
#[derive(Debug, Clone, PartialEq)]
pub struct MisreadTable {
    entries: [MisreadExcess; ENTRY_COUNT],
}

impl MisreadTable {
    /// A table charging nothing anywhere: every base scores at the fitted rate alone.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            entries: [MisreadExcess::default(); ENTRY_COUNT],
        }
    }

    /// A table from its entries, in [`RunEntry::all`]'s order.
    #[must_use]
    pub fn from_entries(entries: [MisreadExcess; ENTRY_COUNT]) -> Self {
        Self { entries }
    }

    /// **The table a read group gets without a fit**: what the parameter pre-pass fits on HG002
    /// chromosome 20 at 30× (`str_context_error_fit.md` §5), regenerated from that fit's
    /// parameters file on 2026-11-02. One Illumina library's numbers, so a read group of another
    /// platform that falls back to them is charged a guess, and its provenance says so.
    #[must_use]
    pub fn defaults() -> Self {
        // (as the run, other wrong), in table order.
        const FITTED_ON_HG002_CHR20_30X: [(f64, f64); ENTRY_COUNT] = [
            (0.002073, 0.000000), // 8 bases, next base
            (0.003439, 0.000000), // 9 bases, next base
            (0.007784, 0.000921), // 10 bases, next base
            (0.012721, 0.001714), // 11 bases, next base
            (0.020489, 0.002889), // 12 bases, next base
            (0.027402, 0.002647), // 13 bases, next base
            (0.037959, 0.004247), // 14 bases, next base
            (0.049556, 0.006237), // 15 bases, next base
            (0.063910, 0.004949), // 16 bases, next base
            (0.072944, 0.006976), // 17 bases, next base
            (0.092829, 0.008675), // 18 bases, next base
            (0.105203, 0.007809), // 19 bases, next base
            (0.105110, 0.007294), // 20 bases, next base
            (0.121211, 0.009084), // 21 bases, next base
            (0.123169, 0.010322), // 22 bases, next base
            (0.117377, 0.009209), // 23+ bases, next base
            (0.000000, 0.001927), // 8–11 bases, 2nd–20th base
            (0.000000, 0.013675), // 12–15 bases, 2nd–20th base
            (0.000000, 0.038420), // 16+ bases, 2nd–20th base
        ];
        let mut entries = [MisreadExcess::default(); ENTRY_COUNT];
        for (entry, (as_the_run, other_wrong)) in entries.iter_mut().zip(FITTED_ON_HG002_CHR20_30X)
        {
            *entry = MisreadExcess {
                as_the_run,
                other_wrong,
            };
        }
        Self { entries }
    }

    /// The excess at one entry.
    #[must_use]
    pub fn at(&self, entry: RunEntry) -> MisreadExcess {
        self.entries[entry.index()]
    }

    /// Every entry's excess, in table order.
    #[must_use]
    pub fn entries(&self) -> &[MisreadExcess; ENTRY_COUNT] {
        &self.entries
    }

    /// Whether every entry is zero, so scoring under it is the flat comparison.
    #[must_use]
    pub fn charges_nothing(&self) -> bool {
        self.entries
            .iter()
            .all(|excess| excess.as_the_run == 0.0 && excess.other_wrong == 0.0)
    }
}

/// [`MisreadTable::defaults`], built once: what a scoring context points at where no read group's
/// own table applies.
pub static DEFAULT_MISREAD_TABLE: std::sync::LazyLock<MisreadTable> =
    std::sync::LazyLock::new(MisreadTable::defaults);

/// **The entry a base belongs to, and the base of the run before it** — what reading the base "as
/// the run" means at a next-base entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContextAfterRun {
    pub entry: RunEntry,
    pub run_base: u8,
}

/// **Each base's context, walking `met` in the order a read met it** — `out[i]` for `met[i]`
/// (`str_context_error_fit.md` §2.3):
///
/// - a run is a maximal stretch of one base, `A`, `C`, `G` or `T`;
/// - a base's candidates are the runs of at least [`SHORTEST_RUN_CHARGED`] bases that ended within
///   the last [`BASES_CHARGED_AFTER_A_RUN`] bases before it;
/// - of those, the longer wins, then the nearer.
///
/// A run cut off by the start of `met` is measured from there, so a caller wanting the length of
/// every run that can reach a base must hand over at least [`BASES_CHARGED_AFTER_A_RUN`] plus 23
/// bases before it, the shortest run of the top length class.
pub fn contexts_along(met: &[u8], out: &mut Vec<Option<ContextAfterRun>>) {
    out.clear();
    out.resize(met.len(), None);
    // The winning candidate's (length, nearness) per base, compared as a tuple.
    let mut winning: Vec<(usize, isize)> = vec![(0, 0); met.len()];
    let mut start = 0;
    while start < met.len() {
        let mut last = start;
        while last + 1 < met.len() && met[last + 1] == met[start] {
            last += 1;
        }
        let length = last - start + 1;
        // A stretch of `N` is no sequence at all, so no run.
        if length >= SHORTEST_RUN_CHARGED && matches!(met[start], b'A' | b'C' | b'G' | b'T') {
            for distance in 1..=BASES_CHARGED_AFTER_A_RUN {
                let at = last + distance;
                if at >= met.len() {
                    break;
                }
                let rank = (length, -(distance as isize));
                if out[at].is_some() && winning[at] >= rank {
                    continue;
                }
                winning[at] = rank;
                out[at] = Some(ContextAfterRun {
                    entry: RunEntry::of(length, distance).expect("a charged length"),
                    run_base: met[start],
                });
            }
        }
        start = last + 1;
    }
}

/// **Whether a run in `bases` is long enough, or one base short of long enough, to be charged,
/// under a table that charges anything** — what decides whether a candidate is scored once per
/// read direction at all: a one-base slip can lengthen such a run to a charged length in the
/// stretching a read is compared with. A table of zeros costs no second scoring.
///
/// `bases` is walked once, so a caller can hand over several slices chained without joining them.
#[must_use]
pub fn a_run_may_reach_a_charged_length<'b>(
    bases: impl IntoIterator<Item = &'b u8>,
    table: &MisreadTable,
) -> bool {
    if table.charges_nothing() {
        return false;
    }
    let (mut previous, mut length) = (None, 0_usize);
    for &base in bases {
        if previous == Some(base) {
            length += 1;
        } else {
            (previous, length) = (Some(base), 1);
        }
        if length + 1 >= SHORTEST_RUN_CHARGED {
            return true;
        }
    }
    false
}

/// **How probable `observed`'s letters are against `candidate`**, both over the same bases, for a
/// read sequenced in `direction`, with `left_flank` and `right_flank` the bases around them, and
/// `table` the read group's misread excess. Bases no charged run precedes are scored at `rate`
/// exactly as the flat comparison scores them.
///
/// `None` when the two differ in length, which the flat comparison answers with zero.
#[must_use]
pub fn letters_after_runs(
    observed: &[u8],
    candidate: &[u8],
    (left_flank, right_flank): (&[u8], &[u8]),
    direction: ReadDirection,
    rate: ErrorRate,
    table: &MisreadTable,
) -> Option<f64> {
    if observed.len() != candidate.len() {
        return None;
    }
    let epsilon = rate.get();
    let match_ln = float::ln(1.0 - epsilon);
    let mismatch_ln = float::ln(epsilon / 3.0);
    // The sequence in the order the read met it, and where the candidate starts in it.
    let (met, first): (Vec<u8>, usize) = match direction {
        ReadDirection::Forward => ([left_flank, candidate].concat(), left_flank.len()),
        ReadDirection::Reverse => (
            right_flank
                .iter()
                .rev()
                .chain(candidate.iter().rev())
                .copied()
                .collect(),
            right_flank.len(),
        ),
    };
    let mut contexts = Vec::new();
    contexts_along(&met, &mut contexts);
    let mut total_ln = 0.0;
    for (at, context) in contexts.iter().enumerate().skip(first) {
        let along = at - first;
        let (candidate_base, read_base) = match direction {
            ReadDirection::Forward => (candidate[along], observed[along]),
            ReadDirection::Reverse => {
                let back = candidate.len() - 1 - along;
                (candidate[back], observed[back])
            }
        };
        total_ln += match context {
            None => {
                if read_base == candidate_base {
                    match_ln
                } else {
                    mismatch_ln
                }
            }
            Some(context) => {
                let excess = table.at(context.entry);
                match context.entry.distance() {
                    DistanceAfterRun::NextBase => {
                        if read_base == candidate_base {
                            float::ln(
                                (1.0 - epsilon - excess.as_the_run - excess.other_wrong)
                                    .max(epsilon),
                            )
                        } else if read_base == context.run_base {
                            float::ln(epsilon / 3.0 + excess.as_the_run)
                        } else {
                            float::ln(epsilon / 3.0 + excess.other_wrong / 2.0)
                        }
                    }
                    DistanceAfterRun::FurtherOn => {
                        if read_base == candidate_base {
                            float::ln((1.0 - epsilon - excess.other_wrong).max(epsilon))
                        } else {
                            float::ln((epsilon + excess.other_wrong) / 3.0)
                        }
                    }
                }
            }
        };
    }
    Some(float::exp(total_ln))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rate() -> ErrorRate {
        ErrorRate::try_new(0.002).unwrap()
    }

    fn contexts(met: &[u8]) -> Vec<Option<ContextAfterRun>> {
        let mut out = Vec::new();
        contexts_along(met, &mut out);
        out
    }

    /// Where no charged run comes first, the comparison is the flat one, bit for bit.
    #[test]
    fn without_a_long_run_it_is_the_flat_comparison() {
        let candidate = b"ACGTACGTAAAACGT";
        let observed = b"ACGTACCTAAAACGT";
        let flat = float::exp(14.0 * float::ln(0.998) + float::ln(0.002 / 3.0));
        for direction in [ReadDirection::Forward, ReadDirection::Reverse] {
            let scored = letters_after_runs(
                observed,
                candidate,
                (b"GGC", b"TCA"),
                direction,
                rate(),
                &MisreadTable::defaults(),
            )
            .unwrap();
            assert!(
                (scored - flat).abs() <= 1e-12 * flat,
                "{scored} against {flat}"
            );
        }
    }

    /// A16 G A10, the `G` read as `A`: about a tenth for a forward read, which crossed the sixteen
    /// `A`s first, and much less for a reverse read, which crossed only the ten.
    #[test]
    fn the_base_after_a_long_run_is_cheap_to_misread_only_from_the_run_s_side() {
        let candidate = [&[b'A'; 16][..], b"G", &[b'A'; 10]].concat();
        let observed = [b'A'; 27];
        let score = |direction| {
            letters_after_runs(
                &observed,
                &candidate,
                (b"", b""),
                direction,
                rate(),
                &MisreadTable::defaults(),
            )
            .unwrap()
        };
        let (forward, reverse) = (score(ReadDirection::Forward), score(ReadDirection::Reverse));
        // The misread: ε/3 + 0.0639 after sixteen against ε/3 + 0.0078 after ten, about 7.6-fold.
        // The bases after the `G` are charged too: the forward read's ten `A`s at 0.0384 each
        // (0.9596¹⁰ = 0.66), the reverse read's sixteen at 0.0019 (0.94) — so about 5.3.
        let ratio = forward / reverse;
        assert!(ratio > 4.8 && ratio < 5.8, "{ratio}");
        let flat = float::exp(26.0 * float::ln(0.998) + float::ln(0.002 / 3.0));
        assert!(
            reverse > flat * 5.0,
            "reverse {reverse} against flat {flat}"
        );
    }

    /// A table charging nothing scores every read exactly as the flat comparison, even after a
    /// long run — what a fitted table of zeros must give.
    #[test]
    fn a_table_of_zeros_is_the_flat_comparison() {
        let candidate = [&[b'A'; 16][..], b"G", &[b'A'; 10]].concat();
        let observed = [b'A'; 27];
        let flat = float::exp(26.0 * float::ln(0.998) + float::ln(0.002 / 3.0));
        for direction in [ReadDirection::Forward, ReadDirection::Reverse] {
            let scored = letters_after_runs(
                &observed,
                &candidate,
                (b"", b""),
                direction,
                rate(),
                &MisreadTable::zero(),
            )
            .unwrap();
            assert!(
                (scored - flat).abs() <= 1e-12 * flat,
                "{scored} against {flat}"
            );
        }
    }

    /// After `GT`×10 nothing is charged: a repeat of two bases is not a run.
    #[test]
    fn a_dinucleotide_repeat_is_not_a_run() {
        assert!(
            contexts(&[&b"GT".repeat(10)[..], b"CAAC"].concat())
                .iter()
                .all(Option::is_none)
        );
    }

    #[test]
    fn a_base_s_entry_follows_the_rule() {
        // A run of 12 then `C`: the `C` is the next base after a 12-base run, the run's base `A`.
        let found = contexts(b"AAAAAAAAAAAACGT");
        assert_eq!(
            found[12],
            Some(ContextAfterRun {
                entry: RunEntry::of(12, 1).unwrap(),
                run_base: b'A',
            })
        );
        assert_eq!(found[13].unwrap().entry, RunEntry::of(12, 2).unwrap());
        assert!(found[..12].iter().all(Option::is_none));
        // Seven bases are not a charged run, and neither are eight `N`s.
        assert!(contexts(b"AAAAAAACGT").iter().all(Option::is_none));
        assert!(contexts(b"NNNNNNNNCGT").iter().all(Option::is_none));
    }

    /// The longest run wins over a nearer, shorter one: eleven bases after `A`×16, a `T`×8 just
    /// ended charges less than the sixteen still do.
    #[test]
    fn the_longest_run_wins_then_the_nearest() {
        let met = [&[b'A'; 16][..], b"CG", &[b'T'; 8][..], b"C"].concat();
        let found = contexts(&met);
        let last = found[met.len() - 1].unwrap();
        assert_eq!(last.entry, RunEntry::of(16, 11).unwrap());
        // Two runs of eight: the nearer one is the next-base entry.
        let met = [&[b'A'; 8][..], b"C", &[b'T'; 8][..], b"G"].concat();
        let last = contexts(&met)[met.len() - 1].unwrap();
        assert_eq!(last.entry, RunEntry::of(8, 1).unwrap());
        assert_eq!(last.run_base, b'T');
    }

    #[test]
    fn every_entry_round_trips_through_its_description() {
        let mut seen = 0;
        for entry in RunEntry::all() {
            let (shortest, _) = entry.lengths();
            let distance = match entry.distance() {
                DistanceAfterRun::NextBase => 1,
                DistanceAfterRun::FurtherOn => 2,
            };
            assert_eq!(RunEntry::of(shortest, distance), Some(entry));
            seen += 1;
        }
        assert_eq!(seen, ENTRY_COUNT);
        assert_eq!(RunEntry::of(7, 1), None);
        assert_eq!(RunEntry::of(16, 21), None);
        assert_eq!(RunEntry::of(40, 1), RunEntry::of(23, 1));
    }

    #[test]
    fn a_run_one_base_short_of_charged_asks_for_the_second_scoring() {
        let table = MisreadTable::defaults();
        assert!(!a_run_may_reach_a_charged_length(b"AAAAAACGT", &table));
        assert!(a_run_may_reach_a_charged_length(b"AAAAAAACGT", &table));
        assert!(!a_run_may_reach_a_charged_length(
            b"GTGTGTGTGTGTGTC",
            &table
        ));
        assert!(!a_run_may_reach_a_charged_length(
            b"AAAAAAAAAAAAC",
            &MisreadTable::zero()
        ));
    }
}
