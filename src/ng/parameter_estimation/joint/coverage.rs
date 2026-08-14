//! Building one sample's coverage-by-window summary as the walk runs.
//!
//! [`CoverageByWindow`](super::census::CoverageByWindow) is the finished object and lives
//! with the records it travels beside; this module is what fills it. Design:
//! `doc/devel/ng/spec/parameter_prepass_joint_records.md` §4.
//!
//! **Why a window and not the records' own depths.** The fit's third class of site — a locus
//! the *sample* carries more copies of than the reference does, collecting two copies' reads
//! at one position — is told apart by the coverage *around* it and never by its own depth.
//! The records hold one binned depth per kept position, and the kept positions are one in a
//! few hundred, so a 500 bp window holds one or two of them. The summary is therefore over
//! **every position the walk had in scope**, kept or not.
//!
//! # The two passes, and why the denominator is not the walk's
//!
//! A window's mean depth is `Σ depth / positions in scope`, and the denominator comes from
//! the **reference and the analysed regions**, not from what the walk emitted. A window whose
//! reads are missing has a low mean depth, which is the truth about it; counting only the
//! positions a read reached would make it read normal. So:
//!
//! 1. [`CoverageAccumulator::observe_reference`] is handed the reference bases of each
//!    analysed span. It fills every window's denominator and its GC content.
//! 2. [`CoverageAccumulator::add_depths`] is handed the walk's per-position depths.
//!
//! # Which windows exist
//!
//! **The grid is a function of the reference and the analysed regions alone**, so two samples'
//! summaries are comparable by construction and neither stores a coordinate — the same
//! property, and the same reason, as the kept positions of
//! [`loci`](super::loci). [`CoverageGrid`] is that derivation, and the fit rebuilds it from
//! inputs it already holds.
//!
//! **Windows the analysed regions do not touch are not windows.** A run restricted to a BED
//! touches a small share of a reference's grid, and giving the untouched ones an entry would
//! put every one of them in
//! [`CoverageByWindow`](super::census::CoverageByWindow)'s short-window list — the list whose
//! whole design assumes it is sparse — and then weight them, at mean depth zero, into any
//! wider mean summed across them.

use crate::ng::parameter_estimation::joint::census::CoverageByWindow;
use crate::ng::parameter_estimation::joint::loci::SelectableRegions;
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};

/// The stored window width. **Fine on purpose**: a window's mean separates one copy from two
/// only once it has collected about 12,000 aligned bases, which is 500 bp at 25× and 5 kb at
/// 2.5×, and summing adjacent windows back to a wider mean is exact while splitting one is
/// impossible (`parameter_prepass_joint_records.md` §4.1).
pub const DEFAULT_WINDOW_BP: Bp = Bp(500);

/// GC bins for the depth-against-GC curve: two percentage points wide.
pub const GC_BINS: usize = 51;

/// Below this many windows a GC bin's own median is scatter, and the whole sample's median
/// stands in for it.
const MIN_WINDOWS_PER_GC_BIN: usize = 50;

/// How much of a window must be in scope before its mean depth is allowed to set the
/// sample's median or a GC bin's.
///
/// A window cut down to a handful of positions has a mean, and that mean is noise; letting it
/// into the median moves the scale every relative coverage is read against.
const MIN_IN_SCOPE_SHARE: f64 = 0.5;

/// Which windows exist, and where a position falls among them.
///
/// Derived from the analysed regions and the window width — never stored, and rebuilt by the
/// fit from the same two inputs the identity already makes every sample agree on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoverageGrid {
    window_bp: Bp,
    /// `(contig, window ordinal)` for every window an analysed span touches, sorted. The
    /// slot a window occupies in the summary is its index here.
    windows: Vec<(u32, u32)>,
}

impl CoverageGrid {
    /// The windows `analysed` touches, in reference order.
    ///
    /// # Panics
    ///
    /// When `window_bp` is zero — a grid of zero-width windows has no positions in it.
    pub fn over(analysed: &SelectableRegions, window_bp: Bp) -> Self {
        assert!(window_bp.get() > 0, "a window has to be at least one base");
        let width = window_bp.get();
        let mut windows = Vec::new();
        for span in analysed.spans() {
            let first = span.start.get() / width;
            let last = span.end.get() / width;
            for ordinal in first..=last {
                windows.push((span.contig.get(), ordinal as u32));
            }
        }
        // Two spans may share a boundary window, and the spans arrive sorted but a contig's
        // may not be contiguous, so both the sort and the dedup do work.
        windows.sort_unstable();
        windows.dedup();
        Self { window_bp, windows }
    }

    pub fn window_bp(&self) -> Bp {
        self.window_bp
    }

    pub fn len(&self) -> usize {
        self.windows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Which slot holds this position, or `None` where the analysed regions never reached it.
    pub fn slot_of(&self, contig: ContigId, position: Position) -> Option<usize> {
        let key = (contig.get(), (position.get() / self.window_bp.get()) as u32);
        self.windows.binary_search(&key).ok()
    }

    /// The reference span this slot covers — its whole window, not the part in scope.
    pub fn span_of(&self, slot: usize) -> GenomeRegion {
        let (contig, ordinal) = self.windows[slot];
        let start = u64::from(ordinal) * self.window_bp.get();
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(start + self.window_bp.get() - 1),
        }
    }

    /// Whether slot `b` is the window immediately after slot `a` on the same contig.
    ///
    /// **What summing adjacent windows may cross.** A run over a BED holds windows that are
    /// neighbours in the summary and megabases apart in the genome, and a wider mean taken
    /// across that pair is a mean over two unrelated stretches.
    pub fn adjacent(&self, a: usize, b: usize) -> bool {
        let (contig_a, ordinal_a) = self.windows[a];
        let (contig_b, ordinal_b) = self.windows[b];
        contig_a == contig_b && ordinal_b == ordinal_a + 1
    }
}

/// One sample's window sums, before they are scaled and packed.
///
/// Fourteen bytes a window while the walk runs, against the finished summary's one: 22 MB on
/// tomato's whole grid and 87 MB on GRCh38's, held by one sample at a time and dropped when
/// it finishes.
pub struct CoverageAccumulator {
    grid: CoverageGrid,
    /// Σ observation depth over the positions in scope.
    depth_sum: Vec<u64>,
    /// Positions in scope — the denominator, from the reference and not from the walk.
    in_scope: Vec<u16>,
    /// Reference `G` or `C`, and `A` or `T`, over the positions in scope.
    gc: Vec<u16>,
    at: Vec<u16>,
}

impl CoverageAccumulator {
    pub fn new(grid: CoverageGrid) -> Self {
        let windows = grid.len();
        Self {
            grid,
            depth_sum: vec![0; windows],
            in_scope: vec![0; windows],
            gc: vec![0; windows],
            at: vec![0; windows],
        }
    }

    pub fn grid(&self) -> &CoverageGrid {
        &self.grid
    }

    /// The reference bases of one analysed span, in order from `region.start`.
    ///
    /// Fills the denominators and the GC counts. Call it once for each analysed span; calling
    /// it twice for the same span counts that span twice.
    ///
    /// # Panics
    ///
    /// When `bases` is not as long as `region`, which would silently shorten a denominator.
    pub fn observe_reference(&mut self, region: GenomeRegion, bases: &[u8]) {
        assert_eq!(
            bases.len() as u64,
            region.len(),
            "a span's bases and its coordinates describe the same stretch"
        );
        for (offset, base) in bases.iter().enumerate() {
            let position = Position(region.start.get() + offset as u64);
            let Some(slot) = self.grid.slot_of(region.contig, position) else {
                continue;
            };
            self.in_scope[slot] = self.in_scope[slot].saturating_add(1);
            match base.to_ascii_uppercase() {
                b'G' | b'C' => self.gc[slot] = self.gc[slot].saturating_add(1),
                b'A' | b'T' => self.at[slot] = self.at[slot].saturating_add(1),
                _ => {}
            }
        }
    }

    /// The walk's depths at consecutive positions from `start`.
    ///
    /// Positions the grid does not hold are dropped — a locus may extend past the analysed
    /// regions that admitted it, and its overhang belongs to no window.
    pub fn add_depths(&mut self, contig: ContigId, start: Position, depths: &[u32]) {
        for (offset, depth) in depths.iter().enumerate() {
            let position = Position(start.get() + offset as u64);
            if let Some(slot) = self.grid.slot_of(contig, position) {
                self.depth_sum[slot] += u64::from(*depth);
            }
        }
    }

    /// This sample's median window depth, over the windows at least half in scope.
    ///
    /// **The one number every relative coverage is read against**, which is why a window cut
    /// down to a sliver is not allowed to set it.
    pub fn median_depth(&self) -> f32 {
        let mut depths: Vec<f32> = self.usable_slots().map(|slot| self.mean(slot)).collect();
        median(&mut depths)
    }

    /// Depth against GC content: the median window depth in each GC bin, falling back to the
    /// sample's own median where a bin holds too few windows to have one.
    ///
    /// **On tomato this curve spans a factor of 1.79** — median window depth runs from 16.2
    /// reads a position at 20% GC to 29.0 at 36% — which is larger than the doubling the
    /// duplicated-site class looks for, so a window near an extreme of GC reads high for a
    /// reason that has nothing to do with copy number.
    pub fn gc_curve(&self) -> Vec<f32> {
        let overall = self.median_depth();
        let mut per_bin: Vec<Vec<f32>> = vec![Vec::new(); GC_BINS];
        for slot in self.usable_slots() {
            if let Some(fraction) = self.gc_fraction(slot) {
                per_bin[gc_bin(fraction)].push(self.mean(slot));
            }
        }
        per_bin
            .iter_mut()
            .map(|depths| {
                if depths.len() >= MIN_WINDOWS_PER_GC_BIN {
                    median(depths)
                } else {
                    overall
                }
            })
            .collect()
    }

    /// This window's mean depth over the positions in scope, in reads a position.
    pub fn mean(&self, slot: usize) -> f32 {
        if self.in_scope[slot] == 0 {
            0.0
        } else {
            (self.depth_sum[slot] as f64 / f64::from(self.in_scope[slot])) as f32
        }
    }

    /// The reference GC fraction over this window's positions in scope, where it has one.
    pub fn gc_fraction(&self, slot: usize) -> Option<f32> {
        let called = u32::from(self.gc[slot]) + u32::from(self.at[slot]);
        (called > 0).then(|| f64::from(self.gc[slot]) as f32 / called as f32)
    }

    /// How many positions of this window the analysed regions hold.
    pub fn in_scope(&self, slot: usize) -> u16 {
        self.in_scope[slot]
    }

    /// The windows whose mean depth is allowed to set a scale.
    fn usable_slots(&self) -> impl Iterator<Item = usize> + '_ {
        let floor = (MIN_IN_SCOPE_SHARE * self.grid.window_bp.get() as f64).ceil() as u32;
        (0..self.grid.len()).filter(move |&slot| u32::from(self.in_scope[slot]) >= floor)
    }

    /// The packed summary the records travel with.
    pub fn finish(self) -> CoverageByWindow {
        let median = self.median_depth();
        let curve = self.gc_curve();
        let means: Vec<f32> = (0..self.grid.len()).map(|slot| self.mean(slot)).collect();
        CoverageByWindow::new(self.grid.window_bp, median, &means, &self.in_scope, curve)
    }
}

fn gc_bin(fraction: f32) -> usize {
    ((f64::from(fraction) * (GC_BINS - 1) as f64).round() as usize).min(GC_BINS - 1)
}

/// The middle value, or zero where there is none. Sorts in place, which is why it takes a
/// mutable slice rather than pretending not to.
fn median(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable_by(|a, b| a.partial_cmp(b).expect("no NaN among window depths"));
    values[values.len() / 2]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn regions(spans: &[(u32, u64, u64)]) -> SelectableRegions {
        SelectableRegions::new(
            spans
                .iter()
                .map(|&(contig, start, end)| GenomeRegion {
                    contig: ContigId(contig),
                    start: Position(start),
                    end: Position(end),
                })
                .collect(),
        )
        .expect("the test's spans are disjoint")
    }

    // ---- the grid ------------------------------------------------------------------

    #[test]
    fn the_grid_holds_only_the_windows_the_analysed_regions_touch() {
        // Two spans a long way apart on one contig: four windows, not the 2,001 between
        // them. **This is what keeps the short-window list sparse** — every untouched
        // window would otherwise land in it at zero positions.
        let grid = CoverageGrid::over(&regions(&[(0, 0, 999), (0, 1_000_000, 1_000_999)]), Bp(500));
        assert_eq!(grid.len(), 4);
        assert_eq!(grid.slot_of(ContigId(0), Position(0)), Some(0));
        assert_eq!(grid.slot_of(ContigId(0), Position(999)), Some(1));
        assert_eq!(grid.slot_of(ContigId(0), Position(500_000)), None);
        assert_eq!(grid.slot_of(ContigId(0), Position(1_000_000)), Some(2));
    }

    #[test]
    fn two_spans_sharing_a_boundary_window_share_its_slot() {
        let grid = CoverageGrid::over(&regions(&[(0, 0, 600), (0, 601, 1_200)]), Bp(500));
        assert_eq!(grid.len(), 3);
        assert_eq!(
            grid.slot_of(ContigId(0), Position(600)),
            grid.slot_of(ContigId(0), Position(601)),
            "both spans reach into the window at 500–999"
        );
    }

    #[test]
    fn adjacency_is_genomic_and_not_positional() {
        let grid = CoverageGrid::over(&regions(&[(0, 0, 999), (0, 1_000_000, 1_000_999)]), Bp(500));
        assert!(grid.adjacent(0, 1), "windows 0 and 1 abut in the genome");
        assert!(
            !grid.adjacent(1, 2),
            "neighbouring slots a megabase apart are not one wider window"
        );
    }

    // ---- what the summary is, and what it is not -------------------------------------

    /// Spec §7.10: **a window's mean depth is not the mean over the kept positions inside
    /// it**, and an implementation that quietly derived one from the other would pass every
    /// other check here.
    #[test]
    fn a_windows_mean_is_over_every_position_and_not_over_a_chosen_one() {
        let grid = CoverageGrid::over(&regions(&[(0, 0, 99)]), Bp(100));
        let mut accumulator = CoverageAccumulator::new(grid);
        accumulator.observe_reference(
            GenomeRegion {
                contig: ContigId(0),
                start: Position(0),
                end: Position(99),
            },
            &[b'A'; 100],
        );
        // One position deep, ninety-nine shallow — the shape a kept position landing on a
        // pile-up produces.
        let mut depths = vec![1_u32; 100];
        depths[7] = 101;
        accumulator.add_depths(ContigId(0), Position(0), &depths);

        assert!((accumulator.mean(0) - 2.0).abs() < 1e-6, "(99 + 101) / 100");
        assert_eq!(
            depths[7], 101,
            "the position a kept locus would have reported"
        );
    }

    #[test]
    fn the_denominator_is_the_reference_and_not_what_the_walk_reached() {
        let grid = CoverageGrid::over(&regions(&[(0, 0, 99)]), Bp(100));
        let mut accumulator = CoverageAccumulator::new(grid);
        accumulator.observe_reference(
            GenomeRegion {
                contig: ContigId(0),
                start: Position(0),
                end: Position(99),
            },
            &[b'A'; 100],
        );
        // Ten positions carrying ten reads each; the other ninety were never reached.
        accumulator.add_depths(ContigId(0), Position(0), &[10; 10]);
        assert!(
            (accumulator.mean(0) - 1.0).abs() < 1e-6,
            "100 reads over 100 positions in scope, not over the 10 that were covered"
        );
    }

    // ---- summing back to a wider window (spec §7.11) ----------------------------------

    #[test]
    fn ten_windows_summed_equal_one_walk_over_the_same_five_kilobases() {
        let span = GenomeRegion {
            contig: ContigId(0),
            start: Position(0),
            end: Position(4_999),
        };
        let grid = CoverageGrid::over(&regions(&[(0, 0, 4_999)]), Bp(500));
        let mut fine = CoverageAccumulator::new(grid);
        fine.observe_reference(span, &[b'A'; 5_000]);
        let depths: Vec<u32> = (0..5_000).map(|i| (i % 17) as u32).collect();
        fine.add_depths(ContigId(0), Position(0), &depths);
        let summary = fine.finish();

        let direct = depths.iter().map(|d| f64::from(*d)).sum::<f64>() / 5_000.0;
        let summed = f64::from(summary.mean_depth_over(0..10));
        assert!(
            (summed - direct).abs() < 0.02,
            "ten 500 bp windows summed give {summed}, one 5 kb walk gives {direct}"
        );
    }

    /// Spec §7.11: **plant a short window** — one the analysed regions cut down — and assert
    /// that weighting it by its own position count changes the answer. A sum that treated
    /// every window as full would agree everywhere else and be wrong only at contig ends and
    /// region edges.
    #[test]
    fn a_short_window_is_weighted_by_its_own_positions() {
        // Two windows: a full one at depth 1, and a 50-base one at depth 7. Seven and not
        // seventy — the stored byte is `round(32 × mean / median)`, so a window above eight
        // times the sample's own median saturates, and a test that read the saturated value
        // back would be measuring the ceiling rather than the weighting.
        let grid = CoverageGrid::over(&regions(&[(0, 0, 549)]), Bp(500));
        let mut accumulator = CoverageAccumulator::new(grid);
        accumulator.observe_reference(
            GenomeRegion {
                contig: ContigId(0),
                start: Position(0),
                end: Position(549),
            },
            &[b'A'; 550],
        );
        accumulator.add_depths(ContigId(0), Position(0), &[1; 500]);
        accumulator.add_depths(ContigId(0), Position(500), &[7; 50]);
        let summary = accumulator.finish();

        assert_eq!(summary.positions(0), 500);
        assert_eq!(summary.positions(1), 50, "the second window is a stub");

        let weighted = summary.mean_depth_over(0..2);
        let as_if_full = (summary.mean_depth(0) + summary.mean_depth(1)) / 2.0;
        assert!(
            (weighted - (500.0 + 350.0) / 550.0).abs() < 0.02,
            "850 reads over 550 positions, got {weighted}"
        );
        assert!(
            (as_if_full - weighted).abs() > 2.0,
            "treating the stub as full gives {as_if_full} against {weighted}; if these agreed \
             the test would not be watching anything"
        );
    }

    // ---- the GC curve ------------------------------------------------------------------

    #[test]
    fn the_gc_curve_falls_back_where_a_bin_is_thin() {
        // Sixty windows at 50% GC and one at 0%: the populated bin gets its own median, the
        // lone window's bin gets the sample's.
        let grid = CoverageGrid::over(&regions(&[(0, 0, 61 * 100 - 1)]), Bp(100));
        let mut accumulator = CoverageAccumulator::new(grid);
        for window in 0..61_u64 {
            let start = window * 100;
            let bases: Vec<u8> = if window == 60 {
                vec![b'A'; 100]
            } else {
                (0..100)
                    .map(|i| if i % 2 == 0 { b'G' } else { b'A' })
                    .collect()
            };
            accumulator.observe_reference(
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(start),
                    end: Position(start + 99),
                },
                &bases,
            );
            accumulator.add_depths(ContigId(0), Position(start), &[10; 100]);
        }
        let curve = accumulator.gc_curve();
        assert!((curve[gc_bin(0.5)] - 10.0).abs() < 1e-6);
        assert!(
            (curve[gc_bin(0.0)] - accumulator.median_depth()).abs() < 1e-6,
            "one window is not a median"
        );
    }

    #[test]
    fn a_window_barely_in_scope_does_not_set_the_median() {
        // Ninety-nine full windows at depth 4, and one 10-base sliver at depth 400.
        let mut spans: Vec<(u32, u64, u64)> =
            (0..99).map(|w| (0, w * 500, w * 500 + 499)).collect();
        spans.push((0, 99 * 500, 99 * 500 + 9));
        let grid = CoverageGrid::over(&regions(&spans), Bp(500));
        let mut accumulator = CoverageAccumulator::new(grid);
        for window in 0..99_u64 {
            let start = window * 500;
            accumulator.observe_reference(
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(start),
                    end: Position(start + 499),
                },
                &[b'A'; 500],
            );
            accumulator.add_depths(ContigId(0), Position(start), &[4; 500]);
        }
        accumulator.observe_reference(
            GenomeRegion {
                contig: ContigId(0),
                start: Position(99 * 500),
                end: Position(99 * 500 + 9),
            },
            &[b'A'; 10],
        );
        accumulator.add_depths(ContigId(0), Position(99 * 500), &[400; 10]);

        assert!((accumulator.median_depth() - 4.0).abs() < 1e-6);
        assert!(
            (accumulator.mean(99) - 400.0).abs() < 1e-6,
            "the sliver keeps its own mean; it is only barred from setting the scale"
        );
    }

    // ---- what travels with the summary --------------------------------------------------

    #[test]
    fn the_stored_width_travels_and_a_different_one_is_a_different_summary() {
        let span = GenomeRegion {
            contig: ContigId(0),
            start: Position(0),
            end: Position(4_999),
        };
        let build = |width: Bp| {
            let mut accumulator =
                CoverageAccumulator::new(CoverageGrid::over(&regions(&[(0, 0, 4_999)]), width));
            accumulator.observe_reference(span, &[b'A'; 5_000]);
            accumulator.add_depths(ContigId(0), Position(0), &[3; 5_000]);
            accumulator.finish()
        };
        let fine = build(Bp(500));
        let coarse = build(Bp(5_000));
        assert_eq!(fine.window_bp(), Bp(500));
        assert_eq!(coarse.window_bp(), Bp(5_000));
        assert_ne!(
            fine.window_bp(),
            coarse.window_bp(),
            "the width is what the identity check compares; summaries on two grids are not \
             comparable and are refused rather than resampled"
        );
    }
}
