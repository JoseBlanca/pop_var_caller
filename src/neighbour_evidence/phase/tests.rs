use super::*;

/// A biallelic site's two haplotypes from each copy's probability of the alternative.
fn site(first_alternative: f64, second_alternative: f64) -> HaplotypePair {
    [
        vec![1.0 - first_alternative, first_alternative],
        vec![1.0 - second_alternative, second_alternative],
    ]
}

/// A path over `n` biallelic sites with copy 0 on the alternative, swapped where `swapped(site)`.
fn path(n: usize, swapped: impl Fn(usize) -> bool) -> Vec<HaplotypePair> {
    (0..n)
        .map(|index| {
            if swapped(index) {
                site(0.05, 0.95)
            } else {
                site(0.95, 0.05)
            }
        })
        .collect()
}

fn layout(n: usize) -> PairLayout {
    PairLayout::new(std::iter::repeat_n(2, n))
}

fn every_site_heterozygous(n: usize) -> Vec<Option<(usize, usize)>> {
    vec![Some((0, 1)); n]
}

fn positions(n: usize) -> Vec<u64> {
    (0..n as u64).map(|index| 1 + index * 1_000).collect()
}

/// **A round equal to the consensus with one stretch swapped is swapped back exactly there**
/// (spec §7.4): a stretch of six marked sites inside the window — more than the two changes of
/// state it costs.
#[test]
fn the_alignment_swaps_a_stretch_back_exactly() {
    let pairs = layout(20);
    let mut phase = SamplePhase::new(&[2; 20], &pairs);
    phase.add_round(&path(20, |_| false));
    let round = path(20, |site| (8..14).contains(&site));
    let swapped = phase.alignment(&round);
    let expected: Vec<bool> = (0..20).map(|site| (8..14).contains(&site)).collect();
    assert_eq!(swapped, expected);
    // Added, the swapped round agrees with the consensus everywhere.
    phase.add_round(&round);
    for index in 0..20 {
        let (first, second) = phase.consensus(index);
        assert!((first[1] - 0.95).abs() < 1e-6 && (second[1] - 0.05).abs() < 1e-6);
    }
}

/// **A single disagreeing site is not worth two changes of state**: it stays "kept".
#[test]
fn the_alignment_keeps_an_isolated_disagreement() {
    let pairs = layout(20);
    let mut phase = SamplePhase::new(&[2; 20], &pairs);
    phase.add_round(&path(20, |_| false));
    assert_eq!(
        phase.alignment(&path(20, |site| site == 10)),
        vec![false; 20]
    );
}

/// One round a path, each path's bits recorded and the path added to the consensus.
fn rounds(n: usize, paths: &[Vec<HaplotypePair>]) -> (PairLayout, SamplePhase) {
    let pairs = layout(n);
    let mut phase = SamplePhase::new(&vec![2; n], &pairs);
    for path in paths {
        phase.record_path(&pairs, path);
        phase.add_round(path);
    }
    (pairs, phase)
}

/// **The link test sees a switch in one path that the alignment absorbs** (spec §2.3): four paths
/// agree, one turns over between sites 9 and 10. Lined up, it agrees with the consensus on both
/// sides; read unaligned, it breaks the link there and nowhere else.
#[test]
fn the_link_test_catches_a_switch_the_alignment_absorbs() {
    let n = 20;
    let mut paths: Vec<Vec<HaplotypePair>> = (0..4).map(|_| path(n, |_| false)).collect();
    paths.push(path(n, |site| site >= 10));
    let (pairs, phase) = rounds(n, &paths);
    let phased = phase.phased(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        300_000.0,
        MIN_PATH_AGREEMENT,
    );
    for (site, written) in phased.iter().enumerate() {
        let written = written.expect("every site is phased");
        let block_start = if site < 10 { 1 } else { 1 + 10 * 1_000 };
        assert_eq!(written.block, block_start, "site {site}");
        assert_eq!((written.first, written.second), (1, 0));
    }
    // With one path of five allowed to disagree, the link holds.
    let phased = phase.phased(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        300_000.0,
        0.75,
    );
    assert!(
        phased
            .iter()
            .all(|written| written.is_some_and(|w| w.block == 1))
    );
}

/// **The consensus must lean by at least [`LEAN_MIN`]** (spec §2.2): 0.75 against 0.25 leans by
/// exactly 0.5 and is phased; 0.74 against 0.25 is not.
#[test]
fn the_lean_minimum_is_inclusive() {
    let n = 3;
    let make = |middle: f64| -> Vec<HaplotypePair> {
        vec![site(0.95, 0.05), site(middle, 0.25), site(0.95, 0.05)]
    };
    for (middle, phased_middle) in [(0.75, true), (0.74, false)] {
        let (pairs, phase) = rounds(n, &[make(middle)]);
        let phased = phase.phased(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            MIN_PATH_AGREEMENT,
        );
        assert_eq!(phased[1].is_some(), phased_middle, "middle {middle}");
    }
}

/// **A homozygote, an uncalled record and a block of one are written unphased** (spec §2.2), and a
/// link longer than `S` ends a block.
#[test]
fn only_heterozygotes_in_blocks_of_two_are_phased() {
    let n = 6;
    let (pairs, phase) = rounds(n, &[path(n, |_| false)]);
    let heterozygotes = vec![
        Some((0, 1)),
        Some((0, 1)),
        None,
        Some((0, 1)),
        Some((0, 1)),
        Some((0, 1)),
    ];
    let positions = vec![1, 1_001, 2_001, 3_001, 50_001, 51_001];
    let phased = phase.phased(&pairs, &heterozygotes, &positions, 30_000.0, 1.0);
    let blocks: Vec<Option<u64>> = phased.iter().map(|w| w.map(|w| w.block)).collect();
    // Sites 0, 1 and 3 link (no uncalled record breaks a block); 3 → 4 is past S; 4 and 5 link.
    assert_eq!(
        blocks,
        vec![Some(1), Some(1), None, Some(1), Some(50_001), Some(50_001)]
    );
    let heterozygotes = vec![Some((0, 1)), None, None, None, Some((0, 1)), None];
    let phased = phase.phased(&pairs, &heterozygotes, &positions, 30_000.0, 1.0);
    assert!(phased.iter().all(Option::is_none), "two blocks of one");
}

/// **A broken link between two records at one position leaves the second unphased**, so no two
/// blocks share a `PS` (spec §2.3).
#[test]
fn a_broken_link_at_one_position_leaves_the_second_unphased() {
    let n = 4;
    let mut paths: Vec<Vec<HaplotypePair>> = (0..4).map(|_| path(n, |_| false)).collect();
    paths.push(path(n, |site| site >= 2));
    let (pairs, phase) = rounds(n, &paths);
    let positions = vec![1, 1_001, 1_001, 2_001];
    let phased = phase.phased(
        &pairs,
        &every_site_heterozygous(n),
        &positions,
        300_000.0,
        1.0,
    );
    assert_eq!(phased[0].map(|w| w.block), Some(1));
    assert_eq!(phased[1].map(|w| w.block), Some(1));
    assert_eq!(phased[2], None);
    // Site 3 follows site 1 across the same switch: a new block, which has one record.
    assert_eq!(phased[3], None);
}

/// **A multi-allelic `1|2`**: at three alleles the genotype's pair is `(1, 2)`; copy 0 leaning to
/// allele 2 writes `2|1`.
#[test]
fn a_multi_allelic_pair_is_phased_on_its_own_alleles() {
    let three =
        |first: [f64; 3], second: [f64; 3]| -> HaplotypePair { [first.to_vec(), second.to_vec()] };
    let path = vec![
        site(0.95, 0.05),
        three([0.02, 0.03, 0.95], [0.03, 0.95, 0.02]),
        site(0.95, 0.05),
    ];
    let pairs = PairLayout::new([2, 3, 2].into_iter());
    let mut phase = SamplePhase::new(&[2, 3, 2], &pairs);
    phase.record_path(&pairs, &path);
    phase.add_round(&path);
    let phased = phase.phased(
        &pairs,
        &[Some((0, 1)), Some((1, 2)), Some((0, 1))],
        &positions(3),
        300_000.0,
        1.0,
    );
    let middle = phased[1].expect("the 1/2 is phased");
    assert_eq!((middle.first, middle.second), (2, 1));
    assert_eq!(phased[0].map(|w| (w.first, w.second)), Some((1, 0)));
}

/// **The pair entries of a site are its pairs `(a, b)`, `a < b`, in order**, one after another.
#[test]
fn the_pair_entries_are_in_order() {
    let pairs = PairLayout::new([2, 4, 3].into_iter());
    assert_eq!(pairs.entry(0, 0, 1, 2), 0);
    let four: Vec<usize> = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        .iter()
        .map(|&(a, b)| pairs.entry(1, a, b, 4))
        .collect();
    assert_eq!(four, vec![1, 2, 3, 4, 5, 6]);
    assert_eq!(pairs.entry(2, 1, 2, 3), 9);
    assert_eq!(pairs.entries(), 10);
}

/// A window's phase as [`SamplePhase::phased`] gives it: `(site, copy 0's allele, block)` for the
/// phased sites of a biallelic window.
fn window_phase(n: usize, phased: &[(usize, usize, u64)]) -> Vec<Option<WrittenPhase>> {
    let mut out = vec![None; n];
    for &(site, first, block) in phased {
        out[site] = Some(WrittenPhase {
            first,
            second: 1 - first,
            block,
        });
    }
    out
}

/// **A window joins the last written block where three shared records agree** (spec §2.4): sites
/// 0–3 are the margin, 4–7 the core. The last block wrote copy 0 = allele 1 at sites 1–3; this
/// window agrees — the core keeps the block's `PS` and orientation — or agrees reversed, and the
/// core is swapped as a whole.
#[test]
fn a_window_joins_where_three_shared_records_agree() {
    let positions = positions(8);
    let in_core = |site: usize| site >= 4;
    let written = [(1, 1), (2, 1), (3, 1)];
    let window = window_phase(
        8,
        &[
            (1, 1, 1_001),
            (2, 1, 1_001),
            (3, 1, 1_001),
            (4, 0, 1_001),
            (6, 1, 1_001),
        ],
    );
    let out = written_in_core(&window, in_core, &positions, Some((77, &written)));
    assert_eq!(out[..4], [None; 4], "the margin is not written");
    assert_eq!(
        out[4],
        Some(WrittenPhase {
            first: 0,
            second: 1,
            block: 77
        })
    );
    assert_eq!(
        out[6],
        Some(WrittenPhase {
            first: 1,
            second: 0,
            block: 77
        })
    );

    let reversed = window_phase(
        8,
        &[
            (1, 0, 1_001),
            (2, 0, 1_001),
            (3, 0, 1_001),
            (4, 0, 1_001),
            (6, 1, 1_001),
        ],
    );
    let out = written_in_core(&reversed, in_core, &positions, Some((77, &written)));
    assert_eq!(
        out[4],
        Some(WrittenPhase {
            first: 1,
            second: 0,
            block: 77
        })
    );
    assert_eq!(
        out[6],
        Some(WrittenPhase {
            first: 0,
            second: 1,
            block: 77
        })
    );
}

/// **No join without three agreeing records in the block that reaches the core**: the windows
/// disagreeing at one record, only two shared records, or the shared records in a block that ends
/// before the core — each core block is then named by its first core record, and a core block of
/// one record is unphased.
#[test]
fn a_window_does_not_join_without_agreement() {
    let positions = positions(8);
    let in_core = |site: usize| site >= 4;
    let written = [(1, 1), (2, 1), (3, 1)];
    let named_by_core = |out: &[Option<WrittenPhase>]| {
        assert_eq!(out[4].map(|phase| phase.block), Some(4_001));
        assert_eq!(out[6].map(|phase| phase.block), Some(4_001));
    };
    let mixed = window_phase(
        8,
        &[
            (1, 1, 1_001),
            (2, 0, 1_001),
            (3, 1, 1_001),
            (4, 0, 1_001),
            (6, 1, 1_001),
        ],
    );
    named_by_core(&written_in_core(
        &mixed,
        in_core,
        &positions,
        Some((77, &written)),
    ));
    let two = window_phase(
        8,
        &[(2, 1, 2_001), (3, 1, 2_001), (4, 0, 2_001), (6, 1, 2_001)],
    );
    named_by_core(&written_in_core(
        &two,
        in_core,
        &positions,
        Some((77, &written)),
    ));
    let ends_before = window_phase(
        8,
        &[
            (1, 1, 1_001),
            (2, 1, 1_001),
            (3, 1, 1_001),
            (4, 0, 4_001),
            (6, 1, 4_001),
        ],
    );
    named_by_core(&written_in_core(
        &ends_before,
        in_core,
        &positions,
        Some((77, &written)),
    ));

    // One core record: phased only as part of a joined block.
    let lone = window_phase(
        8,
        &[(1, 1, 1_001), (2, 1, 1_001), (3, 1, 1_001), (5, 0, 1_001)],
    );
    let joined = written_in_core(&lone, in_core, &positions, Some((77, &written)));
    assert_eq!(joined[5].map(|phase| phase.block), Some(77));
    let alone = written_in_core(&lone, in_core, &positions, None);
    assert_eq!(alone[5], None);
}
