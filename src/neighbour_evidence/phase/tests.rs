use super::*;
use crate::read_phase::fragments::{FragmentAllele, PackedRecordFragments, RecordFragments};

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
    let swapped = phase.align(&round);
    let expected: Vec<bool> = (0..20).map(|site| (8..14).contains(&site)).collect();
    assert_eq!(swapped, expected);
    // Added, the swapped round agrees with the consensus everywhere.
    phase.add_round(&round);
    for index in 0..20 {
        let (mut first, mut second) = (Vec::new(), Vec::new());
        phase.consensus_into(index, &mut first, &mut second);
        assert!((first[1] - 0.95).abs() < 1e-6 && (second[1] - 0.05).abs() < 1e-6);
    }
}

/// **A single disagreeing site is not worth two changes of state**: it stays "kept".
#[test]
fn the_alignment_keeps_an_isolated_disagreement() {
    let pairs = layout(20);
    let mut phase = SamplePhase::new(&[2; 20], &pairs);
    phase.add_round(&path(20, |_| false));
    assert_eq!(phase.align(&path(20, |site| site == 10)), vec![false; 20]);
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
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
    for (site, written) in phased.iter().enumerate() {
        let written = written.expect("every site is phased");
        let block_start = if site < 10 { 1 } else { 1 + 10 * 1_000 };
        assert_eq!(written.phase_set, block_start, "site {site}");
        assert_eq!((written.copy_0_allele, written.copy_1_allele), (1, 0));
    }
    // With one path of five allowed to disagree, the link holds.
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            PathAgreement::try_new(0.75).expect("a share"),
            &SiteReads::default(),
        )
        .0;
    assert!(
        phased
            .iter()
            .all(|written| written.is_some_and(|w| w.phase_set == 1))
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
        let phased = phase
            .phase_heterozygotes(
                &pairs,
                &every_site_heterozygous(n),
                &positions(n),
                300_000.0,
                PathAgreement::ALL,
                &SiteReads::default(),
            )
            .0;
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
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &heterozygotes,
            &positions,
            30_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
    let blocks: Vec<Option<u64>> = phased.iter().map(|w| w.map(|w| w.phase_set)).collect();
    // Sites 0, 1 and 3 link (no uncalled record breaks a block); 3 → 4 is past S; 4 and 5 link.
    assert_eq!(
        blocks,
        vec![Some(1), Some(1), None, Some(1), Some(50_001), Some(50_001)]
    );
    let heterozygotes = vec![Some((0, 1)), None, None, None, Some((0, 1)), None];
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &heterozygotes,
            &positions,
            30_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
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
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions,
            300_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
    assert_eq!(phased[0].map(|w| w.phase_set), Some(1));
    assert_eq!(phased[1].map(|w| w.phase_set), Some(1));
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
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &[Some((0, 1)), Some((1, 2)), Some((0, 1))],
            &positions(3),
            300_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
    let middle = phased[1].expect("the 1/2 is phased");
    assert_eq!((middle.copy_0_allele, middle.copy_1_allele), (2, 1));
    assert_eq!(
        phased[0].map(|w| (w.copy_0_allele, w.copy_1_allele)),
        Some((1, 0))
    );
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
            copy_0_allele: first,
            copy_1_allele: 1 - first,
            phase_set: block,
        });
    }
    out
}

/// The last written block of phase set 77 that wrote copy 0's allele `allele` at each `(site, allele)`
/// of `written` — biallelic records, so copy 1 carries the other.
fn last_block(written: &[(usize, usize)]) -> LastBlock {
    LastBlock {
        phase_set: 77,
        records: written
            .iter()
            .map(|&(record, allele)| JoinRecord {
                record,
                copy_0_allele: allele,
                copy_1_allele: 1 - allele,
            })
            .collect(),
    }
}

/// [`join_and_name_core_blocks`] on an eight-site window whose core is sites 4–7, records being
/// sites.
fn join(window: &[Option<WrittenPhase>], last: Option<&LastBlock>) -> Vec<Option<WrittenPhase>> {
    join_and_name_core_blocks(
        window,
        |site| site >= 4,
        &positions(8),
        last,
        |record| record,
    )
}

/// **A window joins the last written block where three shared records agree** (spec §2.4): sites
/// 0–3 are the margin, 4–7 the core. The last block wrote copy 0 = allele 1 at sites 1–3; this
/// window agrees — the core keeps the block's `PS` and orientation — or agrees reversed, and the
/// core is swapped as a whole.
#[test]
fn a_window_joins_where_three_shared_records_agree() {
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
    let out = join(&window, Some(&last_block(&written)));
    assert_eq!(out[..4], [None; 4], "the margin is not written");
    assert_eq!(
        out[4],
        Some(WrittenPhase {
            copy_0_allele: 0,
            copy_1_allele: 1,
            phase_set: 77
        })
    );
    assert_eq!(
        out[6],
        Some(WrittenPhase {
            copy_0_allele: 1,
            copy_1_allele: 0,
            phase_set: 77
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
    let out = join(&reversed, Some(&last_block(&written)));
    assert_eq!(
        out[4],
        Some(WrittenPhase {
            copy_0_allele: 1,
            copy_1_allele: 0,
            phase_set: 77
        })
    );
    assert_eq!(
        out[6],
        Some(WrittenPhase {
            copy_0_allele: 0,
            copy_1_allele: 1,
            phase_set: 77
        })
    );
}

/// **No join without three agreeing records in the block that reaches the core**: the windows
/// disagreeing at one record, only two shared records, or the shared records in a block that ends
/// before the core — each core block is then named by its first core record, and a core block of
/// one record is unphased.
#[test]
fn a_window_does_not_join_without_agreement() {
    let written = [(1, 1), (2, 1), (3, 1)];
    let named_by_core = |out: &[Option<WrittenPhase>]| {
        assert_eq!(out[4].map(|phase| phase.phase_set), Some(4_001));
        assert_eq!(out[6].map(|phase| phase.phase_set), Some(4_001));
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
    named_by_core(&join(&mixed, Some(&last_block(&written))));
    let two = window_phase(
        8,
        &[(2, 1, 2_001), (3, 1, 2_001), (4, 0, 2_001), (6, 1, 2_001)],
    );
    named_by_core(&join(&two, Some(&last_block(&written))));
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
    named_by_core(&join(&ends_before, Some(&last_block(&written))));

    // One core record: phased only as part of a joined block.
    let lone = window_phase(
        8,
        &[(1, 1, 1_001), (2, 1, 1_001), (3, 1, 1_001), (5, 0, 1_001)],
    );
    let joined = join(&lone, Some(&last_block(&written)));
    assert_eq!(joined[5].map(|phase| phase.phase_set), Some(77));
    let alone = join(&lone, None);
    assert_eq!(alone[5], None);
}

/// **A link breaks where a path agrees in orientation but does not lean** (spec §2.3, code review
/// 2026-11-02): at site 1 the fifth path leans by 0.6·0.6 − 0.4·0.4 = 0.2, below [`LEAN_MIN`].
#[test]
fn a_link_breaks_where_one_path_does_not_lean() {
    let n = 3;
    let mut paths: Vec<Vec<HaplotypePair>> = (0..4).map(|_| path(n, |_| false)).collect();
    paths.push(vec![site(0.95, 0.05), site(0.6, 0.4), site(0.95, 0.05)]);
    let (pairs, phase) = rounds(n, &paths);
    let phased = phase
        .phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            PathAgreement::ALL,
            &SiteReads::default(),
        )
        .0;
    assert!(phased.iter().all(Option::is_none), "{phased:?}");
}

/// **A share asks for its paths rounded up**: with 3 of 5 paths agreeing, 0.7 (3.5 paths, so 4)
/// breaks the link and 0.5 (2.5, so 3) keeps it.
#[test]
fn a_share_asks_for_whole_paths_rounded_up() {
    let n = 2;
    let pairs = layout(n);
    let mut phase = SamplePhase::new(&[2; 2], &pairs);
    for _ in 0..3 {
        let agreeing = path(n, |_| false);
        phase.record_path(&pairs, &agreeing);
        phase.add_round(&agreeing);
    }
    for _ in 0..2 {
        phase.record_path(&pairs, &path(n, |site| site == 1));
    }
    let at = |share: f64| {
        phase
            .phase_heterozygotes(
                &pairs,
                &every_site_heterozygous(n),
                &positions(n),
                300_000.0,
                PathAgreement::try_new(share).expect("a share"),
                &SiteReads::default(),
            )
            .0
    };
    assert!(
        at(0.7).iter().all(Option::is_none),
        "3 of 5 is short of 0.7"
    );
    assert!(at(0.5).iter().all(Option::is_some), "3 of 5 meets 0.5");
}

/// **A share is a number in `(0, 1]`**.
#[test]
fn a_path_agreement_is_a_share() {
    for share in [0.0, -0.1, 1.000_001, f64::NAN, f64::INFINITY] {
        assert!(PathAgreement::try_new(share).is_none(), "{share}");
    }
    for share in [1.0, 0.5, 1e-3] {
        assert_eq!(
            PathAgreement::try_new(share).map(PathAgreement::get),
            Some(share)
        );
    }
    assert_eq!(PathAgreement::default(), PathAgreement::DEFAULT);
}

/// **A tie between "kept" and "swapped" goes to "kept"**: one agreeing and one disagreeing site,
/// each state costing 1.
#[test]
fn the_alignment_resolves_a_tie_to_kept() {
    let pairs = layout(2);
    let mut phase = SamplePhase::new(&[2; 2], &pairs);
    phase.add_round(&path(2, |_| false));
    assert_eq!(phase.align(&path(2, |site| site == 1)), vec![false, false]);
}

/// **The alignment ignores sites where the round does not lean**: a stretch where the round is
/// swapped by only 0.4 is not marked, so it is not swapped back.
#[test]
fn the_alignment_ignores_sites_where_the_round_does_not_lean() {
    let pairs = layout(20);
    let mut phase = SamplePhase::new(&[2; 20], &pairs);
    phase.add_round(&path(20, |_| false));
    let round: Vec<HaplotypePair> = (0..20)
        .map(|index| {
            if (8..14).contains(&index) {
                site(0.3, 0.7)
            } else {
                site(0.95, 0.05)
            }
        })
        .collect();
    assert_eq!(phase.align(&round), vec![false; 20]);
}

/// **A join reads the shared records nearest the core**: the farthest of four disagrees, the
/// three nearest agree, and the window joins.
#[test]
fn a_join_reads_the_records_nearest_the_core() {
    let written = [(0, 0), (1, 1), (2, 1), (3, 1)];
    let window = window_phase(
        8,
        &[
            (0, 1, 1),
            (1, 1, 1),
            (2, 1, 1),
            (3, 1, 1),
            (4, 0, 1),
            (6, 1, 1),
        ],
    );
    let out = join(&window, Some(&last_block(&written)));
    assert_eq!(out[4].map(|phase| phase.phase_set), Some(77));
}

/// **A join reads the block of the window's first phased core record**, not of its first phased
/// record: site 0 is in a block of its own in this window, and the core continues the last block.
#[test]
fn a_join_reads_the_block_of_the_first_core_record() {
    let written = [(0, 1), (1, 1), (2, 1), (3, 1)];
    let window = window_phase(
        8,
        &[
            (0, 1, 1),
            (1, 1, 1_001),
            (2, 1, 1_001),
            (3, 1, 1_001),
            (4, 0, 1_001),
            (6, 1, 1_001),
        ],
    );
    let out = join(&window, Some(&last_block(&written)));
    assert_eq!(out[4].map(|phase| phase.phase_set), Some(77));
}

// ---------------------------------------------------------------- reads in the link test

/// Fragments at each site, as `(fragment, allele)` at Phred 30, and every site's written `GQ`.
fn read_fixture(
    n: usize,
    at: &[(usize, &[(u32, u16)])],
    quality: impl Fn(usize) -> Option<u32>,
) -> (Vec<PackedRecordFragments>, Vec<Option<u32>>) {
    let mut fragments = vec![PackedRecordFragments::default(); n];
    for &(site, entries) in at {
        fragments[site] = RecordFragments::from_entries(
            entries
                .iter()
                .map(|&(fragment, allele)| FragmentAllele {
                    sample: 0,
                    fragment,
                    allele,
                    error_phred: 30,
                })
                .collect(),
        )
        .into();
    }
    (fragments, (0..n).map(quality).collect())
}

fn site_reads<'a>(
    fragments: &'a [PackedRecordFragments],
    qualities: &[Option<u32>],
) -> SiteReads<'a> {
    SiteReads {
        fragments: fragments.iter().map(|at| at.of_sample(0)).collect(),
        genotype_qualities: qualities.to_vec(),
        drawn_without_reads: Vec::new(),
    }
}

/// Each site's written copy-0 allele and phase set, `None` where unphased.
fn written(phased: &[Option<WrittenPhase>]) -> Vec<Option<(usize, u64)>> {
    phased
        .iter()
        .map(|written| written.map(|w| (w.copy_0_allele, w.phase_set)))
        .collect()
}

/// **A kept population link stands where the reads are silent or agree, and breaks where they
/// contradict it at ln 9 or more** (read phase spec §2.7, rows 1 and 2); the read check counts the
/// contradicted link among those reads reached. Five paths agree everywhere, so every link is kept;
/// two fragments say sites 1 and 2 are trans where the consensus has them cis.
#[test]
fn reads_break_a_population_link_they_contradict_and_leave_the_rest() {
    let n = 4;
    let paths: Vec<Vec<HaplotypePair>> = (0..5).map(|_| path(n, |_| false)).collect();
    let (pairs, phase) = rounds(n, &paths);
    let (fragments, qualities) =
        read_fixture(n, &[(1, &[(7, 0), (8, 1)]), (2, &[(7, 1), (8, 0)])], |_| {
            Some(40)
        });
    let (phased, counts) = phase.phase_heterozygotes(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        300_000.0,
        PathAgreement::ALL,
        &site_reads(&fragments, &qualities),
    );
    assert_eq!(
        written(&phased),
        vec![
            Some((1, 1)),
            Some((1, 1)),
            Some((1, 2_001)),
            Some((1, 2_001))
        ]
    );
    assert_eq!(
        counts,
        PassReadCounts {
            population_links_reached: 1,
            population_links_contradicted: 1,
            population_links_broken: 1,
            links_made_by_reads: 0,
            records_attached_by_reads: 0,
        }
    );
}

/// **Where the draw used reads, the read check counts only a link with an end held out of it**
/// (reads-in-draw spec §2.6). The fixture of the test above, with two fragments per link at every
/// pair of consecutive sites: cis at 0–1 and 2–3, trans at 1–2. Every site drawn with reads, the
/// check counts nothing, while the link test still breaks 1–2; site 2 held out, it counts 1–2 and
/// 2–3; site 0 held out, only 0–1.
#[test]
fn the_read_check_counts_only_links_with_an_end_held_out_of_the_draw() {
    let n = 4;
    let paths: Vec<Vec<HaplotypePair>> = (0..5).map(|_| path(n, |_| false)).collect();
    let (pairs, phase) = rounds(n, &paths);
    let (fragments, qualities) = read_fixture(
        n,
        &[
            (0, &[(1, 0), (2, 1)]),
            (1, &[(1, 0), (2, 1), (7, 0), (8, 1)]),
            (2, &[(7, 1), (8, 0), (9, 1), (10, 0)]),
            (3, &[(9, 1), (10, 0)]),
        ],
        |_| Some(40),
    );
    let counts_with = |held_out: Option<usize>| {
        let mut reads = site_reads(&fragments, &qualities);
        reads.drawn_without_reads = (0..n).map(|site| Some(site) == held_out).collect();
        phase
            .phase_heterozygotes(
                &pairs,
                &every_site_heterozygous(n),
                &positions(n),
                300_000.0,
                PathAgreement::ALL,
                &reads,
            )
            .1
    };
    let reached_and_contradicted = |counts: PassReadCounts| {
        (
            counts.population_links_reached,
            counts.population_links_contradicted,
        )
    };
    let none_held_out = counts_with(None);
    assert_eq!(reached_and_contradicted(none_held_out), (0, 0));
    assert_eq!(
        none_held_out.population_links_broken, 1,
        "the link test still acts"
    );
    assert_eq!(reached_and_contradicted(counts_with(Some(2))), (2, 1));
    assert_eq!(reached_and_contradicted(counts_with(Some(0))), (1, 0));
}

/// **Reads link two candidates the population did not, and a link whose relation is the opposite
/// of the consensus's swaps the following block as a whole** (§2.7, row 3). One path of five turns
/// over at site 10, so the population keeps no link across 9–10; two fragments say trans there,
/// where the consensus has the two cis. Sites 10 on are written the other way round, in site 0's
/// phase set. With `GQ` 29 at site 10 (row 4) the reads link nothing.
#[test]
fn reads_link_what_the_population_did_not_and_swap_the_block_after() {
    let n = 20;
    let mut paths: Vec<Vec<HaplotypePair>> = (0..4).map(|_| path(n, |_| false)).collect();
    paths.push(path(n, |site| site >= 10));
    let (pairs, phase) = rounds(n, &paths);
    let at: &[(usize, &[(u32, u16)])] = &[(9, &[(3, 0), (4, 1)]), (10, &[(3, 1), (4, 0)])];

    let (fragments, qualities) = read_fixture(n, at, |_| Some(40));
    let (phased, counts) = phase.phase_heterozygotes(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        300_000.0,
        PathAgreement::ALL,
        &site_reads(&fragments, &qualities),
    );
    let expected: Vec<Option<(usize, u64)>> = (0..n)
        .map(|site| Some((usize::from(site < 10), 1)))
        .collect();
    assert_eq!(written(&phased), expected);
    assert_eq!(counts.links_made_by_reads, 1);

    let (fragments, qualities) = read_fixture(n, at, |site| Some(if site == 10 { 29 } else { 40 }));
    let (phased, counts) = phase.phase_heterozygotes(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        300_000.0,
        PathAgreement::ALL,
        &site_reads(&fragments, &qualities),
    );
    assert_eq!(phased[10].map(|w| w.phase_set), Some(10_001));
    assert_eq!(counts.links_made_by_reads, 0);
}

/// Five sites; the consensus leans at 0, 3 and 4 and not at 1 and 2, every path agreeing.
fn with_two_records_the_consensus_does_not_lean_on() -> (PairLayout, SamplePhase) {
    let n = 5;
    let path: Vec<HaplotypePair> = (0..n)
        .map(|site| {
            if site == 1 || site == 2 {
                self::site(0.5, 0.5)
            } else {
                self::site(0.95, 0.05)
            }
        })
        .collect();
    rounds(n, &vec![path; 5])
}

/// **Heterozygotes the consensus does not lean on are attached by reads, as a chain from the
/// candidate before them** (§2.7): 1 is cis with 0 and 2 trans with 1, so 1 is written as 0 is
/// and 2 the other way round, all in 0's phase set; the chain's last record says nothing against 3.
/// When two fragments say 2 is cis with 3, which the block has trans, the chain is dropped; when
/// only one fragment links 1 and 2, the chain stops at 1.
#[test]
fn a_chain_of_records_the_consensus_does_not_lean_on_is_attached_by_reads() {
    let (pairs, phase) = with_two_records_the_consensus_does_not_lean_on();
    let n = 5;
    let phased_with = |at: &[(usize, &[(u32, u16)])]| {
        let (fragments, qualities) = read_fixture(n, at, |_| Some(40));
        phase.phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            PathAgreement::ALL,
            &site_reads(&fragments, &qualities),
        )
    };

    let (phased, counts) = phased_with(&[
        (0, &[(1, 1), (2, 0)]),
        (1, &[(1, 1), (2, 0), (5, 1), (6, 0)]),
        (2, &[(5, 0), (6, 1)]),
    ]);
    assert_eq!(
        written(&phased),
        vec![
            Some((1, 1)),
            Some((1, 1)),
            Some((0, 1)),
            Some((1, 1)),
            Some((1, 1))
        ]
    );
    assert_eq!(counts.records_attached_by_reads, 2);

    let (phased, counts) = phased_with(&[
        (0, &[(1, 1), (2, 0)]),
        (1, &[(1, 1), (2, 0), (5, 1), (6, 0)]),
        (2, &[(5, 0), (6, 1), (9, 0), (10, 1)]),
        (3, &[(9, 0), (10, 1)]),
    ]);
    assert_eq!(written(&phased)[1..3], [None, None], "{phased:?}");
    assert_eq!(counts.records_attached_by_reads, 0);

    let (phased, counts) = phased_with(&[
        (0, &[(1, 1), (2, 0)]),
        (1, &[(1, 1), (2, 0), (5, 1)]),
        (2, &[(5, 0)]),
    ]);
    assert_eq!(written(&phased)[1..3], [Some((1, 1)), None]);
    assert_eq!(counts.records_attached_by_reads, 1);
}

/// **A heterozygote below `GQ` 30 is never linked by reads**, however well they agree, and one with
/// no `GQ` — a margin record an earlier window wrote — neither: at a wrong heterozygote every
/// fragment shows one allele, and two fragments then agree half the time (§2.5).
#[test]
fn a_heterozygote_below_gq_30_is_never_attached_by_reads() {
    let (pairs, phase) = with_two_records_the_consensus_does_not_lean_on();
    let n = 5;
    for quality in [Some(29), None] {
        let (fragments, qualities) = read_fixture(
            n,
            &[(0, &[(1, 1), (2, 1)]), (1, &[(1, 1), (2, 1)])],
            |site| if site == 1 { quality } else { Some(40) },
        );
        let (phased, counts) = phase.phase_heterozygotes(
            &pairs,
            &every_site_heterozygous(n),
            &positions(n),
            300_000.0,
            PathAgreement::ALL,
            &site_reads(&fragments, &qualities),
        );
        assert_eq!(phased[1], None, "GQ {quality:?}");
        assert_eq!(counts.records_attached_by_reads, 0);
    }
}

// ------------------------------------------- the limits on when reads may phase (review Mi6)

/// The two-records fixture, phased with these fragments, every site at `GQ` 40 unless `quality`
/// says otherwise, and the copying length `segment_bases`.
fn attached_with(
    at: &[(usize, &[(u32, u16)])],
    quality: impl Fn(usize) -> Option<u32>,
    segment_bases: f64,
) -> (Vec<Option<WrittenPhase>>, PassReadCounts) {
    let (pairs, phase) = with_two_records_the_consensus_does_not_lean_on();
    let n = 5;
    let (fragments, qualities) = read_fixture(n, at, quality);
    phase.phase_heterozygotes(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        segment_bases,
        PathAgreement::ALL,
        &site_reads(&fragments, &qualities),
    )
}

/// **Reads make no link beyond the copying length** (§2.7): sites 1,000 bases apart with `S` 500,
/// so the population keeps no link; two good fragments say sites 1 and 2 are cis, and still no
/// link is made. Kills a read-made link without the `S` limit (review mutant p1).
#[test]
fn reads_make_no_link_beyond_the_copying_length() {
    let n = 4;
    let paths: Vec<Vec<HaplotypePair>> = (0..5).map(|_| path(n, |_| false)).collect();
    let (pairs, phase) = rounds(n, &paths);
    let (fragments, qualities) =
        read_fixture(n, &[(1, &[(7, 0), (8, 0)]), (2, &[(7, 0), (8, 0)])], |_| {
            Some(40)
        });
    let (phased, counts) = phase.phase_heterozygotes(
        &pairs,
        &every_site_heterozygous(n),
        &positions(n),
        500.0,
        PathAgreement::ALL,
        &site_reads(&fragments, &qualities),
    );
    assert_eq!(written(&phased), vec![None; n]);
    assert_eq!(counts.links_made_by_reads, 0);
}

/// **A heterozygote beyond the copying length is not attached** (§2.7): with `S` 500, two good
/// fragments over sites 0 and 1, 1,000 bases apart, attach nothing. Kills attaching without the
/// `S` limit (mutant p2).
#[test]
fn a_heterozygote_beyond_the_copying_length_is_not_attached() {
    let (phased, counts) = attached_with(
        &[(0, &[(1, 1), (2, 1)]), (1, &[(1, 1), (2, 1)])],
        |_| Some(40),
        500.0,
    );
    assert_eq!(phased[1], None);
    assert_eq!(counts.records_attached_by_reads, 0);
}

/// **A chain stops at the first heterozygote reads do not link** (§2.7): fragments 1 and 2 cover
/// sites 0 and 2 and not site 1, so the chain from 0 stops at 1 and site 2 is not attached, though
/// the reads link it to 0. Kills the chain stepping over a site it could not link (mutant p3).
#[test]
fn a_chain_stops_at_the_first_heterozygote_reads_do_not_link() {
    let (phased, counts) = attached_with(
        &[(0, &[(1, 1), (2, 1)]), (2, &[(1, 1), (2, 1)])],
        |_| Some(40),
        300_000.0,
    );
    assert_eq!(written(&phased)[1..3], [None, None]);
    assert_eq!(counts.records_attached_by_reads, 0);
}

/// **No chain grows from a candidate below `GQ` 30** (§2.7): site 0 at `GQ` 29, two good
/// fragments over sites 0 and 1; site 1 is not attached. Kills attaching from `p` without its
/// `GQ` check (mutant p4) — every earlier test lowered only the attached site's `GQ`.
#[test]
fn no_chain_grows_from_a_candidate_below_gq_30() {
    let (phased, counts) = attached_with(
        &[(0, &[(1, 1), (2, 1)]), (1, &[(1, 1), (2, 1)])],
        |site| Some(if site == 0 { 29 } else { 40 }),
        300_000.0,
    );
    assert_eq!(phased[1], None);
    assert_eq!(counts.records_attached_by_reads, 0);
}

/// **A chain is tested against `q` only where `p` and `q` share a block** (§2.7): fragments 11 and
/// 12 say sites 0 and 3 are trans where the consensus has them cis, so that population link breaks
/// and 0 and 3 are in two blocks. The chain from 0 — site 1 cis with 0, site 2 trans with 1 — is
/// then kept whatever the reads say between 2 and 3, for both values of `q`'s allele on fragments 9
/// and 10. Kills testing the chain against a `q` in another block (mutant p6), which drops it for
/// one of the two.
#[test]
fn a_chain_is_not_tested_against_a_candidate_in_another_block() {
    for q in [0, 1] {
        let at: &[(usize, &[(u32, u16)])] = &[
            (0, &[(1, 1), (2, 1), (11, 0), (12, 0)]),
            (1, &[(1, 1), (2, 1), (5, 1), (6, 1)]),
            (2, &[(5, 0), (6, 0), (9, 0), (10, 0)]),
            (3, &[(9, q), (10, q), (11, 1), (12, 1)]),
        ];
        let (phased, counts) = attached_with(at, |_| Some(40), 300_000.0);
        assert_eq!(counts.population_links_broken, 1, "q {q}");
        assert!(
            phased[1].is_some() && phased[2].is_some(),
            "q {q}: {phased:?}"
        );
        assert_eq!(counts.records_attached_by_reads, 2, "q {q}");
    }
}
