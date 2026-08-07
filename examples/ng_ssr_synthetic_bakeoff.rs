//! **Synthetic STR delimiter bake-off (step 3).** Reads with *known* truth, run through the real
//! delimiters, designed to populate the four buckets the real-data investigation found:
//!
//!   * **both OK**   — clean, in-frame alleles, non-motif flanks: both should recover the length.
//!   * **favours unit_slip** — homopolymer / short-tract expansions, where flat_gap collapses the
//!     measurement onto the reference length (the documented failure).
//!   * **favours flat_gap** — a flank that *continues the motif* by one base: unit_slip tends to
//!     grab the half-unit into the tract (the +1 "boundary convention" seen at period 2).
//!   * **breaks both** — interrupted / out-of-frame tracts, and long alleles that outrun the read
//!     (partials, the regime a per-read oracle cannot reach on real data).
//!
//! Unlike the real-data oracle, here the truth is *injected*: we build the read from a chosen
//! allele, so its exact tract length and boundary are known, and we can score both delimiters
//! against it — reaching the long-tract regime and defining the boundary unambiguously.
//!
//! ```text
//! ng_ssr_synthetic_bakeoff            # TSV of every scenario to stdout, a summary to stderr
//! ```
//!
//! It calls `BestPathAligner::align` directly: the delimitation *is* the diverging component, and
//! everything around it (read fetch, CIGAR windowing, the tally) is shared between the two aligners
//! and cannot cause divergence. No BAM, no reference file, no I/O.

use std::process::ExitCode;

use pop_var_caller::ng::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner;
use pop_var_caller::ng::alignment::ssr_best_path_unit_slip::SsrUnitSlipAligner;
use pop_var_caller::ng::alignment::{
    PerQualityEmission, ReadBases, RepeatContext, RepeatGeometry, RepeatSpan, StutterModel,
};
use pop_var_caller::ng::locus_generation::ssr::RepeatDelimiter;
use pop_var_caller::ng::types::{Bp, Motif};

/// Two flank bodies (aperiodic, no long homopolymer at the junction). The junction base is
/// overwritten per scenario so it either **breaks** the motif (clean) or **continues** it (the
/// boundary trap).
const BODY_L: &[u8] = b"GATCTTGCAAGCTGGAATCCGTTAC";
const BODY_R: &[u8] = b"CAGTTCACGATCCTAAGGCTTGACG";

/// How many bases of flank each scenario lays down either side of the tract — **the
/// aligner's anchor**, and settable so the bake-off can answer how short it may be.
///
/// It matters because `flank_bp` is capped by the bundle radius
/// (`SsrGeneratorConfig::check_flank_within`), so narrowing the radius to keep more tracts
/// nameable as loci also shortens every read's anchor. ng's default moved 30 → 20 on
/// 2026-08-07 and this is the instrument that says whether that costs recovery: the tract
/// length is *injected*, so a recovered length can be scored against truth — which the
/// real-data survey cannot do, and which is why its "no degradation" reading could not
/// carry the decision on its own.
///
/// Bounded above by the bodies' own 25 bases. The truncation keeps the **junction** end of
/// each body, since that is the base the scenario overwrites to break or continue the motif.
static FLANK_LEN: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(BODY_L.len());

fn flank_len() -> usize {
    FLANK_LEN
        .load(std::sync::atomic::Ordering::Relaxed)
        .min(BODY_L.len())
}

/// The left body, truncated to [`flank_len`] — keeping its **last** bases, so the junction
/// base the caller overwrites stays at the end.
fn body_left() -> Vec<u8> {
    BODY_L[BODY_L.len() - flank_len()..].to_vec()
}

/// The right body, truncated to [`flank_len`] — keeping its **first** bases, so the junction
/// base stays at the start.
fn body_right() -> Vec<u8> {
    BODY_R[..flank_len()].to_vec()
}

/// Primitive motifs, one per period 1..=6.
fn motif_for(period: usize) -> &'static [u8] {
    match period {
        1 => b"A",
        2 => b"AC",
        3 => b"AAG",
        4 => b"AAGG",
        5 => b"AAGGC",
        6 => b"AAGGCC",
        _ => unreachable!(),
    }
}

/// A base not in `avoid` — used to break the motif at a flank junction so the tract boundary is
/// unambiguous.
fn break_base(avoid: &[u8]) -> u8 {
    *b"ACGT".iter().find(|b| !avoid.contains(b)).unwrap()
}

/// Left flank whose last base does not extend the motif leftward.
fn clean_left(motif: &[u8]) -> Vec<u8> {
    let mut f = body_left();
    *f.last_mut().unwrap() = break_base(&[motif[0], motif[motif.len() - 1]]);
    f
}

/// Right flank whose first base does not extend the motif rightward.
fn clean_right(motif: &[u8]) -> Vec<u8> {
    let mut f = body_right();
    f[0] = break_base(&[motif[0], motif[motif.len() - 1]]);
    f
}

/// Right flank that **continues the motif by one base** — the boundary trap: the tract truly ends
/// where we say, but the next base looks like the start of another unit.
fn continuing_right(motif: &[u8]) -> Vec<u8> {
    let mut f = clean_right(motif);
    f[0] = motif[0];
    f
}

fn tract(motif: &[u8], units: usize) -> Vec<u8> {
    motif
        .iter()
        .copied()
        .cycle()
        .take(units * motif.len())
        .collect()
}

/// A synthetic scenario: a reference frame, a read built from a chosen allele, and the truth.
///
/// `period`, `ref_tract_len` and `allele` describe the case rather than drive it — they say what
/// was built, which is what makes a scenario readable at the construction site, and the TSV keys
/// on `family` instead. Kept deliberately: a field removed here has to be re-derived by whoever
/// next reads the scenario list.
#[expect(
    dead_code,
    reason = "period/ref_tract_len/allele describe the scenario at its construction site"
)]
struct Scenario {
    family: &'static str,
    period: usize,
    motif: Vec<u8>,
    frame: Vec<u8>,
    left_flank_len: usize,
    right_flank_len: usize,
    read: Vec<u8>,
    ref_tract_len: usize,
    /// The measurement a correct delimiter must return for a **complete** read; ignored for partials.
    truth_len: usize,
    /// `true` = the read spans the tract (expect a length); `false` = it runs off (expect a partial).
    complete: bool,
    /// A short human tag for the allele, e.g. `+2u`, `+1bp-oof`, `interrupt`.
    allele: String,
}

impl Scenario {
    fn geometry(&self) -> RepeatGeometry {
        RepeatGeometry {
            left_flank_len: Bp(self.left_flank_len as u64),
            right_flank_len: Bp(self.right_flank_len as u64),
            motif: Motif::new(&self.motif).unwrap(),
        }
    }
}

/// Build a complete-read scenario: `read = left + sample_tract + right[..keep_right]`, reference
/// frame carries the full `right`. `right` may continue the motif (boundary trap). `keep_right`
/// limits how much right flank the *read* shows — a small value puts the tract near the read's 3′
/// end, the realistic boundary trap (production admits a complete read on ≥5 bp flank).
#[allow(clippy::too_many_arguments)]
fn complete_keep(
    family: &'static str,
    period: usize,
    ref_units: usize,
    sample_tract: Vec<u8>,
    right: Vec<u8>,
    keep_right: usize,
    allele: String,
) -> Scenario {
    let motif = motif_for(period).to_vec();
    let left = clean_left(&motif);
    let ref_tract = tract(&motif, ref_units);
    let mut frame = left.clone();
    frame.extend_from_slice(&ref_tract);
    frame.extend_from_slice(&right);
    let mut read = left.clone();
    read.extend_from_slice(&sample_tract);
    read.extend_from_slice(&right[..keep_right.min(right.len())]);
    Scenario {
        family,
        period,
        motif,
        left_flank_len: left.len(),
        right_flank_len: right.len(),
        ref_tract_len: ref_tract.len(),
        truth_len: sample_tract.len(),
        complete: true,
        allele,
        frame,
        read,
    }
}

/// A complete-read scenario with the full right flank in the read (the easy case).
fn complete(
    family: &'static str,
    period: usize,
    ref_units: usize,
    sample_tract: Vec<u8>,
    right: Vec<u8>,
    allele: String,
) -> Scenario {
    let keep = right.len();
    complete_keep(family, period, ref_units, sample_tract, right, keep, allele)
}

/// Build a partial-read scenario: the sample allele is long, and the read runs off `reach` bases
/// into the tract (no right flank). A correct delimiter must return a *partial*, not a short length.
fn partial_left(period: usize, ref_units: usize, reach: usize) -> Scenario {
    let motif = motif_for(period).to_vec();
    let left = clean_left(&motif);
    let right = clean_right(&motif);
    let ref_tract = tract(&motif, ref_units);
    let long = tract(&motif, ref_units + 40); // a big expansion, well past the read
    let mut frame = left.clone();
    frame.extend_from_slice(&ref_tract);
    frame.extend_from_slice(&right);
    let mut read = left.clone();
    read.extend_from_slice(&long[..reach.min(long.len())]); // runs off inside the tract
    Scenario {
        family: "PARTIAL",
        period,
        motif,
        left_flank_len: left.len(),
        right_flank_len: right.len(),
        ref_tract_len: ref_tract.len(),
        truth_len: 0,
        complete: false,
        allele: format!("reach{reach}"),
        frame,
        read,
    }
}

/// Build a right-partial scenario: the read begins inside a long tract and runs into the right
/// flank (left border never reached). A correct delimiter must return FromRight, not a spurious
/// short complete.
fn partial_right(period: usize, ref_units: usize, reach: usize) -> Scenario {
    let motif = motif_for(period).to_vec();
    let left = clean_left(&motif);
    let right = clean_right(&motif);
    let ref_tract = tract(&motif, ref_units);
    let long = tract(&motif, ref_units + 40);
    let mut frame = left.clone();
    frame.extend_from_slice(&ref_tract);
    frame.extend_from_slice(&right);
    let start = long.len().saturating_sub(reach.min(long.len()));
    let mut read = long[start..].to_vec();
    read.extend_from_slice(&right);
    Scenario {
        family: "PARTIAL",
        period,
        motif,
        left_flank_len: left.len(),
        right_flank_len: right.len(),
        ref_tract_len: ref_tract.len(),
        truth_len: 0,
        complete: false,
        allele: format!("Rreach{reach}"),
        frame,
        read,
    }
}

/// Build a complete-read scenario with a 1 bp indel injected into the READ's flank (not the tract),
/// `side` = which flank, `insert` = true for an inserted base else a deletion. The reference frame is
/// normal and the allele equals the reference, so the correct measurement is the reference tract
/// length despite the flank indel — a delimiter that cannot place a gap in the flank mismeasures.
fn flank_indel(period: usize, ref_units: usize, right_side: bool, insert: bool) -> Scenario {
    let motif = motif_for(period).to_vec();
    let left = clean_left(&motif);
    let right = clean_right(&motif);
    let ref_tract = tract(&motif, ref_units);
    // The reference frame is normal; its flank lengths are what the geometry reports.
    let mut frame = left.clone();
    frame.extend_from_slice(&ref_tract);
    frame.extend_from_slice(&right);
    // The read carries a 1 bp indel in the middle of one flank (away from the junction).
    let edit = |flank: &[u8], at: usize| -> Vec<u8> {
        let mut f = flank.to_vec();
        if insert {
            f.insert(at, b'T');
        } else {
            f.remove(at);
        }
        f
    };
    let (read_left, read_right) = if right_side {
        (left.clone(), edit(&right, right.len() / 2))
    } else {
        (edit(&left, left.len() / 2), right.clone())
    };
    let mut read = read_left;
    read.extend_from_slice(&ref_tract);
    read.extend_from_slice(&read_right);
    Scenario {
        family: "FLANK_INDEL",
        period,
        motif,
        left_flank_len: left.len(),
        right_flank_len: right.len(),
        ref_tract_len: ref_tract.len(),
        truth_len: ref_tract.len(),
        complete: true,
        allele: format!(
            "{}{}",
            if right_side { "R" } else { "L" },
            if insert { "+1" } else { "-1" }
        ),
        frame,
        read,
    }
}

/// Run one delimiter on `read` against the scenario's frame, returning the measured span. `read` is
/// passed explicitly so a noisy replicate can be scored without rebuilding the scenario.
fn measure_read<A: RepeatDelimiter>(
    aligner: &A,
    scenario: &Scenario,
    read: &[u8],
    stutter: &StutterModel,
) -> RepeatSpan {
    let geometry = scenario.geometry();
    let context = RepeatContext {
        geometry: &geometry,
        stutter,
    };
    let quals = vec![40u8; read.len()];
    let bases = ReadBases::try_new(read, &quals).expect("equal length");
    let mut scratch = A::Scratch::default();
    aligner.align(bases, &scenario.frame, context, &mut scratch)
}

fn measure<A: RepeatDelimiter>(
    aligner: &A,
    scenario: &Scenario,
    stutter: &StutterModel,
) -> RepeatSpan {
    measure_read(aligner, scenario, &scenario.read, stutter)
}

/// A tiny deterministic PRNG (SplitMix64-style) so noisy replicates are reproducible without a
/// dependency and without `Math.random`-style nondeterminism.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn other_base(&mut self, b: u8) -> u8 {
        let bases = [b'A', b'C', b'G', b'T'];
        loop {
            let c = bases[(self.next() % 4) as usize];
            if c != b {
                return c;
            }
        }
    }
}

/// Copy `read`, substituting each base with probability `p` (a confident miscall — quality stays
/// high, the stress case most able to shift a boundary). Length is unchanged, so the truth stands.
fn mutate(read: &[u8], p: f64, rng: &mut Rng) -> Vec<u8> {
    read.iter()
        .map(|&b| if rng.unit() < p { rng.other_base(b) } else { b })
        .collect()
}

/// Did this span correctly measure the scenario? For a complete read, the measured length must equal
/// the truth; for a partial, the span must be a lower bound (FromLeft/FromRight), never a spurious
/// `Between` or `Unanchored`.
fn correct(scenario: &Scenario, span: &RepeatSpan) -> bool {
    if scenario.complete {
        span.measured_length() == Some(scenario.truth_len as u64)
    } else {
        matches!(span, RepeatSpan::FromLeft(_) | RepeatSpan::FromRight(_))
    }
}

/// A short rendering of what a delimiter returned, for the TSV.
#[expect(
    dead_code,
    reason = "kept for the TSV column it was written for, which the current output omits"
)]
fn render(span: &RepeatSpan) -> String {
    match span {
        RepeatSpan::Between(r) => format!("len={}", r.end - r.start),
        RepeatSpan::FromLeft(r) => format!("partialL~{}", r.end - r.start),
        RepeatSpan::FromRight(r) => format!("partialR~{}", r.end - r.start),
        RepeatSpan::Unanchored => "unanchored".to_string(),
    }
}

fn scenarios() -> Vec<Scenario> {
    let mut out = Vec::new();
    for period in 1..=6 {
        let motif = motif_for(period).to_vec();
        // ref tract ~ 24 bp so complete reads span comfortably.
        let ref_units = (24 / period).max(3);
        let clean_r = clean_right(&motif);

        // CLEAN — in-frame stutter/allele, non-motif flanks. The baseline (expect both OK).
        for d in [-2i64, -1, 0, 1, 2] {
            let su = ref_units as i64 + d;
            if su < 2 {
                continue;
            }
            out.push(complete(
                "CLEAN",
                period,
                ref_units,
                tract(&motif, su as usize),
                clean_r.clone(),
                format!("{d:+}u"),
            ));
        }

        // HOMO_EXP — bigger expansions at low period, where flat_gap is prone to collapse to ref.
        if period <= 2 {
            for d in [3i64, 4, 6] {
                out.push(complete(
                    "HOMO_EXP",
                    period,
                    ref_units,
                    tract(&motif, (ref_units as i64 + d) as usize),
                    clean_r.clone(),
                    format!("+{d}u"),
                ));
            }
        }

        // BOUNDARY — allele == reference, but the right flank continues the motif by one base.
        // Truth is the reference tract; a delimiter that grabs the extra base overshoots by 1.
        out.push(complete(
            "BOUNDARY",
            period,
            ref_units,
            tract(&motif, ref_units),
            continuing_right(&motif),
            "flank+1base".to_string(),
        ));

        // BOUND_SHORT — the boundary trap with the tract near the read's 3′ end (only 5 bp of a
        // motif-continuing flank in the read). This is the realistic version of the +1 divergence.
        for d in [0i64, 1] {
            let su = (ref_units as i64 + d) as usize;
            out.push(complete_keep(
                "BOUND_SHORT",
                period,
                ref_units,
                tract(&motif, su),
                continuing_right(&motif),
                5,
                format!("{d:+}u,flank5"),
            ));
        }

        // CLEAN_SHORT — a clean (non-continuing) flank, also only 5 bp in the read: isolates whether
        // a short flank alone (no motif continuation) is enough to split the two delimiters.
        out.push(complete_keep(
            "CLEAN_SHORT",
            period,
            ref_units,
            tract(&motif, ref_units),
            clean_right(&motif),
            5,
            "+0u,flank5".to_string(),
        ));

        // INTERRUPT — one interior base substituted (same length, an interrupted repeat).
        if period >= 2 {
            let mut t = tract(&motif, ref_units);
            let mid = t.len() / 2;
            t[mid] = break_base(&[t[mid]]); // a non-motif base in the interior
            out.push(complete(
                "INTERRUPT",
                period,
                ref_units,
                t,
                clean_r.clone(),
                "1subst".to_string(),
            ));
        }

        // OOF — one extra non-motif base inserted mid-tract (length not a multiple of period).
        if period >= 2 {
            let mut t = tract(&motif, ref_units);
            let mid = t.len() / 2;
            t.insert(mid, break_base(&[motif[0], motif[period - 1]]));
            out.push(complete(
                "OOF",
                period,
                ref_units,
                t,
                clean_r.clone(),
                "+1bp".to_string(),
            ));
        }

        // PARTIAL — a long allele that outruns the read; expect a lower-bound span. Sweep the reach
        // through the hard complete-vs-partial decision zone, both directions.
        for reach in [40usize, 55, 70, 85, 100] {
            out.push(partial_left(period, ref_units, reach));
        }
        for reach in [55usize, 85] {
            out.push(partial_right(period, ref_units, reach));
        }

        // FLANK_INDEL — a real 1 bp indel in the read's flank (not the tract). The tract length is
        // unchanged, so the correct measurement is the reference length — but a delimiter that bans
        // flank gaps cannot represent the indel and may mismeasure. Both sides, insertion + deletion.
        for right_side in [false, true] {
            for insert in [false, true] {
                out.push(flank_indel(period, ref_units, right_side, insert));
            }
        }
    }
    out
}

/// A per-aligner scorecard — the fitness function. All entries are fractions correct (higher better).
/// v2 adds the two axes real-data validation exposed: `partial_noise` (a genuine one-anchored partial
/// must STAY a partial under sequencing error, not collapse to none — the failure that halved
/// anchor_firm's real partials), and `flank_indel` (a real indel in the flank must not derail the
/// tract measurement — the risk of banning flank gaps).
struct Card {
    clean: f64,   // error-free complete reads, no flank indel, measured to the exact length
    partial: f64, // error-free long alleles correctly called a lower bound
    partial_noise: f64, // NOISY long alleles still called a lower bound (mean over rates)
    noise: f64,   // noisy complete reads measured to the exact length (mean over rates)
    flank_indel: f64, // complete reads with a 1 bp flank indel measured to the exact tract length
    composite: f64, // 0 if clean < 0.999; else mean of the four axes above
}

const REPS: u32 = 200;
const RATES: [f64; 4] = [0.005, 0.01, 0.02, 0.04];

fn is_indel(s: &Scenario) -> bool {
    s.family == "FLANK_INDEL"
}

/// Score one delimiter. Noise seeds depend only on the scenario/rate (not the aligner), so every
/// aligner faces the **identical** noisy reads — a paired test.
fn evaluate<A: RepeatDelimiter>(
    aligner: &A,
    scenarios: &[Scenario],
    stutter: &StutterModel,
) -> Card {
    let frac = |it: &mut dyn Iterator<Item = bool>| {
        let (mut ok, mut n) = (0u64, 0u64);
        for b in it {
            ok += b as u64;
            n += 1;
        }
        if n == 0 { 1.0 } else { ok as f64 / n as f64 }
    };
    // Pools: the "clean" baseline excludes the flank-indel family (that is its own axis).
    let clean_pool: Vec<&Scenario> = scenarios
        .iter()
        .filter(|s| s.complete && !is_indel(s))
        .collect();
    let partial_pool: Vec<&Scenario> = scenarios.iter().filter(|s| !s.complete).collect();
    let indel_pool: Vec<&Scenario> = scenarios.iter().filter(|s| is_indel(s)).collect();

    let clean = frac(
        &mut clean_pool
            .iter()
            .map(|s| correct(s, &measure(aligner, s, stutter))),
    );
    let partial = frac(
        &mut partial_pool
            .iter()
            .map(|s| correct(s, &measure(aligner, s, stutter))),
    );
    let flank_indel = frac(
        &mut indel_pool
            .iter()
            .map(|s| correct(s, &measure(aligner, s, stutter))),
    );

    // Noise over a pool of scenarios: mutate each read REPS times at each rate, average correctness.
    let noise_over = |pool: &[&Scenario]| -> f64 {
        let mut acc = 0.0;
        for (ri, &p) in RATES.iter().enumerate() {
            let (mut ok, mut n) = (0u64, 0u64);
            for (si, s) in pool.iter().enumerate() {
                let mut rng = Rng(0xD1B54A32D192ED03 ^ ((ri as u64) << 40) ^ ((si as u64) << 8));
                for _ in 0..REPS {
                    let read = mutate(&s.read, p, &mut rng);
                    ok += correct(s, &measure_read(aligner, s, &read, stutter)) as u64;
                    n += 1;
                }
            }
            acc += ok as f64 / n.max(1) as f64;
        }
        acc / RATES.len() as f64
    };
    let noise = noise_over(&clean_pool);
    let partial_noise = noise_over(&partial_pool);

    let composite = if clean < 0.999 {
        0.0
    } else {
        (partial + partial_noise + noise + flank_indel) / 4.0
    };
    Card {
        clean,
        partial,
        partial_noise,
        noise,
        flank_indel,
        composite,
    }
}

fn main() -> ExitCode {
    // `--flank N` shortens the anchor either side of every tract, so the same scenarios can
    // be scored at each candidate `flank_bp`. Capped at the flank bodies' own 25 bases.
    let mut rest = std::env::args().skip(1);
    while let Some(arg) = rest.next() {
        if arg == "--flank" {
            match rest.next().and_then(|v| v.parse::<usize>().ok()) {
                Some(bases) if (1..=BODY_L.len()).contains(&bases) => {
                    FLANK_LEN.store(bases, std::sync::atomic::Ordering::Relaxed);
                }
                _ => {
                    eprintln!("error: --flank needs 1..={} bases", BODY_L.len());
                    return ExitCode::from(2);
                }
            }
        }
    }
    eprintln!("flank: {} bp either side of the tract", flank_len());

    let flat = SsrFlatGapAligner::new(PerQualityEmission::new());
    let unit = SsrUnitSlipAligner::new(PerQualityEmission::new());
    let stutter = StutterModel::hipstr_shipped();
    let scen = scenarios();

    // A new algorithm registers one line here: (name, evaluate(&aligner, &scen, &stutter)).
    // Winner = highest composite with clean == 1.000. Composite is the mean of partial,
    // partial_noise, noise, and flank_indel — the four axes that matter on real data.
    let anchor_firm = pop_var_caller::ng::alignment::ssr_anchor_firm::SsrAnchorFirmAligner::new(
        PerQualityEmission::new(),
    );
    let noise_robust = pop_var_caller::ng::alignment::ssr_noise_robust::SsrNoiseRobustAligner::new(
        PerQualityEmission::new(),
    );
    let robust_indel = pop_var_caller::ng::alignment::ssr_robust_indel::SsrRobustIndelAligner::new(
        PerQualityEmission::new(),
    );
    let unit_robust = pop_var_caller::ng::alignment::ssr_unit_robust::SsrUnitRobustAligner::new(
        PerQualityEmission::new(),
    );
    let anchor_robust =
        pop_var_caller::ng::alignment::ssr_anchor_robust::SsrAnchorRobustAligner::new(
            PerQualityEmission::new(),
        );
    let cards: Vec<(&str, Card)> = vec![
        ("flat_gap (algo 3)", evaluate(&flat, &scen, &stutter)),
        ("unit_slip (algo 4)", evaluate(&unit, &scen, &stutter)),
        (
            "anchor_firm (algo 5)",
            evaluate(&anchor_firm, &scen, &stutter),
        ),
        (
            "noise_robust (algo 4n)",
            evaluate(&noise_robust, &scen, &stutter),
        ),
        (
            "robust_indel (algo 4r)",
            evaluate(&robust_indel, &scen, &stutter),
        ),
        (
            "unit_robust (algo 4u)",
            evaluate(&unit_robust, &scen, &stutter),
        ),
        (
            "anchor_robust (algo 5r)",
            evaluate(&anchor_robust, &scen, &stutter),
        ),
    ];

    println!(
        "{:<23} {:>6} {:>8} {:>9} {:>7} {:>11} {:>10}",
        "delimiter", "clean", "partial", "p_noise", "noise", "flank_indel", "COMPOSITE"
    );
    for (name, c) in &cards {
        println!(
            "{:<23} {:>6.3} {:>8.3} {:>9.3} {:>7.3} {:>11.3} {:>10.4}",
            name, c.clean, c.partial, c.partial_noise, c.noise, c.flank_indel, c.composite
        );
    }
    println!(
        "\nclean must stay 1.000 (error-free, no-indel complete reads) or COMPOSITE is zeroed.\n\
         partial/p_noise = long alleles kept as a lower bound (error-free / noisy).\n\
         noise = noisy complete reads measured exactly. flank_indel = a 1 bp flank indel tolerated.\n\
         COMPOSITE = mean(partial, p_noise, noise, flank_indel). Scenarios: {} ({} clean, {} partial, {} flank-indel).",
        scen.len(),
        scen.iter().filter(|s| s.complete && !is_indel(s)).count(),
        scen.iter().filter(|s| !s.complete).count(),
        scen.iter().filter(|s| is_indel(s)).count(),
    );
    ExitCode::SUCCESS
}
