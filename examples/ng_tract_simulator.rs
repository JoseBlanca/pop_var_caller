//! **Repeat tracts whose genotypes we chose, sequenced by a slippage we set** — the
//! exact-truth ground the tract QUAL experiment needs and no benchmark can give.
//!
//! ```text
//! ng_tract_simulator <out-dir> [key=value ...]
//! ```
//!
//! # The question this fixture exists to answer
//!
//! `doc/devel/ng/spec/calling_loop_ssr.md` §3.3 names one risk for the site quality a caller
//! writes at a repeat tract: the number leans entirely on the stutter model pricing slip
//! products correctly, and **where the model under-prices them the quality grows confidently
//! wrong**, because every extra read compounds the same mispriced error. On a benchmark that
//! risk cannot be isolated — the reads' true slippage is unknown, and the caller's own model is
//! the only description of it anybody has. Here both are ours: the reads are drawn under a
//! stutter model this program is told, and the caller can be handed that same model or the
//! shipped default. **The gap between the two is the experiment**, and it is settable.
//!
//! # What it writes
//!
//! Into `<out-dir>`, everything a run and its scoring need:
//!
//! - `reference.fa` — one contig of repeat tracts, each wrapped in flanks that carry no repeat
//!   of their own, so the only tandem repeats in the file are the ones put there.
//! - `<sample>.bam` — one file a sample, reads spanning each tract with flank either side.
//!   **Not indexed here**: `samtools index` is the driver's job, so this program writes no
//!   dependency on a tool it cannot check for.
//! - `truth.vcf` — the genotypes drawn, as a VCF: one record per tract that any sample carries
//!   a non-reference allele at, and **no record at a tract every sample is homozygous
//!   reference at**. That absence is the point: a call there is a false positive, and without
//!   tracts that carry nothing there is nothing for a calibration to measure.
//! - `tracts.bed` — `chrom start end period`, the tract ground the scorer restricts to.
//! - `confident.bed` — the whole contig, so the scorer's two masks have the same shape here as
//!   on a real benchmark.
//! - `truth_genotypes.tsv` — every sample's repeat counts at every tract, including the
//!   homozygous-reference ones, for anything that wants the genotype rather than the record.
//! - `truth_haplotypes.tsv` — every sample's two haplotypes at every tract, spelled over the
//!   tract and 40 bases either side, for a scorer that compares sequences whatever records the
//!   caller writes them as.
//! - `slippage_rows.toml` — the stutter model the reads were drawn under, written as the rows
//!   a parameters file states slippage in. Append it to a run's own `.parameters.toml`, with
//!   that file's empty `slippage_by_stratum_and_group` line deleted, and the caller scores
//!   these reads under the model that made them instead of under the shipped default. **That
//!   is the only way the experiment's fitted-against-defaulted split has two sides**, because
//!   no command fits a parameters file yet.
//!
//! # The forward model, in the order it is applied
//!
//! Per sample and tract a diploid genotype is drawn: with probability `variant_rate` the
//! sample carries at least one non-reference allele, whose repeat count sits within
//! `allele_span` of the reference tract's. Then each read picks one of the two alleles, and:
//!
//! 1. **slips** — with probability `slip_share` the read reports a length other than its
//!    allele's. The slip is shorter with probability `shorter_share`, and its size in whole
//!    repeats is one with probability `one_step_share`, two with `(1 - one_step_share) *
//!    one_step_share`, and so on. **`one_step_share` is the share of the slips that moved by
//!    exactly one step, not a per-step multiplier** — the same spelling
//!    [`StutterModel`](pop_var_caller::alignment::StutterModel) uses, and inverting it
//!    makes large slips common while nothing crashes.
//! 2. **mis-reads a base** — every base of the read, in the tract and in its flanks, is
//!    replaced by a different one with probability `substitution_rate`.
//!
//! **Two optional features ride on a non-reference haplotype**, both off by default and drawn
//! only when on, so a fixture without them is the one it always was: a deletion in the right
//! flank, starting a set distance from the tract, and a substituted base in the middle of the
//! tract. Together they are the case a flank indel next to a non-reference tract makes, where
//! the delimiter's junction guard can force a read that crossed the tract to count as partial
//! (found 2026-11-01 at SL4.0ch05:51017179). With the deletion, the read's alignment writes the
//! length change at the tract's start and the deletion where it lies.
//!
//! **Nothing here models a partial-repeat slip, a PCR chimera or a mapping error.** A read is
//! always placed where it came from and always spans its tract whole. So a caller that does
//! badly on this fixture is doing badly on the easiest version of the problem, and one that
//! does well has not thereby been shown to do well on reads.
//!
//! # Settings, as `key=value` after the output directory
//!
//! | key | default | what it does |
//! |---|---|---|
//! | `tracts` | 2000 | how many tracts the contig holds |
//! | `samples` | 1 | how many samples are sequenced |
//! | `depth` | 30 | reads a sample a tract |
//! | `slip_share` | 0.10 | share of reads reporting a length other than their allele's |
//! | `shorter_share` | 0.50 | of those, the share reporting a shorter tract |
//! | `one_step_share` | 0.95 | of those, the share that moved by exactly one repeat |
//! | `substitution_rate` | 0.001 | per-base chance of reading the wrong base |
//! | `variant_rate` | 0.50 | share of (sample, tract) cells carrying a non-reference allele |
//! | `allele_span` | 3 | furthest a drawn allele sits from the reference, in repeats |
//! | `flank_deletion_len` | 0 | bases a non-reference haplotype deletes from the right flank; 0 is none |
//! | `flank_deletion_offset_max` | 12 | furthest the deletion starts from the tract, in bases; tracts cycle 0 to this |
//! | `flank_deletion_share` | 1.0 | share of non-reference haplotypes carrying the deletion |
//! | `impure_share` | 0.0 | share of non-reference haplotypes whose tract has its middle base changed |
//! | `seed` | 20260902 | the draw's seed; the whole fixture is a function of it |
//!
//! The defaults are the shipped stutter model's own numbers — `slip_share` 0.10,
//! `shorter_share` 0.50, `one_step_share` 0.95 — so a run at the defaults is the case where
//! the caller's model is exactly right, and every departure from them is a case where it is
//! not.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// The contig every tract sits on. One, so a scorer's BED needs no contig list.
const CONTIG: &str = "chrSIM";

/// How much flank a tract carries either side of it, in bases.
///
/// **Two jobs, and the larger of the two requirements sets it.** A read has to reach past the
/// tract far enough for the aligner to anchor — twenty bases is ng's own default flank width —
/// and two tracts have to sit far enough apart that the region typing calls them two tracts
/// rather than one cluster with no clean flanks. The second is what makes this 150 rather
/// than 30.
const FLANK_BP: usize = 150;

/// The bases a read reaches past its tract, at least, on each side.
const MIN_ANCHOR_BP: usize = 25;

/// How much the per-read anchor varies, so reads do not all start at one position.
const ANCHOR_JITTER_BP: usize = 20;

/// The motif and reference repeat count of each period this fixture lays down.
///
/// **Every one clears the calling floors** — `SsrSegmentCriteria::default`'s
/// `[8, 6, 6, 6, 5, 4]` copies over periods 1 to 6 — with room to spare, so a drawn allele
/// several repeats short of the reference still leaves a tract the run routes to its repeat
/// path. None reaches the 100 bp cap above which a run types a tandem array as a satellite it
/// will not call.
const TRACT_SHAPES: [(&[u8], u32); 6] = [
    (b"A", 14),
    (b"AC", 11),
    (b"AGC", 9),
    (b"AAGC", 8),
    (b"AAGGC", 7),
    (b"AAGGTC", 6),
];

/// The four bases, in the order every choice here walks them.
const FOUR_BASES: [u8; 4] = *b"ACGT";

/// The largest slip this program will draw, in whole repeats.
///
/// A geometric tail is unbounded and a read reporting a tract twenty repeats from its allele
/// is not slippage. The cap is generous against `one_step_share` 0.95, where a slip of five
/// already happens about one time in a million.
const MAX_SLIP: u32 = 6;

// ---------------------------------------------------------------------------
// The draw
// ---------------------------------------------------------------------------

/// SplitMix64 — dependency-free, so the fixture's determinism rests on this file and not on a
/// crate's version.
struct Draw {
    state: u64,
}

impl Draw {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform draw on `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// A uniform draw on `0..n`.
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n.max(1)
    }
}

// ---------------------------------------------------------------------------
// Settings
// ---------------------------------------------------------------------------

/// Everything the fixture is a function of, with the shipped stutter model's own numbers as
/// the defaults — see the module docstring's table.
#[derive(Debug, Clone, Copy)]
struct Settings {
    tracts: usize,
    samples: usize,
    depth: u32,
    slip_share: f64,
    shorter_share: f64,
    one_step_share: f64,
    substitution_rate: f64,
    variant_rate: f64,
    allele_span: u32,
    flank_deletion_len: usize,
    flank_deletion_offset_max: usize,
    flank_deletion_share: f64,
    impure_share: f64,
    seed: u64,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            tracts: 2000,
            samples: 1,
            depth: 30,
            slip_share: 0.10,
            shorter_share: 0.50,
            one_step_share: 0.95,
            substitution_rate: 0.001,
            variant_rate: 0.50,
            allele_span: 3,
            flank_deletion_len: 0,
            flank_deletion_offset_max: 12,
            flank_deletion_share: 1.0,
            impure_share: 0.0,
            seed: 20_260_902,
        }
    }
}

impl Settings {
    /// Read `key=value` pairs over the defaults, refusing anything unrecognised.
    ///
    /// **A misspelled key is an error rather than a silent default**, because the whole point
    /// of the program is that a run is a function of these numbers and a typo would produce a
    /// fixture nobody asked for while reporting the settings they did.
    fn parse(arguments: &[String]) -> Result<Self, String> {
        let mut settings = Self::default();
        for argument in arguments {
            let (key, value) = argument
                .split_once('=')
                .ok_or_else(|| format!("expected key=value, got `{argument}`"))?;
            let number = |what: &str| -> Result<f64, String> {
                value
                    .parse::<f64>()
                    .map_err(|_| format!("{what} must be a number, got `{value}`"))
            };
            match key {
                "tracts" => settings.tracts = number(key)? as usize,
                "samples" => settings.samples = number(key)? as usize,
                "depth" => settings.depth = number(key)? as u32,
                "slip_share" => settings.slip_share = number(key)?,
                "shorter_share" => settings.shorter_share = number(key)?,
                "one_step_share" => settings.one_step_share = number(key)?,
                "substitution_rate" => settings.substitution_rate = number(key)?,
                "variant_rate" => settings.variant_rate = number(key)?,
                "allele_span" => settings.allele_span = number(key)? as u32,
                "flank_deletion_len" => settings.flank_deletion_len = number(key)? as usize,
                "flank_deletion_offset_max" => {
                    settings.flank_deletion_offset_max = number(key)? as usize
                }
                "flank_deletion_share" => settings.flank_deletion_share = number(key)?,
                "impure_share" => settings.impure_share = number(key)?,
                "seed" => settings.seed = number(key)? as u64,
                other => return Err(format!("unknown setting `{other}`")),
            }
        }
        if settings.tracts == 0 || settings.samples == 0 || settings.depth == 0 {
            return Err("tracts, samples and depth must all be at least one".to_string());
        }
        for (what, share) in [
            ("slip_share", settings.slip_share),
            ("shorter_share", settings.shorter_share),
            ("one_step_share", settings.one_step_share),
            ("substitution_rate", settings.substitution_rate),
            ("variant_rate", settings.variant_rate),
            ("flank_deletion_share", settings.flank_deletion_share),
            ("impure_share", settings.impure_share),
        ] {
            if !(0.0..=1.0).contains(&share) {
                return Err(format!("{what} must be a share in [0, 1], got {share}"));
            }
        }
        if settings.one_step_share <= 0.0 {
            return Err("one_step_share must be above zero, or no slip has a size".to_string());
        }
        // Every read must still carry aligned flank past the deletion, or it stops being a
        // read that crossed the tract.
        if settings.flank_deletion_len > 0
            && settings.flank_deletion_offset_max + settings.flank_deletion_len + 5 > MIN_ANCHOR_BP
        {
            return Err(format!(
                "flank_deletion_offset_max + flank_deletion_len must leave five aligned bases \
                 inside the {MIN_ANCHOR_BP}-base shortest anchor"
            ));
        }
        Ok(settings)
    }
}

// ---------------------------------------------------------------------------
// The reference
// ---------------------------------------------------------------------------

/// One tract on the contig, as the fixture laid it down.
#[derive(Debug, Clone)]
struct Tract {
    motif: Vec<u8>,
    reference_repeats: u32,
    /// One-based, the first base of the repeat itself.
    start: usize,
    /// How many reference bases lie between the tract's last base and the first base a
    /// haplotype's flank deletion takes — zero touches the tract. Cycles through
    /// `0..=flank_deletion_offset_max` tract by tract, each period in turn, so every period meets
    /// every distance.
    flank_deletion_offset: usize,
}

impl Tract {
    fn period(&self) -> usize {
        self.motif.len()
    }

    /// How many bases the reference tract spans.
    fn reference_len(&self) -> usize {
        self.period() * self.reference_repeats as usize
    }

    /// One-based, the last base of the repeat.
    ///
    /// Arithmetic rather than the length of [`Self::reference_bases`], because this is asked
    /// once a read and that one allocates the tract to measure it.
    fn end(&self) -> usize {
        self.start + self.reference_len() - 1
    }

    /// The bases this many repeats of the motif spell.
    fn bases_at(&self, repeats: u32) -> Vec<u8> {
        self.motif.repeat(repeats as usize).to_vec()
    }

    /// The bases this many repeats spell, with the substitution inside the tract that the ch05
    /// case carried when `impure`: the middle base of an allele of `allele_repeats` repeats.
    ///
    /// **The substitution stays where the allele has it when a read slips**, since a slip adds or
    /// removes whole repeats at the tract's end. Placing it in the middle of the slipped length
    /// instead turned every slipped read into a different sequence rather than a stutter copy,
    /// and the caller rightly called it a second allele.
    fn allele_bases(&self, repeats: u32, impure: bool, allele_repeats: u32) -> Vec<u8> {
        let mut bases = self.bases_at(repeats);
        let middle = self.period() * allele_repeats as usize / 2;
        if impure && middle < bases.len() {
            bases[middle] = other_base(bases[middle]);
        }
        bases
    }
}

/// The contig, and where each tract sits on it.
fn build_reference(settings: &Settings, draw: &mut Draw) -> (Vec<u8>, Vec<Tract>) {
    let mut bases = Vec::new();
    let mut tracts = Vec::with_capacity(settings.tracts);
    for index in 0..settings.tracts {
        let (motif, repeats) = TRACT_SHAPES[index % TRACT_SHAPES.len()];
        let motif = motif.to_vec();
        append_aperiodic(&mut bases, FLANK_BP - motif.len(), draw);
        append_spacer(&mut bases, &motif);
        let start = bases.len() + 1; // one-based
        bases.extend(motif.repeat(repeats as usize));
        append_spacer(&mut bases, &motif);
        append_aperiodic(&mut bases, FLANK_BP - motif.len(), draw);
        tracts.push(Tract {
            motif,
            reference_repeats: repeats,
            start,
            flank_deletion_offset: (index / TRACT_SHAPES.len())
                % (settings.flank_deletion_offset_max + 1),
        });
    }
    (bases, tracts)
}

/// **`period` bases that cannot be read as another copy of the motif** — the tract's boundary,
/// written rather than hoped for.
///
/// Every base differs from the motif base it stands opposite, so the neighbouring `period`
/// bases are not an approximate copy either: the repeat catalog admits a tract at 80% purity,
/// so a flank whose last five bases are `AAGGT` against a motif of `AAGGC` **is** a copy as far
/// as the scan is concerned, and the tract it finds is one copy longer and one base earlier
/// than the truth file claims.
///
/// **Measured, on the version that only forced the single junction base:** 5 of 4,000 tracts
/// came back typed as ordinary sequence, because the accidental near-copy moved the catalog's
/// idea of the tract off the injected coordinates. The spacer is what removed them.
///
/// Where a candidate base would complete a third copy of some short unit — the same condition
/// [`append_aperiodic`] avoids — the next base that still differs from the motif is taken.
///
/// # Panics
///
/// If no base satisfies both. It is not known to be reachable, and a fixture that quietly
/// carried an unasked-for repeat would be a wrong measurement rather than a failed run: the
/// driver's own check — the typed regions must equal the injected tracts — would then fail
/// somewhere else, with nothing saying why.
fn append_spacer(bases: &mut Vec<u8>, motif: &[u8]) {
    for motif_base in motif {
        let chosen = FOUR_BASES.iter().copied().find(|candidate| {
            *candidate != *motif_base && !completes_a_third_copy(bases, *candidate)
        });
        let chosen = chosen.unwrap_or_else(|| {
            panic!(
                "no base at position {} both differs from the motif's `{}` and leaves the \
                 flank free of a three-copy repeat",
                bases.len(),
                *motif_base as char,
            )
        });
        bases.push(chosen);
    }
}

/// Append `len` bases carrying no tandem repeat of period 1 to 6 longer than two copies,
/// **checked against everything already on the contig** rather than against the flank alone.
///
/// The base is drawn and then the four are tried in turn from it, so the sequence stays a
/// function of the seed while the constraint holds.
///
/// # Panics
///
/// If all four bases would complete a third copy — for the reason [`append_spacer`] gives.
fn append_aperiodic(bases: &mut Vec<u8>, len: usize, draw: &mut Draw) {
    for _ in 0..len {
        let drawn = draw.below(4) as usize;
        let chosen = (0..FOUR_BASES.len())
            .map(|step| FOUR_BASES[(drawn + step) % FOUR_BASES.len()])
            .find(|candidate| !completes_a_third_copy(bases, *candidate))
            .unwrap_or_else(|| {
                panic!(
                    "every base at position {} would complete a three-copy repeat",
                    bases.len()
                )
            });
        bases.push(chosen);
    }
}

/// Whether appending `base` would make three consecutive copies of any 1-to-6-base unit.
fn completes_a_third_copy(bases: &[u8], base: u8) -> bool {
    for period in 1..=6_usize {
        if bases.len() < period * 3 - 1 {
            continue;
        }
        let end = bases.len();
        let unit = &bases[end + 1 - period..];
        // The candidate is the last base of the third copy, so the unit is the last `period`
        // bases with `base` appended.
        let mut third: Vec<u8> = unit.to_vec();
        third.push(base);
        let second = &bases[end + 1 - 2 * period..end + 1 - period];
        let first = &bases[end + 1 - 3 * period..end + 1 - 2 * period];
        if third == second && second == first {
            return true;
        }
    }
    false
}

fn other_base(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

// ---------------------------------------------------------------------------
// The genotypes and the reads
// ---------------------------------------------------------------------------

/// One haplotype at one tract: its repeat count, whether its tract carries a substitution, and
/// whether its right flank carries the deletion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Haplotype {
    repeats: u32,
    impure: bool,
    flank_deleted: bool,
}

impl Haplotype {
    fn reference(tract: &Tract) -> Self {
        Self {
            repeats: tract.reference_repeats,
            impure: false,
            flank_deleted: false,
        }
    }

    /// What this haplotype reads in place of the reference from the tract's first base to the
    /// one-based `span_end`, which must reach past any flank deletion.
    fn spelled(
        &self,
        tract: &Tract,
        reference: &[u8],
        settings: &Settings,
        span_end: usize,
    ) -> Vec<u8> {
        let mut bases = tract.allele_bases(self.repeats, self.impure, self.repeats);
        let flank = &reference[tract.end()..span_end];
        if self.flank_deleted {
            let from = tract.flank_deletion_offset;
            bases.extend_from_slice(&flank[..from]);
            bases.extend_from_slice(&flank[from + settings.flank_deletion_len..]);
        } else {
            bases.extend_from_slice(flank);
        }
        bases
    }
}

/// One sample's two haplotypes at one tract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Genotype {
    first: Haplotype,
    second: Haplotype,
}

impl Genotype {
    fn is_reference(&self, tract: &Tract) -> bool {
        let reference = Haplotype::reference(tract);
        self.first == reference && self.second == reference
    }
}

/// Draw one sample's genotype at one tract.
///
/// **A non-reference allele is never allowed below one repeat**, because a tract of nothing is
/// not a tract and the caller's ladder has no rung for it.
fn draw_genotype(tract: &Tract, settings: &Settings, draw: &mut Draw) -> Genotype {
    let reference_haplotype = Haplotype::reference(tract);
    let reference = tract.reference_repeats;
    if draw.unit() >= settings.variant_rate {
        return Genotype {
            first: reference_haplotype,
            second: reference_haplotype,
        };
    }
    let span = settings.allele_span.max(1);
    let mut alternative = reference;
    while alternative == reference {
        let step = 1 + draw.below(u64::from(span)) as u32;
        alternative = if draw.unit() < 0.5 {
            reference.saturating_sub(step).max(1)
        } else {
            reference + step
        };
        if alternative == reference {
            // A downward step clamped back onto the reference; try again rather than
            // silently emitting a homozygous-reference cell the truth file calls variant.
            alternative = reference + step;
        }
    }
    // Heterozygous half the time, homozygous for the alternative the other half — a shape
    // that gives both the "one allele hides under the other's stutter" case and the "no
    // reference reads at all" case.
    //
    // The two optional features ride on the alternative haplotype. Each is drawn only when
    // switched on, so a fixture with both off draws exactly what it drew before they existed.
    let alternative = Haplotype {
        repeats: alternative,
        impure: settings.impure_share > 0.0 && draw.unit() < settings.impure_share,
        flank_deleted: settings.flank_deletion_len > 0
            && draw.unit() < settings.flank_deletion_share,
    };
    if draw.unit() < 0.5 {
        Genotype {
            first: reference_haplotype,
            second: alternative,
        }
    } else {
        Genotype {
            first: alternative,
            second: alternative,
        }
    }
}

/// How many repeats one read reports, given the allele it came from.
fn observed_repeats(allele: u32, settings: &Settings, draw: &mut Draw) -> u32 {
    if draw.unit() >= settings.slip_share {
        return allele;
    }
    let mut size = 1_u32;
    while size < MAX_SLIP && draw.unit() >= settings.one_step_share {
        size += 1;
    }
    if draw.unit() < settings.shorter_share {
        allele.saturating_sub(size).max(1)
    } else {
        allele + size
    }
}

/// One read: where it sits on the reference, what it says, and how it aligns.
struct Read {
    reference_start: usize, // one-based
    sequence: Vec<u8>,
    cigar: Vec<(char, usize)>,
    name: String,
}

/// Every read one sample carries at one tract.
fn reads_at(
    tract: &Tract,
    genotype: Genotype,
    reference: &[u8],
    settings: &Settings,
    draw: &mut Draw,
    name_prefix: &str,
    reads: &mut Vec<Read>,
) {
    let reference_tract_len = tract.reference_len();
    for index in 0..settings.depth {
        let left_anchor = MIN_ANCHOR_BP + (draw.below(ANCHOR_JITTER_BP as u64) as usize);
        let right_anchor = MIN_ANCHOR_BP + (draw.below(ANCHOR_JITTER_BP as u64) as usize);
        // The fixture lays down `FLANK_BP` either side, so an anchor can never run off a
        // neighbour's tract; the assertion says so rather than trusting the arithmetic.
        assert!(
            left_anchor <= FLANK_BP && right_anchor <= FLANK_BP,
            "a read's anchor reached past its tract's flank"
        );
        let allele = if draw.unit() < 0.5 {
            genotype.first
        } else {
            genotype.second
        };
        let observed = observed_repeats(allele.repeats, settings, draw);
        let observed_bases = tract.allele_bases(observed, allele.impure, allele.repeats);

        let reference_start = tract.start - left_anchor;
        let mut sequence = Vec::with_capacity(left_anchor + observed_bases.len() + right_anchor);
        sequence.extend_from_slice(&reference[reference_start - 1..tract.start - 1]);
        sequence.extend_from_slice(&observed_bases);
        let right_flank = &reference[tract.end()..tract.end() + right_anchor];
        let deletion = allele.flank_deleted.then_some(
            tract.flank_deletion_offset..tract.flank_deletion_offset + settings.flank_deletion_len,
        );
        match &deletion {
            Some(deleted) => {
                sequence.extend_from_slice(&right_flank[..deleted.start]);
                sequence.extend_from_slice(&right_flank[deleted.end..]);
            }
            None => sequence.extend_from_slice(right_flank),
        }
        for base in &mut sequence {
            if draw.unit() < settings.substitution_rate {
                *base = other_base(*base);
            }
        }

        let observed_len = observed_bases.len();
        // **With a flank deletion the length change goes at the tract's start**, as a
        // left-aligning mapper writes it in a run, so it can never abut the flank deletion as an
        // insertion beside a deletion.
        if let Some(deleted) = deletion {
            let mut cigar = vec![('M', left_anchor)];
            if observed_len > reference_tract_len {
                cigar.push(('I', observed_len - reference_tract_len));
            } else if observed_len < reference_tract_len {
                cigar.push(('D', reference_tract_len - observed_len));
            }
            cigar.push(('M', observed_len.min(reference_tract_len) + deleted.start));
            cigar.push(('D', deleted.len()));
            cigar.push(('M', right_anchor - deleted.end));
            reads.push(Read {
                reference_start,
                sequence,
                cigar: merged(cigar),
                name: format!("{name_prefix}_{}_{index}", tract.start),
            });
            continue;
        }
        let cigar = if observed_len == reference_tract_len {
            vec![('M', left_anchor + reference_tract_len + right_anchor)]
        } else if observed_len > reference_tract_len {
            vec![
                ('M', left_anchor + reference_tract_len),
                ('I', observed_len - reference_tract_len),
                ('M', right_anchor),
            ]
        } else {
            vec![
                ('M', left_anchor + observed_len),
                ('D', reference_tract_len - observed_len),
                ('M', right_anchor),
            ]
        };
        reads.push(Read {
            reference_start,
            sequence,
            cigar,
            name: format!("{name_prefix}_{}_{index}", tract.start),
        });
    }
}

/// A CIGAR with zero-length operations dropped and neighbours of one kind joined.
fn merged(cigar: Vec<(char, usize)>) -> Vec<(char, usize)> {
    let mut out: Vec<(char, usize)> = Vec::with_capacity(cigar.len());
    for (kind, len) in cigar {
        if len == 0 {
            continue;
        }
        match out.last_mut() {
            Some((last, total)) if *last == kind => *total += len,
            _ => out.push((kind, len)),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

fn write_fasta(path: &Path, bases: &[u8]) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(file, ">{CONTIG}")?;
    for chunk in bases.chunks(60) {
        file.write_all(chunk)?;
        writeln!(file)?;
    }
    file.flush()
}

/// A coordinate-sorted single-contig BAM for one sample.
fn write_bam(
    path: &Path,
    contig_len: usize,
    sample: &str,
    reads: &mut [Read],
) -> std::io::Result<()> {
    use bstr::BString;
    use noodles_bam as bam;
    use noodles_core::Position;
    use noodles_sam as sam;
    use sam::alignment::io::Write as _;
    use sam::alignment::record::MappingQuality;
    use sam::alignment::record::cigar::Op;
    use sam::alignment::record::cigar::op::Kind;
    use sam::alignment::record_buf::{QualityScores, Sequence};
    use sam::header::record::value::Map;
    use sam::header::record::value::map::header::Version;
    use sam::header::record::value::map::header::tag::SORT_ORDER;
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    use sam::header::record::value::map::{Header as HeaderMap, ReadGroup, ReferenceSequence};
    use std::num::NonZero;

    reads.sort_by_key(|read| read.reference_start);

    let mut hd = Map::<HeaderMap>::new(Version::new(1, 6));
    hd.other_fields_mut()
        .insert(SORT_ORDER, BString::from("coordinate"));
    let sq = Map::<ReferenceSequence>::new(
        NonZero::new(contig_len).expect("the contig has at least one base"),
    );
    let mut rg = Map::<ReadGroup>::default();
    rg.other_fields_mut().insert(SAMPLE, BString::from(sample));
    let read_group = format!("rg_{sample}");
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(CONTIG.as_bytes().to_vec(), sq)
        .add_read_group(read_group.as_bytes().to_vec(), rg)
        .build();

    let mut writer = bam::io::Writer::new(std::fs::File::create(path)?);
    writer.write_header(&header)?;
    for read in reads.iter() {
        let cigar: sam::alignment::record_buf::Cigar = read
            .cigar
            .iter()
            .map(|(kind, len)| {
                let kind = match kind {
                    'I' => Kind::Insertion,
                    'D' => Kind::Deletion,
                    _ => Kind::Match,
                };
                Op::new(kind, *len)
            })
            .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_name(read.name.as_bytes())
            .set_reference_sequence_id(0)
            .set_flags(sam::alignment::record::Flags::empty())
            .set_mapping_quality(MappingQuality::new(60).expect("60 is a mapping quality"))
            .set_alignment_start(
                Position::try_from(read.reference_start).expect("starts are one-based"),
            )
            .set_cigar(cigar)
            .set_sequence(Sequence::from(read.sequence.clone()))
            .set_quality_scores(QualityScores::from(vec![35_u8; read.sequence.len()]))
            .set_data(
                [(
                    sam::alignment::record::data::field::Tag::READ_GROUP,
                    sam::alignment::record_buf::data::field::Value::String(BString::from(
                        read_group.clone(),
                    )),
                )]
                .into_iter()
                .collect(),
            )
            .build();
        writer.write_alignment_record(&header, &record)?;
    }
    writer.try_finish()
}

/// The truth VCF: one record per tract any sample varies at, in the tract's own coordinates.
///
/// **The REF is the whole reference tract and each ALT a whole alternative tract**, which is
/// the shape ng's own tract records take. Where a haplotype carries the flank deletion, every
/// allele of the record runs on to the deletion's last base, so the deletion is part of the
/// haplotype's allele rather than a second record. Left-alignment is left to `bcftools norm`, which the
/// scorer runs over both sides, so the two meet in the same spelling however either is written
/// here.
fn write_truth_vcf(
    path: &Path,
    tracts: &[Tract],
    genotypes: &[Vec<Genotype>],
    samples: &[String],
    reference: &[u8],
    settings: &Settings,
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(file, "##fileformat=VCFv4.2")?;
    writeln!(file, "##contig=<ID={CONTIG}>")?;
    writeln!(
        file,
        "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Repeat unit length in bases\">"
    )?;
    writeln!(
        file,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        file,
        "##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description=\"Repeat copy number of each \
         called allele, GT order\">"
    )?;
    write!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    )?;
    for sample in samples {
        write!(file, "\t{sample}")?;
    }
    writeln!(file)?;

    for (index, tract) in tracts.iter().enumerate() {
        let reference_haplotype = Haplotype::reference(tract);
        let mut alternatives: Vec<Haplotype> = Vec::new();
        for sample in genotypes {
            for allele in [sample[index].first, sample[index].second] {
                if allele != reference_haplotype && !alternatives.contains(&allele) {
                    alternatives.push(allele);
                }
            }
        }
        if alternatives.is_empty() {
            continue;
        }
        alternatives.sort_unstable_by_key(|one| (one.repeats, one.impure, one.flank_deleted));
        let span_end = record_end(tract, &alternatives, settings);
        let spelled = |haplotype: &Haplotype| {
            String::from_utf8(haplotype.spelled(tract, reference, settings, span_end))
                .expect("ASCII bases")
        };
        let alt_field = alternatives
            .iter()
            .map(spelled)
            .collect::<Vec<_>>()
            .join(",");
        let allele_index = |haplotype: Haplotype| -> usize {
            if haplotype == reference_haplotype {
                0
            } else {
                1 + alternatives
                    .iter()
                    .position(|one| *one == haplotype)
                    .expect("every drawn allele is either the reference or an alternative")
            }
        };
        let mut line = String::new();
        write!(
            line,
            "{CONTIG}\t{}\t.\t{}\t{alt_field}\t100\tPASS\tPERIOD={}\tGT:REPCN",
            tract.start,
            spelled(&reference_haplotype),
            tract.period(),
        )
        .expect("writing to a String");
        for sample in genotypes {
            let genotype = sample[index];
            write!(
                line,
                "\t{}/{}:{},{}",
                allele_index(genotype.first),
                allele_index(genotype.second),
                genotype.first.repeats,
                genotype.second.repeats,
            )
            .expect("writing to a String");
        }
        writeln!(file, "{line}")?;
    }
    file.flush()
}

/// One-based, the last base a tract's truth record covers: the tract's own, or past the flank
/// deletion when any haplotype listed carries it.
fn record_end(tract: &Tract, haplotypes: &[Haplotype], settings: &Settings) -> usize {
    if haplotypes.iter().any(|one| one.flank_deleted) {
        tract.end() + tract.flank_deletion_offset + settings.flank_deletion_len
    } else {
        tract.end()
    }
}

/// How far either side of a tract the haplotype table's window reaches, in bases.
const HAPLOTYPE_WINDOW_BP: usize = 40;

/// **Every sample's two haplotypes at every tract, spelled over a fixed window** — the truth a
/// scorer compares called haplotypes with, whatever records the caller chose to write them as.
/// The window is `HAPLOTYPE_WINDOW_BP` either side of the tract, one-based and inclusive.
fn write_haplotype_table(
    path: &Path,
    tracts: &[Tract],
    genotypes: &[Vec<Genotype>],
    samples: &[String],
    reference: &[u8],
    settings: &Settings,
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "chrom\twindow_start\twindow_end\tperiod\tdeletion_offset\tsample\tvariant\timpure\t\
         flank_deleted\tfirst\tsecond"
    )?;
    for (index, tract) in tracts.iter().enumerate() {
        let window_start = tract.start - HAPLOTYPE_WINDOW_BP;
        let window_end = tract.end() + HAPLOTYPE_WINDOW_BP;
        let left = &reference[window_start - 1..tract.start - 1];
        let spell = |haplotype: &Haplotype| {
            let mut bases = left.to_vec();
            bases.extend(haplotype.spelled(tract, reference, settings, window_end));
            String::from_utf8(bases).expect("ASCII bases")
        };
        for (sample_index, sample) in samples.iter().enumerate() {
            let genotype = genotypes[sample_index][index];
            let carriers = [genotype.first, genotype.second];
            writeln!(
                file,
                "{CONTIG}\t{window_start}\t{window_end}\t{}\t{}\t{sample}\t{}\t{}\t{}\t{}\t{}",
                tract.period(),
                tract.flank_deletion_offset,
                !genotype.is_reference(tract),
                carriers.iter().any(|one| one.impure),
                carriers.iter().any(|one| one.flank_deleted),
                spell(&genotype.first),
                spell(&genotype.second),
            )?;
        }
    }
    file.flush()
}

fn write_tract_bed(path: &Path, tracts: &[Tract]) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    for tract in tracts {
        // BED is half-open and zero-based; the tract's own coordinates are one-based
        // inclusive.
        writeln!(
            file,
            "{CONTIG}\t{}\t{}\t{}",
            tract.start - 1,
            tract.end(),
            tract.period()
        )?;
    }
    file.flush()
}

/// **The slippage this program drew the reads under, as parameters-file rows** — what a fit
/// over these reads would have found if fitting were a thing this caller could do yet.
///
/// `doc/devel/ng/spec/calling_loop_ssr.md` §3.3 asks the QUAL experiment to split its cells by
/// where the model's numbers came from, fitted against defaulted. **No command produces a
/// fitted parameters file** (§3.4 — deferred future work), so on a benchmark every cell is
/// defaulted and the split has one side. Here it has both, because the truth is ours: a run
/// handed these rows is scoring the reads under the model they came from, and a run given
/// `--defaults` is scoring them under the shipped HipSTR numbers. The gap between the two is
/// the axis.
///
/// **A fragment, not a whole file.** The caller writes a complete `.parameters.toml` beside
/// every run's output, and this is the one table to paste into it — the file's own commentary
/// invites exactly that ("There is nothing here to change, only rows to add"), and reusing the
/// caller's writer for the other six groups of numbers means the fixture cannot drift from the
/// file format.
///
/// **One row per `(period, repeat count, slippage group)` the reads can reach.** A tract's
/// scoring context is looked up by the *candidate's* repeat count rather than its reference
/// tract's, so a candidate several repeats either side of the reference needs a row of its own
/// or it falls back to the default and the arm stops being what it claims. The span covers
/// every allele this fixture draws plus every slip it can put on one.
///
/// **The substitution rate rows are deliberately left out.** This fixture's per-base error is
/// the caller's own stated default, so a defaulted cell and a fitted one would carry the same
/// number and the row would say nothing.
fn write_slippage_rows(path: &Path, settings: &Settings) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "# The stutter model `ng_tract_simulator` drew these reads under, as the rows a\n\
         # parameters file states slippage in. Paste them over the\n\
         # run's own `.parameters.toml`: delete its `slippage_by_stratum_and_group = []`\n\
         # line and append this file to the end of it. **The end, not in place** — an\n\
         # array-of-tables closes the table it sits in, so pasting these where the empty\n\
         # array was would make every key after them a key of the last row.\n\
         #\n\
         # share_of_reads_that_slip = {slip}, shorter_share = {shorter}, fall_off = {fall_off:.9}\n\
         # (fall_off is one minus the share of slips that moved by exactly one repeat).",
        slip = settings.slip_share,
        shorter = settings.shorter_share,
        fall_off = 1.0 - settings.one_step_share,
    )?;
    let reach = i64::from(settings.allele_span) + i64::from(MAX_SLIP);
    let mut rows = 0_usize;
    for (motif, reference_repeats) in TRACT_SHAPES {
        let period = motif.len();
        let lowest = (i64::from(reference_repeats) - reach).max(1);
        let highest = i64::from(reference_repeats) + reach;
        for repeats in lowest..=highest {
            writeln!(file, "[[repeat_tracts.slippage_by_stratum_and_group]]")?;
            writeln!(file, "period = {period}")?;
            writeln!(file, "reference_repeats = {repeats}")?;
            writeln!(file, "slippage_group = 0")?;
            writeln!(file, "share_of_reads_that_slip = {}", settings.slip_share)?;
            writeln!(file, "shorter_share = {}", settings.shorter_share)?;
            writeln!(file, "fall_off = {:.9}", 1.0 - settings.one_step_share)?;
            writeln!(
                file,
                "share_of_reads_that_slip_origin = {{ smoothing = \"this_stratum\", \
                 expected_slipped_reads = 0.0 }}"
            )?;
            writeln!(
                file,
                "shorter_share_and_fall_off_origin = {{ expected_slipped_reads = 0.0, \
                 shorter_share_smoothing = \"this_stratum\", \
                 fall_off_smoothing = \"this_stratum\" }}"
            )?;
            writeln!(file)?;
            rows += 1;
        }
    }
    file.flush()?;
    println!("slippage rows     : {rows}");
    Ok(())
}

fn write_confident_bed(path: &Path, contig_len: usize) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "{CONTIG}\t0\t{contig_len}")
}

fn write_genotype_table(
    path: &Path,
    tracts: &[Tract],
    genotypes: &[Vec<Genotype>],
    samples: &[String],
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "chrom\tstart\tend\tperiod\tmotif\treference_repeats\tsample\tfirst\tsecond\tis_reference"
    )?;
    for (index, tract) in tracts.iter().enumerate() {
        for (sample_index, sample) in samples.iter().enumerate() {
            let genotype = genotypes[sample_index][index];
            writeln!(
                file,
                "{CONTIG}\t{}\t{}\t{}\t{}\t{}\t{sample}\t{}\t{}\t{}",
                tract.start - 1,
                tract.end(),
                tract.period(),
                String::from_utf8(tract.motif.clone()).expect("ASCII bases"),
                tract.reference_repeats,
                genotype.first.repeats,
                genotype.second.repeats,
                genotype.is_reference(tract),
            )?;
        }
    }
    file.flush()
}

// ---------------------------------------------------------------------------

fn run(out_dir: &Path, settings: Settings) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(out_dir)?;
    let mut draw = Draw::new(settings.seed);
    let (reference, tracts) = build_reference(&settings, &mut draw);

    let samples: Vec<String> = (0..settings.samples)
        .map(|index| format!("sim{index:03}"))
        .collect();
    let mut genotypes: Vec<Vec<Genotype>> = Vec::with_capacity(samples.len());
    for _ in &samples {
        genotypes.push(
            tracts
                .iter()
                .map(|tract| draw_genotype(tract, &settings, &mut draw))
                .collect(),
        );
    }

    write_fasta(&out_dir.join("reference.fa"), &reference)?;
    write_tract_bed(&out_dir.join("tracts.bed"), &tracts)?;
    write_confident_bed(&out_dir.join("confident.bed"), reference.len())?;
    write_truth_vcf(
        &out_dir.join("truth.vcf"),
        &tracts,
        &genotypes,
        &samples,
        &reference,
        &settings,
    )?;
    write_haplotype_table(
        &out_dir.join("truth_haplotypes.tsv"),
        &tracts,
        &genotypes,
        &samples,
        &reference,
        &settings,
    )?;
    write_slippage_rows(&out_dir.join("slippage_rows.toml"), &settings)?;
    write_genotype_table(
        &out_dir.join("truth_genotypes.tsv"),
        &tracts,
        &genotypes,
        &samples,
    )?;

    let mut variant_cells = 0_usize;
    for (sample_index, sample) in samples.iter().enumerate() {
        let mut reads = Vec::with_capacity(tracts.len() * settings.depth as usize);
        for (tract_index, tract) in tracts.iter().enumerate() {
            let genotype = genotypes[sample_index][tract_index];
            if !genotype.is_reference(tract) {
                variant_cells += 1;
            }
            reads_at(
                tract, genotype, &reference, &settings, &mut draw, sample, &mut reads,
            );
        }
        write_bam(
            &out_dir.join(format!("{sample}.bam")),
            reference.len(),
            sample,
            &mut reads,
        )?;
    }

    let varying_tracts = (0..tracts.len())
        .filter(|index| {
            genotypes
                .iter()
                .any(|sample| !sample[*index].is_reference(&tracts[*index]))
        })
        .count();
    println!("out dir           : {}", out_dir.display());
    println!("contig            : {CONTIG}, {} bases", reference.len());
    println!("tracts            : {}", tracts.len());
    println!("samples           : {}", samples.len());
    println!(
        "reads a sample    : {}",
        tracts.len() * settings.depth as usize
    );
    println!(
        "tracts any sample varies at : {varying_tracts} of {} ({:.1}%)",
        tracts.len(),
        100.0 * varying_tracts as f64 / tracts.len() as f64
    );
    println!(
        "(sample, tract) cells carrying a non-reference allele : {variant_cells} of {} ({:.1}%)",
        tracts.len() * samples.len(),
        100.0 * variant_cells as f64 / (tracts.len() * samples.len()) as f64
    );
    println!("settings          : {settings:?}");
    println!();
    println!("next, from the output directory:");
    println!("  samtools faidx reference.fa");
    println!("  for bam in sim*.bam; do samtools index \"$bam\"; done");
    println!("  pvc repeat-catalog --reference reference.fa");
    println!(
        "  pvc call-from-alignments --reference reference.fa \\\n    \
         --alignment sim000.bam --regions confident.bed --output calls.vcf --defaults"
    );
    Ok(())
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let Some(out_dir) = arguments.first().map(PathBuf::from) else {
        eprintln!("usage: ng_tract_simulator <out-dir> [key=value ...]");
        eprintln!("settings: see this example's module documentation");
        return ExitCode::from(2);
    };
    let settings = match Settings::parse(&arguments[1..]) {
        Ok(settings) => settings,
        Err(message) => {
            eprintln!("{message}");
            return ExitCode::from(2);
        }
    };
    match run(&out_dir, settings) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}
