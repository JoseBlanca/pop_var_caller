//! **A cohort whose STR genotypes we chose, read through libraries that stutter as real ones
//! do** — the simulated benchmark `doc/devel/ng/spec/str_stutter_benchmark.md` specifies, for
//! scoring ng, GATK HaplotypeCaller and HipSTR against the truth the reads came from.
//!
//! ```text
//! ng_str_stutter_benchmark <out-dir> [key=value ...]
//! ```
//!
//! # Why it exists
//!
//! Neither real benchmark can test the STR model's side of ng: GIAB's human libraries are
//! PCR-free and stutter little, and the tomato libraries stutter but have no truth for their STR
//! lengths. Here the genotypes and the stutter are both ours, so every call can be scored
//! (spec §1).
//!
//! # What it writes
//!
//! Into `<out-dir>`:
//!
//! - `reference.fa` — `contigs` contigs of ordinary sequence with STR segments laid into it,
//!   each segment wrapped in 150 bases of flank carrying no repeat of its own.
//! - `<sample>_R1.fastq.gz`, `<sample>_R2.fastq.gz` — paired-end 2 × 150 reads over the whole
//!   genome, BGZF-compressed (any gzip reader takes them). Mapping them is the driver's job.
//! - `truth.vcf` — every variant any sample carries, genotypes phased (the first copy before
//!   the bar). An STR segment is **one record over its window** (the segment and 30 bases either
//!   side, trimmed to the stretch its alleles differ over), so a length change, a substitution
//!   inside the repeat, a SNP beside it and a deletion across its edge on one haplotype are one
//!   allele. Variants away from the segments are one record each.
//! - `segments.bed` — `chrom start end period`, zero-based half-open, one row per STR segment:
//!   the ground the typing must call STR segments, and nothing else.
//! - `segment_alleles.tsv` — every allele drawn at every segment, carried or not, as its whole
//!   window sequence, with its repeat count, how it differs from the reference and its
//!   population frequency.
//! - `segment_genotypes.tsv` — every sample's two alleles at every segment, as indices into
//!   `segment_alleles.tsv`, including the homozygous-reference ones. With the allele table this
//!   is each sample's two truth sequences at every segment.
//! - `libraries.tsv` — each sample's library: its stutter level and its own multiplier.
//! - `stutter_rates.tsv` — the share of molecules that stutter, by level, motif length and
//!   repeat count, before a library's multiplier (§3.2's table, laid out per repeat count).
//! - `confident.bed` — every contig whole.
//!
//! # The forward model, in the order it is applied
//!
//! 1. **The reference**: `segments_per_period` STR segments of each motif length from one to
//!    six, their reference repeat counts spread evenly from the calling floor (8, 6, 6, 6, 5, 4)
//!    to 30 copies for motif lengths one and two and 15 for the rest, shuffled along the genome
//!    and separated by ordinary sequence.
//! 2. **The population**: at each segment two to six alleles differing by whole units, the
//!    non-reference ones sharing a frequency drawn from a Beta(0.5, 1.5); one segment in ten
//!    gains an allele with a substitution inside the repeat, one in ten an allele with a SNP
//!    within 30 bases of an edge, one in thirty an allele with a deletion crossing an edge.
//!    Away from the segments, SNPs at one per thousand bases and indels of one to ten bases at
//!    one per ten thousand, each with its own Beta frequency.
//! 3. **The genotypes**: each sample's two copies drawn from those frequencies, both copies one
//!    draw with probability `inbreeding`.
//! 4. **The libraries**: each sample's stutter level comes from the `design`, and its multiplier
//!    from a log-normal with a standard deviation of 0.4 on the log scale (median one).
//! 5. **The molecules**: fragments of length drawn around 350 (sd 50), placed uniformly over
//!    the sample's two copies of the genome. At every STR segment a molecule overlaps it gains
//!    or loses whole units: it stutters with probability `rate(level, motif length, repeats) ×
//!    multiplier`, shorter with the level's shorter share, by one unit with probability 0.85,
//!    two with 0.85 × 0.15 and so on, at most four. **Stutter acts on the molecule, so both mates
//!    of a pair report the same length**, as both read one PCR product.
//! 6. **The reads**: the first 150 bases of the molecule and the reverse complement of its last
//!    150, from either strand at random; base qualities drawn around Q32 with a tail falling to
//!    about Q22 at the read's end, and each base misread with the probability its quality
//!    states.
//!
//! **The rate**: on long segments it is the level's plateau (spec §3.2); at the calling floor
//! one thirtieth of it; in between it rises evenly on a log scale, reaching the plateau at 16
//! repeats for motif lengths one and two and at 10 for the rest. **Below the calling floor** —
//! an allele shorter than any reference segment the typing admits — it keeps falling at the same
//! slope; the spec gives no number there and the rate is already a thirtieth of the plateau.
//!
//! **Not modelled** (spec §3.3, §7): stutter that depends on the flanks or on an interruption
//! inside the repeat, partial-unit slips, PCR chimeras, stutter that differs between a sample's
//! two alleles, insertion or deletion sequencing errors, GC bias, duplicates, and anything in a
//! real genome's sequence context.
//!
//! # Settings, as `key=value` after the output directory
//!
//! | key | default | what it does |
//! |---|---|---|
//! | `design` | `pcr` | `clean`, `pcr_free`, `pcr`, `heavy` or `one_outlier` (spec §4.1) |
//! | `samples` | 1 | how many samples are sequenced |
//! | `depth` | 30 | reads a position, both mates counted, averaged over the genome |
//! | `segments_per_period` | 500 | STR segments of each motif length |
//! | `genome_bp` | 4000000 | the reference's length, all contigs together |
//! | `contigs` | 4 | how many contigs the reference is cut into |
//! | `inbreeding` | 0 | the chance a sample's two copies at a site are one draw |
//! | `threads` | all | how many samples are sequenced at once |
//! | `seed` | 20261009 | the draw's seed; the whole output is a function of it and the rest |
//!
//! `one_outlier` makes sample `sim000` heavy PCR and every other library PCR-free.
//!
//! The draw and the flank-laying functions are copied from `ng_tract_simulator.rs`, so that the
//! two fixtures can change separately.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};

use pop_var_caller::float;

/// Each mate's length.
const READ_LEN: usize = 150;

/// The fragment length's mean and standard deviation, in bases.
const FRAGMENT_MEAN: f64 = 350.0;
const FRAGMENT_SD: f64 = 50.0;

/// How much flank a segment carries either side of it, spacer included, in bases.
///
/// Enough for two segments never to read as one cluster of repeats, which is what made it
/// 150 in `ng_tract_simulator`.
const FLANK_BP: usize = 150;

/// How far a segment's window reaches into each flank: the stretch an STR allele is written
/// over, and the furthest from a segment's edge that a SNP or a deletion beside it is drawn.
const WINDOW_FLANK_BP: usize = 30;

/// The fewest copies the calling typing admits as an STR segment, motif lengths one to six —
/// `SsrSegmentCriteria::default`. The lowest reference repeat count laid down, and where the
/// stutter rate is a thirtieth of its plateau.
const CALLING_FLOOR: [u32; 6] = [8, 6, 6, 6, 5, 4];

/// The highest reference repeat count laid down, motif lengths one to six (spec §4.2). None
/// reaches the 100-base cap above which the typing calls a tandem array a satellite.
const HIGHEST_REFERENCE_REPEATS: [u32; 6] = [30, 30, 15, 15, 15, 15];

/// The repeat count at which stutter stops rising, motif lengths one to six (spec §3.2).
const PLATEAU_REPEATS: [u32; 6] = [16, 16, 10, 10, 10, 10];

/// How many times less a segment at the calling floor stutters than one at the plateau.
const FLOOR_TO_PLATEAU_FOLD: f64 = 30.0;

/// The share of stutter events that move by exactly one unit; each further unit is this much
/// rarer again's complement (spec §3.1).
const ONE_UNIT_SHARE: f64 = 0.85;

/// The most units one stutter event moves; the geometric tail beyond it is put on it.
const MAX_STUTTER_UNITS: u32 = 4;

/// The highest share of molecules a library may stutter at, after its multiplier: a heavy-PCR
/// homopolymer at 0.30 times a multiplier of 3.5 would otherwise exceed one.
const HIGHEST_STUTTER_SHARE: f64 = 0.95;

/// The standard deviation of a library's multiplier on the log scale (spec §3.3).
const LIBRARY_MULTIPLIER_LOG_SD: f64 = 0.4;

/// Variants away from the segments, per base (spec §4.2).
const SNP_RATE: f64 = 1e-3;
const INDEL_RATE: f64 = 1e-4;

/// The longest indel drawn away from the segments, in bases.
const LONGEST_GENERIC_INDEL: usize = 10;

/// The fewest bases between two variants away from the segments, and between one and the edge
/// of the ordinary sequence it sits in.
const GENERIC_VARIANT_SPACING: usize = 10;

/// The Beta distribution a site's non-reference frequency is drawn from: mean 0.25, with most
/// sites rare and some common.
const NONREF_FREQUENCY_BETA: (f64, f64) = (0.5, 1.5);

/// An STR allele's distance from the reference, in units: one with probability 0.6, then each
/// further unit with probability 0.4 of the last, at most six.
const LENGTH_STEP_CONTINUE: f64 = 0.4;
const LONGEST_LENGTH_STEP: u32 = 6;

/// The share of segments gaining each kind of allele beyond a length change (spec §4.2).
const INTERIOR_SUBSTITUTION_SHARE: f64 = 0.1;
const FLANK_SNP_SHARE: f64 = 0.1;
const EDGE_DELETION_SHARE: f64 = 1.0 / 30.0;

/// The most consecutive copies of a unit the ordinary sequence between segments may carry, by
/// unit length one to six.
///
/// **Low enough that the typing never finds a segment there**: it admits a repeat at 80% purity,
/// so two short runs joined by one mismatch can pass, and a run of four bases is as long as a
/// homopolymer here goes. The driver checks it — the typed STR segments must be the ones in
/// `segments.bed`.
const GROUND_MAX_COPIES: [usize; 6] = [4, 2, 2, 2, 2, 2];

/// The extra bases a molecule is cut with before its stutter is applied, so a removal still
/// leaves it its full length.
const MOLECULE_SLACK: usize = 64;

/// Base qualities: the mean, its spread, and how far the mean falls by the read's last base.
const QUALITY_MEAN: f64 = 32.0;
const QUALITY_SD: f64 = 2.5;
const QUALITY_TAIL_DROP: f64 = 10.0;
const LOWEST_QUALITY: u8 = 2;
const HIGHEST_QUALITY: u8 = 41;

/// The four bases, in the order every choice here walks them.
const FOUR_BASES: [u8; 4] = *b"ACGT";

/// The motifs drawn at each motif length. A motif listed twice is drawn twice as often: A and T
/// homopolymers, and AC and AT dinucleotides, are the common ones in plant and human genomes.
/// Each is primitive — not itself a repeat of a shorter unit.
const MOTIFS: [&[&str]; 6] = [
    &["A", "T", "A", "T", "C", "G"],
    &[
        "AC", "CA", "GT", "TG", "AT", "TA", "AT", "AG", "GA", "CT", "TC",
    ],
    &[
        "AAT", "ATT", "AAG", "CTT", "ATC", "GAT", "AAC", "GTT", "AGC", "GCT", "CAG", "CTG",
    ],
    &[
        "AAAT", "ATTT", "AAAG", "CTTT", "AAGG", "CCTT", "ATCT", "AGAT", "AAAC", "GTTT",
    ],
    &[
        "AAAAT", "ATTTT", "AAAAC", "GTTTT", "AAAAG", "CTTTT", "AAGAG",
    ],
    &["AAAAAT", "ATTTTT", "AAAAAC", "AAAAAG", "AAATTT", "AACAAT"],
];

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

    /// A standard normal draw (Box–Muller).
    fn normal(&mut self) -> f64 {
        let first = 1.0 - self.unit(); // (0, 1], so the logarithm is finite
        let second = self.unit();
        (-2.0 * float::ln(first)).sqrt() * float::cos(2.0 * std::f64::consts::PI * second)
    }

    /// An exponential draw with mean one.
    fn exponential(&mut self) -> f64 {
        -float::ln(1.0 - self.unit())
    }

    /// A Gamma draw with scale one (Marsaglia and Tsang, 2000; shapes under one by the
    /// `U^(1/shape)` boost).
    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let boost = float::powf(1.0 - self.unit(), 1.0 / shape);
            return self.gamma(shape + 1.0) * boost;
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue;
            }
            let v = v * v * v;
            let u = 1.0 - self.unit();
            if float::ln(u) < 0.5 * x * x + d - d * v + d * float::ln(v) {
                return d * v;
            }
        }
    }

    fn beta(&mut self, alpha: f64, beta: f64) -> f64 {
        let x = self.gamma(alpha);
        let y = self.gamma(beta);
        x / (x + y)
    }

    /// An index drawn with probability proportional to its weight.
    fn categorical(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut target = self.unit() * total;
        for (index, weight) in weights.iter().enumerate() {
            if target < *weight {
                return index;
            }
            target -= weight;
        }
        weights.len() - 1
    }

    fn other_base(&mut self, base: u8) -> u8 {
        let index = FOUR_BASES.iter().position(|one| *one == base).unwrap_or(0);
        FOUR_BASES[(index + 1 + self.below(3) as usize) % FOUR_BASES.len()]
    }
}

/// An independent seed for one sample's reads, so samples can be sequenced in any order and on
/// any number of threads and still come out the same.
fn sample_seed(master: u64, sample: usize) -> u64 {
    let mut draw = Draw::new(master ^ (sample as u64 + 1).wrapping_mul(0xD1B5_4A32_D192_ED03));
    draw.next_u64()
}

// ---------------------------------------------------------------------------
// Settings
// ---------------------------------------------------------------------------

/// A library's stutter level (spec §3.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StutterLevel {
    None,
    PcrFree,
    Pcr,
    HeavyPcr,
}

impl StutterLevel {
    /// The share of molecules that stutter on a long segment, by motif length.
    fn plateau(self, period: usize) -> f64 {
        let table = match self {
            Self::None => [0.0; 6],
            Self::PcrFree => [0.065, 0.015, 0.005, 0.003, 0.003, 0.002],
            Self::Pcr => [0.19, 0.15, 0.05, 0.015, 0.015, 0.015],
            Self::HeavyPcr => [0.30, 0.25, 0.08, 0.03, 0.03, 0.03],
        };
        table[period - 1]
    }

    /// Of the molecules that stutter, the share that come out shorter.
    fn shorter_share(self) -> f64 {
        match self {
            Self::None => 0.5,
            Self::PcrFree => 0.70,
            Self::Pcr => 0.78,
            Self::HeavyPcr => 0.80,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::PcrFree => "pcr_free",
            Self::Pcr => "pcr",
            Self::HeavyPcr => "heavy_pcr",
        }
    }
}

/// The share of molecules a library at `level` stutters at a segment of this motif length and
/// repeat count, before the library's own multiplier — the module docstring's rate.
fn stutter_rate(level: StutterLevel, period: usize, repeats: u32) -> f64 {
    let plateau = level.plateau(period);
    let plateau_repeats = f64::from(PLATEAU_REPEATS[period - 1]);
    let floor = f64::from(CALLING_FLOOR[period - 1]);
    let repeats = f64::from(repeats);
    if repeats >= plateau_repeats {
        return plateau;
    }
    let below_plateau = (plateau_repeats - repeats) / (plateau_repeats - floor);
    plateau * float::powf(FLOOR_TO_PLATEAU_FOLD, -below_plateau)
}

/// A cohort design: which level each sample's library is at (spec §4.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Design {
    Clean,
    PcrFree,
    Pcr,
    Heavy,
    OneOutlier,
}

impl Design {
    fn parse(name: &str) -> Result<Self, String> {
        match name {
            "clean" => Ok(Self::Clean),
            "pcr_free" => Ok(Self::PcrFree),
            "pcr" => Ok(Self::Pcr),
            "heavy" => Ok(Self::Heavy),
            "one_outlier" => Ok(Self::OneOutlier),
            other => Err(format!(
                "design must be clean, pcr_free, pcr, heavy or one_outlier, got `{other}`"
            )),
        }
    }

    fn level_of(self, sample: usize) -> StutterLevel {
        match self {
            Self::Clean => StutterLevel::None,
            Self::PcrFree => StutterLevel::PcrFree,
            Self::Pcr => StutterLevel::Pcr,
            Self::Heavy => StutterLevel::HeavyPcr,
            Self::OneOutlier if sample == 0 => StutterLevel::HeavyPcr,
            Self::OneOutlier => StutterLevel::PcrFree,
        }
    }
}

/// Everything the output is a function of — see the module docstring's table.
#[derive(Debug, Clone, Copy)]
struct Settings {
    design: Design,
    samples: usize,
    depth: f64,
    segments_per_period: usize,
    genome_bp: usize,
    contigs: usize,
    inbreeding: f64,
    threads: usize,
    seed: u64,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            design: Design::Pcr,
            samples: 1,
            depth: 30.0,
            segments_per_period: 500,
            genome_bp: 4_000_000,
            contigs: 4,
            inbreeding: 0.0,
            threads: std::thread::available_parallelism().map_or(1, usize::from),
            seed: 20_261_009,
        }
    }
}

impl Settings {
    /// Read `key=value` pairs over the defaults, refusing anything unrecognised: a typo would
    /// otherwise produce a cohort nobody asked for while reporting the settings they did.
    fn parse(arguments: &[String]) -> Result<Self, String> {
        let mut settings = Self::default();
        for argument in arguments {
            let (key, value) = argument
                .split_once('=')
                .ok_or_else(|| format!("expected key=value, got `{argument}`"))?;
            let number = || -> Result<f64, String> {
                value
                    .parse::<f64>()
                    .map_err(|_| format!("{key} must be a number, got `{value}`"))
            };
            match key {
                "design" => settings.design = Design::parse(value)?,
                "samples" => settings.samples = number()? as usize,
                "depth" => settings.depth = number()?,
                "segments_per_period" => settings.segments_per_period = number()? as usize,
                "genome_bp" => settings.genome_bp = number()? as usize,
                "contigs" => settings.contigs = number()? as usize,
                "inbreeding" => settings.inbreeding = number()?,
                "threads" => settings.threads = number()? as usize,
                "seed" => settings.seed = number()? as u64,
                other => return Err(format!("unknown setting `{other}`")),
            }
        }
        if settings.samples == 0 || settings.contigs == 0 || settings.threads == 0 {
            return Err("samples, contigs and threads must all be at least one".to_string());
        }
        if settings.depth.is_nan() || settings.depth <= 0.0 {
            return Err(format!("depth must be above zero, got {}", settings.depth));
        }
        if !(0.0..=1.0).contains(&settings.inbreeding) {
            return Err(format!(
                "inbreeding must be in [0, 1], got {}",
                settings.inbreeding
            ));
        }
        if settings.design == Design::OneOutlier && settings.samples < 2 {
            return Err(
                "one_outlier needs at least two samples: one outlier, the rest clean of it"
                    .to_string(),
            );
        }
        Ok(settings)
    }
}

// ---------------------------------------------------------------------------
// The reference
// ---------------------------------------------------------------------------

struct Contig {
    name: String,
    bases: Vec<u8>,
    /// Zero-based half-open stretches of ordinary sequence between segments — where the
    /// variants away from the segments are drawn.
    ordinary: Vec<(usize, usize)>,
}

/// One STR segment as laid down, and the alleles the population carries at it.
struct Segment {
    contig: usize,
    motif: Vec<u8>,
    reference_repeats: u32,
    /// Zero-based, the first base of the repeat.
    start: usize,
    /// Allele 0 is the reference.
    alleles: Vec<WindowAllele>,
}

impl Segment {
    fn period(&self) -> usize {
        self.motif.len()
    }

    /// Zero-based, one past the repeat's last base.
    fn end(&self) -> usize {
        self.start + self.period() * self.reference_repeats as usize
    }

    fn window_start(&self) -> usize {
        self.start - WINDOW_FLANK_BP
    }

    fn window_end(&self) -> usize {
        self.end() + WINDOW_FLANK_BP
    }
}

fn segment_name(index: usize) -> String {
    format!("seg{index:05}")
}

/// The contigs, and every segment's place on them; the segments carry no alleles yet.
fn build_reference(
    settings: &Settings,
    draw: &mut Draw,
) -> Result<(Vec<Contig>, Vec<Segment>), String> {
    let mut shapes: Vec<(Vec<u8>, u32)> = Vec::new();
    for period in 1..=6_usize {
        let floor = CALLING_FLOOR[period - 1];
        let span = HIGHEST_REFERENCE_REPEATS[period - 1] - floor + 1;
        let motifs = MOTIFS[period - 1];
        for index in 0..settings.segments_per_period {
            let repeats = floor + (index as u32 % span);
            let motif = motifs[draw.below(motifs.len() as u64) as usize];
            shapes.push((motif.as_bytes().to_vec(), repeats));
        }
    }
    // Fisher–Yates, so motif lengths and repeat counts are interleaved along the genome.
    for index in (1..shapes.len()).rev() {
        let other = draw.below(index as u64 + 1) as usize;
        shapes.swap(index, other);
    }

    let segment_bp: usize = shapes
        .iter()
        .map(|(motif, repeats)| 2 * FLANK_BP + motif.len() * *repeats as usize)
        .sum();
    let gaps = shapes.len() + settings.contigs;
    let mean_gap = settings.genome_bp.saturating_sub(segment_bp) / gaps;
    if mean_gap < 2 * READ_LEN {
        return Err(format!(
            "{} segments leave {mean_gap} bases of ordinary sequence between them in {} bases; \
             raise genome_bp or lower segments_per_period",
            shapes.len(),
            settings.genome_bp
        ));
    }

    let mut contigs = Vec::with_capacity(settings.contigs);
    let mut segments = Vec::with_capacity(shapes.len());
    for contig_index in 0..settings.contigs {
        let first = contig_index * shapes.len() / settings.contigs;
        let last = (contig_index + 1) * shapes.len() / settings.contigs;
        let mut bases = Vec::with_capacity(settings.genome_bp / settings.contigs + mean_gap);
        let mut ordinary = Vec::new();
        let mut gap = |bases: &mut Vec<u8>, draw: &mut Draw| {
            let len = (mean_gap as f64 * (0.6 + 0.8 * draw.unit())) as usize;
            let start = bases.len();
            append_ordinary(bases, len, draw);
            ordinary.push((start, bases.len()));
        };
        for (motif, repeats) in &shapes[first..last] {
            gap(&mut bases, draw);
            append_aperiodic(&mut bases, FLANK_BP - motif.len(), draw);
            append_spacer(&mut bases, motif);
            let start = bases.len();
            bases.extend(motif.repeat(*repeats as usize));
            append_spacer(&mut bases, motif);
            append_aperiodic(&mut bases, FLANK_BP - motif.len(), draw);
            segments.push(Segment {
                contig: contig_index,
                motif: motif.clone(),
                reference_repeats: *repeats,
                start,
                alleles: Vec::new(),
            });
        }
        gap(&mut bases, draw);
        contigs.push(Contig {
            name: format!("chr{}", contig_index + 1),
            bases,
            ordinary,
        });
    }
    Ok((contigs, segments))
}

/// **`period` bases that cannot be read as another copy of the motif** — the segment's
/// boundary, written rather than hoped for. Every base differs from the motif base it stands
/// opposite, because the catalog admits a repeat at 80% purity and a near-copy beside the
/// segment would move its edge (`ng_tract_simulator.rs` measured 5 in 4,000 without it).
///
/// # Panics
///
/// If no base both differs from the motif's and leaves the flank free of a three-copy repeat.
fn append_spacer(bases: &mut Vec<u8>, motif: &[u8]) {
    for motif_base in motif {
        let chosen = FOUR_BASES.iter().copied().find(|candidate| {
            *candidate != *motif_base && !completes_copies(bases, *candidate, &[2; 6])
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

/// Append `len` bases carrying no tandem repeat of a unit of one to six bases longer than two
/// copies, checked against everything already on the contig.
fn append_aperiodic(bases: &mut Vec<u8>, len: usize, draw: &mut Draw) {
    append_capped(bases, len, &[2; 6], draw);
}

/// Append `len` bases of the ordinary sequence between segments: homopolymers up to
/// [`GROUND_MAX_COPIES`] long, no unit of two to six bases more than twice in a row.
fn append_ordinary(bases: &mut Vec<u8>, len: usize, draw: &mut Draw) {
    append_capped(bases, len, &GROUND_MAX_COPIES, draw);
}

/// Append `len` drawn bases, none completing more consecutive copies of a unit than
/// `max_copies` allows for its length. The base is drawn and the four tried in turn from it, so
/// the sequence stays a function of the seed while the constraint holds.
///
/// # Panics
///
/// If all four bases would break the cap.
fn append_capped(bases: &mut Vec<u8>, len: usize, max_copies: &[usize; 6], draw: &mut Draw) {
    for _ in 0..len {
        let drawn = draw.below(4) as usize;
        let chosen = (0..FOUR_BASES.len())
            .map(|step| FOUR_BASES[(drawn + step) % FOUR_BASES.len()])
            .find(|candidate| !completes_copies(bases, *candidate, max_copies))
            .unwrap_or_else(|| {
                panic!(
                    "every base at position {} would break the repeat cap",
                    bases.len()
                )
            });
        bases.push(chosen);
    }
}

/// Whether appending `base` would make more consecutive copies of some unit of one to six bases
/// than `max_copies` allows for that unit's length.
fn completes_copies(bases: &[u8], base: u8, max_copies: &[usize; 6]) -> bool {
    (1..=6_usize).any(|period| {
        let copies = max_copies[period - 1] + 1;
        let span = copies * period;
        if bases.len() + 1 < span {
            return false;
        }
        let tail = &bases[bases.len() + 1 - span..];
        let at = |index: usize| if index == span - 1 { base } else { tail[index] };
        (0..span - period).all(|index| at(index) == at(index + period))
    })
}

// ---------------------------------------------------------------------------
// The population
// ---------------------------------------------------------------------------

/// How an STR allele differs from the reference window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum AlleleChange {
    Reference,
    Length,
    InteriorSubstitution,
    FlankSnp,
    EdgeDeletion,
}

impl AlleleChange {
    fn name(self) -> &'static str {
        match self {
            Self::Reference => "reference",
            Self::Length => "length",
            Self::InteriorSubstitution => "substitution_inside",
            Self::FlankSnp => "snp_beside",
            Self::EdgeDeletion => "deletion_across_edge",
        }
    }

    /// Whether the allele is the reference window or a whole-unit length change of it, with no
    /// other difference — the alleles the other kinds are derived from.
    fn is_pure_length(self) -> bool {
        matches!(self, Self::Reference | Self::Length)
    }
}

/// One allele at one segment, written over the segment's whole window.
#[derive(Debug, Clone)]
struct WindowAllele {
    sequence: Vec<u8>,
    /// Whole units of the motif the allele's repeat holds; the stutter rate is looked up by it.
    repeats: u32,
    change: AlleleChange,
    frequency: f64,
    /// Zero-based half-open, within `sequence`: the longest stretch of pure repeat, where
    /// stutter adds or removes units. `None` when it is shorter than two units.
    pure_span: Option<(usize, usize)>,
}

/// Every allele the population carries at one segment, with its frequency.
fn draw_segment_alleles(segment: &Segment, window: &[u8], draw: &mut Draw) -> Vec<WindowAllele> {
    let period = segment.period();
    let left = &window[..WINDOW_FLANK_BP];
    let right = &window[window.len() - WINDOW_FLANK_BP..];
    let spell = |repeats: u32| {
        let mut sequence = left.to_vec();
        sequence.extend(segment.motif.repeat(repeats as usize));
        sequence.extend_from_slice(right);
        sequence
    };

    let wanted = 2 + draw.below(5) as usize; // two to six alleles
    let mut lengths = vec![segment.reference_repeats];
    for _ in 0..50 {
        if lengths.len() == wanted {
            break;
        }
        let mut step = 1;
        while step < LONGEST_LENGTH_STEP && draw.unit() < LENGTH_STEP_CONTINUE {
            step += 1;
        }
        let candidate = if draw.unit() < 0.5 {
            segment
                .reference_repeats
                .checked_sub(step)
                .filter(|repeats| *repeats >= 1)
        } else {
            Some(segment.reference_repeats + step)
        };
        if let Some(candidate) = candidate.filter(|one| !lengths.contains(one)) {
            lengths.push(candidate);
        }
    }

    let nonreference = draw.beta(NONREF_FREQUENCY_BETA.0, NONREF_FREQUENCY_BETA.1);
    let weights: Vec<f64> = (1..lengths.len()).map(|_| draw.exponential()).collect();
    let total: f64 = weights.iter().sum();
    let mut alleles: Vec<WindowAllele> = lengths
        .iter()
        .enumerate()
        .map(|(index, repeats)| WindowAllele {
            sequence: spell(*repeats),
            repeats: *repeats,
            change: if index == 0 {
                AlleleChange::Reference
            } else {
                AlleleChange::Length
            },
            frequency: if index == 0 {
                1.0 - nonreference
            } else {
                nonreference * weights[index - 1] / total
            },
            pure_span: None,
        })
        .collect();

    let pick_donor = |alleles: &[WindowAllele], least_repeats: u32, draw: &mut Draw| {
        let donors: Vec<usize> = (0..alleles.len())
            .filter(|index| {
                alleles[*index].change.is_pure_length() && alleles[*index].repeats >= least_repeats
            })
            .collect();
        (!donors.is_empty()).then(|| donors[draw.below(donors.len() as u64) as usize])
    };

    if draw.unit() < INTERIOR_SUBSTITUTION_SHARE
        && let Some(donor) = pick_donor(&alleles, 3, draw)
    {
        let repeats = alleles[donor].repeats as usize;
        let unit =
            (repeats / 3 + draw.below((repeats / 3).max(1) as u64) as usize).clamp(1, repeats - 2);
        let at = WINDOW_FLANK_BP + unit * period + draw.below(period as u64) as usize;
        let mut sequence = alleles[donor].sequence.clone();
        sequence[at] = draw.other_base(sequence[at]);
        let repeats = alleles[donor].repeats;
        add_derived(
            &mut alleles,
            donor,
            sequence,
            repeats,
            AlleleChange::InteriorSubstitution,
            draw,
        );
    }
    if draw.unit() < FLANK_SNP_SHARE
        && let Some(donor) = pick_donor(&alleles, 1, draw)
    {
        let distance = 1 + draw.below(WINDOW_FLANK_BP as u64) as usize;
        let mut sequence = alleles[donor].sequence.clone();
        let at = if draw.unit() < 0.5 {
            WINDOW_FLANK_BP - distance
        } else {
            sequence.len() - WINDOW_FLANK_BP + distance - 1
        };
        sequence[at] = draw.other_base(sequence[at]);
        let repeats = alleles[donor].repeats;
        add_derived(
            &mut alleles,
            donor,
            sequence,
            repeats,
            AlleleChange::FlankSnp,
            draw,
        );
    }
    if draw.unit() < EDGE_DELETION_SHARE
        && let Some(donor) = pick_donor(&alleles, 2, draw)
    {
        let run_len = period * alleles[donor].repeats as usize;
        let deleted = 2 + draw.below(9) as usize; // two to ten bases
        let in_repeat = 1 + draw.below((deleted - 1).min(run_len - period) as u64) as usize;
        let in_flank = deleted - in_repeat;
        let mut sequence = alleles[donor].sequence.clone();
        let range = if draw.unit() < 0.5 {
            WINDOW_FLANK_BP - in_flank..WINDOW_FLANK_BP + in_repeat
        } else {
            WINDOW_FLANK_BP + run_len - in_repeat..WINDOW_FLANK_BP + run_len + in_flank
        };
        sequence.drain(range);
        let repeats = ((run_len - in_repeat) / period) as u32;
        add_derived(
            &mut alleles,
            donor,
            sequence,
            repeats,
            AlleleChange::EdgeDeletion,
            draw,
        );
    }

    for allele in &mut alleles {
        allele.pure_span = pure_span(&allele.sequence, period);
    }
    alleles
}

/// Add an allele derived from `donor`, taking a share of the donor's frequency so the
/// frequencies still sum to one. A derived allele equal to one already there is dropped.
fn add_derived(
    alleles: &mut Vec<WindowAllele>,
    donor: usize,
    sequence: Vec<u8>,
    repeats: u32,
    change: AlleleChange,
    draw: &mut Draw,
) {
    if alleles.iter().any(|allele| allele.sequence == sequence) {
        return;
    }
    let frequency = alleles[donor].frequency * (0.2 + 0.6 * draw.unit());
    alleles[donor].frequency -= frequency;
    alleles.push(WindowAllele {
        sequence,
        repeats,
        change,
        frequency,
        pure_span: None,
    });
}

/// The longest stretch of `sequence` that repeats with this period, searched in the middle of
/// the window where the repeat sits, or `None` if it is shorter than two units.
fn pure_span(sequence: &[u8], period: usize) -> Option<(usize, usize)> {
    let from = WINDOW_FLANK_BP.saturating_sub(period);
    let to = (sequence.len() + period).saturating_sub(WINDOW_FLANK_BP);
    let to = to.min(sequence.len());
    let mut best: Option<(usize, usize)> = None;
    let mut index = from;
    while index + period < to {
        if sequence[index] != sequence[index + period] {
            index += 1;
            continue;
        }
        let begin = index;
        while index + period < to && sequence[index] == sequence[index + period] {
            index += 1;
        }
        let span = (begin, index + period);
        if best.is_none_or(|(start, end)| span.1 - span.0 > end - start) {
            best = Some(span);
        }
    }
    best.filter(|(start, end)| end - start >= 2 * period)
}

/// One variant away from the segments: a SNP or an indel, in VCF's spelling (an indel's first
/// base is the one before it).
struct GenericVariant {
    contig: usize,
    /// Zero-based, the first base of `reference`.
    position: usize,
    reference: Vec<u8>,
    alternative: Vec<u8>,
    frequency: f64,
}

impl GenericVariant {
    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternative.len() == 1
    }
}

fn draw_generic_variants(contigs: &[Contig], draw: &mut Draw) -> Vec<GenericVariant> {
    let mut variants = Vec::new();
    for (contig_index, contig) in contigs.iter().enumerate() {
        for &(start, end) in &contig.ordinary {
            let mut position = start + GENERIC_VARIANT_SPACING;
            loop {
                position += (draw.exponential() / (SNP_RATE + INDEL_RATE)) as usize;
                if position + LONGEST_GENERIC_INDEL + GENERIC_VARIANT_SPACING >= end {
                    break;
                }
                let bases = &contig.bases;
                let (reference, alternative) = if draw.unit() < INDEL_RATE / (SNP_RATE + INDEL_RATE)
                {
                    let mut len = 1;
                    while len < LONGEST_GENERIC_INDEL && draw.unit() < 0.5 {
                        len += 1;
                    }
                    if draw.unit() < 0.5 {
                        (
                            bases[position..=position + len].to_vec(),
                            vec![bases[position]],
                        )
                    } else {
                        let mut inserted = vec![bases[position]];
                        inserted.extend((0..len).map(|_| FOUR_BASES[draw.below(4) as usize]));
                        (vec![bases[position]], inserted)
                    }
                } else {
                    (
                        vec![bases[position]],
                        vec![draw.other_base(bases[position])],
                    )
                };
                let frequency = draw.beta(NONREF_FREQUENCY_BETA.0, NONREF_FREQUENCY_BETA.1);
                let next = position + reference.len() + GENERIC_VARIANT_SPACING;
                variants.push(GenericVariant {
                    contig: contig_index,
                    position,
                    reference,
                    alternative,
                    frequency,
                });
                position = next;
            }
        }
    }
    variants
}

/// One sample's two copies at a site, as allele indices: both one draw with probability
/// `inbreeding`, two independent draws otherwise.
fn draw_pair(frequencies: &[f64], inbreeding: f64, draw: &mut Draw) -> [u8; 2] {
    let first = draw.categorical(frequencies) as u8;
    if draw.unit() < inbreeding {
        [first, first]
    } else {
        [first, draw.categorical(frequencies) as u8]
    }
}

/// One sample's library.
#[derive(Debug, Clone, Copy)]
struct Library {
    level: StutterLevel,
    multiplier: f64,
}

/// The whole drawn truth: reference, alleles, genotypes and libraries.
struct Cohort {
    contigs: Vec<Contig>,
    segments: Vec<Segment>,
    generic: Vec<GenericVariant>,
    /// `[segment][sample]`, indices into the segment's alleles, first copy first.
    segment_genotypes: Vec<Vec<[u8; 2]>>,
    /// `[variant][sample]`, 1 where the copy carries the alternative.
    generic_genotypes: Vec<Vec<[u8; 2]>>,
    samples: Vec<String>,
    libraries: Vec<Library>,
    /// Per contig, every segment and variant on it in order along the contig.
    events: Vec<Vec<Event>>,
}

#[derive(Debug, Clone, Copy)]
enum Event {
    Segment(usize),
    Generic(usize),
}

fn draw_cohort(settings: &Settings, draw: &mut Draw) -> Result<Cohort, String> {
    let (contigs, mut segments) = build_reference(settings, draw)?;
    for segment in &mut segments {
        let window = &contigs[segment.contig].bases[segment.window_start()..segment.window_end()];
        segment.alleles = draw_segment_alleles(segment, window, draw);
    }
    let generic = draw_generic_variants(&contigs, draw);
    let samples: Vec<String> = (0..settings.samples)
        .map(|index| format!("sim{index:03}"))
        .collect();
    let libraries = (0..settings.samples)
        .map(|sample| Library {
            level: settings.design.level_of(sample),
            multiplier: float::exp(LIBRARY_MULTIPLIER_LOG_SD * draw.normal()),
        })
        .collect();
    let segment_genotypes = segments
        .iter()
        .map(|segment| {
            let frequencies: Vec<f64> = segment.alleles.iter().map(|one| one.frequency).collect();
            (0..settings.samples)
                .map(|_| draw_pair(&frequencies, settings.inbreeding, draw))
                .collect()
        })
        .collect();
    let generic_genotypes = generic
        .iter()
        .map(|variant| {
            let frequencies = [1.0 - variant.frequency, variant.frequency];
            (0..settings.samples)
                .map(|_| draw_pair(&frequencies, settings.inbreeding, draw))
                .collect()
        })
        .collect();

    let mut events: Vec<Vec<(usize, Event)>> = vec![Vec::new(); contigs.len()];
    for (index, segment) in segments.iter().enumerate() {
        events[segment.contig].push((segment.window_start(), Event::Segment(index)));
    }
    for (index, variant) in generic.iter().enumerate() {
        events[variant.contig].push((variant.position, Event::Generic(index)));
    }
    let events = events
        .into_iter()
        .map(|mut on_contig| {
            on_contig.sort_by_key(|(position, _)| *position);
            on_contig.into_iter().map(|(_, event)| event).collect()
        })
        .collect();

    Ok(Cohort {
        contigs,
        segments,
        generic,
        segment_genotypes,
        generic_genotypes,
        samples,
        libraries,
        events,
    })
}

// ---------------------------------------------------------------------------
// The molecules and the reads
// ---------------------------------------------------------------------------

/// Where one copy of the genome can stutter: one segment's pure repeat, in the copy's own
/// coordinates.
#[derive(Debug, Clone, Copy)]
struct StutterSite {
    pure_start: usize,
    pure_end: usize,
    period: usize,
    repeats: u32,
}

/// One copy of one contig in one sample.
struct Haplotype {
    bases: Vec<u8>,
    /// In order along the copy.
    sites: Vec<StutterSite>,
}

fn build_haplotype(cohort: &Cohort, contig: usize, sample: usize, copy: usize) -> Haplotype {
    let reference = &cohort.contigs[contig].bases;
    let mut bases = Vec::with_capacity(reference.len() + 4096);
    let mut sites = Vec::new();
    let mut cursor = 0;
    for event in &cohort.events[contig] {
        match *event {
            Event::Segment(index) => {
                let segment = &cohort.segments[index];
                let allele =
                    &segment.alleles[cohort.segment_genotypes[index][sample][copy] as usize];
                bases.extend_from_slice(&reference[cursor..segment.window_start()]);
                let offset = bases.len();
                bases.extend_from_slice(&allele.sequence);
                if let Some((start, end)) = allele.pure_span {
                    sites.push(StutterSite {
                        pure_start: offset + start,
                        pure_end: offset + end,
                        period: segment.period(),
                        repeats: allele.repeats,
                    });
                }
                cursor = segment.window_end();
            }
            Event::Generic(index) => {
                if cohort.generic_genotypes[index][sample][copy] == 1 {
                    let variant = &cohort.generic[index];
                    bases.extend_from_slice(&reference[cursor..variant.position]);
                    bases.extend_from_slice(&variant.alternative);
                    cursor = variant.position + variant.reference.len();
                }
            }
        }
    }
    bases.extend_from_slice(&reference[cursor..]);
    Haplotype { bases, sites }
}

/// How many units one molecule's repeat moves at one site: zero when it does not stutter,
/// negative when it comes out shorter.
fn stutter_units(library: Library, period: usize, repeats: u32, draw: &mut Draw) -> i32 {
    let share = (stutter_rate(library.level, period, repeats) * library.multiplier)
        .min(HIGHEST_STUTTER_SHARE);
    if draw.unit() >= share {
        return 0;
    }
    let mut size = 1;
    while size < MAX_STUTTER_UNITS && draw.unit() >= ONE_UNIT_SHARE {
        size += 1;
    }
    if draw.unit() < library.level.shorter_share() {
        -(size as i32)
    } else {
        size as i32
    }
}

/// Copy the molecule starting at `start` on `haplotype` into `molecule`, at most `length` bases,
/// with stutter applied at every site it overlaps; returns how many sites stuttered.
///
/// **A stutter edit is placed at the start of the site's pure repeat, or at the molecule's own
/// start if that falls inside the repeat**, and is clamped to the pure stretch the molecule
/// holds: a removal never takes the last unit, an insertion copies a unit of the repeat.
fn cut_molecule(
    haplotype: &Haplotype,
    start: usize,
    length: usize,
    library: Library,
    draw: &mut Draw,
    molecule: &mut Vec<u8>,
) -> usize {
    let end = haplotype.bases.len().min(start + length + MOLECULE_SLACK);
    molecule.clear();
    molecule.extend_from_slice(&haplotype.bases[start..end]);
    let first = haplotype
        .sites
        .partition_point(|site| site.pure_end <= start);
    let overlapping = haplotype.sites[first..]
        .iter()
        .take_while(|site| site.pure_start < start + length)
        .count();
    let mut stuttered = 0;
    // Last site first, so an edit never moves the coordinates of one still to come.
    for site in haplotype.sites[first..first + overlapping].iter().rev() {
        let units = stutter_units(library, site.period, site.repeats, draw);
        if units == 0 {
            continue;
        }
        let at = site.pure_start.max(start) - start;
        let pure_end = site.pure_end.min(end) - start;
        let period = site.period;
        if units > 0 {
            if at + period > pure_end {
                continue;
            }
            let mut unit = [0_u8; 6];
            unit[..period].copy_from_slice(&molecule[at..at + period]);
            let copies = units as usize;
            molecule.splice(
                at..at,
                std::iter::repeat_n(&unit[..period], copies)
                    .flatten()
                    .copied(),
            );
        } else {
            let removable = ((pure_end - at) / period).saturating_sub(1);
            let copies = (units.unsigned_abs() as usize).min(removable);
            if copies == 0 {
                continue;
            }
            molecule.drain(at..at + copies * period);
        }
        stuttered += 1;
    }
    molecule.truncate(length);
    stuttered
}

fn reverse_complement_into(source: &[u8], target: &mut Vec<u8>) {
    target.clear();
    target.extend(source.iter().rev().map(|base| match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => *other,
    }));
}

/// Draw one read's base qualities and misread its bases at the rate they state.
fn sequence_read(
    read: &mut [u8],
    qualities: &mut Vec<u8>,
    error_probability: &[f64],
    draw: &mut Draw,
) {
    qualities.clear();
    let last = (read.len().max(2) - 1) as f64;
    for (index, base) in read.iter_mut().enumerate() {
        let fraction = index as f64 / last;
        let mean = QUALITY_MEAN - QUALITY_TAIL_DROP * fraction * fraction * fraction;
        let quality = (mean + QUALITY_SD * draw.normal())
            .round()
            .clamp(f64::from(LOWEST_QUALITY), f64::from(HIGHEST_QUALITY))
            as u8;
        if draw.unit() < error_probability[quality as usize] {
            *base = draw.other_base(*base);
        }
        qualities.push(quality + 33);
    }
}

/// What sequencing one sample produced.
#[derive(Debug, Default, Clone, Copy)]
struct SampleCounts {
    fragments: usize,
    stutter_events: usize,
}

/// Sequence one sample: build its two copies of every contig, cut `fragments` molecules from
/// them and write both mates of each.
fn sequence_sample(
    cohort: &Cohort,
    sample: usize,
    fragments: usize,
    seed: u64,
    out_dir: &Path,
) -> std::io::Result<SampleCounts> {
    let mut draw = Draw::new(seed);
    let haplotypes: Vec<Haplotype> = (0..2)
        .flat_map(|copy| {
            (0..cohort.contigs.len())
                .map(move |contig| build_haplotype(cohort, contig, sample, copy))
        })
        .collect();
    let total: u64 = haplotypes.iter().map(|one| one.bases.len() as u64).sum();
    let library = cohort.libraries[sample];
    let error_probability: Vec<f64> = (0..=HIGHEST_QUALITY)
        .map(|quality| float::powf(10.0, -f64::from(quality) / 10.0))
        .collect();

    let name = &cohort.samples[sample];
    let open = |mate: u8| -> std::io::Result<noodles_bgzf::io::Writer<std::fs::File>> {
        let path = out_dir.join(format!("{name}_R{mate}.fastq.gz"));
        Ok(noodles_bgzf::io::Writer::new(std::fs::File::create(path)?))
    };
    let mut first_out = open(1)?;
    let mut second_out = open(2)?;

    let mut counts = SampleCounts::default();
    let mut molecule = Vec::new();
    let mut flipped = Vec::new();
    let mut first = Vec::with_capacity(READ_LEN);
    let mut second = Vec::with_capacity(READ_LEN);
    let mut qualities = Vec::with_capacity(READ_LEN);
    while counts.fragments < fragments {
        let length = (FRAGMENT_MEAN + FRAGMENT_SD * draw.normal())
            .round()
            .max(READ_LEN as f64) as usize;
        let mut position = draw.below(total);
        let haplotype = haplotypes
            .iter()
            .find(|one| {
                let len = one.bases.len() as u64;
                if position < len {
                    true
                } else {
                    position -= len;
                    false
                }
            })
            .expect("a position below the total lies in some copy");
        let start = position as usize;
        if start + length > haplotype.bases.len() {
            continue;
        }
        counts.stutter_events +=
            cut_molecule(haplotype, start, length, library, &mut draw, &mut molecule);
        if draw.unit() < 0.5 {
            reverse_complement_into(&molecule, &mut flipped);
            std::mem::swap(&mut molecule, &mut flipped);
        }
        first.clear();
        first.extend_from_slice(&molecule[..READ_LEN]);
        reverse_complement_into(&molecule[molecule.len() - READ_LEN..], &mut second);

        let read_name = format!("{name}:{}", counts.fragments);
        for (read, out) in [(&mut first, &mut first_out), (&mut second, &mut second_out)] {
            sequence_read(read, &mut qualities, &error_probability, &mut draw);
            out.write_all(b"@")?;
            out.write_all(read_name.as_bytes())?;
            out.write_all(b"\n")?;
            out.write_all(read)?;
            out.write_all(b"\n+\n")?;
            out.write_all(&qualities)?;
            out.write_all(b"\n")?;
        }
        counts.fragments += 1;
    }
    first_out.try_finish()?;
    second_out.try_finish()?;
    Ok(counts)
}

// ---------------------------------------------------------------------------
// Writing the truth
// ---------------------------------------------------------------------------

fn create(path: &Path) -> std::io::Result<std::io::BufWriter<std::fs::File>> {
    Ok(std::io::BufWriter::new(std::fs::File::create(path)?))
}

fn write_fasta(path: &Path, contigs: &[Contig]) -> std::io::Result<()> {
    let mut file = create(path)?;
    for contig in contigs {
        writeln!(file, ">{}", contig.name)?;
        for chunk in contig.bases.chunks(60) {
            file.write_all(chunk)?;
            writeln!(file)?;
        }
    }
    file.flush()
}

/// Trim the bases every allele shares at the end, then at the start, keeping at least one base
/// in each — VCF's spelling. Returns the zero-based position of the first base kept.
fn trim_shared_bases(position: usize, alleles: &mut [Vec<u8>]) -> usize {
    while alleles.iter().all(|one| one.len() >= 2)
        && alleles.iter().all(|one| one.last() == alleles[0].last())
    {
        for allele in alleles.iter_mut() {
            allele.pop();
        }
    }
    let mut position = position;
    while alleles.iter().all(|one| one.len() >= 2)
        && alleles.iter().all(|one| one[0] == alleles[0][0])
    {
        for allele in alleles.iter_mut() {
            allele.remove(0);
        }
        position += 1;
    }
    position
}

/// The truth VCF: every site any sample carries a non-reference allele at, genotypes phased.
fn write_truth_vcf(out: &mut impl std::io::Write, cohort: &Cohort) -> std::io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##source=ng_str_stutter_benchmark")?;
    for contig in &cohort.contigs {
        writeln!(
            out,
            "##contig=<ID={},length={}>",
            contig.name,
            contig.bases.len()
        )?;
    }
    writeln!(
        out,
        "##INFO=<ID=SEGMENT,Number=1,Type=String,Description=\"The STR segment whose window \
         this record spells\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"The segment's motif\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=CHANGE,Number=A,Type=String,Description=\"How each alternative differs \
         from the reference window: length, substitution_inside, snp_beside or \
         deletion_across_edge\">"
    )?;
    writeln!(
        out,
        "##INFO=<ID=REPEATS,Number=R,Type=Integer,Description=\"Whole motif units in each \
         allele's repeat\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample in &cohort.samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;

    let mut line = String::new();
    for (contig_index, events) in cohort.events.iter().enumerate() {
        let contig = &cohort.contigs[contig_index];
        for event in events {
            line.clear();
            match *event {
                Event::Segment(index) => {
                    let segment = &cohort.segments[index];
                    let genotypes = &cohort.segment_genotypes[index];
                    let mut carried: Vec<u8> = genotypes.iter().flatten().copied().collect();
                    carried.sort_unstable();
                    carried.dedup();
                    carried.retain(|allele| *allele != 0);
                    if carried.is_empty() {
                        continue;
                    }
                    let mut alleles = vec![segment.alleles[0].sequence.clone()];
                    alleles.extend(
                        carried
                            .iter()
                            .map(|one| segment.alleles[*one as usize].sequence.clone()),
                    );
                    let position = trim_shared_bases(segment.window_start(), &mut alleles);
                    let changes: Vec<&str> = carried
                        .iter()
                        .map(|one| segment.alleles[*one as usize].change.name())
                        .collect();
                    let repeats: Vec<String> = std::iter::once(0)
                        .chain(carried.iter().copied())
                        .map(|one| segment.alleles[one as usize].repeats.to_string())
                        .collect();
                    write!(
                        line,
                        "{}\t{}\t{}\t{}\t{}\t.\tPASS\tSEGMENT={};MOTIF={};CHANGE={};REPEATS={}\tGT",
                        contig.name,
                        position + 1,
                        segment_name(index),
                        String::from_utf8_lossy(&alleles[0]),
                        alleles[1..]
                            .iter()
                            .map(|one| String::from_utf8_lossy(one).into_owned())
                            .collect::<Vec<_>>()
                            .join(","),
                        segment_name(index),
                        String::from_utf8_lossy(&segment.motif),
                        changes.join(","),
                        repeats.join(","),
                    )
                    .expect("writing to a String");
                    let index_in_record = |allele: u8| -> usize {
                        if allele == 0 {
                            0
                        } else {
                            1 + carried
                                .iter()
                                .position(|one| *one == allele)
                                .expect("every carried allele is listed")
                        }
                    };
                    for pair in genotypes {
                        write!(
                            line,
                            "\t{}|{}",
                            index_in_record(pair[0]),
                            index_in_record(pair[1])
                        )
                        .expect("writing to a String");
                    }
                }
                Event::Generic(index) => {
                    let variant = &cohort.generic[index];
                    let genotypes = &cohort.generic_genotypes[index];
                    if genotypes.iter().flatten().all(|copy| *copy == 0) {
                        continue;
                    }
                    write!(
                        line,
                        "{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT",
                        contig.name,
                        variant.position + 1,
                        String::from_utf8_lossy(&variant.reference),
                        String::from_utf8_lossy(&variant.alternative),
                    )
                    .expect("writing to a String");
                    for pair in genotypes {
                        write!(line, "\t{}|{}", pair[0], pair[1]).expect("writing to a String");
                    }
                }
            }
            writeln!(out, "{line}")?;
        }
    }
    out.flush()
}

fn write_tables(out_dir: &Path, cohort: &Cohort) -> std::io::Result<()> {
    let mut bed = create(&out_dir.join("segments.bed"))?;
    let mut alleles = create(&out_dir.join("segment_alleles.tsv"))?;
    let mut genotypes = create(&out_dir.join("segment_genotypes.tsv"))?;
    writeln!(
        alleles,
        "segment\tchrom\twindow_start\twindow_end\tsegment_start\tsegment_end\tperiod\tmotif\t\
         reference_repeats\tallele\trepeats\tchange\tfrequency\tsequence"
    )?;
    writeln!(genotypes, "segment\tsample\tfirst_copy\tsecond_copy")?;
    for (index, segment) in cohort.segments.iter().enumerate() {
        let contig = &cohort.contigs[segment.contig].name;
        writeln!(
            bed,
            "{contig}\t{}\t{}\t{}",
            segment.start,
            segment.end(),
            segment.period()
        )?;
        for (allele_index, allele) in segment.alleles.iter().enumerate() {
            writeln!(
                alleles,
                "{}\t{contig}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{allele_index}\t{}\t{}\t{:.6}\t{}",
                segment_name(index),
                segment.window_start(),
                segment.window_end(),
                segment.start,
                segment.end(),
                segment.period(),
                String::from_utf8_lossy(&segment.motif),
                segment.reference_repeats,
                allele.repeats,
                allele.change.name(),
                allele.frequency,
                String::from_utf8_lossy(&allele.sequence),
            )?;
        }
        for (sample, pair) in cohort.samples.iter().zip(&cohort.segment_genotypes[index]) {
            writeln!(
                genotypes,
                "{}\t{sample}\t{}\t{}",
                segment_name(index),
                pair[0],
                pair[1]
            )?;
        }
    }
    bed.flush()?;
    alleles.flush()?;
    genotypes.flush()?;

    let mut confident = create(&out_dir.join("confident.bed"))?;
    for contig in &cohort.contigs {
        writeln!(confident, "{}\t0\t{}", contig.name, contig.bases.len())?;
    }
    confident.flush()?;

    let mut libraries = create(&out_dir.join("libraries.tsv"))?;
    writeln!(libraries, "sample\tlevel\tmultiplier\tshorter_share")?;
    for (sample, library) in cohort.samples.iter().zip(&cohort.libraries) {
        writeln!(
            libraries,
            "{sample}\t{}\t{:.6}\t{}",
            library.level.name(),
            library.multiplier,
            library.level.shorter_share()
        )?;
    }
    libraries.flush()?;

    let mut rates = create(&out_dir.join("stutter_rates.tsv"))?;
    writeln!(
        rates,
        "level\tperiod\trepeats\tshare_of_molecules_that_stutter"
    )?;
    for level in [
        StutterLevel::PcrFree,
        StutterLevel::Pcr,
        StutterLevel::HeavyPcr,
    ] {
        for period in 1..=6 {
            for repeats in 1..=40 {
                writeln!(
                    rates,
                    "{}\t{period}\t{repeats}\t{:.6}",
                    level.name(),
                    stutter_rate(level, period, repeats)
                )?;
            }
        }
    }
    rates.flush()
}

// ---------------------------------------------------------------------------

fn run(out_dir: &Path, settings: Settings) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(out_dir)?;
    let mut draw = Draw::new(settings.seed);
    let cohort = draw_cohort(&settings, &mut draw)?;

    write_fasta(&out_dir.join("reference.fa"), &cohort.contigs)?;
    write_truth_vcf(&mut create(&out_dir.join("truth.vcf"))?, &cohort)?;
    write_tables(out_dir, &cohort)?;

    let genome_len: usize = cohort.contigs.iter().map(|one| one.bases.len()).sum();
    let fragments = (settings.depth * genome_len as f64 / (2 * READ_LEN) as f64).round() as usize;
    let next = AtomicUsize::new(0);
    let threads = settings.threads.min(settings.samples);
    let counts: Vec<(usize, SampleCounts)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| -> std::io::Result<Vec<(usize, SampleCounts)>> {
                    let mut done = Vec::new();
                    loop {
                        let sample = next.fetch_add(1, Ordering::Relaxed);
                        if sample >= settings.samples {
                            return Ok(done);
                        }
                        let seed = sample_seed(settings.seed, sample);
                        let one = sequence_sample(&cohort, sample, fragments, seed, out_dir)?;
                        done.push((sample, one));
                    }
                })
            })
            .collect();
        workers
            .into_iter()
            .map(|worker| worker.join().expect("a sequencing thread panicked"))
            .collect::<std::io::Result<Vec<_>>>()
            .map(|all| all.into_iter().flatten().collect())
    })?;

    let carried_alleles: usize = cohort
        .segments
        .iter()
        .zip(&cohort.segment_genotypes)
        .map(|(segment, genotypes)| {
            (1..segment.alleles.len())
                .filter(|allele| {
                    genotypes
                        .iter()
                        .flatten()
                        .any(|one| usize::from(*one) == *allele)
                })
                .count()
        })
        .sum();
    let varying_segments = cohort
        .segment_genotypes
        .iter()
        .filter(|genotypes| genotypes.iter().flatten().any(|one| *one != 0))
        .count();
    let variant_cells: usize = cohort
        .segment_genotypes
        .iter()
        .flatten()
        .filter(|pair| pair.iter().any(|one| *one != 0))
        .count();
    let stutter_events: usize = counts.iter().map(|(_, one)| one.stutter_events).sum();
    println!("out dir             : {}", out_dir.display());
    println!(
        "reference           : {} contigs, {genome_len} bases",
        cohort.contigs.len()
    );
    println!("STR segments        : {}", cohort.segments.len());
    println!(
        "segments any sample varies at : {varying_segments} of {}, carrying {carried_alleles} \
         non-reference alleles",
        cohort.segments.len()
    );
    println!(
        "(sample, segment) cells carrying a non-reference allele : {variant_cells} of {}",
        cohort.segments.len() * cohort.samples.len()
    );
    println!(
        "variants away from segments : {} drawn ({} SNPs)",
        cohort.generic.len(),
        cohort.generic.iter().filter(|one| one.is_snp()).count()
    );
    println!(
        "fragments a sample  : {fragments} ({} reads)",
        2 * fragments
    );
    println!("stutter events      : {stutter_events} over all samples");
    println!("settings            : {settings:?}");
    Ok(())
}

fn main() -> ExitCode {
    let arguments: Vec<String> = std::env::args().skip(1).collect();
    let Some(out_dir) = arguments.first().map(PathBuf::from) else {
        eprintln!("usage: ng_str_stutter_benchmark <out-dir> [key=value ...]");
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Whether `observed` successes in `trials` sit within four standard errors of `expected`.
    /// Four, because the tests below check about a hundred shares and three would fail one of
    /// them by chance every few runs of a changed seed.
    fn within_sampling_error(observed: usize, trials: usize, expected: f64) -> bool {
        let share = observed as f64 / trials as f64;
        let standard_error = (expected * (1.0 - expected) / trials as f64).sqrt();
        (share - expected).abs() <= 4.0 * standard_error + 1.0 / trials as f64
    }

    #[test]
    fn the_rate_is_a_thirtieth_at_the_floor_and_the_plateau_from_its_repeat_count() {
        for level in [
            StutterLevel::PcrFree,
            StutterLevel::Pcr,
            StutterLevel::HeavyPcr,
        ] {
            for period in 1..=6 {
                let plateau = level.plateau(period);
                let floor = stutter_rate(level, period, CALLING_FLOOR[period - 1]);
                assert!((floor - plateau / 30.0).abs() < 1e-12, "{level:?} {period}");
                let at = PLATEAU_REPEATS[period - 1];
                assert_eq!(stutter_rate(level, period, at), plateau);
                assert_eq!(stutter_rate(level, period, at + 10), plateau);
                assert!(stutter_rate(level, period, at - 1) < plateau);
            }
        }
        assert_eq!(stutter_rate(StutterLevel::None, 1, 30), 0.0);
    }

    /// **Counted stutter within sampling error of the table**: share stuttered, share shorter
    /// and the size of each event, for every level and motif length at four repeat counts.
    #[test]
    fn counted_stutter_matches_the_table() {
        let trials = 200_000;
        let mut draw = Draw::new(7);
        for level in [
            StutterLevel::PcrFree,
            StutterLevel::Pcr,
            StutterLevel::HeavyPcr,
        ] {
            let library = Library {
                level,
                multiplier: 1.0,
            };
            let mut stuttered_all = 0;
            let mut shorter = 0;
            let mut by_size = [0_usize; MAX_STUTTER_UNITS as usize + 1];
            for period in 1..=6 {
                let floor = CALLING_FLOOR[period - 1];
                let plateau = PLATEAU_REPEATS[period - 1];
                for repeats in [floor, (floor + plateau) / 2, plateau, plateau + 5] {
                    let mut stuttered = 0;
                    for _ in 0..trials {
                        let units = stutter_units(library, period, repeats, &mut draw);
                        if units != 0 {
                            stuttered += 1;
                            shorter += usize::from(units < 0);
                            by_size[units.unsigned_abs() as usize] += 1;
                        }
                    }
                    stuttered_all += stuttered;
                    let expected = stutter_rate(level, period, repeats);
                    assert!(
                        within_sampling_error(stuttered, trials, expected),
                        "{level:?} period {period} repeats {repeats}: {stuttered} of {trials} \
                         stuttered against {expected}"
                    );
                }
            }
            assert!(
                within_sampling_error(shorter, stuttered_all, level.shorter_share()),
                "{level:?}: {shorter} of {stuttered_all} shorter"
            );
            let tail = 1.0 - ONE_UNIT_SHARE;
            for (size, moved) in by_size.iter().enumerate().skip(1) {
                let expected = if size < MAX_STUTTER_UNITS as usize {
                    ONE_UNIT_SHARE * float::powi(tail, size as i32 - 1)
                } else {
                    float::powi(tail, size as i32 - 1)
                };
                assert!(
                    within_sampling_error(*moved, stuttered_all, expected),
                    "{level:?}: {moved} of {stuttered_all} moved {size} units, expected {expected}"
                );
            }
        }
    }

    /// A library's multiplier scales the rate, and the result is capped.
    #[test]
    fn the_multiplier_scales_the_rate() {
        let trials = 200_000;
        let mut draw = Draw::new(11);
        let library = Library {
            level: StutterLevel::Pcr,
            multiplier: 2.0,
        };
        let stuttered = (0..trials)
            .filter(|_| stutter_units(library, 2, 20, &mut draw) != 0)
            .count();
        assert!(within_sampling_error(stuttered, trials, 0.30));
        let capped = Library {
            level: StutterLevel::HeavyPcr,
            multiplier: 5.0,
        };
        let stuttered = (0..trials)
            .filter(|_| stutter_units(capped, 1, 20, &mut draw) != 0)
            .count();
        assert!(within_sampling_error(
            stuttered,
            trials,
            HIGHEST_STUTTER_SHARE
        ));
    }

    /// Through the molecule cutter, as the reads see it: a molecule spanning a pure segment
    /// whole comes out longer or shorter by whole units at the table's rate.
    #[test]
    fn molecules_stutter_at_the_table_rate() {
        let mut draw = Draw::new(3);
        let mut bases = Vec::new();
        append_aperiodic(&mut bases, 120, &mut draw);
        append_spacer(&mut bases, b"AC");
        let pure_start = bases.len();
        bases.extend(b"AC".repeat(18));
        let pure_end = bases.len();
        append_spacer(&mut bases, b"AC");
        append_aperiodic(&mut bases, 120, &mut draw);
        let haplotype = Haplotype {
            sites: vec![StutterSite {
                pure_start,
                pure_end,
                period: 2,
                repeats: 18,
            }],
            bases,
        };
        let library = Library {
            level: StutterLevel::Pcr,
            multiplier: 1.0,
        };
        let trials = 100_000;
        let mut molecule = Vec::new();
        let mut changed = 0;
        let mut shorter = 0;
        for _ in 0..trials {
            let full = haplotype.bases.len() + MOLECULE_SLACK;
            let events = cut_molecule(&haplotype, 0, full, library, &mut draw, &mut molecule);
            let difference = molecule.len() as i64 - haplotype.bases.len() as i64;
            assert_eq!(difference % 2, 0, "a stutter moved a partial unit");
            assert_eq!(events, usize::from(difference != 0));
            changed += usize::from(difference != 0);
            shorter += usize::from(difference < 0);
            assert!(molecule[..pure_start] == haplotype.bases[..pure_start]);
        }
        assert!(within_sampling_error(changed, trials, 0.15), "{changed}");
        assert!(
            within_sampling_error(shorter, changed, 0.78),
            "{shorter} of {changed}"
        );
    }

    fn small_cohort(samples: usize) -> Cohort {
        let settings = Settings {
            samples,
            segments_per_period: 40,
            genome_bp: 300_000,
            contigs: 2,
            inbreeding: 0.3,
            ..Settings::default()
        };
        draw_cohort(&settings, &mut Draw::new(settings.seed)).expect("the settings fit")
    }

    /// **The truth VCF reproduces the truth sequences**: applying its phased records to the
    /// reference gives each sample's two copies of every contig exactly, the copies the reads
    /// were cut from — and so every segment's window as the allele table spells it.
    #[test]
    fn truth_vcf_reproduces_the_truth_sequences() {
        let cohort = small_cohort(6);
        let mut vcf = Vec::new();
        write_truth_vcf(&mut vcf, &cohort).expect("writing to memory");
        let text = String::from_utf8(vcf).expect("ASCII");

        let contig_index = |name: &str| {
            cohort
                .contigs
                .iter()
                .position(|one| one.name == name)
                .expect("a declared contig")
        };
        let mut records: Vec<Vec<&str>> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect())
            .collect();
        assert!(records.len() > 100, "{} records", records.len());
        records
            .sort_by_key(|fields| (contig_index(fields[0]), fields[1].parse::<usize>().unwrap()));
        let mut segment_records = 0;
        for sample in 0..cohort.samples.len() {
            for copy in 0..2 {
                for (contig, reference) in cohort.contigs.iter().enumerate() {
                    let mut rebuilt = Vec::new();
                    let mut cursor = 0;
                    for fields in records.iter().filter(|one| contig_index(one[0]) == contig) {
                        let position = fields[1].parse::<usize>().unwrap() - 1;
                        let alleles: Vec<&str> = std::iter::once(fields[3])
                            .chain(fields[4].split(','))
                            .collect();
                        let called: usize = fields[9 + sample]
                            .split('|')
                            .nth(copy)
                            .unwrap()
                            .parse()
                            .unwrap();
                        assert!(position >= cursor, "records overlap at {position}");
                        if called != 0 {
                            rebuilt.extend_from_slice(&reference.bases[cursor..position]);
                            rebuilt.extend_from_slice(alleles[called].as_bytes());
                            cursor = position + fields[3].len();
                        }
                        segment_records += usize::from(fields[2] != ".");
                    }
                    rebuilt.extend_from_slice(&reference.bases[cursor..]);
                    let built = build_haplotype(&cohort, contig, sample, copy);
                    assert!(
                        rebuilt == built.bases,
                        "sample {sample} copy {copy} {}",
                        reference.name
                    );
                }
            }
        }
        assert!(segment_records > 0);
    }

    /// Every allele drawn is distinct, the frequencies sum to one, each kind of derived allele
    /// turns up, and an allele's pure stretch is its repeat.
    #[test]
    fn alleles_are_distinct_and_their_frequencies_sum_to_one() {
        let cohort = small_cohort(1);
        let mut kinds = std::collections::HashSet::new();
        for segment in &cohort.segments {
            let total: f64 = segment.alleles.iter().map(|one| one.frequency).sum();
            assert!((total - 1.0).abs() < 1e-9);
            for (index, allele) in segment.alleles.iter().enumerate() {
                kinds.insert(allele.change);
                assert!(
                    segment.alleles[..index]
                        .iter()
                        .all(|earlier| earlier.sequence != allele.sequence)
                );
                if allele.change.is_pure_length() && allele.repeats >= 2 {
                    let (start, end) = allele.pure_span.expect("a pure repeat has a span");
                    assert!(end - start >= segment.period() * allele.repeats as usize);
                }
            }
            assert_eq!(segment.alleles[0].change, AlleleChange::Reference);
            let window =
                &cohort.contigs[segment.contig].bases[segment.window_start()..segment.window_end()];
            assert_eq!(segment.alleles[0].sequence, window);
        }
        assert_eq!(kinds.len(), 5, "{kinds:?}");
    }

    #[test]
    fn motifs_are_primitive_and_of_their_length() {
        for (index, motifs) in MOTIFS.iter().enumerate() {
            for motif in *motifs {
                let bytes = motif.as_bytes();
                assert_eq!(bytes.len(), index + 1, "{motif}");
                for unit in 1..bytes.len() {
                    if bytes.len() % unit == 0 {
                        assert!(
                            bytes[..unit].repeat(bytes.len() / unit) != bytes,
                            "{motif} repeats {}",
                            String::from_utf8_lossy(&bytes[..unit])
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn the_ordinary_sequence_keeps_its_repeat_caps() {
        let mut draw = Draw::new(5);
        let mut bases = Vec::new();
        append_ordinary(&mut bases, 50_000, &mut draw);
        let longest_homopolymer = bases
            .chunk_by(|a, b| a == b)
            .map(<[u8]>::len)
            .max()
            .unwrap();
        assert_eq!(longest_homopolymer, GROUND_MAX_COPIES[0]);
        for period in 2..=6 {
            for start in 0..bases.len() - 3 * period {
                let window = &bases[start..start + 3 * period];
                assert!(
                    (0..2 * period).any(|index| window[index] != window[index + period]),
                    "three copies of period {period} at {start}"
                );
            }
        }
    }

    #[test]
    fn settings_refuse_what_they_do_not_know() {
        assert!(Settings::parse(&["sample=3".to_string()]).is_err());
        assert!(Settings::parse(&["design=pcrfree".to_string()]).is_err());
        assert!(Settings::parse(&["design=one_outlier".to_string()]).is_err());
        let settings =
            Settings::parse(&["design=one_outlier".to_string(), "samples=3".to_string()]).unwrap();
        assert_eq!(settings.design.level_of(0), StutterLevel::HeavyPcr);
        assert_eq!(settings.design.level_of(1), StutterLevel::PcrFree);
    }
}
