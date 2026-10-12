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
//!   allele. An STR bundle is likewise one record over its window (the hull and 30 bases either
//!   side), with `BUNDLE`, `MOTIFS` and `MEMBER_REPEATS` in place of `SEGMENT`, `MOTIF` and
//!   `REPEATS`; those three keys are declared in the header only when the file has a bundle.
//!   Variants away from the segments are one record each.
//! - `segments.bed` — `chrom start end period`, zero-based half-open, one row per STR segment:
//!   the ground the typing must call STR segments, and nothing else.
//! - `segment_alleles.tsv` — every allele drawn at every segment, carried or not, as its whole
//!   window sequence, with its repeat count, how it differs from the reference and its
//!   population frequency.
//! - `segment_genotypes.tsv` — every sample's two alleles at every segment, as indices into
//!   `segment_alleles.tsv`, including the homozygous-reference ones. With the allele table this
//!   is each sample's two truth sequences at every segment.
//! - `bundles.bed`, `bundle_alleles.tsv`, `bundle_genotypes.tsv` — the same three for the STR
//!   bundles (below), with every member's motif, span and repeat count. Kept apart from the
//!   segments' files so that those keep their columns and every reader of them keeps reading
//!   segments only. Written with their header alone when `bundles=0`, and then every other file
//!   is byte for byte what the simulator wrote before bundles existed; with bundles laid, the
//!   segments sit elsewhere along the genome, so a default run is a different cohort.
//! - `bundle_runs.tsv` — each bundle's shape (`separated`, `abutting`, `interrupted` or
//!   `abutting_interrupted`) and where its reference stutters, as `motif:offset:length` from the
//!   hull's first base: the format `ng_str_bundle_runs` writes the caller's reference runs in, so
//!   the typing check compares the two. A file of its own, so that with `hard_bundle_share=0` every
//!   other file is byte for byte what the simulator wrote before hard shapes existed.
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
//!    and separated by ordinary sequence. Beside them, `bundles` **STR bundles**: two or three
//!    repeats (*members*) separated by gaps of 1 to 14 bases, so close that none of them has
//!    clean flanks of its own, the whole stretch from the first member's start to the last one's
//!    end (the *hull*) at most 100 bases and wrapped in 150 bases of clean flank as a segment is
//!    (`doc/devel/ng/spec/str_bundle_loci.md` §9). A share `hard_bundle_share` of them are laid in
//!    a **hard shape** instead: repeats that abut, often a repeat running into a homopolymer of its
//!    main base (`(TTTC)4 T(11)`), or repeats split into pure pieces by an interruption — the
//!    shapes six in ten of HG002's bundles have. Each pure piece is then a member.
//! 2. **The population**: at each segment two to six alleles differing by whole units, the
//!    non-reference ones sharing a frequency drawn from a Beta(0.5, 1.5); one segment in ten
//!    gains an allele with a substitution inside the repeat, one in ten an allele with a SNP
//!    within 30 bases of an edge, one in thirty an allele with a deletion crossing an edge.
//!    Away from the segments, SNPs at one per thousand bases and indels of one to ten bases at
//!    one per ten thousand, each with its own Beta frequency. A bundle's alleles are drawn the
//!    same way, the length changes per member (one member at a time, a second one too in one
//!    allele in four), and the derived alleles carry a substitution inside a member, a SNP
//!    beside the hull or a deletion across the hull's outer edge.
//! 3. **The genotypes**: each sample's two copies drawn from those frequencies, both copies one
//!    draw with probability `inbreeding`.
//! 4. **The libraries**: each sample's stutter level comes from the `design`, and its multiplier
//!    from a log-normal with a standard deviation of 0.4 on the log scale (median one).
//! 5. **The molecules**: fragments of length drawn around 350 (sd 50), placed uniformly over
//!    the sample's two copies of the genome. At every STR segment a molecule overlaps it gains
//!    or loses whole units: it stutters with probability `rate(level, motif length, repeats) ×
//!    multiplier`, shorter with the level's shorter share, by one unit with probability 0.85,
//!    two with 0.85 × 0.15 and so on, at most four. **Stutter acts on the molecule, so both mates
//!    of a pair report the same length**, as both read one PCR product. A separated bundle
//!    stutters in every member independently, each at the rate for its own motif length and its
//!    own repeat count in the allele the molecule copies; a hard-shaped one in every **run** of the
//!    allele's sequence — each stretch of pure whole copies of one motif at or above its calling
//!    floor, found by the rule the caller uses for a bundle's reference runs — each at its own rate.
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
//! | `bundles` | 500 | STR bundles, all motif lengths together; `0` lays none and leaves every other output as it was without them |
//! | `hard_bundle_share` | 0.6 | the share of bundles laid in a hard shape (HG002's); `0` lays every bundle as before hard shapes existed, byte for byte |
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

/// The bundles laid by default: one for every six segments, a few hundred, so that a cohort's
/// bundle cells are counted in thousands while the segments stay the bulk of the genome.
const DEFAULT_BUNDLES: usize = 500;

/// The share of bundles with three members; the rest have two. HG002's tandem-repeat regions
/// hold 229 three-member bundles in 2,473 (`str_bundle_loci.md` §1), about one in eleven; this
/// is set higher so that three-member bundles are numerous enough to count.
const THREE_MEMBER_SHARE: f64 = 0.15;

/// How often a bundle member has each motif length, one to six, in hundredths: HG002's shares.
/// Of the 5,067 members of the 2,410 bundles of at most 100 bases that the typing finds in HG002's
/// tandem-repeat regions, 2,540, 1,519, 249, 390, 167 and 202 have motif lengths one to six
/// (measured 2026-10-11, `str_stutter_benchmark.md` §11.1). Bundles are mostly short repeats side by
/// side; drawn uniform over one to six, two members at their floors and a gap of the mean length
/// would already make a hull of 44 bases, against HG002's median of 34.
const BUNDLE_MOTIF_LENGTH_WEIGHTS: [u64; 6] = [50, 30, 5, 8, 3, 4];

/// The share of mono- and dinucleotide members drawn **long**: a repeat count uniform from
/// [`LONG_MEMBER_LEAST_REPEATS`] to the segments' highest, 30. The rest of the members sit near
/// their floor ([`BUNDLE_MEAN_EXTRA_COPIES`]). Chosen so that about one bundle in eight holds a
/// mono- or dinucleotide run of 16 copies or more, as in HG002 (284 of 2,410 bundles of at most 100
/// bases, counting pure runs) — the members where stutter does its damage. Of 10,000 bundles drawn
/// with the default seed, 1,164 hold one.
const LONG_MEMBER_SHARE: f64 = 0.07;

/// The fewest copies of a long member ([`LONG_MEMBER_SHARE`]).
const LONG_MEMBER_LEAST_REPEATS: u32 = 16;

/// The mean number of copies a member that is not long has above its calling floor: a geometric
/// count, so about three in four of these members sit at the floor and one in six a copy above it.
/// Far lower than HG002's (about 2.6 extra pure copies for mono- and dinucleotides) because the
/// simulator cannot lay what shortens HG002's hulls — members that overlap or abut, in 898 of its
/// 2,410 bundles — and its gaps, uniform over 1 to 14 bases, average 7.5 against HG002's 6.2. Chosen
/// so that the median hull comes near HG002's 34 bases: 36 over 10,000 bundles drawn with the
/// default seed, against 37 with a mean of 0.5 (and 8 long members in 100).
const BUNDLE_MEAN_EXTRA_COPIES: f64 = 0.3;

/// The longest gap between two members, in bases: one less than the typing's bundle radius
/// (`DEFAULT_BUNDLE_THRESHOLD`, 15), the widest gap at which two repeats are still one bundle.
const LONGEST_BUNDLE_GAP: usize = 14;

/// The longest hull laid, in bases: the default `--max-str-bundle-len` of
/// `str_bundle_loci.md` §5.2, above which a bundle is typed as a satellite.
const LONGEST_BUNDLE_HULL: usize = 100;

/// The scores [`scanner_repeats`] gives a base that equals the one a motif length before it, and
/// takes from one that does not: copies of `tandem_repeat.rs`'s `DEFAULT_MATCH_REWARD` and
/// `DEFAULT_MISMATCH_PENALTY`, held here so the bundles drawn do not change when the library's do.
/// The test `the_scanner_copy_uses_the_library_scores` fails when the two part.
const SCANNER_MATCH_REWARD: i64 = 2;
const SCANNER_MISMATCH_PENALTY: i64 = 7;

/// The redraws [`draw_bundle_shape`] makes before it gives up, far more than today's rules need; the
/// bound turns a rule that can no longer be met into a message rather than a hang.
const BUNDLE_SHAPE_ATTEMPTS: usize = 100_000;

/// The share of a bundle's length-change alleles that change a second member as well as the
/// first.
const SECOND_MEMBER_CHANGE_SHARE: f64 = 0.25;

/// Mixed into the seed for the bundles' own draw, so that laying bundles takes nothing from the
/// draw everything else comes from, and `bundles=0` leaves that draw untouched.
const BUNDLE_SEED_SALT: u64 = 0x5851_F42D_4C95_7F2D;

/// The share of bundles laid in a **hard shape** — members that abut or carry an interruption —
/// by default: HG002's. Of the 2,410 bundles of at most 100 bases the typing finds in HG002's
/// tandem-repeat regions, 1,447 are not the separated shape (each member exactly periodic, no two
/// members touching or overlapping, no two reference runs abutting); measured 2026-10-12, the
/// classification and its numbers recorded in `str_stutter_benchmark.md` §11.4 (the script that
/// measured it was scratch and is not in the repository).
const DEFAULT_HARD_BUNDLE_SHARE: f64 = 0.6;

/// In a hard shape, the chance that two neighbouring members abut, with no base between them,
/// instead of being separated by a gap, when the first is not a homopolymer.
const ABUTTING_JUNCTION_SHARE: f64 = 0.4;

/// How much less often a junction abuts after a homopolymer. In HG002, 19 of the 191 pairs of
/// abutting reference runs start with a homopolymer, although homopolymers are half of all bundle
/// members.
const AFTER_HOMOPOLYMER_ABUTTING_FACTOR: f64 = 0.1;

/// When two members abut and the first is not a homopolymer, the chance the second is a homopolymer
/// of the first's most frequent base — the `(TTTC)n T(m)` shape. In HG002, 118 of the 191 pairs of
/// abutting reference runs are a repeat and a homopolymer of one of its bases; and a repeat mostly
/// of one base runs on into such a homopolymer in the scanner's reading, which is where most of
/// HG002's nested members come from.
const HOMOPOLYMER_TAIL_SHARE: f64 = 0.8;

/// In a hard shape, the chance that each member is interrupted: split into pure pieces of its motif
/// by a substituted copy, an inserted base or a copy that lost a base.
const INTERRUPTED_MEMBER_SHARE: f64 = 0.45;

/// How many interruptions an interrupted member carries — one, two or three — as weights: HG002's
/// impure bundle members hold 852 with one cluster of bases off their period, 414 with two and 199
/// with three or more.
const INTERRUPTION_COUNT_WEIGHTS: [f64; 3] = [852.0, 414.0, 199.0];

/// The mean number of whole copies, geometric, an interrupted member gains on top of its drawn
/// count, plus one. HG002's impure members are longer than its pure ones — a median of 12 copies
/// for mono- and dinucleotides against 9 and 6 — and a member split at the floor would leave
/// pieces too short for the typing to see as a repeat at all.
const INTERRUPTED_MEMBER_MEAN_EXTRA_COPIES: f64 = 1.0;

/// The most extra copies an interrupted member gains on top of its drawn count. The hull cap
/// ([`LONGEST_BUNDLE_HULL`]) binds long before it; it only bounds the geometric draw.
const INTERRUPTED_MEMBER_MOST_EXTRA_COPIES: u32 = 30;

/// The fewest copies of the short repeat a homopolymer tail follows (`(TTC)4 T9`); it is drawn
/// from this up to its calling floor, as HG002's are.
const TAIL_LEAD_LEAST_COPIES: u32 = 3;

/// The share of hard bundles laid with **no reference run**: every repeat impure, no stretch of
/// pure copies at its calling floor. HG002: 123 of the 1,447 hard bundles
/// (`str_stutter_benchmark.md` §11.4). Without it the simulator laid 6 in 3,000.
const NO_RUN_HARD_BUNDLE_SHARE: f64 = 123.0 / 1447.0;

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
    bundles: usize,
    /// The share of bundles laid in a hard shape ([`draw_hard_bundle_shape`]); the rest are laid
    /// separated ([`draw_bundle_shape`]). `0` lays every bundle as before hard shapes existed.
    hard_bundle_share: f64,
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
            bundles: DEFAULT_BUNDLES,
            hard_bundle_share: DEFAULT_HARD_BUNDLE_SHARE,
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
                "bundles" => settings.bundles = number()? as usize,
                "hard_bundle_share" => settings.hard_bundle_share = number()?,
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
        if !(0.0..=1.0).contains(&settings.hard_bundle_share) {
            return Err(format!(
                "hard_bundle_share must be in [0, 1], got {}",
                settings.hard_bundle_share
            ));
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

/// **An STR bundle's reference sequence, before it is placed**: its members in order along the
/// contig, each a motif and how many whole copies of it the reference holds, and the bases
/// between consecutive members.
///
/// In a hard shape a *member* is one pure piece: an interrupted repeat is several members of one
/// motif, with the interruption as the bases between them, and two abutting repeats have no base
/// between them. So a member is always pure whole copies, which is what a length change moves.
#[derive(Debug, Clone, PartialEq, Eq)]
struct BundleShape {
    /// Each member's motif and reference repeat count, first member first.
    members: Vec<(Vec<u8>, u32)>,
    /// The bases after each member but the last. In the separated shape, one to
    /// [`LONGEST_BUNDLE_GAP`] bases, none of them continuing either neighbour's motif; in a hard
    /// shape also none (two members that abut) or an interruption (a mangled copy, or one stray
    /// base, between two pieces of one repeat).
    gaps: Vec<Vec<u8>>,
    kind: BundleShapeKind,
}

/// Which of the shapes [`draw_bundle_shape`] and [`draw_hard_bundle_shape`] lay a bundle is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum BundleShapeKind {
    /// Pure members separated by gaps, each seen by the typing exactly as drawn: the shape the
    /// simulator laid before hard shapes existed.
    Separated,
    /// Some neighbouring members abut, none is interrupted.
    Abutting,
    /// Some member is interrupted, none abut.
    Interrupted,
    /// Both.
    AbuttingAndInterrupted,
}

impl BundleShapeKind {
    fn name(self) -> &'static str {
        match self {
            Self::Separated => "separated",
            Self::Abutting => "abutting",
            Self::Interrupted => "interrupted",
            Self::AbuttingAndInterrupted => "abutting_interrupted",
        }
    }
}

impl BundleShape {
    fn reference_repeats(&self) -> Vec<u32> {
        self.members.iter().map(|(_, repeats)| *repeats).collect()
    }

    fn periods(&self) -> Vec<usize> {
        self.members.iter().map(|(motif, _)| motif.len()).collect()
    }

    /// The hull spelled with `repeats[m]` whole copies of member `m`'s motif and the reference
    /// gaps between them.
    fn spell(&self, repeats: &[u32]) -> Vec<u8> {
        let mut hull = Vec::new();
        for (index, ((motif, _), copies)) in self.members.iter().zip(repeats).enumerate() {
            hull.extend(motif.repeat(*copies as usize));
            if let Some(gap) = self.gaps.get(index) {
                hull.extend_from_slice(gap);
            }
        }
        hull
    }

    /// Where each member's run lies in [`Self::spell`]'s hull for these repeat counts:
    /// zero-based half-open offsets from the hull's first base.
    fn member_runs(&self, repeats: &[u32]) -> Vec<(usize, usize)> {
        let mut runs = Vec::with_capacity(self.members.len());
        let mut cursor = 0;
        for (index, ((motif, _), copies)) in self.members.iter().zip(repeats).enumerate() {
            let end = cursor + motif.len() * *copies as usize;
            runs.push((cursor, end));
            cursor = end + self.gaps.get(index).map_or(0, Vec::len);
        }
        runs
    }

    fn hull_len(&self) -> usize {
        let members: usize = self
            .members
            .iter()
            .map(|(motif, repeats)| motif.len() * *repeats as usize)
            .sum();
        members + self.gaps.iter().map(Vec::len).sum::<usize>()
    }

    fn first_motif(&self) -> &[u8] {
        &self.members[0].0
    }

    fn last_motif(&self) -> &[u8] {
        &self.members[self.members.len() - 1].0
    }
}

/// One STR bundle as laid down, and the alleles the population carries at it.
struct Bundle {
    contig: usize,
    shape: BundleShape,
    /// Zero-based, the first base of the first member — the hull's start.
    start: usize,
    /// Allele 0 is the reference.
    alleles: Vec<BundleAllele>,
}

impl Bundle {
    /// Zero-based, one past the last member's last base.
    fn end(&self) -> usize {
        self.start + self.shape.hull_len()
    }

    fn window_start(&self) -> usize {
        self.start - WINDOW_FLANK_BP
    }

    fn window_end(&self) -> usize {
        self.end() + WINDOW_FLANK_BP
    }

    /// Each member's reference span on the contig, zero-based half-open.
    fn member_spans(&self) -> Vec<(usize, usize)> {
        self.shape
            .member_runs(&self.shape.reference_repeats())
            .into_iter()
            .map(|(start, end)| (self.start + start, self.start + end))
            .collect()
    }

    /// The members' motifs, comma-separated, first member first.
    fn motifs(&self) -> String {
        self.shape
            .members
            .iter()
            .map(|(motif, _)| String::from_utf8_lossy(motif).into_owned())
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn bundle_name(index: usize) -> String {
    format!("bun{index:05}")
}

/// Join numbers with `separator` — a bundle's per-member values in its tables and record.
fn joined<T: std::fmt::Display>(values: impl IntoIterator<Item = T>, separator: &str) -> String {
    values
        .into_iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(separator)
}

/// One STR bundle's shape: members, motifs, repeat counts and gap bases, drawn until every rule
/// below holds — a redraw, so a rule shapes the distribution rather than bending one draw.
///
/// - **Two members, or three** with probability [`THREE_MEMBER_SHARE`], drawn once; the rules
///   below redraw the members and gaps, never their number. Each member is drawn by
///   [`draw_bundle_member`]: motif lengths in HG002's proportions, most members near their
///   calling floor and some mono- and dinucleotides long.
/// - **Gaps of 1 to [`LONGEST_BUNDLE_GAP`] bases**, each drawn by [`draw_gap`] so it continues
///   neither neighbour's motif.
/// - **A hull of at most [`LONGEST_BUNDLE_HULL`] bases.** Long members are therefore rarer in
///   bundles than in segments, most of all for motif lengths four to six.
/// - **The typing's scanner reports exactly the members** ([`scanner_repeats`]): nothing joins
///   two members into one repeat, nothing stretches a member past its drawn span, and nothing
///   in a gap clears a copy floor. Two members of one motif length — whatever their motifs,
///   since the scanner compares each base with the one a motif length before it — are read as
///   one repeat when the gap is too short to break the run; this is the rule that catches it.
///
/// # Panics
///
/// If no shape keeps every rule in [`BUNDLE_SHAPE_ATTEMPTS`] redraws.
fn draw_bundle_shape(draw: &mut Draw) -> BundleShape {
    // Drawn once, outside the redraw, so the share holds although three members are refused
    // more often than two.
    let member_count = if draw.unit() < THREE_MEMBER_SHARE {
        3
    } else {
        2
    };
    for _ in 0..BUNDLE_SHAPE_ATTEMPTS {
        let members: Vec<(Vec<u8>, u32)> = (0..member_count)
            .map(|_| draw_bundle_member(draw))
            .collect();
        let gap_lengths: Vec<usize> = (1..member_count)
            .map(|_| 1 + draw.below(LONGEST_BUNDLE_GAP as u64) as usize)
            .collect();
        let members_bp: usize = members
            .iter()
            .map(|(motif, repeats)| motif.len() * *repeats as usize)
            .sum();
        if members_bp + gap_lengths.iter().sum::<usize>() > LONGEST_BUNDLE_HULL {
            continue;
        }

        let mut hull = members[0].0.repeat(members[0].1 as usize);
        let mut gaps = Vec::with_capacity(gap_lengths.len());
        for (index, len) in gap_lengths.iter().enumerate() {
            let (left, right) = (&members[index].0, &members[index + 1].0);
            let Some(gap) = draw_gap(&hull, left, right, *len, draw) else {
                break;
            };
            hull.extend_from_slice(&gap);
            hull.extend(right.repeat(members[index + 1].1 as usize));
            gaps.push(gap);
        }
        if gaps.len() < gap_lengths.len() {
            continue;
        }
        let shape = BundleShape {
            members,
            gaps,
            kind: BundleShapeKind::Separated,
        };
        if scanner_reports_exactly_the_members(&shape) {
            return shape;
        }
    }
    panic!(
        "no {member_count}-member bundle shape kept every rule in {BUNDLE_SHAPE_ATTEMPTS} draws: \
         the hull cap, the gap rules or the scanner copy can no longer be met together"
    );
}

/// One bundle member before the shape's rules are applied: its motif and reference repeat count.
///
/// - **Motif length** from [`BUNDLE_MOTIF_LENGTH_WEIGHTS`], HG002's shares; the motif from
///   [`MOTIFS`], as a segment's is.
/// - **Repeat count**: for a mono- or dinucleotide, with probability [`LONG_MEMBER_SHARE`],
///   uniform from [`LONG_MEMBER_LEAST_REPEATS`] to [`HIGHEST_REFERENCE_REPEATS`]; otherwise the
///   calling floor plus a geometric number of extra copies with mean [`BUNDLE_MEAN_EXTRA_COPIES`],
///   held at the highest. Never below the floor, so the typing counts each member as a repeat.
fn draw_bundle_member(draw: &mut Draw) -> (Vec<u8>, u32) {
    let total_weight: u64 = BUNDLE_MOTIF_LENGTH_WEIGHTS.iter().sum();
    let mut ticket = draw.below(total_weight);
    let period = 1 + BUNDLE_MOTIF_LENGTH_WEIGHTS
        .iter()
        .position(|weight| {
            if ticket < *weight {
                return true;
            }
            ticket -= weight;
            false
        })
        .expect("the ticket is below the weights' sum");
    let motifs = MOTIFS[period - 1];
    let motif = motifs[draw.below(motifs.len() as u64) as usize];
    (motif.as_bytes().to_vec(), draw_member_repeats(period, draw))
}

/// A bundle member's reference repeat count for its motif length, as [`draw_bundle_member`]
/// states it: long with probability [`LONG_MEMBER_SHARE`] for a mono- or dinucleotide, otherwise
/// the floor plus a geometric number of extra copies.
fn draw_member_repeats(period: usize, draw: &mut Draw) -> u32 {
    let floor = CALLING_FLOOR[period - 1];
    let highest = HIGHEST_REFERENCE_REPEATS[period - 1];
    if period <= 2 && draw.unit() < LONG_MEMBER_SHARE {
        let span = highest - LONG_MEMBER_LEAST_REPEATS + 1;
        LONG_MEMBER_LEAST_REPEATS + draw.below(u64::from(span)) as u32
    } else {
        floor + draw_geometric(BUNDLE_MEAN_EXTRA_COPIES, highest - floor, draw)
    }
}

/// A geometric count with this mean, held at `most`: each further one with probability
/// `mean / (1 + mean)`.
fn draw_geometric(mean: f64, most: u32, draw: &mut Draw) -> u32 {
    let another = mean / (1.0 + mean);
    let mut count = 0;
    while count < most && draw.unit() < another {
        count += 1;
    }
    count
}

/// **One STR bundle in a hard shape**: repeats that abut, or that carry an interruption — the
/// shapes about six in ten of HG002's bundles have ([`DEFAULT_HARD_BUNDLE_SHARE`]), which
/// [`draw_bundle_shape`] does not lay. Drawn as that function draws, then hardened:
///
/// - **Two or three repeats**, each by [`draw_bundle_member`], as in the separated shape.
/// - **Each neighbouring pair abuts** with probability [`ABUTTING_JUNCTION_SHARE`], a tenth of it
///   after a homopolymer ([`AFTER_HOMOPOLYMER_ABUTTING_FACTOR`]); otherwise a gap
///   is drawn between them as in the separated shape. When they abut, the second becomes a
///   homopolymer of the first's most frequent base with probability [`HOMOPOLYMER_TAIL_SHARE`] — the
///   `(TTTC)n T(m)` shape, and the first is then drawn short, three copies up to its floor, as
///   HG002's are. Two abutting repeats of one motif would be one longer repeat, and are redrawn.
/// - **Each repeat is interrupted** with probability [`INTERRUPTED_MEMBER_SHARE`]: lengthened, then
///   cut into pure pieces of its motif by one to three interruptions ([`draw_interrupted_member`]).
/// - **At least one of the two happens**, or the shape is redrawn: a hard shape is never the
///   separated one by construction.
/// - **With probability [`NO_RUN_HARD_BUNDLE_SHARE`], a bundle with no run**: every repeat is
///   interrupted, and the shape is redrawn until no stretch of pure copies in the hull reaches
///   its calling floor ([`reference_runs`] finds nothing) — a cluster of impure repeats the
///   caller calls with no stutter modelled, as 123 of HG002's 1,447 hard bundles are.
/// - **A hull of at most [`LONGEST_BUNDLE_HULL`] bases.**
/// - **The typing must see one bundle over exactly the laid hull**
///   ([`scanner_reports_one_bundle_over_the_hull`]). Which repeats it reports inside is left to it:
///   they may nest, overlap or join two pieces, as on HG002.
///
/// Each pure piece is one of the returned shape's members, so a length change moves one pure
/// stretch, and the stutter sites are found from the sequence (`stutter_sites_of`), not the members.
///
/// # Panics
///
/// If no shape keeps every rule in [`BUNDLE_SHAPE_ATTEMPTS`] redraws.
fn draw_hard_bundle_shape(draw: &mut Draw) -> BundleShape {
    let repeat_count = if draw.unit() < THREE_MEMBER_SHARE {
        3
    } else {
        2
    };
    let without_runs = draw.unit() < NO_RUN_HARD_BUNDLE_SHARE;
    'attempt: for _ in 0..BUNDLE_SHAPE_ATTEMPTS {
        let mut repeats: Vec<(Vec<u8>, u32)> = (0..repeat_count)
            .map(|_| draw_bundle_member(draw))
            .collect();
        // Junction by junction, so a homopolymer tail made at one junction lowers the chance
        // that the next one abuts, as for any homopolymer.
        let mut abutting: Vec<bool> = Vec::with_capacity(repeat_count - 1);
        for index in 0..repeat_count - 1 {
            let left_is_homopolymer = repeats[index].0.len() == 1;
            let share = if left_is_homopolymer {
                ABUTTING_JUNCTION_SHARE * AFTER_HOMOPOLYMER_ABUTTING_FACTOR
            } else {
                ABUTTING_JUNCTION_SHARE
            };
            let abuts = draw.unit() < share;
            abutting.push(abuts);
            if abuts && !left_is_homopolymer && draw.unit() < HOMOPOLYMER_TAIL_SHARE {
                let left = &repeats[index].0;
                let base = majority_base(left);
                // HG002's tails follow a short repeat — `(AAAC)4 A8`, `(TTC)4 T9` — that the
                // scanner reads by running on into the tail: three copies up to its floor.
                let floor = CALLING_FLOOR[left.len() - 1];
                let span = floor - TAIL_LEAD_LEAST_COPIES + 1;
                repeats[index].1 = TAIL_LEAD_LEAST_COPIES + draw.below(u64::from(span)) as u32;
                repeats[index + 1] = (vec![base], draw_member_repeats(1, draw));
            }
        }
        let interruptions: Vec<usize> = repeats
            .iter()
            .map(|_| {
                // A bundle with no run has every repeat interrupted: a pure repeat is a run.
                if without_runs || draw.unit() < INTERRUPTED_MEMBER_SHARE {
                    1 + draw.categorical(&INTERRUPTION_COUNT_WEIGHTS)
                } else {
                    0
                }
            })
            .collect();
        let abuts_somewhere = abutting.contains(&true);
        let interrupted_somewhere = interruptions.iter().any(|count| *count > 0);
        if !abuts_somewhere && !interrupted_somewhere {
            continue;
        }
        for (index, abuts) in abutting.iter().enumerate() {
            if *abuts && is_rotation_of(&repeats[index].0, &repeats[index + 1].0) {
                continue 'attempt;
            }
        }

        let mut members: Vec<(Vec<u8>, u32)> = Vec::new();
        let mut gaps: Vec<Vec<u8>> = Vec::new();
        let mut hull: Vec<u8> = Vec::new();
        for (index, ((motif, repeats), count)) in repeats.iter().zip(&interruptions).enumerate() {
            if index > 0 {
                if abutting[index - 1] {
                    gaps.push(Vec::new());
                } else {
                    let len = 1 + draw.below(LONGEST_BUNDLE_GAP as u64) as usize;
                    let left = &members.last().expect("a member before every gap").0;
                    let Some(gap) = draw_gap(&hull, left, motif, len, draw) else {
                        continue 'attempt;
                    };
                    hull.extend_from_slice(&gap);
                    gaps.push(gap);
                }
            }
            let (pieces, between) = draw_interrupted_member(motif, *repeats, *count, draw);
            for (piece, copies) in pieces.into_iter().enumerate() {
                if piece > 0 {
                    hull.extend_from_slice(&between[piece - 1]);
                    gaps.push(between[piece - 1].clone());
                }
                hull.extend(motif.repeat(copies as usize));
                members.push((motif.clone(), copies));
            }
            if hull.len() > LONGEST_BUNDLE_HULL {
                continue 'attempt;
            }
        }
        if without_runs && !reference_runs(&hull).is_empty() {
            continue 'attempt;
        }
        let kind = match (abuts_somewhere, interrupted_somewhere) {
            (true, true) => BundleShapeKind::AbuttingAndInterrupted,
            (true, false) => BundleShapeKind::Abutting,
            (false, true) => BundleShapeKind::Interrupted,
            (false, false) => {
                unreachable!("a shape that neither abuts nor is interrupted is redrawn")
            }
        };
        let shape = BundleShape {
            members,
            gaps,
            kind,
        };
        if scanner_reports_one_bundle_over_the_hull(&shape)
            && !scanner_reads_it_as_separated(&shape)
        {
            return shape;
        }
    }
    panic!(
        "no hard {repeat_count}-repeat bundle shape kept every rule in {BUNDLE_SHAPE_ATTEMPTS} \
         draws: the hull cap, the gap rules or the one-bundle rule can no longer be met together"
    );
}

/// One repeat of `repeats` copies of `motif`, cut by `interruptions` interruptions (none leaves it
/// whole): the whole copies of each pure piece, in order, and the bases between consecutive
/// pieces. With any interruption the repeat first gains `1 +` a geometric number of copies of mean
/// [`INTERRUPTED_MEMBER_MEAN_EXTRA_COPIES`], so its pieces stay long enough to be seen; the cut
/// points are drawn among its copies, each piece keeping at least one. Each interruption is, with
/// these chances:
///
/// - **a substituted copy** (one in two): a copy of the motif with one base changed;
/// - **a stray base** (one in four): one base, not the motif's first, inserted between two copies;
/// - **a copy that lost a base** (one in four): a copy of the motif with one base removed — for a
///   homopolymer, which would only be one copy fewer, a substituted copy instead.
///
/// The pieces keep the motif's phase: each starts on the motif's first base.
fn draw_interrupted_member(
    motif: &[u8],
    repeats: u32,
    interruptions: usize,
    draw: &mut Draw,
) -> (Vec<u32>, Vec<Vec<u8>>) {
    if interruptions == 0 {
        return (vec![repeats], Vec::new());
    }
    let extra = 1 + draw_geometric(
        INTERRUPTED_MEMBER_MEAN_EXTRA_COPIES,
        INTERRUPTED_MEMBER_MOST_EXTRA_COPIES,
        draw,
    );
    let total = (repeats + extra).max(interruptions as u32 + 1);
    let mut cuts: Vec<u32> = Vec::with_capacity(interruptions);
    while cuts.len() < interruptions {
        let cut = 1 + draw.below(u64::from(total - 1)) as u32;
        if !cuts.contains(&cut) {
            cuts.push(cut);
        }
    }
    cuts.sort_unstable();
    let mut pieces = Vec::with_capacity(interruptions + 1);
    let mut from = 0;
    for cut in cuts.iter().copied().chain(std::iter::once(total)) {
        pieces.push(cut - from);
        from = cut;
    }
    let between = (0..interruptions)
        .map(|_| {
            let kind = draw.below(4);
            let at = draw.below(motif.len() as u64) as usize;
            match kind {
                2 => vec![draw.other_base(motif[0])],
                3 if motif.len() > 1 => {
                    let mut copy = motif.to_vec();
                    copy.remove(at);
                    copy
                }
                _ => {
                    let mut copy = motif.to_vec();
                    copy[at] = draw.other_base(copy[at]);
                    copy
                }
            }
        })
        .collect();
    (pieces, between)
}

/// The base a motif holds most often, the earliest on a tie: `T` for `TTTC`, `A` for `AC`.
fn majority_base(motif: &[u8]) -> u8 {
    let count = |base: u8| motif.iter().filter(|one| **one == base).count();
    let mut best = motif[0];
    for base in motif {
        if count(*base) > count(best) {
            best = *base;
        }
    }
    best
}

/// Whether two motifs of one length are rotations of each other — one repeat read from another
/// base, so that two abutting runs of them would be a single run.
fn is_rotation_of(first: &[u8], second: &[u8]) -> bool {
    first.len() == second.len()
        && [first, first]
            .concat()
            .windows(second.len())
            .any(|w| w == second)
}

/// Whether the repeats the typing's scanner finds in this shape, with a spacer either side as it
/// will be laid, make **one bundle over exactly the hull**: at least two repeats, each within the
/// bundle radius of the reach of those before it (the typing's chaining, `segment_criteria.rs`,
/// `joins_cluster`), the first starting on the hull's first base and the furthest ending on its
/// last. The hull is what the typing calls the bundle (`str_bundle_loci.md` §5.1), so the bundle
/// laid is the bundle typed.
fn scanner_reports_one_bundle_over_the_hull(shape: &BundleShape) -> bool {
    let (probe, offset, hull_end) = scanner_probe(shape);
    let repeats = scanner_repeats(&probe);
    let Some(((first_start, first_end, _), rest)) = repeats.split_first() else {
        return false;
    };
    if rest.is_empty() || *first_start != offset {
        return false;
    }
    let bundle_radius = LONGEST_BUNDLE_GAP + 1;
    let mut reach = *first_end;
    for (start, end, _) in rest {
        if start.saturating_sub(reach) >= bundle_radius {
            return false;
        }
        reach = reach.max(*end);
    }
    reach == hull_end
}

/// Whether the typing's scanner reads this shape as the separated one — every repeat it reports
/// pure (each base equal to the one a period before it) and starting past the end of every repeat
/// before it — so that the typing could not tell it from a shape [`draw_bundle_shape`] lays. A hard
/// shape read so is redrawn: an interruption the scanner splits a repeat at, leaving two pure
/// repeats a few bases apart, is a separated bundle to the typing, and HG002's share of hard
/// bundles ([`DEFAULT_HARD_BUNDLE_SHARE`]) counts what the typing sees.
fn scanner_reads_it_as_separated(shape: &BundleShape) -> bool {
    let (probe, _, _) = scanner_probe(shape);
    let repeats = scanner_repeats(&probe);
    let pure = repeats.iter().all(|(start, end, period)| {
        (start + period..*end).all(|at| probe[at] == probe[at - period])
    });
    pure && repeats.windows(2).all(|pair| pair[1].0 > pair[0].1)
}

/// A shape's hull with the spacer either side as [`build_reference`] lays it, for the scanner
/// copy: the bases, and where the hull starts and ends in them.
fn scanner_probe(shape: &BundleShape) -> (Vec<u8>, usize, usize) {
    let mut probe = Vec::new();
    append_spacer(&mut probe, shape.first_motif());
    let offset = probe.len();
    probe.extend(shape.spell(&shape.reference_repeats()));
    let hull_end = probe.len();
    append_spacer(&mut probe, shape.last_motif());
    (probe, offset, hull_end)
}

/// `len` gap bases to follow `hull` (which ends with a whole copy of `left`) and precede a run
/// of `right`, or `None` if some position has no base that keeps every rule:
///
/// - **the first bases continue no copy of `left`**: base `i` differs from `left[i]` for each
///   `i` below `left`'s length — the spacer's rule ([`append_spacer`]) applied after the run;
/// - **the last bases continue no copy of `right` backwards**: the base `j` places before the
///   right member differs from the base a copy of `right` would have there;
/// - **no three copies of any unit** of one to six bases end on it, as in a clean flank.
///
/// The base is drawn and the four tried in turn from it, as [`append_capped`] does.
fn draw_gap(
    hull: &[u8],
    left: &[u8],
    right: &[u8],
    len: usize,
    draw: &mut Draw,
) -> Option<Vec<u8>> {
    let mut so_far = hull.to_vec();
    let mut gap = Vec::with_capacity(len);
    for index in 0..len {
        let from_end = len - index;
        let continues_left = |base: u8| index < left.len() && base == left[index];
        let continues_right =
            |base: u8| from_end <= right.len() && base == right[right.len() - from_end];
        let drawn = draw.below(4) as usize;
        let chosen = (0..FOUR_BASES.len())
            .map(|step| FOUR_BASES[(drawn + step) % FOUR_BASES.len()])
            .find(|base| {
                !continues_left(*base)
                    && !continues_right(*base)
                    && !completes_copies(&so_far, *base, &[2; 6])
            })?;
        so_far.push(chosen);
        gap.push(chosen);
    }
    Some(gap)
}

/// Whether the repeats the typing's scanner finds in this bundle, with a spacer either side as
/// it will be laid, are exactly its members, each over its drawn span.
fn scanner_reports_exactly_the_members(shape: &BundleShape) -> bool {
    let mut probe = Vec::new();
    append_spacer(&mut probe, shape.first_motif());
    let offset = probe.len();
    probe.extend(shape.spell(&shape.reference_repeats()));
    append_spacer(&mut probe, shape.last_motif());
    let expected: Vec<(usize, usize, usize)> = shape
        .member_runs(&shape.reference_repeats())
        .into_iter()
        .zip(shape.periods())
        .map(|((start, end), period)| (offset + start, offset + end, period))
        .collect();
    scanner_repeats(&probe) == expected
}

/// **The repeats the region typing would see in `bases`**, as `(start, end, period)`,
/// zero-based half-open, in start order — this file's own copy of the typing's arithmetic, so
/// the fixture does not change when the caller's code does:
///
/// 1. For each motif length `p` from one to six, every base scores
///    [`SCANNER_MATCH_REWARD`] if it equals the base `p` before it and loses
///    [`SCANNER_MISMATCH_PENALTY`] otherwise; each maximal scoring stretch is one repeat
///    (`tandem_repeat.rs`, `find_tandem_repeats`), kept when it holds two copies and its
///    first `p` bases are not themselves a repeat of a shorter unit.
/// 2. A repeat at least half covered by one of a shorter motif length that divides its own is
///    the same repeat seen at a multiple of its period, and is dropped.
/// 3. What is left must clear the calling copy floor ([`CALLING_FLOOR`]) to count
///    (`segment_criteria.rs`, `prefilter`).
fn scanner_repeats(bases: &[u8]) -> Vec<(usize, usize, usize)> {
    let mut found: Vec<(usize, usize, usize)> = Vec::new();
    for period in 1..=6_usize {
        if period >= bases.len() {
            continue;
        }
        let scores: Vec<i64> = (period..bases.len())
            .map(|index| {
                if bases[index] == bases[index - period] {
                    SCANNER_MATCH_REWARD
                } else {
                    -SCANNER_MISMATCH_PENALTY
                }
            })
            .collect();
        let mut stretches = Vec::new();
        maximal_scoring_stretches(&scores, 0, &mut stretches);
        for (first, last) in stretches {
            let (start, end) = (first, last + period + 1);
            if (end - start) / period >= 2 && is_primitive(&bases[start..start + period]) {
                found.push((start, end, period));
            }
        }
    }
    found.sort_by_key(|(start, _, period)| (*period, *start));
    let mut kept: Vec<(usize, usize, usize)> = Vec::with_capacity(found.len());
    for (start, end, period) in found {
        let len = end - start;
        let seen_shorter = kept.iter().any(|(other_start, other_end, other_period)| {
            *other_period < period
                && period.is_multiple_of(*other_period)
                && 2 * (*other_end)
                    .min(end)
                    .saturating_sub((*other_start).max(start))
                    >= len
        });
        if !seen_shorter {
            kept.push((start, end, period));
        }
    }
    kept.retain(|(start, end, period)| {
        (end - start) / period >= CALLING_FLOOR[period - 1] as usize
    });
    kept.sort_unstable();
    kept
}

/// The maximal scoring stretches of `scores`, as inclusive index pairs offset by `offset`:
/// the highest-scoring stretch, then the same on what lies left and right of it — the
/// recursive definition the Ruzzo–Tompa algorithm the scanner runs computes in one pass.
fn maximal_scoring_stretches(scores: &[i64], offset: usize, out: &mut Vec<(usize, usize)>) {
    let mut best: Option<(usize, usize, i64)> = None;
    let mut run_start = 0;
    let mut run_score = 0_i64;
    for (index, score) in scores.iter().enumerate() {
        if run_score <= 0 {
            run_start = index;
            run_score = 0;
        }
        run_score += score;
        if run_score > 0 && best.is_none_or(|(_, _, highest)| run_score > highest) {
            best = Some((run_start, index, run_score));
        }
    }
    let Some((first, last, _)) = best else {
        return;
    };
    maximal_scoring_stretches(&scores[..first], offset, out);
    out.push((offset + first, offset + last));
    maximal_scoring_stretches(&scores[last + 1..], offset + last + 1, out);
}

/// Whether `unit` is not itself a whole number of copies of a shorter unit.
fn is_primitive(unit: &[u8]) -> bool {
    (1..unit.len())
        .filter(|shorter| unit.len().is_multiple_of(*shorter))
        .all(|shorter| unit[shorter..] != unit[..unit.len() - shorter])
}

/// What [`build_reference`] lays at one place along a contig.
enum PlannedFeature {
    /// An STR segment: its motif and reference repeat count.
    Segment(Vec<u8>, u32),
    Bundle(BundleShape),
}

impl PlannedFeature {
    /// The feature's bases with its two flanks.
    fn footprint(&self) -> usize {
        let core = match self {
            Self::Segment(motif, repeats) => motif.len() * *repeats as usize,
            Self::Bundle(shape) => shape.hull_len(),
        };
        2 * FLANK_BP + core
    }
}

/// The reference as [`build_reference`] lays it: the contigs, and every segment and bundle on
/// them, none carrying alleles yet.
struct LaidReference {
    contigs: Vec<Contig>,
    segments: Vec<Segment>,
    bundles: Vec<Bundle>,
}

/// The contigs, and every segment's and bundle's place on them; neither carries alleles yet.
fn build_reference(
    settings: &Settings,
    bundle_shapes: Vec<BundleShape>,
    draw: &mut Draw,
) -> Result<LaidReference, String> {
    let mut shapes: Vec<PlannedFeature> = Vec::new();
    for period in 1..=6_usize {
        let floor = CALLING_FLOOR[period - 1];
        let span = HIGHEST_REFERENCE_REPEATS[period - 1] - floor + 1;
        let motifs = MOTIFS[period - 1];
        for index in 0..settings.segments_per_period {
            let repeats = floor + (index as u32 % span);
            let motif = motifs[draw.below(motifs.len() as u64) as usize];
            shapes.push(PlannedFeature::Segment(motif.as_bytes().to_vec(), repeats));
        }
    }
    shapes.extend(bundle_shapes.into_iter().map(PlannedFeature::Bundle));
    // Fisher–Yates, so motif lengths and repeat counts are interleaved along the genome.
    for index in (1..shapes.len()).rev() {
        let other = draw.below(index as u64 + 1) as usize;
        shapes.swap(index, other);
    }

    let segment_bp: usize = shapes.iter().map(PlannedFeature::footprint).sum();
    let gaps = shapes.len() + settings.contigs;
    let mean_gap = settings.genome_bp.saturating_sub(segment_bp) / gaps;
    if mean_gap < 2 * READ_LEN {
        return Err(format!(
            "{} segments and bundles leave {mean_gap} bases of ordinary sequence between them in \
             {} bases; raise genome_bp or lower segments_per_period or bundles",
            shapes.len(),
            settings.genome_bp
        ));
    }

    let mut contigs = Vec::with_capacity(settings.contigs);
    let mut segments = Vec::with_capacity(shapes.len());
    let mut bundles = Vec::new();
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
        for shape in &shapes[first..last] {
            gap(&mut bases, draw);
            match shape {
                PlannedFeature::Segment(motif, repeats) => {
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
                PlannedFeature::Bundle(shape) => {
                    let (first_motif, last_motif) = (shape.first_motif(), shape.last_motif());
                    append_aperiodic(&mut bases, FLANK_BP - first_motif.len(), draw);
                    append_spacer(&mut bases, first_motif);
                    let start = bases.len();
                    bases.extend(shape.spell(&shape.reference_repeats()));
                    append_spacer(&mut bases, last_motif);
                    append_aperiodic(&mut bases, FLANK_BP - last_motif.len(), draw);
                    bundles.push(Bundle {
                        contig: contig_index,
                        shape: shape.clone(),
                        start,
                        alleles: Vec::new(),
                    });
                }
            }
        }
        gap(&mut bases, draw);
        contigs.push(Contig {
            name: format!("chr{}", contig_index + 1),
            bases,
            ordinary,
        });
    }
    Ok(LaidReference {
        contigs,
        segments,
        bundles,
    })
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
    longest_periodic_stretch(sequence, from, to.min(sequence.len()), period)
}

/// The longest stretch of `sequence[from..to]` that repeats with this period, or `None` if it
/// is shorter than two units.
fn longest_periodic_stretch(
    sequence: &[u8],
    from: usize,
    to: usize,
    period: usize,
) -> Option<(usize, usize)> {
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

/// One allele at one STR bundle, written over the bundle's whole window.
#[derive(Debug, Clone)]
struct BundleAllele {
    sequence: Vec<u8>,
    /// Whole units of each member's motif the allele holds, first member first; each member's
    /// stutter rate is looked up by its own.
    repeats: Vec<u32>,
    change: AlleleChange,
    frequency: f64,
    /// Zero-based half-open, within `sequence`: where each member's run lies.
    member_runs: Vec<(usize, usize)>,
    /// Within `sequence`, per member: the longest stretch of the member's run that is pure
    /// repeat, where stutter adds or removes units of that member. `None` when it is shorter
    /// than two units.
    pure_spans: Vec<Option<(usize, usize)>>,
    /// Where a molecule copying this allele can stutter, within `sequence`, in order. In the
    /// separated shape, each member's pure span (one site per member, as before hard shapes
    /// existed); in a hard shape, the runs the caller's own rule finds in the allele's sequence
    /// ([`stutter_sites_of`]), so stutter acts on the pure stretches and not on how the shape was
    /// drawn.
    stutter_sites: Vec<StutterSite>,
}

/// Every allele the population carries at one bundle, with its frequency — the bundle's
/// counterpart of [`draw_segment_alleles`], with the same counts and shares:
///
/// - **two to six alleles**, the reference and length changes. A length change moves one
///   member by a whole number of its units (one most often, as for a segment), and with
///   probability [`SECOND_MEMBER_CHANGE_SHARE`] a second member as well; no member falls
///   below one copy.
/// - **derived alleles**, each from a donor among the pure length alleles: a substitution
///   inside a member of at least three copies, a SNP within 30 bases of the hull's edge, a
///   deletion of two to ten bases across the hull's left or right edge, cutting into the
///   first or the last member.
fn draw_bundle_alleles(bundle: &Bundle, window: &[u8], draw: &mut Draw) -> Vec<BundleAllele> {
    let shape = &bundle.shape;
    let left = &window[..WINDOW_FLANK_BP];
    let right = &window[window.len() - WINDOW_FLANK_BP..];
    let periods = shape.periods();
    let pure_allele = |repeats: Vec<u32>, change: AlleleChange, frequency: f64| {
        let mut sequence = left.to_vec();
        sequence.extend(shape.spell(&repeats));
        sequence.extend_from_slice(right);
        let member_runs = shape
            .member_runs(&repeats)
            .into_iter()
            .map(|(start, end)| (WINDOW_FLANK_BP + start, WINDOW_FLANK_BP + end))
            .collect();
        BundleAllele {
            sequence,
            repeats,
            change,
            frequency,
            member_runs,
            pure_spans: Vec::new(),
            stutter_sites: Vec::new(),
        }
    };
    let step = |draw: &mut Draw| {
        let mut step = 1;
        while step < LONGEST_LENGTH_STEP && draw.unit() < LENGTH_STEP_CONTINUE {
            step += 1;
        }
        step
    };
    let moved = |repeats: u32, draw: &mut Draw| -> Option<u32> {
        let units = step(draw);
        if draw.unit() < 0.5 {
            repeats.checked_sub(units).filter(|one| *one >= 1)
        } else {
            Some(repeats + units)
        }
    };

    let wanted = 2 + draw.below(5) as usize;
    let mut lengths: Vec<Vec<u32>> = vec![shape.reference_repeats()];
    for _ in 0..50 {
        if lengths.len() == wanted {
            break;
        }
        let mut repeats = shape.reference_repeats();
        let first = draw.below(repeats.len() as u64) as usize;
        let Some(count) = moved(repeats[first], draw) else {
            continue;
        };
        repeats[first] = count;
        if draw.unit() < SECOND_MEMBER_CHANGE_SHARE {
            let second =
                (first + 1 + draw.below(repeats.len() as u64 - 1) as usize) % repeats.len();
            if let Some(count) = moved(repeats[second], draw) {
                repeats[second] = count;
            }
        }
        if !lengths.contains(&repeats) {
            lengths.push(repeats);
        }
    }

    let nonreference = draw.beta(NONREF_FREQUENCY_BETA.0, NONREF_FREQUENCY_BETA.1);
    let weights: Vec<f64> = (1..lengths.len()).map(|_| draw.exponential()).collect();
    let total: f64 = weights.iter().sum();
    let mut alleles: Vec<BundleAllele> = lengths
        .into_iter()
        .enumerate()
        .map(|(index, repeats)| {
            if index == 0 {
                pure_allele(repeats, AlleleChange::Reference, 1.0 - nonreference)
            } else {
                let frequency = nonreference * weights[index - 1] / total;
                pure_allele(repeats, AlleleChange::Length, frequency)
            }
        })
        .collect();

    let pure_donors = |alleles: &[BundleAllele]| -> Vec<usize> {
        (0..alleles.len())
            .filter(|index| alleles[*index].change.is_pure_length())
            .collect()
    };

    if draw.unit() < INTERIOR_SUBSTITUTION_SHARE {
        let candidates: Vec<(usize, usize)> = pure_donors(&alleles)
            .into_iter()
            .flat_map(|donor| (0..periods.len()).map(move |member| (donor, member)))
            .filter(|(donor, member)| alleles[*donor].repeats[*member] >= 3)
            .collect();
        if !candidates.is_empty() {
            let (donor, member) = candidates[draw.below(candidates.len() as u64) as usize];
            let period = periods[member];
            let repeats = alleles[donor].repeats[member] as usize;
            let unit = (repeats / 3 + draw.below((repeats / 3).max(1) as u64) as usize)
                .clamp(1, repeats - 2);
            let at = alleles[donor].member_runs[member].0
                + unit * period
                + draw.below(period as u64) as usize;
            let mut derived = alleles[donor].clone();
            derived.sequence[at] = draw.other_base(derived.sequence[at]);
            derived.change = AlleleChange::InteriorSubstitution;
            add_derived_bundle_allele(&mut alleles, donor, derived, draw);
        }
    }
    if draw.unit() < FLANK_SNP_SHARE {
        let donors = pure_donors(&alleles);
        let donor = donors[draw.below(donors.len() as u64) as usize];
        let distance = 1 + draw.below(WINDOW_FLANK_BP as u64) as usize;
        let mut derived = alleles[donor].clone();
        let at = if draw.unit() < 0.5 {
            WINDOW_FLANK_BP - distance
        } else {
            derived.sequence.len() - WINDOW_FLANK_BP + distance - 1
        };
        derived.sequence[at] = draw.other_base(derived.sequence[at]);
        derived.change = AlleleChange::FlankSnp;
        add_derived_bundle_allele(&mut alleles, donor, derived, draw);
    }
    if draw.unit() < EDGE_DELETION_SHARE {
        let from_left = draw.unit() < 0.5;
        let member = if from_left { 0 } else { periods.len() - 1 };
        let donors: Vec<usize> = pure_donors(&alleles)
            .into_iter()
            .filter(|donor| alleles[*donor].repeats[member] >= 2)
            .collect();
        if !donors.is_empty() {
            let donor = donors[draw.below(donors.len() as u64) as usize];
            let period = periods[member];
            let (run_start, run_end) = alleles[donor].member_runs[member];
            let run_len = run_end - run_start;
            let deleted = 2 + draw.below(9) as usize; // two to ten bases
            let in_repeat = 1 + draw.below((deleted - 1).min(run_len - period) as u64) as usize;
            let in_flank = deleted - in_repeat;
            let mut derived = alleles[donor].clone();
            if from_left {
                derived
                    .sequence
                    .drain(run_start - in_flank..run_start + in_repeat);
                // The first run now starts where the deletion did; every later one moves left
                // by all the deleted bases.
                for (index, run) in derived.member_runs.iter_mut().enumerate() {
                    let start = if index == 0 {
                        run_start - in_flank
                    } else {
                        run.0 - deleted
                    };
                    *run = (start, run.1 - deleted);
                }
            } else {
                derived
                    .sequence
                    .drain(run_end - in_repeat..run_end + in_flank);
                derived.member_runs[member].1 -= in_repeat;
            }
            derived.repeats[member] = ((run_len - in_repeat) / period) as u32;
            derived.change = AlleleChange::EdgeDeletion;
            add_derived_bundle_allele(&mut alleles, donor, derived, draw);
        }
    }

    for allele in &mut alleles {
        allele.pure_spans = allele
            .member_runs
            .iter()
            .zip(&periods)
            .map(|((start, end), period)| {
                longest_periodic_stretch(&allele.sequence, *start, *end, *period)
            })
            .collect();
        allele.stutter_sites = if shape.kind == BundleShapeKind::Separated {
            allele
                .pure_spans
                .iter()
                .zip(&periods)
                .zip(&allele.repeats)
                .filter_map(|((span, period), repeats)| {
                    span.map(|(pure_start, pure_end)| StutterSite {
                        pure_start,
                        pure_end,
                        period: *period,
                        repeats: *repeats,
                    })
                })
                .collect()
        } else {
            stutter_sites_of(&allele.sequence)
        };
    }
    alleles
}

/// **Where a molecule stutters in a hard-shaped bundle allele**: one site per run of `sequence` —
/// the allele's window — as the caller finds a hull's reference runs ([`reference_runs`]), each
/// stuttering at its own motif length and its whole copies. A pure stretch below its calling
/// floor does not stutter, which is what the caller assumes and what keeps an interrupted repeat's
/// short pieces quiet; a piece a length change lifts to the floor does.
///
/// In the reference allele the window's flanks hold no run (no unit of one to six bases more than
/// twice in a row, and a spacer beside the hull that continues no copy), so its sites are the
/// hull's runs — what the caller finds. In a derived allele a SNP beside the hull or a deletion
/// across its edge can lengthen a run into the flank; the molecule then slips on the longer
/// stretch, as a real one would.
fn stutter_sites_of(sequence: &[u8]) -> Vec<StutterSite> {
    reference_runs(sequence)
        .into_iter()
        .map(|(start, length, period)| StutterSite {
            pure_start: start,
            pure_end: start + length,
            period,
            repeats: (length / period) as u32,
        })
        .collect()
}

/// **The runs of `bases`** as `(start, length, period)`, in order — this file's own copy of the
/// rule the caller uses for a bundle's reference runs (`region_typing/reference_runs.rs`,
/// `find_reference_runs`; `str_bundle_loci.md` §6.1), at the calling floors, so the fixture does
/// not change when the caller's code does:
///
/// 1. for each period from one to six, every maximal stretch that repeats with that period, of at
///    least one copy and one more base, whose first `period` bases are not a repeat of a shorter
///    unit;
/// 2. kept when its whole copies reach [`CALLING_FLOOR`];
/// 3. chosen longest first (then earlier, then shorter period); a stretch overlapping a chosen
///    run gives up the shared bases — cut back to whole copies where it lost its end, re-read from
///    its new first base where it lost its start — and each remaining piece is kept if it still
///    reaches the floor.
///
/// The test `the_run_rule_copy_agrees_with_the_library` fails when the two part.
fn reference_runs(bases: &[u8]) -> Vec<(usize, usize, usize)> {
    let mut stretches: Vec<(usize, usize, usize)> = Vec::new();
    for period in 1..=6_usize {
        let floor = CALLING_FLOOR[period - 1] as usize;
        let last_comparable = bases.len().saturating_sub(period);
        let mut at = 0;
        while at < last_comparable {
            if bases[at] != bases[at + period] {
                at += 1;
                continue;
            }
            let start = at;
            while at < last_comparable && bases[at] == bases[at + period] {
                at += 1;
            }
            let length = at + period - start;
            if is_primitive(&bases[start..start + period]) && length / period >= floor {
                stretches.push((start, length, period));
            }
        }
    }
    stretches.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)).then(a.2.cmp(&b.2)));
    let mut chosen: Vec<(usize, usize, usize)> = Vec::new();
    for (start, length, period) in stretches {
        let end = start + length;
        let mut covered: Vec<(usize, usize)> = chosen
            .iter()
            .filter(|(other, other_length, _)| *other < end && other + other_length > start)
            .map(|(other, other_length, _)| (*other, other + other_length))
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
        for (piece_start, piece_end) in pieces {
            let piece_length = if piece_end < end {
                (piece_end - piece_start) / period * period
            } else {
                piece_end - piece_start
            };
            if piece_length / period >= CALLING_FLOOR[period - 1] as usize {
                chosen.push((piece_start, piece_length, period));
            }
        }
    }
    chosen.sort_unstable();
    chosen
}

/// [`add_derived`] for a bundle: add `derived`, which came from `donor`, with a share of the
/// donor's frequency — unless an allele with its sequence is already there.
fn add_derived_bundle_allele(
    alleles: &mut Vec<BundleAllele>,
    donor: usize,
    mut derived: BundleAllele,
    draw: &mut Draw,
) {
    if alleles
        .iter()
        .any(|allele| allele.sequence == derived.sequence)
    {
        return;
    }
    derived.frequency = alleles[donor].frequency * (0.2 + 0.6 * draw.unit());
    alleles[donor].frequency -= derived.frequency;
    alleles.push(derived);
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
    bundles: Vec<Bundle>,
    /// `[bundle][sample]`, indices into the bundle's alleles, first copy first.
    bundle_genotypes: Vec<Vec<[u8; 2]>>,
    samples: Vec<String>,
    libraries: Vec<Library>,
    /// Per contig, every segment, bundle and variant on it in order along the contig.
    events: Vec<Vec<Event>>,
}

#[derive(Debug, Clone, Copy)]
enum Event {
    Segment(usize),
    Generic(usize),
    Bundle(usize),
}

/// The whole truth. **Everything about the bundles is drawn from a second draw**, seeded from
/// the first's seed and [`BUNDLE_SEED_SALT`]: their shapes, alleles and genotypes. The first
/// draw lays them along the contigs with the segments, and gives them their flanks, but with
/// `bundles=0` it makes exactly the draws it made before bundles existed.
fn draw_cohort(settings: &Settings, draw: &mut Draw) -> Result<Cohort, String> {
    let mut bundle_draw = Draw::new(settings.seed ^ BUNDLE_SEED_SALT);
    // At a share of 0 no draw is spent choosing, so the bundles are the ones drawn before hard
    // shapes existed.
    let bundle_shapes: Vec<BundleShape> = (0..settings.bundles)
        .map(|_| {
            if settings.hard_bundle_share > 0.0 && bundle_draw.unit() < settings.hard_bundle_share {
                draw_hard_bundle_shape(&mut bundle_draw)
            } else {
                draw_bundle_shape(&mut bundle_draw)
            }
        })
        .collect();
    let LaidReference {
        contigs,
        mut segments,
        mut bundles,
    } = build_reference(settings, bundle_shapes, draw)?;
    for segment in &mut segments {
        let window = &contigs[segment.contig].bases[segment.window_start()..segment.window_end()];
        segment.alleles = draw_segment_alleles(segment, window, draw);
    }
    for bundle in &mut bundles {
        let window = &contigs[bundle.contig].bases[bundle.window_start()..bundle.window_end()];
        bundle.alleles = draw_bundle_alleles(bundle, window, &mut bundle_draw);
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
    let bundle_genotypes = bundles
        .iter()
        .map(|bundle| {
            let frequencies: Vec<f64> = bundle.alleles.iter().map(|one| one.frequency).collect();
            (0..settings.samples)
                .map(|_| draw_pair(&frequencies, settings.inbreeding, &mut bundle_draw))
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
    for (index, bundle) in bundles.iter().enumerate() {
        events[bundle.contig].push((bundle.window_start(), Event::Bundle(index)));
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
        bundles,
        bundle_genotypes,
        samples,
        libraries,
        events,
    })
}

// ---------------------------------------------------------------------------
// The molecules and the reads
// ---------------------------------------------------------------------------

/// Where one copy of the genome can stutter: one segment's or one bundle member's pure repeat,
/// in the copy's own coordinates.
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
            Event::Bundle(index) => {
                let bundle = &cohort.bundles[index];
                let allele = &bundle.alleles[cohort.bundle_genotypes[index][sample][copy] as usize];
                bases.extend_from_slice(&reference[cursor..bundle.window_start()]);
                let offset = bases.len();
                bases.extend_from_slice(&allele.sequence);
                // Each site stutters at its own motif length and count.
                sites.extend(allele.stutter_sites.iter().map(|site| StutterSite {
                    pure_start: offset + site.pure_start,
                    pure_end: offset + site.pure_end,
                    ..*site
                }));
                cursor = bundle.window_end();
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

/// One STR segment's or bundle's truth record before it is written: the window and every allele
/// drawn there, with each sample's two copies.
struct WindowRecord<'a> {
    contig: &'a str,
    /// Zero-based, the window's first base.
    window_start: usize,
    id: &'a str,
    /// Every allele's window sequence, the reference first.
    sequences: &'a [&'a [u8]],
    /// `[sample]`, indices into `sequences`, first copy first.
    genotypes: &'a [[u8; 2]],
}

impl WindowRecord<'_> {
    /// Write the record into `line` — the reference and the alleles some sample carries,
    /// trimmed to the stretch they differ over, `info(carried)` as its INFO, then each sample's
    /// phased genotype — or write nothing and return `false` when no sample carries a
    /// non-reference allele. `carried` is the carried alleles' indices, ascending, as they are
    /// listed in ALT.
    fn write(&self, line: &mut String, info: impl FnOnce(&[u8]) -> String) -> bool {
        let mut carried: Vec<u8> = self.genotypes.iter().flatten().copied().collect();
        carried.sort_unstable();
        carried.dedup();
        carried.retain(|allele| *allele != 0);
        if carried.is_empty() {
            return false;
        }
        let mut alleles = vec![self.sequences[0].to_vec()];
        alleles.extend(
            carried
                .iter()
                .map(|one| self.sequences[*one as usize].to_vec()),
        );
        let position = trim_shared_bases(self.window_start, &mut alleles);
        write!(
            line,
            "{}\t{}\t{}\t{}\t{}\t.\tPASS\t{}\tGT",
            self.contig,
            position + 1,
            self.id,
            String::from_utf8_lossy(&alleles[0]),
            alleles[1..]
                .iter()
                .map(|one| String::from_utf8_lossy(one).into_owned())
                .collect::<Vec<_>>()
                .join(","),
            info(&carried),
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
        for pair in self.genotypes {
            write!(
                line,
                "\t{}|{}",
                index_in_record(pair[0]),
                index_in_record(pair[1])
            )
            .expect("writing to a String");
        }
        true
    }
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
    if !cohort.bundles.is_empty() {
        writeln!(
            out,
            "##INFO=<ID=BUNDLE,Number=1,Type=String,Description=\"The STR bundle whose window \
             this record spells\">"
        )?;
        writeln!(
            out,
            "##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"The bundle's members' motifs, \
             first member first\">"
        )?;
        writeln!(
            out,
            "##INFO=<ID=MEMBER_REPEATS,Number=R,Type=String,Description=\"Whole motif units in \
             each member of each allele, joined by _ in member order\">"
        )?;
    }
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
                    let sequences: Vec<&[u8]> = segment
                        .alleles
                        .iter()
                        .map(|one| one.sequence.as_slice())
                        .collect();
                    let info = |carried: &[u8]| {
                        let changes: Vec<&str> = carried
                            .iter()
                            .map(|one| segment.alleles[*one as usize].change.name())
                            .collect();
                        let repeats: Vec<String> = std::iter::once(0)
                            .chain(carried.iter().copied())
                            .map(|one| segment.alleles[one as usize].repeats.to_string())
                            .collect();
                        format!(
                            "SEGMENT={};MOTIF={};CHANGE={};REPEATS={}",
                            segment_name(index),
                            String::from_utf8_lossy(&segment.motif),
                            changes.join(","),
                            repeats.join(","),
                        )
                    };
                    let record = WindowRecord {
                        contig: &contig.name,
                        window_start: segment.window_start(),
                        id: &segment_name(index),
                        sequences: &sequences,
                        genotypes: &cohort.segment_genotypes[index],
                    };
                    if !record.write(&mut line, info) {
                        continue;
                    }
                }
                Event::Bundle(index) => {
                    let bundle = &cohort.bundles[index];
                    let sequences: Vec<&[u8]> = bundle
                        .alleles
                        .iter()
                        .map(|one| one.sequence.as_slice())
                        .collect();
                    let info = |carried: &[u8]| {
                        let changes: Vec<&str> = carried
                            .iter()
                            .map(|one| bundle.alleles[*one as usize].change.name())
                            .collect();
                        let repeats: Vec<String> = std::iter::once(0)
                            .chain(carried.iter().copied())
                            .map(|one| joined(&bundle.alleles[one as usize].repeats, "_"))
                            .collect();
                        format!(
                            "BUNDLE={};MOTIFS={};CHANGE={};MEMBER_REPEATS={}",
                            bundle_name(index),
                            bundle.motifs(),
                            changes.join(","),
                            repeats.join(","),
                        )
                    };
                    let record = WindowRecord {
                        contig: &contig.name,
                        window_start: bundle.window_start(),
                        id: &bundle_name(index),
                        sequences: &sequences,
                        genotypes: &cohort.bundle_genotypes[index],
                    };
                    if !record.write(&mut line, info) {
                        continue;
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

/// The bundles' three tables, the segments' counterparts. Per-member values are listed in member
/// order, comma-separated; coordinates are zero-based half-open.
///
/// - `bundles.bed` — `chrom start end name periods motifs reference_repeats member_starts
///   member_ends`: the hull, then each member's motif length, motif, reference repeat count and
///   span.
/// - `bundle_alleles.tsv` — one row per allele, as `segment_alleles.tsv`, with the hull for the
///   segment's span, the per-member lists for its period, motif and repeat counts, and `repeats`
///   the allele's own count in each member.
/// - `bundle_genotypes.tsv` — each sample's two alleles, as indices into `bundle_alleles.tsv`.
/// - `bundle_runs.tsv` — `bundle chrom hull_start hull_end shape runs`: the shape's kind
///   ([`BundleShapeKind::name`]) and where the reference stutters, each site as
///   `motif:offset:length` from the hull's first base, joined by `,` (`.` for none) — the
///   format `ng_str_bundle_runs` writes the caller's reference runs in, so the typing check can
///   compare the two. A file of its own, so the three above are unchanged by it.
fn write_bundle_tables(out_dir: &Path, cohort: &Cohort) -> std::io::Result<()> {
    let mut bed = create(&out_dir.join("bundles.bed"))?;
    let mut alleles = create(&out_dir.join("bundle_alleles.tsv"))?;
    let mut genotypes = create(&out_dir.join("bundle_genotypes.tsv"))?;
    let mut runs = create(&out_dir.join("bundle_runs.tsv"))?;
    writeln!(runs, "bundle\tchrom\thull_start\thull_end\tshape\truns")?;
    writeln!(
        bed,
        "#chrom\tstart\tend\tname\tperiods\tmotifs\treference_repeats\tmember_starts\tmember_ends"
    )?;
    writeln!(
        alleles,
        "bundle\tchrom\twindow_start\twindow_end\thull_start\thull_end\tmembers\tperiods\tmotifs\t\
         reference_repeats\tallele\trepeats\tchange\tfrequency\tsequence"
    )?;
    writeln!(genotypes, "bundle\tsample\tfirst_copy\tsecond_copy")?;
    for (index, bundle) in cohort.bundles.iter().enumerate() {
        let contig = &cohort.contigs[bundle.contig].name;
        let name = bundle_name(index);
        let periods = joined(bundle.shape.periods(), ",");
        let motifs = bundle.motifs();
        let reference_repeats = joined(bundle.shape.reference_repeats(), ",");
        let spans = bundle.member_spans();
        writeln!(
            bed,
            "{contig}\t{}\t{}\t{name}\t{periods}\t{motifs}\t{reference_repeats}\t{}\t{}",
            bundle.start,
            bundle.end(),
            joined(spans.iter().map(|(start, _)| start), ","),
            joined(spans.iter().map(|(_, end)| end), ","),
        )?;
        for (allele_index, allele) in bundle.alleles.iter().enumerate() {
            writeln!(
                alleles,
                "{name}\t{contig}\t{}\t{}\t{}\t{}\t{}\t{periods}\t{motifs}\t{reference_repeats}\t\
                 {allele_index}\t{}\t{}\t{:.6}\t{}",
                bundle.window_start(),
                bundle.window_end(),
                bundle.start,
                bundle.end(),
                spans.len(),
                joined(&allele.repeats, ","),
                allele.change.name(),
                allele.frequency,
                String::from_utf8_lossy(&allele.sequence),
            )?;
        }
        for (sample, pair) in cohort.samples.iter().zip(&cohort.bundle_genotypes[index]) {
            writeln!(genotypes, "{name}\t{sample}\t{}\t{}", pair[0], pair[1])?;
        }
        let reference = &bundle.alleles[0];
        let sites: Vec<String> = reference
            .stutter_sites
            .iter()
            .map(|site| {
                format!(
                    "{}:{}:{}",
                    String::from_utf8_lossy(
                        &reference.sequence[site.pure_start..site.pure_start + site.period]
                    ),
                    site.pure_start
                        .checked_sub(WINDOW_FLANK_BP)
                        .expect("a reference stutter site lies inside the hull"),
                    site.pure_end - site.pure_start
                )
            })
            .collect();
        writeln!(
            runs,
            "{name}\t{contig}\t{}\t{}\t{}\t{}",
            bundle.start,
            bundle.end(),
            bundle.shape.kind.name(),
            if sites.is_empty() {
                ".".to_string()
            } else {
                sites.join(",")
            }
        )?;
    }
    bed.flush()?;
    alleles.flush()?;
    genotypes.flush()?;
    runs.flush()
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
    write_bundle_tables(out_dir, cohort)?;

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
    let three_member_bundles = cohort
        .bundles
        .iter()
        .filter(|bundle| bundle.shape.members.len() == 3)
        .count();
    let hard_bundles = cohort
        .bundles
        .iter()
        .filter(|bundle| bundle.shape.kind != BundleShapeKind::Separated)
        .count();
    let varying_bundles = cohort
        .bundle_genotypes
        .iter()
        .filter(|genotypes| genotypes.iter().flatten().any(|one| *one != 0))
        .count();
    println!(
        "STR bundles         : {} ({three_member_bundles} with three members, {hard_bundles} in a \
         hard shape), {varying_bundles} varying in some sample",
        cohort.bundles.len()
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

    fn small_cohort(samples: usize, bundles: usize) -> Cohort {
        let settings = Settings {
            samples,
            segments_per_period: 40,
            bundles,
            // Hard shapes on, as a default run lays them.
            hard_bundle_share: DEFAULT_HARD_BUNDLE_SHARE,
            genome_bp: 300_000,
            contigs: 2,
            inbreeding: 0.3,
            ..Settings::default()
        };
        draw_cohort(&settings, &mut Draw::new(settings.seed)).expect("the settings fit")
    }

    /// FNV-1a, 64 bits, continuing from `hash`.
    fn fnv(bytes: &[u8], hash: u64) -> u64 {
        bytes.iter().fold(hash, |hash, byte| {
            (hash ^ u64::from(*byte)).wrapping_mul(0x0100_0000_01b3)
        })
    }

    /// **With no bundles, the cohort is the one drawn before bundles existed**: the truth VCF and
    /// the reference's names and bases hash to what the build before bundles wrote for these
    /// settings (`samples=6 segments_per_period=40 genome_bp=300000 contigs=2 inbreeding=0.3`,
    /// the default seed; 33,623 bytes of VCF). The bundles' draws come from their own stream, so
    /// nothing about the segments may move when there are none.
    #[test]
    fn no_bundles_draws_the_cohort_drawn_before_bundles() {
        let cohort = small_cohort(6, 0);
        let mut vcf = Vec::new();
        write_truth_vcf(&mut vcf, &cohort).expect("writing to memory");
        assert_eq!(vcf.len(), 33_623);
        let hash = cohort
            .contigs
            .iter()
            .fold(fnv(&vcf, 0xcbf2_9ce4_8422_2325), |hash, contig| {
                fnv(&contig.bases, fnv(contig.name.as_bytes(), hash))
            });
        assert_eq!(hash, 0x899e_4139_8adc_6a9f);
    }

    /// Every bundle shape keeps its rules: two or three members, each at or above its calling
    /// floor; gaps of 1 to 14 bases, none continuing either neighbour's motif; a hull of at most
    /// 100 bases; and the scanner's arithmetic sees exactly the members.
    #[test]
    fn bundle_shapes_keep_their_rules() {
        let mut draw = Draw::new(17);
        let mut three_members = 0;
        let mut one_motif_length_pairs = 0;
        for _ in 0..2_000 {
            let shape = draw_bundle_shape(&mut draw);
            assert!((2..=3).contains(&shape.members.len()));
            three_members += usize::from(shape.members.len() == 3);
            assert_eq!(shape.gaps.len(), shape.members.len() - 1);
            assert!(shape.hull_len() <= LONGEST_BUNDLE_HULL, "{shape:?}");
            for (motif, repeats) in &shape.members {
                assert!(*repeats >= CALLING_FLOOR[motif.len() - 1], "{shape:?}");
            }
            for (index, gap) in shape.gaps.iter().enumerate() {
                assert!((1..=LONGEST_BUNDLE_GAP).contains(&gap.len()), "{shape:?}");
                let (left, right) = (&shape.members[index].0, &shape.members[index + 1].0);
                one_motif_length_pairs += usize::from(left.len() == right.len());
                for (offset, base) in gap.iter().enumerate().take(left.len()) {
                    assert_ne!(
                        *base, left[offset],
                        "gap continues the left motif: {shape:?}"
                    );
                }
                for (from_end, base) in gap.iter().rev().enumerate().take(right.len()) {
                    assert_ne!(
                        *base,
                        right[right.len() - 1 - from_end],
                        "gap continues the right motif backwards: {shape:?}"
                    );
                }
            }
            assert!(scanner_reports_exactly_the_members(&shape));
        }
        assert!(
            within_sampling_error(three_members, 2_000, THREE_MEMBER_SHARE),
            "{three_members} of 2000 with three members"
        );
        assert!(one_motif_length_pairs > 100, "{one_motif_length_pairs}");
    }

    /// A bundle member's draw is the one stated: motif lengths in HG002's shares, the long share of
    /// mono- and dinucleotides between 16 and 30 copies, and every other member at or above its
    /// floor with [`BUNDLE_MEAN_EXTRA_COPIES`] extra copies on average.
    #[test]
    fn bundle_members_follow_the_stated_draw() {
        let mut draw = Draw::new(23);
        let trials = 20_000;
        let mut by_motif_length = [0_usize; 6];
        let (mut short_motif_members, mut long_members) = (0, 0);
        let (mut near_floor_members, mut extra_copies) = (0, 0);
        for _ in 0..trials {
            let (motif, repeats) = draw_bundle_member(&mut draw);
            let period = motif.len();
            by_motif_length[period - 1] += 1;
            let floor = CALLING_FLOOR[period - 1];
            assert!((floor..=HIGHEST_REFERENCE_REPEATS[period - 1]).contains(&repeats));
            if period <= 2 {
                short_motif_members += 1;
            }
            // Below 16 copies a mono- or dinucleotide can only be a near-floor draw: the floors
            // are 8 and 6 and a geometric count of mean 0.3 reaches 16 too rarely to matter.
            if period <= 2 && repeats >= LONG_MEMBER_LEAST_REPEATS {
                long_members += 1;
            } else {
                near_floor_members += 1;
                extra_copies += (repeats - floor) as usize;
            }
        }
        let total_weight: u64 = BUNDLE_MOTIF_LENGTH_WEIGHTS.iter().sum();
        for (count, weight) in by_motif_length.iter().zip(BUNDLE_MOTIF_LENGTH_WEIGHTS) {
            assert!(
                within_sampling_error(*count, trials, weight as f64 / total_weight as f64),
                "{by_motif_length:?}"
            );
        }
        assert!(
            within_sampling_error(long_members, short_motif_members, LONG_MEMBER_SHARE),
            "{long_members} long of {short_motif_members}"
        );
        let mean_extra = extra_copies as f64 / near_floor_members as f64;
        assert!(
            (mean_extra - BUNDLE_MEAN_EXTRA_COPIES).abs() < 0.03,
            "mean extra copies {mean_extra}"
        );
    }

    /// The scanner copy reads two runs of one motif length a base apart as one repeat — the case
    /// the shape rule exists to refuse — two runs of different lengths as two, and stretches a run
    /// over a gap base that continues its motif.
    #[test]
    fn the_scanner_copy_joins_close_runs_of_one_motif_length() {
        let mut joined = b"AC".repeat(10);
        joined.push(b'G');
        joined.extend(b"AC".repeat(10));
        assert_eq!(scanner_repeats(&joined), vec![(0, 41, 2)]);

        let mut apart = b"AC".repeat(10);
        apart.extend_from_slice(b"GC");
        apart.extend(b"AAT".repeat(6));
        assert_eq!(scanner_repeats(&apart), vec![(0, 20, 2), (22, 40, 3)]);
        // A gap base that continues a motif stretches its run: `T` before `AAT` is one more
        // base of it.
        let mut stretched = b"AC".repeat(10);
        stretched.extend_from_slice(b"GT");
        stretched.extend(b"AAT".repeat(6));
        assert_eq!(scanner_repeats(&stretched), vec![(0, 20, 2), (21, 40, 3)]);

        assert!(is_primitive(b"AAAAAC"));
        assert!(!is_primitive(b"ATAT"));
    }

    /// A bundle's alleles: distinct, frequencies summing to one, the reference the window itself,
    /// each pure allele's member runs exactly its counts of its motifs, and every kind of derived
    /// allele turning up somewhere.
    #[test]
    fn bundle_alleles_spell_their_member_counts() {
        let cohort = small_cohort(1, 120);
        let mut kinds = std::collections::HashSet::new();
        for bundle in &cohort.bundles {
            let total: f64 = bundle.alleles.iter().map(|one| one.frequency).sum();
            assert!((total - 1.0).abs() < 1e-9);
            let window =
                &cohort.contigs[bundle.contig].bases[bundle.window_start()..bundle.window_end()];
            assert_eq!(bundle.alleles[0].sequence, window);
            assert_eq!(bundle.alleles[0].repeats, bundle.shape.reference_repeats());
            for (index, allele) in bundle.alleles.iter().enumerate() {
                kinds.insert(allele.change);
                assert!(
                    bundle.alleles[..index]
                        .iter()
                        .all(|earlier| earlier.sequence != allele.sequence)
                );
                assert!(allele.repeats.iter().all(|count| *count >= 1));
                if allele.change.is_pure_length() {
                    for ((run, (motif, _)), count) in allele
                        .member_runs
                        .iter()
                        .zip(&bundle.shape.members)
                        .zip(&allele.repeats)
                    {
                        assert_eq!(
                            allele.sequence[run.0..run.1],
                            motif.repeat(*count as usize)[..]
                        );
                    }
                }
            }
        }
        assert_eq!(kinds.len(), 5, "{kinds:?}");
    }

    /// **A separated bundle stutters in every member, a hard one in its runs**: each copy of a
    /// separated bundle carries one stutter site per member with a pure run, at that member's motif
    /// length and its count in the allele the copy carries; a hard-shaped one carries one per run
    /// of the allele's sequence ([`stutter_sites_of`]).
    #[test]
    fn a_bundle_stutters_in_its_members_or_its_runs() {
        let cohort = small_cohort(3, 60);
        for sample in 0..3 {
            for copy in 0..2 {
                for contig in 0..cohort.contigs.len() {
                    let haplotype = build_haplotype(&cohort, contig, sample, copy);
                    let mut sites = haplotype.sites.iter();
                    for event in &cohort.events[contig] {
                        let expected: Vec<(usize, u32)> = match *event {
                            Event::Segment(index) => {
                                let segment = &cohort.segments[index];
                                let allele = &segment.alleles
                                    [cohort.segment_genotypes[index][sample][copy] as usize];
                                allele
                                    .pure_span
                                    .map(|_| (segment.period(), allele.repeats))
                                    .into_iter()
                                    .collect()
                            }
                            Event::Bundle(index) => {
                                let bundle = &cohort.bundles[index];
                                let allele = &bundle.alleles
                                    [cohort.bundle_genotypes[index][sample][copy] as usize];
                                assert_eq!(allele.pure_spans.len(), bundle.shape.members.len());
                                if bundle.shape.kind == BundleShapeKind::Separated {
                                    allele
                                        .pure_spans
                                        .iter()
                                        .zip(bundle.shape.periods())
                                        .zip(&allele.repeats)
                                        .filter(|((span, _), _)| span.is_some())
                                        .map(|((_, period), repeats)| (period, *repeats))
                                        .collect()
                                } else {
                                    stutter_sites_of(&allele.sequence)
                                        .iter()
                                        .map(|site| (site.period, site.repeats))
                                        .collect()
                                }
                            }
                            Event::Generic(_) => Vec::new(),
                        };
                        for (period, repeats) in expected {
                            let site = sites.next().expect("a site for every member");
                            assert_eq!((site.period, site.repeats), (period, repeats));
                        }
                    }
                    assert!(sites.next().is_none());
                }
            }
        }
    }

    /// **The truth VCF reproduces the truth sequences**: applying its phased records to the
    /// reference gives each sample's two copies of every contig exactly, the copies the reads
    /// were cut from — and so every segment's and every bundle's window as the allele tables spell it.
    #[test]
    fn truth_vcf_reproduces_the_truth_sequences() {
        let cohort = small_cohort(6, 40);
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
        let mut bundle_records = 0;
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
                        segment_records += usize::from(fields[2].starts_with("seg"));
                        bundle_records += usize::from(fields[2].starts_with("bun"));
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
        assert!(bundle_records > 0);
    }

    /// Every allele drawn is distinct, the frequencies sum to one, each kind of derived allele
    /// turns up, and an allele's pure stretch is its repeat.
    #[test]
    fn alleles_are_distinct_and_their_frequencies_sum_to_one() {
        let cohort = small_cohort(1, 0);
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
        assert_eq!(Settings::default().bundles, DEFAULT_BUNDLES);
        assert_eq!(
            Settings::parse(&["bundles=0".to_string()]).unwrap().bundles,
            0
        );
    }

    /// **Every bundle allele's member runs lie on its members' bases, derived alleles included**:
    /// between consecutive runs sit exactly the reference gap bases, so a length change in an
    /// earlier member moves the later runs with it; an untouched member's run is its motif `count`
    /// times; the member a deletion across the edge cut is the tail (left edge) or the head (right
    /// edge) of a run of its motif holding `count` whole copies; a substitution changes one base of
    /// one run. And every member's stretch where stutter acts lies inside its run and repeats with
    /// the member's motif length.
    #[test]
    fn every_bundle_allele_keeps_its_member_runs_on_its_members() {
        let cohort = small_cohort(1, 150);
        let mut edge_deletions = [0_usize; 2];
        for bundle in &cohort.bundles {
            let shape = &bundle.shape;
            let last = shape.members.len() - 1;
            for allele in &bundle.alleles {
                let runs = &allele.member_runs;
                for (index, gap) in shape.gaps.iter().enumerate() {
                    assert_eq!(&allele.sequence[runs[index].1..runs[index + 1].0], &gap[..]);
                }
                let mut changed_bases = 0;
                for (member, ((motif, _), count)) in
                    shape.members.iter().zip(&allele.repeats).enumerate()
                {
                    let period = motif.len();
                    let run = &allele.sequence[runs[member].0..runs[member].1];
                    let whole = motif.repeat(*count as usize);
                    let cut_left = allele.change == AlleleChange::EdgeDeletion
                        && member == 0
                        && runs[0].0 < WINDOW_FLANK_BP;
                    let cut_right = allele.change == AlleleChange::EdgeDeletion
                        && member == last
                        && allele.sequence.len() - runs[last].1 < WINDOW_FLANK_BP;
                    if cut_left || cut_right {
                        edge_deletions[usize::from(cut_right)] += 1;
                        assert_eq!(run.len() / period, *count as usize);
                        assert!((period..run.len()).all(|at| run[at] == run[at - period]));
                        if cut_left {
                            assert_eq!(&run[run.len() - period..], &motif[..]);
                        } else {
                            assert_eq!(&run[..period], &motif[..]);
                        }
                    } else {
                        assert_eq!(run.len(), whole.len(), "{:?}", allele.change);
                        changed_bases += run.iter().zip(&whole).filter(|(a, b)| a != b).count();
                    }
                    if let Some((start, end)) = allele.pure_spans[member] {
                        assert!(runs[member].0 <= start && end <= runs[member].1);
                        assert!(end - start >= 2 * period);
                        let stretch = &allele.sequence[start..end];
                        assert!(
                            (period..stretch.len()).all(|at| stretch[at] == stretch[at - period])
                        );
                    }
                }
                let expected = usize::from(allele.change == AlleleChange::InteriorSubstitution);
                assert_eq!(changed_bases, expected, "{:?}", allele.change);
            }
        }
        assert!(
            edge_deletions[0] > 0 && edge_deletions[1] > 0,
            "{edge_deletions:?}"
        );
    }

    /// **The scanner copy finds what the library's scanner finds**, after the calling copy floor,
    /// on bundle-like sequences whose gaps follow no rule — so the copy cannot drift from
    /// `tandem_repeat::find_tandem_repeats` without this test saying so. The copy stays the
    /// fixture's own; this only reports when the two part.
    #[test]
    fn the_scanner_copy_agrees_with_the_library_scanner() {
        use pop_var_caller::tandem_repeat::{PeriodRange, ScanParams, find_tandem_repeats};
        let periods = PeriodRange::new(1, 6).expect("1..=6 is a valid range");
        let mut draw = Draw::new(29);
        let mut differing = 0;
        for _ in 0..3_000 {
            let mut bases = Vec::new();
            append_aperiodic(&mut bases, 20, &mut draw);
            for _ in 0..2 + draw.below(2) {
                let period = 1 + draw.below(6) as usize;
                let motifs = MOTIFS[period - 1];
                let motif = motifs[draw.below(motifs.len() as u64) as usize].as_bytes();
                bases.extend(motif.repeat(2 + draw.below(14) as usize));
                for _ in 0..draw.below(8) {
                    bases.push(FOUR_BASES[draw.below(4) as usize]);
                }
            }
            append_aperiodic(&mut bases, 20, &mut draw);
            let library: Vec<(usize, usize, usize)> =
                find_tandem_repeats(&bases, periods, &ScanParams::default())
                    .into_iter()
                    .map(|one| {
                        (
                            one.start as usize,
                            one.end as usize,
                            usize::from(one.period),
                        )
                    })
                    .filter(|(start, end, period)| {
                        (end - start) / period >= CALLING_FLOOR[period - 1] as usize
                    })
                    .collect();
            let mut library = library;
            library.sort_unstable();
            differing += usize::from(scanner_repeats(&bases) != library);
        }
        assert_eq!(
            differing, 0,
            "{differing} of 3000 sequences scanned differently"
        );
    }

    /// **A zero hard share lays the bundles laid before hard shapes existed**: the truth VCF and the
    /// reference hash to a recorded value for these settings (`samples=6 segments_per_period=40
    /// bundles=40 genome_bp=300000 contigs=2 inbreeding=0.3`, the default seed). The value was
    /// taken from this build at share 0, after every output file of a full run on these settings
    /// was compared byte for byte with the build before hard shapes and found equal
    /// (`str_stutter_benchmark.md` §11.4); the comparison is the evidence, the hash guards it.
    #[test]
    fn a_zero_hard_share_draws_the_bundles_drawn_before_hard_shapes() {
        let settings = Settings {
            samples: 6,
            segments_per_period: 40,
            bundles: 40,
            hard_bundle_share: 0.0,
            genome_bp: 300_000,
            contigs: 2,
            inbreeding: 0.3,
            ..Settings::default()
        };
        let cohort =
            draw_cohort(&settings, &mut Draw::new(settings.seed)).expect("the settings fit");
        assert!(
            cohort
                .bundles
                .iter()
                .all(|bundle| bundle.shape.kind == BundleShapeKind::Separated)
        );
        let mut vcf = Vec::new();
        write_truth_vcf(&mut vcf, &cohort).expect("writing to memory");
        let hash = cohort
            .contigs
            .iter()
            .fold(fnv(&vcf, 0xcbf2_9ce4_8422_2325), |hash, contig| {
                fnv(&contig.bases, fnv(contig.name.as_bytes(), hash))
            });
        assert_eq!(hash, ZERO_HARD_SHARE_HASH);
    }

    /// The hash [`a_zero_hard_share_draws_the_bundles_drawn_before_hard_shapes`] expects.
    const ZERO_HARD_SHARE_HASH: u64 = 0xe542_83ee_1512_eb5b;

    /// Every hard shape keeps its rules: not the separated kind; every member at least one copy; a
    /// hull of at most 100 bases; two abutting members never of one motif; at least one abutting
    /// junction or interruption, as its kind says; and the scanner's arithmetic — which agrees with
    /// the library's scanner on these hulls too — sees one bundle over exactly the hull. And the
    /// kinds all turn up, with members nested or overlapping in a good share of them, as on HG002.
    #[test]
    fn hard_shapes_keep_their_rules() {
        use pop_var_caller::tandem_repeat::{PeriodRange, ScanParams, find_tandem_repeats};
        let periods = PeriodRange::new(1, 6).expect("1..=6 is a valid range");
        let mut draw = Draw::new(31);
        let mut kinds = std::collections::HashMap::new();
        let mut nested_or_overlapping = 0;
        let mut without_runs = 0;
        let shapes = 2_000;
        for _ in 0..shapes {
            let shape = draw_hard_bundle_shape(&mut draw);
            *kinds.entry(shape.kind).or_insert(0) += 1;
            assert_ne!(shape.kind, BundleShapeKind::Separated);
            assert!(shape.hull_len() <= LONGEST_BUNDLE_HULL, "{shape:?}");
            assert_eq!(shape.gaps.len(), shape.members.len() - 1);
            assert!(shape.members.iter().all(|(_, copies)| *copies >= 1));
            let abutting = shape.gaps.iter().filter(|gap| gap.is_empty()).count();
            for (index, gap) in shape.gaps.iter().enumerate() {
                if gap.is_empty() {
                    assert!(!is_rotation_of(
                        &shape.members[index].0,
                        &shape.members[index + 1].0
                    ));
                }
            }
            let interrupted = shape
                .gaps
                .iter()
                .enumerate()
                .filter(|(index, gap)| {
                    !gap.is_empty() && shape.members[*index].0 == shape.members[index + 1].0
                })
                .count();
            match shape.kind {
                // Two repeats of one motif a gap apart look interrupted here, so an abutting
                // shape is only asked to abut.
                BundleShapeKind::Abutting => assert!(abutting > 0, "{shape:?}"),
                BundleShapeKind::Interrupted => {
                    assert!(abutting == 0 && interrupted > 0, "{shape:?}")
                }
                BundleShapeKind::AbuttingAndInterrupted => {
                    assert!(abutting > 0 && interrupted > 0, "{shape:?}")
                }
                BundleShapeKind::Separated => unreachable!("asserted above"),
            }
            without_runs +=
                usize::from(reference_runs(&shape.spell(&shape.reference_repeats())).is_empty());
            assert!(
                scanner_reports_one_bundle_over_the_hull(&shape),
                "{shape:?}"
            );

            let mut probe = Vec::new();
            append_spacer(&mut probe, shape.first_motif());
            probe.extend(shape.spell(&shape.reference_repeats()));
            append_spacer(&mut probe, shape.last_motif());
            let mut library: Vec<(usize, usize, usize)> =
                find_tandem_repeats(&probe, periods, &ScanParams::default())
                    .into_iter()
                    .map(|one| {
                        (
                            one.start as usize,
                            one.end as usize,
                            usize::from(one.period),
                        )
                    })
                    .filter(|(start, end, period)| {
                        (end - start) / period >= CALLING_FLOOR[period - 1] as usize
                    })
                    .collect();
            library.sort_unstable();
            let found = scanner_repeats(&probe);
            assert_eq!(found, library, "{shape:?}");
            nested_or_overlapping +=
                usize::from(found.windows(2).any(|pair| pair[1].0 < pair[0].1));
        }
        assert_eq!(kinds.len(), 3, "{kinds:?}");
        assert!(
            nested_or_overlapping > shapes / 5,
            "{nested_or_overlapping}"
        );
        // NO_RUN_HARD_BUNDLE_SHARE is about 1 in 12: 2,000 draws give about 170.
        assert!(
            (shapes / 20..shapes / 8).contains(&without_runs),
            "{without_runs} of {shapes} hard shapes have no run"
        );
    }

    /// **The run rule's copy finds what the caller finds**: on HG002's two examples
    /// (`str_bundle_loci.md` §5.1), on hard hulls, on every allele window of a cohort's bundles
    /// (flanks, length changes, substitutions and edge deletions included — what
    /// [`stutter_sites_of`] is called on) and on bundle-like sequences whose gaps follow no rule,
    /// [`reference_runs`] equals `find_reference_runs` at the default floors.
    #[test]
    fn the_run_rule_copy_agrees_with_the_library() {
        use pop_var_caller::region_typing::reference_runs::find_reference_runs;
        use pop_var_caller::region_typing::segment_criteria::MinCopies;
        let library = |bases: &[u8]| -> Vec<(usize, usize, usize)> {
            find_reference_runs(bases, &MinCopies::default())
                .into_iter()
                .map(|run| (run.offset() as usize, run.length() as usize, run.period()))
                .collect()
        };
        let ct_at_t = [b"CT".repeat(20), b"AT".repeat(8), b"T".repeat(8)].concat();
        assert_eq!(
            reference_runs(&ct_at_t),
            vec![(0, 40, 2), (40, 16, 2), (56, 8, 1)]
        );
        assert_eq!(
            reference_runs(b"AAAAAAAAATAAAAATAAAATAAAAATAAA"),
            vec![(0, 9, 1)]
        );
        let mut draw = Draw::new(37);
        let mut differing = 0;
        for _ in 0..1_000 {
            let shape = draw_hard_bundle_shape(&mut draw);
            let hull = shape.spell(&shape.reference_repeats());
            differing += usize::from(reference_runs(&hull) != library(&hull));
        }
        for _ in 0..2_000 {
            let mut bases = Vec::new();
            for _ in 0..2 + draw.below(3) {
                let period = 1 + draw.below(6) as usize;
                let motifs = MOTIFS[period - 1];
                let motif = motifs[draw.below(motifs.len() as u64) as usize].as_bytes();
                bases.extend(motif.repeat(2 + draw.below(14) as usize));
                for _ in 0..draw.below(4) {
                    bases.push(FOUR_BASES[draw.below(4) as usize]);
                }
            }
            differing += usize::from(reference_runs(&bases) != library(&bases));
        }
        let cohort = small_cohort(4, 150);
        let mut windows = 0;
        for allele in cohort.bundles.iter().flat_map(|bundle| &bundle.alleles) {
            windows += 1;
            differing += usize::from(reference_runs(&allele.sequence) != library(&allele.sequence));
        }
        assert!(windows > 300, "{windows} allele windows");
        assert_eq!(differing, 0, "{differing} sequences found different runs");
    }

    /// **A hard bundle's truth names its stutter sites**: in every allele of every hard-shaped
    /// bundle, each site is pure whole copies of its period at or above its calling floor and the
    /// sites do not overlap; the reference's are exactly the caller's runs of the hull
    /// (`find_reference_runs`), moved by the window's flank; and some bundles have none.
    #[test]
    fn a_hard_bundle_stutters_in_pure_stretches_at_the_floor() {
        use pop_var_caller::region_typing::reference_runs::find_reference_runs;
        use pop_var_caller::region_typing::segment_criteria::MinCopies;
        let cohort = small_cohort(1, 150);
        let mut hard = 0;
        let mut sites = 0;
        let mut without_sites = 0;
        for bundle in &cohort.bundles {
            if bundle.shape.kind == BundleShapeKind::Separated {
                continue;
            }
            hard += 1;
            let reference = &bundle.alleles[0];
            let hull =
                &reference.sequence[WINDOW_FLANK_BP..WINDOW_FLANK_BP + bundle.shape.hull_len()];
            let callers: Vec<(usize, usize, usize)> =
                find_reference_runs(hull, &MinCopies::default())
                    .into_iter()
                    .map(|run| {
                        let start = WINDOW_FLANK_BP + run.offset() as usize;
                        (start, start + run.length() as usize, run.period())
                    })
                    .collect();
            let laid: Vec<(usize, usize, usize)> = reference
                .stutter_sites
                .iter()
                .map(|site| (site.pure_start, site.pure_end, site.period))
                .collect();
            assert_eq!(laid, callers, "{:?}", bundle.shape);
            sites += laid.len();
            without_sites += usize::from(laid.is_empty());
            for (index, allele) in bundle.alleles.iter().enumerate() {
                let mut previous_end = 0;
                for site in &allele.stutter_sites {
                    assert!(site.pure_start >= previous_end);
                    previous_end = site.pure_end;
                    let stretch = &allele.sequence[site.pure_start..site.pure_end];
                    assert!(
                        (site.period..stretch.len())
                            .all(|at| stretch[at] == stretch[at - site.period])
                    );
                    assert!(site.repeats >= CALLING_FLOOR[site.period - 1]);
                    assert_eq!(site.repeats as usize, stretch.len() / site.period);
                    if index == 0 {
                        assert!(site.pure_start >= WINDOW_FLANK_BP);
                        assert!(site.pure_end <= WINDOW_FLANK_BP + bundle.shape.hull_len());
                    }
                }
            }
        }
        assert!(hard > 50, "{hard} hard bundles");
        assert!(
            sites > hard,
            "{sites} reference sites in {hard} hard bundles"
        );
        assert!(without_sites > 0, "no hard bundle without a run");
    }

    /// The scanner copy's scores and the widest gap laid are the library's, so a retuned scanner or
    /// bundle radius in `src/` shows up here rather than as different bundles.
    #[test]
    fn the_scanner_copy_uses_the_library_scores() {
        use pop_var_caller::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD;
        use pop_var_caller::tandem_repeat::{DEFAULT_MATCH_REWARD, DEFAULT_MISMATCH_PENALTY};
        assert_eq!(SCANNER_MATCH_REWARD, i64::from(DEFAULT_MATCH_REWARD));
        assert_eq!(
            SCANNER_MISMATCH_PENALTY,
            i64::from(DEFAULT_MISMATCH_PENALTY)
        );
        assert_eq!(LONGEST_BUNDLE_GAP as u64, DEFAULT_BUNDLE_THRESHOLD - 1);
    }
}
