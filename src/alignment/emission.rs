//! Emission — scoring one read base against one reference base.
//!
//! This is a **component, not a variant**: per-base quality and a flat error rate are two
//! configurations of the same algorithm, which is what makes comparing them a swap rather
//! than a fork (spec §3). Every aligner in this module takes its emission model as a type
//! parameter, never behind `dyn` — this is called once per matrix cell, so a virtual call
//! per cell is not affordable (arch §4). The `Sized` supertrait on [`Emission`] makes that
//! decision a compile error rather than a convention.
//!
//! Two implementations land here: [`PerQualityEmission`], which reads the read's own
//! quality score, and [`FlatEmission`], which ignores it in favour of one rate. The first
//! is the default, on the reasoning that the read already carries a confidence for every
//! base and throwing it away to align the read would discard information we paid for
//! (spec §4.1); the second is the quality-blind end of the comparison.
//!
//! ## The scores must not move by a bit
//!
//! The repeat-aware aligner built on [`PerQualityEmission`]'s table breaks ties on its values,
//! so these are not merely "close enough" numbers: a table that moves the last bit of one entry
//! can move a measured repeat by a byte. The table's bits are written out rather than computed,
//! so they are the same on every platform, and `per_quality_table_matches_the_dindel_model`
//! checks them against the formula.

use crate::float;
use crate::types::{BaseQual, DomainError};

/// What a base is worth at one quality: the score if it agrees with the reference, and
/// the score if it does not. Both are natural logarithms.
///
/// This exists so quality can be resolved **once per matrix row** rather than once per
/// cell — see [`Emission::scores_for`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaseScores {
    /// Score for a read base equal to the reference base, `ln(1 − ε)`.
    pub match_ln: f64,
    /// Score for a read base differing from the reference base, `ln(ε / 3)`.
    pub mismatch_ln: f64,
}

impl BaseScores {
    /// Pick the match or mismatch score by comparing the two bases.
    ///
    /// **Comparison is raw byte equality** — see [`Emission::scores_for`] for the
    /// precondition that puts on the caller.
    #[inline]
    #[must_use]
    pub fn pick(&self, read_base: u8, reference_base: u8) -> f64 {
        if read_base == reference_base {
            self.match_ln
        } else {
            self.mismatch_ln
        }
    }
}

/// Scores one read base against one reference base, in log space.
///
/// # Contract: pure and total
///
/// An implementation is a pure function of its arguments and its own constructor state —
/// no hidden mutation, so a score never depends on call order or thread count, which the
/// cohort's byte-identity guarantee rests on.
///
/// **Total** is the sharp end: *no* argument may produce a non-finite score. In
/// particular a quality-zero base must not score `-inf`. Quality zero means an error
/// probability of 1, so the naive `ln(1 − ε)` is `ln(0)`, and a single `-inf` annihilates
/// **every path** through that cell — the whole read becomes unalignable because of one
/// worthless base. Production floors it precisely to stop that (spec §4.2), and both
/// implementations here inherit the floor. `emission_is_finite_at_every_quality` pins it.
///
/// A `+inf` is worse still, and is the reason totality is stated unconditionally rather
/// than only for well-formed input: a positive log-probability makes the event it scores
/// *infinitely preferred*, so one would not produce a slightly wrong alignment but a
/// uniformly wrong one, with nothing in the output to show for it.
///
/// # Precondition: the caller canonicalizes the bases
///
/// Bases are compared by **raw byte equality**, exactly as production does. Two
/// consequences the caller owns, because nothing here can check them without paying for
/// it on every cell:
///
/// - **Case matters.** A soft-masked reference spells its bases in lower case, and `b'a'`
///   is not `b'A'` — a soft-masked stretch would score as a mismatch at every base.
///   Soft-masking marks repeats, which is precisely where the repeat-aware aligner works,
///   so this is not a remote hazard. ng offers both shapes deliberately:
///   [`RefSeq::fetch`] canonicalizes (upper-cases ACGT, folds everything else to `N`)
///   while [`RawRefSeq`] returns bases verbatim, soft-mask intact. **Feed this the
///   canonical form.**
/// - **`N` is a base like any other here.** An `N` read base against an `N` reference
///   base scores a full-confidence *match*. Neither this component nor production's
///   treats ambiguity specially; if that should change it is a scoring-model decision for
///   the aligner steps, not a silent fix here.
///
/// Both behaviours are pinned by `emission_compares_bases_by_raw_byte_equality` so they
/// are recorded rather than discovered.
///
/// [`RefSeq::fetch`]: crate::ref_seq::RefSeq::fetch
/// [`RawRefSeq`]: crate::ref_seq::RawRefSeq
pub trait Emission: Sized {
    /// Resolve the two scores for one quality, so a caller can hoist the lookup out of
    /// its inner loop.
    ///
    /// **This is the primary method, and the shape is deliberate.** Arch §6 left open
    /// "whether quality arrives per call or as a pre-resolved row", to be settled once
    /// two implementations existed; they now do, and the answer is per row. A quality
    /// belongs to a *read base*, so it is constant along a whole row of the alignment
    /// matrix while the reference base varies along it — production's own loop resolves
    /// it once per row and leaves a compare-and-select in the inner loop
    /// ([src/ssr/pileup/alignment.rs](../../../ssr/pileup/alignment.rs)). Billing a table
    /// lookup per *cell* would put work in the hottest loop in the module that the
    /// structure of the problem does not require.
    ///
    /// No measurement backs the choice — nothing is built yet to measure. It rests on
    /// matching the loop shape of the code this module must reproduce byte for byte.
    fn scores_for(&self, quality: BaseQual) -> BaseScores;

    /// Score an **inserted** read base, which has no reference base to score against.
    ///
    /// The two implementations have **no common source** for this value, so it is worth
    /// stating where each comes from. Production's per-quality path scores an inserted
    /// base against a uniform base composition, `ln(1/4)`. Its flat path has no emission
    /// for an inserted base at all — the value it uses there is a *transition* cost
    /// (`gap = eps`, [src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)), not
    /// an emission. So [`FlatEmission`]'s value is a **decision, not a port** (arch §2.3);
    /// see its documentation for the reasoning.
    fn insert_ln(&self) -> f64;

    /// Score `read_base` against `reference_base` at one quality — the convenience form
    /// of [`Self::scores_for`], for callers not in an inner loop.
    #[inline]
    fn emit_ln(&self, read_base: u8, reference_base: u8, quality: BaseQual) -> f64 {
        self.scores_for(quality).pick(read_base, reference_base)
    }
}

/// The smallest positive value a floored log may take the logarithm of. Using
/// [`f64::MIN_POSITIVE`] rather than an arbitrary epsilon keeps the floored probability as
/// close to the true zero as the type allows, so a floored base is *merely* wildly
/// improbable (about `-708` in log space) instead of impossible.
const PROBABILITY_FLOOR: f64 = f64::MIN_POSITIVE;

/// An inserted read base has no reference base to be compared with, so it is scored
/// against a uniform base composition: any of the four bases, `1/4`. Production's
/// per-quality path uses exactly this (`INS_EMIT_LN`).
///
/// A `const` rather than a lazily-computed static: this is read on the hot path, and the
/// module's whole argument against `dyn` is that per-cell overhead is not affordable — an
/// atomic initialisation check would be the same kind of cost. `uniform_base_ln_is_ln_of_a_quarter`
/// pins the literal against `float::ln(0.25)`.
const UNIFORM_BASE_LN: f64 = -1.386_294_361_119_890_6;

/// Per-Phred-quality emission scores, indexed by the raw quality byte — all 256 of them,
/// so no quality a BAM can hold can index out of range.
///
/// The Dindel base-quality model: with `ε = 10^(−Q/10)`, a match scores `ln(1 − ε)` and a
/// mismatch `ln(ε / 3)`, the error split evenly across the three other bases, each probability
/// floored at [`PROBABILITY_FLOOR`] before its logarithm. The floor binds only at quality zero,
/// where `1 − ε` is exactly zero; the mismatch floor never binds (the smallest `ε / 3` any `u8`
/// quality gives is about `1.05e-26`, at Q255), which
/// `mismatch_floor_never_binds_over_the_quality_domain` asserts.
///
/// **Written out as bits, not computed, so every platform scores alike.** `f64::powf` and
/// `f64::ln` call the platform's maths library, and those are not required to round alike: at
/// Q4 macOS's gives a match score two units in the last place away from glibc's
/// (`0xbfe03ee1794fdd29` against `…27`). One such difference was enough for a repeat-tract
/// delimitation to break a tie the other way on macOS — the tract measured one byte longer than
/// on Linux (`alignment::delimit_parity`, seed `0x5eed0001`, case 287). These are the bits
/// glibc produced on aarch64 Linux with Rust 1.98, the values every run had used on Linux
/// before; `per_quality_table_matches_the_dindel_model` re-derives them from the formula through
/// [`crate::float`], to within rounding. They stay glibc's bits rather than libm's: already the same on
/// every platform, and changing them would move measured tracts for no gain.
static PER_QUALITY_LN: [BaseScores; 256] = [
    scores(0xc086232bdd7abcd2, 0xbff193ea7aad030b), // Q0
    scores(0xbff94db76c25f264, 0xbff5430e069e140f), // Q1
    scores(0xbfefe62362284804, 0xbff8f231928f2513), // Q2
    scores(0xbfe641bc893e144c, 0xbffca1551e803618), // Q3
    scores(0xbfe03ee1794fdd27, 0xc000283c5538a38e), // Q4
    scores(0xbfd8540e7db6ff76, 0xc001ffce1b312c10), // Q5
    scores(0xbfd2835eb6f3768d, 0xc003d75fe129b492), // Q6
    scores(0xbfcc7c916d63a729, 0xc005aef1a7223d14), // Q7
    scores(0xbfc6165572aff83d, 0xc00786836d1ac596), // Q8
    scores(0xbfc138ffa686c7d3, 0xc0095e1533134e19), // Q9
    scores(0xbfbaf8e8210a415c, 0xc00b35a6f90bd69b), // Q10
    scores(0xbfb5301b5c35244b, 0xc00d0d38bf045f1d), // Q11
    scores(0xbfb0af412e6f7610, 0xc00ee4ca84fce79f), // Q12
    scores(0xbfaa537efbd9b517, 0xc0105e2e257ab811), // Q13
    scores(0xbfa4ccc792ac9150, 0xc01149f70876fc51), // Q14
    scores(0xbfa073cfd310c471, 0xc01235bfeb734093), // Q15
    scores(0xbf9a0cdf371d5393, 0xc0132188ce6f84d4), // Q16
    scores(0xbf94a3588615a7bd, 0xc0140d51b16bc915), // Q17
    scores(0xbf905bfa6fc24d66, 0xc014f91a94680d56), // Q18
    scores(0xbf89f24b3db557eb, 0xc015e4e377645197), // Q19
    scores(0xbf8495453e6fd4bc, 0xc016d0ac5a6095d8), // Q20
    scores(0xbf8055322655bf0f, 0xc017bc753d5cda19), // Q21
    scores(0xbf79ed071c97f2ba, 0xc018a83e20591e5b), // Q22
    scores(0xbf74948af0a9857d, 0xc01994070355629b), // Q23
    scores(0xbf7056c9ab3d0327, 0xc01a7fcfe651a6dc), // Q24
    scores(0xbf69f248eb25b279, 0xc01b6b98c94deb1d), // Q25
    scores(0xbf649a6f506e2570, 0xc01c5761ac4a2f5f), // Q26
    scores(0xbf605c8c62163836, 0xc01d432a8f4673a0), // Q27
    scores(0xbf59fccc17b83728, 0xc01e2ef37242b7e0), // Q28
    scores(0xbf54a3a4773eb613, 0xc01f1abc553efc22), // Q29
    scores(0xbf5064670d979b73, 0xc02003429c1da031), // Q30
    scores(0xbf4a09f4bd4ebbc3, 0xc02079270d9bc252), // Q31
    scores(0xbf44ae86395c96fa, 0xc020ef0b7f19e473), // Q32
    scores(0xbf406d5130d1fa42, 0xc02164eff0980693), // Q33
    scores(0xbf3a1875b8ae5a61, 0xc021dad4621628b4), // Q34
    scores(0xbf34ba42b4ef63d1, 0xc02250b8d3944ad4), // Q35
    scores(0xbf3076c686f7678e, 0xc022c69d45126cf5), // Q36
    scores(0xbf2a27a84952ef96, 0xc0233c81b6908f15), // Q37
    scores(0xbf24c670c82e428c, 0xc023b266280eb136), // Q38
    scores(0xbf208084c54d40d0, 0xc024284a998cd356), // Q39
    scores(0xbf1a3738d2cf1cc2, 0xc0249e2f0b0af577), // Q40
    scores(0xbf14d2dbb7c60f02, 0xc02514137c891797), // Q41
    scores(0xbf108a6aa997a821, 0xc02589f7ee0739b8), // Q42
    scores(0xbf0a46fd609d729e, 0xc025ffdc5f855bd8), // Q43
    scores(0xbf04df690d1e7575, 0xc02675c0d1037df9), // Q44
    scores(0xbf00946782bfa4ba, 0xc026eba54281a019), // Q45
    scores(0xbefa56e0e44eec72, 0xc0276189b3ffc23a), // Q46
    scores(0xbef4ec0b805b93be, 0xc027d76e257de45b), // Q47
    scores(0xbef09e72f086842d, 0xc0284d5296fc067b), // Q48
    scores(0xbeea66d8cd676e03, 0xc028c337087a289c), // Q49
    scores(0xbee4f8bc681df714, 0xc029391b79f84abc), // Q50
    scores(0xbee0a888bfb56e51, 0xc029aeffeb766cdc), // Q51
    scores(0xbeda76dfd05491c0, 0xc02a24e45cf48efe), // Q52
    scores(0xbed505786e135f35, 0xc02a9ac8ce72b11e), // Q53
    scores(0xbed0b2a6d608dde0, 0xc02b10ad3ff0d33f), // Q54
    scores(0xbeca86f347105d68, 0xc02b8691b16ef55f), // Q55
    scores(0xbec5123de7675c5d, 0xc02bfc7622ed177f), // Q56
    scores(0xbec0bccc26f4b525, 0xc02c725a946b39a0), // Q57
    scores(0xbeba9711dff4578d, 0xc02ce83f05e95bc0), // Q58
    scores(0xbeb51f0c0005d7d2, 0xc02d5e2377677de2), // Q59
    scores(0xbeb0c6f82d74d230, 0xc02dd407e8e5a002), // Q60
    scores(0xbeaaa73af4594ad2, 0xc02e49ec5a63c222), // Q61
    scores(0xbea52be24fbea82e, 0xc02ebfd0cbe1e443), // Q62
    scores(0xbea0d12aa8840e2a, 0xc02f35b53d600663), // Q63
    scores(0xbe9ab76e3378b960, 0xc02fab99aede2884), // Q64
    scores(0xbe9538c0a48b8597, 0xc03010bf102e2552), // Q65
    scores(0xbe90db6379650e1f, 0xc0304bb148ed3662), // Q66
    scores(0xbe8ac7ab77d2a8d4, 0xc03086a381ac4773), // Q67
    scores(0xbe8545a6e7c8053e, 0xc030c195ba6b5883), // Q68
    scores(0xbe80e5a2929824dd, 0xc030fc87f32a6993), // Q69
    scores(0xbe7ad7f2b1049b9f, 0xc031377a2be97aa3), // Q70
    scores(0xbe755295103538be, 0xc031726c64a88bb4), // Q71
    scores(0xbe70efe7eef6ee83, 0xc031ad5e9d679cc4), // Q72
    scores(0xbe6ae843db50020b, 0xc031e850d626add4), // Q73
    scores(0xbe655f8b1c2341eb, 0xc03223430ee5bee5), // Q74
    scores(0xbe60fa338e80ebe9, 0xc0325e3547a4cff5), // Q75
    scores(0xbe5af89ef5aee37c, 0xc03299278063e105), // Q76
    scores(0xbe556c890b95f8ff, 0xc032d419b922f216), // Q77
    scores(0xbe5104857243339b, 0xc0330f0bf1e20325), // Q78
    scores(0xbe4b090402dae72a, 0xc03349fe2aa11436), // Q79
    scores(0xbe45798ee5cd2b2a, 0xc03384f063602546), // Q80
    scores(0xbe410edd9d22fa4c, 0xc033bfe29c1f3656), // Q81
    scores(0xbe3b1973096f3066, 0xc033fad4d4de4766), // Q82
    scores(0xbe35869ca8e7ae3e, 0xc03435c70d9d5877), // Q83
    scores(0xbe31193c10922e3b, 0xc03470b9465c6988), // Q84
    scores(0xbe2b29ec10b877aa, 0xc034abab7f1b7a98), // Q85
    scores(0xbe2593b26074641f, 0xc034e69db7da8ba8), // Q86
    scores(0xbe2123a0d0497014, 0xc035218ff0999cb8), // Q87
    scores(0xbe1b3a6f205cac19, 0xc0355c822958adc9), // Q88
    scores(0xbe15a0d0003a78e5, 0xc03597746217bed9), // Q89
    scores(0xbe112e0be024e4bc, 0xc035d2669ad6cfe9), // Q90
    scores(0xbe0b4afc402e8e73, 0xc0360d58d395e0f9), // Q91
    scores(0xbe05adf5c01d6008, 0xc036484b0c54f209), // Q92
    scores(0xbe01387d401288d2, 0xc036833d4514031a), // Q93
    scores(0xbdfb5b938017638f, 0xc036be2f7dd3142a), // Q94
    scores(0xbdf5bb23000ec1e4, 0xc036f921b692253a), // Q95
    scores(0xbdf142f500094fb0, 0xc0373413ef51364a), // Q96
    scores(0xbdeb6c35000bc004, 0xc0376f062810475a), // Q97
    scores(0xbde5c859000769ee, 0xc037a9f860cf586b), // Q98
    scores(0xbde14d730004ad83, 0xc037e4ea998e697c), // Q99
    scores(0xbddb7ce00005e728, 0xc0381fdcd24d7a8c), // Q100
    scores(0xbdd5d5960003b97a, 0xc0385acf0b0c8b9c), // Q101
    scores(0xbdd157f80002599a, 0xc03895c143cb9cac), // Q102
    scores(0xbdcb8d940002f72c, 0xc038d0b37c8aadbd), // Q103
    scores(0xbdc5e2dc0001df01, 0xc0390ba5b549becd), // Q104
    scores(0xbdc1628400012e3b, 0xc0394697ee08cfdd), // Q105
    scores(0xbdbb9e5800017d64, 0xc039818a26c7e0ed), // Q106
    scores(0xbdb5f0280000f0a4, 0xc039bc7c5f86f1fd), // Q107
    scores(0xbdb16d10000097d5, 0xc039f76e9846030e), // Q108
    scores(0xbdabaf200000bf9a, 0xc03a3260d105141e), // Q109
    scores(0xbda5fd80000078e5, 0xc03a6d5309c4252e), // Q110
    scores(0xbda177b000004c47, 0xc03aa8454283363e), // Q111
    scores(0xbd9bc00000006042, 0xc03ae3377b42474e), // Q112
    scores(0xbd960ae000003cbc, 0xc03b1e29b4015860), // Q113
    scores(0xbd91824000002652, 0xc03b591becc06970), // Q114
    scores(0xbd8bd0c00000305b, 0xc03b940e257f7a80), // Q115
    scores(0xbd86184000001e83, 0xc03bcf005e3e8b90), // Q116
    scores(0xbd818d0000001340, 0xc03c09f296fd9ca0), // Q117
    scores(0xbd7be1800000184b, 0xc03c44e4cfbcadb1), // Q118
    scores(0xbd76258000000f54, 0xc03c7fd7087bbec1), // Q119
    scores(0xbd719780000009ac, 0xc03cbac9413acfd1), // Q120
    scores(0xbd6bf30000000c35, 0xc03cf5bb79f9e0e1), // Q121
    scores(0xbd663300000007b3, 0xc03d30adb2b8f1f1), // Q122
    scores(0xbd61a200000004dc, 0xc03d6b9feb780302), // Q123
    scores(0xbd5c040000000622, 0xc03da69224371412), // Q124
    scores(0xbd564000000003de, 0xc03de1845cf62522), // Q125
    scores(0xbd51ae0000000271, 0xc03e1c7695b53632), // Q126
    scores(0xbd4c140000000314, 0xc03e5768ce744742), // Q127
    scores(0xbd465000000001f2, 0xc03e925b07335854), // Q128
    scores(0xbd41b8000000013a, 0xc03ecd4d3ff26964), // Q129
    scores(0xbd3c28000000018c, 0xc03f083f78b17a74), // Q130
    scores(0xbd365800000000fa, 0xc03f4331b1708b84), // Q131
    scores(0xbd31c0000000009e, 0xc03f7e23ea2f9c94), // Q132
    scores(0xbd2c3000000000c7, 0xc03fb91622eeada5), // Q133
    scores(0xbd2670000000007e, 0xc03ff4085badbeb5), // Q134
    scores(0xbd21d0000000004f, 0xc040177d4a3667e2), // Q135
    scores(0xbd1c400000000064, 0xc04034f66695f06a), // Q136
    scores(0xbd1680000000003f, 0xc040526f82f578f3), // Q137
    scores(0xbd11e00000000028, 0xc0406fe89f55017b), // Q138
    scores(0xbd0c400000000032, 0xc0408d61bbb48a03), // Q139
    scores(0xbd06800000000020, 0xc040aadad814128b), // Q140
    scores(0xbd02000000000014, 0xc040c853f4739b13), // Q141
    scores(0xbcfc800000000019, 0xc040e5cd10d3239b), // Q142
    scores(0xbcf6800000000010, 0xc04103462d32ac24), // Q143
    scores(0xbcf200000000000a, 0xc04120bf499234ac), // Q144
    scores(0xbcec00000000000c, 0xc0413e3865f1bd34), // Q145
    scores(0xbce7000000000008, 0xc0415bb1825145bc), // Q146
    scores(0xbce2000000000005, 0xc041792a9eb0ce44), // Q147
    scores(0xbcdc000000000006, 0xc04196a3bb1056cc), // Q148
    scores(0xbcd6000000000004, 0xc041b41cd76fdf54), // Q149
    scores(0xbcd2000000000003, 0xc041d195f3cf67dc), // Q150
    scores(0xbccc000000000003, 0xc041ef0f102ef065), // Q151
    scores(0xbcc8000000000002, 0xc0420c882c8e78ed), // Q152
    scores(0xbcc4000000000002, 0xc0422a0148ee0175), // Q153
    scores(0xbcc0000000000001, 0xc042477a654d89fd), // Q154
    scores(0xbcb8000000000001, 0xc04264f381ad1285), // Q155
    scores(0xbcb0000000000001, 0xc042826c9e0c9b0d), // Q156
    scores(0xbcb0000000000001, 0xc0429fe5ba6c2395), // Q157
    scores(0xbca0000000000000, 0xc042bd5ed6cbac1e), // Q158
    scores(0xbca0000000000000, 0xc042dad7f32b34a6), // Q159
    scores(0xbca0000000000000, 0xc042f8510f8abd2e), // Q160
    scores(0xbca0000000000000, 0xc04315ca2bea45b6), // Q161
    scores(0xbca0000000000000, 0xc04333434849ce3e), // Q162
    scores(0x0000000000000000, 0xc04350bc64a956c6), // Q163
    scores(0x0000000000000000, 0xc0436e358108df4e), // Q164
    scores(0x0000000000000000, 0xc0438bae9d6867d7), // Q165
    scores(0x0000000000000000, 0xc043a927b9c7f05f), // Q166
    scores(0x0000000000000000, 0xc043c6a0d62778e7), // Q167
    scores(0x0000000000000000, 0xc043e419f287016f), // Q168
    scores(0x0000000000000000, 0xc04401930ee689f7), // Q169
    scores(0x0000000000000000, 0xc0441f0c2b46127f), // Q170
    scores(0x0000000000000000, 0xc0443c8547a59b08), // Q171
    scores(0x0000000000000000, 0xc04459fe6405238f), // Q172
    scores(0x0000000000000000, 0xc04477778064ac18), // Q173
    scores(0x0000000000000000, 0xc04494f09cc4349f), // Q174
    scores(0x0000000000000000, 0xc044b269b923bd28), // Q175
    scores(0x0000000000000000, 0xc044cfe2d58345b0), // Q176
    scores(0x0000000000000000, 0xc044ed5bf1e2ce38), // Q177
    scores(0x0000000000000000, 0xc0450ad50e4256c1), // Q178
    scores(0x0000000000000000, 0xc045284e2aa1df48), // Q179
    scores(0x0000000000000000, 0xc04545c7470167d1), // Q180
    scores(0x0000000000000000, 0xc04563406360f059), // Q181
    scores(0x0000000000000000, 0xc04580b97fc078e1), // Q182
    scores(0x0000000000000000, 0xc0459e329c200169), // Q183
    scores(0x0000000000000000, 0xc045bbabb87f89f1), // Q184
    scores(0x0000000000000000, 0xc045d924d4df1279), // Q185
    scores(0x0000000000000000, 0xc045f69df13e9b02), // Q186
    scores(0x0000000000000000, 0xc04614170d9e2389), // Q187
    scores(0x0000000000000000, 0xc046319029fdac12), // Q188
    scores(0x0000000000000000, 0xc0464f09465d3499), // Q189
    scores(0x0000000000000000, 0xc0466c8262bcbd22), // Q190
    scores(0x0000000000000000, 0xc04689fb7f1c45aa), // Q191
    scores(0x0000000000000000, 0xc046a7749b7bce32), // Q192
    scores(0x0000000000000000, 0xc046c4edb7db56ba), // Q193
    scores(0x0000000000000000, 0xc046e266d43adf42), // Q194
    scores(0x0000000000000000, 0xc046ffdff09a67cb), // Q195
    scores(0x0000000000000000, 0xc0471d590cf9f053), // Q196
    scores(0x0000000000000000, 0xc0473ad2295978db), // Q197
    scores(0x0000000000000000, 0xc047584b45b90163), // Q198
    scores(0x0000000000000000, 0xc04775c4621889eb), // Q199
    scores(0x0000000000000000, 0xc047933d7e781273), // Q200
    scores(0x0000000000000000, 0xc047b0b69ad79afc), // Q201
    scores(0x0000000000000000, 0xc047ce2fb7372383), // Q202
    scores(0x0000000000000000, 0xc047eba8d396ac0c), // Q203
    scores(0x0000000000000000, 0xc0480921eff63493), // Q204
    scores(0x0000000000000000, 0xc048269b0c55bd1c), // Q205
    scores(0x0000000000000000, 0xc048441428b545a4), // Q206
    scores(0x0000000000000000, 0xc048618d4514ce2c), // Q207
    scores(0x0000000000000000, 0xc0487f06617456b5), // Q208
    scores(0x0000000000000000, 0xc0489c7f7dd3df3c), // Q209
    scores(0x0000000000000000, 0xc048b9f89a3367c5), // Q210
    scores(0x0000000000000000, 0xc048d771b692f04d), // Q211
    scores(0x0000000000000000, 0xc048f4ead2f278d5), // Q212
    scores(0x0000000000000000, 0xc0491263ef52015d), // Q213
    scores(0x0000000000000000, 0xc0492fdd0bb189e5), // Q214
    scores(0x0000000000000000, 0xc0494d562811126d), // Q215
    scores(0x0000000000000000, 0xc0496acf44709af6), // Q216
    scores(0x0000000000000000, 0xc049884860d0237d), // Q217
    scores(0x0000000000000000, 0xc049a5c17d2fac06), // Q218
    scores(0x0000000000000000, 0xc049c33a998f348d), // Q219
    scores(0x0000000000000000, 0xc049e0b3b5eebd16), // Q220
    scores(0x0000000000000000, 0xc049fe2cd24e459e), // Q221
    scores(0x0000000000000000, 0xc04a1ba5eeadce26), // Q222
    scores(0x0000000000000000, 0xc04a391f0b0d56af), // Q223
    scores(0x0000000000000000, 0xc04a5698276cdf36), // Q224
    scores(0x0000000000000000, 0xc04a741143cc67bf), // Q225
    scores(0x0000000000000000, 0xc04a918a602bf047), // Q226
    scores(0x0000000000000000, 0xc04aaf037c8b78cf), // Q227
    scores(0x0000000000000000, 0xc04acc7c98eb0157), // Q228
    scores(0x0000000000000000, 0xc04ae9f5b54a89df), // Q229
    scores(0x0000000000000000, 0xc04b076ed1aa1267), // Q230
    scores(0x0000000000000000, 0xc04b24e7ee099af0), // Q231
    scores(0x0000000000000000, 0xc04b42610a692377), // Q232
    scores(0x0000000000000000, 0xc04b5fda26c8ac00), // Q233
    scores(0x0000000000000000, 0xc04b7d5343283487), // Q234
    scores(0x0000000000000000, 0xc04b9acc5f87bd10), // Q235
    scores(0x0000000000000000, 0xc04bb8457be74599), // Q236
    scores(0x0000000000000000, 0xc04bd5be9846ce20), // Q237
    scores(0x0000000000000000, 0xc04bf337b4a656a9), // Q238
    scores(0x0000000000000000, 0xc04c10b0d105df30), // Q239
    scores(0x0000000000000000, 0xc04c2e29ed6567b9), // Q240
    scores(0x0000000000000000, 0xc04c4ba309c4f041), // Q241
    scores(0x0000000000000000, 0xc04c691c262478c9), // Q242
    scores(0x0000000000000000, 0xc04c869542840151), // Q243
    scores(0x0000000000000000, 0xc04ca40e5ee389d9), // Q244
    scores(0x0000000000000000, 0xc04cc1877b431261), // Q245
    scores(0x0000000000000000, 0xc04cdf0097a29aea), // Q246
    scores(0x0000000000000000, 0xc04cfc79b4022371), // Q247
    scores(0x0000000000000000, 0xc04d19f2d061abfa), // Q248
    scores(0x0000000000000000, 0xc04d376becc13481), // Q249
    scores(0x0000000000000000, 0xc04d54e50920bd0a), // Q250
    scores(0x0000000000000000, 0xc04d725e25804593), // Q251
    scores(0x0000000000000000, 0xc04d8fd741dfce1a), // Q252
    scores(0x0000000000000000, 0xc04dad505e3f56a3), // Q253
    scores(0x0000000000000000, 0xc04dcac97a9edf2a), // Q254
    scores(0x0000000000000000, 0xc04de84296fe67b3), // Q255
];

/// One row of [`PER_QUALITY_LN`], from the bit patterns of its two scores.
const fn scores(match_ln_bits: u64, mismatch_ln_bits: u64) -> BaseScores {
    BaseScores {
        match_ln: f64::from_bits(match_ln_bits),
        mismatch_ln: f64::from_bits(mismatch_ln_bits),
    }
}

/// Scores each base with the read's own quality — the default, and production's model.
///
/// Holds a borrow of the shared table rather than reaching for it per call, so holding one
/// costs a pointer. That is what lets an aligner take it by value as a type parameter.
///
/// **Why the borrow is a field.** [`Self::scores_for`] is called from inside the delimiters'
/// innermost DP loop. While the table was a `LazyLock`, a deref there was an acquire-load of the
/// once-guard plus a *live call site* to the initialiser — `cargo asm` on `classify::delimit`
/// showed seven of each inside one function, and the per-row constants were re-read from the
/// stack on every cell — and holding the borrow turned that into a register-resident pointer.
/// The table is a plain `static` now, so the guard is gone either way; the borrow stays because
/// it is still the cheapest thing to carry.
///
/// `Default` is written out rather than derived, deliberately: a derived `Default` is the
/// one construction path that would keep compiling if this gained a field, silently
/// zero-filling it, whereas [`Self::new`] would fail to compile and say so.
#[derive(Debug, Clone, Copy)]
pub struct PerQualityEmission {
    scores: &'static [BaseScores; 256],
}

impl PerQualityEmission {
    /// Build the per-quality model. The table itself is shared.
    #[must_use]
    pub fn new() -> Self {
        Self {
            scores: &PER_QUALITY_LN,
        }
    }
}

impl Default for PerQualityEmission {
    fn default() -> Self {
        Self::new()
    }
}

impl Emission for PerQualityEmission {
    #[inline]
    fn scores_for(&self, quality: BaseQual) -> BaseScores {
        self.scores[quality.get() as usize]
    }

    #[inline]
    fn insert_ln(&self) -> f64 {
        UNIFORM_BASE_LN
    }
}

/// Scores every base at one fixed error rate, ignoring the read's qualities — the
/// quality-blind end of the emission comparison (spec §3).
///
/// The rate is the sample group's per-base substitution rate, so it is the experiment's
/// configuration rather than a per-read fact; that is why it is constructor state while
/// the repeat geometry travels per call (arch §4).
///
/// Both scores are computed **once, at construction**, because a flat model's scores do
/// not depend on anything the matrix varies — which also keeps the per-cell work to one
/// comparison and one field read. There is deliberately **no `Default`**: a default error
/// rate would be exactly the kind of behaviourally-significant hidden value that has to be
/// visible at the call site.
#[derive(Debug, Clone, Copy)]
pub struct FlatEmission {
    scores: BaseScores,
}

impl FlatEmission {
    /// Build a flat model from a per-base error rate, `ε`, matching production's flat
    /// path (`align_subst`,
    /// [src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)): a match scores
    /// `1 − ε` and a mismatch `ε / 3`, converted to log space here.
    ///
    /// # Errors
    ///
    /// [`DomainError::ErrorRate`] if `error_rate` is not a finite probability in `[0, 1]`.
    ///
    /// # Why this one *is* checked, when nothing else in the module is
    ///
    /// Arch §3 bans `Result` across this module, and the reason is hot-path cost: a
    /// fallible alignment would push error handling onto the caller's hottest loop for
    /// cases that are answers, not failures. **That argument does not reach this
    /// constructor.** The rate is the run's fixed configuration, consumed *once* into two
    /// `f64` fields — checking it costs nothing per matrix cell, which is the only place
    /// the ban was defending. Arch §3's own escape clause names the remedy as "a checked
    /// constructor on the context type", and this type is that context type.
    ///
    /// The previous shape — a `debug_assert!` plus a release clamp — was defensible while
    /// the only callers were tests, but it left the guarantee compiled out of the build the
    /// project actually runs. Now the check is real in every profile. The clamping
    /// behaviour survives in [`Self::from_unchecked_rate`], which the totality test uses to
    /// exercise what a violated precondition would once have produced.
    ///
    /// The `PROBABILITY_FLOOR` at each endpoint is a separate thing and is **not** a check:
    /// it is ported model behaviour, binding only at `ε = 0` (which would make every
    /// mismatch impossible) and `ε = 1` (every match impossible). Within the accepted range
    /// it changes nothing.
    ///
    /// [`MismatchFraction`]: crate::types::MismatchFraction
    pub fn try_new(error_rate: f64) -> Result<Self, DomainError> {
        if !error_rate.is_finite() || !(0.0..=1.0).contains(&error_rate) {
            return Err(DomainError::ErrorRate(error_rate));
        }
        Ok(Self::from_unchecked_rate(error_rate))
    }

    /// Everything [`Self::new`] does **except** the debug assertion — that is, exactly
    /// what a release build runs.
    ///
    /// It exists to be testable. The totality contract has to hold once the assertion is
    /// compiled out, but a test cannot reach that path through `new`, because the
    /// assertion fires in the test profile. Testing the clamp through this function is
    /// the only way to assert the release behaviour without conditionally compiling the
    /// test out of the build anyone actually runs.
    fn from_unchecked_rate(error_rate: f64) -> Self {
        let error_rate = if error_rate.is_nan() {
            1.0
        } else {
            error_rate.clamp(0.0, 1.0)
        };
        Self {
            scores: BaseScores {
                match_ln: float::ln((1.0 - error_rate).max(PROBABILITY_FLOOR)),
                mismatch_ln: float::ln((error_rate / 3.0).max(PROBABILITY_FLOOR)),
            },
        }
    }
}

impl Emission for FlatEmission {
    #[inline]
    fn scores_for(&self, _quality: BaseQual) -> BaseScores {
        self.scores
    }

    /// `ln(1/4)` — the same value [`PerQualityEmission`] uses, and **a decision rather
    /// than a port** (arch §2.3).
    ///
    /// The reasoning is that an inserted base has no reference base in *either* model, so
    /// there is nothing for `ε` to describe: `ε` is the chance of misreading a base whose
    /// true identity the reference supplies, and an insertion has no such base. Scoring
    /// it against a uniform composition is the same modelling assumption in both, so the
    /// two agree here and differ only where they are meant to — matched bases.
    ///
    /// The alternative would be to follow production's flat path, whose unequal-length
    /// route charges `gap = eps` for a base absorbed in the flank slop. That is a
    /// **transition** cost, not an emission, and putting it in this slot would price the
    /// same event twice once the aligner's own gap model is applied on top.
    #[inline]
    fn insert_ln(&self) -> f64 {
        UNIFORM_BASE_LN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Phred quality to error probability, spelled out independently of the table under
    /// test so the test cannot pass by sharing the implementation's mistake.
    fn error_probability(quality: u8) -> f64 {
        float::powf(10.0, -f64::from(quality) / 10.0)
    }

    /// **The written-out table is the Dindel model.** Re-derived from the formula through
    /// [`crate::float`], spelled independently of the table. The table holds the bits glibc
    /// produced, and libm rounds some entries' last place differently, so every entry must agree
    /// to within `4 · f64::EPSILON` of its magnitude, or of 1.0 where the magnitude is smaller —
    /// high qualities' match scores are near zero but are computed from `1 − ε`, whose rounding
    /// step is about `1.1e-16` in absolute terms. That still catches a wrong row (adjacent rows'
    /// mismatch scores differ by about 0.23), a transposed pair or a changed floor. The exact bits,
    /// which the repeat-aware aligner breaks ties on, are fixed by the literal itself and are the
    /// same on every platform.
    #[test]
    fn per_quality_table_matches_the_dindel_model() {
        let emission = PerQualityEmission::new();
        for quality in 0..=u8::MAX {
            let error_rate = error_probability(quality);
            let expected = BaseScores {
                match_ln: float::ln((1.0 - error_rate).max(f64::MIN_POSITIVE)),
                mismatch_ln: float::ln((error_rate / 3.0).max(f64::MIN_POSITIVE)),
            };
            let written = emission.scores_for(BaseQual(quality));
            for (name, got, want) in [
                ("match", written.match_ln, expected.match_ln),
                ("mismatch", written.mismatch_ln, expected.mismatch_ln),
            ] {
                assert!(
                    (got - want).abs() <= 4.0 * f64::EPSILON * want.abs().max(1.0),
                    "table's {name} score diverged from the Dindel model at Q{quality}: \
                     {got:e} against {want:e}"
                );
            }
        }
    }

    /// The published Dindel model at qualities where the arithmetic is checkable by hand:
    /// Q10 → ε = 0.1, Q20 → 0.01, Q30 → 0.001. Tolerance is appropriate here — these are
    /// decimal literals a human can verify, not the bit-level contract, which
    /// `per_quality_table_matches_the_dindel_model` carries.
    #[test]
    fn per_quality_emission_reproduces_the_dindel_model() {
        let emission = PerQualityEmission::new();
        for (quality, error_rate) in [(10u8, 0.1f64), (20, 0.01), (30, 0.001)] {
            assert!(
                (emission.emit_ln(b'A', b'A', BaseQual(quality)) - float::ln(1.0 - error_rate))
                    .abs()
                    < 1e-12,
                "match score at Q{quality}"
            );
            assert!(
                (emission.emit_ln(b'A', b'C', BaseQual(quality)) - float::ln(error_rate / 3.0))
                    .abs()
                    < 1e-12,
                "mismatch score at Q{quality}"
            );
        }
    }

    /// **The table's exact bits.** The tolerance test above cannot see an entry moved by a few units
    /// in the last place, and the repeat-aware aligner breaks ties on exactly those units. A digest
    /// over all 512 bit patterns — FNV-1a over the 64-bit words, match then mismatch, Q0 to Q255 —
    /// recorded from the table as `75722336` wrote it, fails on any single changed bit.
    #[test]
    fn per_quality_table_keeps_its_written_bits() {
        let emission = PerQualityEmission::new();
        let digest = (0..=u8::MAX)
            .map(|quality| emission.scores_for(BaseQual(quality)))
            .flat_map(|scores| [scores.match_ln.to_bits(), scores.mismatch_ln.to_bits()])
            .fold(0xcbf2_9ce4_8422_2325_u64, |hash, word| {
                (hash ^ word).wrapping_mul(0x0000_0100_0000_01b3)
            });
        assert_eq!(digest, WRITTEN_TABLE_DIGEST, "digest {digest:#018x}");
    }

    const WRITTEN_TABLE_DIGEST: u64 = 0x0454_3206_1128_dbc0;

    /// **The floor, which is the reason this is a table.** At Q0 the error probability is
    /// 1, so an unfloored `ln(1 − ε)` is `ln(0) = -inf`, and one such cell annihilates
    /// every path through it — the read stops being alignable because of a single
    /// worthless base. This test fails if the floor is dropped.
    #[test]
    fn per_quality_emission_floors_the_quality_zero_match_instead_of_annihilating_it() {
        let emission = PerQualityEmission::new();
        let match_ln = emission.emit_ln(b'A', b'A', BaseQual(0));

        assert!(
            match_ln.is_finite(),
            "a Q0 match must be finite, was {match_ln}"
        );
        // Checked against the value independently, not against `PROBABILITY_FLOOR` —
        // otherwise editing the constant would move both sides of the assertion. The
        // floor is documented as the smallest positive normal f64, ~2.2250738585072014e-308.
        assert_eq!(PROBABILITY_FLOOR, 2.225_073_858_507_201_4e-308);
        assert_eq!(match_ln, float::ln(2.225_073_858_507_201_4e-308));
        assert!((match_ln - -708.396_418_532_264_1).abs() < 1e-9);
        // A Q0 mismatch (ε/3 = 1/3) stays an ordinary, unfloored score.
        assert!((emission.emit_ln(b'A', b'C', BaseQual(0)) - float::ln(1.0 / 3.0)).abs() < 1e-12);
    }

    /// The port floors the mismatch term where production does not. That is safe only
    /// because the floor cannot bind there — assert it rather than assume it, since the
    /// claim is what makes the two tables bit-identical.
    #[test]
    fn mismatch_floor_never_binds_over_the_quality_domain() {
        for quality in 0..=u8::MAX {
            let smallest = error_probability(quality) / 3.0;
            assert!(
                smallest > PROBABILITY_FLOOR,
                "the mismatch floor would bind at Q{quality} ({smallest:e})"
            );
        }
        // The tightest case, at the top of the domain, with a wide margin.
        assert!(error_probability(u8::MAX) / 3.0 > 1e-27);
    }

    /// Totality, over the whole domain rather than at sampled points: every one of the
    /// 256 quality bytes a BAM can hold must give a finite score for both outcomes, in
    /// both implementations. This is the contract the trait states.
    #[test]
    fn emission_is_finite_at_every_quality() {
        let per_quality = PerQualityEmission::new();
        let flat = FlatEmission::try_new(0.01).expect("a valid test error rate");
        for quality in 0..=u8::MAX {
            for score in [
                per_quality.emit_ln(b'A', b'A', BaseQual(quality)),
                per_quality.emit_ln(b'A', b'C', BaseQual(quality)),
                flat.emit_ln(b'A', b'A', BaseQual(quality)),
                flat.emit_ln(b'A', b'C', BaseQual(quality)),
            ] {
                assert!(score.is_finite(), "non-finite score at Q{quality}");
                assert!(score <= 0.0, "a log-probability above zero at Q{quality}");
            }
        }
        assert!(per_quality.insert_ln().is_finite());
        assert!(flat.insert_ln().is_finite());
    }

    /// **The check is now real in every build profile**, which is the point of `try_new`:
    /// the previous `debug_assert!` was compiled out of the release build this project
    /// actually runs, so the guarantee existed only in tests.
    #[test]
    fn try_new_rejects_a_rate_that_is_not_a_probability() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.5, 1.5, 1e300] {
            assert!(
                matches!(FlatEmission::try_new(bad), Err(DomainError::ErrorRate(_))),
                "rate {bad} was not rejected"
            );
        }
        // Both endpoints are legal — degenerate, not invalid, and the floors make them safe.
        assert!(FlatEmission::try_new(0.0).is_ok());
        assert!(FlatEmission::try_new(1.0).is_ok());
        assert!(FlatEmission::try_new(0.01).is_ok());
    }

    /// **Totality must survive a violated precondition**, because the debug assertion
    /// that guards it is compiled out of the release build this project runs. Each of
    /// these rates would, without the clamp, produce a score that inverts the model
    /// rather than merely skewing it: `+inf` mismatches, or a match score above zero.
    ///
    /// Goes through `from_unchecked_rate` rather than `new` because `new`'s debug
    /// assertion fires in the test profile — the path under test is precisely the one
    /// that runs when that assertion is gone.
    #[test]
    fn flat_emission_stays_total_for_rates_outside_the_contract() {
        for error_rate in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            -0.5,
            1.5,
            f64::MAX,
        ] {
            let flat = FlatEmission::from_unchecked_rate(error_rate);
            for score in [
                flat.emit_ln(b'A', b'A', BaseQual(30)),
                flat.emit_ln(b'A', b'C', BaseQual(30)),
            ] {
                assert!(score.is_finite(), "non-finite score for rate {error_rate}");
                assert!(
                    score <= 0.0,
                    "log-probability above zero ({score}) for rate {error_rate}"
                );
            }
        }
    }

    /// A match must outscore a mismatch at every usable quality. If that ever inverts,
    /// the aligner prefers disagreeing with the reference and every result is nonsense.
    ///
    /// **Q0 and Q1 are genuine exceptions, and the boundary is asserted rather than
    /// assumed.** A match beats a mismatch exactly when `1 − ε > ε / 3`, i.e. when
    /// `ε < 0.75`, which is `Q > 1.249`. So the model says a mismatch is the *likelier*
    /// reading at Q0 and Q1 — correctly: at ε ≈ 0.79 the base is very nearly noise, and
    /// there are three ways to disagree against one way to agree. This is a property of
    /// the Dindel model being ported, not a defect, but it is exactly the kind of edge a
    /// caller would assume away, so both sides of the boundary are pinned here.
    #[test]
    fn a_match_outscores_a_mismatch_above_the_quality_one_crossover() {
        let emission = PerQualityEmission::new();

        for quality in 2..=u8::MAX {
            assert!(
                emission.emit_ln(b'A', b'A', BaseQual(quality))
                    > emission.emit_ln(b'A', b'C', BaseQual(quality)),
                "match did not outscore mismatch at Q{quality}"
            );
        }

        for quality in [0u8, 1] {
            assert!(
                emission.emit_ln(b'A', b'A', BaseQual(quality))
                    < emission.emit_ln(b'A', b'C', BaseQual(quality)),
                "expected the mismatch to win below the crossover, at Q{quality}"
            );
        }
    }

    /// The same ordering for the flat model. Its `ε` is a **constructor parameter**, so
    /// unlike the quality table the crossover at `ε = 0.75` is reachable by
    /// configuration — and without this test, swapping the two structurally identical
    /// lines that build `match_ln` and `mismatch_ln` would pass the whole suite.
    #[test]
    fn flat_emission_orders_match_above_mismatch_below_the_crossover() {
        for error_rate in [0.0, 0.001, 0.01, 0.1, 0.5, 0.74] {
            let flat = FlatEmission::try_new(error_rate).expect("a valid test error rate");
            assert!(
                flat.emit_ln(b'A', b'A', BaseQual(30)) > flat.emit_ln(b'A', b'C', BaseQual(30)),
                "match did not outscore mismatch at ε = {error_rate}"
            );
        }
        // Above the crossover the order reverses, for the same reason it does at Q0/Q1.
        for error_rate in [0.76, 0.9, 1.0] {
            let flat = FlatEmission::try_new(error_rate).expect("a valid test error rate");
            assert!(
                flat.emit_ln(b'A', b'A', BaseQual(30)) < flat.emit_ln(b'A', b'C', BaseQual(30)),
                "expected the mismatch to win at ε = {error_rate}"
            );
        }
    }

    #[test]
    fn flat_emission_scores_every_quality_alike() {
        let flat = FlatEmission::try_new(0.02).expect("a valid test error rate");
        let at_low = flat.emit_ln(b'A', b'A', BaseQual(2));
        let at_high = flat.emit_ln(b'A', b'A', BaseQual(60));
        assert_eq!(at_low, at_high);
        assert!((at_low - float::ln(0.98)).abs() < 1e-12);
        assert!((flat.emit_ln(b'A', b'G', BaseQual(2)) - float::ln(0.02 / 3.0)).abs() < 1e-12);
    }

    /// The flat model's degenerate endpoints are floored rather than annihilating, the
    /// same discipline the quality table applies at Q0.
    #[test]
    fn flat_emission_floors_both_degenerate_error_rates() {
        let perfect = FlatEmission::try_new(0.0).expect("a valid test error rate");
        assert_eq!(perfect.emit_ln(b'A', b'A', BaseQual(30)), 0.0); // ln(1)
        assert_eq!(
            perfect.emit_ln(b'A', b'C', BaseQual(30)),
            float::ln(PROBABILITY_FLOOR)
        );

        let hopeless = FlatEmission::try_new(1.0).expect("a valid test error rate");
        assert_eq!(
            hopeless.emit_ln(b'A', b'A', BaseQual(30)),
            float::ln(PROBABILITY_FLOOR)
        );
        assert!(hopeless.emit_ln(b'A', b'C', BaseQual(30)).is_finite());
    }

    /// Both implementations score an inserted base against a uniform base composition.
    /// They agree deliberately — see `FlatEmission::insert_ln` for why this is a decision
    /// rather than two ports that happened to coincide.
    #[test]
    fn both_implementations_score_an_inserted_base_against_a_uniform_composition() {
        assert_eq!(PerQualityEmission::new().insert_ln(), UNIFORM_BASE_LN);
        assert_eq!(
            FlatEmission::try_new(0.01)
                .expect("a valid test error rate")
                .insert_ln(),
            UNIFORM_BASE_LN
        );
    }

    /// The literal is written out so it costs nothing on the hot path; this is what keeps
    /// it honest.
    #[test]
    fn uniform_base_ln_is_ln_of_a_quarter() {
        assert_eq!(UNIFORM_BASE_LN, float::ln(0.25));
    }

    /// **The caller's precondition, pinned so it is recorded rather than discovered.**
    /// Bases are compared by raw byte equality, so a soft-masked (lower-case) reference
    /// mismatches everywhere, and `N` against `N` is a full-confidence match. Both follow
    /// production; neither is fixed here, because doing so would be a scoring-model change
    /// smuggled into a component.
    #[test]
    fn emission_compares_bases_by_raw_byte_equality() {
        let emission = PerQualityEmission::new();
        let quality = BaseQual(30);
        let matched = emission.emit_ln(b'A', b'A', quality);
        let mismatched = emission.emit_ln(b'A', b'C', quality);

        // Soft-masked reference: same base, different case, scored as a mismatch.
        assert_eq!(emission.emit_ln(b'A', b'a', quality), mismatched);
        assert_eq!(emission.emit_ln(b'a', b'a', quality), matched);

        // `N` is not special: it matches itself at full confidence.
        assert_eq!(emission.emit_ln(b'N', b'N', quality), matched);
        assert_eq!(emission.emit_ln(b'N', b'A', quality), mismatched);
    }

    /// Emission depends on whether the bases agree, never on which bases they are — the
    /// model has no per-base composition beyond that.
    #[test]
    fn emission_depends_only_on_whether_the_bases_agree() {
        let emission = PerQualityEmission::new();
        let quality = BaseQual(25);
        for (read_base, reference_base) in [(b'A', b'A'), (b'C', b'C'), (b'G', b'G'), (b'T', b'T')]
        {
            assert_eq!(
                emission.emit_ln(read_base, reference_base, quality),
                emission.emit_ln(b'A', b'A', quality)
            );
        }
        for (read_base, reference_base) in [(b'A', b'C'), (b'C', b'A'), (b'G', b'T'), (b'T', b'N')]
        {
            assert_eq!(
                emission.emit_ln(read_base, reference_base, quality),
                emission.emit_ln(b'A', b'C', quality)
            );
        }
    }

    /// The row-resolved and per-call forms must agree — `emit_ln` is a provided method
    /// over `scores_for`, and a caller that hoists must get the same scores as one that
    /// does not.
    #[test]
    fn row_resolved_scores_agree_with_the_per_call_form() {
        let per_quality = PerQualityEmission::new();
        let flat = FlatEmission::try_new(0.01).expect("a valid test error rate");
        for quality in 0..=u8::MAX {
            let quality = BaseQual(quality);
            let row = per_quality.scores_for(quality);
            assert_eq!(
                row.pick(b'A', b'A'),
                per_quality.emit_ln(b'A', b'A', quality)
            );
            assert_eq!(
                row.pick(b'A', b'C'),
                per_quality.emit_ln(b'A', b'C', quality)
            );

            let flat_row = flat.scores_for(quality);
            assert_eq!(flat_row.pick(b'A', b'A'), flat.emit_ln(b'A', b'A', quality));
            assert_eq!(flat_row.pick(b'A', b'C'), flat.emit_ln(b'A', b'C', quality));
        }
    }
}
