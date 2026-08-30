//! **The two derived bindings of spec §3.1, in the spelling `[fitted_from]` carries them in** —
//! the reference's content digest and the census the fit read.
//!
//! **Two of the four and not all four.** The other two are names: the sample list and the
//! read-group table are read straight off the run's own `ReadGroups` in `from_run_parameters`,
//! and there is nothing to derive. These two are derived values, and until this step neither was
//! derived anywhere.
//!
//! The *shape* of that section landed with the rest of the file (step A1), the writer fills it
//! (B1) and the reader reads it back (C2). What was missing until here is the step from a run's
//! own inputs to the two bindings that are not names: `of_run` took the reference digest as
//! **text**, with its own documentation saying "nothing here can check that it is one", and took
//! a census identity nothing in this tree could build. **The census identity is minted here and
//! the reference digest is spelled here** — `ReferenceDigest` is computed upstream, over the
//! reference itself — so that the value a run writes and the value a later run compares against
//! come out of one function rather than out of two callers agreeing.
//!
//! # Why a digest and not the value
//!
//! Every term of the census is written as a digest rather than as what it holds. **Seven of the
//! twelve arrive that way**, because the census made that choice one level down: the selection
//! terms hold two other modules' whole configuration, and their only use is an equality
//! ([`SelectionTermsDigest`](crate::ng::parameter_estimation::joint::census::SelectionTermsDigest)).
//! **The other five are digested here, and that is this file's choice rather than an inherited
//! one** — the census file writes the per-stratum locus counts, the read cap and the depth cap
//! as values (`census_file.rs`, `encode_header`). One rule across all twelve is what is bought,
//! and what it costs is that a file can say *whether* it matches and not *what it was built at*.
//!
//! **Twelve terms and not one digest over them**, for the reason the shape's own documentation
//! gives: a mismatch has to be *named*, because every one of these fails the same way — silently
//! — and "the terms differ" is not something anyone can act on.
//!
//! # The names are the census's, and so is their order
//!
//! A term is named in the words
//! [`RecordingTerms::first_disagreement`](crate::ng::parameter_estimation::joint::census::RecordingTerms::first_disagreement)
//! uses, because that is the sentence the fit already prints when two samples disagree, and a
//! second vocabulary for the same twelve values would mean a reader met one word here and
//! another there. Seven arrive already named, from the selection's own table; the other five are
//! written out below.
//!
//! **The order is that function's checking order, and it is load-bearing rather than tidy**: it
//! reports the *first* value two censuses disagree on, so where two differ, the order decides
//! which one a run names. `every_term_is_named_as_the_census_names_it` moves one value at a time
//! and `the_terms_are_in_the_order_the_census_checks_them` moves two, which is the only way to
//! see an order at all.

use md5::{Digest, Md5};

use super::{CensusIdentity, CensusTerm};
use crate::ng::parameter_estimation::joint::census::RecordingTerms;

#[cfg(test)]
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
#[cfg(test)]
use crate::ng::parameter_estimation::joint::census::{
    DepthCap, DepthLadderDigest, ReadCap, SelectionTermsDigest,
};
#[cfg(test)]
use crate::ng::parameter_estimation::joint::loci::{
    BlockDigest, CatalogBuildSettings, CensusLociDigest, ReferenceDigest, RegionSetDigest,
    SelectionTerms,
};
#[cfg(test)]
use crate::ng::repeat_catalog::{StrRepeatCriteria, StratumCounts};
#[cfg(test)]
use crate::ng::tandem_repeat::ScanParams;
#[cfg(test)]
use crate::ng::types::{Bp, ContigId};

impl CensusIdentity {
    /// **The census a fit ran under, as the file names it** — one term a value, digested.
    ///
    /// The order is [`RecordingTerms::first_disagreement`]'s own checking order: the seven
    /// selection values first, in the order that type compares them, then the five that say
    /// what came back and in what units.
    #[must_use]
    pub fn of(terms: &RecordingTerms) -> Self {
        // **Destructured without `..` on purpose**, which is this struct's own convention
        // upstream: a value added to `RecordingTerms` stops this compiling rather than quietly
        // dropping out of the identity, and a value that drops out lets a file fitted under
        // other terms read back as this run's own without a word.
        let RecordingTerms {
            selection,
            kept_loci,
            ssr_stratum_counts,
            read_cap,
            depth_ladder,
            depth_cap,
        } = terms;

        let mut named = Vec::with_capacity(selection.fields().len() + 5);
        for (term, digest) in selection.fields() {
            named.push(CensusTerm {
                term: (*term).to_owned(),
                digest: hex_digest(digest),
            });
        }
        // **The whole digest and every block**, because that is what this value's own equality
        // compares: two censuses agreeing on the whole and differing in a block are two
        // censuses the fit refuses to pool, so an identity built from the whole alone would
        // call them the same.
        named.push(CensusTerm {
            term: "the loci actually kept".to_owned(),
            digest: a_digest_over(|hasher| {
                hasher.update(kept_loci.whole());
                for block in kept_loci.blocks() {
                    hasher.update(block.contig.get().to_le_bytes());
                    hasher.update(block.megabase.to_le_bytes());
                    hasher.update(block.digest.to_le_bytes());
                }
            }),
        });
        named.push(CensusTerm {
            term: "per-stratum locus counts".to_owned(),
            digest: a_digest_over(|hasher| {
                for ((period, reference_repeats), loci) in ssr_stratum_counts.iter_sorted() {
                    hasher.update(period.to_le_bytes());
                    hasher.update(reference_repeats.to_le_bytes());
                    hasher.update(loci.to_le_bytes());
                }
            }),
        });
        // **The two caps are digested like everything else, though each is one small integer.**
        // Writing them as numbers would be friendlier to read, and this file's own rule is that
        // a key a person can read is a key a person can edit — which is exactly wrong for a
        // binding, whose whole use is an equality nobody should be able to satisfy by typing.
        // One rule for all twelve terms, and `[fitted_from]`'s note says they are not editable.
        named.push(CensusTerm {
            term: "per-locus read cap".to_owned(),
            digest: a_digest_over(|hasher| hasher.update(read_cap.0.to_le_bytes())),
        });
        named.push(CensusTerm {
            term: "depth ladder edges".to_owned(),
            digest: hex_digest(&depth_ladder.0),
        });
        named.push(CensusTerm {
            term: "per-position depth cap".to_owned(),
            digest: a_digest_over(|hasher| hasher.update(depth_cap.get().to_le_bytes())),
        });

        Self { terms: named }
    }
}

/// A digest of whatever is fed to it, in the file's spelling.
fn a_digest_over(feed: impl FnOnce(&mut Md5)) -> String {
    let mut hasher = Md5::new();
    feed(&mut hasher);
    let digest: [u8; 16] = hasher.finalize().into();
    hex_digest(&digest)
}

/// **The file's spelling of a 16-byte digest**: 32 characters of lower-case hex.
///
/// The one place either binding's text is produced, so the string a run writes and the string a
/// later run compares it against cannot be two spellings of the same bytes.
pub(super) fn hex_digest(digest: &[u8; 16]) -> String {
    use std::fmt::Write as _;

    digest
        .iter()
        .fold(String::with_capacity(32), |mut hex, byte| {
            write!(hex, "{byte:02x}").expect("a string never fails");
            hex
        })
}

// ---------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------

/// **A run's selection terms, as a test needs them.**
///
/// **No two of the seven render alike**, so a digest taken over the wrong field is visible.
/// `catalog_built_under.criteria` is deliberately *not* `StrRepeatCriteria::default()`: it and
/// `ssr_criteria` are the same type and mean different things — what the catalog file was built
/// at, against what this run asked that file for — and a fixture giving both the default would
/// see nothing if a digest read one for the other.
#[cfg(test)]
fn a_runs_selection_terms() -> SelectionTerms {
    let built_under = StrRepeatCriteria {
        min_flank_bp: Bp(StrRepeatCriteria::default().min_flank_bp.get() + 5),
        ..StrRepeatCriteria::default()
    };
    SelectionTerms {
        seed: 42,
        reference: ReferenceDigest([7; 16]),
        analysed_regions: RegionSetDigest([9; 16]),
        catalog_built_under: CatalogBuildSettings {
            criteria: built_under,
            scan: ScanParams::default(),
            tool_version: "0.1.0".to_string(),
        },
        ssr_criteria: StrRepeatCriteria::default(),
        generic_target: 2_000_000,
        ssr_cap: 1_000,
    }
}

/// The digest a walk takes over every kept locus at once.
#[cfg(test)]
const THE_WHOLE_A_WALK_DIGESTED: [u8; 16] = [3; 16];

/// The per-megabase blocks a walk kept, **two of them on different contigs**, so a digest that
/// drops the contig — leaving two blocks that differ only in which chromosome they are on — has
/// something to miss.
#[cfg(test)]
fn the_blocks_a_walk_kept() -> Vec<BlockDigest> {
    vec![
        BlockDigest {
            contig: ContigId(1),
            megabase: 4,
            digest: 0x0102_0304_0506_0708,
        },
        BlockDigest {
            contig: ContigId(2),
            megabase: 7,
            digest: 0x1112_1314_1516_1718,
        },
    ]
}

/// The kept loci as a walk would report them — a whole digest **and** its blocks, which is what
/// this value's own equality compares.
#[cfg(test)]
fn the_loci_a_walk_kept() -> CensusLociDigest {
    CensusLociDigest::from_parts(THE_WHOLE_A_WALK_DIGESTED, the_blocks_a_walk_kept())
}

/// **The per-stratum locus counts a walk produced — two strata, and non-empty on purpose.**
///
/// An empty table would let every edit below pass for the wrong reason: going from no bytes to
/// some bytes moves a digest whatever of the entry survives, so a digest that dropped the period
/// or the count would still look alive.
#[cfg(test)]
fn the_strata_a_walk_counted() -> StratumCounts {
    StratumCounts::from_counted([((2, 6), 3), ((3, 5), 7)])
}

/// **A census's recording terms, as a test needs them.**
///
/// `pub(super)` because the module's shared file fixture is built from the identity these mint,
/// with the digests replaced — see `a_census_a_run_could_have_fitted_under` in `mod.rs`.
#[cfg(test)]
pub(super) fn a_censuss_recording_terms() -> RecordingTerms {
    RecordingTerms {
        selection: SelectionTermsDigest::of(&a_runs_selection_terms()),
        kept_loci: the_loci_a_walk_kept(),
        ssr_stratum_counts: the_strata_a_walk_counted(),
        read_cap: ReadCap(100),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
        depth_cap: DepthCap::new(124),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One of the twelve values moved, and what moved.
    ///
    /// **Two arguments rather than one**, so that edits *compose*: seven of the twelve live
    /// inside `SelectionTerms` and reach `RecordingTerms` only through
    /// `SelectionTermsDigest::of`, so a closure that rebuilt the digest from a fresh selection
    /// would discard whatever an earlier closure had done. `moved` applies them in order and
    /// digests once at the end, which is what lets a test move *two* values.
    type OneValueMoved = (&'static str, fn(&mut SelectionTerms, &mut RecordingTerms));

    /// **The twelve values `RecordingTerms::first_disagreement` checks, in its order**, each
    /// with an edit that moves it and nothing else.
    fn the_twelve_edits() -> Vec<OneValueMoved> {
        vec![
            ("the selection seed", |selection, _| selection.seed = 43),
            ("the reference", |selection, _| {
                selection.reference = ReferenceDigest([8; 16]);
            }),
            ("the analysed regions", |selection, _| {
                selection.analysed_regions = RegionSetDigest([10; 16]);
            }),
            ("what the catalog was built at", |selection, _| {
                selection.catalog_built_under.tool_version = "0.2.0".to_string();
            }),
            ("what this run asked the catalog for", |selection, _| {
                selection.ssr_criteria.max_str_len_bp =
                    Bp(selection.ssr_criteria.max_str_len_bp.get() + 1);
            }),
            ("the generic target", |selection, _| {
                selection.generic_target = 3_000_000;
            }),
            ("the per-stratum cap", |selection, _| {
                selection.ssr_cap = 2_000
            }),
            // **A block and not the whole**, which is the half an identity built from
            // `kept_loci.whole()` alone would miss: the census compares both.
            ("one megabase of the kept loci", |_, terms| {
                let mut blocks = the_blocks_a_walk_kept();
                blocks[0].megabase += 1;
                terms.kept_loci = CensusLociDigest::from_parts(THE_WHOLE_A_WALK_DIGESTED, blocks);
            }),
            ("one stratum's locus count", |_, terms| {
                terms.ssr_stratum_counts = StratumCounts::from_counted([((2, 6), 4), ((3, 5), 7)]);
            }),
            ("the read cap", |_, terms| terms.read_cap = ReadCap(101)),
            ("the depth ladder", |_, terms| {
                terms.depth_ladder = DepthLadderDigest([1; 16]);
            }),
            ("the depth cap", |_, terms| {
                terms.depth_cap = DepthCap::new(123)
            }),
        ]
    }

    /// A census with the named edits applied, in the order given.
    fn moved(edits: &[OneValueMoved], which: &[usize]) -> RecordingTerms {
        let mut selection = a_runs_selection_terms();
        let mut terms = a_censuss_recording_terms();
        for &at in which {
            (edits[at].1)(&mut selection, &mut terms);
        }
        terms.selection = SelectionTermsDigest::of(&selection);
        terms
    }

    /// Which terms' digests differ between two identities, by name.
    ///
    /// Asserts as it walks that the two carry the same names in the same order, so a caller
    /// cannot read a permutation as a difference.
    fn terms_that_differ<'a>(mine: &'a CensusIdentity, theirs: &CensusIdentity) -> Vec<&'a str> {
        assert_eq!(mine.terms.len(), theirs.terms.len());
        mine.terms
            .iter()
            .zip(&theirs.terms)
            .inspect(|(before, after)| assert_eq!(before.term, after.term))
            .filter(|(before, after)| before.digest != after.digest)
            .map(|(before, _)| before.term.as_str())
            .collect()
    }

    #[test]
    fn a_digest_is_thirty_two_characters_of_lower_case_hex() {
        for term in CensusIdentity::of(&a_censuss_recording_terms()).terms {
            assert_eq!(
                term.digest.len(),
                32,
                "{} was digested to {:?}",
                term.term,
                term.digest
            );
            assert!(
                term.digest
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
                "{} was digested to {:?}",
                term.term,
                term.digest
            );
        }
    }

    /// **Every value the census can name a disagreement on is a term this identity carries, and
    /// moving that value moves that term's digest and no other.**
    ///
    /// One edit at a time, and three things asserted together: the census itself reports the
    /// disagreement under some name, this identity carries a term of exactly that name, and it
    /// is the **only** term whose digest moved. So a name typed wrongly here fails, and so does
    /// a value this identity forgets to digest — which is the failure that would let a file
    /// fitted under other terms be read back as this run's own.
    #[test]
    fn every_term_is_named_as_the_census_names_it() {
        let edits = the_twelve_edits();
        assert_eq!(
            edits.len(),
            12,
            "the census refuses to pool across twelve values"
        );

        let mine = a_censuss_recording_terms();
        let identity = CensusIdentity::of(&mine);
        for (at, (what_moved, _)) in edits.iter().enumerate() {
            let theirs = moved(&edits, &[at]);
            let named = mine.first_disagreement(&theirs).unwrap_or_else(|| {
                panic!("the census reports no disagreement after {what_moved} moved")
            });
            let moved_terms = terms_that_differ(&identity, &CensusIdentity::of(&theirs));
            assert_eq!(
                moved_terms,
                vec![named],
                "after {what_moved} moved, the census names {named:?} and the identity moved \
                 {moved_terms:?}"
            );
        }
    }

    /// **The twelve are written in the order the census checks them**, which one edit at a time
    /// cannot see: any permutation passes the test above.
    ///
    /// `first_disagreement` reports the *first* value two censuses differ on, so with two moved
    /// it names the earlier of the two in its own order. Every one of the 66 pairs is tried, and
    /// the identity's own earlier-moved term must be the one the census names — which holds only
    /// if the two orders are the same one.
    ///
    /// **This is what decides which term a run reports** when a re-fitted census has drifted in
    /// more than one place, and it is the only test that can fail on an order.
    #[test]
    fn the_terms_are_in_the_order_the_census_checks_them() {
        let edits = the_twelve_edits();
        let mine = a_censuss_recording_terms();
        let identity = CensusIdentity::of(&mine);

        let mut pairs = 0;
        for earlier in 0..edits.len() {
            for later in (earlier + 1)..edits.len() {
                let theirs = moved(&edits, &[earlier, later]);
                let named = mine
                    .first_disagreement(&theirs)
                    .expect("two values moved and the census sees neither");
                let moved_terms = terms_that_differ(&identity, &CensusIdentity::of(&theirs));
                assert_eq!(
                    moved_terms.len(),
                    2,
                    "moving {:?} and {:?} moved {moved_terms:?}",
                    edits[earlier].0,
                    edits[later].0
                );
                assert_eq!(
                    moved_terms[0], named,
                    "moving {:?} and {:?}, the census names {named:?} and this identity has \
                     {:?} first",
                    edits[earlier].0, edits[later].0, moved_terms[0]
                );
                pairs += 1;
            }
        }
        assert_eq!(pairs, 66, "twelve values make sixty-six pairs");
    }

    /// **Two of the twelve values are not scalars, and every part of each one has to reach the
    /// digest** — which moving the value as a whole cannot show.
    ///
    /// The kept loci are a whole digest and a list of blocks, each block a contig, a megabase and
    /// a digest; the per-stratum counts are a list of (period, reference repeats) keys and their
    /// counts. `every_term_is_named_as_the_census_names_it` moves *one* part of each, so a
    /// version of [`CensusIdentity::of`] that dropped the whole digest, or a block's contig, or a
    /// stratum's period, passed it — **seven such mutants survived the suite before this test
    /// existed**. The census's own equality compares all of them, so a part that does not reach
    /// the digest is a pair of censuses the fit refuses to pool and this file calls the same:
    /// two runs whose kept loci differ only in which chromosome a megabase sits on would mint
    /// byte-identical identities, and nothing would demote.
    ///
    /// Each part is moved alone, and two things asserted: the census still names the term that
    /// *contains* the part, and that term is the only one whose digest moved.
    #[test]
    fn every_part_of_a_composite_value_reaches_its_term() {
        /// The term the part belongs to, what was moved, and the move.
        type OnePartMoved = (&'static str, &'static str, fn(&mut RecordingTerms));

        let parts: Vec<OnePartMoved> = vec![
            (
                "the loci actually kept",
                "the digest over every locus, with the blocks held still",
                |terms| {
                    terms.kept_loci =
                        CensusLociDigest::from_parts([4; 16], the_blocks_a_walk_kept());
                },
            ),
            (
                "the loci actually kept",
                "which contig a block covers, with its megabase and digest held still",
                |terms| {
                    let mut blocks = the_blocks_a_walk_kept();
                    blocks[0].contig = ContigId(9);
                    terms.kept_loci =
                        CensusLociDigest::from_parts(THE_WHOLE_A_WALK_DIGESTED, blocks);
                },
            ),
            ("the loci actually kept", "a block's own digest", |terms| {
                let mut blocks = the_blocks_a_walk_kept();
                blocks[0].digest ^= 1;
                terms.kept_loci = CensusLociDigest::from_parts(THE_WHOLE_A_WALK_DIGESTED, blocks);
            }),
            (
                "the loci actually kept",
                "one more block, with every existing one held still",
                |terms| {
                    let mut blocks = the_blocks_a_walk_kept();
                    blocks.push(BlockDigest {
                        contig: ContigId(3),
                        megabase: 0,
                        digest: 0x2122_2324_2526_2728,
                    });
                    terms.kept_loci =
                        CensusLociDigest::from_parts(THE_WHOLE_A_WALK_DIGESTED, blocks);
                },
            ),
            (
                // **The moved period keeps the stratum where it was in the sorted order**, and
                // that is the whole of this case. `iter_sorted` sorts on (period, repeats), so
                // a period moved from 2 to 9 sends this stratum past its neighbour and the
                // bytes after it move whether the period itself is digested or not — a version
                // that dropped the period passed such a test. At period 1 the order is what it
                // was, every other byte is what it was, and only the period differs.
                "per-stratum locus counts",
                "a stratum's motif period, with its repeats, count and sorted place held still",
                |terms| {
                    terms.ssr_stratum_counts =
                        StratumCounts::from_counted([((1, 6), 3), ((3, 5), 7)]);
                },
            ),
            (
                "per-stratum locus counts",
                "a stratum's reference repeat count, with its period and count held still",
                |terms| {
                    terms.ssr_stratum_counts =
                        StratumCounts::from_counted([((2, 9), 3), ((3, 5), 7)]);
                },
            ),
            (
                "per-stratum locus counts",
                "how many loci one stratum holds, with both its keys held still",
                |terms| {
                    terms.ssr_stratum_counts =
                        StratumCounts::from_counted([((2, 6), 4), ((3, 5), 7)]);
                },
            ),
            (
                "per-stratum locus counts",
                "one more stratum, with the grand total held still",
                |terms| {
                    terms.ssr_stratum_counts =
                        StratumCounts::from_counted([((2, 6), 3), ((3, 5), 4), ((4, 4), 3)]);
                },
            ),
        ];

        let mine = a_censuss_recording_terms();
        let identity = CensusIdentity::of(&mine);
        for (term, what_moved, move_it) in parts {
            let mut theirs = a_censuss_recording_terms();
            move_it(&mut theirs);
            let named = mine.first_disagreement(&theirs).unwrap_or_else(|| {
                panic!("the census reports no disagreement after {what_moved} moved")
            });
            assert_eq!(
                named, term,
                "moving {what_moved}, the census names {named:?} rather than {term:?}"
            );
            let moved_terms = terms_that_differ(&identity, &CensusIdentity::of(&theirs));
            assert_eq!(
                moved_terms,
                vec![term],
                "after {what_moved} moved, the identity moved {moved_terms:?}"
            );
        }
    }

    /// The twelve are twelve, and none of them is written twice.
    #[test]
    fn the_identity_names_each_value_once() {
        let identity = CensusIdentity::of(&a_censuss_recording_terms());
        let mut names: Vec<&str> = identity.terms.iter().map(|t| t.term.as_str()).collect();
        assert_eq!(names.len(), 12, "{names:?}");
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), 12, "a term is written twice: {names:?}");
    }

    /// **A digest is the bytes and not their rendering** — the two ends of the byte range, in
    /// order, so a formatter that dropped a leading zero or printed upper case fails.
    #[test]
    fn a_digest_spells_every_byte_as_two_lower_case_characters() {
        assert_eq!(
            hex_digest(&[
                0x00, 0x0f, 0xa0, 0xff, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
            ]),
            "000fa0ff0102030405060708090a0b0c"
        );
    }
}
