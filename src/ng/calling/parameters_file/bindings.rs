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

use std::collections::BTreeMap;

use super::{CensusIdentity, CensusTerm, ParametersFile, ParametersFileError, ReadGroupRow};
use crate::ng::parameter_estimation::joint::census::RecordingTerms;
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::read::input::read_groups::{ReadGroup, ReadGroups};

#[cfg(test)]
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
#[cfg(test)]
use crate::ng::parameter_estimation::joint::census::{
    DepthCap, DepthLadderDigest, ReadCap, SelectionTermsDigest,
};
#[cfg(test)]
use crate::ng::parameter_estimation::joint::loci::{
    BlockDigest, CatalogBuildSettings, CensusLociDigest, RegionSetDigest, SelectionTerms,
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
// The three bindings that refuse (spec §6)
// ---------------------------------------------------------------------

impl ParametersFile {
    /// **Refuse a file whose numbers were fitted from inputs that are not this run's.**
    ///
    /// Spec §6's first three bindings: the reference's content digest, the sample list in order
    /// by name, and the read-group table. Each refuses, because a file that fails one of them
    /// cannot be *interpreted* against this run — a file fitted against another assembly has its
    /// repeat strata cut at other tract lengths, one listing other samples has its inbreeding
    /// coefficients against other plants, and one whose read-group table does not cover the
    /// run's leaves a library with no calibration and no contamination row, which surfaces as a
    /// panic at whichever locus first carries one of that library's reads.
    ///
    /// **The fourth binding — the census — is not here**, because it does not refuse: a file
    /// fitted from a different census of this same cohort is still interpretable, merely less
    /// warranted, and what a run does with one is step D3's.
    ///
    /// # It is handed the run, not a second file
    ///
    /// The two arguments are the same two `of_run` writes from, and deliberately so: the writer
    /// and this check read one pair of inputs, so a file this run wrote is a file this run
    /// accepts. **The sample list is not among them** — it is derived from `read_groups` exactly
    /// as `of_run` derives it, because a run's samples *are* its read-group table's samples in
    /// first-seen order, and a second argument for it would be a second thing to keep in step.
    ///
    /// # The two axes are compared differently, and neither choice is free
    ///
    /// **Samples by position.** Every per-sample row of the file is read *by name* into this list
    /// and handed to calling as a *position*, so a list holding the right names in another order
    /// is not a cohort that happens to be shuffled — it is a run that gives each plant its
    /// neighbour's inbreeding coefficient and batch, with nothing downstream able to see it.
    /// `validate` refuses the same disagreement inside the file for the same reason.
    ///
    /// **Read groups by their number.** Row *order* in `fitted_from.read_groups` carries no
    /// meaning anywhere in this module — `validate` sorts the ids before checking they are dense,
    /// the projection reads the table only for its length, and every other section joins on the
    /// `read_group` key — so a file whose rows are written in another order is the same file.
    /// **Comparing these two tables positionally refused one**: two lanes of one plant swapped
    /// leaves the file's first-seen sample order unchanged, so it validates, projects, and is the
    /// file this run would have written.
    ///
    /// # What it does not check
    ///
    /// **Anything about the file's agreement with itself**, which is [`Self::validate`]'s: that
    /// the read-group ids run `0..n` with no gap, that the sample list is the file's own
    /// read-group table's first-seen order, that every keyed table covers its axis. A caller
    /// reading a file into a run wants both, and step D3 is where they are run together. **This
    /// one assumes none of them** — it is a public entry point, and a file that has not been
    /// through `validate` is compared just as soundly, because every join here is a lookup
    /// rather than a position.
    ///
    /// **⚑ Two things follow from that, and both are why the two run together.** A file whose
    /// own sample list disagrees with its own read-group table is refused here *as though the
    /// run were wrong*, where `validate` names the file precisely; and a file with no samples
    /// and no read groups matches a run with none, which `validate` refuses and this cannot see.
    /// Neither is reachable once `validate` has run first.
    ///
    /// # Errors
    ///
    /// [`ParametersFileError::FittedFromOtherInputs`], naming a key of the file and both values.
    /// **The first difference stops the walk**, in spec §6's own order — reference, then samples,
    /// then read groups. What that order decides is which refusal a run mismatched in more than
    /// one place hears about first, and the reference leads because its consequence is the worst
    /// of the three: a plausible VCF whose repeat strata were cut on another assembly, where the
    /// other two go missing loudly.
    pub fn refuse_if_not_this_runs_inputs(
        &self,
        reference: &ReferenceDigest,
        read_groups: &ReadGroups,
    ) -> Result<(), ParametersFileError> {
        let bound = &self.fitted_from;

        let this_runs_reference = hex_digest(&reference.0);
        if bound.reference_digest != this_runs_reference {
            return Err(fitted_from_other_inputs(
                "fitted_from.reference_digest",
                &bound.reference_digest,
                &this_runs_reference,
            ));
        }

        // **The run's samples are its read-group table's, in first-seen order**, which is what
        // `of_run` writes and what `validate` holds the file's own list to.
        let this_runs_samples: Vec<&str> = read_groups
            .read_groups_per_sample()
            .iter()
            .map(|of_sample| of_sample.sample.as_ref())
            .collect();
        for (at, (in_the_file, in_the_run)) in bound
            .samples
            .iter()
            .zip(&this_runs_samples)
            .enumerate()
            .map(|(at, (mine, theirs))| (at, (mine.as_str(), *theirs)))
        {
            if in_the_file != in_the_run {
                return Err(fitted_from_other_inputs(
                    format!("fitted_from.samples[{at}]"),
                    format!("{in_the_file:?}"),
                    format!("{in_the_run:?}"),
                ));
            }
        }
        // **`zip` stops at the shorter, so this is the other half of the walk above** and not a
        // tidier restatement of it: a cohort that gained or lost a plant agrees on every position
        // the two lists share. **It names the plants rather than counting them**, which is what
        // spec §6 asks for — a reader told *2 against 3* has to diff two lists by eye, and one of
        // them is not written down anywhere they can look.
        if bound.samples.len() != this_runs_samples.len() {
            return Err(fitted_from_other_inputs(
                "fitted_from.samples",
                a_list_of(
                    "samples",
                    bound.samples.iter().map(|name| format!("{name:?}")),
                ),
                a_list_of(
                    "samples",
                    this_runs_samples.iter().map(|name| format!("{name:?}")),
                ),
            ));
        }

        // **Joined on the read group's own number, never on the row's place in the table.** The
        // ids the file carries are `0..n` once `validate` has passed and the run's are `0..n` by
        // construction, so equal id *sets* is the whole of spec §6's gap check — and the join is
        // still sound on a file `validate` has not seen, since it looks each id up rather than
        // assuming any of them.
        let this_runs_lanes: BTreeMap<u32, &ReadGroup> = read_groups
            .iter()
            .map(|(id, declared)| (id.get(), declared))
            .collect();
        let the_files_lanes: BTreeMap<u32, &ReadGroupRow> = bound
            .read_groups
            .iter()
            .map(|row| (row.read_group, row))
            .collect();
        if the_files_lanes.keys().ne(this_runs_lanes.keys()) {
            // **Each lane is its number *and* its `@RG ID`.** The `@RG ID` is what the reader's
            // alignment files call it and the number is what every other section of this file
            // joins on, and either alone gives a message that reads as a contradiction: a file
            // numbered `0, 1, 3` for this run's three lanes has the same three names as the run,
            // so names alone print two identical lists beside "these differ".
            return Err(fitted_from_other_inputs(
                "fitted_from.read_groups",
                a_list_of(
                    "read groups",
                    the_files_lanes
                        .values()
                        .map(|row| format!("{} {:?}", row.read_group, row.declared_id)),
                ),
                a_list_of(
                    "read groups",
                    this_runs_lanes
                        .iter()
                        .map(|(number, lane)| format!("{number} {:?}", lane.id)),
                ),
            ));
        }
        for (number, (row, declared)) in the_files_lanes
            .values()
            .zip(this_runs_lanes.values())
            .map(|(row, lane)| (row.read_group, (*row, *lane)))
        {
            // **The three names a reader joins a row to a lane by.** The number itself is what
            // the two tables were joined on, so it cannot differ here.
            for (key, in_the_file, in_the_run) in [
                (
                    "declared_id",
                    format!("{:?}", row.declared_id),
                    format!("{:?}", declared.id),
                ),
                (
                    "library",
                    format!("{:?}", row.library),
                    format!("{:?}", declared.library.value),
                ),
                (
                    "sample",
                    format!("{:?}", row.sample),
                    format!("{:?}", declared.sample),
                ),
            ] {
                if in_the_file != in_the_run {
                    return Err(fitted_from_other_inputs(
                        format!("fitted_from.read_groups[read_group = {number}].{key}"),
                        in_the_file,
                        in_the_run,
                    ));
                }
            }
        }

        Ok(())
    }
}

/// How many names a refusal spells out before it says how many more there are.
///
/// **Spec §9 commits this file to 3,000 samples**, and a refusal that printed three thousand
/// names is one nobody reads. Five is enough to see which cohort is which and short enough to
/// sit on one line beside its twin.
const NAMES_BEFORE_A_TALLY: usize = 5;

/// A list as a refusal spells it: how many there are, what they are, then the first few.
///
/// **The count first**, because the two lists a refusal prints differ in length by construction
/// and that is the fact a reader acts on; the entries are what says *which*. **And the noun**,
/// because the two lists sit inside one sentence and a bare number twice over is a sentence a
/// reader has to parse rather than read.
fn a_list_of(what: &str, entries: impl Iterator<Item = String>) -> String {
    let entries: Vec<String> = entries.collect();
    let shown = entries
        .iter()
        .take(NAMES_BEFORE_A_TALLY)
        .cloned()
        .collect::<Vec<String>>()
        .join(", ");
    let and_more = entries.len().saturating_sub(NAMES_BEFORE_A_TALLY);
    let tail = if and_more == 0 {
        String::new()
    } else {
        format!(" and {and_more} more")
    };
    format!("{} {what} ({shown}{tail})", entries.len())
}

/// One binding that does not hold, named with both values.
fn fitted_from_other_inputs(
    field: impl Into<String>,
    in_the_file: impl Into<String>,
    in_the_run: impl Into<String>,
) -> ParametersFileError {
    ParametersFileError::FittedFromOtherInputs {
        field: field.into(),
        in_the_file: in_the_file.into(),
        in_the_run: in_the_run.into(),
    }
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
    use super::super::tests::a_file_using_every_shape;
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

    // -----------------------------------------------------------------
    // The three refusals (spec §6, §13 test 4)
    // -----------------------------------------------------------------

    /// The run the module's shared file fixture was fitted by: its reference, and its lanes.
    fn the_run_the_fixture_was_fitted_by() -> (ReferenceDigest, ReadGroups) {
        (
            super::super::tests::THE_REFERENCE_A_RUN_FITTED_AGAINST,
            ReadGroups::of_lanes(&[
                ("HWI.3", "TS-1", "lib3"),
                ("HWI.4", "TS-1", "lib4"),
                ("HWI.5", "Ailsa ‘Craig’ \"×2\"", "lib5"),
            ]),
        )
    }

    /// The three parts of the refusal a run gets, or a panic naming what it got instead.
    fn refused(
        file: &ParametersFile,
        reference: &ReferenceDigest,
        read_groups: &ReadGroups,
    ) -> (String, String, String) {
        match file.refuse_if_not_this_runs_inputs(reference, read_groups) {
            Err(ParametersFileError::FittedFromOtherInputs {
                field,
                in_the_file,
                in_the_run,
            }) => (field, in_the_file, in_the_run),
            Err(other) => panic!("refused, but as {other}"),
            Ok(()) => panic!("accepted"),
        }
    }

    /// **The file a run wrote is the file that run accepts** — the case every refusal below is a
    /// departure from, and the one a green suite would still pass without.
    #[test]
    fn the_run_that_wrote_the_file_is_not_refused() {
        let (reference, read_groups) = the_run_the_fixture_was_fitted_by();
        a_file_using_every_shape()
            .refuse_if_not_this_runs_inputs(&reference, &read_groups)
            .expect("the fixture is the file this run's inputs would write");
    }

    /// **Spec §6's first binding.** A file fitted against another assembly gives a plausible VCF
    /// with its repeat strata cut at other tract lengths, so it is refused rather than demoted.
    #[test]
    fn a_file_fitted_against_another_reference_is_refused() {
        let (_, read_groups) = the_run_the_fixture_was_fitted_by();
        let another = ReferenceDigest([0xab; 16]);
        let (field, in_the_file, in_the_run) =
            refused(&a_file_using_every_shape(), &another, &read_groups);

        assert_eq!(field, "fitted_from.reference_digest");
        // **Both values, which is what the ruling of 2026-08-30 asks for.** A reader with three
        // copies of one assembly on disk cannot act on a field name alone.
        assert_eq!(in_the_file, "0123456789abcdef0123456789abcdef");
        assert_eq!(in_the_run, "ab".repeat(16));
    }

    /// **Spec §6's second binding**, and it names the samples.
    #[test]
    fn a_file_listing_samples_the_run_does_not_have_is_refused() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let another_cohort = ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", "TS-9", "lib5"),
        ]);
        let (field, in_the_file, in_the_run) =
            refused(&a_file_using_every_shape(), &reference, &another_cohort);

        assert_eq!(field, "fitted_from.samples[1]");
        assert_eq!(in_the_file, "\"Ailsa ‘Craig’ \\\"×2\\\"\"");
        assert_eq!(in_the_run, "\"TS-9\"");
    }

    /// **A cohort of the same plants in another order is refused too**, and that is the point of
    /// comparing by position: every per-sample row is read by name into this list and handed to
    /// calling as a position, so a run whose order differs gives each plant its neighbour's
    /// inbreeding coefficient. Nothing downstream could see it.
    #[test]
    fn the_same_samples_in_another_order_are_refused() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let reordered = ReadGroups::of_lanes(&[
            ("HWI.5", "Ailsa ‘Craig’ \"×2\"", "lib5"),
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
        ]);
        let (field, in_the_file, in_the_run) =
            refused(&a_file_using_every_shape(), &reference, &reordered);
        assert_eq!(field, "fitted_from.samples[0]");
        assert_eq!(in_the_file, "\"TS-1\"");
        assert_eq!(in_the_run, "\"Ailsa ‘Craig’ \\\"×2\\\"\"");
    }

    /// A run that has a plant the file never saw, with the file's own list a prefix of it.
    #[test]
    fn a_cohort_with_one_more_sample_is_refused_by_its_count() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let one_more = ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", "Ailsa ‘Craig’ \"×2\"", "lib5"),
            ("HWI.6", "TS-4", "lib6"),
        ]);
        let (field, in_the_file, in_the_run) =
            refused(&a_file_using_every_shape(), &reference, &one_more);

        // **It names the plants, which is what spec §6 asks for** — a reader told *2 against
        // 3* has to diff two lists by eye, and the run's is not written down anywhere.
        assert_eq!(field, "fitted_from.samples");
        assert_eq!(
            in_the_file,
            "2 samples (\"TS-1\", \"Ailsa ‘Craig’ \\\"×2\\\"\")"
        );
        assert_eq!(
            in_the_run,
            "3 samples (\"TS-1\", \"Ailsa ‘Craig’ \\\"×2\\\"\", \"TS-4\")"
        );
    }

    /// **Spec §6's third binding.** A run with a library the file's table does not cover is the
    /// gap §6 names: that library's calibration and contamination row are simply absent, and the
    /// symptom without this check is a panic at whichever locus first carries one of its reads.
    #[test]
    fn a_run_whose_read_groups_the_file_does_not_cover_is_refused() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let a_fourth_lane = ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", "Ailsa ‘Craig’ \"×2\"", "lib5"),
            ("HWI.6", "Ailsa ‘Craig’ \"×2\"", "lib6"),
        ]);
        let (field, in_the_file, in_the_run) =
            refused(&a_file_using_every_shape(), &reference, &a_fourth_lane);

        // The samples still agree — both cohorts are the same two plants — so the walk reaches
        // the read-group table, which is what this test is about. **The lane is named by its
        // `@RG ID`**: the run's fourth is `HWI.6`, and the file has no row for it.
        assert_eq!(field, "fitted_from.read_groups");
        assert_eq!(
            in_the_file,
            "3 read groups (0 \"HWI.3\", 1 \"HWI.4\", 2 \"HWI.5\")"
        );
        assert_eq!(
            in_the_run,
            "4 read groups (0 \"HWI.3\", 1 \"HWI.4\", 2 \"HWI.5\", 3 \"HWI.6\")"
        );
    }

    /// **Each of the three names on a read-group row is compared, and both values are
    /// asserted.**
    ///
    /// The three are what a reader joins a row to a lane by. **The values and not only the
    /// field**: a review swapped `in_the_file` and `in_the_run` in the library arm and the whole
    /// module stayed green, and the message that mutant produces tells a geneticist the file
    /// says `lib9` where the file says `lib4`.
    #[test]
    fn a_read_group_row_that_differs_in_any_name_is_refused() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let awkward = "Ailsa ‘Craig’ \"×2\"";
        for (key, in_the_files_row, in_the_runs_lane, lanes) in [
            (
                "declared_id",
                "\"HWI.4\"",
                "\"HWI.9\"",
                [
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.9", "TS-1", "lib4"),
                    ("HWI.5", awkward, "lib5"),
                ],
            ),
            (
                "library",
                "\"lib4\"",
                "\"lib9\"",
                [
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.4", "TS-1", "lib9"),
                    ("HWI.5", awkward, "lib5"),
                ],
            ),
            (
                // The plant a lane belongs to, moved without moving the sample *list*: both
                // cohorts still hold these two plants in this order, and only which lane is
                // whose has changed. Nothing before this check can see it.
                "sample",
                "\"TS-1\"",
                "\"Ailsa ‘Craig’ \\\"×2\\\"\"",
                [
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.4", awkward, "lib4"),
                    ("HWI.5", awkward, "lib5"),
                ],
            ),
        ] {
            let (field, in_the_file, in_the_run) = refused(
                &a_file_using_every_shape(),
                &reference,
                &ReadGroups::of_lanes(&lanes),
            );
            assert_eq!(
                field,
                format!("fitted_from.read_groups[read_group = 1].{key}")
            );
            assert_eq!(in_the_file, in_the_files_row);
            assert_eq!(in_the_run, in_the_runs_lane);
        }
    }

    /// **A file whose rows are numbered otherwise is a file for other read groups**, however
    /// well its names match: every other section of the file joins on `read_group`, so a table
    /// numbered `0, 1, 3` files this run's lanes under numbers the run does not use.
    #[test]
    fn a_read_group_table_numbered_otherwise_is_refused() {
        let (reference, read_groups) = the_run_the_fixture_was_fitted_by();
        let mut file = a_file_using_every_shape();
        file.fitted_from.read_groups[2].read_group = 3;
        let (field, in_the_file, in_the_run) = refused(&file, &reference, &read_groups);

        // **The number and the `@RG ID` together**, which is what keeps this message from
        // printing two identical lists: the three lanes are the same three lanes, and the file
        // has filed the last of them under a number the run does not use.
        assert_eq!(field, "fitted_from.read_groups");
        assert_eq!(
            in_the_file,
            "3 read groups (0 \"HWI.3\", 1 \"HWI.4\", 3 \"HWI.5\")"
        );
        assert_eq!(
            in_the_run,
            "3 read groups (0 \"HWI.3\", 1 \"HWI.4\", 2 \"HWI.5\")"
        );
    }

    /// **⚑ A file whose rows are written in another order is the same file, and is accepted.**
    ///
    /// Row order in `fitted_from.read_groups` carries no meaning anywhere in this module:
    /// `validate` sorts the ids before checking they are dense, the projection reads the table
    /// only for its length, and every other section joins on the `read_group` key. **The first
    /// draft of this check joined the two tables positionally and refused such a file** — two
    /// lanes of one plant swapped, which leaves the file's own first-seen sample order unchanged,
    /// so it validates and projects and is the file this run would have written.
    #[test]
    fn a_file_whose_rows_are_written_in_another_order_is_the_same_file() {
        let (reference, read_groups) = the_run_the_fixture_was_fitted_by();
        let mut file = a_file_using_every_shape();
        file.fitted_from.read_groups.swap(0, 1);

        file.validate()
            .expect("row order is not the file's meaning");
        file.refuse_if_not_this_runs_inputs(&reference, &read_groups)
            .expect("and it is still this run's file");
    }

    /// **The three are checked in spec §6's own order**, which is what decides which refusal a
    /// run mismatched in more than one place hears about.
    ///
    /// The reference leads because its consequence is the worst of the three — a plausible VCF
    /// whose repeat strata were cut on another assembly — where the other two go missing loudly.
    #[test]
    fn the_reference_is_the_first_thing_checked() {
        let another_cohort = ReadGroups::of_lanes(&[("HWI.9", "TS-9", "lib9")]);
        let (field, ..) = refused(
            &a_file_using_every_shape(),
            &ReferenceDigest([0xab; 16]),
            &another_cohort,
        );
        assert_eq!(field, "fitted_from.reference_digest");
    }

    /// **Every binding refusal names a key the file actually contains** — the same guarantee
    /// `validate`'s own refusals carry, so a reader meets one vocabulary and can find what they
    /// are told about by searching the file they have in front of them.
    ///
    /// **The alternative was prose**, in the shape `Freshness`'s `"the pileup's header"` uses one
    /// level down. It was dropped because two of the three field names it produced —
    /// *the sample in position 3*, *the library of read group 2* — appear nowhere in a produced
    /// file, so the promise their own doc comment made was false for two cases in three.
    #[test]
    fn every_refusal_names_a_key_the_file_contains() {
        let file = a_file_using_every_shape();
        let text = file.to_toml();
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let awkward = "Ailsa ‘Craig’ \"×2\"";

        let elsewhere: Vec<(ReferenceDigest, ReadGroups)> = vec![
            (ReferenceDigest([0xab; 16]), ReadGroups::of_lanes(&[])),
            (
                reference,
                ReadGroups::of_lanes(&[
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.4", "TS-1", "lib4"),
                    ("HWI.5", "TS-9", "lib5"),
                ]),
            ),
            (
                reference,
                ReadGroups::of_lanes(&[
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.4", "TS-1", "lib4"),
                    ("HWI.5", awkward, "lib5"),
                    ("HWI.6", awkward, "lib6"),
                ]),
            ),
            (
                reference,
                ReadGroups::of_lanes(&[
                    ("HWI.3", "TS-1", "lib3"),
                    ("HWI.9", "TS-1", "lib4"),
                    ("HWI.5", awkward, "lib5"),
                ]),
            ),
        ];
        assert_eq!(
            elsewhere.len(),
            4,
            "one run for each shape of refusal this check can raise"
        );
        for (reference, lanes) in &elsewhere {
            let (field, ..) = refused(&file, reference, lanes);
            // **Every segment, not just the last**, and a segment carrying a row key is checked
            // as its key name — the same walk `validate`'s own version of this test makes.
            for segment in field.split('.') {
                let key = segment.split(['[', '(']).next().unwrap_or(segment);
                assert!(
                    !key.is_empty() && text.contains(key),
                    "the refusal names {field}, whose segment {key:?} is not a key in the file"
                );
            }
        }
    }

    /// **A cohort of 3,000 is the top of the committed range, and a refusal is still one line.**
    ///
    /// Spec §9 prices this file at 3,000 samples; a refusal that spelled three thousand names is
    /// one nobody reads. The names are what says *which* cohort and the count is what a reader
    /// acts on, so both are kept and the list is cut.
    #[test]
    fn a_refusal_over_a_cohort_of_thousands_still_fits_on_a_line() {
        let (reference, _) = the_run_the_fixture_was_fitted_by();
        let named: Vec<(String, String, String)> = (0..3_000)
            .map(|plant| {
                (
                    format!("HWI.{plant}"),
                    format!("TS-{plant}"),
                    format!("lib{plant}"),
                )
            })
            .collect();
        let lanes: Vec<(&str, &str, &str)> = named
            .iter()
            .map(|(id, sample, library)| (id.as_str(), sample.as_str(), library.as_str()))
            .collect();

        let (field, in_the_file, in_the_run) = refused(
            &a_file_using_every_shape(),
            &reference,
            &ReadGroups::of_lanes(&lanes),
        );

        assert_eq!(field, "fitted_from.samples[0]");
        assert_eq!(in_the_file, "\"TS-1\"");
        assert_eq!(in_the_run, "\"TS-0\"");

        // And where the two lists agree on every shared position, the tally rather than the list.
        let mut of_three_thousand = a_file_using_every_shape();
        of_three_thousand.fitted_from.samples =
            named.iter().map(|(_, sample, _)| sample.clone()).collect();
        let (field, in_the_file, in_the_run) = refused(
            &of_three_thousand,
            &reference,
            &ReadGroups::of_lanes(&lanes[..2_999]),
        );
        assert_eq!(field, "fitted_from.samples");
        assert_eq!(
            in_the_file,
            "3000 samples (\"TS-0\", \"TS-1\", \"TS-2\", \"TS-3\", \"TS-4\" and 2995 more)"
        );
        assert_eq!(
            in_the_run,
            "2999 samples (\"TS-0\", \"TS-1\", \"TS-2\", \"TS-3\", \"TS-4\" and 2994 more)"
        );
        // **Asserted rather than described**, so the claim cannot go stale: a cohort at the top
        // of spec §9's committed range still names five plants and a tally, on one line.
        assert_eq!(in_the_file.len(), 67);
    }

    /// A binding refusal has no line and no parser behind it, and both accessors say so.
    #[test]
    fn a_binding_refusal_has_no_line_in_the_file_to_send_anyone_to() {
        let (_, read_groups) = the_run_the_fixture_was_fitted_by();
        let refusal = a_file_using_every_shape()
            .refuse_if_not_this_runs_inputs(&ReferenceDigest([0xab; 16]), &read_groups)
            .expect_err("another reference");

        assert_eq!(refusal.line(), None);
        assert_eq!(refusal.rendered_by_the_parser(), refusal.to_string());
        assert!(
            refusal
                .to_string()
                .contains("0123456789abcdef0123456789abcdef")
                && refusal.to_string().contains(&"ab".repeat(16)),
            "both values reach the message a run prints: {refusal}"
        );
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
