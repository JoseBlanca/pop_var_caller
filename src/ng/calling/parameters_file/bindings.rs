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

use super::{
    BaseQualityCalibrationRow, CensusIdentity, CensusTerm, InbreedingRow, ParametersFile,
    ParametersFileError, ReadGroupRow, RepeatRouting, RepeatTracts, RunParametersFromFile,
    StatedConstants, SubstitutionRateRow, Warrant,
};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::census::RecordingTerms;
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::read::input::read_groups::{ReadGroup, ReadGroups};
use crate::ng::region_typing::segment_criteria::SsrSegmentCriteria;
use crate::ng::repeat_catalog::StrRepeatCriteria;

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
use crate::ng::repeat_catalog::StratumCounts;
#[cfg(test)]
use crate::ng::tandem_repeat::ScanParams;
#[cfg(test)]
use crate::ng::types::{Bp, ContigId};

impl CensusIdentity {
    /// **The run that has no census at all** — no terms, because there was no store of evidence.
    ///
    /// **Two runs reach it and both are real** (spec §7, `run_streaming.md` §2): the defaults run,
    /// which fitted nothing, and any direct-mode run, which reads its evidence from the alignment
    /// files and builds no psp and no census. Spec §7 says such a run writes its parameters file
    /// like any other, so `[fitted_from.census]` has to be able to say *there was none*.
    ///
    /// **An empty list of terms rather than an absent section**, and the difference is what a
    /// later reader does with it. [`ParametersFile::census_disagreement`] already treats a term
    /// one identity has and the other does not as a disagreement, so a psp-mode run reading this
    /// file finds one at the first term and demotes — which is the right answer, since none of
    /// these numbers was fitted under *its* census. An absent section would have to be given that
    /// meaning separately, in a second place, and could come to disagree with this one.
    ///
    /// **The demotion costs such a file nothing**, and that is worth knowing before it looks
    /// alarming: demotion is [`weaker_of`](crate::ng::parameter_estimation::Provenance::weaker_of)
    /// against `Supplied`, and every number a run with no census wrote is already `supplied` or
    /// `defaulted`.
    #[must_use]
    pub fn of_a_run_with_no_census() -> Self {
        Self { terms: Vec::new() }
    }

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
// What the run counted as a repeat (spec §3.9)
// ---------------------------------------------------------------------

impl RepeatRouting {
    /// **What a run asked the repeat catalog for, as the file spells it.**
    ///
    /// Destructured without `..` on purpose, the same convention the census identity keeps
    /// upstream: an axis added to the criteria stops this compiling rather than quietly
    /// dropping out of the record, and an axis that drops out lets two runs that routed
    /// differently write identical files.
    #[must_use]
    pub fn of(criteria: &StrRepeatCriteria) -> Self {
        let StrRepeatCriteria {
            classification,
            min_flank_bp,
            max_str_len_bp,
        } = criteria;
        let SsrSegmentCriteria {
            periods,
            min_copies,
            min_purity,
            min_score,
            bundle_threshold,
        } = classification;
        Self {
            min_copies: std::array::from_fn(|index| {
                min_copies.for_period(u8::try_from(index + 1).expect("six periods fit a u8"))
            }),
            min_period: periods.min(),
            max_period: periods.max(),
            max_str_len: max_str_len_bp.get(),
            min_purity: *min_purity,
            min_flank_bp: min_flank_bp.get(),
            min_score: *min_score,
            bundle_threshold: *bundle_threshold,
        }
    }
}

impl ParametersFile {
    /// **Which routing threshold this file's run and `asked` first differ on** — `None` where
    /// they routed the same ground, and `None` where the file does not say.
    ///
    /// **Nothing here refuses and nothing is demoted**, which is what separates this from the
    /// census comparison beside it (spec §3.9, the owner's ruling of 2026-09-02): a file fitted
    /// under another census carries numbers fitted elsewhere, where a file written by a run that
    /// routed differently carries numbers that are as warranted as they ever were. Only the
    /// ground they are applied to has moved, and the caller's job is to say so.
    ///
    /// **A file with no `[repeat_routing]` answers `None`**, because it makes no claim to
    /// disagree with — spec §5's rule that absence is not a value.
    ///
    /// The axes are compared in the order the file writes them, and the name returned is the
    /// flag a person would move, so that a caller can quote it without a second table.
    #[must_use]
    pub fn routing_disagreement(&self, asked: &StrRepeatCriteria) -> Option<&'static str> {
        let recorded = self.repeat_routing.as_ref()?;
        let mine = RepeatRouting::of(asked);
        // Exhaustive on purpose: an axis added to the record must be answered for here, or two
        // runs that differ only on it would compare equal.
        let RepeatRouting {
            min_copies,
            min_period,
            max_period,
            max_str_len,
            min_purity,
            min_flank_bp,
            min_score,
            bundle_threshold,
        } = recorded;
        if *min_copies != mine.min_copies {
            return Some("--min-copies");
        }
        if *min_period != mine.min_period {
            return Some("--min-period");
        }
        if *max_period != mine.max_period {
            return Some("--max-period");
        }
        if *max_str_len != mine.max_str_len {
            return Some("--max-str-len");
        }
        // Bit equality, not a tolerance: this is a recorded input rather than a fitted number,
        // and a run either typed the same value or typed another one.
        if min_purity.to_bits() != mine.min_purity.to_bits() {
            return Some("--min-purity");
        }
        if *min_flank_bp != mine.min_flank_bp {
            return Some("the flank floor");
        }
        if *min_score != mine.min_score {
            return Some("the scanner score floor");
        }
        if *bundle_threshold != mine.bundle_threshold {
            return Some("the bundling distance");
        }
        None
    }
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

// ---------------------------------------------------------------------
// The fourth binding, which demotes (spec §2.1, §6, §13 test 5)
// ---------------------------------------------------------------------

impl ParametersFile {
    /// **Read this file into the parameters *this run* scores with.**
    ///
    /// The whole of spec §6 in one door, in the order that gives the better message:
    ///
    /// 1. [`Self::validate`] — is this a parameters file that means anything at all? A file whose
    ///    own sample list disagrees with its own read-group table is the file's fault, and
    ///    running the bindings first would report it as the *run's*.
    /// 2. [`Self::refuse_if_not_this_runs_inputs`] — the three bindings that refuse.
    /// 3. the census, which **demotes rather than refusing**: the numbers are still numbers, and
    ///    a file fitted from another census of this same cohort is interpretable, merely less
    ///    warranted.
    /// 4. [`Self::to_run_parameters`], on the file as demotion left it.
    ///
    /// # `census: None` — the run that has no census to compare against
    ///
    /// **Direct mode has none** (`run_streaming.md` §2): it reads its evidence from the alignment
    /// files, builds no psp and runs no fit, and it is *the* mode this file format exists for —
    /// the parameters file is that mode's user-facing input. So the argument is an `Option`, and
    /// the question is what `None` should do with the fourth binding.
    ///
    /// **`None` keeps the file's warrants**, and spec §2.1 settles it rather than leaving it to
    /// taste. Demoting on every read was considered there and rejected *because it breaks the
    /// two-mode oracle*: the same cohort called in direct mode from a file and in psp mode from
    /// the fit in memory must report the same warrants for identical genotypes. Direct mode is
    /// exactly the mode with no census — so demoting whenever there is nothing to compare
    /// against **is** demoting on every read, under another name, and would break that oracle at
    /// every locus.
    ///
    /// **What `None` gives up is said plainly**: a file fitted under another census of this
    /// cohort reads back into a direct-mode run with its `fitted_here` warrants intact, where the
    /// same file in psp mode would be demoted. That is a difference in what the run *reports*
    /// and never in what it *computes* — §2 is explicit that consumers combine warrants and do
    /// not branch on them — and the alternative trades it for a difference in what two modes
    /// report about the same call, which is the larger of the two.
    ///
    /// # Why the demotion happens to the *file* and not to the parameters
    ///
    /// Spec §2.1's trap for the coder is that **demotion is per-file, not per-number** — there is
    /// no state in which the binding leaves some of a file's numbers fitted and others not. The
    /// shortest way to be that is to demote the file, whose warrants are five public fields, and
    /// then project it once. Demoting afterwards would mean reaching into `RunParameters`,
    /// `StratumFits` and `Estimate` separately and getting all five right, which is five chances
    /// to leave one behind.
    ///
    /// **⚑ Validating first is load-bearing for more than the message.**
    /// [`Self::demoted_to_no_better_than_supplied`] is public and moves warrants, and one of the
    /// warrants it moves is one `validate` has an opinion about: a `fitted_here` outlier weight
    /// is a state no run can mean and `validate` refuses it, but demotion turns it into a
    /// `supplied` one, which is legal. So demoting an unvalidated file *launders* it —
    /// `file.demoted_to_no_better_than_supplied().to_run_parameters()` accepts what
    /// `file.to_run_parameters()` turns down. This door cannot, because nothing is demoted until
    /// `validate` has passed.
    ///
    /// # Errors
    ///
    /// Whatever `validate`, the bindings or the projection refuse — in that order.
    pub fn to_run_parameters_for(
        &self,
        reference: &ReferenceDigest,
        read_groups: &ReadGroups,
        census: Option<&CensusIdentity>,
    ) -> Result<ParametersForThisRun, ParametersFileError> {
        self.validate()?;
        self.refuse_if_not_this_runs_inputs(reference, read_groups)?;

        // **An identity naming no terms is the same fact as no identity at all**, and both
        // reach here: the writer spells *this run had no census* as
        // `CensusIdentity::of_a_run_with_no_census()` and the reader spells it `None`. A driver
        // holding one identity — the natural shape, since `of_run` takes one by value — would
        // otherwise hand this the writer's spelling and get the demotion the argument above
        // exists to avoid.
        let census = census.filter(|census| !census.terms.is_empty());
        let agreement = match census.map(|census| self.census_disagreement(census)) {
            None => CensusAgreement::NothingToCompareAgainst,
            Some(None) => CensusAgreement::TheSameCensus,
            Some(Some(term)) => CensusAgreement::FittedUnderAnother(term),
        };
        // **Two costs here, both once per run and both deliberate.** `validate` runs a second
        // time inside the projection, on a file that has already passed it — both are public
        // entry points and neither may assume the other ran. And the demotion copies the file,
        // whose largest axis spec §9 prices at up to 62 MB at 3,000 samples; the copy is
        // transient and dropped as soon as the projection is done, and demoting the file rather
        // than the assembled parameters is what makes the demotion provably per-file.
        let from_file = if agreement.demoted_the_file() {
            self.demoted_to_no_better_than_supplied()
                .to_run_parameters()?
        } else {
            self.to_run_parameters()?
        };
        Ok(ParametersForThisRun {
            from_file,
            census: agreement,
        })
    }

    /// **Which recording term the file's census and this run's first differ on** — `None` where
    /// they are the same census.
    ///
    /// The same question, and the same answer, as
    /// [`RecordingTerms::first_disagreement`](crate::ng::parameter_estimation::joint::census::RecordingTerms::first_disagreement)
    /// asks of two samples one level down, which is why [`CensusIdentity::of`] mints its terms in
    /// that function's own words and order.
    ///
    /// **A term one identity has and the other does not is a disagreement too**, and it is the
    /// shape a file written by another build produces: the census can grow a thirteenth value,
    /// and a file that names twelve was fitted under terms this build cannot even compare.
    /// Whichever list is longer names it, since only that one has a name to give.
    #[must_use]
    pub fn census_disagreement(&self, census: &CensusIdentity) -> Option<String> {
        for (mine, theirs) in self.fitted_from.census.terms.iter().zip(&census.terms) {
            if mine.term != theirs.term || mine.digest != theirs.digest {
                return Some(mine.term.clone());
            }
        }
        // **The shorter list runs out first**, so whichever identity is longer names the term.
        // The file's own is preferred where it has one, because that is the value a reader can
        // look at.
        self.fitted_from
            .census
            .terms
            .get(census.terms.len())
            .or_else(|| census.terms.get(self.fitted_from.census.terms.len()))
            .map(|extra| extra.term.clone())
    }

    /// **Every number in the file, no better warranted than `supplied`** — spec §2.1's demotion.
    ///
    /// **Named for what it does and not for what §2.1 says**, because the two differ by exactly
    /// the thing below: *demoted to supplied* would assert the sentence this doc exists to
    /// correct, at every call site, to a reader who never opens it.
    ///
    /// A number the fit called `fitted_here` was fitted from *some* cohort's data; read into a
    /// run over a different one it is a number somebody handed over. §2.1 keeps the warrant
    /// where the file's binding matches and demotes the whole file where it does not.
    ///
    /// # It is `weaker_of` and not an assignment, and for one number that matters
    ///
    /// `Provenance` ranks `Supplied` **above** `Defaulted` — a number the run was handed says
    /// nothing about this data, and a stated constant says less than nothing. So assigning
    /// `Supplied` would *promote* every defaulted number — every one of them a claim that
    /// somebody chose a value nobody chose. **The repeat-tract outlier weight is the strongest
    /// case rather than the only one**: it has no fitted state at all and its two legal warrants
    /// are `supplied` and `defaulted`-at-the-project's-own-0.01, so the promotion there is one
    /// `validate` would accept and no reader could see. A defaulted calibration multiplier of
    /// exactly 1.0 promoted to `supplied` says the same false thing more quietly.
    /// [`Provenance::weaker_of`] is a no-op for every already-weaker number and the right answer
    /// for the rest.
    ///
    /// **So *every warrant is `Supplied`* is not true of a demoted file, and cannot be** — spec
    /// §13's fifth test says it and the ladder says otherwise. What is true is that no warrant is
    /// stronger than `Supplied` and none was promoted.
    ///
    /// # What is not demoted, because it is not a warrant
    ///
    /// Five numbers in this file carry a `Warrant` and the demotion moves all five. **The rest
    /// carry other vocabularies, and none of those has a *handed over* state**: a slippage
    /// number says it came off the stratum's own fit, its period's curve, or a blend; the prior
    /// seed says which rung it came from; a contamination fraction says which reads it was
    /// fitted from.
    ///
    /// **⚑ So a demoted file still says *this run's own* about numbers that are not this run's,
    /// and that is a defect rather than a nicety.** `SeedRung::FittedCurve` reads "both moments
    /// came off **the run's own** fitted population curve", and after a demotion the run that
    /// fitted it is a different run. Nothing in the file says otherwise.
    ///
    /// **This is open and it is the owner's**, recorded in `PROJECT_STATUS.md`, which offers
    /// three ways out and recommends the one D3 did *not* take — refusing such a file like the
    /// other three bindings. D3 builds what the plan and §2.1 describe; if the owner takes that
    /// recommendation, this method and the door above it go.
    #[must_use]
    pub fn demoted_to_no_better_than_supplied(&self) -> Self {
        let mut demoted = self.clone();
        // **Destructured without `..` on purpose**, which is this module's convention where a
        // walk has to reach everything ([`CensusIdentity::of`] does it over `RecordingTerms`).
        // A section added to the file, or a key added to either section that holds a warranted
        // number outside a row, stops this compiling rather than quietly keeping its warrant —
        // and a number the demotion forgets is exactly the per-number exemption spec §2.1 says
        // does not exist.
        //
        // **The rows are destructured too**, for the same reason and because nothing else would
        // notice: `every_warrant_in`, in this module's tests, pushes one warrant a row, so a
        // second warranted key on a row type would leave its count unmoved.
        let Self {
            format_version: _,
            ploidy: _,
            fitted_from: _,
            base_quality_calibration,
            contamination: _,
            sequencing_batches: _,
            inbreeding,
            ordinary_site_prior: _,
            repeat_tracts,
            stated_constants,
            // Carries no warrant: it is what the run typed, not what a fit found, so there is
            // nothing here for a demotion to weaken.
            repeat_routing: _,
        } = &mut demoted;
        for BaseQualityCalibrationRow {
            read_group: _,
            error_probability_multiplier,
        } in &mut base_quality_calibration.by_read_group
        {
            error_probability_multiplier.warrant =
                no_better_than_supplied(error_probability_multiplier.warrant);
        }
        for InbreedingRow {
            sample: _,
            inbreeding_coefficient,
        } in &mut inbreeding.by_sample
        {
            inbreeding_coefficient.warrant =
                no_better_than_supplied(inbreeding_coefficient.warrant);
        }
        let RepeatTracts {
            fallback_length_spectrum_concentration,
            slippage_group_by_read_group: _,
            slippage_by_stratum_and_group: _,
            length_spectrum_by_stratum: _,
            length_spectrum_by_period: _,
            substitution_rate_by_stratum,
        } = repeat_tracts;
        fallback_length_spectrum_concentration.warrant =
            no_better_than_supplied(fallback_length_spectrum_concentration.warrant);
        for SubstitutionRateRow {
            read_group: _,
            period: _,
            reference_repeats: _,
            ploidy: _,
            rate,
        } in substitution_rate_by_stratum
        {
            rate.warrant = no_better_than_supplied(rate.warrant);
        }
        let StatedConstants {
            repeat_tract_outlier_weight,
            repeat_tract_junk_decay_per_unit,
        } = stated_constants;
        repeat_tract_outlier_weight.warrant =
            no_better_than_supplied(repeat_tract_outlier_weight.warrant);
        repeat_tract_junk_decay_per_unit.warrant =
            no_better_than_supplied(repeat_tract_junk_decay_per_unit.warrant);
        demoted
    }
}

/// **What a run got when it read a parameters file for itself** — the numbers, and whether it had
/// to demote them.
#[derive(Debug)]
pub struct ParametersForThisRun {
    /// What calling scores with, and the two things `RunParameters` does not keep.
    ///
    /// **Named `from_file` rather than `parameters`** so that reaching the run's own parameters
    /// through it reads as `from_file.parameters` rather than as one word twice.
    pub from_file: RunParametersFromFile,
    /// **What comparing this run's census against the file's found** — including the case where
    /// there was nothing to compare against.
    pub census: CensusAgreement,
}

/// **What comparing a run's census against a file's found** — spec §6's fourth binding, which
/// demotes rather than refusing.
///
/// **Three states rather than two, and the third is why.** Until step F1 this was
/// `Option<String>`, `None` documented as *they agree*. A run with no census of its own — every
/// direct-mode run (`run_streaming.md` §2) — makes no comparison at all, and folding that into
/// `None` tells whatever reports it that the file matched this run's census when no census
/// existed. The demotion is what a run says about its own numbers (§13 test 5), so the difference
/// has to survive as far as the report.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum CensusAgreement {
    /// The run's census and the file's name the same terms with the same digests.
    TheSameCensus,
    /// **This run has no census**, so the fourth binding was not checked and nothing was demoted.
    /// The file's warrants stand — spec §2.1's grounds for not demoting on every read.
    NothingToCompareAgainst,
    /// They differ, at this term, in the census's own words — so every number in the file was
    /// demoted to no better than `supplied` (§2.1).
    FittedUnderAnother(String),
}

impl CensusAgreement {
    /// **Which term the two censuses first differ on**, and `None` where they do not differ or
    /// were never compared.
    ///
    /// **A convenience for a caller that only wants to print the term**, and deliberately not the
    /// representation: it collapses the two states this type exists to keep apart.
    #[must_use]
    pub fn term_they_differ_on(&self) -> Option<&str> {
        match self {
            Self::TheSameCensus | Self::NothingToCompareAgainst => None,
            Self::FittedUnderAnother(term) => Some(term),
        }
    }

    /// Whether the file's numbers were demoted because of this comparison.
    #[must_use]
    pub fn demoted_the_file(&self) -> bool {
        matches!(self, Self::FittedUnderAnother(_))
    }
}

/// One warrant, no better founded than *somebody handed this over*.
///
/// **Through [`Provenance`] rather than by matching on [`Warrant`]**, so the ladder that decides
/// which of two warrants is weaker stays in the one place that documents it.
fn no_better_than_supplied(warrant: Warrant) -> Warrant {
    Provenance::from(warrant)
        .weaker_of(Provenance::Supplied)
        .into()
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
    pub(super) fn the_run_the_fixture_was_fitted_by() -> (ReferenceDigest, ReadGroups) {
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

#[cfg(test)]
mod the_fourth_binding_demotes {
    //! **Step D3: a file fitted from another census of this cohort is used, and every number in
    //! it says so.**
    //!
    //! Spec §6's fourth binding is the one that does not refuse. The numbers are still numbers —
    //! a census is a store of *evidence*, and two censuses of one cohort differ in which loci
    //! were kept or at what depth, not in what a plant's genome is — so §2.1 keeps the file and
    //! demotes it wholesale.
    //!
    //! **Its failure is silent, which is why the test is two assertions and not one.** Warrants
    //! change what a run *reports* and never what it *computes* (spec §2), so a demotion that
    //! did not happen gives identical genotypes and a run that overstates every one of them, and
    //! a demotion that reached too far gives identical genotypes and a run that understates
    //! them. Only the pair — **the same answers, and no warrant left above `supplied`** — can
    //! tell those apart.

    use super::super::WarrantedValue;
    use super::super::tests::a_file_using_every_shape;
    use super::*;
    use crate::ng::types::ReadGroupId;

    /// The run the fixture was fitted by: its reference and lanes, from the sibling module that
    /// already builds them, plus the census the file itself names.
    fn this_run() -> (ReferenceDigest, ReadGroups, CensusIdentity) {
        let (reference, read_groups) = super::tests::the_run_the_fixture_was_fitted_by();
        (
            reference,
            read_groups,
            a_file_using_every_shape().fitted_from.census,
        )
    }

    /// The same census with one term's digest moved — a second census of the same cohort.
    fn a_census_of_the_same_cohort_recorded_otherwise(term: &str) -> CensusIdentity {
        let (.., mut census) = this_run();
        let moved = census
            .terms
            .iter_mut()
            .find(|carried| carried.term == term)
            .unwrap_or_else(|| panic!("the census names {term:?}"));
        moved.digest = "ff".repeat(16);
        census
    }

    fn read_for(file: &ParametersFile, census: &CensusIdentity) -> ParametersForThisRun {
        let (reference, read_groups, _) = this_run();
        file.to_run_parameters_for(&reference, &read_groups, Some(census))
            .unwrap_or_else(|error| panic!("{error}"))
    }

    /// **A file fitted under this run's own census keeps every warrant it was written with.**
    ///
    /// This is what spec §2.1 protects and the reason demotion is not unconditional: the same
    /// cohort called in psp mode from the fit in memory and in direct mode from the file that
    /// fit wrote must report the same warrants, or the two-mode oracle has to be told to ignore
    /// a difference that is real.
    #[test]
    fn a_file_fitted_under_this_runs_own_census_is_not_demoted() {
        let file = a_file_using_every_shape();
        let (_, _, census) = this_run();
        let read = read_for(&file, &census);

        assert_eq!(read.census, CensusAgreement::TheSameCensus);
        assert_eq!(
            read.from_file.parameters.calibration_by_read_group()[0].provenance,
            Provenance::FittedHere,
            "the file says this multiplier was fitted, and this run's census is the one it was \
             fitted under"
        );
        assert_eq!(
            read.from_file.inbreeding_by_sample[0].provenance,
            Provenance::FittedHere
        );
    }

    /// **A file fitted under another census of this cohort is used, and the run can say which
    /// term differed** — in the census's own words, which is why `CensusIdentity::of` mints them
    /// in `first_disagreement`'s vocabulary.
    #[test]
    fn a_file_fitted_under_another_census_is_used_and_names_the_term() {
        let file = a_file_using_every_shape();
        for term in [
            "the loci actually kept",
            "per-position depth cap",
            "selection seed",
        ] {
            let read = read_for(&file, &a_census_of_the_same_cohort_recorded_otherwise(term));
            assert_eq!(read.census.term_they_differ_on(), Some(term));
        }
    }

    /// **A run with no census keeps the file's warrants**, and spec §2.1 is what settles it
    /// rather than taste.
    ///
    /// Direct mode has no census (`run_streaming.md` §2) and is the mode this file format exists
    /// for. Demoting whenever there is nothing to compare against would be demoting on every read
    /// in that mode — which §2.1 considered and rejected, because the two-mode oracle requires
    /// the same cohort called from the file and from the fit in memory to report the same
    /// warrants for identical genotypes.
    ///
    /// **The three refusals still fire**, which is the half of §6 a missing census does not
    /// touch: they are about the reference, the samples and the read groups, none of which
    /// direct mode is missing.
    #[test]
    fn a_run_with_no_census_keeps_the_files_warrants_and_still_refuses_the_other_bindings() {
        let file = a_file_using_every_shape();
        let (reference, read_groups, census) = this_run();

        let read = file
            .to_run_parameters_for(&reference, &read_groups, None)
            .expect("a file this run's inputs match");
        // **Not `TheSameCensus`** — nothing was compared, and a report that said the censuses
        // matched would be claiming a check that never ran.
        assert_eq!(read.census, CensusAgreement::NothingToCompareAgainst);
        assert_eq!(
            read.from_file.inbreeding_by_sample[0].provenance,
            Provenance::FittedHere,
            "nothing was demoted"
        );

        // **And the same file under a census that disagrees *is* demoted**, so the `None` arm is
        // not merely a census that happens to match: the two arms give different answers on one
        // file.
        let demoted = file
            .to_run_parameters_for(
                &reference,
                &read_groups,
                Some(&a_census_of_the_same_cohort_recorded_otherwise(
                    "the loci actually kept",
                )),
            )
            .expect("a mismatched census demotes rather than refusing");
        assert_eq!(
            demoted.from_file.inbreeding_by_sample[0].provenance,
            Provenance::Supplied
        );
        assert!(demoted.census.demoted_the_file());

        // **The writer's spelling of the same fact gives the same answer.** A driver holding one
        // `CensusIdentity` would otherwise hand this `Some(of_a_run_with_no_census())` and get
        // every number demoted — the outcome the `None` arm exists to avoid.
        let as_the_writer_spells_it = file
            .to_run_parameters_for(
                &reference,
                &read_groups,
                Some(&CensusIdentity::of_a_run_with_no_census()),
            )
            .expect("a run with no census, said the other way");
        assert_eq!(
            as_the_writer_spells_it.census,
            CensusAgreement::NothingToCompareAgainst
        );
        assert_eq!(
            as_the_writer_spells_it.from_file.inbreeding_by_sample[0].provenance,
            Provenance::FittedHere
        );

        // The reference binding refuses whether or not there is a census to compare.
        let _ = census;
        assert!(
            file.to_run_parameters_for(&ReferenceDigest([0xab; 16]), &read_groups, None)
                .is_err(),
            "a run with no census still refuses a file fitted against another reference"
        );
    }

    /// **A term renamed while its digest stands still is a disagreement**, which comparing
    /// digests alone would miss.
    ///
    /// The twelve terms are named in the census's own words and compared in its own order
    /// (`CensusIdentity::of`), and a build that renames one while its value is unchanged has
    /// recorded its evidence under terms this build cannot line up against its own. Dropping the
    /// name half of the comparison left all 168 tests green.
    #[test]
    fn a_term_renamed_with_its_digest_unmoved_disagrees() {
        let mut file = a_file_using_every_shape();
        let (_, _, census) = this_run();
        let renamed = &mut file.fitted_from.census.terms[3];
        let digest = renamed.digest.clone();
        renamed.term = "what the catalog was built at".to_owned();

        assert_eq!(
            file.fitted_from.census.terms[3].digest, digest,
            "only the name moved"
        );
        assert_eq!(
            file.census_disagreement(&census).as_deref(),
            Some("what the catalog was built at"),
            "and the file's own name for it is what the run is told"
        );
    }

    /// **A census this build cannot even compare is a disagreement**, and it is the shape a file
    /// written by another build produces — one that recorded a thirteenth term, or twelve where
    /// this build knows thirteen.
    #[test]
    fn a_census_naming_a_different_number_of_terms_disagrees() {
        let mut file = a_file_using_every_shape();
        let (_, _, census) = this_run();

        let dropped = file.fitted_from.census.terms.pop().expect("twelve terms");
        assert_eq!(
            file.census_disagreement(&census).as_deref(),
            Some(dropped.term.as_str()),
            "the run's census names a term the file does not"
        );

        file.fitted_from.census.terms.push(dropped);
        file.fitted_from.census.terms.push(CensusTerm {
            term: "something a later build records".to_owned(),
            digest: "0e".repeat(16),
        });
        assert_eq!(
            file.census_disagreement(&census).as_deref(),
            Some("something a later build records")
        );
    }

    /// **⚑ The door applies the demotion**, which is the whole of what D3 composes and the one
    /// thing every other test here can be green without.
    ///
    /// The warrant assertions below and in `no_warrant_survives_the_demotion_stronger_than_
    /// supplied` are made against `demoted_to_no_better_than_supplied` *called directly*, and
    /// the one that goes through the door looks at the **agreeing** case. So a door that noticed
    /// the disagreement, reported it, and then projected the file undemoted left all 168 tests
    /// green — which is exactly this step's own silent failure: identical genotypes and a run
    /// that overstates every warrant it prints.
    #[test]
    fn the_door_demotes_and_not_only_the_method() {
        let read = read_for(
            &a_file_using_every_shape(),
            &a_census_of_the_same_cohort_recorded_otherwise("depth ladder edges"),
        );
        assert!(read.census.demoted_the_file());
        assert_eq!(
            read.from_file.parameters.calibration_by_read_group()[0].provenance,
            Provenance::Supplied,
            "the file says this multiplier was fitted here, and this run's evidence is not what \
             it was fitted from"
        );
        assert_eq!(
            read.from_file.inbreeding_by_sample[0].provenance,
            Provenance::Supplied
        );
        assert_eq!(
            read.from_file
                .parameters
                .ssr_slippage_fits()
                .stated_concentration_warrant(),
            Provenance::Supplied,
            "the bottom rung too, which is the number the base commit made carriable"
        );
        assert_eq!(
            read.from_file
                .parameters
                .repeat_tract_outlier_weight()
                .provenance(),
            Provenance::Defaulted,
            "and the one the demotion must not promote"
        );
    }

    /// **The demotion changes no number a locus is scored against.** Half of §13's fifth test,
    /// and the half that would pass on a demotion that did nothing — which is why the other half
    /// is above and below.
    #[test]
    fn the_demotion_changes_no_number_a_locus_is_scored_against() {
        let file = a_file_using_every_shape();
        let (_, _, census) = this_run();
        let kept = read_for(&file, &census).from_file.parameters;
        let demoted = read_for(
            &file,
            &a_census_of_the_same_cohort_recorded_otherwise("depth ladder edges"),
        )
        .from_file
        .parameters;

        assert_eq!(demoted.ploidy(), kept.ploidy());
        assert_eq!(demoted.read_group_count(), kept.read_group_count());
        assert_eq!(demoted.prior_seed(), kept.prior_seed());
        assert_eq!(
            demoted.contamination_by_read_group(),
            kept.contamination_by_read_group()
        );
        assert_eq!(
            demoted.inbreeding_coefficient_by_sample(),
            kept.inbreeding_coefficient_by_sample()
        );
        // **The batching is a number a locus is scored against**, not bookkeeping: it is the
        // population a contaminating read is drawn from.
        assert_eq!(demoted.sequencing_batches(), kept.sequencing_batches());
        assert_eq!(
            demoted.repeat_tract_outlier_weight().value(),
            kept.repeat_tract_outlier_weight().value()
        );
        for (was, is) in kept
            .calibration_by_read_group()
            .iter()
            .zip(demoted.calibration_by_read_group())
        {
            assert_eq!(is.scale, was.scale, "the multiplier a read is charged");
        }
        for ((key, was), (also, is)) in kept
            .ssr_substitution_rate()
            .zip(demoted.ssr_substitution_rate())
        {
            assert_eq!(key, also);
            assert_eq!(is.value, was.value, "the substitution rate inside a tract");
        }

        // **The slippage numbers and the length spectra are untouched by the demotion**, so they
        // compare whole rather than field by field: neither carries a `Provenance`, and what
        // they carry instead — which curve a number came off — is not a warrant and has no
        // *handed over* state to move to.
        for row in &file.repeat_tracts.slippage_by_stratum_and_group {
            for group in 0..file.fitted_from.read_groups.len() {
                let id = ReadGroupId(u32::try_from(group).expect("three lanes"));
                assert_eq!(
                    demoted
                        .ssr_slippage_fits()
                        .at(id, row.period, row.reference_repeats),
                    kept.ssr_slippage_fits()
                        .at(id, row.period, row.reference_repeats),
                );
            }
            let was = kept
                .ssr_slippage_fits()
                .length_spectrum_at(row.period, row.reference_repeats);
            let is = demoted
                .ssr_slippage_fits()
                .length_spectrum_at(row.period, row.reference_repeats);
            assert_eq!(is.rung(), was.rung());
            assert_eq!(is.concentration(), was.concentration());
            assert_eq!(is.fitted_weights(), was.fitted_weights());
        }
    }

    /// **⚑ No warrant is left above `supplied`, and none was promoted — which is *not* "every
    /// warrant `Supplied`", and cannot be.**
    ///
    /// Spec §13's fifth test says *same calls, every warrant `Supplied`*. The ladder says
    /// otherwise: `Provenance` ranks `Supplied` **above** `Defaulted`, so assigning `Supplied`
    /// would *promote* every defaulted number. The owner's ruling of 2026-08-30 is that the
    /// demotion is `weaker_of(file's warrant, Supplied)`, which is a no-op below `Supplied` and
    /// the right answer above it. So the property that holds is this one.
    #[test]
    fn no_warrant_survives_the_demotion_stronger_than_supplied() {
        let file = a_file_using_every_shape();
        let demoted = file.demoted_to_no_better_than_supplied();

        let before: Vec<(&str, Warrant)> = every_warrant_in(&file);
        let after: Vec<(&str, Warrant)> = every_warrant_in(&demoted);
        assert_eq!(before.len(), after.len());
        assert!(
            before
                .iter()
                .any(|(_, warrant)| *warrant == Warrant::FittedHere),
            "the fixture has something to demote"
        );
        assert!(
            before
                .iter()
                .any(|(_, warrant)| *warrant == Warrant::Defaulted),
            "and something that must not be promoted"
        );

        for ((what, was), (also, is)) in before.iter().zip(&after) {
            assert_eq!(what, also);
            let (was, is) = (Provenance::from(*was), Provenance::from(*is));
            assert_eq!(
                is,
                was.weaker_of(Provenance::Supplied),
                "{what} went from {was:?} to {is:?}"
            );
            assert!(
                is == Provenance::Supplied || is == Provenance::Defaulted,
                "{what} came out of the demotion as {is:?}"
            );
        }
    }

    /// **The one number where assigning `Supplied` rather than taking the weaker would be
    /// visible**: the repeat-tract outlier weight has no fitted state at all, and its two legal
    /// warrants are `supplied` and `defaulted` at the project's own 0.01. A demoted file that
    /// called it `supplied` would claim somebody chose that number when nobody did — and
    /// `validate` refuses `supplied`'s opposite, so the wrongness would be silent.
    #[test]
    fn the_outlier_weight_the_project_guessed_is_not_promoted_to_one_somebody_chose() {
        let file = a_file_using_every_shape();
        assert_eq!(
            file.stated_constants.repeat_tract_outlier_weight.warrant,
            Warrant::Defaulted
        );
        assert_eq!(
            file.demoted_to_no_better_than_supplied()
                .stated_constants
                .repeat_tract_outlier_weight
                .warrant,
            Warrant::Defaulted,
            "the demotion takes the weaker of the two, and defaulted is already weaker"
        );

        // And one a person did choose stays chosen.
        let mut chosen = file;
        chosen.stated_constants.repeat_tract_outlier_weight = WarrantedValue {
            value: 0.02,
            warrant: Warrant::Supplied,
            observations: None,
        };
        assert_eq!(
            chosen
                .demoted_to_no_better_than_supplied()
                .stated_constants
                .repeat_tract_outlier_weight
                .warrant,
            Warrant::Supplied
        );
    }

    /// **The junk decay is not promoted either** — the same property as the outlier weight's
    /// test above, for the second stated constant nothing fits.
    #[test]
    fn the_junk_decay_the_project_guessed_is_not_promoted_to_one_somebody_chose() {
        let file = a_file_using_every_shape();
        assert_eq!(
            file.stated_constants
                .repeat_tract_junk_decay_per_unit
                .warrant,
            Warrant::Defaulted
        );
        assert_eq!(
            file.demoted_to_no_better_than_supplied()
                .stated_constants
                .repeat_tract_junk_decay_per_unit
                .warrant,
            Warrant::Defaulted,
            "the demotion takes the weaker of the two, and defaulted is already weaker"
        );

        // And one a person did choose stays chosen.
        let mut chosen = file;
        chosen.stated_constants.repeat_tract_junk_decay_per_unit = WarrantedValue {
            value: 0.5,
            warrant: Warrant::Supplied,
            observations: None,
        };
        assert_eq!(
            chosen
                .demoted_to_no_better_than_supplied()
                .stated_constants
                .repeat_tract_junk_decay_per_unit
                .warrant,
            Warrant::Supplied
        );
    }

    /// **A demoted file is still a file this caller accepts**, which is not free: `validate`
    /// keys two of its rules on the bottom rung's warrant, and the demotion moves that warrant.
    ///
    /// Only one of the two rules is one the demotion can reach — `fitted_here` becomes
    /// `supplied`, and `defaulted` the demotion leaves alone, which is the whole of the
    /// `weaker_of` ruling. **The rule that had to change was the one this test is a guard on**:
    /// until the base commit, `validate` required `fitted_here` exactly where any stratum was
    /// fitted, and so refused every file the demotion produces.
    #[test]
    fn a_demoted_file_still_validates_and_still_projects() {
        let file = a_file_using_every_shape().demoted_to_no_better_than_supplied();
        file.validate()
            .expect("a demoted file still means something");
        file.to_run_parameters()
            .expect("and still projects to a run");
    }

    /// **The three refusals still refuse through this door**, and a file that means nothing is
    /// named by `validate` rather than blamed on the run.
    #[test]
    fn the_door_runs_validate_first_and_then_the_three_refusals() {
        let (reference, read_groups, census) = this_run();

        // A file at odds with itself: `validate`'s message, not a binding's.
        let mut at_odds = a_file_using_every_shape();
        at_odds.fitted_from.samples.reverse();
        let error = at_odds
            .to_run_parameters_for(&reference, &read_groups, Some(&census))
            .expect_err("a file whose sample list is not its own table's");
        assert!(
            matches!(error, ParametersFileError::Meaningless { .. }),
            "the file is what is wrong, and validate names it: {error}"
        );

        // And a good file against another run: a binding's message.
        let error = a_file_using_every_shape()
            .to_run_parameters_for(&ReferenceDigest([0xab; 16]), &read_groups, Some(&census))
            .expect_err("another reference");
        assert!(
            matches!(error, ParametersFileError::FittedFromOtherInputs { .. }),
            "{error}"
        );
    }

    /// Every warranted number in the file, named, in one order.
    fn every_warrant_in(file: &ParametersFile) -> Vec<(&str, Warrant)> {
        let mut warrants = vec![];
        for row in &file.base_quality_calibration.by_read_group {
            warrants.push(("a calibration", row.error_probability_multiplier.warrant));
        }
        for row in &file.inbreeding.by_sample {
            warrants.push((
                "an inbreeding coefficient",
                row.inbreeding_coefficient.warrant,
            ));
        }
        warrants.push((
            "the bottom rung",
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .warrant,
        ));
        for row in &file.repeat_tracts.substitution_rate_by_stratum {
            warrants.push(("a substitution rate", row.rate.warrant));
        }
        warrants.push((
            "the outlier weight",
            file.stated_constants.repeat_tract_outlier_weight.warrant,
        ));
        // **Five kinds, and the count is asserted so a sixth cannot be added upstream without
        // this walk being extended.** A warranted number the demotion forgets is exactly the
        // per-number exemption spec §2.1 says does not exist.
        assert_eq!(warrants.len(), 3 + 2 + 1 + 1 + 1);
        warrants
    }
}

#[cfg(test)]
mod what_the_run_counted_as_a_repeat {
    use super::super::tests::a_file_using_every_shape;
    use super::super::{ParametersFile, RepeatRouting};
    use crate::ng::region_typing::TypedRegionConfig;
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::types::Bp;

    /// The routing the fixture records — the catalog's own storage floors.
    fn the_fixtures_routing() -> StrRepeatCriteria {
        StrRepeatCriteria::default()
    }

    /// **Every axis of the criteria reaches the record.** A record that dropped one would let two
    /// runs that routed differently write identical files, and the difference the whole section
    /// exists to make visible would be invisible.
    #[test]
    fn the_record_carries_every_axis_the_catalog_was_asked_on() {
        let asked = StrRepeatCriteria::from(&TypedRegionConfig::default());
        let recorded = RepeatRouting::of(&asked);

        assert_eq!(
            recorded.min_copies,
            [8, 6, 6, 6, 5, 4],
            "ng's measured stutter onsets, one per period",
        );
        assert_eq!(recorded.min_period, asked.classification.periods.min());
        assert_eq!(recorded.max_period, asked.classification.periods.max());
        assert_eq!(recorded.max_str_len, asked.max_str_len_bp.get());
        assert_eq!(recorded.min_purity, asked.classification.min_purity);
        assert_eq!(recorded.min_flank_bp, asked.min_flank_bp.get());
        assert_eq!(recorded.min_score, asked.classification.min_score);
        assert_eq!(
            recorded.bundle_threshold,
            asked.classification.bundle_threshold
        );
    }

    /// **A run that routed as the file's run did finds no disagreement**, which is the control
    /// the eight cases below are read against.
    #[test]
    fn the_same_routing_is_no_disagreement() {
        let file = a_file_using_every_shape();
        assert_eq!(file.routing_disagreement(&the_fixtures_routing()), None);
    }

    /// **Each axis is named by the flag that moves it**, and each is tested on its own: a
    /// comparison that stopped at the first field, or that compared the record against itself,
    /// would pass with only one of these.
    ///
    /// The three axes with no flag are named in words rather than left unnamed, because a
    /// difference there is a catalog to rebuild and a person still has to be told which one.
    #[test]
    fn each_axis_that_differs_names_itself() {
        let file = a_file_using_every_shape();
        let with = |change: fn(&mut StrRepeatCriteria)| {
            let mut criteria = the_fixtures_routing();
            change(&mut criteria);
            file.routing_disagreement(&criteria)
        };

        assert_eq!(
            with(|c| c.classification.min_copies =
                crate::ng::region_typing::segment_criteria::MinCopies::default()),
            Some("--min-copies"),
        );
        assert_eq!(
            with(|c| c.classification.periods =
                crate::ng::tandem_repeat::PeriodRange::new(2, 6).expect("a range")),
            Some("--min-period"),
        );
        assert_eq!(
            with(|c| c.classification.periods =
                crate::ng::tandem_repeat::PeriodRange::new(1, 4).expect("a range")),
            Some("--max-period"),
        );
        assert_eq!(with(|c| c.max_str_len_bp = Bp(100)), Some("--max-str-len"));
        assert_eq!(
            with(|c| c.classification.min_purity = 0.9),
            Some("--min-purity"),
        );
        assert_eq!(with(|c| c.min_flank_bp = Bp(30)), Some("the flank floor"));
        assert_eq!(
            with(|c| c.classification.min_score = 30),
            Some("the scanner score floor"),
        );
        assert_eq!(
            with(|c| c.classification.bundle_threshold = 20),
            Some("the bundling distance"),
        );
    }

    /// **A file that does not say what it routed with disagrees with nothing** — spec §5's rule
    /// that absence is not a value, and here the difference between *this file made no claim* and
    /// *this file claims the defaults*. A file written by a build older than the section, or by
    /// hand, is the case.
    #[test]
    fn a_file_with_no_routing_record_makes_no_claim_to_disagree_with() {
        let mut file = a_file_using_every_shape();
        file.repeat_routing = None;
        assert_eq!(
            file.routing_disagreement(&StrRepeatCriteria::from(&TypedRegionConfig::default())),
            None,
            "a run must not be told a file disagrees when the file said nothing",
        );
    }

    /// **A routing difference does not touch a single warrant**, which is what separates it from
    /// the census mismatch beside it: those numbers were fitted elsewhere, these were not, and
    /// only the ground they are applied to has moved.
    #[test]
    fn a_routing_difference_leaves_every_number_as_warranted_as_it_was() {
        let file = a_file_using_every_shape();
        let mut routed_otherwise = the_fixtures_routing();
        routed_otherwise.max_str_len_bp = Bp(100);
        assert!(file.routing_disagreement(&routed_otherwise).is_some());

        // The comparison is a question, not a step: nothing about the file moves when it is
        // asked, so a caller that reports the difference and calls on is calling on the same
        // numbers it would have used had it never asked.
        let untouched: &ParametersFile = &file;
        assert_eq!(untouched, &a_file_using_every_shape());
    }
}
