//! **What a psp's census is, judged before a fit trusts it** — the verdict and nothing that
//! reaches it.
//!
//! A psp carries its sample's census as the file's closing payload
//! (`doc/devel/ng/spec/psp_census_pair.md` §3). A census is a cache, so a psp can carry one that
//! this build cannot read, or one written against a different set of loci from the one the run in
//! hand rebuilds. `estimate-parameters` will refuse such a cohort and name every sample in it
//! (spec §4, §4.1, plan step C3); `regenerate-census` will rebuild exactly those (spec §8, plan
//! step D1). Both ask the same question, and [`CensusVerdict`] is the answer. Neither command
//! reads it yet.

use crate::ng::parameter_estimation::joint::census_file::VERSION;

/// **What a run should do with the census in one psp**, in the words spec §4.2 uses.
///
/// **A verdict and not an action**: the same answers serve a command that refuses
/// ([`estimate-parameters`](crate::pop_var_caller_exp::estimate_parameters)) and one that
/// rebuilds (`regenerate-census`, plan step D1), and which of those happens is the command's
/// business.
///
/// **The four causes cost three different reads to reach.** [`NoCensus`](Self::NoCensus) is in
/// the psp's footer, which was read when the file was opened. [`AnotherFormat`](Self::AnotherFormat)
/// and [`NotACensus`](Self::NotACensus) are one short read at the trailer's front.
/// [`AnotherSelection`](Self::AnotherSelection) needs the run's reference read and its selection
/// rebuilt — once a run, and then one digest compared a sample. A command that has read only the
/// footer and the trailer's front reports what those two can see and stops there.
///
/// **An older census is named as an older census, never as damage.** A version word this build
/// does not know reaches the user as *malformed* today
/// ([`decode_census`](crate::ng::parameter_estimation::joint::census_file::decode_census)), which
/// sends someone looking for a corrupted file when what happened is that this build changed what
/// a census holds (spec §4.2). Plan step B2 reads the version before that check, so it reaches
/// them as [`AnotherFormat`](Self::AnotherFormat) instead.
///
/// **What is not a verdict: a census that decodes badly past its version word.** A trailer whose
/// magic and version are this build's and whose sections are truncated is refused by the reader
/// that meets it, as a [`CensusError`](crate::ng::parameter_estimation::joint::census::CensusError),
/// and no cheap read can tell it from a whole one. Regenerating fixes it, but nothing here says
/// so, because nothing here has looked.
///
/// **Grouping a cohort's verdicts by cause** — spec §4.1's report — groups on the whole value, so
/// two samples carrying censuses of two different old versions are two causes and not one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CensusVerdict {
    /// Nothing the judgement that produced this verdict found wrong.
    ///
    /// **How much that covers is the producer's to say and not this variant's.** The three causes
    /// of spec §4.2 cost three different reads, so a judgement made from the psp's head has
    /// looked at two of them and cannot have looked at the third. `Fresh` from such a judgement
    /// means *nothing the head can see is wrong*; spec §8's rule for skipping a psp in
    /// `regenerate-census` is all three causes, so a skip needs the selection checked too.
    Fresh,

    /// The psp carries no census at all: its trailer is empty.
    ///
    /// A psp written before the census moved into the trailer, or one whose trailer an
    /// [`append`](crate::ng::psp::PspWriter::append) discarded (spec §3.4).
    NoCensus,

    /// The version word at the front of the psp's census is not the one this build reads.
    ///
    /// **Ordinarily that means another version of this program wrote it**, which is what the
    /// message says; a version word corrupted behind an intact magic looks the same from here.
    /// Either way the census cannot be decoded and regenerating it is the fix.
    AnotherFormat {
        /// The version word at the front of this psp's census.
        ///
        /// **The version this build reads is not stored beside it**: it is [`VERSION`], a
        /// constant, so no verdict can carry a pair of numbers that disagrees with the build
        /// printing it.
        version_in_the_psp: u16,
    },

    /// The psp's trailer holds bytes, and they are not a census: they do not begin with a
    /// census's magic, or they end before the version word.
    ///
    /// **A trailer of zero bytes is [`NoCensus`](Self::NoCensus) and not this**, so that the two
    /// causes partition the psps rather than overlap on the empty one.
    ///
    /// **Not one of spec §4.2's three causes, and it is here because the read that reaches them
    /// can meet this instead.** §4.2 asks for the version word to be read before
    /// [`decode_census`](crate::ng::parameter_estimation::joint::census_file::decode_census) can
    /// call an old census malformed; reading that word means first knowing the bytes are a
    /// census, and bytes that are not are damage rather than an old format. Regenerating the
    /// census rewrites exactly those bytes and touches nothing else
    /// ([`replace_trailer`](crate::ng::psp::replace_trailer)), so it is still the repair to try.
    NotACensus,

    /// The psp carries a census this build reads, written against a different set of loci from
    /// the one this run rebuilds (spec §4.2, third row).
    ///
    /// A build whose selection constants changed, or a run under a different reference or
    /// catalog. It is the one cause that cannot be reached without building the selection.
    AnotherSelection,
}

impl std::fmt::Display for CensusVerdict {
    /// **Written to be read after the sample's name and its psp's path** — spec §4.3's report is
    /// one line a sample, and this is that line's last column.
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fresh => write!(out, "carries a census this build reads"),
            Self::NoCensus => write!(out, "carries no census"),
            Self::AnotherFormat { version_in_the_psp } => write!(
                out,
                "carries a census built by {} version of this program: it is version \
                 {version_in_the_psp} and this build reads version {VERSION}",
                match *version_in_the_psp < VERSION {
                    true => "an older",
                    false => "a newer",
                }
            ),
            Self::NotACensus => write!(out, "carries a trailer that is not a census"),
            Self::AnotherSelection => {
                write!(out, "carries a census built under another selection")
            }
        }
    }
}

impl CensusVerdict {
    /// The verdict for a census whose version word is `version_in_the_psp`, or `None` when that
    /// is the version this build reads.
    ///
    /// **`None` and not [`Fresh`](Self::Fresh), because a version word is one of spec §4.2's
    /// three causes and freshness is all three.** A constructor that answered `Fresh` here would
    /// hand a caller the whole judgement in exchange for one comparison — and `regenerate-census`
    /// skips a psp on that word (spec §8).
    #[must_use]
    pub fn of_a_version_word(version_in_the_psp: u16) -> Option<Self> {
        match version_in_the_psp == VERSION {
            true => None,
            false => Some(Self::AnotherFormat { version_in_the_psp }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The cause a user reads is the whole of what a refusal gives them to act on**, so each
    /// one is pinned rather than left to whatever a `Debug` would print.
    #[test]
    fn each_verdict_says_what_the_psp_carries() {
        assert_eq!(
            CensusVerdict::Fresh.to_string(),
            "carries a census this build reads"
        );
        assert_eq!(CensusVerdict::NoCensus.to_string(), "carries no census");
        assert_eq!(
            CensusVerdict::NotACensus.to_string(),
            "carries a trailer that is not a census"
        );
        assert_eq!(
            CensusVerdict::AnotherSelection.to_string(),
            "carries a census built under another selection"
        );
    }

    /// **A version this build does not read is named twice over** — what the file holds and what
    /// this binary reads — because *the wrong version* alone does not say which of the two the
    /// user should change.
    ///
    /// **The expected string names the found version as a literal and this build's from the
    /// constant**, so that a message printing one of them in both places fails whichever of the
    /// two it kept.
    #[test]
    fn a_census_of_another_format_names_both_versions() {
        assert_eq!(
            CensusVerdict::AnotherFormat {
                version_in_the_psp: 1
            }
            .to_string(),
            format!(
                "carries a census built by an older version of this program: it is version 1 \
                 and this build reads version {VERSION}"
            )
        );
    }

    /// **A census from a build newer than this one is not an older census**, and saying so would
    /// send the reader to rebuild psps that a newer binary already reads.
    #[test]
    fn a_census_from_a_newer_build_is_named_as_newer() {
        assert_eq!(
            CensusVerdict::AnotherFormat {
                version_in_the_psp: VERSION + 1
            }
            .to_string(),
            format!(
                "carries a census built by a newer version of this program: it is version {} \
                 and this build reads version {VERSION}",
                VERSION + 1
            )
        );
    }

    /// **The version word alone answers only one of spec §4.2's three questions**, so the version
    /// this build reads is `None` — nothing to report — rather than a verdict of fresh.
    #[test]
    fn this_builds_own_version_is_nothing_to_report() {
        assert_eq!(CensusVerdict::of_a_version_word(VERSION), None);
        assert_eq!(
            CensusVerdict::of_a_version_word(VERSION - 1),
            Some(CensusVerdict::AnotherFormat {
                version_in_the_psp: VERSION - 1,
            })
        );
        assert_eq!(
            CensusVerdict::of_a_version_word(VERSION + 1),
            Some(CensusVerdict::AnotherFormat {
                version_in_the_psp: VERSION + 1,
            })
        );
    }
}
