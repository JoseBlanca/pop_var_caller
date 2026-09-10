//! **What a psp's census is, judged before a fit trusts it** — the verdict, the two cheap reads
//! that reach two of its causes, and a cohort judged whole.
//!
//! A psp carries its sample's census as the file's closing payload
//! (`doc/devel/ng/spec/psp_census_pair.md` §3). A census is a cache, so a psp can carry one that
//! this build cannot read, or one written against a different set of loci from the one the run in
//! hand rebuilds. `estimate-parameters` will refuse such a cohort and name every sample in it
//! (spec §4, §4.1, plan step C3); `regenerate-census` will rebuild exactly those (spec §8, plan
//! step D1). Both ask the same question, and [`CensusVerdict`] is the answer. Neither command
//! reads it yet.
//!
//! **A cohort is judged in one pass and every sample in it is named**
//! ([`what_the_heads_say_about_every_census_in_a_cohort`]), because the wait a refusal saves is
//! per sample: regenerating one census is a quarter of an hour (spec §4), so a run that reported
//! the first stale psp and stopped would cost that wait once a stale sample, one after another,
//! to learn a job that fits in one message.

use std::path::{Path, PathBuf};

use crate::ng::parameter_estimation::joint::census_file::{
    BYTES_THAT_NAME_THE_VERSION, VERSION, version_word_of,
};
use crate::ng::psp::{PspReadError, PspReader};

use super::psp_caller::OpenPspCohort;

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

/// **What the footer and the front of the trailer say about a psp's census** — two of spec
/// §4.2's three causes, for one short read or none.
///
/// The footer says how long the trailer is, and it was read when the file was opened: a length of
/// zero is [`NoCensus`](CensusVerdict::NoCensus), and no byte of trailer is read. Otherwise the
/// first
/// [`BYTES_THAT_NAME_THE_VERSION`] bytes of the trailer are read — one seek and one read — and
/// they either fail to be a census at all ([`NotACensus`](CensusVerdict::NotACensus)), name a
/// version this build does not have ([`AnotherFormat`](CensusVerdict::AnotherFormat)), or name
/// this build's own.
///
/// **[`Fresh`](CensusVerdict::Fresh) here means only that**: nothing these two reads can see is
/// wrong. The third cause — a census built under another selection — needs the run's reference
/// read and its selection rebuilt, and this function has done neither. A command that may only
/// skip a psp once all three causes are ruled out (`regenerate-census`, spec §8) checks the
/// selection as well; a command that refuses before the reference is read
/// (`estimate-parameters`, spec §4.2) is what this alone serves.
///
/// **The whole trailer is never read**, and
/// [`trailer_bytes_read`](crate::ng::psp::trailer_bytes_read) is what says so: a census in a
/// whole-genome psp is tens of megabytes (spec §3.1) and a cohort is thousands of psps, which is
/// the memory the census's own reader is lazy to avoid (spec §5).
///
/// # Errors
///
/// [`PspReadError`] when the seek or the read fails — a file truncated or replaced after it was
/// opened, or an I/O fault. **A footer pointing past the file's own end is not one of these**:
/// `PspReader::open` proves the trailer ends exactly where the footer begins, so such a psp never
/// reaches here. **A psp that cannot be read is not a stale psp**: regenerating its census would
/// not fix it, so it is not a verdict.
pub fn what_the_footer_and_the_trailers_head_say_about_a_census(
    psp: &mut PspReader,
) -> Result<CensusVerdict, PspReadError> {
    if psp.footer().trailer_bytes == 0 {
        return Ok(CensusVerdict::NoCensus);
    }
    let head = psp.trailer_head(BYTES_THAT_NAME_THE_VERSION)?;
    let Some(version) = version_word_of(&head) else {
        return Ok(CensusVerdict::NotACensus);
    };
    Ok(CensusVerdict::of_a_version_word(version).unwrap_or(CensusVerdict::Fresh))
}

/// **What one sample's psp says about its census, and which sample and which file said it** —
/// one row of the report `estimate-parameters` refuses a cohort with (spec §4.3).
///
/// **The individual and the file are both here because neither names the other.** A psp's header
/// carries the sample it holds and not where it was stored, and a file's name is whatever the walk
/// was told to write — the two agree only by convention. What the user acts on is the path, which
/// is what they copy into `regenerate-census`; what tells them which of their accessions it is, is
/// the name in the header.
#[derive(Debug)]
pub struct JudgedPsp {
    /// The individual whose psp this is, out of the file's own header.
    pub sample: String,
    /// The file, as the run was given it.
    pub psp: PathBuf,
    /// What its census turned out to be — or the failure that means it has no verdict.
    ///
    /// **A psp that cannot be *read* is not a psp with a stale census**, and the two must not
    /// arrive as one: regenerating the census of a truncated or vanished file would not fix it,
    /// and the report that tells a user which samples to regenerate would be sending them at a
    /// file whose fault is elsewhere (spec §4.2's causes are all about what the trailer holds).
    pub verdict: Result<CensusVerdict, PspReadError>,
}

/// **Every psp of an opened cohort judged, one verdict a sample, in the run's sample order** —
/// spec §4.1.
///
/// **Nothing here stops early, and that is the whole of it.** Before the census moved into the
/// psp, the cohort opener returned at the first census it could not check, so a run over a cohort
/// with three stale ones reported one, the user regenerated it, and the next run reported the
/// second. Regenerating a census is a quarter of an hour a sample (spec §4), so what that costs
/// is three rounds of that wait to learn a job that could have been named in one. **A file that
/// will not read does not stop it either**: its row carries the read failure and the psps after
/// it are still judged.
///
/// **Two of spec §4.2's three causes**, because that is what
/// [`what_the_footer_and_the_trailers_head_say_about_a_census`] reaches: a verdict of
/// [`Fresh`](CensusVerdict::Fresh) here means nothing the footer and the trailer's front can see
/// is wrong, and a census built under another selection looks exactly like a fresh one until the
/// run's reference is read and its selection rebuilt. That is what lets this be taken before the
/// reference is opened, which is what makes a missing census an immediate refusal rather than one
/// that arrives after minutes of work (spec §4.2).
///
/// **The cost is at most one seek and ten bytes a sample, and nothing at all for a psp with no
/// census**, and no census is decoded — see
/// [`trailer_bytes_read`](crate::ng::psp::trailer_bytes_read). A thousand-sample cohort's
/// censuses are tens of gigabytes (spec §5) and this reads ten kilobytes of them.
#[must_use]
pub fn what_the_heads_say_about_every_census_in_a_cohort(
    cohort: &mut OpenPspCohort,
) -> Vec<JudgedPsp> {
    cohort
        .each_psp_with_its_path()
        .map(|(path, psp)| what_one_psps_head_says_about_its_census(path, psp))
        .collect()
}

/// One psp judged, named by the individual in its header and by `path`.
fn what_one_psps_head_says_about_its_census(path: &Path, psp: &mut PspReader) -> JudgedPsp {
    // **Cloned before the judgement, not after**: reading the trailer's front takes the reader
    // mutably, so the name has to be out of the header first.
    let sample = psp.header().sample.clone();
    JudgedPsp {
        sample,
        psp: path.to_path_buf(),
        verdict: what_the_footer_and_the_trailers_head_say_about_a_census(psp),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use tempfile::TempDir;

    use crate::ng::parameter_estimation::joint::census_file::tests_support::a_census_this_build_wrote;
    use crate::ng::psp::PspWriter;
    use crate::ng::psp::writer::tests_support::a_header;

    /// A finished psp holding no records, sealed with `trailer`.
    ///
    /// **Sealed through the writer**, so the footer this judgement reads is the one a walk would
    /// have written rather than a hand-built one.
    ///
    /// `sample` is the individual the header names, and the file is written to
    /// `<sample>.psp` — the ordinary case, and what a test judging one psp wants.
    fn a_psp_sealed_with(dir: &TempDir, sample: &str, trailer: &[u8]) -> PathBuf {
        a_psp_of_sample_written_to(dir, sample, sample, trailer)
    }

    /// The same, with the individual the header names and the file's stem given separately.
    ///
    /// **The two are separated because nothing outside a convention ties them.** A walk writes
    /// whatever file it was told to; a cohort's psps can as easily be `psp-1.psp` as
    /// `SRR7279481.psp`, and a report that read the individual off the path would be right on
    /// every fixture that names them alike and wrong on real files.
    fn a_psp_of_sample_written_to(
        dir: &TempDir,
        sample: &str,
        stem: &str,
        trailer: &[u8],
    ) -> PathBuf {
        let mut header = a_header(1_000);
        header.sample = sample.to_string();
        // **Each sample's lanes are named after it**, for the cohort fixture's sake: an
        // `@RG ID` is unique across a whole cohort, and the shared header fixture names every
        // sample's lanes alike, so a cohort built from it is refused before a census is judged.
        for lane in &mut header.read_groups {
            lane.id = format!("{sample}-{}", lane.id);
        }
        let path = dir.path().join(format!("{stem}.psp"));
        let writer = PspWriter::create(&path, header).expect("the header writes");
        let _ = writer.finish(trailer).expect("the file seals");
        path
    }

    fn the_verdict_on(path: &Path) -> CensusVerdict {
        let mut psp = PspReader::open(path).expect("a finished psp opens");
        what_the_footer_and_the_trailers_head_say_about_a_census(&mut psp).expect("its head reads")
    }

    /// **A psp carrying the census this build writes passes the head's two checks.**
    ///
    /// The fixture's census is the one `census_file`'s own tests round-trip, written by
    /// `write_census` — so this also says the judgement and the writer agree about where a
    /// census's version word is.
    #[test]
    fn a_psp_carrying_this_builds_census_passes_what_the_head_can_check() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        assert!(
            census.len() > BYTES_THAT_NAME_THE_VERSION,
            "a fixture whose whole census is the head would not show the head is enough",
        );
        let psp = a_psp_sealed_with(&dir, "alpha", &census);

        assert_eq!(the_verdict_on(&psp), CensusVerdict::Fresh);
    }

    /// **An empty trailer is answered out of the footer, which opening already read.**
    ///
    /// Every psp this tree wrote before the census moved into the trailer is this shape, and so
    /// is one whose trailer an append discarded (spec §3.4).
    #[test]
    fn a_psp_whose_trailer_is_empty_carries_no_census() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let psp = a_psp_sealed_with(&dir, "alpha", b"");

        assert_eq!(the_verdict_on(&psp), CensusVerdict::NoCensus);
    }

    /// **A census of another version is named as one, not as damage** (spec §4.2), and the
    /// version it names is the one in the file.
    ///
    /// The fixture is a real census of this build with its version word overwritten, so
    /// everything but that word is what a census actually looks like — a hand-built ten bytes
    /// would pass whatever the judgement did with the rest.
    ///
    /// **Both directions**, because a build older than the psp is as ordinary as a psp older
    /// than the build: one person upgrades, another has not yet.
    #[test]
    fn a_census_written_by_another_build_names_its_version() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        for version in [VERSION - 1, VERSION + 1] {
            let census = a_census_of_version(version);
            let psp = a_psp_sealed_with(&dir, &format!("of-version-{version}"), &census);

            assert_eq!(
                the_verdict_on(&psp),
                CensusVerdict::AnotherFormat {
                    version_in_the_psp: version,
                },
            );
        }
    }

    /// **Judging a psp reads the front of its census and not its census**, which is the whole
    /// reason this is two short reads rather than `trailer()`.
    ///
    /// **Only the byte count tells the two apart** — every verdict above is the same whichever
    /// of them ran — and at a thousand samples the difference is ten kilobytes against tens of
    /// gigabytes (spec §5).
    #[test]
    fn judging_a_psp_reads_the_front_of_its_census_and_not_its_census() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let path = a_psp_sealed_with(&dir, "alpha", &census);
        let mut psp = PspReader::open(&path).expect("a finished psp opens");

        crate::ng::psp::reset_trailer_bytes_read();
        let verdict = what_the_footer_and_the_trailers_head_say_about_a_census(&mut psp)
            .expect("its head reads");

        assert_eq!(verdict, CensusVerdict::Fresh);
        assert_eq!(
            crate::ng::psp::trailer_bytes_read(),
            BYTES_THAT_NAME_THE_VERSION as u64,
            "the census in this psp is {} bytes",
            census.len(),
        );
    }

    /// **A psp with no census is judged without taking a byte of trailer**: the footer that says
    /// the trailer is empty was read when the file was opened.
    #[test]
    fn judging_a_psp_that_has_no_census_takes_no_trailer_bytes() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = a_psp_sealed_with(&dir, "alpha", b"");
        let mut psp = PspReader::open(&path).expect("a finished psp opens");

        crate::ng::psp::reset_trailer_bytes_read();
        let verdict = what_the_footer_and_the_trailers_head_say_about_a_census(&mut psp)
            .expect("its footer was read at open");

        assert_eq!(verdict, CensusVerdict::NoCensus);
        assert_eq!(crate::ng::psp::trailer_bytes_read(), 0);
    }

    /// **A trailer holding something that is not a census is damage, and says so.**
    ///
    /// The payload is the one psp's own tests seal a file with — a trailer that was legal before
    /// the census moved into it, and spec `psp_file_format.md` §3.4 keeps the payload the
    /// writer's business, so a psp meeting this reader with other bytes in it is possible.
    #[test]
    fn a_trailer_that_is_not_a_census_is_not_reported_as_an_old_one() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let psp = a_psp_sealed_with(&dir, "alpha", b"a per-sample summary");

        assert_eq!(the_verdict_on(&psp), CensusVerdict::NotACensus);
    }

    /// **A trailer that stops inside the version word is not a census either**, where a reader
    /// that took whatever bytes were there would read one byte of a census and one of the
    /// footer as a version.
    #[test]
    fn a_trailer_ending_inside_the_version_word_is_not_a_census() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let psp = a_psp_sealed_with(&dir, "alpha", &census[..BYTES_THAT_NAME_THE_VERSION - 1]);

        assert_eq!(the_verdict_on(&psp), CensusVerdict::NotACensus);
    }

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

    /// A cohort of psps, each sealed with the trailer named beside its sample, opened together
    /// in the order given.
    ///
    /// The headers are otherwise the one fixture, so the cohort agrees about its ground, its
    /// catalog and its repeat-tract criteria and reaches the judgement rather than a refusal.
    ///
    /// **No file here is named after the sample in it** — psp *i* of the list is written to
    /// `psp-<i>.psp`. A cohort whose stems were its sample names could not tell a report that
    /// reads the individual out of the header from one that reads it off the path, and those are
    /// two different answers on any cohort a walk named by run accession or by lane.
    fn a_cohort_sealed_with(dir: &TempDir, psps: &[(&str, &[u8])]) -> OpenPspCohort {
        let paths: Vec<PathBuf> = psps
            .iter()
            .enumerate()
            .map(|(index, (sample, trailer))| {
                a_psp_of_sample_written_to(dir, sample, &format!("psp-{index}"), trailer)
            })
            .collect();
        OpenPspCohort::open(&paths).expect("the fixtures agree about everything but their censuses")
    }

    /// The file psp `index` of a cohort fixture was written to.
    fn the_cohorts_file(dir: &TempDir, index: usize) -> PathBuf {
        dir.path().join(format!("psp-{index}.psp"))
    }

    /// A census of this build with its version word overwritten — a psp another build wrote.
    fn a_census_of_version(version: u16) -> Vec<u8> {
        let mut census = a_census_this_build_wrote();
        let word_at = BYTES_THAT_NAME_THE_VERSION - size_of::<u16>();
        census[word_at..BYTES_THAT_NAME_THE_VERSION].copy_from_slice(&version.to_le_bytes());
        census
    }

    /// The individual, the file and the verdict of each row, for a cohort every psp of which
    /// could be read.
    fn each_sample_file_and_verdict(rows: &[JudgedPsp]) -> Vec<(&str, &Path, CensusVerdict)> {
        rows.iter()
            .map(|row| {
                let verdict = *row
                    .verdict
                    .as_ref()
                    .unwrap_or_else(|error| panic!("{} reads: {error}", row.sample));
                (row.sample.as_str(), row.psp.as_path(), verdict)
            })
            .collect()
    }

    /// **Three stale psps in a cohort of five are all named, and in the order the run was given
    /// them** — spec §4.1, and the whole reason this is a pass over the cohort rather than the
    /// search the census cohort opener does today.
    ///
    /// **The first psp is fresh and the last is not**, so a judgement that stopped at the first
    /// stale one would come back with two rows where five are asserted, and one that judged only
    /// the first psp would come back with one.
    ///
    /// **The three are stale for three different causes**, because grouping the report by cause
    /// (spec §4.3) is worth nothing if the pass reports whatever the first stale psp was.
    ///
    /// **The samples are not in alphabetical order**, so a pass that sorted them — or that
    /// walked a map keyed by name — gives a different sequence from the one asserted.
    ///
    /// **Each row's individual and file are asserted together**, and no file here is named after
    /// the sample in it, so a pass that read the individual off the path fails on every row.
    #[test]
    fn every_stale_psp_of_a_cohort_is_named_and_not_only_the_first() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let older = a_census_of_version(VERSION - 1);
        let mut cohort = a_cohort_sealed_with(
            &dir,
            &[
                ("delta", &census),
                ("alpha", b""),
                ("echo", &census),
                ("bravo", &older),
                ("charlie", b"a per-sample summary"),
            ],
        );

        let rows = what_the_heads_say_about_every_census_in_a_cohort(&mut cohort);

        assert_eq!(
            each_sample_file_and_verdict(&rows),
            vec![
                (
                    "delta",
                    the_cohorts_file(&dir, 0).as_path(),
                    CensusVerdict::Fresh
                ),
                (
                    "alpha",
                    the_cohorts_file(&dir, 1).as_path(),
                    CensusVerdict::NoCensus
                ),
                (
                    "echo",
                    the_cohorts_file(&dir, 2).as_path(),
                    CensusVerdict::Fresh
                ),
                (
                    "bravo",
                    the_cohorts_file(&dir, 3).as_path(),
                    CensusVerdict::AnotherFormat {
                        version_in_the_psp: VERSION - 1,
                    }
                ),
                (
                    "charlie",
                    the_cohorts_file(&dir, 4).as_path(),
                    CensusVerdict::NotACensus
                ),
            ],
        );
    }

    /// **A cohort of one stale psp is one row** — not none, and not a refusal.
    ///
    /// One sample is the low end of the range this caller is built for (`CLAUDE.md`), and a pass
    /// written as *report the samples after the first* or *report nothing when there is nothing
    /// to compare* is a pass that says a single-sample run has no work to do when it has all of
    /// it.
    #[test]
    fn a_cohort_of_one_stale_psp_is_one_row() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let mut cohort = a_cohort_sealed_with(&dir, &[("alpha", b"")]);

        let rows = what_the_heads_say_about_every_census_in_a_cohort(&mut cohort);

        assert_eq!(
            each_sample_file_and_verdict(&rows),
            vec![(
                "alpha",
                the_cohorts_file(&dir, 0).as_path(),
                CensusVerdict::NoCensus
            )],
        );
    }

    /// **Each row names the file the run was given**, which is what a user copies into
    /// `regenerate-census` (spec §4.3) — and two psps of one cohort are two files, so a row that
    /// carried the cohort's first path, or none, would still name the right samples.
    #[test]
    fn each_row_names_the_psp_the_run_was_given() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let mut cohort = a_cohort_sealed_with(&dir, &[("alpha", b""), ("bravo", &census)]);

        let rows = what_the_heads_say_about_every_census_in_a_cohort(&mut cohort);

        assert_eq!(
            rows.iter().map(|row| row.psp.clone()).collect::<Vec<_>>(),
            vec![the_cohorts_file(&dir, 0), the_cohorts_file(&dir, 1)],
        );
    }

    /// **A psp that cannot be read is not a stale psp, and it does not stop the rest being
    /// judged.**
    ///
    /// Regenerating the census of a file that will not read would not fix it, so its row carries
    /// the read failure rather than a verdict — and the psps after it in the cohort are judged
    /// all the same, which is what a pass written with `?` would not do: it would come back with
    /// one error where five rows are asserted, reproducing for a truncated file exactly the
    /// early return spec §4.1 exists to remove.
    ///
    /// **The file is emptied after the cohort is open**, which is when a psp becomes unreadable
    /// in a run: opening read its header and footer, and the trailer those point at is gone by
    /// the time the judgement seeks to it.
    #[test]
    fn a_psp_that_will_not_read_has_no_verdict_and_does_not_stop_the_cohort() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let newer = a_census_of_version(VERSION + 1);
        let mut cohort = a_cohort_sealed_with(
            &dir,
            &[
                ("alpha", &census),
                ("bravo", &census),
                ("charlie", b""),
                ("delta", &census),
                ("echo", &newer),
            ],
        );
        std::fs::File::create(the_cohorts_file(&dir, 1)).expect("the file empties");

        let rows = what_the_heads_say_about_every_census_in_a_cohort(&mut cohort);

        let names: Vec<&str> = rows.iter().map(|row| row.sample.as_str()).collect();
        assert_eq!(names, vec!["alpha", "bravo", "charlie", "delta", "echo"]);
        assert!(
            rows[1].verdict.is_err(),
            "an emptied psp has nothing at the trailer offset its footer named, and the \
             judgement came back as {:?}",
            rows[1].verdict,
        );
        assert_eq!(
            *rows[2].verdict.as_ref().expect("charlie reads"),
            CensusVerdict::NoCensus,
            "the psps after the one that would not read are judged all the same",
        );
        assert_eq!(
            *rows[4].verdict.as_ref().expect("echo reads"),
            CensusVerdict::AnotherFormat {
                version_in_the_psp: VERSION + 1,
            },
        );
    }

    /// **A cohort is judged for one seek and ten bytes a sample, and no census is decoded** —
    /// the property that makes a refusal cheap enough to take before the reference is read
    /// (spec §4.2), and the only thing that tells this pass from one that took each trailer
    /// whole.
    ///
    /// Of the four samples here, three carry a census and one carries none, so the whole pass
    /// costs 30 bytes against those three censuses' own size, which the message names; at a
    /// thousand whole-genome samples it is ten kilobytes against tens of gigabytes (spec §5).
    ///
    /// **The psp with no census is the first, and the two assertions catch different faults.**
    /// The byte count catches a pass that reads whole trailers, and also one that never reaches
    /// the last psp: measured, a pass ending one psp early reads 20 bytes where 30 are asserted,
    /// which it would not if the census-less psp were the one it skipped. The row count is what
    /// catches a pass that judges every psp and then drops one — measured at 3 rows against 4,
    /// with the byte count still 30.
    #[test]
    fn judging_a_cohort_reads_the_front_of_each_census_and_no_census() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let census = a_census_this_build_wrote();
        let mut cohort = a_cohort_sealed_with(
            &dir,
            &[
                ("alpha", b""),
                ("bravo", &census),
                ("charlie", &census),
                ("delta", &census),
            ],
        );

        crate::ng::psp::reset_trailer_bytes_read();
        let rows = what_the_heads_say_about_every_census_in_a_cohort(&mut cohort);

        assert_eq!(
            crate::ng::psp::trailer_bytes_read(),
            3 * BYTES_THAT_NAME_THE_VERSION as u64,
            "one psp carries no census and three carry one of {} bytes each",
            census.len(),
        );
        assert_eq!(rows.len(), 4);
    }
}
