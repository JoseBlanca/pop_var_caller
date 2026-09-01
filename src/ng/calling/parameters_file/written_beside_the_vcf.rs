//! **Where a run's parameters file goes, and the write that puts it there** —
//! spec §7, *"every run writes the parameters file it used, beside its VCF, whatever the numbers
//! came from"*.
//!
//! # Why this is a function and not three lines at a call site
//!
//! Spec §7 makes writing **unconditional**, on the stated grounds that *the run that most needs
//! its parameters recorded is the one whose operator did not think to ask*. A rule that holds
//! only where somebody remembered to write it is not unconditional, and the place a run driver
//! forgets is not the writing but the **naming**: a file whose name a run chooses ad hoc is a
//! file the next run overwrites, or one nobody can pair with its VCF.
//!
//! So the path is derived from the VCF's, once, here, and the writing is one call.
//!
//! # The naming rule, and the one thing it is careful about
//!
//! **The VCF's own name with its compression and format suffixes taken off, then
//! `.parameters.toml`.** `calls.vcf.gz` gives `calls.parameters.toml`; so does `calls.vcf`, and
//! `calls.bcf`; a path with none of those suffixes simply gains the new one, so `calls` gives
//! `calls.parameters.toml`.
//!
//! **Two VCFs with different names in one directory keep two parameters files**, which is the
//! reason the stem is used rather than a fixed name like `parameters.toml`. **Two spellings of
//! one stem do not**: `calls.vcf.gz` and `calls.bcf` are one cohort written twice and share
//! `calls.parameters.toml`, which is right where the two hold the same calls and wrong if they
//! ever do not.
//!
//! **What the rule does not give is a way to tell two files apart by their contents** — the file
//! records no run date, no caller version and no command line, so two copies in a directory are
//! distinguishable only by their names. That is recorded rather than solved: what a run stamps
//! into its own output is the command surface's, beside the rest of `pop_var_caller_exp`'s
//! subcommands (spec §11).
//!
//! **The argument is a VCF *file's* path.** Hand it a directory and the parameters file lands
//! beside that directory rather than in it; hand it a path with no file name at all (`/`, `.`,
//! the empty path) and the name is the hidden `.parameters.toml`. Neither is checked, because
//! neither is a state a run driver holding its own output path can be in.
//!
//! **Range (`CLAUDE.md`).** `to_toml` builds the whole file as one `String` before any of it is
//! written, and spec §9 prices the substitution-rate axis at up to 62 MB at 3,000 samples — C4
//! re-measured 185 bytes a row, putting it nearer 79 MB. So this is a single allocation of that
//! size, taken **after** the last locus is called, where §9's memory paragraph prices only the
//! read side. Nothing at the 63 accessions of the tomato cohort; it is what breaks first at the
//! top of the committed range, and it is recorded rather than fixed because streaming the writer
//! is a change to `to_toml`'s shape rather than to this one.

use std::ffi::OsStr;
use std::io;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use super::ParametersFile;

/// The suffix a parameters file always ends in.
const PARAMETERS_SUFFIX: &str = ".parameters.toml";

/// **Where the parameters file for `vcf` goes** — a pure function of the VCF's path, so that a
/// caller can name the file it is about to write, or find one a previous run wrote, without
/// writing anything.
///
/// See this module's header for the rule and what it buys.
#[must_use]
pub fn beside_the_vcf(vcf: &Path) -> PathBuf {
    // **Never through a `String`.** A file name on both platforms this project builds for is
    // arbitrary bytes, and `to_string_lossy` maps every invalid byte to U+FFFD — so
    // `/data/\xff.vcf` and `/data/\xfe.vcf` produced one parameters file between them, each
    // named after neither VCF, and the second run overwrote the first's. `file_stem` and
    // `extension` work on `OsStr`, so the name that comes out is the name that went in.
    //
    // **Compression first, then format**, because that is the order they are applied in:
    // `calls.vcf.gz` is a `calls.vcf` that was compressed.
    let mut stem = vcf;
    let uncompressed;
    if matches_any(stem.extension(), &["gz", "bgz"]) {
        uncompressed = stem.with_extension("");
        stem = &uncompressed;
    }
    let unformatted;
    if matches_any(stem.extension(), &["vcf", "bcf"]) {
        unformatted = stem.with_extension("");
        stem = &unformatted;
    }

    let mut name = stem.file_name().unwrap_or_default().to_os_string();
    name.push(PARAMETERS_SUFFIX);
    vcf.with_file_name(name)
}

/// Whether a path's extension is one of these, compared as bytes and case-sensitively — the way
/// the file systems this project runs on compare them.
fn matches_any(extension: Option<&OsStr>, any_of: &[&str]) -> bool {
    extension.is_some_and(|extension| any_of.iter().any(|candidate| extension == *candidate))
}

/// Where the bytes go before they are renamed into place — the destination's own name with
/// `.tmp` after it, which is the spelling the VCF's sink uses for the same purpose.
fn in_flight_path_for(at: &Path) -> PathBuf {
    let mut name = at.as_os_str().to_os_string();
    name.push(".tmp");
    PathBuf::from(name)
}

impl ParametersFile {
    /// **Write this file beside `vcf`**, and say where it went.
    ///
    /// The whole of spec §7 that can exist before a run driver does: the path rule and the write.
    /// **When it is called is the driver's**, and the one thing the driver has to decide is the
    /// order of the two writes — see [`Self::of_run`]'s own note on what a panic there would cost
    /// after the last locus has been called.
    ///
    /// **Written whole or not at all.** The text goes to a temporary file in the same directory
    /// and is renamed over the destination, which is atomic on both platforms this project builds
    /// for. `fs::write` truncates first, so a write that fails part-way — a full disk is the
    /// ordinary way — would leave a half-written parameters file beside a *complete* VCF, after
    /// every locus of the run had been called. A truncated file that still parses is the worse
    /// half of that: it looks like an answer.
    ///
    /// **It also replaces a file that may be this run's own input.** Spec §7 tells a user to copy
    /// the file their run wrote and edit a line, and a re-run whose supplied file and whose VCF
    /// share a stem writes over the file it was handed. Whether a driver should refuse that is
    /// the driver's; the rename at least makes the replacement whole.
    ///
    /// # The temporary is created the way the VCF's is, and the mode is the reason
    ///
    /// **`File::create` rather than a `NamedTempFile`, changed 2026-09-01.** A named temporary is
    /// created mode `0600` and keeps it through the rename, so under an ordinary `umask 022` a
    /// run wrote its VCF world-readable and the file saying what those calls rest on readable by
    /// nobody but the person who launched it — measured on a real run: `-rw-r--r--` beside
    /// `-rw-------`. On a shared directory that defeats §7 for everyone else in the group, which
    /// is most of who §7 is for. `File::create` is what the VCF's own sink uses
    /// ([`vcf::writer`](crate::ng::vcf::writer)), so the two files now get their mode from the
    /// same rule — the process's `umask` — rather than from which crate created them.
    ///
    /// **Nothing about being written whole changes**: the bytes still land under a sibling name
    /// and are renamed over the destination in one step. What is given up is a temporary that
    /// removes itself on a failed write, so a write that fails part-way leaves
    /// `<name>.tmp` beside the destination — the same thing the VCF's sink leaves, and visible
    /// rather than mistakable for the answer.
    ///
    /// # Errors
    ///
    /// [`io::Error`] from creating, writing or renaming the temporary file — there is nothing
    /// this can add to *permission denied on `/data/calls.parameters.toml`* that a caller does
    /// not already have.
    pub fn write_beside_the_vcf(&self, vcf: &Path) -> io::Result<PathBuf> {
        let at = beside_the_vcf(vcf);
        let in_flight = in_flight_path_for(&at);
        let mut written = std::fs::File::create(&in_flight)?;
        written.write_all(self.to_toml().as_bytes())?;
        written.sync_all()?;
        drop(written);
        std::fs::rename(&in_flight, &at)?;
        Ok(at)
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::a_file_using_every_shape;
    use super::beside_the_vcf;
    use std::path::Path;

    /// Every spelling of a VCF's name this project can produce, and the one name they all give.
    #[test]
    fn the_parameters_file_takes_the_vcfs_own_name() {
        for vcf in [
            "/data/run7/calls.vcf.gz",
            "/data/run7/calls.vcf",
            "/data/run7/calls.bcf",
            "/data/run7/calls.vcf.bgz",
            // No recognised suffix at all: the name gains the new one rather than losing a piece
            // of itself.
            "/data/run7/calls",
        ] {
            assert_eq!(
                beside_the_vcf(Path::new(vcf)),
                Path::new("/data/run7/calls.parameters.toml"),
                "from {vcf}"
            );
        }
    }

    /// **The directory is the VCF's, not the working directory**, which is the whole of *beside*
    /// — a run whose output goes somewhere else must not drop its parameters where it was
    /// launched from.
    #[test]
    fn it_lands_in_the_vcfs_own_directory() {
        assert_eq!(
            beside_the_vcf(Path::new("../elsewhere/tomato.vcf.gz")),
            Path::new("../elsewhere/tomato.parameters.toml")
        );
        assert_eq!(
            beside_the_vcf(Path::new("tomato.vcf.gz")),
            Path::new("tomato.parameters.toml")
        );
    }

    /// **Two cohorts called into one directory keep two parameters files.** A fixed name would
    /// have the second run overwrite the first's, silently, and leave a VCF beside parameters
    /// that are not its own.
    #[test]
    fn two_vcfs_in_one_directory_do_not_share_a_parameters_file() {
        assert_ne!(
            beside_the_vcf(Path::new("/data/tomato.vcf.gz")),
            beside_the_vcf(Path::new("/data/potato.vcf.gz"))
        );
    }

    /// **A name that merely contains the suffixes keeps them**, since only a trailing one is a
    /// suffix: `calls.vcf.backup` is not a VCF this rule should unwrap.
    #[test]
    fn only_a_trailing_suffix_is_stripped() {
        assert_eq!(
            beside_the_vcf(Path::new("/data/calls.vcf.backup")),
            Path::new("/data/calls.vcf.backup.parameters.toml")
        );
        assert_eq!(
            beside_the_vcf(Path::new("/data/gz.vcf.gz")),
            Path::new("/data/gz.parameters.toml")
        );
    }

    /// **Two VCFs whose names are different bytes get two parameters files**, even where neither
    /// name is text. Through `to_string_lossy` both collapsed to one U+FFFD name, so the second
    /// run overwrote the first's parameters and neither file was named after its VCF.
    #[cfg(unix)]
    #[test]
    fn two_names_that_are_not_text_do_not_collide() {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt as _;

        let one = Path::new(OsStr::from_bytes(b"/data/\xff\xfe.vcf"));
        let other = Path::new(OsStr::from_bytes(b"/data/\xfe\xff.vcf"));
        assert_ne!(beside_the_vcf(one), beside_the_vcf(other));
        assert_eq!(
            beside_the_vcf(one),
            Path::new(OsStr::from_bytes(b"/data/\xff\xfe.parameters.toml")),
            "the VCF's own bytes, with the suffix swapped"
        );
    }

    /// **A name that is only a suffix keeps it.** `.vcf` is a hidden file called `.vcf`, not an
    /// empty name with a `.vcf` extension, and `Path::extension` agrees — so it gains the new
    /// suffix rather than losing itself and becoming the hidden `.parameters.toml`.
    #[test]
    fn a_name_that_is_only_a_suffix_is_not_eaten() {
        assert_eq!(
            beside_the_vcf(Path::new("/data/.vcf")),
            Path::new("/data/.vcf.parameters.toml")
        );
    }

    /// The write puts the file where [`beside_the_vcf`] says, and what lands there reads back as
    /// the same file — so a run really is reproducible from its own output (spec §7).
    #[test]
    fn what_is_written_is_at_that_path_and_reads_back() {
        let directory = tempfile::tempdir().expect("a temporary directory");
        let vcf = directory.path().join("calls.vcf.gz");
        let file = a_file_using_every_shape();

        let at = file.write_beside_the_vcf(&vcf).expect("the write succeeds");

        assert_eq!(at, beside_the_vcf(&vcf));
        let text = std::fs::read_to_string(&at).expect("the file is there");
        assert_eq!(
            crate::ng::calling::parameters_file::ParametersFile::from_toml(&text)
                .expect("what a run wrote is what its reader reads"),
            file
        );
    }

    /// **The parameters file is as readable as the VCF beside it.**
    ///
    /// A named temporary is created mode `0600` and keeps it through the rename, so a run under
    /// an ordinary `umask 022` wrote its calls world-readable and the file saying what those
    /// calls rest on readable by nobody else — which defeats spec §7 for everyone in the group
    /// but the person who launched the run. Both files now take their mode from the process's
    /// `umask`, so the test is that the two agree rather than that either is a fixed number.
    #[cfg(unix)]
    #[test]
    fn the_parameters_file_is_as_readable_as_a_file_the_run_creates_beside_it() {
        use std::os::unix::fs::PermissionsExt as _;

        let directory = tempfile::tempdir().expect("a temporary directory");
        let vcf = directory.path().join("calls.vcf.gz");
        // What `File::create` gives under this process's umask — the same call the VCF's own
        // sink makes, so this is the mode the two files are being asked to share.
        std::fs::File::create(&vcf).expect("a stand-in for the VCF");
        let alongside = std::fs::metadata(&vcf)
            .expect("the stand-in is there")
            .permissions()
            .mode()
            & 0o777;

        let at = a_file_using_every_shape()
            .write_beside_the_vcf(&vcf)
            .expect("the write succeeds");

        let written = std::fs::metadata(&at)
            .expect("the parameters are there")
            .permissions()
            .mode()
            & 0o777;
        assert_eq!(
            written, alongside,
            "the parameters file is mode {written:o} and a file created beside it is {alongside:o}",
        );
        assert!(
            written & 0o044 != 0 || alongside & 0o044 == 0,
            "under a umask that lets others read the VCF, they can read the parameters too",
        );
    }

    /// **A directory that does not exist is an error and not a panic**, which is the state a
    /// mistyped `--out` reaches — the run has already called every locus by then.
    #[test]
    fn a_path_that_cannot_be_written_comes_back_as_an_error() {
        let directory = tempfile::tempdir().expect("a temporary directory");
        let vcf = directory.path().join("no-such-directory").join("calls.vcf");
        assert!(
            a_file_using_every_shape()
                .write_beside_the_vcf(&vcf)
                .is_err()
        );
    }
}
