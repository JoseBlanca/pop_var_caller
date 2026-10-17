//! **Where a run's messages go: the terminal and the run's log file**
//! (`doc/devel/ng/spec/command_line.md` §7).
//!
//! Code anywhere in the crate states an event and its importance with `tracing` —
//! `info!("…")`, `debug!("…")`, `warn!("…")` — and says nothing about where it goes. The binary
//! calls [`install`] once, after parsing its command line, and that decides:
//!
//! | Level   | Terminal (stderr)          | Log file |
//! |---------|----------------------------|----------|
//! | `error` | yes, as `error: …`         | yes      |
//! | `warn`  | yes, as `warning: …`       | yes      |
//! | `info`  | yes, unless `--quiet`      | yes      |
//! | `debug` | no                         | yes      |
//! | `trace` | no                         | only with `RUST_LOG=trace` |
//!
//! **The terminal is short and the log is complete.** A subcommand that keeps a log opens it with
//! [`open_log`] before its first event worth keeping, and the binary closes it with [`close_log`]
//! after the run's last word, which on a failure is the error itself.
//!
//! **A progress line is an `info` event with the target [`PROGRESS`].** On a terminal it is drawn
//! in place, each update over the last; the next ordinary message leaves the last update standing
//! and starts below it. Redirected — a cluster job, `2> run.err` — each update is a line of its
//! own, with no control characters.

use std::fmt::{self, Write as _};
use std::fs::{File, OpenOptions};
use std::io::{self, IsTerminal, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};

use thiserror::Error;
use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::layer::{Context, Layer, SubscriberExt};
use tracing_subscriber::registry::Registry;

use super::provenance::{current_command_line, rfc3339_now};

/// The target a progress event carries, so the terminal can draw it in place:
/// `info!(target: PROGRESS, "calling: …")`.
pub const PROGRESS: &str = "pvc::progress";

/// The run's log, while one is open.
static LOG: Mutex<Option<RunLog>> = Mutex::new(None);

/// Every option of the subcommand that was run, as [`install`] was handed them.
static SETTINGS: OnceLock<Vec<Setting>> = OnceLock::new();

/// One option of the subcommand that was run, and the value the run used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Setting {
    /// The option as it is typed, without its dashes: `max-locus-span`.
    pub name: String,
    /// The value, as text; several values are separated by spaces. `None` where the option was
    /// neither typed nor has a default.
    pub value: Option<String>,
    /// Whether the value was typed, rather than taken from the option's default.
    pub typed: bool,
}

/// **Send every event to the terminal and the log**, from now until the process ends.
///
/// `quiet` keeps the terminal to warnings and errors. `settings` are the options of the
/// subcommand being run, which every log opened afterwards records in its header.
///
/// Called once, by the binary. **A second call changes nothing**: a process has one global
/// subscriber, and tests that run several subcommands in one process reach here through none of
/// them.
pub fn install(quiet: bool, settings: Vec<Setting>) {
    let _ = SETTINGS.set(settings);
    let terminal_level = if quiet { Level::WARN } else { Level::INFO };
    let log_level = log_level_from_the_environment();
    let terminal = TerminalLayer {
        most_detailed: terminal_level,
        interactive: io::stderr().is_terminal(),
        progress_line_open: AtomicBool::new(false),
    };
    let log = LogLayer {
        most_detailed: log_level,
    };
    // **One filter for the whole subscriber, at the more detailed of the two.** A layer that
    // filters by itself in `enabled` would disable the event for both; this filter only spares
    // the work of building events neither destination keeps.
    let subscriber = Registry::default()
        .with(LevelFilter::from_level(log_level.max(terminal_level)))
        .with(terminal)
        .with(log);
    let _ = tracing::subscriber::set_global_default(subscriber);
}

/// `debug`, or `trace` when `RUST_LOG` asks for it. Only the level is read: the log's content is
/// set by this module, and `RUST_LOG` exists for a developer who wants everything.
fn log_level_from_the_environment() -> Level {
    match std::env::var("RUST_LOG") {
        Ok(asked) if asked.eq_ignore_ascii_case("trace") => Level::TRACE,
        _ => Level::DEBUG,
    }
}

/// Where a log goes when `--log` does not say: `name` in the directory `output` is written to.
#[must_use]
pub fn beside(output: &Path, name: &str) -> PathBuf {
    match output.parent() {
        Some(directory) if !directory.as_os_str().is_empty() => directory.join(name),
        _ => PathBuf::from(name),
    }
}

/// Why a log could not be opened.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum LogFileError {
    /// A log is already there, and the run was not asked to replace it.
    #[error(
        "a log is already at {} — give --force to replace it together with the output, or \
         --log to write this run's log elsewhere",
        path.display()
    )]
    AlreadyThere {
        /// The log.
        path: PathBuf,
    },

    /// The file could not be created.
    #[error("the log {} could not be created", path.display())]
    NotCreated {
        /// The log.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: io::Error,
    },
}

/// **Open the run's log at `path` and write its header**: what ran, where, when, and every
/// setting's value, including the ones left at their default.
///
/// An existing file is refused unless `replace`, which a subcommand passes as its `--force`: the
/// log belongs to the output it sits beside and is replaced together with it. **A log whose run
/// failed is replaced without `--force`**: it belongs to no output, and refusing it would make
/// the corrected command fail on the mistake it corrects. A log that ends any other way — a run
/// that succeeded, one still running, one that was killed — is refused.
///
/// **The returned handle writes to this file and no other**, which is what a run uses for what it
/// writes itself — the parameters block. Events from anywhere in the crate go to whichever log is
/// open, which in the binary is only ever this one.
///
/// # Errors
///
/// [`LogFileError::AlreadyThere`] and [`LogFileError::NotCreated`].
pub fn open_log(path: &Path, replace: bool) -> Result<RunLog, LogFileError> {
    let mut options = OpenOptions::new();
    options.write(true);
    if replace {
        options.create(true).truncate(true);
    } else {
        // **Refused by the open itself**, not by a check before it, so two runs started at once
        // cannot both decide the file is free.
        options.create_new(true);
    }
    let opened = match options.open(path) {
        Err(source)
            if source.kind() == io::ErrorKind::AlreadyExists && left_by_a_failed_run(path) =>
        {
            OpenOptions::new().write(true).truncate(true).open(path)
        }
        opened => opened,
    };
    let file = opened.map_err(|source| {
        if source.kind() == io::ErrorKind::AlreadyExists {
            LogFileError::AlreadyThere {
                path: path.to_path_buf(),
            }
        } else {
            LogFileError::NotCreated {
                path: path.to_path_buf(),
                source,
            }
        }
    })?;
    let log = RunLog {
        file: Arc::new(Mutex::new(file)),
        path: path.to_path_buf(),
        opened: Instant::now(),
    };
    *lock_the_log() = Some(log.clone());
    write_the_header();
    Ok(log)
}

/// The words [`close_log`] ends a log with, before the exit status.
const FINISHED_WITH_EXIT_STATUS: &str = " finished: exit status ";

/// Whether the log at `path` ends with its run's failure: a last line from [`close_log`] with an
/// exit status other than 0.
fn left_by_a_failed_run(path: &Path) -> bool {
    let Ok(text) = std::fs::read_to_string(path) else {
        return false;
    };
    let Some(last) = text.lines().rev().find(|line| !line.trim().is_empty()) else {
        return false;
    };
    last.split_once(FINISHED_WITH_EXIT_STATUS)
        .is_some_and(|(_, status)| !status.starts_with("0,"))
}

/// The path of the log this run is writing, while one is open.
#[must_use]
pub fn log_path() -> Option<PathBuf> {
    lock_the_log().as_ref().map(|log| log.path.clone())
}

/// **Close the run's log with how the run ended**: the exit status, how long it took, and the
/// process's peak memory. Called by the binary after the run's last word — on a failure, after
/// the error, so the error is in the log too.
///
/// Nothing is lost if this is never called: every line is written to the file as it is logged.
pub fn close_log(exit_status: i32) {
    let Some(log) = lock_the_log().take() else {
        return;
    };
    let memory = match peak_resident_bytes() {
        Some(bytes) => format!("peak memory {}", size_in_bytes(bytes)),
        None => "peak memory not measured on this platform".to_owned(),
    };
    log.write_line(&format!(
        "{} DEBUG{FINISHED_WITH_EXIT_STATUS}{exit_status}, after {}, {memory}\n",
        rfc3339_now(),
        elapsed_text(log.opened.elapsed()),
    ));
}

/// **Close the run's log after a panic**: the panic's message as the log's last event, then how
/// the run ended, with the exit status a panicking Rust process has. Called from the binary's
/// panic hook, because a release build aborts on a panic and nothing after it runs — without
/// this the log would end at its last event, which reads as a run that was killed.
///
/// Written to the file only: the terminal gets Rust's own panic message.
pub fn close_log_after_a_panic(message: &str) {
    if let Some(log) = lock_the_log().clone() {
        log.write_line(&format!("{} ERROR {message}\n", rfc3339_now()));
    }
    close_log(101);
}

/// The log's lock. **A panic while it was held does not stop logging**: a line half-written by a
/// panicking thread is worth less than every line after it.
fn lock_the_log() -> std::sync::MutexGuard<'static, Option<RunLog>> {
    LOG.lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// What ran, where and when, then every setting — the first lines of every log.
fn write_the_header() {
    // **The binary's path, beside the version**, because the version alone does not tell two
    // development builds apart. A git commit would, but a build script that recorded it would make
    // every commit recompile the whole crate.
    match std::env::current_exe() {
        Ok(binary) => tracing::debug!(
            "pvc {}, the binary {}",
            env!("CARGO_PKG_VERSION"),
            binary.display()
        ),
        Err(_) => tracing::debug!("pvc {}", env!("CARGO_PKG_VERSION")),
    }
    tracing::debug!("command line: {}", current_command_line());
    match std::env::current_dir() {
        Ok(directory) => tracing::debug!("working directory: {}", directory.display()),
        Err(error) => tracing::debug!("working directory: unknown ({error})"),
    }
    tracing::debug!(
        "host: {}",
        rustix::system::uname().nodename().to_string_lossy()
    );
    tracing::debug!("threads: {}", rayon::current_num_threads());
    for setting in SETTINGS.get().map_or(&[][..], Vec::as_slice) {
        let value = setting.value.as_deref().unwrap_or("(not given)");
        let origin = if setting.typed || setting.value.is_none() {
            ""
        } else {
            " (default)"
        };
        tracing::debug!("setting: {} = {value}{origin}", setting.name);
    }
}

/// The first line of the parameters block in a log: a TOML comment, so the block from this line to
/// [`PARAMETERS_BLOCK_ENDS`] is a parameters file as it stands.
pub const PARAMETERS_BLOCK_BEGINS: &str = "# ---- parameters file begins ----";

/// The last line of the parameters block in a log.
pub const PARAMETERS_BLOCK_ENDS: &str = "# ---- parameters file ends ----";

/// **A run's log: the file, its path, and when it was opened.** Cloned into the global slot the
/// log destination writes events to, and held by the run for what it writes itself.
#[derive(Clone)]
pub struct RunLog {
    file: Arc<Mutex<File>>,
    path: PathBuf,
    opened: Instant,
}

impl RunLog {
    /// Where the log is.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// **The numbers the run scores with, as a parameters file**, between two marker lines that
    /// are TOML comments (`doc/devel/ng/spec/command_line.md` §7.3). Copied out of the log, the
    /// block is a file `--parameters` reads; [`parameters_in`] does the copying.
    pub fn record_parameters(&self, parameters_file_text: &str) {
        let mut block = format!(
            "{} DEBUG parameters: the numbers this run scores with, as a parameters file between \
             the two marker lines below; copy them out to pass them back with --parameters\n\
             {PARAMETERS_BLOCK_BEGINS}\n",
            rfc3339_now()
        );
        block.push_str(parameters_file_text);
        if !parameters_file_text.ends_with('\n') {
            block.push('\n');
        }
        block.push_str(PARAMETERS_BLOCK_ENDS);
        block.push('\n');
        self.write_line(&block);
    }

    /// Append `text` as it is. **One write**, so a run killed part-way leaves whole lines.
    fn write_line(&self, text: &str) {
        let mut file = self
            .file
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let _ = file.write_all(text.as_bytes());
    }
}

/// **The parameters block of a log, as a parameters file's text**: the lines between
/// [`PARAMETERS_BLOCK_BEGINS`] and [`PARAMETERS_BLOCK_ENDS`], or `None` when the log holds none.
#[must_use]
pub fn parameters_in(log_text: &str) -> Option<String> {
    // **Whole lines only**: the line that introduces the block may name it, and a marker matched
    // inside a sentence would start the block there.
    let mut lines = log_text.split_inclusive('\n');
    lines.find(|line| line.trim_end() == PARAMETERS_BLOCK_BEGINS)?;
    let mut block = String::new();
    for line in lines {
        if line.trim_end() == PARAMETERS_BLOCK_ENDS {
            return Some(block);
        }
        block.push_str(line);
    }
    None
}

/// The process's peak resident memory, in bytes — Linux's `VmHWM`, the high-water mark for the
/// life of the process; `None` elsewhere.
fn peak_resident_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let kib: u64 = status
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))?
        .trim()
        .strip_suffix("kB")?
        .trim()
        .parse()
        .ok()?;
    Some(kib * 1024)
}

/// A size in bytes as a person reads it: `3.2 GiB`, `850 MiB`.
fn size_in_bytes(bytes: u64) -> String {
    // Precision loss is irrelevant: a size to read, to one decimal.
    #[allow(clippy::cast_precision_loss)]
    let as_float = bytes as f64;
    if bytes >= 1 << 30 {
        format!("{:.1} GiB", as_float / f64::from(1u32 << 30))
    } else {
        format!("{:.0} MiB", as_float / f64::from(1u32 << 20))
    }
}

/// The terminal's destination: stderr, short.
struct TerminalLayer {
    most_detailed: Level,
    /// Whether stderr is a terminal, so a progress line can be drawn in place.
    interactive: bool,
    /// Whether the last thing written was a progress line with no line end after it.
    progress_line_open: AtomicBool,
}

impl<S: Subscriber> Layer<S> for TerminalLayer {
    fn on_event(&self, event: &Event<'_>, _context: Context<'_, S>) {
        let metadata = event.metadata();
        if *metadata.level() > self.most_detailed {
            return;
        }
        let message = message_of(event);
        let mut stderr = io::stderr().lock();
        if metadata.target() == PROGRESS && self.interactive {
            // `\r` back to the line's start, `ESC[2K` to clear what the last update left.
            let _ = write!(stderr, "\r\x1b[2K{message}");
            let _ = stderr.flush();
            self.progress_line_open.store(true, Ordering::Relaxed);
            return;
        }
        // **The last progress update stays on the screen** as the state the run had reached, and
        // the message starts on the line below it.
        if self.progress_line_open.swap(false, Ordering::Relaxed) {
            let _ = writeln!(stderr);
        }
        let _ = writeln!(stderr, "{}{message}", terminal_prefix(*metadata.level()));
    }
}

/// `error: ` and `warning: ` before the messages that need one; nothing before the rest.
fn terminal_prefix(level: Level) -> &'static str {
    match level {
        Level::ERROR => "error: ",
        Level::WARN => "warning: ",
        _ => "",
    }
}

/// The log's destination: the open log file, complete, each line stamped with its time.
struct LogLayer {
    most_detailed: Level,
}

impl<S: Subscriber> Layer<S> for LogLayer {
    fn on_event(&self, event: &Event<'_>, _context: Context<'_, S>) {
        let metadata = event.metadata();
        if *metadata.level() > self.most_detailed {
            return;
        }
        let Some(log) = lock_the_log().clone() else {
            return;
        };
        // **One write a line, unbuffered**, so a run that is killed leaves every line it logged.
        log.write_line(&format!(
            "{} {:<5} {}\n",
            rfc3339_now(),
            metadata.level(),
            message_of(event)
        ));
    }
}

/// An event's message, followed by any other fields it carries as `name=value`.
fn message_of(event: &Event<'_>) -> String {
    let mut collected = MessageAndFields::default();
    event.record(&mut collected);
    collected.text
}

/// Gathers an event's fields into one line of text.
#[derive(Default)]
struct MessageAndFields {
    text: String,
}

impl Visit for MessageAndFields {
    fn record_debug(&mut self, field: &Field, value: &dyn fmt::Debug) {
        if field.name() == "message" {
            // The message is the event's text, and `format_args!` prints it as written.
            let _ = write!(self.text, "{value:?}");
        } else {
            let _ = write!(self.text, " {}={value:?}", field.name());
        }
    }

    fn record_str(&mut self, field: &Field, value: &str) {
        if field.name() == "message" {
            self.text.push_str(value);
        } else {
            let _ = write!(self.text, " {}={value}", field.name());
        }
    }
}

/// `1234567` as `1,234,567`.
#[must_use]
pub fn with_thousands(number: u64) -> String {
    let digits = number.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

/// A duration as a person reads it: `1h12m`, `12m05s` or `42s` — to the second, since a
/// progress line comes every few minutes at the least.
#[must_use]
pub fn elapsed_text(took: Duration) -> String {
    let seconds = took.as_secs();
    let (hours, minutes, seconds) = (seconds / 3600, seconds / 60 % 60, seconds % 60);
    if hours > 0 {
        format!("{hours}h{minutes:02}m")
    } else if minutes > 0 {
        format!("{minutes}m{seconds:02}s")
    } else {
        format!("{seconds}s")
    }
}

/// A stretch of genome as a person reads its size: `782 Mb`, `12.4 Mb`, `8.0 kb`, `950 bases`.
#[must_use]
pub fn size_of_ground(bases: u64) -> String {
    // Precision loss is irrelevant: this is a size to read, rounded to three figures.
    #[allow(clippy::cast_precision_loss)]
    let as_float = bases as f64;
    match bases {
        0..1_000 => format!("{bases} bases"),
        1_000..1_000_000 => format!("{:.1} kb", as_float / 1e3),
        1_000_000..100_000_000 => format!("{:.1} Mb", as_float / 1e6),
        _ => format!("{:.0} Mb", as_float / 1e6),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_log_with_no_directory_goes_to_the_working_directory() {
        assert_eq!(
            beside(Path::new("calls.vcf.gz"), "pvc_calls.log"),
            PathBuf::from("pvc_calls.log")
        );
        assert_eq!(
            beside(Path::new("out/calls.vcf.gz"), "pvc_calls.log"),
            PathBuf::from("out/pvc_calls.log")
        );
    }

    #[test]
    fn ground_is_sized_in_the_unit_a_person_reads() {
        assert_eq!(size_of_ground(950), "950 bases");
        assert_eq!(size_of_ground(8_000), "8.0 kb");
        assert_eq!(size_of_ground(12_400_000), "12.4 Mb");
        assert_eq!(size_of_ground(782_000_000), "782 Mb");
    }

    #[test]
    fn only_warnings_and_errors_carry_a_prefix() {
        assert_eq!(terminal_prefix(Level::ERROR), "error: ");
        assert_eq!(terminal_prefix(Level::WARN), "warning: ");
        assert_eq!(terminal_prefix(Level::INFO), "");
    }

    /// The refusal comes from the open itself, so it names the file and leaves it untouched.
    #[test]
    fn an_existing_log_is_refused_unless_it_is_to_be_replaced() {
        std::fs::create_dir_all("tmp").unwrap();
        let directory = tempfile::tempdir_in("tmp").unwrap();
        let path = directory.path().join("pvc_calls.log");
        std::fs::write(&path, "an earlier run's log\n").unwrap();

        let refused = open_log(&path, false);
        assert!(matches!(refused, Err(LogFileError::AlreadyThere { .. })));
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "an earlier run's log\n"
        );
    }

    /// **A failed run's log is replaced without `--force`**; a finished or unfinished one is not.
    #[test]
    fn only_a_failed_runs_log_is_replaced_without_force() {
        std::fs::create_dir_all("tmp").unwrap();
        let directory = tempfile::tempdir_in("tmp").unwrap();
        let path = directory.path().join("pvc_calls.log");
        for (earlier, replaced) in [
            (
                "t DEBUG finished: exit status 1, after 3s, peak memory 1 MiB\n",
                true,
            ),
            (
                "t DEBUG finished: exit status 0, after 3s, peak memory 1 MiB\n",
                false,
            ),
            ("t INFO  calling: 40%\n", false),
        ] {
            std::fs::write(&path, earlier).unwrap();
            assert_eq!(left_by_a_failed_run(&path), replaced, "{earlier}");
        }
    }

    /// **The block a run writes is the text it was handed**, found again between its markers.
    #[test]
    fn the_parameters_block_reads_back_as_written() {
        std::fs::create_dir_all("tmp").unwrap();
        let directory = tempfile::tempdir_in("tmp").unwrap();
        let path = directory.path().join("pvc_calls.log");
        let log = open_log(&path, false).expect("a new log");
        let toml = "format = 3\n[reference]\nname = \"x\"\n";
        log.record_parameters(toml);

        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(parameters_in(&text).as_deref(), Some(toml));
        assert_eq!(parameters_in("a log with no block\n"), None);
    }

    #[test]
    fn sizes_are_read_in_mib_and_gib() {
        assert_eq!(size_in_bytes(850 << 20), "850 MiB");
        assert_eq!(size_in_bytes(3 << 30), "3.0 GiB");
    }
}
