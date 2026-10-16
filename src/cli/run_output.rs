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
use std::sync::{Mutex, OnceLock};
use std::time::Duration;

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
static LOG: Mutex<Option<OpenLog>> = Mutex::new(None);

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
/// log belongs to the output it sits beside and is replaced together with it.
///
/// # Errors
///
/// [`LogFileError::AlreadyThere`] and [`LogFileError::NotCreated`].
pub fn open_log(path: &Path, replace: bool) -> Result<(), LogFileError> {
    let mut options = OpenOptions::new();
    options.write(true);
    if replace {
        options.create(true).truncate(true);
    } else {
        // **Refused by the open itself**, not by a check before it, so two runs started at once
        // cannot both decide the file is free.
        options.create_new(true);
    }
    let file = options.open(path).map_err(|source| {
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
    *lock_the_log() = Some(OpenLog {
        file,
        path: path.to_path_buf(),
    });
    write_the_header();
    Ok(())
}

/// The path of the log this run is writing, while one is open.
#[must_use]
pub fn log_path() -> Option<PathBuf> {
    lock_the_log().as_ref().map(|log| log.path.clone())
}

/// **Close the run's log**, after the run's last word. Nothing is lost if this is never called:
/// every line is written to the file as it is logged.
pub fn close_log() {
    if let Some(mut log) = lock_the_log().take() {
        let _ = log.file.flush();
    }
}

/// The log's lock. **A panic while it was held does not stop logging**: a line half-written by a
/// panicking thread is worth less than every line after it.
fn lock_the_log() -> std::sync::MutexGuard<'static, Option<OpenLog>> {
    LOG.lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// What ran, where and when, then every setting — the first lines of every log.
fn write_the_header() {
    tracing::debug!("pvc {}", env!("CARGO_PKG_VERSION"));
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

/// An open log file and where it is.
struct OpenLog {
    file: File,
    path: PathBuf,
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
        let mut log = lock_the_log();
        let Some(open) = log.as_mut() else {
            return;
        };
        let line = format!(
            "{} {:<5} {}\n",
            rfc3339_now(),
            metadata.level(),
            message_of(event)
        );
        // **One write a line, unbuffered**, so a run that is killed leaves every line it logged.
        let _ = open.file.write_all(line.as_bytes());
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
}
