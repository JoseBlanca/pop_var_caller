//! **Every parameter of the SNP/indel fit, after every pass** — a diagnostic, off unless
//! [`PATH_VARIABLE`] names a file.
//!
//! What it is for: telling a fit that is still converging from one that has stopped moving,
//! and seeing *which* parameter keeps it from meeting its stopping rule. One row per
//! parameter per pass, tab-separated:
//!
//! ```text
//! start  pass  log_likelihood  parameter  value
//! ```
//!
//! `log_likelihood` is the pass's own — the one the parameters *entering* the pass produce —
//! and `value` is where the maximisation after that pass left the parameter. Values are
//! printed with seventeen significant digits so a relative change of one in a million survives.
//! Nothing here changes what the fit computes.
//!
//! **Once a fit has chosen its best start, it writes what it returns**, under that start and
//! `pass` one after its last: one row a parameter with the value the fit returns — which is not
//! always the last pass's, since a start whose last accelerated step is refused at the pass limit
//! returns the step before it — and one row a parameter with its standard error, named
//! `standard_error:` and the parameter's name, NaN where it has none. Their `log_likelihood` is
//! the final pass's, at the returned values. So each pass's move can be read in units of the
//! parameter's own error, against the values the fit returned.
//!
//! The file is one a process: a process that fits twice writes both fits' rows under the same
//! start numbers. `estimate-parameters` fits once.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

/// The variable naming the file the trace is written to.
pub const PATH_VARIABLE: &str = "PVC_JOINT_FIT_TRACE";

/// The open file, or `None` where the variable was not set. Resolved once per process.
fn sink() -> Option<&'static Mutex<BufWriter<File>>> {
    static SINK: OnceLock<Option<Mutex<BufWriter<File>>>> = OnceLock::new();
    SINK.get_or_init(|| {
        let path = std::env::var_os(PATH_VARIABLE)?;
        let mut file = BufWriter::new(File::create(&path).unwrap_or_else(|why| {
            panic!(
                "{PATH_VARIABLE} names {} and it could not be created: {why}",
                std::path::Path::new(&path).display()
            )
        }));
        writeln!(file, "start\tpass\tlog_likelihood\tparameter\tvalue")
            .expect("the trace header could not be written");
        Some(Mutex::new(file))
    })
    .as_ref()
}

/// One row of the trace as a test captures it: start, pass, log-likelihood, parameter, value.
#[cfg(test)]
pub(super) type TraceRow = (usize, u32, f64, String, f64);

#[cfg(test)]
thread_local! {
    /// The rows this thread's fit writes, while a test captures them ([`captured`]).
    static CAPTURED: std::cell::RefCell<Option<Vec<TraceRow>>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `fit` with this thread's trace captured instead of written, and return its rows with its
/// result. The fit writes its rows from the thread that calls it — each pass's from the
/// alternation, the returned values' once it has chosen — so a capture on that thread sees them all.
#[cfg(test)]
pub(super) fn captured<R>(fit: impl FnOnce() -> R) -> (R, Vec<TraceRow>) {
    CAPTURED.with(|rows| *rows.borrow_mut() = Some(Vec::new()));
    let result = fit();
    let rows = CAPTURED.with(|rows| rows.borrow_mut().take().unwrap_or_default());
    (result, rows)
}

/// Whether a trace is being written — asked before the rows are built, so a run without one
/// builds nothing.
pub(super) fn is_on() -> bool {
    #[cfg(test)]
    if CAPTURED.with(|rows| rows.borrow().is_some()) {
        return true;
    }
    sink().is_some()
}

/// Write one pass's rows: every `(name, value)` of the parameters as the pass left them.
pub(super) fn write_the_pass(
    start: usize,
    pass: u32,
    log_likelihood: f64,
    parameters: &[(String, f64)],
) {
    #[cfg(test)]
    if CAPTURED
        .with(|rows| {
            rows.borrow_mut().as_mut().map(|rows| {
                rows.extend(
                    parameters
                        .iter()
                        .map(|(name, value)| (start, pass, log_likelihood, name.clone(), *value)),
                );
            })
        })
        .is_some()
    {
        return;
    }
    let Some(sink) = sink() else {
        return;
    };
    let mut file = sink
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    for (name, value) in parameters {
        writeln!(
            file,
            "{start}\t{pass}\t{log_likelihood:.17e}\t{name}\t{value:.17e}"
        )
        .expect("a trace row could not be written");
    }
    file.flush().expect("the trace could not be flushed");
}
