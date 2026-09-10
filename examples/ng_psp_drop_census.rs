//! **Empty a psp's trailer, so that `regenerate-census` has something to rebuild** — the arming
//! step of plan step E1's whole-file oracle
//! (`doc/devel/ng/impl_plan/psp_census_pair.md`, step D4).
//!
//! ```text
//! cargo run --release --example ng_psp_drop_census -- <a.psp> [<b.psp> ...]
//! ```
//!
//! # Why this exists
//!
//! The oracle is: copy a walked psp, rebuild its census with `regenerate-census`, and the copy is
//! the walked file again **byte for byte, whole** — header, blocks, index, trailer, footer. On
//! fixtures that is `a_regenerated_psp_is_the_walked_one_byte_for_byte`; on real reads it is
//! `scripts/ng_fit_stage_end_to_end.sh`.
//!
//! **A copy handed straight to the command is skipped**, because a psp whose census is the one
//! this run would write is left alone (`psp_census_pair.md` §8) — so the comparison would be the
//! walk's own bytes against themselves and would pass whatever the command did. The copy has to
//! be owed a rebuild first, and the way to owe one is to have no census at all.
//!
//! **No shipped subcommand empties a trailer**, and a shell script has no business seeking into a
//! psp: the file has to be cut at the trailer's offset and a footer re-encoded behind the cut with
//! its trailer length set to zero, which is `replace_trailer`'s job and not `dd`'s. This is that
//! one library call with a command line around it. It is a test instrument, not part of the
//! pipeline.
//!
//! # What it prints, and when it fails
//!
//! One line per file — the psp's path and how many trailer bytes went — and it **exits 1 if a
//! psp's trailer was already empty**, naming it. That is not tidiness: a psp with no census is one
//! the command was going to rebuild anyway, so a caller that ignored it would go on to compare a
//! file that was never armed and call the pass an oracle.
//!
//! **A failure part-way through says which psps are already empty**, because they are, and because
//! `replace_trailer` can leave the one it failed on with no footer at all — a psp nothing will
//! open until it is written again.

use std::path::PathBuf;
use std::process::ExitCode;

use pop_var_caller::ng::psp::{PspReader, replace_trailer};

fn main() -> ExitCode {
    let psps: Vec<PathBuf> = std::env::args_os().skip(1).map(PathBuf::from).collect();
    if psps.is_empty() {
        eprintln!(
            "usage: ng_psp_drop_census <a.psp> [<b.psp> ...]\n\
             empties each psp's trailer, so that regenerate-census is owed a rebuild of it."
        );
        return ExitCode::from(2);
    }
    match drop_the_census_from(&psps) {
        Ok(()) => ExitCode::SUCCESS,
        Err(complaint) => {
            eprintln!("{complaint}");
            ExitCode::FAILURE
        }
    }
}

/// Empty each psp's trailer, reporting what each one held.
///
/// **The whole list is read before anything is written**, so that a run naming one unreadable file
/// leaves every psp as it found it rather than half of them emptied.
fn drop_the_census_from(psps: &[PathBuf]) -> Result<(), String> {
    let mut held = Vec::with_capacity(psps.len());
    for psp in psps {
        let bytes = PspReader::open(psp)
            .map_err(|source| format!("{} does not open as a psp: {source}", psp.display()))?
            .footer()
            .trailer_bytes;
        if bytes == 0 {
            return Err(format!(
                "{} already has an empty trailer, so emptying it would arm nothing — whatever \
                 rebuilt it next was going to be asked for that rebuild anyway",
                psp.display(),
            ));
        }
        held.push(bytes);
    }
    let mut emptied: Vec<&PathBuf> = Vec::with_capacity(psps.len());
    for (psp, bytes) in psps.iter().zip(held) {
        // **The failure's own words**, which say whether the file is as it was or cut short and
        // footerless — the difference between running this again and re-walking the sample.
        if let Err(failure) = replace_trailer(psp, b"") {
            return Err(format!(
                "{}: {failure} — {}{}",
                psp.display(),
                failure.source,
                already_emptied(&emptied),
            ));
        }
        emptied.push(psp);
        println!("{}: {bytes} bytes of census dropped", psp.display());
    }
    Ok(())
}

/// The psps this run already emptied, for a failure message to carry.
///
/// **Empty when there are none**, so a run that fails on its first psp does not claim to have
/// changed anything.
fn already_emptied(emptied: &[&PathBuf]) -> String {
    if emptied.is_empty() {
        return String::new();
    }
    let names: Vec<String> = emptied
        .iter()
        .map(|psp| psp.display().to_string())
        .collect();
    format!("; already emptied: {}", names.join(", "))
}
