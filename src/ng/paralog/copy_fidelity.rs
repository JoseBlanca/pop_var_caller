//! **These files are still production's, verbatim — asserted, not claimed.**
//!
//! Five files are guarded today, against four production originals:
//!
//! | ng's copy | lines | production's original |
//! |---|---|---|
//! | `coverage_model.rs` | 1,158 past its module header | `src/paralog/coverage_model.rs` |
//! | `locus_score.rs` | 797 | `src/paralog/locus_score.rs` |
//! | `prior.rs` | 524 | `src/paralog/prior.rs` |
//! | `model_params.rs` | 236 | the items of `src/paralog/mod.rs` |
//! | `calibration.rs` | 83 | a **span** of `src/var_calling/paralog_filter/calibrate.rs` |
//!
//! `locus_score.rs` differs from its original on three declared lines and `calibration.rs`
//! on seven; the other three differ on none. ng owns this file outright. It exists so the copies can
//! be checked from outside themselves — a check written *inside* `coverage_model.rs`
//! would be lines production's file does not have, so writing it would break the very
//! identity it asserts.
//!
//! # Why a textual check and not a green suite
//!
//! Every copy here arrives with production's own tests transcribed alongside it, so both
//! trees stay green through a divergence — that is exactly what makes the tests useless
//! as a fidelity statement. A `diff` is the cheapest and most durable one available, and
//! it survives a future re-copy. The port's whole claim is that it agrees with what it
//! was ported from (`doc/devel/ng/spec/hidden_paralog_filter.md` §1.2: *no new
//! statistics*), and this is what says so.
//!
//! # Where each copy's production content begins
//!
//! [`CopyBegins`] names the two shapes. Most copies are whole files, and ng appends a
//! note to the end of production's module header saying the file is a copy; production's
//! header must survive **verbatim, line for line, at the top**, or a rewritten header
//! would pass while a whole section went missing. `model_params.rs` is the other shape:
//! production keeps those items in its `mod.rs`, whose declarations ng cannot share, so
//! the comparison starts at the first item's `///` on each side instead.
//!
//! Past that point the two sides are compared **line for line, and then byte for byte**.
//! The second comparison is not redundant: `str::lines()` strips a trailing `\r` and
//! cannot see a missing final newline, so a line-ending normalisation on a re-copy would
//! leave two files `git` reports as different and a line compare reports as identical.
//!
//! # Where a span copy ends, and the one thing that is normalised
//!
//! `calibration.rs` copies a *span* of a larger file rather than the whole of it, and stops
//! at a declared marker — the first line **past** the copy. **One boundary is normalised,
//! and it is the only exception to "byte for byte" anywhere in this file**: production
//! separates its items with blank lines, so a span ends with one, and `cargo fmt` strips a
//! trailing blank line from ng's copy. A span's end is therefore trimmed on both sides to no
//! trailing blank lines and one final newline. It is a boundary rather than content, and a
//! whole-file copy is not normalised at all.
//!
//! A declared marker that production no longer has is **the guard failing**, not a passing
//! comparison: without that the extraction would run to the end of production's file and the
//! comparison would fail somewhere unrelated, blaming the copy.
//!
//! # The lines a copy may change: three sanctioned kinds
//!
//! A copy that names `crate::paralog::` reaches back into production from inside ng — at
//! compile time if it is a `use`, and at *read* time if it is a rustdoc link, which sends a
//! reader of ng's own API documentation into the frozen tree. `locus_score.rs` has three
//! such lines: one import in its transcribed test module, and two link definitions naming
//! `SingleCopyCoverageModel`, of which **one is inside production's module header**. Each is
//! declared as a [`Repoint`] and the substitutions are applied to the whole of production's
//! file, header included, so a header line can be repointed without breaking the rule that
//! production's header survives verbatim.
//!
//! `calibration.rs` has **seven** such lines, and two of them are the other two kinds that
//! [`WhyRepointed`] names:
//!
//! - a **path into production**, as above;
//! - a **rustdoc link to a production item ng does not copy** — `[`calibrate`]`, the
//!   spill-streaming driver ng leaves behind — which would otherwise resolve to nothing;
//! - a **`pub(crate)` widened to `pub`**, on five items. ng re-exports them from
//!   `src/ng/paralog/mod.rs`, so that `crate::ng::paralog::ParalogCalibration` resolves the
//!   way `crate::ng::paralog::ParalogPrior` does, and a `pub use` cannot re-export a
//!   `pub(crate)` item. **This is ng widening its own copy, not production** — a different
//!   rule from the one `src/ng/mod.rs` states about widening a *production* item so a parity
//!   test can see it, and with a different subject.
//!
//! **Every rule is exact rather than loose.** Production's side must appear **exactly once**,
//! ng's side must be **exactly** the declared replacement at that same line, and the
//! indentation and the line ending must match. Each kind further requires that ng's line be
//! production's line with **one named substring replaced and nothing else changed** — a
//! check on the declaration itself, which the content comparison cannot make, because it
//! compares production's file against that very declaration. So a *new* path edit, a dropped
//! declaration, and a declaration that quietly changes a constant all fail.
//!
//! # Releasing a file
//!
//! A file this port deliberately changes is **deleted from `guarded_copies` in the commit
//! that changes it**, not commented out, and the release is recorded in the table below
//! and in that file's own module header. Switching the guard off at the first divergence
//! would throw away a working check on every file that is still a copy.
//!
//! | released at | file | why |
//! |---|---|---|
//! | — | — | nothing released yet |
//!
//! # What is deliberately *not* copied
//!
//! One thing, and it is left out on purpose rather than missing. `calibration.rs`'s span
//! stops before production's `cohort_inbreeding`, the helper that fills a slice with one
//! inbreeding coefficient repeated for every sample. Production settled on that because the per-individual proxy it had was
//! contaminated by divergence; ng's parameters file already carries a coefficient fitted per
//! sample, and the copied scorer takes it as a per-sample slice
//! (`doc/devel/ng/spec/hidden_paralog_filter.md` §3.2). It is the one piece of the filter ng
//! **replaces** rather than ports.
//!
//! [`the_item_the_span_stops_before_is_still_productions`] pins that: the span's end marker
//! is that helper's doc comment, so a rename would otherwise leave the guard green while the
//! copy silently grew.
//!
//! **A deleted entry and a dropped one look the same**, which is why
//! [`every_file_in_the_module_is_guarded_ng_s_own_or_released`] makes the directory
//! listing the authority: a file in neither list fails the build, so the guard cannot go
//! quietly empty.
//!
//! **One release is already scheduled.** `coverage_model.rs` fits from production's
//! `CoverageByGcHistogram` because ng's own copy of that type is on the unmerged
//! `ng-window-coverage` branch. When it lands, the `use` line changes and this test fails
//! — which is the point: the swap becomes a commit that says what it did, rather than a
//! line nobody noticed.

/// Where a copied file's production content begins on each side.
#[derive(Clone, Copy)]
enum CopyBegins {
    /// After the leading `//!` run. ng appends its note to production's header, so the
    /// header itself is compared as a prefix and the note is the one permitted addition.
    AfterTheModuleHeader,
    /// At the first item's `///` doc comment — for a copy whose production original is a
    /// `mod.rs` whose declarations ng does not share.
    AtTheFirstItemDoc,
}

/// Why a copy is allowed to change one line. The mechanism sanctions *these* three reasons
/// and nothing else, and each check is exact, so it cannot be used to wave through a changed
/// constant.
#[derive(Clone, Copy, PartialEq)]
enum WhyRepointed {
    /// A path that would otherwise name production from inside ng — at compile time if it is
    /// a `use`, at read time if it is a rustdoc link. ng's line must name `crate::ng::paralog`
    /// where production's names `crate::paralog`.
    PathIntoProduction,
    /// A documentation reference to a production item ng deliberately does not copy, which
    /// would otherwise be a link that resolves to nothing. `link` is the rustdoc link
    /// production's line carries; `replaced_by` is the plain words ng puts in its place. The
    /// rest of the line must be untouched.
    DocLinkToAnItemNgDoesNotCopy {
        link: &'static str,
        replaced_by: &'static str,
    },
    /// A `pub(crate)` item that has to be `pub` in ng, because production keeps it inside a
    /// private module and ng's `paralog` is a public one whose consumers — the run wiring at
    /// plan step B1 — sit outside it. **The two lines must be identical but for that one
    /// keyword**: nothing else about the item may move, and no other visibility change is
    /// sanctioned.
    VisibilityNgsPublicModuleNeeds,
}

/// One line a copy is allowed to change: production's text, the text ng puts in its place,
/// and which of the three sanctioned reasons it is.
struct Repoint {
    production_line: &'static str,
    ng_line: &'static str,
    why: WhyRepointed,
}

/// One file ng copied from `src/paralog/`, paired with the original it must equal.
///
/// A named struct rather than a tuple because two of its fields are same-typed file
/// contents: transposed, a tuple still *passes* — line equality is symmetric — and the
/// mistake surfaces only on the day the guard fires, with every message naming the wrong
/// side.
struct GuardedCopy {
    /// The copy's file name, used to name it in a failure message.
    file_name: &'static str,
    /// The original's path, so a message names the file that must not be edited rather than
    /// a directory. Not every copy comes from `src/paralog/`: `calibration.rs` takes a span
    /// of `src/var_calling/paralog_filter/calibrate.rs`.
    production_path: &'static str,
    /// The whole of the production original, read at compile time.
    production_source: &'static str,
    /// The whole of ng's copy, read at compile time.
    ng_source: &'static str,
    /// Where the copied content starts on each side.
    begins: CopyBegins,
    /// The first line **past** the copy in production's original, for a copy that takes a
    /// span of a larger file. `None` where the copy runs to the end. Compared on the trimmed
    /// line, and a marker that never matches is a guard failure, not a pass.
    ends_before: Option<&'static str>,
    /// The sanctioned substitutions, each of which must occur exactly once.
    repoints: &'static [Repoint],
}

/// Every file ng copied into this module, with the original it must equal.
///
/// `include_str!` resolves relative to this file, so both sides are pinned at compile
/// time: a deleted or moved file is a build error, not a silently skipped case.
fn guarded_copies() -> Vec<GuardedCopy> {
    vec![
        GuardedCopy {
            file_name: "coverage_model.rs",
            production_path: "src/paralog/coverage_model.rs",
            production_source: include_str!("../../paralog/coverage_model.rs"),
            ng_source: include_str!("coverage_model.rs"),
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            repoints: &[],
        },
        GuardedCopy {
            file_name: "locus_score.rs",
            production_path: "src/paralog/locus_score.rs",
            production_source: include_str!("../../paralog/locus_score.rs"),
            ng_source: include_str!("locus_score.rs"),
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            // Three paths into production. The import would build ng's test fixture from
            // production's types (and does not compile, since ng's are distinct); the two
            // link definitions would send a reader of ng's own docs into the frozen tree.
            // The first of those two is a `//!` line inside production's module header,
            // which is why substitutions are applied to the whole file rather than only to
            // the content past the header.
            repoints: &[
                Repoint {
                    production_line: "use crate::paralog::{GridSpec, SfsPriorSpec};",
                    ng_line: "use crate::ng::paralog::{GridSpec, SfsPriorSpec};",
                    why: WhyRepointed::PathIntoProduction,
                },
                Repoint {
                    production_line: "//! [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel",
                    ng_line: "//! [`SingleCopyCoverageModel`]: crate::ng::paralog::SingleCopyCoverageModel",
                    why: WhyRepointed::PathIntoProduction,
                },
                Repoint {
                    production_line: "/// [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel",
                    ng_line: "/// [`SingleCopyCoverageModel`]: crate::ng::paralog::SingleCopyCoverageModel",
                    why: WhyRepointed::PathIntoProduction,
                },
            ],
        },
        GuardedCopy {
            file_name: "model_params.rs",
            production_path: "src/paralog/mod.rs",
            production_source: include_str!("../../paralog/mod.rs"),
            ng_source: include_str!("model_params.rs"),
            begins: CopyBegins::AtTheFirstItemDoc,
            ends_before: None,
            repoints: &[],
        },
        GuardedCopy {
            file_name: "prior.rs",
            production_path: "src/paralog/prior.rs",
            production_source: include_str!("../../paralog/prior.rs"),
            ng_source: include_str!("prior.rs"),
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            repoints: &[],
        },
        // **The one span copy.** Production keeps the calibration in the same file as the
        // spill-streaming driver ng does not want, so this takes that file from its first
        // item down to but not including the cohort inbreeding coefficient — the type ng
        // replaces with a per-sample one rather than porting.
        GuardedCopy {
            file_name: "calibration.rs",
            production_path: "src/var_calling/paralog_filter/calibrate.rs",
            production_source: include_str!("../../var_calling/paralog_filter/calibrate.rs"),
            ng_source: include_str!("calibration.rs"),
            begins: CopyBegins::AtTheFirstItemDoc,
            ends_before: Some(
                "/// The cohort inbreeding coefficient `F`, one value for every sample.",
            ),
            repoints: &[
                Repoint {
                    production_line: "/// match [`crate::paralog::prior::DEFAULT_EM_START`] (`0.03`); they are",
                    ng_line: "/// match [`crate::ng::paralog::prior::DEFAULT_EM_START`] (`0.03`); they are",
                    why: WhyRepointed::PathIntoProduction,
                },
                Repoint {
                    production_line: "/// Tuning knobs for [`calibrate`]. [`Default`] is the production configuration.",
                    ng_line: "/// Tuning knobs for the calibration. [`Default`] is the production configuration.",
                    why: WhyRepointed::DocLinkToAnItemNgDoesNotCopy {
                        link: "[`calibrate`]",
                        replaced_by: "the calibration",
                    },
                },
                // Production keeps these inside a private module; ng's `paralog` is public
                // and the run wiring that will use them (plan step B1) sits outside it, so
                // `pub(crate)` here would read as dead code and the module would export a
                // calibration nothing can name.
                Repoint {
                    production_line: "pub(crate) const DEFAULT_FALLBACK_PARALOG_PRIOR: f64 = 0.03;",
                    ng_line: "pub const DEFAULT_FALLBACK_PARALOG_PRIOR: f64 = 0.03;",
                    why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
                },
                Repoint {
                    production_line: "pub(crate) struct CalibrationConfig {",
                    ng_line: "pub struct CalibrationConfig {",
                    why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
                },
                Repoint {
                    production_line: "pub(crate) struct ParalogCalibration {",
                    ng_line: "pub struct ParalogCalibration {",
                    why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
                },
                Repoint {
                    production_line: "pub(crate) fn flags(&self, lr: f64) -> bool {",
                    ng_line: "pub fn flags(&self, lr: f64) -> bool {",
                    why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
                },
                Repoint {
                    production_line: "pub(crate) fn posterior(&self, lr: f64) -> Option<f64> {",
                    ng_line: "pub fn posterior(&self, lr: f64) -> Option<f64> {",
                    why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
                },
            ],
        },
    ]
}

/// How many leading lines to drop before the copied content starts.
fn lines_before_the_copy(source: &str, begins: CopyBegins) -> usize {
    match begins {
        CopyBegins::AfterTheModuleHeader => {
            source.lines().take_while(|l| l.starts_with("//!")).count()
        }
        CopyBegins::AtTheFirstItemDoc => source
            .lines()
            .take_while(|l| !l.starts_with("/// "))
            .count(),
    }
}

/// The same text with any **whole blank lines** at its end dropped, and one final newline.
///
/// The span's end boundary, and the only place either side is normalised. Production
/// separates its items with blank lines, so a span ends with one; ng's copy cannot keep it,
/// because `cargo fmt` strips a trailing blank line and the copy would then fail on every
/// run. **Only whole blank lines go** — trailing spaces or a `\r` on the last *content* line
/// are kept, because those are exactly the byte-level differences the comparison exists to
/// catch.
fn without_trailing_blank_lines(text: &str) -> String {
    let mut lines: Vec<&str> = text.split_inclusive('\n').collect();
    while lines.last().is_some_and(|line| line.trim().is_empty()) {
        lines.pop();
    }
    let mut kept: String = lines.concat();
    if !kept.ends_with('\n') {
        kept.push('\n');
    }
    kept
}

/// The copied content itself, with its line endings and its final newline intact.
///
/// `ends_before` stops the extraction at the first line whose trimmed text equals it — the
/// first line *past* the copy. A marker that is declared and never found is a guard failure
/// and is reported as one; without that, a renamed marker would silently extend the copy to
/// the end of production's file, and the comparison would fail somewhere unrelated.
fn the_copied_content(
    source: &str,
    begins: CopyBegins,
    ends_before: Option<&str>,
) -> Result<String, ()> {
    let after_the_start = source
        .split_inclusive('\n')
        .skip(lines_before_the_copy(source, begins));
    let Some(marker) = ends_before else {
        return Ok(after_the_start.collect());
    };
    // **The marker must occur exactly once past the copy's start.** A marker that comes to
    // match earlier — a doc sentence production repeats, say — would truncate the span and
    // the comparison would then blame ng's copy for lines it never claimed to have.
    let occurrences = source
        .split_inclusive('\n')
        .skip(lines_before_the_copy(source, begins))
        .filter(|line| line.trim() == marker)
        .count();
    if occurrences != 1 {
        return Err(());
    }
    let after_the_start = source
        .split_inclusive('\n')
        .skip(lines_before_the_copy(source, begins));
    let mut content = String::new();
    for line in after_the_start {
        if line.trim() == marker {
            // **The blank lines between the span and what follows belong to neither.**
            // Production's file separates its items with them; ng's copy cannot keep one at
            // its end, because `cargo fmt` strips a trailing blank line and the copy would
            // then fail the guard on every run. So the span's end is normalised to no blank
            // lines and one final newline — the only place either side is normalised at all,
            // and it is a boundary rather than content.
            return Ok(without_trailing_blank_lines(&content));
        }
        content.push_str(line);
    }
    Err(())
}

/// **Production's whole file with each sanctioned path substitution applied**, so that every
/// comparison below — the module-header prefix included — can stay exact.
///
/// Applied to the whole file rather than to the copied content alone because one of
/// `locus_score.rs`'s three repointed lines is a `//!` line *inside* production's module
/// header, and the header is compared as a byte-verbatim prefix.
///
/// The indentation and the line ending are taken from production's own line, so a copy that
/// re-indented an import or changed its ending still fails. `Err` names a substitution that
/// did not occur exactly once — the guard reporting that *it* is broken, not that the copy
/// drifted.
fn production_with_repoints_applied(copy: &GuardedCopy) -> Result<String, String> {
    for repoint in copy.repoints {
        // **Every kind is exact: ng's line must be production's with one named substring
        // replaced.** Testing only that the two sides *contain* something would sanction any
        // other difference on the same line — a doc comment quoting a constant production
        // does not have, or an extra item on a `use` line — and the content comparison
        // cannot cross-check it, because it compares against the declaration itself.
        let sanctioned = match repoint.why {
            WhyRepointed::PathIntoProduction => {
                repoint.production_line.contains("crate::paralog")
                    && repoint
                        .production_line
                        .replace("crate::paralog", "crate::ng::paralog")
                        == repoint.ng_line
            }
            WhyRepointed::DocLinkToAnItemNgDoesNotCopy { link, replaced_by } => {
                link.starts_with("[`")
                    && !replaced_by.contains('[')
                    && repoint.production_line.contains(link)
                    && repoint.production_line.replace(link, replaced_by) == repoint.ng_line
            }
            WhyRepointed::VisibilityNgsPublicModuleNeeds => {
                repoint.production_line.contains("pub(crate)")
                    && repoint.production_line.replace("pub(crate)", "pub") == repoint.ng_line
            }
        };
        if !sanctioned {
            return Err(format!(
                "{}: THE GUARD COULD NOT RUN — `{}` → `{}` is not one of the three \
                 sanctioned kinds of substitution. A repoint turns a `crate::paralog` path \
                 into a `crate::ng::paralog` one, replaces a documentation link to a \
                 production item ng does not copy with plain words, or widens a `pub(crate)` \
                 to a `pub` — and changes nothing else on the line. It may not change what \
                 the copy computes. Fix the declaration in copy_fidelity.rs; ng's copy is \
                 not the file at fault",
                copy.file_name, repoint.production_line, repoint.ng_line,
            ));
        }
    }
    let content = copy.production_source;
    let mut rewritten = String::with_capacity(content.len());
    let mut applied = vec![0usize; copy.repoints.len()];
    for line in content.split_inclusive('\n') {
        let body = line.trim_end_matches(['\n', '\r']);
        let ending = &line[body.len()..];
        let indent_len = body.len() - body.trim_start().len();
        match copy
            .repoints
            .iter()
            .position(|r| body.trim_start() == r.production_line)
        {
            Some(index) => {
                applied[index] += 1;
                rewritten.push_str(&body[..indent_len]);
                rewritten.push_str(copy.repoints[index].ng_line);
                rewritten.push_str(ending);
            }
            None => rewritten.push_str(line),
        }
    }
    for (index, times) in applied.iter().enumerate() {
        if *times != 1 {
            return Err(format!(
                "{}: THE GUARD COULD NOT RUN — the sanctioned repoint of `{}` occurs {} \
                 times in the production original, not once. Either production changed \
                 the line, or a second occurrence appeared. Fix the `Repoint` \
                 declaration in copy_fidelity.rs; do not touch ng's copy until it is \
                 right, because until then nothing is checking it",
                copy.file_name, copy.repoints[index].production_line, times
            ));
        }
    }
    Ok(rewritten)
}

/// **How ng's copy differs from production's original, or `None` if it does not.**
///
/// `Some(_)` is the whole message a failing guard prints. It is extracted from the
/// `#[test]` so that each way of rejecting can be driven with synthetic strings
/// (`the_comparison_accepts_an_appended_note_and_rejects_every_other_edit`): with only the
/// real, currently-identical file pair to run on, every rejecting branch would be
/// exercised in its accepting direction alone, and a refactor that turned one of them into
/// `assert!(true)` would leave the suite green.
fn how_the_copy_differs(copy: &GuardedCopy) -> Option<String> {
    let GuardedCopy {
        file_name,
        production_path,
        production_source: _,
        ng_source,
        begins,
        ends_before,
        repoints: _,
    } = *copy;

    // Every comparison below runs against production's file with the sanctioned
    // substitutions already made, header included.
    let repointed = match production_with_repoints_applied(copy) {
        Ok(source) => source,
        Err(the_guard_is_broken) => return Some(the_guard_is_broken),
    };
    let production_source = repointed.as_str();

    if let CopyBegins::AfterTheModuleHeader = begins {
        let production_header = lines_before_the_copy(production_source, begins);
        let ng_header = lines_before_the_copy(ng_source, begins);
        if ng_header < production_header {
            return Some(format!(
                "{file_name}: ng's module header is {ng_header} lines to production's \
                 {production_header} — ng appends its note to production's header, it does \
                 not write over it"
            ));
        }
        for (index, (ng_line, production_line)) in ng_source
            .split_inclusive('\n')
            .zip(production_source.split_inclusive('\n'))
            .take(production_header)
            .enumerate()
        {
            if ng_line != production_line {
                return Some(format!(
                    "{file_name}, module-header line {}: production's header must survive \
                     verbatim; ng's note is appended to it, not written over it.\n  \
                     ng's copy:  {ng_line:?}\n  the original: {production_line:?}",
                    index + 1
                ));
            }
        }
    }

    // A copy that takes a span needs its end marker found on **both** sides; ng's own copy
    // ends where the span ends, so the marker is absent there and the extraction runs to the
    // file's end, which is right.
    let Ok(production_body) = the_copied_content(production_source, begins, ends_before) else {
        return Some(format!(
            "{file_name}: THE GUARD COULD NOT RUN — the line this copy is declared to end \
             before, `{}`, is not in {production_path} any more. Without it the copy would \
             silently extend to the end of production's file. Fix `ends_before` in \
             copy_fidelity.rs; nothing is checking ng's copy until it is right",
            ends_before.unwrap_or_default(),
        ));
    };
    let ng_body = {
        let whole = the_copied_content(ng_source, begins, None)
            .expect("no end marker was asked for, so extraction cannot fail");
        // Normalised the same way as production's span, and only where there is one.
        match ends_before {
            None => whole,
            Some(_) => without_trailing_blank_lines(&whole),
        }
    };

    // **A guard that extracted nothing would compare nothing and print `ok`.** With
    // `AtTheFirstItemDoc` that is one refactor away: production's `mod.rs` losing every
    // `/// ` line — an item doc rewritten as `/**`, or reordered behind a `#[doc = …]` —
    // empties both sides, and two empty strings are equal.
    if production_body.trim().is_empty() {
        return Some(format!(
            "{file_name}: THE GUARD COULD NOT RUN — no copied content was found in the \
             production original. `CopyBegins` did not find where the copy starts, so the \
             comparison had nothing to compare and would otherwise have passed"
        ));
    }

    // The per-line loop runs **before** the length check, so that a deletion is still
    // reported with the line number where the two sides part rather than as a bare count.
    let ng_lines: Vec<&str> = ng_body.lines().collect();
    let production_lines: Vec<&str> = production_body.lines().collect();
    for (index, (ng_line, production_line)) in
        ng_lines.iter().zip(production_lines.iter()).enumerate()
    {
        if ng_line != production_line {
            return Some(format!(
                "{file_name}, line {} of the copied content — ng's copy and \
                 {production_path} have diverged. This is a VERBATIM COPY and it must not \
                 be edited: not to \
                 satisfy clippy, not to tidy, not to rename. Production is frozen and this \
                 copy is the baseline the port is measured against \
                 (doc/devel/ng/spec/hidden_paralog_filter.md §1.2). Whichever side moved, \
                 revert it. If the line is a path into production, check first that its \
                 `Repoint` is still declared here — a deleted one lands in this message, and \
                 reverting ng's line would put production's path back. If the divergence is \
                 deliberate, release the file from `guarded_copies` in the commit that \
                 introduces it.\n  \
                 ng's copy:  {ng_line:?}\n  the original: {production_line:?}",
                index + 1
            ));
        }
    }
    if ng_lines.len() != production_lines.len() {
        return Some(format!(
            "{file_name}: ng's copy has {} lines of copied content to production's {}, and \
             every line they share is identical — one was added or deleted at the end",
            ng_lines.len(),
            production_lines.len(),
        ));
    }

    // Line for line is not byte for byte: `lines()` drops a trailing `\r` and cannot see a
    // missing final newline, and this module's header claims the stronger property.
    (ng_body != production_body).then(|| {
        format!(
            "{file_name}: a byte-level difference the line comparison cannot see — a line \
             ending, or the final newline. ng's copy is {} bytes of copied content to \
             production's {}",
            ng_body.len(),
            production_body.len(),
        )
    })
}

/// **Every guarded file is still production's, line for line and byte for byte.**
///
/// The failure message names the original and says the file must not be edited, at the
/// same moment a "do not edit" banner would have been read — and a failing build cannot be
/// skimmed past, where a comment can.
#[test]
fn the_copies_are_still_productions() {
    for copy in guarded_copies() {
        if let Some(difference) = how_the_copy_differs(&copy) {
            panic!("{difference}");
        }
    }
}

/// **The comparison's own rejections, pinned on synthetic files.**
///
/// The real pair is identical, so without these every rejecting branch would run only in
/// its accepting direction and a weakening of any of them would be invisible.
#[test]
fn the_comparison_accepts_an_appended_note_and_rejects_every_other_edit() {
    let production = "//! A\n//! B\n\ncode1\ncode2\n";
    let whole_file = |ours: &'static str| {
        how_the_copy_differs(&GuardedCopy {
            file_name: "f.rs",
            production_path: "a synthetic original",
            production_source: production,
            ng_source: ours,
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            repoints: &[],
        })
    };

    assert!(whole_file(production).is_none(), "an identical file");
    assert!(
        whole_file("//! A\n//! B\n//! ng's note\n\ncode1\ncode2\n").is_none(),
        "a note appended to production's header is the one permitted difference",
    );

    assert!(
        whole_file("//! A\n//! REWRITTEN\n\ncode1\ncode2\n").is_some(),
        "a rewritten production header line",
    );
    assert!(
        whole_file("//! ng's note\n//! A\n//! B\n\ncode1\ncode2\n").is_some(),
        "a note prepended rather than appended",
    );
    assert!(
        whole_file("//! A\n\ncode1\ncode2\n").is_some(),
        "a dropped header line",
    );
    assert!(
        whole_file("//! A\n//! B\n\ncode1\ncodeX\n").is_some(),
        "a drifted body line",
    );
    assert!(
        whole_file("//! A\n//! B\n\ncode1\n").is_some(),
        "a deleted body line",
    );
    assert!(
        whole_file("//! A\n//! B\n\ncode1\ncode2\ncode3\n").is_some(),
        "an added body line",
    );
    assert!(
        whole_file("//! A\n//! B\n\ncode1\ncode2").is_some(),
        "a stripped final newline — the line comparison cannot see it, the byte one can",
    );
    assert!(
        whole_file("//! A\n//! B\n\ncode1\r\ncode2\n").is_some(),
        "a CRLF ending on a body line — likewise",
    );
    assert!(
        whole_file("//! A\n//! B\r\n\ncode1\ncode2\n").is_some(),
        "a CRLF ending on a header line: the header prefix is compared by bytes too, or a \
         re-copy that normalised line endings would pass in the one region ng may extend",
    );

    // The other shape: the comparison starts at the first item's `///`, so whatever each
    // side puts above it is its own.
    let items_only = |theirs: &'static str, ours: &'static str| {
        how_the_copy_differs(&GuardedCopy {
            file_name: "f.rs",
            production_path: "a synthetic original",
            production_source: theirs,
            ng_source: ours,
            begins: CopyBegins::AtTheFirstItemDoc,
            ends_before: None,
            repoints: &[],
        })
    };
    assert!(
        items_only(
            "//! their header\npub mod thing;\n\n/// One\npub const A: u8 = 1;\n",
            "//! ng's header\n\n/// One\npub const A: u8 = 1;\n",
        )
        .is_none(),
        "different declarations above the first item are each side's own",
    );
    assert!(
        items_only(
            "//! their header\npub mod thing;\n\n/// One\npub const A: u8 = 1;\n",
            "//! ng's header\n\n/// One\npub const A: u8 = 2;\n",
        )
        .is_some(),
        "a drifted item",
    );
}

/// **A sanctioned substitution is applied where declared and nowhere else.**
///
/// Without this, the whole `Repoint` path runs only in its accepting direction, on the one
/// real file pair that has repoints — so deleting the "exactly once" check, or loosening the
/// exact line match to a prefix, would leave the suite green. Two of its rejections are
/// caught incidentally by the content comparison; the third, a sanctioned line that comes to
/// occur twice, is caught by nothing else.
#[test]
fn a_repoint_is_applied_where_declared_and_nowhere_else() {
    let production = "//! A\n//! [`T`]: crate::paralog::T\n\nuse crate::paralog::T;\ncode\n";
    let with = |ours: &'static str, repoints: &'static [Repoint]| {
        how_the_copy_differs(&GuardedCopy {
            file_name: "f.rs",
            production_path: "a synthetic original",
            production_source: production,
            ng_source: ours,
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            repoints,
        })
    };
    const BOTH: &[Repoint] = &[
        Repoint {
            production_line: "use crate::paralog::T;",
            ng_line: "use crate::ng::paralog::T;",
            why: WhyRepointed::PathIntoProduction,
        },
        Repoint {
            production_line: "//! [`T`]: crate::paralog::T",
            ng_line: "//! [`T`]: crate::ng::paralog::T",
            why: WhyRepointed::PathIntoProduction,
        },
    ];

    assert!(
        with(
            "//! A\n//! [`T`]: crate::ng::paralog::T\n\nuse crate::ng::paralog::T;\ncode\n",
            BOTH,
        )
        .is_none(),
        "both declared substitutions — one in the module header, one in the body — are the \
         copy's to make",
    );
    assert!(
        with(
            "//! A\n//! [`T`]: crate::paralog::T\n\nuse crate::paralog::T;\ncode\n",
            BOTH,
        )
        .is_some(),
        "a declared substitution the copy did not make",
    );
    assert!(
        with(
            "//! A\n//! [`T`]: crate::ng::paralog::T\n\nuse crate::ng::paralog::T;\ncodeX\n",
            BOTH,
        )
        .is_some(),
        "a body line edited beside the sanctioned ones",
    );
    assert!(
        with(
            "//! A\n//! [`T`]: crate::paralog::T\n\nuse crate::ng::paralog::T;\ncode\n",
            &[],
        )
        .is_some(),
        "an undeclared path edit — the whole point of declaring them",
    );

    // The structural rejection: the declared line is not production's any more, so nothing
    // is checking ng's copy and the message has to say that rather than blame the copy.
    const STALE: &[Repoint] = &[Repoint {
        production_line: "use crate::paralog::Renamed;",
        ng_line: "use crate::ng::paralog::Renamed;",
        why: WhyRepointed::PathIntoProduction,
    }];
    let the_guard_is_broken = with("//! A\n//! [`T`]: crate::paralog::T\n\ncode\n", STALE)
        .expect("a repoint that matches nothing must be reported");
    assert!(
        the_guard_is_broken.contains("THE GUARD COULD NOT RUN"),
        "a repoint that no longer matches says the guard is off, not that the copy drifted; \
         got {the_guard_is_broken:?}",
    );
}

/// **A copy whose content the guard cannot find fails, rather than passing on nothing.**
///
/// `AtTheFirstItemDoc` looks for the first `/// ` line; if production's file has none — an
/// item doc rewritten as `/**`, or reordered behind a `#[doc = …]` — both sides extract the
/// empty string and every comparison below holds. That is the guard failing *open*, which is
/// worse than either of its other outcomes because nothing distinguishes it from a pass.
#[test]
fn a_guard_that_can_find_no_copied_content_fails_rather_than_passing() {
    let difference = how_the_copy_differs(&GuardedCopy {
        file_name: "f.rs",
        production_path: "a synthetic original",
        production_source: "//! header\npub mod thing;\n",
        ng_source: "//! ng's header\n",
        begins: CopyBegins::AtTheFirstItemDoc,
        ends_before: None,
        repoints: &[],
    })
    .expect("no copied content must be reported, not silently accepted");
    assert!(
        difference.contains("THE GUARD COULD NOT RUN"),
        "got {difference:?}"
    );
}

/// **A span copy stops where it is declared to, and fails if it cannot find that line.**
///
/// One copy in this module takes a span of a larger production file rather than the whole of
/// it. Two things can go wrong that a whole-file copy cannot: the span could run past its
/// end and compare content ng never copied, and the marker could be renamed in production,
/// after which the extraction would silently run to the end of the file and the comparison
/// would fail somewhere unrelated. The second is a guard failure and says so.
#[test]
fn a_span_copy_stops_at_its_marker_and_fails_when_the_marker_is_gone() {
    let production =
        "//! header\nuse a;\n\n/// One\npub const A: u8 = 1;\n\n/// Two\npub const B: u8 = 2;\n";
    let span = |ours: &'static str, ends_before: Option<&'static str>| {
        how_the_copy_differs(&GuardedCopy {
            file_name: "f.rs",
            production_path: "a synthetic original",
            production_source: production,
            ng_source: ours,
            begins: CopyBegins::AtTheFirstItemDoc,
            ends_before,
            repoints: &[],
        })
    };

    assert!(
        span(
            "//! ng's header\nuse a;\n\n/// One\npub const A: u8 = 1;\n",
            Some("/// Two"),
        )
        .is_none(),
        "a copy of the first item alone, stopping before the second — the blank line between \
         them belongs to neither side",
    );
    assert!(
        span(
            "//! ng's header\n\n/// One\npub const A: u8 = 1;\n\n/// Two\npub const B: u8 = 2;\n",
            Some("/// Two"),
        )
        .is_some(),
        "a copy that ran past the declared end",
    );

    let the_guard_is_broken = span(
        "//! ng's header\n\n/// One\npub const A: u8 = 1;\n",
        Some("/// A line production does not have"),
    )
    .expect("an end marker that matches nothing must be reported");
    assert!(
        the_guard_is_broken.contains("THE GUARD COULD NOT RUN"),
        "a renamed end marker says the guard is off, not that the copy drifted; got \
         {the_guard_is_broken:?}",
    );
}

/// **A declaration that is not one of the three sanctioned kinds is refused, for every kind.**
///
/// The sanction check is the one thing the content comparison cannot cross-check: the
/// comparison rewrites production's file *using* the declaration, so a declaration that
/// smuggles an edit agrees with a copy that carries the same edit. Until this test existed,
/// the check had no case in its rejecting direction at all, and two of its three kinds were
/// containment tests that let an arbitrary edit through. Measured, before the fix: declaring
/// `locus_score.rs`'s import repoint with an ng-side line reading
/// `use crate::ng::paralog::{GridSpec, SfsPriorSpec}; const SMUGGLED: f64 = 0.05;` and putting
/// that line in ng's copy left the guard green.
#[test]
fn a_declaration_that_smuggles_an_edit_is_refused() {
    let production = "//! header\n\nuse crate::paralog::T;\n/// Tuning knobs for [`gone`].\npub(crate) struct S {\n";
    let refuse = |repoints: &'static [Repoint]| {
        how_the_copy_differs(&GuardedCopy {
            file_name: "f.rs",
            production_path: "a synthetic original",
            production_source: production,
            // ng's copy carries exactly what the declaration says, so only the declaration
            // itself can be what refuses this.
            ng_source: {
                let mut ours = String::from("//! header\n//! ng's note\n\n");
                for line in production.lines().skip(2) {
                    let replaced = repoints
                        .iter()
                        .find(|r| r.production_line == line)
                        .map_or(line, |r| r.ng_line);
                    ours.push_str(replaced);
                    ours.push('\n');
                }
                Box::leak(ours.into_boxed_str())
            },
            begins: CopyBegins::AfterTheModuleHeader,
            ends_before: None,
            repoints,
        })
    };

    // Each kind, well formed: accepted.
    assert!(
        refuse(&[
            Repoint {
                production_line: "use crate::paralog::T;",
                ng_line: "use crate::ng::paralog::T;",
                why: WhyRepointed::PathIntoProduction,
            },
            Repoint {
                production_line: "/// Tuning knobs for [`gone`].",
                ng_line: "/// Tuning knobs for the thing.",
                why: WhyRepointed::DocLinkToAnItemNgDoesNotCopy {
                    link: "[`gone`]",
                    replaced_by: "the thing",
                },
            },
            Repoint {
                production_line: "pub(crate) struct S {",
                ng_line: "pub struct S {",
                why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
            },
        ])
        .is_none(),
        "three well-formed declarations, one of each kind",
    );

    // Each kind, smuggling something else on the same line: refused, and named as a guard
    // failure rather than as a drifted copy.
    for (what, repoints) in [
        (
            "a path substitution that also declares a constant",
            &[Repoint {
                production_line: "use crate::paralog::T;",
                ng_line: "use crate::ng::paralog::T; const SMUGGLED: f64 = 0.05;",
                why: WhyRepointed::PathIntoProduction,
            }][..],
        ),
        (
            "a doc-link substitution that also rewrites the sentence",
            &[Repoint {
                production_line: "/// Tuning knobs for [`gone`].",
                ng_line: "/// The fallback prior is 0.05, not 0.03.",
                why: WhyRepointed::DocLinkToAnItemNgDoesNotCopy {
                    link: "[`gone`]",
                    replaced_by: "the thing",
                },
            }][..],
        ),
        (
            "a visibility widening that also renames the item",
            &[Repoint {
                production_line: "pub(crate) struct S {",
                ng_line: "pub struct Renamed {",
                why: WhyRepointed::VisibilityNgsPublicModuleNeeds,
            }][..],
        ),
    ] {
        let refused = refuse(repoints).unwrap_or_else(|| {
            panic!(
                "{what} must be refused; the declaration is the one thing the content \
                    comparison cannot cross-check"
            )
        });
        assert!(
            refused.contains("THE GUARD COULD NOT RUN"),
            "{what}: the message must say the guard is off, not that ng's copy drifted; got \
             {refused:?}",
        );
    }
}

/// **The item `calibration.rs`'s span stops before is still the one it means to stop before.**
///
/// The span ends at `cohort_inbreeding`'s doc comment, and that helper is the one piece of
/// the filter ng replaces rather than ports (see this module's header). Nothing else ties
/// the marker to the item: rename `cohort_inbreeding` while keeping its first doc sentence
/// and the guard stays green over a copy that now means something different. This asserts
/// the declaration production still makes.
#[test]
fn the_item_the_span_stops_before_is_still_productions() {
    let calibrate = include_str!("../../var_calling/paralog_filter/calibrate.rs");
    assert!(
        calibrate.contains("fn cohort_inbreeding("),
        "src/var_calling/paralog_filter/calibrate.rs no longer declares `cohort_inbreeding`. \
         calibration.rs's span is declared to stop before that helper's doc comment because \
         it is the piece ng replaces with a per-sample coefficient; if it has been renamed or \
         removed, the span's `ends_before` marker and this module's header both have to be \
         re-decided rather than adjusted"
    );
}

/// **Every file in this module is guarded, ng's own, or explicitly released.**
///
/// `guarded_copies` is hand-maintained and the release protocol *deletes* from it, so a
/// dropped entry and a deliberate release are indistinguishable — and a guard that covers
/// nothing still prints `ok`. This makes the directory listing the authority instead: a
/// file in none of the three lists fails the build.
#[test]
fn every_file_in_the_module_is_guarded_ng_s_own_or_released() {
    // Read off `guarded_copies` rather than restated beside it: a second hand-maintained
    // list is the very thing this test exists to remove, and a length check — which is all
    // two lists can cheaply agree on — passes when one name is swapped in only one of them.
    let guarded: Vec<&str> = guarded_copies().iter().map(|c| c.file_name).collect();
    const NG_OWN: [&str; 3] = ["mod.rs", "copy_fidelity.rs", "production_parity.rs"];

    // Grows as files are released; keep it in step with the table in this module's header.
    const RELEASED: [&str; 0] = [];

    let module_directory = concat!(env!("CARGO_MANIFEST_DIR"), "/src/ng/paralog");
    let mut guarded_files_on_disk = 0;
    for entry in std::fs::read_dir(module_directory).expect("the module directory is readable") {
        let file_name = entry.expect("a directory entry").file_name();
        let file_name = file_name.to_string_lossy().into_owned();
        if !file_name.ends_with(".rs") {
            continue;
        }
        assert!(
            guarded.contains(&file_name.as_str())
                || NG_OWN.contains(&file_name.as_str())
                || RELEASED.contains(&file_name.as_str()),
            "{file_name} is in src/ng/paralog/ but is neither guarded, ng's own, nor \
             released — add it to `guarded_copies` and to `GUARDED`, or say which of the \
             other two it is",
        );
        guarded_files_on_disk += usize::from(guarded.contains(&file_name.as_str()));
    }
    assert_eq!(
        guarded_files_on_disk,
        guarded.len(),
        "every file `guarded_copies` names must be on disk in src/ng/paralog/: {guarded:?}",
    );
}
