//! **These files are still production's, verbatim — asserted, not claimed.**
//!
//! Three files are guarded today: `coverage_model.rs`, 1,158 lines past its module header,
//! against `src/paralog/coverage_model.rs`; `locus_score.rs`, 797 lines, against
//! `src/paralog/locus_score.rs`; and `model_params.rs`, 236 lines, against the items of
//! `src/paralog/mod.rs`. `locus_score.rs` differs from its original on three declared lines
//! and no others; the two whole-file copies differ on none. ng owns this file outright. It exists so the copies can
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
//! # The lines a copy may change: paths into `src/paralog/`
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
//! **The rule is exact rather than loose.** Production's side must appear **exactly once**,
//! ng's side must be **exactly** the declared replacement at that same line, and the
//! indentation and the line ending must match. A repoint may only turn a `crate::paralog`
//! path into a `crate::ng::paralog` one — asserted, so that the mechanism cannot be used to
//! sanction a changed constant. So a *new* path edit, or a dropped declaration, still fails;
//! only the sanctioned substitutions pass, and each is written down beside the file it
//! applies to.
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

/// One line a copy is allowed to change, because leaving it would point ng's code back
/// into production: production's text, and the text ng puts in its place.
struct Repoint {
    production_line: &'static str,
    ng_line: &'static str,
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
    /// The whole of the production original, read at compile time.
    production_source: &'static str,
    /// The whole of ng's copy, read at compile time.
    ng_source: &'static str,
    /// Where the copied content starts on each side.
    begins: CopyBegins,
    /// The sanctioned path substitutions, each of which must occur exactly once.
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
            production_source: include_str!("../../paralog/coverage_model.rs"),
            ng_source: include_str!("coverage_model.rs"),
            begins: CopyBegins::AfterTheModuleHeader,
            repoints: &[],
        },
        GuardedCopy {
            file_name: "locus_score.rs",
            production_source: include_str!("../../paralog/locus_score.rs"),
            ng_source: include_str!("locus_score.rs"),
            begins: CopyBegins::AfterTheModuleHeader,
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
                },
                Repoint {
                    production_line: "//! [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel",
                    ng_line: "//! [`SingleCopyCoverageModel`]: crate::ng::paralog::SingleCopyCoverageModel",
                },
                Repoint {
                    production_line: "/// [`SingleCopyCoverageModel`]: crate::paralog::SingleCopyCoverageModel",
                    ng_line: "/// [`SingleCopyCoverageModel`]: crate::ng::paralog::SingleCopyCoverageModel",
                },
            ],
        },
        GuardedCopy {
            file_name: "model_params.rs",
            production_source: include_str!("../../paralog/mod.rs"),
            ng_source: include_str!("model_params.rs"),
            begins: CopyBegins::AtTheFirstItemDoc,
            repoints: &[],
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

/// The copied content itself, with its line endings and its final newline intact.
fn the_copied_content(source: &str, begins: CopyBegins) -> String {
    source
        .split_inclusive('\n')
        .skip(lines_before_the_copy(source, begins))
        .collect()
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
        assert!(
            repoint.production_line.contains("crate::paralog")
                && repoint.ng_line.contains("crate::ng::paralog"),
            "{}: a `Repoint` may only turn a `crate::paralog` path into a \
             `crate::ng::paralog` one, and `{}` → `{}` does not — the mechanism sanctions \
             paths, not changes of behaviour",
            copy.file_name,
            repoint.production_line,
            repoint.ng_line,
        );
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
        production_source: _,
        ng_source,
        begins,
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
                     ng's copy:    {ng_line:?}\n  src/paralog/: {production_line:?}",
                    index + 1
                ));
            }
        }
    }

    let production_body = the_copied_content(production_source, begins);
    let ng_body = the_copied_content(ng_source, begins);

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
                "{file_name}, line {} of the copied content — ng's copy and src/paralog/ \
                 have diverged. This is a VERBATIM COPY and it must not be edited: not to \
                 satisfy clippy, not to tidy, not to rename. Production is frozen and this \
                 copy is the baseline the port is measured against \
                 (doc/devel/ng/spec/hidden_paralog_filter.md §1.2). Whichever side moved, \
                 revert it. If the line is a path into src/paralog/, check first that its \
                 `Repoint` is still declared here — a deleted one lands in this message, and \
                 reverting ng's line would put production's path back. If the divergence is \
                 deliberate, release the file from `guarded_copies` in the commit that \
                 introduces it.\n  \
                 ng's copy:    {ng_line:?}\n  src/paralog/: {production_line:?}",
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
            production_source: production,
            ng_source: ours,
            begins: CopyBegins::AfterTheModuleHeader,
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
            production_source: theirs,
            ng_source: ours,
            begins: CopyBegins::AtTheFirstItemDoc,
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
            production_source: production,
            ng_source: ours,
            begins: CopyBegins::AfterTheModuleHeader,
            repoints,
        })
    };
    const BOTH: &[Repoint] = &[
        Repoint {
            production_line: "use crate::paralog::T;",
            ng_line: "use crate::ng::paralog::T;",
        },
        Repoint {
            production_line: "//! [`T`]: crate::paralog::T",
            ng_line: "//! [`T`]: crate::ng::paralog::T",
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
        production_source: "//! header\npub mod thing;\n",
        ng_source: "//! ng's header\n",
        begins: CopyBegins::AtTheFirstItemDoc,
        repoints: &[],
    })
    .expect("no copied content must be reported, not silently accepted");
    assert!(
        difference.contains("THE GUARD COULD NOT RUN"),
        "got {difference:?}"
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
