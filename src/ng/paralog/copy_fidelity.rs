//! **These files are still production's, verbatim — asserted, not claimed.**
//!
//! Two files are guarded today: `coverage_model.rs`, 1,158 lines past its module header,
//! against `src/paralog/coverage_model.rs`; and `model_params.rs`, 236 lines, against the
//! items of `src/paralog/mod.rs`. ng owns this file outright. It exists so the copies can
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
        },
        GuardedCopy {
            file_name: "model_params.rs",
            production_source: include_str!("../../paralog/mod.rs"),
            ng_source: include_str!("model_params.rs"),
            begins: CopyBegins::AtTheFirstItemDoc,
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
        production_source,
        ng_source,
        begins,
    } = *copy;

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
                 revert it — or, if the divergence is deliberate, release the file from \
                 `guarded_copies` in the commit that introduces it.\n  \
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

/// **Every file in this module is guarded, ng's own, or explicitly released.**
///
/// `guarded_copies` is hand-maintained and the release protocol *deletes* from it, so a
/// dropped entry and a deliberate release are indistinguishable — and a guard that covers
/// nothing still prints `ok`. This makes the directory listing the authority instead: a
/// file in none of the three lists fails the build.
#[test]
fn every_file_in_the_module_is_guarded_ng_s_own_or_released() {
    const GUARDED: [&str; 2] = ["coverage_model.rs", "model_params.rs"];
    const NG_OWN: [&str; 2] = ["mod.rs", "copy_fidelity.rs"];
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
            GUARDED.contains(&file_name.as_str())
                || NG_OWN.contains(&file_name.as_str())
                || RELEASED.contains(&file_name.as_str()),
            "{file_name} is in src/ng/paralog/ but is neither guarded, ng's own, nor \
             released — add it to `guarded_copies` and to `GUARDED`, or say which of the \
             other two it is",
        );
        guarded_files_on_disk += usize::from(GUARDED.contains(&file_name.as_str()));
    }
    assert_eq!(
        guarded_files_on_disk,
        GUARDED.len(),
        "every file named in `GUARDED` must be on disk in src/ng/paralog/",
    );
    assert_eq!(
        guarded_copies().len(),
        GUARDED.len(),
        "`GUARDED` and `guarded_copies` must name the same set of files",
    );
}
