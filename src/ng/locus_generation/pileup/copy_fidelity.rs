//! **The copy is production's, and this is what says so.**
//!
//! ng owns this file outright — it is not one of the eight copies, and it exists so
//! the copies can be checked from outside themselves (a test living inside
//! `tests.rs` cannot compare `tests.rs` against its original, because it is part of
//! it). Keeping it here also keeps the boundary a **file** boundary: everything in
//! `pileup/` is either a verbatim copy or ng's own, and no file is both.
//!
//! # Why the copies carry no "do not edit" banner (owner, 2026-07-29)
//!
//! A review asked for a per-file banner, on a real worry: someone arrives at
//! `genome_walk.rs` from a grep hit or a stack trace, tidies something, and the baseline
//! this whole plan rests on is quietly gone. **Declined, deliberately.** A banner would
//! widen the plan's sanctioned edit set, and it would cost `genome_walk.rs` and `errors.rs`
//! their byte-identity with their originals — `diff` exit 0, the cheapest and most durable
//! fidelity statement available, and one that survives a future re-copy.
//!
//! This test is the answer instead: it fails the build, at the same moment the banner would
//! have been read, with a message that names the production original and says the file must
//! not be edited. A failing build is a stronger signal than a comment, and it cannot be
//! skimmed past.

/// **"These eight files are production's, verbatim" — asserted, not claimed.**
///
/// That property is what the whole port rests on: it is the baseline every later
/// change is measured against, and `locus_generation_pileup.md` §3's rule is
/// *transcribe first, change second*. Until this test existed the only evidence was a
/// `diff` run by hand into a gitignored scratch directory, so it protected the property
/// exactly once — and nothing would have failed if either tree drifted. Both suites are
/// independent copies, so both stay green through a divergence; that is precisely why
/// the check has to be textual.
///
/// The sanctioned divergences are stripped from **ng's side only**, and each is
/// stripped by an exact match rather than a prefix or a pattern:
///
/// - `read_group: ReadGroupId(0),` — the field ng's `PreparedRead` adds, in the test
///   fixtures. Not making that type `#[non_exhaustive]` is what forces these.
/// - `use crate::ng::types::ReadGroupId;` — the import those lines require.
/// - the module doc note at the head of this file, which is the only prose ng adds.
///
/// The two module-path repoints are **lifted out of both sides and matched by presence
/// rather than by position**, because their position is not ours to control: repointing
/// `crate::pileup::walker::decompose` to `crate::ng::…::decompose` changes where rustfmt
/// sorts the line, so `cigar_cursor.rs`'s test imports come out in a different order
/// from production's and no line-for-line comparison can survive it. Each repoint must
/// still appear exactly once on each side, so a *new* path edit — or a dropped one —
/// fails here and has to be sanctioned deliberately. Everything else matches line for
/// line, in order.
///
/// *(That reordering is how this test earned its keep on its first run: it is a real
/// divergence, introduced by A3's `cargo fmt` and not recorded anywhere until now.)*
///
/// `driver.rs` → `genome_walk.rs` is the pairing that cannot be inferred from the
/// filenames, which is the other reason this belongs in the tree rather than in a
/// scratch file.
#[test]
fn the_eight_copies_are_still_productions() {
    /// One ng line that is allowed to have no counterpart in production's file.
    fn is_ng_addition(line: &str) -> bool {
        let line = line.trim();
        line == "read_group: ReadGroupId(0)," || line == "use crate::ng::types::ReadGroupId;"
    }

    /// The path repoints, as `(production line, ng line)`.
    const REPOINTS: &[(&str, &str)] = &[
        (
            "use crate::pileup::walker::tests::MockFasta;",
            "use crate::ng::locus_generation::pileup::tests::MockFasta;",
        ),
        (
            "use crate::pileup::walker::decompose::decompose;",
            "use crate::ng::locus_generation::pileup::decompose::decompose;",
        ),
    ];

    // (label, production source, ng copy). `include_str!` resolves relative to this
    // file, so both sides are pinned at compile time — a deleted file is a build error,
    // not a silently skipped case.
    let pairs: [(&str, &str, &str); 8] = [
        (
            "driver.rs -> genome_walk.rs",
            include_str!("../../../pileup/walker/driver.rs"),
            include_str!("genome_walk.rs"),
        ),
        (
            "open_record.rs",
            include_str!("../../../pileup/walker/open_record.rs"),
            include_str!("open_record.rs"),
        ),
        (
            "cigar_cursor.rs",
            include_str!("../../../pileup/walker/cigar_cursor.rs"),
            include_str!("cigar_cursor.rs"),
        ),
        (
            "decompose.rs",
            include_str!("../../../pileup/walker/decompose.rs"),
            include_str!("decompose.rs"),
        ),
        (
            "active_read_set.rs",
            include_str!("../../../pileup/walker/active_read_set.rs"),
            include_str!("active_read_set.rs"),
        ),
        (
            "chain_id_allocator.rs",
            include_str!("../../../pileup/walker/chain_id_allocator.rs"),
            include_str!("chain_id_allocator.rs"),
        ),
        (
            "errors.rs",
            include_str!("../../../pileup/walker/errors.rs"),
            include_str!("errors.rs"),
        ),
        (
            "tests.rs",
            include_str!("../../../pileup/walker/tests.rs"),
            include_str!("tests.rs"),
        ),
    ];

    for (label, production, ours) in pairs {
        // The module doc note. ng may only **append** to production's header, never
        // rewrite it: the leading `//!` run on ng's side must start with production's,
        // line for line, and only the lines past it are ng's own. Skipping a prefix by
        // count instead would let a rewritten header pass while a suffix was dropped.
        let head = |source: &str| source.lines().take_while(|l| l.starts_with("//!")).count();
        let (their_head, our_head) = (head(production), head(ours));
        assert!(
            our_head >= their_head
                && ours
                    .lines()
                    .take(their_head)
                    .eq(production.lines().take(their_head)),
            "{label}: production's module header must survive verbatim; ng's note is \
             appended to it, not written over it",
        );

        // Each repoint must be present exactly once on its own side, and is then lifted
        // out so the remaining lines can be compared in order. Checked before the
        // comparison, so a *dropped* repoint fails as a missing line rather than as a
        // confusing offset further down.
        let mut theirs: Vec<&str> = Vec::new();
        for line in production.lines().skip(their_head) {
            match REPOINTS.iter().find(|(theirs, _)| *theirs == line.trim()) {
                Some((_, ours_line)) => assert_eq!(
                    ours.lines().filter(|l| l.trim() == *ours_line).count(),
                    1,
                    "{label}: production's `{}` should be repointed to `{ours_line}` \
                     exactly once in ng's copy",
                    line.trim(),
                ),
                None => theirs.push(line),
            }
        }
        let mine: Vec<&str> = ours
            .lines()
            .skip(our_head)
            .filter(|line| !is_ng_addition(line))
            .filter(|line| !REPOINTS.iter().any(|(_, ours)| *ours == line.trim()))
            .collect();

        assert_eq!(
            mine.len(),
            theirs.len(),
            "{label}: ng's copy has {} lines to production's {} after stripping the \
             sanctioned additions — the copy has drifted",
            mine.len(),
            theirs.len(),
        );
        for (number, (ours_line, their_line)) in mine.iter().zip(theirs.iter()).enumerate() {
            assert_eq!(
                ours_line,
                their_line,
                "{label}, line {} of the stripped file: this is a VERBATIM COPY of \
                 src/pileup/walker/, and it must not be edited — not to satisfy clippy, \
                 not to tidy, not to rename. Production is frozen and this copy is the \
                 baseline the whole port is measured against \
                 (doc/devel/ng/spec/locus_generation_pileup.md §3). If the divergence is \
                 deliberate, it belongs in this test's sanctioned set, in the commit that \
                 introduces it",
                number + 1,
            );
        }
    }
}
