//! **The copy is production's, and this is what says so.**
//!
//! ng owns this file outright — it is not one of the eight copies, and it exists so
//! the copies can be checked from outside themselves (a test living inside
//! `tests.rs` cannot compare `tests.rs` against its original, because it is part of
//! it). Keeping it here also keeps the boundary a **file** boundary: every file in
//! `pileup/` is either ng's own, a *still*-verbatim copy, or a copy this plan has
//! **released** and begun changing.
//!
//! That third state is the one the wording used to deny — "no file is both" was
//! written when there were only two — and A0 created it three lines below, in
//! `genome_walk.rs`, `open_record.rs` and `errors.rs`. A released file is no longer
//! guarded here and no longer claims to be production's; the release table below is
//! the record of which are which.
//!
//! # Narrowed file by file, never retired wholesale
//!
//! Plan 3 changes the copies on purpose, and each step releases **the file it
//! changes** from the set below — in the same commit, and recorded in that file's
//! own module header. Switching the guard off at the first divergence would throw
//! away a working check on every file that is still a copy, which is most of them.
//!
//! | released at | files | why |
//! |---|---|---|
//! | A0 | `driver.rs` → `genome_walk.rs`, `open_record.rs`, `errors.rs` | the reference is reached through ng's `RefSeq`, not production's `MultiChromRefFetcher` |
//! | B2 | `tests.rs` | the walk emits ng's `SampleLocusObservations`, so the inherited suite is adapted to it (spec §12) |
//! | D2 | `active_read_set.rs` | a per-read "ever contributed" flag, so `reads_silent_over_footprint` can be counted — the read spec §13's accounting cannot otherwise place (spec §6) |
//! | perf 2026-08-05 | `cigar_cursor.rs` | the cursor states when its two queries coincide, so the fold can reuse the events it already has (−5.8 % of the walk at ~130×, −4.0 % at 30×) |
//! | depth cap 2026-08-05 | `chain_id_allocator.rs` | the two ceilings that made the walk refuse reads on ordinary data — 4,096 held reads and 10,000 pending mates — plus the counter that says how much headroom the second one has |
//!
//! # The pin is lifted, and this file's job is now narrower
//!
//! The owner released the last constraint on 2026-08-05: *"You might copy any file from
//! production and then change it. ng is destined to replace production in the near
//! future."* So there is no file this test may not lose. What it still does, for as long
//! as one entry remains, is say **which** files are untouched — and `decompose.rs` is the
//! last of the eight. When it goes, this test goes with it rather than becoming a
//! zero-length array that asserts nothing.
//!
//! **The bar for releasing a file is now a measured gain, not a plan step** (owner,
//! 2026-08-05: *"You can copy the files from production into ng and modify them there.
//! No problem with ng diverging more and more with production."*). Two other edits to
//! this same file were measured in the same session and **reverted** rather than
//! released — a smaller returned-event buffer (−0.21 % at 130×, nothing at 30×) and a
//! narrower event type (nothing at either) — because a permanent divergence is not worth
//! a change that cannot be seen.
//!
//! Everything still listed in `PAIRS` is byte-for-byte production's, modulo the
//! sanctioned additions below.
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

/// **"These files are still production's, verbatim" — asserted, not claimed.**
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
#[test]
fn the_remaining_copies_are_still_productions() {
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
    //
    // **A released file is deleted from this array, not commented out**, so a `git log`
    // on it shows the whole set that was ever guarded and the array itself shows only
    // what is guarded *now*. The release table in this module's header is the record.
    let pairs: [(&str, &str, &str); 1] = [(
        "decompose.rs",
        include_str!("../../../pileup/walker/decompose.rs"),
        include_str!("decompose.rs"),
    )];

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
