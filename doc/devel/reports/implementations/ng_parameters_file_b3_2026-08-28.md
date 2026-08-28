# ng parameters file — B3: the provenance comments

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone B, step B3 — the last step of Milestone B
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §4, §8
**Code:** [src/ng/calling/parameters_file/to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs), and the file it now produces at [testdata/every_shape_as_written.toml](../../../../src/ng/calling/parameters_file/testdata/every_shape_as_written.toml)

---

## 1. Plan

Each defaulted value gets a comment beside it saying where the default came from. **This is why
the format is TOML** — spec §4 chose it over JSON on the grounds that "what a person needs beside a
number is where it came from and what moving it costs", and that in JSON an annotation has to
become another field, which then has to be parsed, validated and kept in step with the value it
describes.

## 2. What the file now carries

**162 lines where B2's was 81, and 81 of them are explanation.** Only **two** of those lines are
per-row — one note, on the one value the fixture defaults inside a table. The comments cost a fixed
number of lines rather than one a row, so a run of 3,000 samples whose coefficients were all fitted
writes 3,000 rows and no notes.

Three kinds:

- **A note above each section**, saying what the section is in the reader's language and what its
  absences mean. `[repeat_tracts]` carries the most, **30 lines** — it is the largest section, its
  vocabulary is the fit's rather than a geneticist's, and both of the previous reader's unanswered
  guesses lived there. `[contamination]` is next at 10, because spec §5's three states meet there
  and B2's review found them to be the three a reader cannot tell apart.
- **A note above each defaulted value**, saying where that default came from — the **five** the
  file can express are in one `origins` module, so a number's origin and the sentence a reader is
  shown cannot drift apart, and so that step E1 has one list to reconcile against rather than five
  call sites.
- **A note at the top of the file** giving the one rule for editing it.

**Notes go above the thing they describe, not after it.** A trailing comment lengthens a line that
is already the longest thing in the file, and TOML permits a comment inside a multi-line array — so
a per-row note sits on the lines above its row, where a reader meets it first. Generated notes are
wrapped to a stated width of **80 characters all in** — the `# ` and a row's indent are part of it
— counted in characters and not bytes. **No comment line in the produced file exceeds 80**, against
14 data lines that do.

## 3. Two widenings of B3's contract, both from the review of B2, both recorded rather than assumed

B3's contract reads "each defaulted value gets a comment beside it saying where the default came
from". Two comments the file now carries are not that, and both close a hazard the review measured:

- **The offset convention on the length spectra.** `shares_by_repeat_offset = [0.1, 0.8, 0.1]` runs
  `-span ..= +span` from the reference tract length, and that convention lived only in the Rust.
  The review called it "the only edit in the file that is wrong without being invalid": a user who
  reads the array as starting at zero — the natural reading of an array called *by offset* —
  produces a file that parses, deserialises, and shifts every length prior in that stratum by one
  repeat unit.
- **The rule for editing a value at all.** Spec §1.2 goal 3's own worked example is raising one
  library's error rate; the review ran it and found that after the edit the row still says
  `warrant = "fitted_here", observations = { reads = 812344 }`, so the file asserts a hand-typed
  number was fitted from 812,344 reads. **Goal 2 is what a goal-3 edit silently breaks.** The
  file's opening note now says: change the warrant to `supplied` and delete the `observations`.

Both are the same machinery and one comment each. Neither is a shape change.

## 4. ⛦ What B3 cannot deliver, and it is spec §8's headline

**§8's third bullet asks for a comment this file has nowhere to put.** The per-(stratum × slippage
group) slippage numbers are to be defaulted from the GIAB HG002 alignments, and §8 says the default
"is marked `Defaulted` like the rest, and its origin — which alignments, at what depth, on which
date — is written into the file as a comment beside it, so a user can see what they are inheriting
without opening this document".

**A slippage number carries no `Warrant`.** It carries a *smoothing origin* — its stratum's own
fit, its period's curve, or a blend — which is a different vocabulary with no word for *defaulted*.
So there is no state in which a slippage number is defaulted, and nothing for §8's comment to
attach to. This is the same gap the design review of B1 raised as its Major and that
`PROJECT_STATUS.md` has carried open since Checkpoint A: **spec §2.1's wholesale demotion has
nowhere to write itself for these numbers either.** One decision fixes both.

The second half of §8's bullet is not blocked by the shape at all: **the measurement does not
exist** (§12 question 1), so no run can produce a defaulted slippage number today whatever the
shape says. The `origins` module records both halves where the next implementer will meet them.

## 5. Changes made

| file | change |
|---|---|
| `parameters_file/to_toml.rs` | `note`, `note_lines`, `scalar_with_note`, `one_a_line_with_notes` and `wrapped`; the `origins` module with the five defaults the file can express; a note on every section; five tests |
| `parameters_file/testdata/every_shape_as_written.toml` | regenerated — 81 → 161 lines |

## 6. Tests

**Five new; the module's suite went 58 → 63 passing**, and three of the five that were already there were strengthened by the review.

- `every_defaulted_number_says_where_its_default_came_from` — and that the fitted rows around it
  carry no note, so the comment is about the default rather than decoration.
- `a_run_that_defaulted_nothing_writes_no_per_row_notes` — the cost does not scale with the cohort.
- `the_comments_change_what_a_reader_learns_and_not_what_it_reads` — **the one that matters most.**
  A TOML comment runs to the end of its line, so a note that landed inside a row rather than above
  it would silently truncate the document. Stripping every comment line must leave a file that
  reads back as the same value, and so must leaving them in.
- `no_comment_line_is_longer_than_the_prose_it_carries` — prose that needs scrolling is prose
  nobody reads. It caught four of my own section notes at 83 and 84 characters.
- `a_note_wraps_by_characters_at_the_width_it_states` — added after the mutation pass, which found
  that moving the wrap boundary by one and counting bytes instead of characters both survived
  everything: no note in the file happens to sit on the boundary, and none carries a byte that is
  not a character.

## 7. Validation

Run in the dev container, `./scripts/dev.sh`:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean, exit 0 |
| `cargo test --lib ng::calling::parameters_file` | **63 passed, 0 failed, 2 ignored** |
| `cargo test --lib` | **4,986 passed, 0 failed, 13 ignored** |

`cargo test --all-targets --all-features` is not the gate: pre-existing panic in
`benches/psp_writer_perf.rs:386`, verified on clean `main`.

## 8. What the reviews found, and it was all in the prose

Two agents: a mutation pass over the comment machinery, and **the file's reader again** — the same
geneticist exercise B2's review ran on the uncommented file, now re-run to ask whether the comments
answer the questions it had.

**Four comments were wrong, and one answered the previous reader's exact question with the wrong
option.** The contamination note said "how much of each library's **DNA** came from somebody else"
where the quantity is a share of that read group's **reads** and the grain is the lane — and the
spec's own reason is that index hopping happens on a flowcell rather than in a tube, which mislabels
reads and does not put another individual's DNA in the tube. A reader following that note thinks in
libraries, edits one row, and believes a four-lane library is uncorrected when three of its lanes
still are. Also wrong: "each value carries a `warrant`" (eight of the file's numbers do, and the
file's one editing rule rested on it); a batching note whose justification was refuted by the three
rows printed beneath it; and one word, `rung`, given an authoritative four-way definition for the
prior's ladder while a second `rung` meaning something else sits twenty lines below it.

**The mutation pass ran 34 and found five survivors**, four of them blind spots rather than defects:
two origin strings that no fixture ever defaults, so a text put beside the wrong quantity was
invisible; and the wrapper's boundary and its unit, since no note happens to sit at 78 characters or
to carry a byte that is not a character. **The fifth was a real hole** — the flat concentration's
note fired whatever the warrant said, and the test for "nothing defaulted, so nothing explains a
default" had only listed two of the five texts. All five now fail a test.

**And the mutation pass established that the step's stated hazard is genuinely guarded.** A newline
inside a note — a comment landing inside a row rather than above it — fails seven tests, because
`WarrantedValue` refuses unknown fields and requires three, so a truncated inline table fails to
parse rather than deserialising short.

## 9. One defect this step made and caught

**My own edit script ate the Rust line continuations** in the `origins` strings, so the first
produced file carried runs of ten spaces mid-sentence — "used exactly as they came, because
          no usable error rate". It reached the golden file and was caught by reading that file
rather than by any test. The strings are now `concat!` of whole lines, which cannot do it again.
Recorded because the mechanism generalises: a generated artefact is only as good as the last
reading of it, and no assertion in this module would have failed on the mangled text.
