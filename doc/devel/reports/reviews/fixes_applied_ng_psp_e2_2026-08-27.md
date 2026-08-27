# Fixes applied — ng psp E2 (the review of `4a1d789f`)

*2026-08-27. Answers [the review](ng_psp_e2_2026-08-27.md). Step E2 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. Three of the four Majors are one test, and all three agents found it

**`a_read_that_comes_back_across_a_block_boundary_arrives_again` contained neither behaviour its
docstring named.** Its fixture was one read pair whose mates fell either side of the cut, so the
live set was **empty at the boundary** — the test asserted so itself — and it read the second block
with a *fresh* reader, which has no state to carry. So neither the restatement of a read still
covering, nor a reader that failed to reset, could show up there.

Measured: of 22 mutations one agent ran, this test caught exactly **one**, and fourteen other tests
caught that one too. Making the writer's `start_block` a no-op, making the reader's a no-op, and
injecting the step's own headline defect all left it green.

**Fixed** with a fixture that has a read still covering at the cut, a pass with the **same** reader
carried across the boundary, an assertion on the new block's first arrival list so the restatement
is visible rather than inferred, and the fresh-reader pass kept beside it and required to agree.
Measured on the strengthened test: it now kills the writer's forgotten reset, the reader's, and a
reader that keeps a file-scope memory of every id it has named — **three where it killed none.**

## 2. The other Majors

| what | what was done |
|---|---|
| **The sortedness property asserted an invariant in a regime where it cannot break.** The live set can only stop being sorted when an arrival sorts *below* an id already live, and uniform random bytes never build a set for one to sort under: the interleaving arm was entered **0 times over 600 cases**, and a mutation making the merge a plain append left the test green | the bytes are a **real record's, damaged**, and the reader is seeded with a live set first. Measured: **113 of 600 damaged streams are now accepted** and have their set checked, where uniform bytes were refused at the first count almost every time. ⚠ It still does not reach the interleaving arm — zero entries either way — so the docstring now says which tests do hold that arm rather than claiming this one does |
| **"A shrink that removed the re-entry would fail rather than pass quietly"** was not what the property test's assertion said: both sides of the arrival count come from the same records, so a draw with no re-entry satisfies it exactly | a `prop_assume!` throws those draws away, so the docstring is true. And the generated records are now **cut into two blocks** at a point the generator chooses, because two of the three re-entry tests used one block and E3 stresses that seam — the count then carries one restatement per read still covering at the cut |
| **`LiveSetChanges::departed()` had no test that could fail** — replacing its body with `&[]` left all 228 tests green, because the only test reading it asserted the list was *empty* | `changes_departed_names_the_reads_that_stopped_covering` asserts a non-empty list. E4's residual arithmetic is specified against these counts |

## 3. The Minors worth naming

- **The fixture generator could not tell two of its four arguments apart.** Measured: transposing
  the mate length with the gap between the mates gives the same 800 identifiers, the same 660
  covering two stretches and the same 1,460 stretches, because a mate starts at the sum of the two
  either way — so "30-base mates, a 40-base hole" was a claim nothing checked. The four are a
  named struct at the call site now, and one assertion **measures** the mate length: the pair that
  started at record 0 must cover records 0 to 29 and no further. Swapping them fails it.
- **The three headline figures were stated in a comment and asserted nowhere**, so drift to, say,
  81 % and 1,300 stretches would have passed with the comment silently wrong. 800, 660 and 1,460
  are assertions now, beside the fraction that states the intent. That also closes a circularity
  the review named: `stretches_of` and `derive_changes` are two implementations of the same set
  difference, so the arrival equality is a differential between them rather than a check against a
  number derived without the code.
- **The "more than four fifths" threshold appeared nowhere in the assertion it governed** —
  `× 10 > × 8`. It is a named constant carrying the fixture's measured 82.5 % and the corpora's
  83 % and 91 %.
- **The departure-count bound was relaxable by one** with all 228 tests green: the guard test used
  nine departures against a set of two, far from the boundary. `the_departure_count_may_equal_the_live_set_and_no_more`
  pins both sides — a whole set departing at once is ordinary, one more is damage named at the
  count.
- **The justification for "more than twice" argued the wrong way.** "A read straddling two walked
  regions is named twice" gives, under an allocator that never reuses an id, **two identifiers of
  one stretch each** — not one identifier of two. Only the spliced-alignment clause supports the
  test, and that is what the docstring says now.
- **Six hand-rolled write-and-read loops, three added by this step**, each a place the block
  restart could be forgotten independently. `round_trip` returns the arrival count now, and the
  three new tests read it off.
- Names: `alone` for a reader, where the same word already named a `Vec` of live sets;
  `all_of_them` and `all_the_stretches` as pronouns for counts whose nouns the failure message had
  to supply; `live_last_time`/`now` in the test helper where the production code says
  `previously_live`/`now_live`. And `pair_span - 1` is `saturating_sub`, so the helper cannot
  underflow if E3 reuses it with zero-length mates.

## 4. ⚠ A near miss of my own, and it is the one the plan warns about

Applying one of these fixes by replacing a **slice of the file between two markers**, I cut from
the sortedness test's docstring to the next `    }\n}\n` — which was the end of the enclosing
`proptest!` block, not the end of that test. **Three tests were silently deleted**:
`a_set_that_empties_and_fills_again_reads_back`,
`an_id_that_departs_and_arrives_in_one_record_is_damage` and `the_edges_of_the_range_read_back`.

The suite stayed green, because deleted tests do not fail. What caught it was the **count**: 29
down to 26. All three are restored and the count is back.

This is verbatim the hazard `plan-driven-implementation` records — *"a green suite is evidence
about the tree that was compiled, not about the tree in `git add`'s hands"* — arriving by deletion
rather than by an unreverted mutation. The habit that caught it is cheap and is now worth stating:
**after any edit that moves more than a few lines, compare the test count, not only the colour.**

## 5. What the review confirmed rather than found

- **3,000,000 fuzzed inputs across two independent harnesses**, no panic, no hang, and no
  allocation driven by a declared count. One built its corpus from paired-end coverage and refused
  to run unless more than half of each stream's identifiers re-enter; largest single case 30,390
  bytes against a 1 MiB ceiling, widest 923 µs against 50 ms. **Both harnesses were shown able to
  fail** — a planted panic, a planted 120 ms case and a planted 8 MiB allocation were all reported.
- **The step's central oracle holds the thing it is for.** Closing the hole between the mates —
  the fixture losing its re-entry — is caught by `most_reads_go_live_twice…` and by **no other test
  in the suite**. And the re-entry paths E1's review found dead are genuinely reached now:
  instrumentation counted the interleaving arm entered **1,320 times** by that one oracle.
- **Every number in the report, the commit message and the status entry re-derived by two agents
  independently, and all hold.**
- **The diff matched its stated intent exactly**: every hunk lands after `#[cfg(test)]`.
- **`errors`, `defaults` and `module_structure` are honestly clean**, checked over the whole file
  rather than the diff.
- **One mutation changed no behaviour and is reported as that**: applying arrivals before
  departures in the reader. They commute, because an arrival already live is refused, so the
  arriving set is disjoint from the pre-move live set.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::chain_ids` | 29 passed; 0 failed |
| `cargo test --lib ng::psp` | 230 passed; 0 failed |
| `cargo test --lib` | 4,767 passed; 0 failed; 14 ignored |

**Every defect the review found surviving, re-injected against the strengthened tests:**

| defect | tests failed |
|---|---|
| the writer does not restart at a block | 4 — **including the boundary test, which it survived** |
| the reader does not restart at a block | 3 — **same** |
| the reader keeps a file-scope memory of every read it has named | 6 — **same** |
| `LiveSetChanges::departed()` returns nothing | 1 |
| the departure-count bound is relaxed by one | 1 |
| `apply_arrivals` appends instead of merging | 4 |
| the fixture's mate length and gap are transposed | 1 |

Seven for seven.
