# ng psp store — F2 fixes applied

*2026-08-28. Applies [the F2 review](ng_psp_f2_2026-08-28.md) to the tree at `07f6273f`,
branch `ng-psp-encoding`. Every finding is addressed; none is deferred.*

---

## 1. Finding by finding

| | finding | what was done |
|---|---|---|
| **M1** | no accepted footer ever had a trailer holding anything | `a_footer()` carries **1,288 trailer bytes**, and the round-trip set keeps an empty trailer as a second case rather than the only one. The widest-value fixture gained a non-zero trailer too. Both wrong abut rules now fail: folding `trailer_bytes` into the comparison fails 3 tests, checking only when the trailer is empty fails 1. The helper's doc records what it carried and what that cost. |
| **M2** | the destructure forces a field to be mentioned, not written | `the_footer_constant_is_the_width_of_the_fields_it_stands_for` **destructures** instead of reading fields off a value, which is the one place the constant is tied to the field *set*. Verified by the reviewer: with a seventh field added, that test stops compiling, and the only minimal repair that compiles turns it red at 48 against 52. |
| Minor | messages did not carry the damage instruction | both arithmetic refusals now end "this psp is damaged" / "so this one is damaged". |
| Minor | the two overflow checks could report each other's section | `the_first_section_to_overflow_is_the_one_named` — a footer where **both** overflow, asserting the index is named, which pins the order as well as the labels. |
| Minor | nothing tested that arbitrary bytes never panic | `no_forty_eight_bytes_make_the_decoder_panic` — 20,000 deterministic draws, half carrying the real magic so the arithmetic past the magic check is actually reached. The generator is three lines and seeded, so a failure is reproducible from the seed. |
| Minor | the byte-coverage test compared against the wrong thing | it compares each flip against **what the pristine bytes decode to**, not against the footer they came from. The mutation it could not see before — masking the checksum's top byte — now fails 3 tests. |
| Minor | `FileSection` public without `#[non_exhaustive]` | added, matching `IndexEntryField` next door. |
| Minor | `SectionRunsOffTheEnd` named a check the function cannot make | → `SectionEndIsPastAnyFile`. It fires on `u64` overflow, and the doc twelve lines above says nothing here knows the file's length; F4 has to map this variant and the two readings led to different mappings. |
| Minor | *abut* appears in no design document | → `IndexDoesNotEndWhereTheTrailerBegins`, which is the spec's own phrasing, and the message says it in words. |
| Minor | `FileSection::Index` spelled it three ways | → `BlockIndex`, matching the spec's vocabulary table and the `Display` that already rendered it so. |
| Minor | a file shorter than 48 bytes is unrepresentable here | said on the function, routed to `open`. |
| Minor | the two `expect`s were unmarked | **kept, and marked `PANIC-FREE` with the reason.** The reviewer's argument is right and worth recording: this function's argument is a fixed-size `[u8; FOOTER_BYTES]`, so every window has a length the compiler knows — unlike `index.rs`, where the length was a runtime fact and the same shape was replaced after a mutation turned eight tests into panics. |

## 2. The three wrong numbers

**N1 is the one that matters, because it was a wrong *mechanism*, not a miscount.** The report
and the commit message said moving the magic check after the offset arithmetic is caught by three
tests. It is caught by **one**. The mutation that produced the 3 deleted the check outright,
making it a duplicate of the row above it — so a defect table with ten rows described nine
defects. Re-run with the check genuinely relocated: one test.

**This is the third time in Milestone F that a mutation of mine was not the defect I labelled
it.** The first was an unreachable arm appended after the arm that already matched; the second a
`Vec::new()` that over-allocates nothing either; this one a deletion labelled as a move. All
three were reported before I looked at what they actually did. **The habit that catches it is
asserting the anchor matched and then reading the mutant, not the label** — and it caught two
more no-op mutations of mine during this very fix pass, where a replacement string had not
matched at all and both the "fix" and its "verification" were silently nothing.

N2: four names re-exported, not five. N3: the width is asserted four ways, not three. Both
corrected in the F2 report.

## 3. Verification

Seven defects injected into the fixed tree, **seven caught**, each with its anchor asserted to
have matched before the run:

| defect | caught by |
|---|---|
| the abut rule folds the trailer's length in | 3 |
| the abut rule runs only for an empty trailer | 1 |
| the overflow checks report each other's section | 1 |
| the magic is checked after the arithmetic (genuinely moved) | 1 |
| the decoder masks off the checksum's top byte | 3 |
| a seventh field is added to `Footer` | does not compile — which is the guard |
| *(the ten from F2's own table re-verified by the review: nine scored exactly as claimed)* | |

Gate, in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — **library 4,824 passed, 0
  failed**, against 4,822 at `07f6273f`. `ng::psp::footer` is **16 tests**, against 14. The 21
  example failures are the known pre-existing breakage.
