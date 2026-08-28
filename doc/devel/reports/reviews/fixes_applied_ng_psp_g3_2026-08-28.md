# ng psp store — G3 fixes applied

*2026-08-28, on top of `b9c05193`. Answers [the G3 review](ng_psp_g3_2026-08-28.md).*

---

## 1. The two Blockers

**B1 — the trailer's offset had no lower bound that meant anything.** `replace_trailer` now reads
the header for its length alone and refuses a footer whose `trailer_offset` falls inside it. The
probe file — a 3,742-byte psp with `trailer_offset = 4` — is refused with every byte still there,
and `a_footer_that_puts_the_trailer_inside_the_header_is_refused` is that probe.

**This reverses the step's own §2.1 decision, and the reversal is the point.** Spec §6.7's table
lists two refusals for this operation and the argument was that reading more would cost. Reading
the header costs one small parse and is what stands between a file with nonsense offsets and a
fresh footer blessing it. **It also earns a refusal class §6.7 does not list for this operation**:
a file written by a newer format now comes back as `UnsupportedVersion` rather than being
rewritten. That is the safe answer, and **spec §6.7's table should gain the row** — recorded here
rather than edited, because the table is the owner's.

**B2 — the file is truncated at the trailer's offset *before* anything is written.** From that
instant until the new footer lands the file has no footer, so an interruption leaves something
every reader refuses, which is what the ⚠ claimed and did not deliver.
`a_rewrite_stopped_part_way_leaves_a_file_no_reader_accepts` walks every stopping point at four
payload lengths and asserts that exactly one — the complete write — opens.

⚠ **That test builds the torn files itself rather than interrupting a real write**, so what pins
the *ordering* is a behavioural test: deleting the leading truncation makes
`a_trailer_of_any_length_replaces_the_one_before_it` fail on the shorter payload. The two together
are the guarantee; neither alone is.

## 2. Finding by finding

| finding | what was done |
|---|---|
| **M1** the one-sided sections fixture | replaced by `sections_that_stop_short_or_run_past_the_footer_are_both_refused`, nudging both ways and asserting the file is untouched each time — the reader's own shape, which is now also literally the same rule (M4) |
| **M2** the 47-byte foreign file | split into `a_file_shorter_than_a_footer_is_refused_as_incomplete` and `a_foreign_file_longer_than_a_footer_is_told_apart_from_a_killed_run`, each asserting the variant and the file's bytes. **The second changes §2.4's decision**: sharing `open`'s footer read brought `open`'s distinction with it, and a foreign file now comes back as `NotAnNgPsp` — *you handed me the wrong file* and *rebuild this one* are different instructions |
| **M3** the unheld causes | assertions added where the variants are actually reached — in `writer.rs`'s own tests, which is where they are raised. `RecordRefused`'s cause is matched to `BlockWriteError::ContigOutOfOrder`; `WouldNotBeReadable`'s to `NotReadable::BlockIndex` with its `reason`; and the lost-records case is asserted to carry **no** cause, which is the other half of the rule. `ManifestRefusal::Compressor` keeps its variant with a ⚠ saying nothing reaches it today and why |
| **M4** the copied footer read | `read_and_check_the_footer` is now one function in `reader.rs` that `PspReader::open` and `replace_trailer` both call. The copied reasoning is gone with the copy |
| **M5** retry versus rebuild | **documented, not built.** The doc names which failures leave the file byte-identical and which leave it unreadable. A field a caller could match on is what the review asked for; the caller — `run_streaming.md`'s cohort driver — does not exist yet, and the field would touch every `Io` site in the module |
| **Minor** the doubled chain | the five `BlockReadError` / `BlockCutRuleError` messages that interpolated `{source}` while also marking it `#[source]` no longer do. This was deferred at G2 as pre-existing; G3 is where it became visible in a chain the step itself created |
| **Minor** the wrong trim mechanism | gone with the comment, replaced by the truncate-first ⚠ and its measured torn-state counts |
| **Minor** `retrailer.rs` | renamed `trailer.rs` |
| **Minor** arch §1 | `trailer.rs` is in the tree |
| **The two wrong numbers** | corrected in the G3 implementation report, with what was actually run beside what the label claimed |

## 3. Verification

In the container, on the tree being committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib ng::psp` — **373 passed, 0 failed**, against 370 before the fixes and 361
  before the step.
- `cargo test --lib` — **4,910 passed, 0 failed, 14 ignored**.

**Seven defects injected against the fixed tree, seven caught**:

| defect | caught by |
|---|---|
| the truncation happens last rather than first | 1 |
| the rewrite starts at the index rather than the trailer | 4 |
| the footer keeps the old trailer length | 4 |
| the sections rule accepts sections that stop short | 3 |
| the trailer-past-the-header check is dropped | 1 |
| a file shorter than a footer saturates rather than refusing | 3 |
| the index decoder's own account is dropped from the cause | 1 |

The sections rule and the short-file rule are now caught by tests in **both** modules, which is
what sharing the function bought.
