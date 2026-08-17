# Fixes applied — ng cohort merge B1

*2026-08-17, against [the review](ng_cohort_merge_b1_2026-08-17.md) of plan step B1
([impl report](../implementations/ng_cohort_merge_b1_2026-08-17.md)). Every finding is answered
below; the code is [build.rs](../../../../src/ng/run/cohort_merge/build.rs) and one doc paragraph
in [close.rs](../../../../src/ng/run/cohort_merge/close.rs).*

## The shape of what changed

**Two arguments that were only meaningful together became one handle, and that closed three findings
at once.** `project_into(member, sequence, buffer)` let any member of the locus be paired with any
other member's sequence — a well-formed allele, no panic (M1). It is now

```rust
reference.placing(member).project_into(sequence, &mut buffer)
```

where `placing` checks the member's placement once and `MemberPlacement` carries the offset and the
width. The handle also owns `projectable_sequences()`, the member's complete sequences, so the loop
a caller writes cannot reach a partial (M2) and the per-member checks are no longer redone per
sequence.

**Nine tests were added, and each kills a mutant the review found surviving.** Re-run after the fixes,
each of the five survivors and both new assertions now fails exactly one test:

| mutation | result |
|---|---|
| `span_of(member.region)` → `member.region.len()` in the projection | FAILED, 58 passed / 1 failed |
| `NOT_COVERED` `0` → `b'N'` | FAILED, 58 / 1 |
| drop the projection's contig assertion | FAILED, 58 / 1 |
| drop the projection's reach assertion | FAILED, 58 / 1 |
| gather only each sample's first observation | FAILED, 58 / 1 |
| drop the new overlap-disagreement assertion | FAILED, 58 / 1 |
| drop the new `Verdict::Build` assertion | FAILED, 58 / 1 |

(`./scripts/dev.sh cargo test --lib ng::run::cohort_merge`, one mutation at a time from a pristine
copy; driver `tmp/mutate_b1.sh`.)

## Findings, one by one

| id | what it said | what was done |
|---|---|---|
| **M1** | a sequence is not tied to its member | **Applied.** `LocusReferenceBases::placing(member) -> MemberPlacement`, with `project_into(sequence, buffer)` on the handle. |
| **M2** | the `Partial` panic fires on ordinary data | **Applied, differently.** The newtype the review proposed lives in `locus_generation`, one module out, and would ripple; the handle's `projectable_sequences()` closes the same hole where the API is — a caller looping through it cannot reach a partial. The assertion stays as the backstop for a caller that reaches into `observations` itself, and both paths are tested. |
| **M3** | the ceiling test pins the gather, not the projection | **Applied.** The test now projects the member ending *at* the ceiling as well, which is the one `span_of` and `len()` disagree about; the doc says which member does the work and why the other cannot. |
| **M4** | nothing pins the sentinel's one property | **Applied.** `a_reference_containing_n_bases_is_gathered_like_any_other`, and the constant's doc now says why `N` in particular would be the wrong sentinel. |
| **M5** | the projection's contig guard is untested | **Applied.** `projecting_a_member_from_another_contig_is_refused`, on the `placing` path. |
| **M6** | no sample ever has two observations | **Applied.** `one_samples_two_observations_both_reach_the_gather`, gathering `ACGTA` from one sample's 10–11 and 12–14. |
| **M7** | the reference-width check is the producer's guarantee, not the pairing's | **Applied as documentation**, which is what the finding asked for: `over`'s doc now separates it from the other three and instructs the psp step to make it a `RunError` beside `ObservationExceedsReachCeiling`. Nothing to move today — `RunError` does not exist yet. |
| **Mi1** | the uncovered message names no position | **Applied.** It now names the position and the offset; the test asserts `leave position 13 uncovered`. |
| **Mi2** | the sample reads as a name | **Applied.** "sample index {}". Naming the accession needs the run's sample table, which `SampleMembers` does not carry — left for the caller objects. |
| **Mi3** | `over` accepts any verdict and allocates first | **Applied.** `assert_eq!(locus.verdict, Verdict::Build, …)` is the first statement, which states spec §3.2 in code and is what bounds the allocation; the doc's cost claim now rests on it. |
| **Mi4** | `expect`s without `// PANIC-FREE:`, one with a false invariant | **Applied.** Three comments; the first no longer claims a contig bounds a region's span — it claims what is true, that `usize` is 64 bits here. |
| **Mi5** | overlapping members disagreeing is silent | **Applied.** The gather now writes base by base and refuses a disagreement, naming the position and both bases. Costs one comparison per base of a span the verdict bounds. This replaces the report's §7 "recorded rather than checked". |
| **Mi6** | the projection's reach guard is untested | **Applied**, with the matching start-before test. |
| **Mi7** | `bases()`'s reference-allele claim is untested | **Applied.** `a_member_that_matched_the_reference_projects_to_the_locus_reference`, at a non-zero offset. |
| **Mi8** | prose contradicts the test's own assertion | **Applied.** "neither member alone carries all five of those bases" — and the two fixture-only length assertions the nits flagged are gone with it. |
| **Mi9** | the one-base test cannot fail for its stated reason | **Applied.** It now asserts the gathered reference too. |
| **Mi10** | `close.rs`'s `span_of` doc invites the surviving mutation | **Applied.** The exception is removed: nothing in the module reaches for `len()`, observations included. |
| **Mi11** | `member` names two types | **Applied.** The gather's outer binding is `sample_members`; `member` means one `SampleLocusObservations` everywhere. |
| **Mi12** | `LocusReference` reads as a borrow | **Applied.** Renamed `LocusReferenceBases`, matching the field it is gathered from. |
| **Nits** | eleven | **Applied**, except the two below. |

**Two nits not applied, with reasons.**

- **`NOT_COVERED` renamed to `UNCOVERED_BYTE`.** Kept: its two uses read as a state
  (`vec![NOT_COVERED; span]`, `*slot == NOT_COVERED`), and the doc now says what byte it is and why
  `N` would not do, which is what the nit was reaching for.
- **A shared `assert_inside(outer, inner, …)` helper** for the contig and reach checks, now written
  twice. The reviewer filed it as a note-not-a-finding at two occurrences; with the handle, the
  second pair runs once per member rather than once per sequence, and the messages differ in what
  they name. Left for whoever adds the third.

## Claims corrected

Both came from the claim-verification pass, and both were mine about my own work:

- **"the failure spec §3.3 rejects clamping for"** — §3.3 rejects clamping for three other reasons.
  "An allele no molecule carried" is the *next* bullet, about splitting an event. The doc comment no
  longer attributes it.
- **"Four checks are release-level assertions"** — there are eight assertion sites, and the gloss
  ("they fire on a caller pairing a locus with somebody else's members") never fitted the witness
  check. The report now counts them and separates the two kinds.

## Validation

In the container, after the fixes:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean (`Finished dev profile … in 3.38s`).
- `cargo test --lib ng::run::cohort_merge` — `ok. 59 passed; 0 failed` (50 at review time, 38 before
  the step; 21 of the 59 are `build.rs`'s).
- `cargo test --lib` — see the commit.
