# Fixes applied — ng parameters file, A3

**Date:** 2026-08-28
**Review:** [ng_parameters_file_a3_2026-08-28.md](ng_parameters_file_a3_2026-08-28.md)

---

## 1. The Blocker, and the proof the fix works

Part 5 of `each_of_the_five_states_is_a_missing_key_and_not_a_value` asserted that no row of a
`Vec` matched a key the test had just `pop`ped out of that same `Vec` — a property of `Vec::pop`,
which no implementation of this module can make fail. **Parts 4 and 5 now assert against the
emitted document**, on the gap the fixture already has: three strata crossed with two slippage
groups is six pairs of which three carry a row, and one of the three strata carries a length
spectrum.

**Measured, in the container, from a checksummed copy restored to the same checksum
(`5289bbe9…`):** re-running the review's own mutation — a zero row for the (stratum × slippage
group) pair that put no reads — against the fixed test:

| | before the fix | after |
|---|---|---|
| `each_of_the_five_states_is_a_missing_key_and_not_a_value` | **passed** | **fails**, on `the pair with no reads has no row at all, rather than a row of zeros` |
| `the_whole_shape_emits_the_documented_toml` | failed | fails |

Before the fix the golden file was the only thing standing between the caller and a densified
slippage axis — and an author regenerating it would have accepted the change silently.

## 2. Findings table

| # | severity | status | note |
|---|---|---|---|
| B1 | Blocker | **Applied** | parts 4 and 5 read the emitted document; proof above |
| M1 | Major | **Applied** as documentation + a test | both unwanted spellings named in one paragraph, with C2 as the refuser |
| M2 | Major | **Applied with adaptation** | documented and pinned rather than made unspellable — see §3 |
| M3 | Major | **Raised at Checkpoint A** | the spec and the code disagree; not this step's to settle |
| M4 | Major | **Raised at Checkpoint A** | as above |
| m1 | Minor | **Applied** | C2 owns it, not C1 — C1 is parsing and both shapes parse |
| m2 | Minor | **Applied** | folded into the B1 fix |
| m3 | Minor | **Applied** | `replacen`, and the comment says why |
| m4 | Minor | **Applied** | `by_sample`'s emptiness named, with C2 as the refuser, matching its sibling |
| Nits | Nit | **Applied**, both | |

## 3. The one adaptation

M2 asked for `NonZeroU64` on the two evidence counts, which would make *a measurement carrying the
evidence of not being measured* unwritable. **Not taken.** A contamination fit that returns an
estimate with no evidence behind it would then have no file it could be written to at all, and
whether that can happen is a question about the estimator rather than about this shape — step C4's
round trip on a real fit is what answers it. Instead the state is documented where the type is, the
refusal is assigned to step C2, and
`the_shape_accepts_two_things_step_c2_must_refuse` pins that it is accepted today, so C2's landing
inverts a failing test rather than adding one nobody remembered to write.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **14 passed, 0 failed, 1 ignored** |
| `cargo test --lib` | **4,934 passed, 0 failed, 12 ignored** |
| `cargo doc --no-deps` | zero unresolved links in this module |
