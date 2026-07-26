# ng read preparation — A2+B1+B2: `LeftAlignPreparer`, the whole v1 transform

**Date:** 2026-07-26 · **Plan:** [read_preparation.md](../../ng/impl_plan/read_preparation.md) steps A2,
B1, B2 (**merged — see §7**) · **Spec:** [read_preparation.md](../../ng/spec/read_preparation.md) §5,
§6, §7 · **Arch:** [read_preparation.md](../../ng/arch/read_preparation.md) §2, §3

Combined implementation + review + fixes report for one loop iteration.

## 1. Plan

Land step 2's v1 implementation: the preparer type and its scratch, the indel predicate that decides
whether a read needs the reference at all, and the transform itself — fetch the window, round-trip
the CIGAR through an `AlignmentNormalizer`, and build production's `PreparedRead`.

## 2. Assumptions and choices the plan left open

- **A zero-reference-span read is passed through unnormalized.** An all-insertion CIGAR consumes no
  reference, so there is nothing to align against and a zero-length window would be a question with
  no answer. Production guards the same case — `fetch_raw_slice` returns `None` when `ref_span == 0`,
  which skips `F3` — so this keeps the two in step. Taken from the oracle, not invented.
- **The window-length check is a `debug_assert_eq!`, not a runtime error.** The arch doc asks for
  "check the length, do not trust it"; `RefSeq::fetch_into`'s contract already guarantees it writes
  exactly `length` on success or returns `Err`, so a runtime check would defend against a broken
  `RefSeq` impl rather than a reachable state. The assertion documents the invariant where the silent
  failure would occur.

## 3. Changes made

All in **[src/ng/read/left_align.rs](../../../../src/ng/read/left_align.rs)**:

- `cigar_has_indel` — **exactly production's `is_indel`** (`Insertion | Deletion`). It is what decides
  whether a read fetches at all, so anything narrower would silently skip normalization for a read
  production would shift.
- `LeftAlignScratch { reference_window }` — the per-worker buffer. Reads with no indel never touch it.
- `LeftAlignPreparer<R: RefSeq, N: AlignmentNormalizer = DefaultAlignmentNormalizer>` — the reference
  and the normalizer as **visible type parameters**, with `with_default_normalizer` for callers and
  `new` for the bake-off. Deliberately not `Clone`/`Copy`: no real reference accessor is either, and
  the derive would imply an ownership model the arch doc rejects (one preparer per worker).
- `canonicalize` — the transform: skip on zero span, fetch the **uppercased** window, assert it is
  full, move the CIGAR into an `Alignment { reference_offset: 0, .. }`, normalize, move it back.
- `into_prepared` — delegates the whole field build to production's `prepare_passthrough`, per spec
  §11's leaning (call it; port only if its `qual.clone()` profiles).

## 4. Tests added (9)

| test | what it proves |
|---|---|
| `a_read_with_no_indel_never_fetches_a_reference_window` | **the conditional fetch, pinned** — the reference used is one whose `fetch_into` *panics*, so this asserts the fetch is not reached, which inspecting the output could never show |
| `an_indel_in_a_homopolymer_is_moved_to_its_leftmost_spelling` | the transform works: `TAAAAG`/`TAAAG`, `M4 D1 M1` → `M1 D1 M4`, `alignment_start` unmoved, bases untouched |
| `the_indel_predicate_is_exactly_production_s` | the whole **negative** set (`SoftClip`, `Skip`, `HardClip`, `Padding`, `SeqMatch`, `SeqMismatch`) is not an indel |
| `a_failed_fetch_is_an_error_not_a_decline` | `Err`, never `Ok(None)` — the §7 split |
| `an_all_insertion_read_is_passed_through_without_fetching` | the zero-span guard, again against the panicking reference |
| `base_qualities_pass_through_uncapped` | v1's whole quality handling (BAQ deferred) |
| `a_fresh_scratch_holds_no_window` | `Default` pre-allocates nothing — it decides nothing |
| `the_default_constructor_holds_its_reference`, `an_alternative_normalizer_can_be_named` | both constructors; the second keeps the bake-off surface open |

## 5. Validation

In the container (`./scripts/dev.sh`), after `cargo fmt`:

- `cargo fmt --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib ng::read::left_align` — **9 passed; 0 failed.**
- `cargo test --all-targets --all-features` — see the commit message for the recorded result.

## 6. Review and fixes applied

Reviewed inline against the `rust-code-review` checklists (mechanism agreed with the owner). **No
Blocker, no Major.** Two Minor, both applied:

1. `canonicalize` had landed in a **second, redundant inherent `impl` block**; merged into the one
   above it.
2. `LeftAlignScratch.reference_window` was `pub(super)` where private suffices — the tests are a child
   module and reach it regardless. Visibility minimised.

## 7. Deviation: three plan steps merged into one commit

**What happened.** The plan splits A2 (types), B1 (predicate + no-indel path) and B2 (indel path),
with B2 marked "own commit, do not bundle" because its failure mode is silent. Committing A2 alone
fails `cargo clippy -- -D warnings`:

```
error: field `reference_window` is never read
error: field `normalizer` is never read
```

Every field whose first reader lives in a later step is a hard error under `-D warnings`, and once
the test-only accessor became `#[cfg(test)]` that was **all three** fields. Keeping the plan's
boundaries would have required `#[allow(dead_code)]` on the entire struct — which is a louder signal
that the boundary is wrong than that the lint is, and the skill is explicit that checks are not to be
relaxed to get a step committed.

**What was preserved.** B2's isolation requirement exists so `git bisect` lands on the transform
rather than on an unrelated neighbour. This commit contains **only step-2 code**, so it still does.
What is lost is the finer types/transform split, which was never the risky boundary.

**The alternative, if the literal split is wanted:** three `#[allow(dead_code)]` attributes in A2,
removed by B2. Recorded here so the choice can be reversed cheaply.

## 8. Tradeoffs and follow-ups

- **An indel-bearing read whose window runs past the contig end ends the run** (`Err`), where
  production skips `F3` and emits the read un-left-aligned. That is the deliberate divergence spec §7
  records; **C1's parity fixture must exclude such reads** or assert ng's abort.
- **Left-alignment is not yet proven against production** — that is C1. What is proven here is that
  the transform moves an indel, keeps the placement, and leaves bases and qualities alone.
- The `OPEN:` items from arch §6 are untouched: `Ok(None)` carries no reason, and call-vs-port of
  `prepare_passthrough` stands at "call it".
