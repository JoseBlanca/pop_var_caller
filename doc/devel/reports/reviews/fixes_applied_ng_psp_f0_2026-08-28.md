# ng psp store — F0 fixes applied

*2026-08-28. Applies [the F0 review](ng_psp_f0_2026-08-28.md) to the tree at `10793f41`,
branch `ng-psp-encoding`. Every finding is addressed; none is deferred.*

---

## 1. What changed, and why the change got smaller

**The key in the file is exactly what the owner ruled on. Everything else moved.**

Three renamings and one deletion carry most of it:

- **`FieldCardinality` → `FieldShape`, `OneValue`/`AList` → `Scalar`/`List`, wire key
  `cardinality` → `shape`, spellings `one-value`/`a-list` → `scalar`/`list`** (M1, M2). These
  are production's own type name, variant names and wire tokens for this exact concept
  ([`src/psp/registry.rs:60`](../../../../src/psp/registry.rs)). The word `cardinality` goes
  back to meaning what production means by it — *how often a field appears* — which ng's
  manifest does not carry.
- **`FieldSpec` loses the field it gained** (Mi1). Shape is derived: `FieldSpec::shape()` defers
  to `FieldEncoding::shape()`. Two consequences: a `FieldSpec` can no longer hold a
  contradiction, so the writer cannot produce one rather than merely refusing to; and the check
  moves to the decode path, which is the only side where two accounts exist. **`FieldSpec::new`
  is gone with it, and its twenty-four call sites are struct literals again** — that churn was
  only ever needed because a field had been added.

The diff against `10793f41` is smaller than `10793f41` itself was.

## 2. Finding by finding

| | finding | what was done |
|---|---|---|
| **B1** | nothing states what an encoding's shape *is* | **`every_encoding_lays_down_the_shape_its_bytes_have`** names all eight answers against literals. The review's own mutation — `SignedVarint`, `FixedWidthInteger` and `IeeeFloat` moved to the list arm, which left the suite green at 255 — is now caught, **by this test and no other**. `FieldEncoding::shape`'s doc says why it exists. |
| **M1** | `cardinality` is production's word for the other half | renamed throughout to `shape`, with production's `scalar`/`list` tokens. `FieldShape`'s doc names production's `Cardinality` explicitly and says the two must not be confused, giving the `.psp`-extension collision as the reason. |
| **M2** | the doc contradicts the code | rewritten around *what one appearance looks like*. The doc now uses `mapq-sum` as its worked example: a scalar, of which a five-observation record holds five. |
| **M3** | hand-written array, fake completeness test | `FieldShape`, `ALL_SHAPES` and `spelled()` are **generated from one source** by a `field_shapes!` macro, so a variant cannot reach one without reaching all three. The replacement test checks what generation cannot — that no two shapes share a spelling, and that every spelling in the list parses back. |
| **Mi1** | stored where it should be derived | field dropped; `FieldSpec::shape()` derives; `check_declared_shape` runs on decode only. |
| **Mi2** | the fixture omits both list-shaped encodings | `a_manifest()` now declares **all eight**, and its doc records that it carried six and omitted exactly the two that mattered. `the_header_text_names_each_fields_shape` asserts both spellings appear and that the per-field line count holds. |
| **Mi3** | only pinned against a spelling resembling neither | the test now refuses **eight** spellings: `per-observation`, `Scalar`, `LIST`, `" list "`, `"list "`, `s`, `scalars`, `""`. Case folding and trimming are both caught. |
| **Mi4** | the message can be inverted with its tests green | both tests assert the **whole sentence**, not `contains` of its ingredients. The message says `field "x" declares shape "scalar", but its encoding "chain-id-list" lays down "list"`. |
| **Mi5** | passes with `#[serde(default)]` added | renamed to `a_field_that_declares_no_shape_is_refused_by_the_parser` and asserts `missing field \`shape\`` and that the message points at `[[manifest.field]]` — the line number a defaulted key would lose. |
| **Mi6** | the fixture name contains the refused value | the field is named `a-field`, and the assertion is the exact message. |
| **Mi7** | a field name can forge a manifest line | `check_manifest` refuses a field name holding whitespace or a control character, mirroring the contig-name rule and citing it. `a_field_name_that_could_forge_a_line_in_the_header_is_refused` covers the newline forgery, a space, a tab and a NUL. |
| **Mi8** | `from_manifest` ignores shape, silently | its doc now says it needs no shape check: `Header::decode` has already refused a disagreeing file, and the walk steps over an unknown field by its encoding — the same source the shape derives from. |
| **Mi9** | the report cites untracked scratch scripts | §4 of the F0 report describes each mutation precisely enough to reapply, and the filenames are gone. |
| **Mi10** | `#[non_exhaustive]`; `pub` on an internal method | `#[non_exhaustive]` kept — inert inside the crate, and `FieldEncoding` carries it, so dropping it on the type read beside it would be the surprise. `shape()` stays `pub` on `FieldEncoding` because `FieldSpec::shape()` is public and defers to it. |
| **Nit** | `KNOWN_FIELD_COUNT` says "the head's four" | corrected to five. |

## 3. Verification

Six defects injected into the fixed tree, **six caught** — each applied to a clean copy, the psp
suite run, the failing tests counted (the script is regenerable from these descriptions alone):

| defect | caught by |
|---|---|
| `SignedVarint`, `FixedWidthInteger`, `IeeeFloat` moved from the scalar arm of `FieldEncoding::shape` to the list arm | 1 — `every_encoding_lays_down_the_shape_its_bytes_have` |
| `shape_of`'s comparison becomes `eq_ignore_ascii_case` on a trimmed value | 1 — `an_unknown_shape_is_refused_and_the_message_lists_the_known_ones` |
| `check_declared_shape`'s message swaps its declared and laid-down operands | 2 |
| the field-name whitespace rule is replaced by `if false` | 1 — `a_field_name_that_could_forge_a_line_in_the_header_is_refused` |
| `WireHeader::from` writes the literal `"scalar"` for every field | 12 |
| `check_declared_shape` compares the declared shape with itself | 2 |

**⚠ The first attempt at the Blocker's mutation was a no-op and was reported as a survivor
before I looked at it.** It appended a `#[allow(unreachable_patterns)]` arm *after* the arm that
already matched those three encodings, so the earlier arm won and nothing changed. Rewriting it
as a *move* rather than an addition produced the real mutant. **A mutation that does not compile
into a behaviour change is not evidence about a test** — the dispatch rule that says so
(*prove it changed behaviour before recording a survivor*) applies to the author too.

Gate, in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — **library 4,792 passed, 0
  failed**; `ng::psp::` **255** of them, against 253 at `10793f41` and 246 before F0. The 21
  example failures are the known pre-existing `ref.fa.repeats.parquet` breakage.

## 4. Not done, and where it went

- **`ALL_ENCODINGS`' identical drift flaw and its fake completeness test** — pre-existing,
  Major, demonstrated by the reviewer. `FieldEncoding`'s variants carry parameters, so the
  one-source generation the shape half now has is a different macro and a change of its own.
  Recorded in the review's *Out of scope observations*.
- **Whether the manifest should also carry how often a field appears** — the owner's, raised at
  Checkpoint F.
