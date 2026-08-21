# Applying the B1 + B2 review — ng calling foundations

*2026-08-21. Branch `ng-calling-foundations`. Input:
[`ng_calling_b12_2026-08-21.md`](../reviews/ng_calling_b12_2026-08-21.md). Every finding in that
report is accounted for below.*

## Findings table

| id | severity | decision | status |
|---|---|---|---|
| M1 — the overflow guard is untested and cites another caller's cap | Major | Apply | **Applied** |
| M2 — the copies constructor refuses only the empty vector | Major | Apply | **Applied** |
| M3 — the parallel-length contract is enforceable nowhere | Major | Apply | **Applied with adaptation** — the constructor, not B3 (§ below) |
| M4 — the emptiness guard checks the harmless property | Major | Apply | **Applied** (subsumed by M3) |
| Mi1 — `LocusKind`'s payload variant untested | Minor | Apply | **Applied** |
| Mi2 — the out-of-range test uses a one-allele table | Minor | Apply | **Applied**, and then again (§ below) |
| Mi3 — a zero-length reference allele is accepted | Minor | Apply | **Applied** |
| Mi4 — no property test for the admit/lookup round-trip | Minor | Apply | **Applied** |
| Mi5 — `alleles()` leaks the storage and hides the count | Minor | Apply | **Applied** (`iter()` + `len()`) |
| Mi6 — the prune has no route through the API | Minor | Apply | **Applied** (documentation; the method is the loop plan's) |
| Mi7 — the module doc's headline claim is inverted | Minor | Apply | **Applied** |
| Mi8 — "three things it drives" against a four-item list | Minor | Apply | **Applied** |
| Mi9 — `push` names the `Vec` operation, not the domain event | Minor | Apply | **Applied** (`admit`) |
| Mi10 — `get` is the vaguest name of the three accessors | Minor | Apply | **Applied** (`bases_of`) |
| Mi11 — the cohort sum is named with the spec's per-sample term | Minor | Apply | **Applied** (documentation) |
| Minor — four terms undefined before they carry arguments | Minor | Apply | **Applied** |
| Minor — `&[Box<[u8]>]` called a "flat view" | Minor | Apply | **Applied** |
| Minor — `kind` public while `alleles` is private | Minor | Apply | **Applied** (private + `kind()`) |
| Nits — list maintenance, repeated argument, single-allele fixture, overstated comment | Nit | Apply | **Applied** |
| Out of scope — `LocusKind` belongs in `types.rs` | — | Defer | **Deferred** (§ follow-ups) |
| Out of scope — `ng/mod.rs`'s partial `pub use` list | — | Dispute | **Won't fix** — judged correct by the reviewer |
| Out of scope — `ng/mod.rs`'s five-clause module doc | — | Defer | **Deferred** |

## M1 — the overflow guard

The `expect` stays: exceeding what an `AlleleId` can name is a caller bug, not a data condition,
and the alternative spelling `len() as u16` wraps at 65,536 to `AlleleId(0)` — the reference — which
is silent corruption rather than a failure. What changed is the doc and the test.

**The doc's reachability claim was wrong and is gone.** It cited production's
`DEFAULT_MAX_ALLELES_PER_RECORD` and `MAX_ALLELES_PER_VAR_CAP`, which belong to a caller ng shares
no code with; the reviewer's `grep` found no allele cap anywhere in `src/ng/`. It now says so: ng
has no cap yet, step 6 is where one will live, and until then this check is the only thing between
a pathological locus and an id that wraps onto the reference.

**Verified.** With the checked conversion replaced by `len() as u16`:

```
test ng::calling::tests::admitting_past_the_widest_table_an_id_can_name_is_refused - should panic ... FAILED
test result: FAILED. 10 passed; 1 failed; 0 ignored; 0 measured; 3962 filtered out
```

## M2 and M3 — the copies constructor, and where the length check lives

`ExpectedAlleleCopies::new` now takes the allele table and checks two things:

```rust
pub fn new(copies: Vec<f64>, alleles: &CandidateAlleles) -> Self
```

one entry per allele, and every entry finite and at or above zero.

**The two reviewers disagreed about where the length check belongs, and the constructor won.**
`errors` wanted it here; `reliability` wanted `len()` on both types and one `assert_eq!` in B3's
`LocusInference`, on the ground that a constructor taking `&CandidateAlleles` "makes the copies
type depend on the allele table for no other reason".

That objection does not survive its own premise: **the dependency is the type's definition, not an
addition to it.** The doc comment says entry *i* is allele *i* and that the two must stay the same
length; taking the table writes down a relationship already asserted in prose. Enforcing it at
construction makes the mismatch *unrepresentable* rather than *detected*, which is the project's
stated design rule and the standard the A2 and B2 reviews have both been applying. B3's check would
catch only pairings B3 assembles; a consumer handed a bare `ExpectedAlleleCopies` from anywhere else
would be back where it started. Both halves of `reliability`'s remedy were taken anyway: `len()` and
`is_empty()` are on both types, because consumers want them and B3 reads them.

The emptiness `assert!` is gone, subsumed: a table always holds at least the reference, so a copies
vector matching its length cannot be empty.

**Verified.** A negative and a `NaN` are each pinned by their own `#[should_panic]` test. The
mutations the review found surviving — a clamp at zero, a scrub of non-finites — are now
*unreachable*, because the validation rejects those inputs before any repair could act on them; that
makes them no-ops rather than survivors, and reporting them as survivors would be wrong. The
mutation that *is* still reachable is the same repair applied **before** the validation, which is the
realistic slip:

```
test ng::calling::tests::expected_copies_reject_a_negative_count - should panic ... FAILED
test ng::calling::tests::expected_copies_reject_a_count_that_is_not_a_number - should panic ... FAILED
test result: FAILED. 10 passed; 2 failed; 0 ignored; 0 measured; 3962 filtered out
```

## Mi2 — the review's own fix was not enough, and the second attempt is measured

The reviewer's remedy for the narrowed-id mutation was a three-allele fixture. **It does not kill
it.** Applied as given, with the id truncated to a `u8` inside the lookup:

```
test result: ok. 11 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out
```

The mutation needs a table wider than 255 alleles to differ at all, which the reviewer said in the
finding's own honest caveat while proposing a fixture three alleles wide. So the test that landed
builds a 300-allele table in which **each allele spells its own id in bases**, and checks that ids
either side of 256 resolve to themselves:

```
test ng::calling::tests::an_id_above_a_bytes_worth_resolves_to_its_own_allele ... FAILED
assertion `left == right` failed: allele 256 did not resolve to itself
test result: FAILED. 11 passed; 1 failed; 0 ignored; 0 measured; 3962 filtered out
```

The three-allele arm was kept as well — it is where an id one past the end of a real table is
asked for — but it is the wide one that carries the finding.

## Mi1, Mi3, Mi4 and the naming changes

- The kind test now builds all three `LocusKind` variants, and checks the `Ssr` payload in full —
  motif and both flanks. **Verified:** a mutation emptying the left flank inside `new` now fails
  that test where it previously survived the whole suite.
- `new` refuses a reference allele spelling no bases, with a `#[should_panic]` test. ng's own
  `reference_bases` is documented as never empty, and this repository's VCF encoder does not reject
  an empty `REF`, so the alternative is an unparseable record rather than a crash.
- The property test `every_id_admit_mints_resolves_back_to_the_allele_that_was_admitted` went into
  the module, in the file style step A1 established: arbitrary admission sequences up to 40 alleles,
  asserting ids are dense from 1, each resolves to its own bases, the reference never moves, and the
  first id past the end resolves to nothing.

Renames and documentation: `push` → **`admit`**, the architecture's own word for how a discovered
sequence enters the table; `get` → **`bases_of`**, which says what comes back; `alleles()` →
**`iter()`** yielding `&[u8]` plus **`len()`**, so the storage stays private and the count has a
name; `kind` is private with a `kind()` accessor, so a table cannot be re-routed to another read
model after its alleles were chosen by the right one.

The module doc's inverted claim is corrected — nothing in the tree forbade the sideways imports,
which is precisely why three architecture documents each had to write a no-import rule by hand, and
the one folder is what makes those rules unnecessary. Its opening now says what the module does
before naming the step numbers, and the four items are four. *Row builder* and *the prior's seed*
are defined where they carry the argument for `kind` existing. `ExpectedAlleleCopies`' doc opens by
distinguishing the cohort sum from a sample's own copies, because the prior's leave-one-out term
subtracts the second from the first and the two have the same shape.

## The two not applied

- **`LocusKind` moving to `src/ng/types.rs`** — deferred. The reviewer established it is not a
  back-reference (the rule is peer-to-stage, and both are pipeline stages in pipeline order) but
  that ten files across four stages consume it and `module_layout.md` already lists it among the
  shared vocabulary. Pre-existing; this diff adds the tenth consumer. Moving it is a mechanical
  change across ten files and does not belong inside a step that adds two types.
- **The calling types in `ng/mod.rs`'s `pub use` list** — won't fix, and the reviewer agreed: that
  list re-exports 10 of about 22 public `types.rs` items and no step module re-exports anything, so
  it is a partial convenience list rather than a contract.

## Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.89s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 12 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --all-targets --all-features` | 101 | lib `test result: ok. 3963 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 730.22s`; every integration-test binary ok; the run then hits the same **pre-existing** `benches/psp_writer_perf.rs:386` panic A1's review records |

Twelve tests where the review saw six. Each mutation run above restored the file from a copy
afterwards, and the `proptest-regressions` seed lines the mutant runs left behind were reverted with
`git checkout`.

## Follow-ups this run created

1. **The prune.** No method drops an allele, and that is deliberate: dropping allele *k* renumbers
   every id above it, so it cannot be `Vec::retain` handed out free — every `AlleleId` and
   `Genotype` minted before it goes stale. The method has to return the remapping, and it belongs
   with the step that needs one ([`calling_loop.md`](../../ng/impl_plan/calling_loop.md)). Now said
   in the type's doc comment rather than only in a report.
2. **`LocusKind` to `src/ng/types.rs`** — see above.
3. **`src/ng/mod.rs`'s module doc is a five-clause "and" chain** that each step extends. Worth a
   rewrite when a step has reason to touch it for its own sake.
