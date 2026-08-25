# ng candidate alleles — A2: review and the fixes applied

*Review report, 2026-08-24. Branch `ng-candidate-alleles`. Scope: the working-tree diff of step
A2 of [`candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md) — 391 added lines in one
file, **the output vocabulary, no logic**. Reviewed at `52646056` + the step patch, four agents
each in its own worktree. Fixes applied in the same pass; this file carries both.*

**Verdict: approve with changes, all applied.** 1 Blocker, 6 Major and 17 Minor as filed —
**four distinct defects above Minor** once the convergent filings are merged. The Blocker and one
of the Majors are the same shape as A1's: *a wrong number that produces plausible output instead
of a crash.* The other two are new: a hole an assertion could not see, and a method that would
have gone on silently working after the repeat-tract path adds a field to the type it clears.

---

## 1. Which categories ran

Four agents: **reliability**, **naming**, **idiomatic + errors**, and a **design-fidelity** pass
written for this step — does the code match arch §2.2–§2.4 field by field, and can B1, B2, C1,
C2, C3 and D1 each be built on what A2 declares. The last one is not a standard checklist; it
exists because A2's whole job is to be built on, and the cheapest place to find a missing
capability is before the fold is written.

| agent | Blocker | Major | Minor |
|---|---|---|---|
| reliability | 1 | 1 | 3 |
| design fidelity | — | 2 | 3 |
| idiomatic + errors | — | 3 | 6 |
| naming | — | — | 5 |

**Seventeen mutations by the reliability agent alone** (2 survived, 1 changed no behaviour), and
the design-fidelity agent proved its two findings by *writing* the code — a full `generic.rs`
implementing plan steps B1 through C3 against the A2 vocabulary, compiled, then discarded.

## 2. What was actually wrong

### 2.1 `alternative_allele_count` had no test at all — Blocker

Deleting its `- 1` passed all 13 tests. **No test in the module constructed a `LocusSelection`**,
so nothing exercised the one method on it.

Every answer moves by one: 1 → 2 on a two-allele table, and **0 → 1 at a locus that selected down
to the reference alone**, which spec §6.2 measures at 27.4% of tomato loci. The doc comment names
the consumer — the genotype prior divides its alternative concentration by this number — so an
answer one too high dilutes every real allele at every locus, and nothing panics.

**Fixed** by `the_alternative_count_excludes_the_reference`, asserting at one, two and six
alleles, plus a `selection_over` fixture that builds a real `LocusSelection`. Re-running the
mutation now fails that test and only that test.

### 2.2 The remapping's fourth hole: two merge alleles onto one candidate id — Major

Found independently by reliability and design fidelity, each by writing a passing test against
the unmodified code. `admitted` guarded the merge index three ways and the candidate id not at
all, so `admit(1, AlleleId(1))` followed by `admit(2, AlleleId(1))` was accepted: **both indices
in range, each written once, and `candidate_for(1) == candidate_for(2)`.**

**No bounds check can see this**, and it is worse than the off-by-one the plan isolates C1 for.
The evidence hand-off re-keys the merge's rows through this map, so two different sequences' reads
land on one candidate and the read likelihood scores two alleles as one — with an ordinary-looking
genotype coming out. `AlleleId(7)` where two alleles had been admitted was accepted too, leaving a
candidate table entry with no evidence behind it.

**Fixed** by carrying the admission count in the type and asserting `candidate.get() ==
num_admitted`. That is not an extra constraint invented here: every id comes from
`CandidateAlleles::admit`, which returns the table's previous length, so the next dense id is what
a correct caller always passes. It also turns `num_admitted()` from an O(table_len) scan into a
field read. Two tests, and removing the check fails both.

**The check that shipped is not the one either reviewer proposed, and the difference matters for
the path neither of them was reviewing.** One suggested bounding the id by the table length plus a
`debug_assert!` scan for aliasing, and warned — correctly — that a `candidate <= table_index`
bound would be **wrong for the repeat-tract path**, which admits by ladder rung rather than in
merge-table order (`arch/candidate_alleles_ssr.md` §3). The check applied here relates the id to
*how many alleles have been admitted so far* and never to `table_index`, so it is independent of
the order the caller visits the table in: it says only "you passed the id `admit` just handed
you". It holds for both paths, catches aliasing exactly rather than probabilistically, and costs
one comparison instead of a scan.

### 2.3 `reset_for` would not notice a new buffer — Major

`SelectionScratch::reset_for` cleared its two fields by name through `self`. The design-fidelity
agent added a third field, ran `cargo check --lib`, and got a clean build with only a `never read`
warning — one that vanishes the moment `ssr.rs` reads it. **And this is not hypothetical:**
`arch/candidate_alleles_ssr.md` §5 commits, as a decided item, to adding ladder buffers to this
exact struct.

So the repeat-tract path would have added a buffer, compiled, and carried the previous locus's
values into the next — which is the single failure `reset_for`'s own doc comment says it exists
to prevent.

**Fixed** by destructuring: `let Self { per_allele, ranked } = self;`. A field added later is now
`error[E0027]: pattern does not mention field ...` at exactly the line that must handle it.

### 2.4 `LocusSelection` had no door, so its invariant had no enforcer — Major

Three agents raised this. The type has four `pub` fields and no constructor, so C1, C3,
`select_ssr` and every test each write their own struct literal, and the invariant the doc comment
states — `unmatched` runs parallel to `CohortObservation::per_sample` — is enforced nowhere.

Plan step A2 asks for the invariants "in doc comments **and** asserted where they are built".
Deferring the *assertion* to C1/C3 is right, since nothing at A2 builds one. What was missing is
the place to put it.

**Fixed** by `LocusSelection::new`, taking `covering_samples` as a plain `usize` — the lighter of
the two shapes the reviewer offered, and the one that keeps this module from importing the
merge's assembled locus type. It asserts both halves: one leftover per covering sample, and
`remap.num_admitted() == alleles.len()`, so an admitted allele always has bases and a candidate
always has evidence. Two tests.

**The fields stay `pub` because arch §2.4 declares them so**, which means the constructor can be
bypassed. It still earns its place — it is what `select_generic` and `select_ssr` will call, so
the checks run on every value a run produces. **Making the fields private, so this is the only
door, is raised at Checkpoint A.**

### 2.5 `AlleleSummary::reached_the_bar` was a second copy of one fact — Minor

Two agents, and the reasoning is the same both times: the fold raises `samples_clearing_the_bar`
in the same branch that would have set the flag, so the two can only ever agree, and the only
thing a stored flag can do is disagree with the count. It was also a field the architecture does
not declare.

**Removed**, replaced by `cleared_the_bar()`, derived. The reference is not asked the question at
all — it is admitted before any sample's evidence is read (spec §6.1), so C1 seeds it structurally
rather than reading a flag that says it passed.

### 2.6 The judgment call the plan flagged, answered

A2 left the leftover's reads and mass off `AlleleSummary`, where arch §2.4 declares them. **The
design-fidelity agent was asked to check that reasoning and confirmed it, with a third argument
the author had not made:**

- `AlleleSummary` is per allele with no sample axis, while `unmatched` is per covering sample, so
  no per-allele total can produce the output;
- survival is not known during B1's pass, so the number could not be filled anyway;
- **and a cohort total is a sum in allele-major order where C3's oracle demands the per-sample
  rows' own sum** — so the bitwise check C3 is isolated for would fail *by construction* if the
  total were the source.

No reader wants it: the ranking uses the share, the samples clearing and `cohort_reads`, and arch
§6 already records that a future `q_sum` bar reserves nothing here. **Arch §2.4 and plan step B2
both owe an edit**, listed in §4 below.

### 2.7 Smaller things, all applied

- **`admitted` → `admit`**, matching `CandidateAlleles::admit`, the crate's verb for the same
  event. It read as a query on a mutating method.
- **`all_dropped` → `with_all_dropped`**, which cannot be misread as a predicate.
- **`num_alternatives` → `alternative_allele_count`**, which is how the crate already spells this
  exact quantity in `genotype_prior/seed_generic.rs` — the divisor the doc comment points at.
- **`SelectionScratch` is no longer `Clone`.** A cloned per-worker buffer pays the allocation the
  type exists to avoid, and carries the previous locus's fold — the one state that must never be
  read.
- **The allocation-contract sentence was not true of the code beneath it.** It said selection
  allocates nothing per locus beyond the surviving table; `AlleleRemap::with_all_dropped`
  allocates per locus and is not in the scratch. Now stated as what it is — the remapping is
  output, not a working buffer — with a pointer to arch §6's open item on the bitset.
- **The double-admission message now names both ids**, not just the index.
- **Three new panic tests**: the write-side range check, the dense-id check from both sides.
- `LocusSelection::verdict` was the only undocumented field in the diff.
- The `reset_for` doc named three wrong versions its grow test catches; the reviewer found a
  **fourth** — reset rows in place, grow but never shrink — that only the *shrink* test catches.
  Both tests now say what each one catches alone, because a reader who tries the three obvious
  mutants concludes the pair is redundant.

## 2a. One finding taken only in part

**`Clone` was proposed for removal from three types and was removed from one.** The reviewer
stripped it from `SelectionScratch`, `LocusSelection` and `AlleleRemap` and the suite still passed,
which shows none is *used* today.

`SelectionScratch` is removed and the reason is behavioural: a cloned per-worker buffer pays the
allocation the type exists to avoid and carries the previous locus's fold, which is the one state
that must never be read. The other two are kept. They are output values a consumer may reasonably
copy — the calling loop's own plan has not been written yet — and "no caller today" is a weak
argument for narrowing a type whose first consumers arrive two milestones from now. Recorded as a
disagreement rather than silently dropped; re-raise it when the loop's input edge exists and can
be asked whether it clones.

## 3. Checked and found sound

- **All four measured figures CHECKED-CORRECT** against the specs, including the one most likely
  to be misattributed: "98.6% of repeat tracts at 5× against 0.2% at 300×" matches
  `candidate_alleles_ssr.md`'s own table, and the attribution to a *cohort-summed* depth rule is
  what that document says at §6.
- **`alternative_allele_count` cannot underflow.** `CandidateAlleles` has one constructor, which
  pushes the reference, and one mutator, which only pushes, with private fields — so `len() >= 1`
  is structural rather than conventional. Verified in `calling/mod.rs`, not assumed.
- **Milestones B through D all build on these types with nothing missing.** C3 needs no third data
  structure. B2's bases come from `CohortObservation::alleles` via the table indices in `ranked`,
  and the borrow that could have failed — sorting `ranked` while the comparator reads `per_allele`
  — compiles under edition 2024. Nothing here blocks the repeat-tract path.
- **`AlleleRemap`'s half of the parallelism invariant is genuinely enforced**: the length is fixed
  at construction and both accessors assert in release.

## 4. Raised at Checkpoint A, not taken here

Carried forward from A1: the `support` naming, the cap's `u16`, `new_const`'s name, and the GATK
six-alternates question. A2 adds:

1. **Arch §2.4's `AlleleSummary` owes an edit** — the leftover's reads and mass are not fields of
   it and cannot be, for the three reasons in §2.6. **Plan step B2 asks for the same two fields**
   and needs the matching edit.
2. **`SelectionScratch::ranked`** is a bare participle on a struct field; the rename edits arch
   §2.4.
3. **The merge-table index crosses `AlleleRemap`'s public surface as a bare `usize`** while the
   other index space is an `AlleleId` — on the very type whose reason to exist is that the two are
   different numbers. Filed at Medium confidence with the counter-argument attached: the upstream
   producer `SupportedAllele::allele` is already a bare `usize`, so this is a crate-edge decision.
4. **Making `LocusSelection`'s fields private**, so `new` is the only door (§2.4).
5. **`spec/candidate_alleles.md` mis-cites `ng_step_interfaces.md` §3** for the `Admission` sketch,
   which is in that document's §2.

Two obligations recorded rather than raised, because they belong to steps that are still to come:
`Truncated { dropped: u16 }` over an uncapped `usize` merge table needs a saturating conversion at
C2; and spec §8 names three caller bugs that must assert, of which A2 lands one — the other two, a
non-finite `q_sum` and a sample with rows but no reads, belong to B1 and C3.

## 5. Validation after the fixes

All in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean; and `cargo clippy --lib`
  alone, which is where `dead_code` fires and `--tests` hides it — also clean.
- `cargo doc --lib --no-deps` — 23 unresolved intra-doc links, all pre-existing on `main`, none in
  this file.
- `cargo test --lib allele_candidates` — 20 passed (13 before the fixes, 7 added).
- **Mutation re-check on the fixed tree:** removing `- 1` from the alternative count fails 1 test;
  removing the dense-id assertion fails 2. Each applied from a file backup and reverted by
  restoring it, with the restore verified before the next run.

`cargo clippy --all-targets --all-features` remains red on `main` with 14 errors in five benches
and examples, none in `src/`, so this step is gated on `--lib`, `--lib --tests` and `cargo doc`.

## 6. One thing worth keeping from how this review ran

**The design-fidelity agent's most useful output was code, not prose.** Both of its findings were
proved by writing something — a third scratch field that compiled when it should not have, and a
complete `generic.rs` covering plan steps B1 through C3 that demonstrated the vocabulary is
sufficient. Asking "can the next four steps be built on this?" before writing them cost one agent
and closed a defect that would have surfaced in the repeat-tract path, on another branch, months
later.
