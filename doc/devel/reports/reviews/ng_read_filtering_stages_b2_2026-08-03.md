# Code Review: ng_read_filtering_stages_b2 (and Milestone B close-out)
**Date:** 2026-08-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B2 — `DecodedContainer::fill_raw_read` takes a whole `NoodlesRawAlignedRead` and
sets both its halves — plus an audit of Milestone B, this being its checkpoint
**Status:** Approve-with-changes (all applied)

---

## 1. Scope

- **What was reviewed:** the B2 working-tree diff (169 patch lines, two source files), exported
  as `tmp/review_2026-08-03_ng-read-filtering-stages-b2/b2.patch` and re-applied onto a detached
  `d5fb526`.
- **In-scope files:** `src/ng/read/input/aligned_reads_reader/{container,cram}.rs` and the impl
  report; for the close-out, the plan, spec, arch, `PROJECT_STATUS.md` and `8cf6f03..HEAD`.
- **Categories dispatched:** `reliability` (B2's whole subject is a read-group stamp, and a
  wrong one is silent) and `refactor_safety` + `extras` + the milestone close-out. Two agents
  for a two-file diff is the proportion the change warranted.

## 2. Verdict

**Approve-with-changes.** The change is correct and the diff matches its intent exactly. What
the review found is that **B2's justification for adding no tests did not hold**: the one
mutation the impl report ran was the one that could not survive, and two real gaps sat behind
it. Both are now closed.

## 3. Execution status

Both agents reproduced the orchestrator's numbers. Findings labeled "Needs verification": **0** —
the reliability agent ran seven mutations, each `grep`-confirmed present before its run and the
tree `git diff --stat`-confirmed back to the two-file patch after.

## 4. Findings

### Major

**M1: no CRAM fixture combines more than one container with more than one read group**
**Categories:** reliability · **Confidence:** High (mutation-verified, missing test written and run)

The gap is a hole between two fixtures. `indexed_cram_declaring` is the only multi-read-group
CRAM and holds three records — one container, since noodles writes 10,240 per container.
`multi_container_cram` is the only multi-container CRAM, declares a single `@RG`, and its cursor
tests open it as `ReadGroupResolution::Sole`, an arm that short-circuits per-record resolution
entirely. **So no test reached the per-record read-group arm past a container boundary at all.**

Overriding every read's group with a hard-coded `ReadGroupId(0)` from the second container
onwards left **`1542 passed; 0 failed`**.

**Failure scenario:** the one `cram.rs`'s own doc names as the worst on this arm — a library's
reads silently attributed to another read group, corrupting the per-library error model with no
error raised. The tomato acceptance dumps are CRAM and multi-container, but they compare loci
and observations, not read-group attribution, so byte-identity does not cover it either.

**Fix (applied):** `multi_container_cram_two_read_groups` and
`a_read_past_the_first_container_carries_its_own_read_group`. Verified to fail under the
mutation at record 10241 and pass on the shipped code.

**M2: the read group's *value* rested on one test, and the second test credited for it threw the value away**
**Categories:** reliability · **Confidence:** High (mutation-verified)

`read_group = None` fires both tests the impl report named — but only because
`RawAlignedRead::decode` refuses an unstamped buffer, so it proves *presence*, not correctness.
The two wrong-but-valid mutants (`Some(ReadGroupId(0))`, `Some(self.index[0].owner)`) each died
to **one** test. `a_shared_cram_serves_each_open_only_its_own_reads` collected
`(qname, read_group)` pairs and discarded the group with `.map(|(qname, _)| qname)`.

**Fix (applied):** assert the pairs. One line; both mutants now die twice.

### Minor

**Mi1** — `fill_raw_read`'s doc promises it sets "both halves" while the body reached the two
fields by name, so a third field would compile silently and leave the doc vouching for something
false. **Applied:** exhaustive destructure.
**Mi2** — `container.rs` has **no test module at all**, so the function B2 changed has no direct
test; uncovered classes: an unnamed record, empty spans, and the clear-and-refill claim its own
doc makes. **Deferred**, recorded — pre-existing coverage of pre-existing behaviour, and a new
test module is its own work. Worth taking before C2.
**Mi3** — `out.data_mut().clear()` has a doc-comment invariant naming its own silent failure
mode, and deleting it leaves the suite green. Currently unreachable (every caller passes a fresh
buffer), so latent rather than live. **Deferred with Mi2.**
**Mi4** — the impl report's §4 conclusion generalised from the single weakest mutant.
**Applied:** §4 now carries the seven-row log.

### From the close-out

**Mi5** — **the deviation recorded at B1 was recorded inaccurately, and it is one the owner is
ruling on.** Both `PROJECT_STATUS.md` and the B1 reports said the new error variant went
"against spec §1's *adds no new error*". **Spec §1 contains no such sentence**: it says "change
the meaning of any error", which adding a variant does not do, and the only "No new error type"
statement is **arch §4**, scoped to `ReadFilterError` — a different enum. "The spec forbade it"
and "the spec is silent" call for different decisions. **Applied:** corrected in all three
places.
**Mi6** — the spec still said **"no code yet"** after five code commits; its §1 bullet described
a probe loop B1 deleted and cited a line that is now an unrelated function. **Applied.**
**Mi7** — §11's reuse map — **the table Milestone C executes against** — had drifted by up to
120 lines (`RecordSource` 366 → 338, `ReadFilter::next` 895 → 776). **Applied.**
**Mi8** — `with_validated_contigs`' shipped doc still asserted the comparison "proves strictly
more", the exact claim B1's review overturned and which produced B1's Blocker. **Applied.**

### Verified with no finding

Ordering inside `fill_raw_read` (no `?`, no early return, no observable half-filled state); the
BAM and in-memory arms untouched and still clearing, with the module contract still accurate;
and — **for the first time in this plan** — every *factual* claim in the impl report checked
out, including "exactly one caller" for the deleted accessor, verified against the base commit,
with no test whose property could have been lost.

## 5. Milestone B close-out

**Milestone B delivers what the plan asked.** Checkpoint B's three requirements are met: the
dumps and the anchor unchanged, and the walk probe's `seconds` measured before and after — with
B1 correctly recording the measurement rather than adopting the estimate.

**Nothing in `8cf6f03..HEAD` is unaccounted for.** Eight commits, all mapping to named units:
A1–A3, the three Checkpoint-A follow-ups, B1, B2.

**Deviation (b) — spec §9 Q2's cost arithmetic — is recorded accurately**, and a second agent
confirmed it independently in `ref_seq.rs`: the 52 µs is per-`WindowedRefSeq`-construction
("34 µs **per accessor** to clone"), while the per-contig-open figure is a different number
entirely. Spec §9 Q2 multiplies the per-accessor figure by 2,580.

**Deviation (a) was recorded inaccurately** — see Mi5, now corrected.

## 6. Missing tests

M1's fixture and test are applied. Mi2's three (`fill_raw_read_sets_no_name_for_a_record_that_
has_none`, `fill_raw_read_leaves_no_tail_of_a_longer_previous_read`,
`fill_raw_read_clears_auxiliary_tags_left_by_a_previous_fill`) are specified in
`tmp/review_2026-08-03_ng-read-filtering-stages-b2/reliability.md` and deferred with their
reason.

## 7. What's good

- **The reviewer proved the container-boundary gap was not an equivalent mutant** by writing the
  missing fixture and test, running them green on the shipped code, and red under the mutation —
  rather than asserting the gap from the mutation alone.
- **The close-out checked a claim the *reviewer's own brief* asserted** — that B1 deviated from
  spec §1 — and found it false. That is the review catching the orchestrator, which is what the
  fan-out is for.
- **Seven mutations on a 169-line diff.** The proportion was right: this diff is small and its
  subject is silent-failure-prone.

## 8. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --lib
cargo test --lib ng::read::input::open_bam::tests::a_read_past_the_first_container_carries_its_own_read_group
```
Plus the four acceptance dumps against the `8cf6f03` baseline — the two tomato ones are the
CRAM path this step touches.
