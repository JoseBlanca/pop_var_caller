# Handoff — the census lives inside the psp, Milestone B, step B4 and Checkpoint B

**Written 2026-09-09**, at commit `a628d7ad` on branch `census-vs-psp-perf`, in the worktree
`/Users/jose/devel/pop_var_caller-census-vs-psp`. **The working tree is clean and every step so far
is committed.** This exists so a fresh session can pick the plan up without re-deriving anything.

---

## What is left in this session's scope

**One step, B4, and then Checkpoint B — a hard pause.** The plan
([`psp_census_pair.md`](../../ng/impl_plan/psp_census_pair.md)) says:

> ☐ **B4 — a cohort judged whole.** Every psp of an opened cohort judged with B2, no early return;
> one verdict a sample, in the order given. A test with three stale psps in five asserts all three
> are named.
> *Depends:* B2, B3. *Source:* spec §4.1.

After committing it, stop and hand back. **Do not start Milestone C.**

---

## The three things B4 has to get right, already established

1. **No early return, and the reason.** Spec §4.1 exists because today's census cohort opener
   returns at the first census it cannot check, so a cohort with three stale ones reports one.
   B4's whole point is that one run tells the user the whole regeneration job.
2. **A psp that cannot be *read* is not a stale psp.** B2's judgement returns
   `Result<CensusVerdict, PspReadError>`. B2's review named the trap explicitly: a B4 written with
   `?` reproduces for read failures exactly the bug §4.1 exists to fix. **Hold a result a sample**
   — or fold the read failure into its own report row — rather than propagating the first.
   Its test should have an unreadable psp *and* several stale ones, not only stale ones.
3. **The cohort is `OpenPspCohort`, whose `psps` field is private** and whose readers need `&mut`
   for a trailer read. Either add an accessor or put the loop in `psp_caller.rs`; the judgement
   itself lives in `src/ng/run/census_freshness.rs` and is
   `what_the_footer_and_the_trailers_head_say_about_a_census(&mut PspReader)`.

---

## How the work is run — the loop, unchanged from A and B

Per step: **implement → review → apply fixes → commit**, one step at a time, never batched.

**The review arrangement, and why it departs from the `rust-code-review` skill.** The skill asks for
one agent per category, each in its own worktree; the build tree is 17 GB and the disk has 129 GB
free, so that does not fit. Instead: **one read-only agent over grouped categories, forbidden to
edit any file or run `cargo`, and asked to *name* mutations rather than run them.** The orchestrator
runs the ones it names, one at a time, reverting each from a backup and proving the restore with
`diff` before the next.

**That is where the value has been.** Across A and B, every step had at least one test that looked
sound and was green under the mutation it was written to catch:

- B1: the message printed the version this build reads from a field, and both tests set that field
  to this build's version — a message ignoring the field passed all four.
- B2: **every test passed with the whole trailer read instead of its front**, which is the property
  the step exists for. The fix was an instrument, not a test: `PspReader` now counts the trailer
  bytes it hands out.
- B3: nothing pinned which refusal wins when two settings differ at once.

---

## The gate, and it is a set and not a count

**The tree is not green and that is not yours to fix.** Baseline, measured on the untouched tree at
`b87c297d` and unchanged since:

| gate | baseline |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,682 lib tests**, 15 ignored, 20 of 21 targets green; the one failure is `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors of 5 kinds in 6 files** |
| `cargo check --all-targets --keep-going` | **4 examples** do not compile (1, 4, 8, 11 errors) |
| `cargo fmt --check` | **4 files** |

**At `a628d7ad` the suite is 6,700 lib tests** — the baseline's 6,682 plus B1's 4, B2's 10 and B3's
4 — with everything else identical to the table.

**Compare lists, never counts.** `tmp/b/gate.sh <label>` runs all four and prints clippy's error
*kinds* and their project-file locations, the targets `check` cannot compile, and the files `fmt`
would rewrite. That is deliberate: a check written from the baseline's own failures cannot see a
failure it has never seen before, and this session's B2 introduced exactly that — a new clippy kind
(`&PathBuf` where `&Path` does) and a new unformatted file, both invisible to a count.

**`cargo fmt` reformats the whole tree**, the four baseline files included. Run it, then
`git checkout -- examples/ng_call_cohort_end_to_end.rs examples/ng_call_from_psps_cost.rs
examples/ng_census_locus_spans.rs src/ng/psp/block.rs` before committing.

**⚠ And never run `git checkout --` on a file holding uncommitted work.** This session did, on
`psp_caller.rs` mid-step, and lost B3's whole implementation; it had to be retyped from the
conversation. If you need to undo a mutation, restore from the backup copy the mutation script made.

---

## Building

macOS with Apple's `container`: **every `cargo` goes through this worktree's own script**,
`/Users/jose/devel/pop_var_caller-census-vs-psp/scripts/dev.sh cargo …`. Run `python3` and `sed -i`
through it too. The shell's working directory does not persist between commands — use absolute
paths.

---

## What landed in Milestone B, in one line each

- **B1** (`66427170`) — `CensusVerdict`: fresh; no census; a census of another format, naming the
  version in the psp; a trailer that is not a census; a census built under another selection. With
  its `Display`, which is the last column of the refusal's one line a sample.
- **B2** (`89083b0b`) — `what_the_footer_and_the_trailers_head_say_about_a_census`: two of spec
  §4.2's three causes for one short read or none. Plus `PspReader::trailer_head`,
  `census_file::version_word_of`, and the trailer-byte counter that makes "it does not read the
  census" a measurement instead of a claim.
- **B3** (`a628d7ad`) — the cohort opener compares the catalog and the repeat-tract criteria as well
  as the ground, and refuses naming both samples and the field.

Reports: [`ng_psp_census_pair_milestone_b_2026-09-09.md`](ng_psp_census_pair_milestone_b_2026-09-09.md)
and [`../reviews/psp_census_pair_milestone_b_2026-09-09.md`](../reviews/psp_census_pair_milestone_b_2026-09-09.md).
Append B4 to both.

---

## Two things to raise at Checkpoint B, not to act on

1. **The spec has a factual error in the sentence that argues for B3.** `psp_census_pair.md` §6 says
   `SegmentationInputs::first_difference` "is called only from its own tests". It is called from
   `PspVariantCaller::open` and was before this branch, which means a calling run *already* refused a
   cohort typed two ways — naming one sample and the run rather than the pair. **The decision is
   unaffected** (the check still belongs in the shared opener, which is what gives the fit and
   `generate-census` a refusal they had no way to make), so the spec was left as it is and this is
   the owner's to rule on.
2. **A census that decodes badly past its version word has no verdict.** Magic and version right,
   sections truncated: no cheap read can tell it from a whole one, so it reaches the user as the
   census reader's own error rather than as *regenerate this*. Recorded in `CensusVerdict`'s own
   doc; nobody has asked for more.

---

## The commit convention

`<type>(ng): <step-id> — <title>`, ending with

```
Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
```

Flip the step's `☐` to `✅` in the plan in the same commit. Commit the code, the tests, both report
sections and the plan together.
