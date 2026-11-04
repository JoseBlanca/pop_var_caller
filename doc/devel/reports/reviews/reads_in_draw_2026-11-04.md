# Code Review: reads_in_draw
**Date:** 2026-11-04
**Reviewer:** rust-code-review skill (single reviewer, no fan-out)
**Scope:** reads in the neighbour-evidence pass's haplotype draw, branch `reads-in-draw`, diff `480c6341..1e5c94ab`
**Status:** Approve-with-changes; every finding disposed 2026-11-04 (§11): M1 fixed by holding records out of the draw, Mi1 accepted by the owner

---

## 1. Scope

- **What was reviewed:** the commit `1e5c94ab`, which builds
  [`neighbour_evidence_reads_in_draw.md`](../../ng/spec/neighbour_evidence_reads_in_draw.md): along
  each path a round draws, a sample's read fragments orient its consecutive heterozygous records
  through a two-state chain, and the joint's heterozygous mass moves onto the drawn order.
- **Reviewed against:** commit `1e5c94ab`.
- **In-scope files:** `src/neighbour_evidence/{model,refine,random,stage}.rs`,
  `src/neighbour_evidence/model/tests.rs`, `src/neighbour_evidence/stage/tests.rs`, the spec, and the
  implementation report `doc/devel/reports/implementations/reads_in_draw_2026-11-04.md`. Read as
  direct callees and consumers: `src/read_phase/link.rs` (`score_read_link`, `Heterozygote`),
  `src/read_phase/fragments.rs`, `src/neighbour_evidence/phase.rs` (`phase_heterozygotes`,
  `link_to_next`), `src/cli/calling_run.rs` (the read check's log line).
- **Deliberately out of scope:** the link test itself (reviewed 2026-11-03), `vendor/`.
- **How it was reviewed.** At the owner's standing request no git worktrees were made and no
  sub-agents dispatched: one reviewer read the code in place and applied every checklist that the
  skill's triage selects (reliability, errors, naming, defaults, idiomatic, refactor_safety,
  module_structure, float_portability, smells, extras). Mutations were run on a copy of the crate
  under `tmp/probe_sign/`, built into `tmp/probe_target/`; no tracked file was changed.
  `PROJECT_STATUS.md` has no block for this branch, so it was not updated.

## 2. Verdict

**Approve-with-changes.** No blocker. The chain's arithmetic is right: its sign agrees with
`score_read_link`'s *cis*, it cannot overflow at the read cap, a group of one draws nothing, and
output without linked fragments is bit-identical to before. Two major findings. The first is a
consequence of the design rather than a coding slip: the run log's *read check*, which the read phase
spec presents as an estimate of wrong population links "that does not come from the model", now
measures links the reads themselves made, and on the dense 64-founder cohort it fell from 125
contradicted of 3,339 to 1 of 8,034 while the true wrong-link rate is still 2.0%. The second is that
no test pins the sign of `L` where `refine.rs` hands it to the chain.

## 3. Execution status

Run at `1e5c94ab`, release profile, native cargo on macOS arm64:

- `cargo fmt --check` — exit 0, no diff.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 101, five errors, all in files
  this commit does not touch: `src/parameter_estimation/joint/ssr_fit/tract_exclusion.rs:80`, `:85`
  (`doc list item without indentation`), `:316` (negated comparison on a partially ordered type),
  `:590` (`manual Range::contains`), and `src/run/census_fit.rs:629` (`too many arguments (8/7)`).
  None in the scope.
- `cargo test --release --lib neighbour_evidence` — `test result: ok. 106 passed; 0 failed`.
- Not run: `cargo doc`, `cargo audit`, the integration tests.
- Mutation probes, on the copy under `tmp/probe_sign/`:
  1. `L`'s sign flipped where `refine.rs` hands it to the chain — `cargo test --release --lib`:
     `test result: FAILED. 5185 passed; 2 failed`. The two failures are the pre-existing debug-only
     `should_panic` tests (`an_offset_past_the_reference_panics`,
     `lgamma_panics_on_non_positive_in_debug`), which fail in every release build. **The mutation
     survives the whole suite** (M2).
  2. `orient_by_reads` made a no-op inside `draw_haplotypes` — `cargo test --release --lib
     neighbour_evidence`: `105 passed; 1 failed`, the one failure
     `the_pass_with_reads_writes_the_same_bytes_at_any_thread_count`.
     `reads_change_the_haplotypes_and_not_the_walk` passes (Mi2).

## 4. Open questions and assumptions

1. **What should the read check report now?** (M1.) It can be relabelled as no longer independent of
   the model, restricted to links the chain did not see, or computed from paths drawn without reads.
   Recommended: relabel now (one log line and one spec paragraph), and decide the rest when the
   tomato read check of spec §6 point 6 is run — that comparison is the one this finding confounds.
2. **Should the draw have the link test's `GQ` guard?** (Mi1.) The guard exists because a false
   heterozygote gives a confident link of random sign; the draw's only guard is a heterozygous mass
   of one half on the path. The measurements show fewer wrong links overall, so the recommendation is
   to keep the code and record the decision in spec §5.

## 5. Top 3 priorities

1. **M1 — the read check's meaning.** The log line still says the check "does not come from the
   model"; on the cohort this change is for, it now reports about 1 contradicted link in 8,000.
2. **M2 — a test of the sign of `L` through `refine.rs`.** Flipping it there is a one-character change
   that the whole lib suite passes; a flipped sign would put every read-linked pair on the wrong
   copies in every round.
3. **Mi2 — `reads_change_the_haplotypes_and_not_the_walk` cannot fail on a no-op draw.**

## 6. Findings

### Blocker

None.

### Major

**M1: [src/cli/calling_run.rs](../../../../src/cli/calling_run.rs#L516) and
[src/neighbour_evidence/phase.rs](../../../../src/neighbour_evidence/phase.rs#L483) — the read check
now measures links the reads themselves made, and still says it does not come from the model.**
*Confidence: High.* The read check (read phase spec §2.8) counts the population links that reads
reach at `|L|` ≥ ln 199 with `GQ` ≥ 30 at both ends, and how many of those the reads contradict. The
log line calls it "a check of the population links that does not come from the model". With reads in
the draw, two such records are usually consecutive in a sample's chain, so every one of the 36 paths
is drawn with a factor of at least `e^{±2.6}` towards the reads' relation; the consensus follows the
reads, and the reads then "confirm" it. Measured from the implementation report's own run logs
(`tmp/reads_in_draw/*/calls_{off,0.5}.log`, 50 samples at 3×, SNPs 1 in 100 bases):

| founders | contradicted, draw off | contradicted, draw on | true wrong links, draw on |
|---|---|---|---|
| 64 | 125 of 3,339 (3.7%) | 1 of 8,034 | 916 of 45,443 (2.0%) |
| 16 | 7 of 4,933 | 0 of 8,627 | 30 of 57,168 |

**Failure scenario:** someone reads the tomato run's log after this merge and sees about 1
contradicted link in 8,000, and takes the population phase to be that good; spec §6 point 6 plans
exactly this comparison "with and without the change". **Why it matters:** the check is the only
run-level estimate of phase error on real data, where there is no truth. **Fix:** at least, change
the log line and spec §2.8 to say the check is no longer independent of the model wherever reads
enter the draw, and note in the reads-in-draw spec §2.5 that the read check *does* change. To keep an
independent figure, count it only at links whose two records were never in one read-linked group of
that sample (the chain knows which pairs it linked; a per-sample bitset over sites, set in
`orient_by_reads`, would carry it to `link_to_next`). *Category: extras (diff matches stated intent),
reliability.*

**M2: [src/neighbour_evidence/refine.rs](../../../../src/neighbour_evidence/refine.rs#L929) — no test
pins the sign of `L` where it reaches the chain.** *Confidence: High (mutation run).* The model tests
hand the chain synthetic links through their own closures; the stage tests check only that output
with reads differs from output without, which a flipped sign also satisfies; the `GQ` test runs with
the draw off. Replacing `then_some(found.log_odds)` with `then_some(-found.log_odds)` passes all
5,185 lib tests that passed before. **Failure scenario:** a later refactor that changes which allele
`orientable_pair` hands first, or wraps `L` for another caller with the other convention, writes
every read-linked pair in the opposite relation to its reads, and nothing fails. **Fix:** a
`refine.rs` unit test with two heterozygous records whose templates are alike, the sample's fragments
showing allele 0 at both (cis, `L` ≈ +7.8 for two fragments), and the drawn haplotypes asserted to put
allele 0 on one copy at both records in nearly every round; and the mirror case with trans. Run it
against the mutation above. *Category: reliability.*

### Minor

**Mi1: [src/neighbour_evidence/model.rs](../../../../src/neighbour_evidence/model.rs#L709) — the draw
orients records the link test would not, and the spec does not say so.** *Confidence: High.* The link
test makes a read link only with `GQ` ≥ 30 at both ends, because a false heterozygote gives a
confident link of random sign (`link.rs`, `GQ_MIN`). The chain's only gate is a heterozygous mass of
one half on the path. The report found the consequence (a record phased at `GQ` 14) and changed
`a_column_reads_newly_phase_has_gq_30_or_more` to run with the draw off, calling it something "the
population is allowed to do". **Failure scenario:** at 3×, a true homozygote whose path gives it
heterozygous mass 0.6 and one fragment to a neighbouring heterozygote gets a 98-in-100 orientation
from that fragment, which other samples then copy. The measured cohorts show fewer wrong links
overall, so this is a documentation gap, not a defect shown. **Fix:** add it to spec §5's traps and
record the owner's ruling there; no code change unless the tomato or HG002 runs of spec §6 point 6
show it. *Category: defaults.*

**Mi2: [src/neighbour_evidence/model/tests.rs](../../../../src/neighbour_evidence/model/tests.rs#L665)
— `reads_change_the_haplotypes_and_not_the_walk` cannot fail when the draw does nothing.**
*Confidence: High (mutation run).* Its doc says the haplotypes "change, and they do". Every site has
the same templates and likelihoods, so without reads every site's copy-one margin is the same value,
and "every site has the same `one[0] > 0.5`" holds whether or not the chain ran. With
`orient_by_reads` made a no-op the test passes. **Fix:** assert that the with-reads haplotypes
differ from the without-reads ones, and that at every heterozygous site `one[0]` is above the
without-reads value by most of the heterozygous mass:
```rust
assert_ne!(with, without, "the reads change the haplotypes");
assert!(with.iter().zip(&without).all(|([one, _], [before, _])| (one[0] - before[0]).abs() > 0.2));
```
*Category: reliability.*

**Mi3: [doc/devel/reports/implementations/reads_in_draw_2026-11-04.md](../implementations/reads_in_draw_2026-11-04.md)
— the report leaves out how often the link test's read rows still fire, which spec §6 point 3 asked
for.** *Confidence: High.* The numbers are in its own logs: on the dense 64-founder cohort, links made
by reads fell from 1,997 to 93, population links broken by reads from 480 to 12, and heterozygotes
attached by reads from 837 to 145. That is the evidence that the draw has taken over most of the link
test's work, and it bears on M1. **Fix:** add the three counts for both dense cohorts to the report.
*Category: extras.*

### Nits

- `model.rs:825` — the guard `if index + 1 < group.len()` repeats what `drawn.get(index + 1)` already
  decides; `drawn` has the group's length. Drop the guard, or index `links` and `drawn` directly.
- `refine.rs:82` — `HET_MASS_MIN`'s doc says one half gave "the most heterozygotes phased at both 16
  and 64 founders"; the report's table gives the phased count at 64 founders only.
- `stage.rs:74`, `refine.rs:265` — a public field and a public builder that exist only for one test.
  `pub(crate)` would keep them out of the library's API.
- Spec §1.5 and §2.2 still point at `copies_along` (`model.rs:677`), now a `#[cfg(test)]` wrapper; the
  joint is built in `joint_along` (`model.rs:893`).
- `orient_by_reads` allocates `orientable` and `links` for every path, and the link closure runs two
  binary searches per site per call, each site looked up twice. The report measures the whole step at
  about 3%, so this is for a later performance pass, not now.

## 7. Unconfirmed

- [src/neighbour_evidence/stage/tests.rs](../../../../src/neighbour_evidence/stage/tests.rs#L493) —
  `fragments_a_line_carries_reach_the_pass_and_break_a_link_they_contradict` now runs with reads in
  the draw. Its doc says the link test breaks the link. If the two records were consecutive in the
  sample's chain, the draw would put them in trans and they would stay in one phase set, so the
  assertion would fail; since it passes, the draw probably never links them (an orientable record with
  no fragments likely lies between them). That was not checked. If it holds, the test is fine but
  depends on a fixture detail; running it with `join_stage_drawing(.., false)` would make its premise
  explicit.

## 8. Out of scope observations

- Clippy at `-D warnings` fails on five lints in `tract_exclusion.rs` and `census_fit.rs` (section 3),
  outside this commit.

## 9. What was checked and holds

- **Sign convention.** `orient_group`'s *straight* is copy one carrying the lower allele of each pair
  (`joint[a·A + b]`, `a < b`); `same` gives `e^{+L/2}` to equal orientations
  ([model.rs:806](../../../../src/neighbour_evidence/model.rs#L806)). `score_read_link` counts a
  fragment as cis when "`xᵢ` is the lower allele exactly when `xⱼ` is"
  ([link.rs:211](../../../../src/read_phase/link.rs#L211)), so positive `L` favours equal
  orientations, as it should. The pair handed over is `orientable_pair`'s `(a, b)` with `a < b`, and
  `Heterozygote::of` orders it again, so a `1/2` heterozygote reads lower allele 1, not 0.
- **No overflow or underflow.** One fragment weighs at most `ln(0.9802 / 0.0198)` = 3.90 nats (the
  error floor 0.01), and a sample keeps at most 32 read alleles a record, so `|L|` ≤ 124.8 and
  `e^{|L|/2}` ≤ 1.3 × 10²⁷. `exp` overflows only past `|L|` = 1,419, about 364 fragments at the
  largest weight. The forward pass normalises at every site, so nothing accumulates along a group,
  and its total is positive because an orientable site's own weights sum to at least one half and
  the link factor is at least `e^{−62.4}`.
- **Groups of one.** `orient_by_reads` calls `orient_group` only when `index − start > 1`; a lone
  orientable site keeps its joint and reads no uniform. Traced for a linked pair followed by a lone
  site, and tested for all-lone sites by `an_unlinked_site_keeps_its_joint`.
- **Bit identity without links.** The old `copies_along` summed `joint[i] / total` into each margin;
  `joint_along` now stores `joint[i] / total` and `copies_of` sums the stored values in the same
  order — the same IEEE operations. With no link, `orient_by_reads` writes nothing and reads no
  uniform. A window without fragments skips the step entirely
  ([refine.rs:937](../../../../src/neighbour_evidence/refine.rs#L937)).
- **Determinism and stream independence.** `StreamPurpose::Orientation` is a fifth purpose code, keyed
  by window, sample and round; the own path's orientations are drawn before the extra paths' in
  `walk_with_draw`, and neither the `Path` nor the `ExtraPaths` stream is read by the new step.
  `the_pass_with_reads_writes_the_same_bytes_at_any_thread_count` passes at one and eight threads.
- **Hidden records.** Filtered before `orientable_pair`, by binary search in the sample's sorted
  hidden list; `hidden_and_homozygous_sites_stay_out_of_a_group` would fail if they joined, since the
  hidden site's joint would lose one of its two orders.
- **Float portability.** The only transcendental call is `crate::float::exp`; sums run over two
  elements in a fixed order.

## 10. Commands to re-verify

- `cargo test --release --lib neighbour_evidence`.
- The M2 mutation: `then_some(found.log_odds)` → `then_some(-found.log_odds)` at
  `src/neighbour_evidence/refine.rs:929`; the test of section 6, M2, must fail.
- The Mi2 mutation: make `orient_by_reads` return at once; `reads_change_the_haplotypes_and_not_the_walk`
  must fail once fixed.
- The read check counts: `grep -o "read check.*attached by reads" tmp/reads_in_draw/f64_F0_n50_d3_s0.01/calls_{off,0.5}.log`.

## 11. Dispositions (2026-11-04)

| finding | disposition |
|---|---|
| M1, the read check is no longer independent | **Fixed, on the owner's ruling.** One in twenty of each sample's core records with fragments is held out of the draw, and the check counts only links with an end held out (spec §2.6). On the 64-founder dense cohort it reads 9 of 474 (1.9%) against a true 2.1%; the implementation report has the cost. |
| M2, no test pins the sign of `L` | **Fixed.** The closure's body is now `refine::draw_link`; `the_draws_link_is_positive_for_cis` checks cis, trans and a `1/2` heterozygote, and fails under the §10 mutation. |
| Mi1, no `GQ` guard in the draw | **Accepted by the owner.** No guard; the reasons are in spec §5. |
| Mi2, a test that passes on a draw doing nothing | **Fixed.** It now asserts the haplotypes with reads differ from those without. |
| Mi3, the link test's read rows | **Fixed.** A table in the implementation report. |
| Nit, redundant guard in `orient_group` | Fixed. |
| Nit, `HET_MASS_MIN` doc claims a count the report lacked | Fixed: the report's threshold table now gives the 16-founder phased counts. |
| Nit, test-only switch is public | Fixed: `pub(crate)`. |
| Nit, stale `copies_along` pointers in the spec | Fixed. |
| Nit, per-path allocations in `orient_by_reads` | Not changed: the measured cost is 3% of a run at the densest cohort, and the vectors are a few entries a group. |
