# ng candidate alleles — B1: review and the fixes applied

*2026-08-24. Step B1 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md).
Reviewed at `823c7b77` plus the step's working-tree diff, six agents in six isolated worktrees.
Implementation report: [`../implementations/ng_candidate_alleles_b1_2026-08-24.md`](../implementations/ng_candidate_alleles_b1_2026-08-24.md).*

---

## 1. Which categories ran

`reliability`, `errors`, `naming`, `idiomatic`, `smells`, and a sixth that is not a checklist:
**design fidelity**, asked whether the steps that follow could be built on what this one declared,
and told to answer by writing them. Each agent got its own worktree, detached at `823c7b77` with
the step's patch applied, and was told to mutation-test rather than review by reading.

`defaults`, `module_structure`, `unsafe_concurrency` and `tooling` were not dispatched: the step
adds two private functions to one existing file, has no public API, no configuration, no
concurrency primitive and no `Cargo.toml` change.

**The fan-out earned its cost in one specific way.** The two findings that mattered most were
*not* wrong code — the fold computed the right answer on every input either agent could build.
They were **tests that could not fail**, and only mutation found them. This is the third milestone
running in which that is the headline, which is now a pattern rather than an accident.

## 2. What was actually wrong

### 2.1 Nothing made the *share* half of the rule decide — Blocker

The admission rule is `max(2 reads, share × the sample's compared reads)`, and its two halves bind
at opposite ends of the depth range. **In all seven of the step's fixtures the floor decided**, so
the share term was never load-bearing. Two mutants proved it: replacing
`min_allele_support.reached_by(pooled, compared)` with `pooled >= floor.get()` — the share deleted
outright — left the suite green, and so did asking the rule against the allele's own reads instead
of the sample's.

The step's own test helper documents the hazard and then walks into it. Its doc comment says the
shipped 5-in-100 share "is inert below 41 compared reads and every fixture below is a handful of
reads", and the author's response was to raise the *share* in the fixtures rather than the
*depth* — which moves the number without making the term decide, because `ceil(0.5 × 9) = 5` is
still reachable by a floor of 2.

**Why it matters:** the share is the only half of the rule doing any work at high coverage. On the
GIAB trio at 300× a 5-in-100 share keeps 2,308 of the merge's 15,474 alternatives where a
count-only bar keeps 10,793 (spec §3.3). A fold that shipped applying the floor alone would admit
sequencing error as a candidate allele at depth, dilute the genotype prior's alternative
concentration across it, and crash nothing.

**Fixed** by `the_share_refuses_what_the_floor_would_admit`: one sample with 100 compared reads and
3 on the alternative, against 2 reads or 5 in 100 — the floor would admit it, the share asks 5, so
it is refused; at 5 reads it is admitted, so the fixture is a boundary and not a blanket refusal at
depth. It also pins the denominator: asked against the allele's own 3 reads, `ceil(0.05 × 3) = 1`
and the allele would clear. Both mutants now fail it.

### 2.2 `samples_clearing_the_bar` was asserted as 1 everywhere and never above — Blocker

Every assertion on that field expected exactly `1`, so replacing `+= 1` with `= 1` left the suite
green. No fixture had two samples clearing the rule on one allele — the only two-sample fixture was
built so the second sample deliberately *fails*.

**Why it matters:** that count is the cap's first tie-break, and spec §4.1 makes it the *deciding*
key at low coverage — "at 3 reads every admitted allele has a share near 0.67, the first key ties,
and how many samples showed it decides, which is the only signal there is at that depth". A count
stuck at one makes the first two ranking keys tie together and drops the ranking through to the
cohort read total, which is production's ranking and the one spec §4.1 exists not to be: at a
thousand samples it truncates the private alleles first.

**Fixed** by `every_sample_that_cleared_the_rule_is_counted`.

### 2.3 A sample with rows and no reads was skipped where spec §8 says assert — Major

The fold guarded a zero denominator with `continue`. Spec §8 names three caller bugs this step
must answer with assertions held in release, and one of them verbatim is "a compared-read count of
zero on a sample that has rows". Two agents found it independently, and the errors agent showed the
cost is larger than the wording: **because the skip fired before the row loop, such a sample was
never checked by the two conditions the fold does assert.** It demonstrated that with a probe — a
zero-read sample carrying a row naming allele 9 of a two-allele table passed without panicking.

The guard also conflated two states. A sample with **no rows at all** is legitimate and reachable:
`per_sample` holds the samples that *covered* the locus, so a sample whose reads all stopped inside
it has partials and nothing else.

**Fixed** by separating them: an empty `supported` is stepped over with a comment saying why, and
rows that carry no reads now assert. Two tests pin the pair, and disabling the assertion fails one
of them.

### 2.4 The doc comment said the reference is not asked the rule; the code asks it — Major

The fold's documentation read "the reference is folded like any other allele and is not asked to
pass anything", while the code ran `reached_by` for allele 0 like any other. Nothing pinned the
reference's row beyond its read total, so guarding that branch with `allele != 0` also left the
suite green.

Inert today, because C1 seeds the reference structurally. Not inert as documentation: a C1 author
who loops `cleared_the_bar()` over the whole table, trusting that sentence, either double-seeds the
reference or drops it — and no test in the module would have noticed either way.

**Fixed** in both directions: the sentence now says the rule *is* asked of the reference and that
nothing downstream reads the answer, and `the_reference_row_is_folded_like_every_other_allele`
asserts all three of its fields.

### 2.5 The share was maximised in one direction only — Major

Every fixture with two samples on one allele put the larger share first. A first-wins assignment —
`if best == 0.0 { best = share }` — therefore survived the whole suite, while making the cap's
first ranking key depend on the order the samples happen to be walked in. Spec §8 requires the
output to be byte-identical at any worker count.

**Fixed** by `the_largest_share_wins_when_it_arrives_last`, the mirror of the existing fixture. The
assignment is now `.max(share)`, which has no first/last asymmetry to get wrong.

### 2.6 A sample listed twice was folded twice, in silence — Minor

Rows out of allele order panicked; a repeated entry in `per_sample` did not. It is the same failure
one level up — one sample's evidence lifts its allele's cohort total and clears the rule as two
samples — and the fold now asserts ascending sample order beside the other two.

### 2.7 Smaller things, all applied

- **The one measured number in the new prose was wrong.** The pooling argument said a per-row rule
  would be stricter on "the samples that carry two libraries — 157 of 1,707". `read_groups.md` §1
  says 157 carry *more than one*: 133 carry two, 20 carry three, and four carry 7, 16, 16 and 42.
  The population is also samples with more than one *read group*, which is the axis the rows are
  keyed on. Corrected, with the breakdown.
- **A comment cited the architecture section that says the opposite.** The zero-denominator guard
  justified itself by arch §2.5, whose next sentence is "a `NaN` share is a caller bug and
  **asserts in the fold**". Moot now that the fold asserts.
- **Neither panic message named the locus.** With `panic = "abort"` the message is the whole
  post-mortem, so a caller bug at one locus among millions could not be re-run. All three
  assertions now carry `observation.region`.
- **`#[allow(dead_code)]` on `compared_reads_of` was redundant** — rustc treats an item carrying an
  `allow` as a live root, so the one on `summarise_alleles` already covered its callee, and the
  attribute's stated reason claimed a diagnostic that does not exist. Removed, and the surviving
  allow now says it covers both. *Verified the hard way during the fix: dropping **both** turns
  `cargo clippy --lib` red with two errors, which is how the redundancy was confirmed rather than
  assumed.*
- **The `row` test helper took `read_group` and `num_reads` as adjacent bare `u32`s**, so a
  transposition would still have been well-formed. Split into `row(allele, num_reads, q_sum)` for
  the one-library case and `row_from_group(allele, ReadGroupId, num_reads, q_sum)` for the two-lane
  fixtures, which is also the shape that says which fixtures are *about* read groups.
- **The parameter was named `bar`** where the config field it will be passed from is
  `min_allele_support`. Renamed, and the test helper with it.
- Four coverage gaps the reliability agent named and the fix took: a table of one allele, a locus
  no sample covers, `compared_reads_of` on a sample with no rows, and that sample folded before a
  carrier.

## 3. Checked and found sound

- **`chunk_by` as the pooling tool.** Three agents judged it against arch §3.1's nomination of
  `SampleSupport::pooled_support_for` and all three preferred what is written: the method answers
  one allele, so reaching a sample's distinct alleles through it costs a scan of the whole table
  per sample, and it rebuilds all six of the merge's quality moments where the rule reads one. The
  reason is now in the doc comment, which is what was missing.
- **`get_mut(..).unwrap_or_else(|| panic!(..))`** is the idiomatic form when the message needs
  formatting; `expect(&format!(..))` is what `clippy::expect_fun_call` rejects.
- **The fold as a free function rather than a method on `SelectionScratch`** — the scratch's methods
  are about buffer lifetime, the fold's subject is the locus, and B2 and C1 are free functions over
  the same buffer.
- **The saturating arithmetic and the `f64` share**: every bound unreachable, no `as` cast, no
  precision loss.
- **`reset_for`'s destructuring did what A2 built it to do.** The design-fidelity agent added a
  field to `SelectionScratch` while prototyping C3 and the build failed at that method with
  `error[E0027]`, rather than silently carrying one locus's values into the next.

## 4. Raised, not applied — decisions for the owner

Four things the review surfaced that this step must not decide on its own. They are put to the
owner in the session's Checkpoint B message; none blocks the commit, because the code as committed
follows the documents as written.

1. **The cap's ranking is fed by samples the rule refused.** The share is maximised over *every*
   sample, so a sample with one compared read contributes a share of 1.0 to any allele it shows
   while clearing nothing. The design-fidelity agent built the case: at a cap of two alleles, an
   allele earned at 20 of one sample's 100 reads is cut in favour of one earned at 10 of another's
   100, because a third sample showed the second allele once. Spec §4.1's wording licenses the
   code — "maximised over samples" — and its own prediction about behaviour at 3 reads only holds
   if the maximum is taken over the samples that cleared the rule. **Recommendation: restrict the
   maximum to samples that cleared the rule, and correct §4.1's sentence.** It is one line of the
   fold and one of the spec.
2. **Spec §1.3 defines "compared reads" twice, and the two differ.** Once as
   `SampleLocusObservations::reads_compared_with_reference`, once as "the sum of a sample's rows in
   the merge's table". They differ by `reads_removed_as_evidence`, and the plan's B1 oracle cites
   the first while the code implements the second. **The code is right** — arch §3.1 says the same,
   and a read the merge withheld from the table cannot be in a numerator either — but the spec
   sentence should lose its first half.
3. **The measurement harness and the shipped fold now disagree by design.**
   `examples/ng_candidate_selection_probe.rs` asks the rule per `(allele, read group)` row; the
   fold pools the rows first. The probe's published figures were therefore taken under a stricter
   rule for the samples carrying more than one read group, so they are a **lower bound** on what
   the shipped rule admits. Milestone D deletes the probe's copy and reconciles them; one sentence
   in the plan's D1 entry would record why a difference is expected there.
4. **B1 computes every sample's denominator and discards it.** C3 needs it again to ask the rule
   per `(sample, allele)` over the alleles the cap cut. The agent built both versions — recomputing,
   and carrying one `Vec<u32>` in `SelectionScratch` — and they agree to the bit. **Not applied
   here**, because the cap binds at 23 of 53,935 tomato loci and none of the trio's, so the second
   computation is needed at one locus in 2,300; the buffer would be a field with no reader for two
   steps. Recorded so C3 chooses deliberately.

## 5. Validation after the fixes

All in the container, on the tree that was committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. Run separately because `dead_code`
  fires there and `--tests` hides it, which is how the redundant allow was diagnosed.
- `cargo doc --lib --no-deps` — completes; 35 diagnostics, all pre-existing unresolved intra-doc
  links and redundant explicit link targets, none of them in `allele_candidates`.
- `cargo test --lib` — **4,236 passed, 0 failed, 14 ignored**, in 42.5 s, against 4,219 at
  `823c7b77`. The suite is 46 seconds rather than fifteen minutes since `main`'s `c8e297d3`, which
  is merged into this branch.

`cargo clippy --all-targets` is red on `main` with 14 pre-existing errors in five benches and
examples, none in `src/` and none touched here, which is why the gate is `--lib --tests` plus
`--lib` alone.

**Mutation re-check — twelve mutations, twelve killed**, each by the test written for it:

| mutation | test that fails |
|---|---|
| the share dropped, floor alone | `the_share_refuses_what_the_floor_would_admit` |
| the rule asked against the allele's own reads | `the_share_refuses_what_the_floor_would_admit` |
| `samples_clearing_the_bar` stuck at one | `every_sample_that_cleared_the_rule_is_counted` |
| the first sample's share wins | `the_largest_share_wins_when_it_arrives_last` |
| the reference is not asked the rule | `the_reference_row_is_folded_like_every_other_allele` |
| the zero-denominator assertion disabled | `a_sample_with_rows_and_no_reads_is_refused` |
| the sample-order assertion disabled | `a_sample_listed_twice_is_refused` |
| the no-rows skip made a stop | `a_sample_with_only_partial_reads_is_stepped_over_and_the_next_sample_is_folded` |
| the denominator swallows the partials | that test and `the_denominator_is_..._and_nothing_else` |
| the larger read-group row wins over the sum | `one_samples_two_read_groups_are_one_sample` and the oracle |
| the allele-order assertion disabled | `rows_out_of_allele_order_are_refused` |
| the scratch grows instead of resetting | `folding_a_second_locus_...` and `a_locus_no_sample_covers_...` |

Eight of the twelve were survivors before the fixes.

## 6. One thing worth keeping from how this review ran

**The design-fidelity agent's deliverable was code.** It wrote B2, C1, C2 and C3 against this
step's fold and compiled them — 41 tests green in its worktree, the step's own 32 unchanged plus
nine of its own, including the plan's C1 round-trip oracle and C3's second oracle verbatim. Both of
its substantive findings (item 1 and item 4 of §4 above) were *produced* by writing the later step,
not by reading the current one. That is the second milestone running where the highest-value agent
was the one asked to build forward rather than to check backward.

**And the two Blockers were both invisible to reading.** Every agent that read the fold judged it
correct, and it *is* correct. What was wrong was that four of its properties were unpinned, and the
only way to see that was to break each one and watch the suite stay green.
