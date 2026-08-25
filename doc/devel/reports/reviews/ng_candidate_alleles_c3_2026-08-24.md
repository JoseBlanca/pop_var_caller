# ng candidate alleles — C3: review and the fixes applied

*2026-08-24. Step C3 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
Milestone C's last step. Reviewed at `1041e30e` plus the step's working-tree diff, three agents in
three isolated worktrees.*

---

## 1. What was actually wrong

### 1.1 Nothing separated the earned count from the running pool — Blocker

Setting `earned_reads_cut_by_the_cap = leftover.num_reads` — the whole pool so far rather than this
allele's reads — **left all 89 tests green.** All four fixtures that produce a non-zero count gave
the affected sample **exactly one** dropped allele, where the two totals are the same number.

A sample carrying one error read on a rule-dropped allele *and* six on a cap-cut allele it earned
separates them: 7 against 6. **Fixed**, with the rule-dropped allele placed *before* the cap-cut one
in the merge table, because after it the running total and the right answer coincide again.

### 1.2 A shifted leftover was invisible — Blocker

Both alignment tests were written at C1, when every leftover was a zero, and both fixtures still
produced all-default leftovers. A fill that skips samples with no support rows and pads the tail
keeps the length right, slips past `LocusSelection::new`'s assertion, and **passed all 89** while
sliding the third sample's leftover onto the second — which is the exact failure the second test's
own doc comment claims to guard against. **Fixed** by giving the samples either side of a
partial-only sample different error mass, so position is observable in the values.

### 1.3 The denominator was untested against the plausible wrong one — Blocker

The existing test separates "this sample's compared reads" from "the allele's own reads". It does
not separate it from "the compared reads left *after* selection", and narrowing the denominator that
way **passed all 89**: every fixture was at a share of zero or at a depth where both agree.

**The wrong denominator is systematically smaller** — it is the right one minus what selection
dropped — so it asks less of every sample and no-calls samples the rule never meant to touch. A
sample at 6 reads in 16 separates them: the rule asks 8, the narrowed denominator asks 5. **Fixed.**

### 1.4 Spec §8's third assertion was never implemented — Major

Spec §8 and arch §1 both name **three** caller bugs this module holds as assertions in release. Two
were implemented. The third — a non-finite quality sum — was not, and **C3 is the first step that
reads `q_sum` at all**, so it is C3's to own. Run in release, a `NaN` on a dropped row produces a
`NaN` pool with no panic, and that `NaN` then enters every genotype's data likelihood so none is
preferred: the locus comes out called with nothing chosen and nothing failed. **Implemented, with a
test.**

### 1.5 The 400-sample test could not tell the cap rule from the pool rule — Major

At that fixture every sample with a non-empty pool was exactly a sample the cap cut, so the pool
rule produced the same count of 395. The test exists precisely because spec §4.1's reassurance was
measured at 63 accessions — and it held only *how many*, not *which*. **Fixed** by giving every
sample one read of a sequence nobody earns: all 400 now have a non-empty pool while only 395 lose
an allele. The pool rule now fails it.

### 1.6 A refactor accident, and two false claims of my own

**Lifting `one_run_per_allele` out of the fold put it between `summarise_alleles`'s doc comment and
`summarise_alleles`**, so a 70-line block including a three-item `# Panics` list came to document a
20-line helper that raises one of them, and the module's central fold had no doc comment at all.
Caught by the naming agent, fixed by moving the helper above the block.

**Two claims in the new prose were wrong.** The pool oracle said a bit-for-bit equality pinned the
addition order — its four masses are exact binary fractions, so every ordering gives exactly −5.25;
the test pins the *set* of rows and the value, not the order, and now says so. And
`one_run_per_allele`'s `# Panics` was unconditional where the iterator is lazy: a caller that stops
early can walk past a disorder it never reaches. Both current callers exhaust it, which is what
makes the guarantee hold, and the doc now says that rather than implying more.

### 1.7 Vocabulary

The module used "the bar" and "the admission rule" interchangeably across six steps while the
identifiers say *bar*, and C3's prose said only "rule" — including in an assertion message silently
reworded from B1's. The message is restored, and the module header now states the equivalence once
and, more importantly, states what is **not** interchangeable: the bar decides which sequences are
worth calling over and the cap decides how many a locus may carry, and which of the two dropped an
allele is what decides whether a sample keeps its genotype.

## 2. Checked and found sound

- **The leftover's arithmetic is spec §5 and §5.1 exactly.** The design-fidelity agent checked the
  cap-versus-rule keying by assertion over 24 synthetic cohort configurations rather than one
  fixture — `capped || earned == 0` and `earned <= num_reads` never fired from 1 to 2,000 samples
  and 3 to 300 reads a position.
- **`select_generic` is complete for Checkpoint C's wording**, walked clause by clause through spec
  §3, §4, §5 and §6 — the bar, the supplied numerator, the cap and its four ranking keys, the pool
  and both counts, the reference invariant, both verdicts, reference-only as first class. What
  remains is on the plan's own out-of-scope list or in Milestone D.
- **`LocusSelection`'s surface is enough to build arch §3.2's hand-off**, compiled: nothing beyond
  `alleles()`, `verdict()`, `unmatched()`, `remap().candidate_for()` and `UnmatchedSupport`'s three
  fields was needed. One note for the calling loop's plan: arch types `supported` as a borrowed
  slice, but the re-keyed rows are new values, so the loop needs an arena.
- The fold's behaviour is unchanged by the refactor; a sample whose every allele was dropped keeps
  its genotype; saturation is correct; and a cap-cut allele nobody earned is unreachable.

## 3. Raised, not applied — for Checkpoint C

**The range numbers, which are the substantive ones.** All synthetic — the shipped measurement
harness still carries its own copy of the rule until D1:

- **The cap stops being a safety valve well before 400 samples.** It binds at essentially every
  locus from 400 on, with merge tables reaching 145 to 1,953 alleles. Spec §4.1's "measured, that is
  rare: 23 of 53,935 tomato loci" is a fact about 63 accessions and does not carry.
- **Samples emitted as missing:** between 6 in 10,000 and 9 in 1,000 at 63 samples and 3 reads a
  position; between 5 in 1,000 and 4 in 100 at 400 samples and 3 reads; **between 1 sample in 20 and
  1 in 8 at 400 samples and 30 reads.**
- **The cost peaks in the middle of the depth range** — 4.2% at 3 reads, 13.3% at 30, 10.4% at 300
  (400 samples) — because at 300 the 5-in-100 share removes alternatives before the cap sees them.
- **The pool at 3 reads a position is 3.3–3.6% of a sample's reads**, against 0.75–1.0% at 30 and
  300. Spec §5 measured 0.36 in 100 on 63 tomato accessions; a heterozygote drawing 1 read of 3
  fails a floor of 2 and the whole allele goes into the pool.

Also: **arch §3.2 is stale** — it says `read_likelihoods.md` §2.1 "declares `unmatched_q_sum` alone"
and that the count must travel with it, where that document now declares
`genotype_must_be_missing: bool`. And the items C1 and C2 raised stand.

## 4. Validation after the fixes

- `cargo fmt --check` clean; both `clippy` gates clean;
- `cargo test --lib` **4,288 passed, 0 failed, 14 ignored** in 42.9 s, against 4,276 at `1041e30e`.

**Ten mutations, ten killed**, four of them survivors before the fixes: the earned count as the
running pool (now fails 6), the denominator narrowed to the surviving alleles (1), the
missing-genotype rule keyed on the pool (8), and a leftover fill that skips row-less samples and
pads (1).

## 5. What six steps of this have taught

**Every Blocker on this plan — eleven of them — was a test that could not fail, and not one was
wrong code.** The fold, the ranking, the admission pass, the cap and the leftover each computed the
right answer on every input any reviewer could build.

They also had one recurring shape, and it is now precise enough to check for: **a fixture built at a
size, depth or cohort where the term under test is not the term that decides.** Shallow enough that
the rule's floor decided and its share never did. Single-sample enough that a cohort sum matched a
per-sample rule, and a cohort total matched a within-sample share. One dropped allele per sample, so
a per-allele count matched a running total. One question would have caught all eleven: **at this
size, would the simplest wrong rule give the same answer here?**
