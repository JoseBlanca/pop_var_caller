# ng genotype prior — E2: the shape both consumers read, and the support it is not

*Implementation report, 2026-08-23. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step E2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone E, on top of E1
(`fcf8cdb6`).*

## 1. What it is

One public function and five tests. `fill_seed_share_per_candidate` fills a caller's buffer with
the share of the genotype prior's starting mass each candidate at a repeat tract carries — the
same shape [E1](ng_calling_prior_e1_2026-08-23.md)'s `fill_ssr_seed` scales into a concentration,
before it is scaled, normalised to sum to one.

It exists so the prior's belief about which repeat lengths are common has **one implementation**
behind both its consumers (arch §5): the genotype prior itself, and the read likelihood's
contamination term, which needs a distribution over the lengths that contaminating DNA might
carry and has no measured one to use.

## 2. The finding: it is per candidate, and the term it feeds is per observed length

The mixture the export serves is `(1 − λ − c)·copy-mixture + λ·uniform + c·seed(o)`, where `o` is
an **observation** — a read. This buffer has one entry per **candidate**. Those coincide only
when every candidate is a distinct length and every read lands on one of them, and three cases
break that:

| case | what it does to the export |
|---|---|
| two candidates of one length, each taking the rung's full share | the modal *length* arrives with **0.8** of the mass where the geometry says **0.667**, on a tract with the mode spelled twice and one length above, at the fallback decay |
| a read at a length the candidate prune dropped | no entry — while the mixture's sibling uniform term is spread over every length the stutter model can reach from a candidate, a strictly larger support |
| a censored read | no length at all, only a lower bound |

The first is deliberate as a concentration and open as spec Q3; it is only wrong when the same
numbers are read as a claim about lengths. **None of the three is this function's to settle.** The
step's response is to name the buffer for what it holds — the arch sketch called it
`seed_length_distribution` — to record all three with their sizes as `OPEN:` in arch §5, and to
put the same marker on `SsrContamination::length_distribution` in
[`read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §4.1, which is the field that has to
say which support it means. The alternative was to leave the likelihood step to discover it at
the point of use.

## 3. Departure from arch §5, recorded there

**Named `fill_seed_share_per_candidate`**, not `seed_length_distribution`. Two reasons: the
module's `fill_*` convention for buffer-fillers, and — the one that matters — *length* is the
inaccurate word, for the reason above. Arch §5 now carries the new signature and the `OPEN:`.

## 4. What the review changed

One reviewer, in its own worktree, given the gate output. Three should-fixes and seven nits, all
applied:

- **The support mismatch above**, which the first draft's doc comment asserted away by calling the
  buffer a length distribution.
- **"Its cost is one pass" was three**, and contradicted the sibling function's own doc sixty
  lines up, which correctly says four for the same code plus a scaling loop. Worse, the sentence
  omitted that a loop calling `fill_ssr_seed` at the same locus **already holds this shape** — it
  is the concentration over its own total on a seed, and exactly the buffer the refusal hands
  back on a refusal. The doc now says both, and says who the export is for.
- **"Once per locus and frozen" does not survive a discovery round**, which appends candidates
  mid-locus; a frozen candidate-parallel buffer is then one entry short and the length assertion
  cannot fire, because by construction it is not refilled. Faithful to what
  `read_likelihoods.md` said, and still wrong, so both documents now say it.
- **The wrong-mode trap was not carried onto the public export.** Measured: passing the reference
  allele's repeat count as the cohort's mode moves the mode's share from **0.711 to 0.108**, a
  factor of 6.6, while the reference's rises from 0.089 to 0.862. Now documented and pinned.
- **The proportionality test was pinned at one spread and one decay**, and a mutation that made
  the two paths diverge only at candidates three or more repeats from the mode survived it —
  because that fixture's widest offset was two, and the sums-to-one test renormalises. It now
  sweeps four spreads and five decays, each at half of whatever that combination's ceiling can
  hold so nothing is refused.
- **"The two are one implementation" has an exception**: a single-candidate locus short-circuits
  in `fill_ssr_seed` and never reaches the shared code; the two agree there only because
  `ALPHA_REF` is 1. Now stated.
- **The stand-in's second half was unsupported.** "False of an adapter or a foreign organism" is
  not a live failure mode — the contamination fraction is fitted from low-level alternative reads
  concentrated on the panel's own segregating alleles, so a foreign organism does not raise it,
  and the sibling spec prices the case at zero: reads unlike the cohort's mode fall to the outlier
  floor, which is where they go today. Dropped.
- **The three documents still naming a function that does not exist** — arch §5, arch
  `read_likelihoods.md` §4.1, and this plan's own E2 line, still unticked. All three corrected.
- **`mod.rs`'s front page said "two functions"** and now has three public ones. Noted there.

**Mutation testing.** Seven mutations. Five caught. Two survived and were proved to change
behaviour: the divergence-beyond-offset-three above, and the wrong-mode substitution — both are
now covered by the strengthened tests. One further mutation changed no behaviour at all
(returning the Simpson index rather than discarding it), which is a signature widening no test
could observe and so not a coverage gap.

## 5. What the tests pin

| test | what it holds |
|---|---|
| `the_shared_export_sums_to_one` | four spreads at five decays, worst departure from one under 1e-15 |
| `the_shared_export_is_the_seed_divided_by_its_total` | the same implementation stands behind both consumers, swept so a divergence in the tail cannot hide |
| `two_spellings_of_one_length_carry_that_lengths_share_twice` | the first of the three support gaps, with its size — 0.8 against 0.667 |
| `the_mode_is_the_cohorts_and_nothing_here_can_check_it` | the wrong-mode trap, with its size — a factor of 6.6 on the mode's share |
| `a_mis_sized_buffer_is_refused_by_the_shared_export_too` | the length check, which the export inherits rather than makes |

## 6. Gates

Green in the container: `cargo fmt --check`, `cargo clippy --lib --tests --all-features -D
warnings`, `cargo test --lib`.
