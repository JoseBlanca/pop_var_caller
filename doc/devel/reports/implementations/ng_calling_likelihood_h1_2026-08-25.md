# ng read likelihoods — H1: the STR row, and three tests that could not fail

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step H1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md).*

## 1. What it is

`src/ng/calling/likelihood/ssr.rs` — **one sample's log-likelihood for every candidate genotype at
one repeat tract**. A read arrived one of two ways: copied from one of this individual's own copies
of the tract, or from somewhere else entirely. For a genotype `g` with copy counts `k_a` over a
ploidy `P` (spec §2.1, §4.5):

```text
log Lg(g)  =  Σ_o  n_o · log[ (1 − λ) · Σ_a (k_a / P) · Lr(o | a)  +  λ · U ]
```

`Lr(o | a)` is the emission built in milestones F and G — `emission` for a read that spanned the
tract, `censored_emission` for one that ran off its own end. `λ` is `DEFAULT_OUTLIER_WEIGHT = 0.01`,
**inherited from production and declared inherited**, and `U` is uniform over the tract lengths the
model can reach.

**Contamination is not here.** Spec §4.5.1's third term, and the per-locus computation of how many
lengths `λ` is spread over, are H2's; both arrive as parameters for now.

## 2. What a row costs, and how that is pinned

**`observations × candidates` emission calls, not `× genotypes`.** Each emission is computed once
into `SsrRowScratch`'s cache and read by every genotype carrying that candidate. Spec §8 calls that
the design rather than an optimisation.

An instrumented model counts the calls, and the fixture is chosen so the three plausible costs are
three different numbers: at three observations, three candidates and a diploid — six genotypes, nine
carried-allele slots — **the design costs 9, recomputing per genotype would cost 18, and
recomputing per carried allele would cost 27.**

**The cache is indexed by an observation's position in the whole slice, partials included**, because
both of the evidence type's filters enumerate the whole slice and then filter. A dense counter
inside a filtered loop addresses the wrong row for every observation above a partial — measured, it
returns an all-`NaN` row, and until the review it passed every test here.

## 3. What the review found, and what it changed

**Three category agents, each in its own worktree**, returning **three Blockers and six Majors**.
Two of the Blockers and much of the rest were the same shape: *the code is right and nothing could
have told you*.

### The three Blockers were all fixtures that could not fail

- **Nothing pinned `k / P`.** Deleting the division by the ploidy — so the copy weights sum to the
  ploidy instead of to one — passed all four original tests while moving every entry of the row by
  7–8 nats, about a factor of a thousand in likelihood. The tests pinned *orderings*, and every
  ordering survives a weighting wrong by a constant factor. Spec §12's seventh test names a
  hand-computed biallelic diploid as its oracle and only the ploidy-4 half had been written; the
  hand calculation now exists and agrees **to the bit**.
- **Every scoring context in the file was interchangeable.** The fixture gave all candidates the
  same slippage row and the same substitution rate, so hoisting the per-candidate lookup out of the
  candidate loop — the one thing spec §4.4 forbids by name — left every row bit-identical. The
  fixture now keys both parameters to the candidate's **repeat count**, across two read groups, so
  the two axes of the lookup each have something to fail against.
- **No fixture held a partial observation**, so `censored_emission` was never reached through the
  row at all. Routing the partial loop to `emission` passed; so did the cache-index defect above.

### The Majors

- **`outlier_weight` was the one input the row did not check**, and both ends of the range produce a
  number rather than a crash: at λ = 1.5, two of three genotypes came back `NaN` and the third a
  plausible −12.74; at λ = 0 a read nothing explains takes every genotype to `−∞`, whose differences
  are `NaN` — the collapse spec §4.5's junk term exists to prevent. Now asserted strictly inside 0
  and 1, with a `should_panic` test at each end.
- **The order-independence test never permuted the candidates**, despite its name and spec §12's
  eighth test. It does now — and the genotypes have to be matched by their copy counts rather than
  by mirroring the index, because at three candidates the reversal maps genotype indices
  `0→5, 1→4, 2→2, 3→3, 4→1, 5→0`.
- **The copy-share table was built twice**, in the very step that moved `MAX_PLOIDY_COPIES` into a
  shared home. Both rows now call one `copy_shares`, which also carries the ploidy check.
- **Three doc comments stated arithmetic that did not hold**: the naive emission cost given as 18
  under a rule that yields 27; a measurement quoted at four candidates that came from the
  three-candidate fixture; and `mod.rs`'s claim that the cache buys "a factor of ten at six
  candidates and a diploid", which is 3.5 — a figure spec §8 records having already corrected once,
  and which `generic.rs` already carries correctly.

### What the independent check confirmed

An agent transcribed spec §2.1 into an oracle that calls the emission model itself, walks
genotype-major instead of observation-major, and shares no code with the row, then compared **6,192
cells** across ploidies 1/2/4, candidate counts 1–4, six observation mixes, four outlier weights and
three reachable-length counts. **Worst disagreement: zero units in the last place — every entry
bit-identical.** The four edge cases the row's own tests did not reach all behave: an empty
observation list gives exactly `0.0`; one candidate gives one genotype, bit-identical; an
observation nothing explains gives exactly `n·log(λU)` in every genotype; ploidy 1 is bit-identical.

## 4. Spec §12's sixth test asked for something no implementation can give

**"The junk term cancels for a read nothing explains ... bit-for-bit."** It cannot: the term is added
to each genotype's running total, and `(a + k) − (b + k)` is not `a − b` in floating point however
carefully `k` is computed.

**The unit matters more than the number, and getting it wrong is how this test was first written.**
Counting units in the last place *of the difference between two genotypes* measures the true
rounding error scaled by `|entry| / |separation|`, and that ratio is set by how many junk reads there
are — the same fixture reports **16 units at 3 junk reads and 3,072 at 300**, with nothing about the
row having changed. A first correction to this document quoted 64 and was measuring exactly that
artefact. Measured relative to **the entries' own magnitude**, which is where the rounding happens,
the worst disagreement over every cell the sweep reaches is **one `f64::EPSILON`** and stays there
under the same sweep — the same shape of bound spec §12's eighth test means by *"the same relative
bound as test 9"*.

**Spec §12's sixth test is corrected to that**, in the same style as items 8 and 9, which already
carry this correction for the same reason. The property is real; the unit was wrong.

*(A two-candidate fixture agrees bitwise, and the first draft of this test used one — exactly the
luck item 9 warns about in its own words: "a single fixture will often agree to the last bit and
that is luck rather than arithmetic".)*

## 5. Two departures from the architecture's sketch, and the reason for each

- **`candidates` arrives as built `SsrCandidate`s, not as `CandidateAlleles`.** A candidate's repeat
  count is not derivable from its bases — an interrupted tract's byte length divided by the period
  is not how many repeats it holds — and F1 settled that this type *consumes* a measurement and
  never makes one.
- **The four things the locus contributes are one `SsrLocusParameters` value**, not four arguments.
  The nine-argument form the architecture sketches trips clippy's argument-count limit on the newly
  merged 1.98 toolchain, and the four travel together anyway: they are built once per locus and H2
  adds a fifth to them.

## 6. Validation

Run in the dev container on this worktree, at the **1.98 pin merged from `main` today**:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features --tests -- -D warnings` — clean.
- `cargo test --lib` — **4,520 passed, 0 failed, 14 ignored**; **207 in `ng::calling::likelihood`**,
  of which **13 are the row's**.

**Four mutations re-run against the repaired tests, all killed**: the copy weights dropping `/ P`
(22 failures), the context lookup hoisted out of the candidate loop (2), a partial scored as a
complete read (1), and the cache indexed by a dense counter inside the filtered loop (1). Each
restore was verified by checksum before the next ran.

## 7. What is still open

- **Contamination and the per-locus reachable-length count** are H2's, as planned.
- **`emission` scores 22 length changes the support calls unreachable** (G1's report §9) — unchanged
  here, and it reaches this row through the emission cache like any other score.
