# ng genotype prior — F1: the comparator, and the arm nothing could name

*Implementation report, 2026-08-23. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step F1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone F, on top of `0ea9e459`.*

## 1. What it is

`PlugInWrightPrior` — the comparator implementation of the step-8 seam. With it, **step 8 is
complete as a set of pure functions**: two implementations behind one trait.

Both turn a concentration into one log-probability per candidate genotype, and both apply the same
two-branch inbreeding mixture. They differ in one line of algebra: the default averages the
genotype probability over every frequency the concentration finds plausible; the comparator
collapses the concentration to a single frequency, `α_a / Σα`, and evaluates there.

Genotype probability is quadratic in the frequency, so `E[p²] = p̄² + Var(p)`: **plugging in
undercounts both homozygotes by exactly the variance and hands the heterozygotes twice that**
(spec §2.2). `Var(p)` is how badly the frequency is pinned down, so the gap vanishes in a large
cohort and dominates at one sample — the corner this caller commits to supporting.

## 2. Why it is kept

One measurement. On the GIAB trio, each sample called on its own at 5×, swapping the comparator
for the default took SNP genotype accuracy at true variants from **83.6% to 94.6%**, and the sites
where a sample carrying two copies of the variant was called heterozygous from **214 to 8**, with
the emitted variant set byte-identical (spec §2.2).

**And the trap it must not fall into.** That gain is not "marginalize" — it is the starting
concentration (spec §2.3). Production's plug-in path regularised its frequency estimate with a
reference pseudocount of 10; marginalizing over *that* gives 20:1 odds on a heterozygote where the
plug-in gave 22:1 — the same wrong answer, computed more expensively. So the comparator runs on the
same seed and supplies no pseudocount of its own, and
`the_row_is_hardy_weinberg_at_the_handed_concentration_and_nothing_else` builds its closed form
from the entries of the buffer that was passed in. Injecting the pseudocount fails that test and
three others.

## 3. What the tests pin, with sizes

At a single sample's seed `(1, θ)` at tomato's fitted diversity, `Var(p) = 2.9955e-4`:

| | default | comparator | difference |
|---|---|---|---|
| hom-ref | 0.99910 | 0.99880 | −Var |
| het | 5.9946e-4 | 1.1986e-3 | +2·Var |
| hom-alt | 2.9991e-4 | 3.5957e-7 | −Var |

The identity is exact rather than asymptotic, so it is asserted at 1e-9 and the worst residual is
5.3e-15. The comparator is **834 times less willing** to call a sample homozygous for the
alternative allele — asserted against `1 + Var/p̄²` rather than against the number.

**The oracle is production's own `wright_genotype_log_priors`**, which is legitimate because the
mixture *is* Wright's formula at two alleles and two copies: `(1 − F)p² + Fp = p² + Fpq`. Worst
disagreement over the sweep, 3.55e-15 nats. The two implementations converge as the frequency
becomes certain — 6.73, 2.86, 0.154, 0.00166 nats at `Σα` scaled by 1, 10², 10⁴, 10⁶ — and the
last step closes the gap 93-fold where the first closes it 2.4-fold, because the gap is only
proportional to `Var` once `Var ≪ p̄²`, that is once `Σα` is past about 1,700 here.

## 4. What the reviews changed

Two reviewers, each in its own worktree, each given the gate output.

**The finding that mattered: nothing tested the seam through the path a run uses.** Every test
named the comparator through its own file. Re-export the default under the comparator's name — a
one-line edit — and all 132 tests still passed while a run selecting between the two arms got the
same prior twice and reported that the two agree. **That is the bake-off's own failure mode, and
the number that would have vanished is the factor of 834 above.** There is now a test that holds
both implementations as `Box<dyn GenotypePriorModel>` through the folder's public names and
asserts the rows differ.

**A deferral that arrived unbuilt.** Step B2 deferred a `name()` method on the trait to "F1, where
the second implementation arrives", and recorded that deriving `Debug` would do in the meantime so
a `&dyn GenotypePriorModel` "can at least be printed". It cannot — the trait has no `Debug`
supertrait, so `Box<dyn GenotypePriorModel>` does not implement it and the stand-in never
compiled. `name()` is now on the trait. The seed's provenance is no substitute: `SeedRegime` and
`SpectrumMatch` describe the input the two implementations *share*, so neither can tell two runs
over one seed apart.

**A pin that could not see an unwritten row.** Every comparison in the file folds its departures
with `f64::max`, which ignores `NaN` — so a row whose entries were never written scored a
departure of zero and passed. Proved: filling only the first genotype left the module's own named
pin green. The row builder now refuses a non-finite entry, which closes all three sites at once.

**Six wrong numbers in my own doc comments**, the same failure as every preceding step, every one
a claim about my own fixture: the table's default hom-alt held `Var` instead of the probability
(2.995e-4 for 2.9991e-4); "833 times" is 834.08 — and it was the one assertion in a block whose
own comment says "against their closed forms rather than against a band"; the pseudocount
paragraph named the fully-reference genotype as moving 9.2 nats when it moves 0.0067, the 9.2
belonging to the genotypes carrying no reference copy; "nine thousand times this tolerance" is
nine trillion; and the comparator's heterozygote is 1.1986e-3, not 1.1985e-3.

**Two claims false at an endpoint.** The Wright equivalence is unqualified and fails at `F = 1`
exactly, where the two floor different things — production floors the finished probability and
lands at −690.78, this floors the weight and adds `ln 2pq`, landing 6.73 nats lower. At `F = 0.99`
they still agree to 9e-16. And the claim that the two implementations are "pinned there against
each other" at `F = 1` describes a test that would fail legitimately: they still differ by
`ln((Σα + 1)/Σα)`, 0.69 nats, because complete inbreeding scales the random-mating branch rather
than removing it.

**A floor with no test.** Deleting the floor on the frequencies survived the whole module. It is
defensive rather than live — it needs a concentration total near 1e300 — but every value check in
this folder is a `debug_assert!` that release compiles out, so in a shipping build it is the only
thing between a zero concentration entry and the `−∞` the seam's contract forbids by name. Now
covered, and the mutation was re-run against the new test: it is the only test that catches it.

**Mutation testing:** seventeen mutations across the two reviews. Three survivors, all proved to
change behaviour, all now covered. One changed no behaviour at all — removing the zero-count guard
— verified over 63 shapes as bit-identical, which is exactly what the code comment claims.

## 5. What Checkpoint F cannot claim, and it is not this step's to fix

The checkpoint's words are "both seam impls run under one trait; the recipe can select either".
The first half is true and now tested. **The second half describes a recipe that does not exist and
is not in any plan's scope.** A reviewer wrote the selection and it works — the trait is
object-safe, `Box<dyn>` and `Arc<dyn … + Send + Sync>` both compile and run, the latter across a
thread boundary — but the measurement the comparator is kept for needs three things and only one
of them is code: the calling loop (planned), candidate selection (unspecced — its own plan says
the spec "needs writing before it needs planning"), and a run-time selection mechanism, which
exists nowhere. The two plans agree with each other about this; what they agree on is that nobody
owns the run.

## 6. Two items raised for the milestone review

- **At `Σα = 9e15` it is the default that degrades, not the comparator.**
  `lgamma(Σα + 2) − lgamma(Σα)` cancels to 64.0 where the true value is 73.46 — one unit in the
  last place of `lgamma(9e15)`. The repeat-tract seed can reach that total: step E1's own
  `the_largest_total_any_shape_can_ask_for_is_finite` pins the worst at `2^53 − 1`. So Milestone E
  and Milestone B meet here, and neither noticed.
- **The recipe's field needs `+ Send + Sync`.** `arch/ng_step_interfaces.md` sketches
  `Box<dyn GenotypePriorModel>`, which will not cross a worker boundary in a chunk-parallel loop.
  Both unit structs are trivially both, so it costs nothing — but the loop plan should say it
  rather than let it be found at the call site.

## 7. Gates

Green in the container: `cargo fmt --check`, `cargo clippy --lib --tests --all-features -D
warnings`, `cargo test --lib`.
