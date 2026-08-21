# Applying the B3 review — ng calling foundations

*2026-08-21. Branch `ng-calling-foundations`. Input:
[`ng_calling_b3_2026-08-21.md`](../reviews/ng_calling_b3_2026-08-21.md). Every finding in that
report is accounted for below.*

## Findings table

| id | severity | decision | status |
|---|---|---|---|
| M1 — the copies/alleles check is tested in one direction only | Major | Apply | **Applied** |
| M2 — `passes == 1` is never built successfully | Major | Apply | **Applied** |
| M3 — the seed marker is pinned on a SNP/indel fixture | Major | Apply | **Applied**, and the guard the review verified came with it |
| M4 — the invariant is `&mut`-reachable after construction | Major | Apply | **Applied** (two private fields, two accessors) |
| M5 — a positional constructor for a check a literal skips | Major | Apply with adaptation | **Applied** — M4's remedy subsumes the bypass; the `#[allow]` gained a reason |
| Mi1 — an empty `per_sample` is neither rejected nor pinned | Minor | Apply | **Applied** |
| Mi2 — the module inventory is stale | Minor | Apply | **Applied** |
| Mi3 — a citation to a section that does not exist | Minor | Apply | **Applied** (`§8` → `§1.4`) |
| Mi4 — `weakest_provenance` names an ordering that does not exist | Minor | Apply | **Applied** (documentation) + follow-up |
| Mi5 — a test doc claims a case its check cannot reach | Minor | Apply | **Applied** |
| Mi6 — the seed marker's ceiling stated as settled where the arch marks it provisional | Minor | Apply | **Applied** |
| Nit — `new` has no summary line before `# Panics` | Nit | Apply | **Applied** |
| Nit — `region.start <= region.end` unchecked | Nit | Apply | **Applied** |
| Nit — the `50` literal will need `DEFAULT_MAX_PASSES` | Nit | Defer | **Deferred** to the loop plan |
| Nit — `passes` reads as a collection; `seed_diversity_unreachable` parses oddly | Nit | Dispute | **Won't fix** — the arch fixes both names |
| Recommendation — `converged` implies `passes >= 2` | — | Dispute | **Won't fix**, on the spec (§ below) |
| Out of scope — `Provenance` has no ordering | — | Defer | **Deferred** (blocks the loop step) |
| Out of scope — `arch/read_likelihoods.md` §8 does not exist | — | Defer | **Deferred** (a design-doc fix) |

## M4 and M5 — two private fields, and what that bought

`alleles` and `cohort_expected_copies` are now private, read through two `#[inline]` accessors. The
other six fields stay public.

**This was cheap because only two of the eight fields participate in an invariant.** The review's
own probe is why it was necessary: `alleles` was public and `CandidateAlleles::admit` takes
`&mut self`, so

```
PROBE-MUTATE: after new(), admit() widened alleles to 3 against 2 copies; passes now 0
```

— a value that passed every check on the way in could be broken afterwards, from a call site that
reads as ordinary use of the discovery round. A struct literal skipping `new` is visible in a diff;
this is not.

It also closes M5's bypass half for every caller outside `calling/`, since a struct with a private
field cannot be written as a literal there. Inside `calling/`, Rust's privacy reaches descendant
modules, so the future inference sub-module can still write one — and that is the right shape: the
loop is the type's builder, and the consumers are the ones who must not be able to break it.

**The eight-argument shape stays**, and the `#[allow(clippy::too_many_arguments)]` now carries a
`reason`: the architecture fixes this type's fields as a flat list of eight, so grouping them to
satisfy the lint would be a design change rather than a refactor. The review's alternative — a
struct literal plus a `checked(self)` method — was not taken because it requires all eight fields
public, which is what M4 is about. **The residual transposition hazard is real but measured**: the
review found that the fixtures set the two flags to opposite values, so a swap inside the
constructor fails two tests, not zero.

## M1, M2, M3 and Mi1 — the four checks and their tests

`new` now makes five assertions where it made two:

| check | why |
|---|---|
| copies are one entry per allele, **either direction** | the prune shrinks the table, so an un-re-cut copies vector is *wider*, and every consumer indexes by `AlleleId` |
| `per_sample` is not empty | a cohort has at least one sample and every sample is called, so an empty list is a dropped result |
| `passes > 0` | zero is a counter never incremented |
| the seed marker is not set on a SNP/indel locus | that marker belongs to the STR prior's seed |
| `region.start <= region.end` | the region reaches the output's position column |

**M1's test covers both directions**, using `catch_unwind` around a small builder so one test can
assert that three alleles refuse two copy counts *and* that two alleles refuse three.

**M2's test builds a locus at `passes: 1` and asserts it exists.** The previous fixture at that
value panicked on the earlier copies check first, so the accepting side of the bound was never
observed.

**M3's fixture moved to the repeat path.** A new `str_two_allele_locus()` helper builds a
`LocusKind::Ssr` table with a motif and both flanks, and the capped-locus test now carries the seed
marker there. A second test asserts a SNP/indel locus is refused for carrying it — which is the
guard the review verified lands clean once the mis-specified fixture is fixed.

**Verified against the two mutations the review found surviving:**

```
--- assert_eq! weakened to  copies.len() >= alleles.len()
test ng::calling::tests::a_locus_cannot_carry_copies_of_a_width_its_alleles_do_not_have ... FAILED
test result: FAILED. 20 passed; 1 failed; 0 ignored; 0 measured; 3962 filtered out

--- passes > 0 tightened to passes > 1
test ng::calling::tests::a_locus_that_settled_on_its_first_pass_is_a_locus ... FAILED
test result: FAILED. 20 passed; 1 failed; 0 ignored; 0 measured; 3962 filtered out
```

## The disputed recommendation: `converged` does **not** imply two passes

The reliability agent recommended asserting that a converged locus took at least two passes,
reasoning that convergence is a comparison between two passes and that a locus reporting
`converged: true, passes: 1` would mean the loop compared against a reused scratch buffer's stale
contents. It filed this as a recommendation and an author question rather than a finding, which was
right — and the answer is no.

`spec/calling_em_loop.md` §2's pseudocode puts an **initialisation E-step outside the frequency
loop**: "one E-step with NO prior — reads only — then sum: the cohort's expected copies", and §3
gives the reason — the prior needs the cohort's expected copies, which do not exist until some pass
has produced them, so the first pass runs on the reads alone. The frequency loop's first pass
therefore has a freshly computed previous estimate to be compared against, not stale scratch. A
locus that converges on it means the prior barely moved the copies away from what the reads alone
said, which is a real outcome — and one to expect exactly where this caller is weakest, at a cohort
of one with three reads a position.

The `passes` field's doc now records this, so the next reader does not re-derive the argument, and
M2's test pins the accepting side of the bound.

## The documentation fixes

The module header's inventory names all four types. `weakest_provenance`'s citation is `§1.4` of
the read-likelihoods architecture, not the `§8` that does not exist there. That field's doc also
now says what the review found: **`Provenance` defines no ordering**, so "weakest" is not yet
computable — four names are not a scale, and where a *supplied* value sits against a *fitted* one is
open. The step that first has to compare two of them must settle it.

The seed marker's doc says the ceiling is **provisional**, which is what the architecture says,
rather than stating it as a settled rule; and it now says outright that `new` refuses the marker on
the SNP/indel path, rather than only that it is never set there.

`new` opens with a summary line. The cohort-of-one test says "where this caller is weakest" rather
than "the hardest case", since `CLAUDE.md`'s claim is about a single **low-coverage** sample and
this fixture fixes no depth.

## Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.85s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 21 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --all-targets --all-features` | 101 | lib `test result: ok. 3972 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 663.05s`; every integration-test binary ok; the run then hits the same **pre-existing** `benches/psp_writer_perf.rs:386` panic A1's review records |

Twenty-one tests in the module where the review saw seventeen. Each mutation run above restored the
file from a copy afterwards.

## Follow-ups this run created

1. **`Provenance` needs an ordering** before anything can compute a "weakest" one
   (`src/ng/parameter_estimation/mod.rs`). It blocks the loop step, not this one.
2. **`arch/read_likelihoods.md` §8 is cited in two design documents and does not exist** — that
   document runs §0 to §7 and the rule is its §1.4. A design-doc fix, out of this loop's remit.
3. **The `50` in the capped-locus test** should read `DEFAULT_MAX_PASSES` once
   `CallingLoopConfig` lands ([`calling_loop.md`](../../ng/impl_plan/calling_loop.md)), so the two
   cannot drift.
4. **`converged` against the pass cap** stays unchecked here: a capped locus should report the cap,
   but the cap is run configuration this type does not see.
