# Review: the portable-float A1 call-site inventory

*Reviews `doc/devel/reports/implementations/portable_float_A1_call_site_inventory_2026-09-13.md`
(step A1 of `doc/devel/implementation_plans/portable_float.md`), in the worktree
`pop_var_caller-portable-float`, branch `portable-float`. 2026-09-13. Nothing under `src/` was
changed and no build was run. The reviewer's scratch files are in `tmp/review_A1/`.*

## Verdict

**The inventory can be used for the decision it was written for.** The count is complete and
exactly right: an independent search, written separately from the report's, finds the same 202
shipped calls and 232 test calls, line for line. Every hot call the report names is in the loop
it says, and none of the per-run, per-table or per-locus calls I traced sits in a hotter loop
than the report says. No finding below is a Blocker.

There are four Fix findings. One call the report counts as reachable from the binary is not:
the reachable count should be 170, not 171. The first-pass command printed in §1.1 does not
reproduce its own count. The summary's description of the fit's innermost loops multiplies
together loop dimensions that belong to different calls. And one per-table class is too cold
at repeat-tract loci with more than 16 candidate alleles. There are also seven Nits. Two of the
Nits matter for step A2's weighting even though they correct nothing: several of the per-read
calls could be moved out of their loops, because their arguments do not change inside them.

## What was checked, and held

1. **Completeness.** A Python script, `tmp/review_A1/classify.py`, strips comments and string
   literals, then matches the method form, the path form (`f64::exp)`) and `libm::`. It also
   decides test code by brace-matching every `#[test]`, `#[cfg(test)]` and
   `#[cfg(any(test, feature = "bench-fixtures"))]` item, and it follows `mod x;` declarations
   made under those attributes. Its shipped list, as `file:line  function  count`, is
   byte-identical to the report's `tmp/portable_float_A1_sites.tsv`
   (`diff review_A1/mine_ship.tsv review_A1/theirs.tsv` → `IDENTICAL`). Totals per function
   are the same as §1.3. My search list was wider than the report's: it added `gamma`,
   `sin_cos`, `erf` and others. The only extra matches were 8 `self.gamma(…)` calls, which are
   a method of the project's own random-number generator and sit in test and bench-fixture code.
   A second search for these names not followed by `(` found only the one `.map(f64::exp)` the
   report already counts. No `mul_add`, `f32` transcendental call, `num_traits::Float` or
   other maths crate appears. The vendored `vendor/noodles-cram` has no such call either.
   The 44 files with a match and the 34 files with a shipped match agree with §1.2 and §1.3.
2. **Citations.** I checked about 150 cited `file:line` references for the stated call and
   the stated enclosing function (`tmp/review_A1/cite.sh`). That includes every per-read,
   per-genotype, per-base and per-fit-iteration row, the loop-bound lines, the call-chain
   lines and the constant-definition lines. All hold, apart from the items listed below.
3. **Hot-call classes.** I traced these sites through their callers:
   - `fit.rs:2098` and `:2132`, which sit in `one_position`, called per census position from
     `expectation_pass`;
   - `ssr_fit.rs:2073`, `:2332` and `:2388`: `Scorer::score` runs `ln_tract` over 256 points
     per golden-section evaluation, and `dirichlet_points` runs 12 sticks × 256 points × up to
     60 bisection steps whenever the spectrum or the concentration changes;
   - `generic.rs:810`, `:830`, `:903` and `:912`, reached per observation from
     `assemble_genotype_log_likelihood_row` and `score_partials`;
   - `ssr_emission.rs:629` together with `emission.rs:558-559`: `substitution_probability`
     builds a new `FlatEmission` on every call;
   - `ssr_unit_robust.rs:278-288`: `SlipCosts::from_model` runs once per `delimit` call, and
     the generator calls `delimit` once or twice per kept read (`locus_generation/ssr.rs:932,
     942, 2519`).

   All classes are right. The slippage-refit round and the discovery round, which would
   multiply the repeat-tract emission cost, are both off in the shipped configuration
   (`summarise_condition.rs` around line 2440, `max_rounds` 0 and discovery `Off`).
4. **Colder classes.** I spot-checked 14 of them, and all hold apart from Finding 4:
   - `stutter.rs:486`, `quality/mod.rs:191`, the `artifact_correction.rs` tests, `genetics.rs:54`
     and the `CountPriorMemo` around it, `genetics.rs:86-88`, `locus_score.rs:226, 246-266,
     454, 466, 479-480, 581`, `paralog/calibration.rs:109-110`;
   - `calibration.rs:238`, `depth_bins.rs:233/299` (every `for_census` caller builds a writer,
     a digest or a configuration, never a position), `census_moments.rs:410`,
     `contamination.rs:613, 880, 1311`, and the share and slippage curves.
5. **"Not reachable" claims.** I checked 10 groups: the five sibling delimiters (their types
   are named only under `#[cfg(test)]` at `locus_generation/ssr.rs:24` and in `examples/`),
   both marginal aligners, `PlugInWrightPrior`, `ReadGroupCalibration::log_scale`,
   `MintedErrorTotals::geometric_mean` and `JointFitConfig::coverage_odds`. All hold.
6. **Constants from literals.** I confirmed 9 of them:
   - `SsrGenerator::new` always stores `StutterModel::hipstr_shipped()`
     (`locus_generation/ssr.rs:2331`), and no setter exists. The shares
     (`stutter.rs:310-315`) pass unchanged through the clamps, so the four slip costs take
     `ln 0.88`, `ln 0.0475` twice, and `ln 0.05`.
   - The flat-gap literals (`:183, :196, :201`).
   - The paralog defaults (`finish.rs:211`, `model_params.rs:14, 27, 45`; 7 carrier
     configurations × 40 points = 280 cells).
   - The depth ladder, `(124/8)^(1/11)`.
   - `CANDIDATE_ALTERNATIVES = 3` and the multiplicity values 1 to 3 (`fit.rs:2047, 2055`).
   - The share-curve fallback (`ssr_fit.rs:1799, 1805`).
   - `QUADRATURE_POINTS = 256`, `TAU.ln()`, and `log_factorial`'s small integers.
7. **Input-range lines.** I checked 12 of them, and all bound what they claim apart from
   Nits 9 and 10:
   - `emission.rs:535-537, 139`, `generic.rs:26, 31, 157-161`, `mod.rs:154, 207, 498-499`;
   - `ssr.rs:155, 254, 696-697`, `prior.rs:45-56, 316`, `ssr_fit.rs:2241, 2348, 2451`;
   - `fit.rs:2820, 2758-2759, 2944`, `stutter.rs:287, 857, 868, 772`.
8. **Summary numbers.** The function × class table re-derives exactly from the working file,
   weighted by the count column. The same is true of "12 per-read", "22 per-genotype",
   "70 per-fit-iteration, 3 from the paralog prior", "30 not reachable, 1 gated", "15 calls
   reach `libm::lgamma`" and every row of the per-file table. The only exception is the
   reachability change in Finding 1.

## Findings

### 1. Fix — `depth_bins.rs:241` is not reachable from the binary

The report marks `ladder_tops`'s `powi` as reachable and per-run (§2.7, §4, and the working file
row `src/parameter_estimation/depth_bins.rs:241 powi 1 per-run yes yes`). `ladder_tops` has one
shipped caller, `DepthBinEdges::new` (`depth_bins.rs:270`), and `new` has one shipped caller,
`impl Default for DepthBinEdges` (`:201-204`). Every call of either is in a test or an example:

```
$ rg -n 'DepthBinEdges::new\(\)|DepthBinEdges::default\(\)' src
src/parameter_estimation/depth_bins.rs:463:  let edges = DepthBinEdges::new();      # tests start at :443
… (all further hits at :493–:781, inside the same test module)
$ rg -n 'DepthBinEdges' examples | rg 'new\(\)|default\(\)'
examples/ng_minted_error_means.rs:294, examples/ng_joint_duplicated_drifting.rs:208
```

No struct holding a `DepthBinEdges` derives `Default`. The two fields of that type are
`CensusWriter.edges` (`census.rs:2258`) and `JointFitConfig.edges` (`fit.rs:554`), and the
second has a hand-written `Default` that calls `for_census` (`fit.rs:614`). The shipped run
builds its ladders only through `for_census`, which calls `widening_ratio` (`:233`) and its own
`powi` (`:299`). Both of those stay reachable.

**Correction.**
- Mark `depth_bins.rs:241` *not reachable*.
- Change the reachable totals from 171 to 170 and the unreachable ones from 30 to 31, in the
  Words list, the Summary sentence and the per-file table, where `depth_bins.rs` becomes
  3 shipped, 2 reachable.
- In the function × class table, `powi` per-run becomes "2 (1)", the `powi` total "9 (7)", the
  per-run total "46 (45)" and the grand total "202 (170)".
- In the §4 row, *reachable* becomes "no".

The decision is unchanged: the call is per-run either way.

### 2. Fix — the §1.1 comment-filter command does not print 434

The second command in §1.1 is annotated `# 434 remain`. As printed, without `-o`, it counts
matching **lines**:

```
$ rg -n --no-heading -e "(\.|f(32|64)::)($F)\b(\(|\))|libm::" src | rg -v '^[^:]+:[0-9]+:\s*//' | wc -l
381
```

The occurrence count needs a third stage:

```
$ … | rg -v '^[^:]+:[0-9]+:\s*//' | rg -o -e "(\.|f(32|64)::)($F)\b(\(|\))|libm::[a-z_0-9]+" | wc -l
434
```

The 434 itself is right, and my independent count agrees with it.

**Correction.** Add the `| rg -o -e … | wc -l` stage to the printed command, or annotate
it "381 lines, 434 occurrences".

### 3. Fix — the Summary merges the loop dimensions of different calls

The summary's per-fit-iteration bullet says the innermost loops "run per census position per
sample per quadrature node per pass (`joint/fit.rs:2098-2456`), per repeat tract per sample per
genotype per quadrature point (`joint/ssr_fit.rs:1948-2073`)". **No single transcendental call
runs over all of those dimensions at once.** The sample × node loop in `fit.rs:2115-2127` and
the per-sample loop in `ssr_fit.rs:2031-2069` only multiply and add. The calls sit either side of
those loops:

- `fit.rs:2098` (`exp`) runs per position × class × candidate × **sample** × 3 genotypes, with
  no loop over nodes.
- `fit.rs:2132` (`ln`) runs per position × class × candidate × **node** (16), with no loop over
  samples.
- `ssr_fit.rs:1948` (`ln`) runs per tract × **sample × genotype** × read group × non-empty
  bucket, once per `refresh`, with no loop over quadrature points.
- `ssr_fit.rs:2073` and `:2332` run per tract × **256 points** per evaluation, with no loop
  over samples or genotypes.

The §2.8 and §2.9 tables state each of these correctly. Only the summary wording multiplies the
dimensions together, which overstates the count of the `ln` at `:2132` by a factor of the
cohort size, and the count of `ssr_fit.rs:2073` by samples × genotypes. That errs on the
cautious side, but step A2 is meant to weight functions by how often they run, so the summary
should not be the source of that weight.

**Correction.** Rewrite the bullet per call, for example: "`exp` per position × sample ×
candidate × 3 (`fit.rs:2098`); `ln` per position × candidate × 16 nodes (`fit.rs:2132`); `ln`
per tract × sample × 91 genotypes × bucket per slippage change (`ssr_fit.rs:1948`); `exp` per
tract × 256 points per evaluation (`ssr_fit.rs:2332`)".

### 4. Fix — `genotype_table.rs:535` runs per locus at tracts with more than 16 candidates

The report classes `log_factorial`'s `ln` as per-table: "built once per ploidy and allele count
and cached per thread (`:173-211`)". The cache only holds tables up to 16 alleles and ploidy 8:

```
src/calling/genotype_table.rs:57:  pub const MAX_CACHED_ALLELE_COUNT: usize = 16;
src/calling/genotype_table.rs:162-163:
    if copies > MAX_CACHED_PLOIDY || allele_count > MAX_CACHED_ALLELE_COUNT {
        return Arc::new(Self::build_uncached(ploidy, allele_count));
```

Every locus builds its table at `summarise_condition.rs:2315`
(`GenotypeTable::build(parameters.ploidy(), candidates.len())`). The repeat-tract candidate cap
ships at 32 (`calling/allele_candidates/ssr.rs:477-478`,
`MaxCandidateAlleles::new_or_panic(32)`), and the comment on the cache bound
(`genotype_table.rs:54-56`) was written against the ordinary cap of 6. So a tract with 17 to 32
candidates rebuilds its table at every locus. `build_uncached` computes one multinomial
coefficient per genotype (`:247`). For a diploid genotype that is `ln 2` once for the ploidy,
plus once more at a homozygote, so up to 528 genotypes × 2 `ln` calls per such locus.

**Correction.** Class the site "per-table, or per-locus × genotypes at repeat tracts past
16 candidates". The constant column, *small integer*, stays right, and the decision is
unchanged, because the argument is always `ln 2` at diploid.

### 5. Nit — `digamma` runs four times per Newton step, and its argument stays at or below 100

§2.8's row for `fit.rs:2950` says "up to 50 Newton steps × 3 calls". `fit.rs:2747-2748` makes
four: `digamma(a)`, `digamma(ab)` twice, and `digamma(b)`. §3.1 says *x* is "below about 108
from the second step on". After the clamp at `:2758-2759` both shapes are at most 50, so
`ab ≤ 100`. The recurrence at `:2944` only raises arguments below 8, so *x* ≤ 100.

**Correction.** Change "× 3 calls" to "× 4 calls" and "below about 108" to "at most 100".

### 6. Nit — several per-read arguments do not change inside their loops

This corrects nothing, but it bears on step A2's weighting and on how much a slower `ln` would
actually cost. Among the 12 reachable per-read calls:

- `emission.rs:558-559` compute `ln(1 − ε)` and `ln(ε/3)` from `context.substitution_rate`,
  which is fixed per (read group, candidate) at a locus. They are recomputed because
  `substitution_probability` builds a fresh `FlatEmission` on each call
  (`ssr_emission.rs:617`).
- `generic.rs:810` and `:912` take `(1 − c)·k/P + c·q`. Where a read group has no fitted
  contamination (*c* = 0), the argument is `k/P`, one of *P* values fixed for the whole run.
- `stutter.rs:783`'s `powi` takes a per-context share raised to one of at most 10 exponents.
- `ssr_unit_robust.rs:278-288` are already listed as constants.

Two hot calls outside the per-read class have the same property:

- `ssr_fit.rs:1948`: the probability depends only on (read group, genotype pair, bucket)
  (`ssr_fit.rs:1946-1947`), so a `refresh` could take at most groups × 91 × 17 distinct
  logarithms rather than one per tract × sample.
- `locus_score.rs:490`: σ is a sample's single-copy depth SD times one of 8 fixed factors, all
  set once per run.

**Correction, optional.** A column or a sentence saying "argument fixed per locus / per run"
would let step A2 separate the calls whose cost can be removed by moving them out of the loop
from those that genuinely vary per iteration.

### 7. Nit — partial reads multiply the repeat-tract emission calls by the number of length changes

§2.1 and §2.2 give `ssr_emission.rs:629` and `stutter.rs:783` as "per observation per
candidate per slip placement". For a read that ran off the tract, `censored_emission`
(`ssr_emission.rs:378`) loops over every reachable length change (`:439`) and calls
`letters_over` for each (`:448`). Its pure-tract branch calls
`probability_at_least_this_much_longer` instead (`:400-401`), which also sums the geometric. Both calls therefore also run per reachable length change at
partial observations. The class stays per-read.

**Correction.** Add "and per reachable length change for partial reads" to both rows.

### 8. Nit — one of the two examples named for `coverage_odds` does not fill it

§2.8 says `examples/ng_joint_duplicated_in_fit.rs` and `examples/ng_joint_records_walk.rs` fill
`JointFitConfig::coverage_odds`. The second passes an empty vector:
`examples/ng_joint_records_walk.rs:302: let coverage_odds = Vec::new();`. The
"not reachable" conclusion is unaffected.

### 9. Nit — the concentration floor cites the constant, not the floors

§3.1 bounds the Dirichlet concentration α at "≥ 1e-12 (`genetics.rs:107`)". Line 107 defines
`MIN_ALT_CONCENTRATION`. The floors themselves are applied where the seed is built:
`seed_generic.rs:417, 425`, `seed_ssr.rs:183` and `genetics.rs:147`. `fill_sample_concentration`
then adds a non-negative leave-one-out count (`genotype_prior/mod.rs:781`). The bound holds.

**Correction.** Cite the four floor lines.

### 10. Nit — the rescale peak at `quality/mod.rs:588` can exceed 1

§3.1 gives the argument as "(0, 1]-ish after each rescale". The logarithm is taken of the peak
*before* the rescale. That peak is a sum of up to `ploidy + 1` products of an entry at most 1
and a weight at most 1 (`quality/mod.rs:560-573`), so it lies in (0, ploidy + 1].

**Correction.** Change the range to "(0, ploidy + 1]".

### 11. Nit — the table cited for the per-quality scores is at `emission.rs:170`

The Summary cites `src/alignment/emission.rs:161-163` for the "tables built from written-out
bits". Those lines are the doc comment explaining why the bits are written out. The table,
`static PER_QUALITY_LN`, starts at `:170`.

## Scratch files

`tmp/review_A1/classify.py` (independent search and test classification),
`tmp/review_A1/mine.tsv`, `tmp/review_A1/mine_ship.tsv`, `tmp/review_A1/cite.sh` (citation
printer), `tmp/review_A1/broad.py` (search for path-form uses without a parenthesis).
