# ng step 4, the SNP/indel path — Milestone E: the four fits

**Date:** 2026-08-09. **Branch:** `ng-parameter-estimation`, `543cc9fe` → `4ea17d46`, nine
commits. **Plan:** [`parameter_prepass_generic.md`](../../ng/impl_plan/parameter_prepass_generic.md),
steps E1–E4, all ✅. **At Checkpoint E.**

Milestone E is where step 4 stops being machinery and starts producing numbers. Milestones
A–D built the vocabulary, the cell table, the path from a locus to a cell, and the fitting
mathematics — and proved all of it before a single locus was read. E turns that into four
fitted parameters: a per-read-group error rate, the sample's genotype frequencies at each
ploidy, its inbreeding coefficient, and the provenance of each.

**Nothing here reads a locus either.** Every fit is proven against a table filled cell by
cell from a known truth, which is this plan's own rule: an accumulator bug and a fit bug
cannot then hide each other. The locus stream reaches a fit in Milestone F.

---

## What was built

### E1 — the per-read-group error rate (`25ee9a1c`, fixes `7f730da6`)

**E1 is not the standalone profile fit the plan and arch §5.2 describe.** The owner settled
that on 2026-08-07, from the rule *if we measured something and it worked OK, build that*:
both arms of `examples/ng_multilib_key_harness.rs` fit `ε` with
`fit_eps_on_read_group(space, freqs)`, which scores each candidate rate at the frequencies it
is handed and **never re-climbs them**. So E1 is the `ε` half of an alternation.

The profile scan gained a **sibling**, not a mode argument — a scan at fixed frequencies is
not a profile likelihood, and a function named for one that sometimes is not one is a name
that has stopped being true.

- `fitting/profile_scan.rs` → **`fitting/ladder_scan.rs`**, for the same reason: it now holds
  two scans and its name said one.
- **`scan_ladder`** is the rung loop both run — the per-ploidy plan, the declared-width
  check, the value checks that name the rung, the tie rule and the rail flag. Only the
  scoring step differs, and it is a closure. `ploidy_plans` was later split out of it so each
  scan can check its own arguments once, before the first rung, rather than at each of 161.
- **`fit_by_fixed_frequency_scan`** scores every rung at frequencies handed in, one set per
  ploidy, climbing nothing. `weighted_log_likelihood` became visible outside its file: the
  scan is that function once per rung.
- **`FixedFrequencyScanResult`** rather than `ScanResult`, resolving the wrinkle the owner
  left open. It carries the winning rung's **index** — the coupled fit's stopping rule is
  stated in rungs — its parameters, its score and the rail flag, and **no genotype
  frequencies**: the only ones it could report are the caller's own argument echoed back.
- **`generic/read_group_error_rate.rs`** — `fit_read_group_error_rates`, one rate per group
  from that group's own table, gathering every ploidy the group covered into one scan.

### E2 — the coupled loop (`2d9b8095`, fixes `12203d92`, `08ec113c`) — *own commit*

One iteration is: the genotype frequencies, climbed on the whole-sample table at the rates
the previous round produced; then each read group's rate, scanned on its own table **at those
frequencies and without re-climbing them**. It stops when every read group's winning rung is
the one it had last round.

- `fit_coupled` takes the accumulator; `fit_coupled_from_tables` takes the two tables and is
  what the tests drive; `fit_by_alternation` takes the start and the cap as parameters, so
  what they cost is a test rather than a recompile.
- An iterate's score is the whole-sample table's likelihood **at that iterate's own rates and
  frequencies** — one objective on one table, which is what makes "best-scoring" a defined
  comparison. Neither block's own score is: step 1's belongs to the previous rates and
  step 2's to a different table.
- `DepthAltHistogram::total_reads` — a library's share `w_g` is its reads over the sample's,
  **not its sites over the sample's**.
- `fit_mixture_weights` became `pub(crate)` at last, its consumer being step 1 of this loop;
  its `# Examples` doctest became a named test rather than being deleted.
- `generic/expected_counts.rs` — the infinite-genome table generator, moved out of E1's test
  module and shared.

### E3 — the runs model (`5b9b2d5c`, fixes `72787673`) — *own commit*

A two-state hidden Markov model over 100 kb windows. Forward–backward in log space,
restarting at every contig boundary; each state's three genotype frequencies re-estimated
freely from its own window posteriors; both transition rates fitted; the ordering constraint
applied by relabelling at the end.

- `F` is the coverage-weighted posterior occupancy, weighted by `total_covered_positions()`
  and **not by loci**.
- The emission is a sum over the window's cells — exactly `weighted_log_likelihood` at a
  state's frequencies, so the two cannot drift.
- **Absent windows are in the chain**, from window zero to the last that received a site on
  each contig. The floor is checked against the windows that **hold sites**, so a
  region-restricted run cannot pass it on padding.
- `resolution_at` interpolates research note §3.6's eight-seed means in log–log.
- `InbreedingStatesNotSeparated` rather than a zero when no start left posterior mass on both
  states.

### E4 — the fallback ladder (`e1ad9290`, fixes `4ea17d46`)

`resolve_error_rates` walks the ladder per read group — **fitted here → borrowed → supplied →
defaulted** — with `MIN_SITES_TO_FIT` gating the first rung. A borrowed value is the plain
mean of the qualifying groups' rates; its `observations` is the **lenders'** reads, because
this group's own reads are precisely the ones that were not enough.
`take_supplied_inbreeding` is separate because a supplied `F` changes what the accumulator
*stores*. Wired into `into_coupled_fit`, so `CoupledFit` stops marking a rate fitted from 500
sites as `FittedHere`.

---

## What it is proven against

| step | the claim | how |
|---|---|---|
| E1 | two libraries four Phred apart each recover their own rate | rungs 80 and 64 of 161; answering either from the other's table lands sixteen rungs out, and pooling lands on rung 70, asserted |
| E1 | the handed frequencies decide the answer | at **six** reads a site the fit moves 80 → 84 and stays interior |
| E2 | from three times the true rates the fixed point is the truth | rungs 61 and 45 in, 80 and 64 out, frequencies inside 1%, converged in two rounds |
| E2 | two different starts reach the same fixed point | the default rung and three-times-the-truth agree |
| E2 | at one read group the alternation reaches the profile scan's answer | asserted on the rate and on the frequencies |
| E2 | the whole-sample table's frequencies reach the fitted rates | two samples with **byte-identical** read-group tables and whole-sample tables that disagree about the individual return different rates |
| E3 | `F` recovers a drawn genome's **realised** autozygous fraction | within 0.02 at four nominal levels from 0.05 to 0.60, where the realised fractions miss their nominal values by up to 0.027 |
| E3 | a floor of false heterozygotes at five times the real rate does not move it | four levels, each scored against its own realised fraction |
| E3 | starts sharing one separation return a failed search, not a zero | `InbreedingStatesNotSeparated` where the default nine land within 0.05 |
| E4 | a thin library comes out `Borrowed`, a supplied one `Supplied`, the rest `Defaulted` | at the fit a caller reads, not only in the ladder's own unit test |

**Suite 3,108 → 3,187** (+79), 0 failed, 5 ignored. **`ng::parameter_estimation` 206 → 285**
(+79). Doctests in the module 1 → 0, by design. `cargo fmt --check`, `cargo clippy
--all-targets --all-features -- -D warnings` and `cargo doc --no-deps --lib` at its
12-unresolved-link pre-existing baseline, none in this module, before every commit.

Both research harnesses were re-run to completion at the start: the multi-library one's
coupled fixed point is the truth in every world, and the inbreeding one's `F` posterior is
0.3073 against a realised 0.3073 at true 0.30, and 0.2634 against 0.2629 at the 3×
false-heterozygote floor — the two figures the design docs quote.

---

## Recorded deviations from the design

All licensed by the architecture's "signatures are illustrative" preamble; none is a
deviation from a measurement.

1. **`fit_coupled(sample, accumulators, ladder, supplied)`** — arch §5.2 gives one parameter.
   The sample name is what the error variants require; the ladder is taken so a caller
   alternating twenty times builds it once; `supplied` is E4's.
2. **`fit_coupled_from_tables`** exists beside it, taking the two tables, so the alternation
   can be proven without a locus stream — this plan's own rule. `fit_coupled` has its own
   test even so, because gutting it left all 3,147 tests green.
3. **`fit_inbreeding(sample, windows, noise, ploidy, outside_rates, starts)`** — arch §5.3
   gives `(windows, error_rate, ploidy)`. It takes a `SampleLibraryNoise` rather than a bare
   rate map, because the shares travel with the rates and the multi-library scoring rule needs
   both; and the sample's own `SampleRates`, because that is where each start's outside state
   begins.
4. **`RunsModelStarts` cannot express what the harness's starts actually vary.** Its three
   `(enter, leave)` pairs imply mean run lengths of 10.5, 50 and 83 windows as well as three
   inside fractions; the type carries only the fractions. The rates are derived here from
   `implied_f` and a stated `MEAN_RUN_WINDOWS_AT_START = 20`. Both rates are fitted, so this
   changes how many passes a start takes and not what it converges to — confirmed by moving
   the constant to 5 and to 200 with every test green.
5. **The runs model fits per-*window* transition rates and converts to per-base for
   reporting**, where arch §5.3 says they are fitted per base. The harness fits per window
   too; it is the arch text that is now wrong about the direction.
6. **E3's tolerance is 0.02, not the plan's four decimal places.** Four decimals come from
   8,004 windows of 100,000 sites; 3,600 windows of 400 sites hold about a fifth as many
   heterozygotes each. What is checked is that the estimator recovers a **drawn** fraction
   across a five-fold range of it.
7. **`fitting/profile_scan.rs` → `fitting/ladder_scan.rs`**, and `ReadGroupErrorRate` →
   `ReadGroupErrorRateFit`, `alternate` → `fit_by_alternation`, `Iterate` → `ScoredIterate`,
   `assemble` → `into_coupled_fit`, `rung_nearest` → `nearest_rung`, `kept` →
   `scoring_output`, `supplied_inbreeding` → `take_supplied_inbreeding`.

---

## What Milestone E does **not** establish

- **Nothing has read a locus.** F1 wires the stream; F2 refits from a directly-filled
  accumulator at ploidy 2 and 4 and at 3 reads and 300×; F3 runs the identities on both
  cohorts.
- **The `Supplied` rung has no source.** `GenericEstimationConfig` has no supplied-rates
  field and does not exist yet (F1). E4 added the parameter; nothing can fill it.
- **`take_supplied_inbreeding` has no production caller** — F1's.
- **Keep-the-best-scoring-iterate is separated from keep-the-last on one world only**, and
  keep-the-best-scoring-*start* is separated from keep-the-first on none: on every drawn
  genome all nine starts land on the same answer with the same score. The rule is defensive
  rather than demonstrated, and the comment claiming a collapsed fit can outscore a real one
  was wrong — measured, collapsed starts score 2,165 nats **worse**.
- **The contig restarts in the runs model change no answer** on any fixture, and the reviewer
  proved why: the cross-boundary transition term is normalised by the wrong contig's total
  and underflows to exactly zero, and an absent window has zero covered positions so it cannot
  move `F` whatever the chain believes about it.
- **Nothing asserts a log-likelihood *value* in the runs model**, so summing it at every
  window instead of once per contig — a 300-fold error — survives.
