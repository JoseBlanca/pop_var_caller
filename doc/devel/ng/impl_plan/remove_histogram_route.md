# Removing the whole-genome histogram route

**Owner's decision, 2026-09-11.** ng keeps one route to its parameters — the census, fitted over
every sample at the same positions — and the per-sample whole-genome histogram route goes, the
runs-of-homozygosity estimator of the inbreeding coefficient with it.

## Why, in one paragraph

No shipped command has ever run the histogram route: `generate-psps`, `estimate-parameters`,
`call-from-psps` and `regenerate-census` all read censuses, and the join that would hand histogram
results to the caller (`RunParameters::from_prepass`) is called only by its own unit tests. The one
parameter the census could not produce was the inbreeding coefficient, because that estimator reads
*which windows of the genome lie in a run of homozygosity* and a census window holds a quarter of
one heterozygote at the shipped budget. Raising the budget to six million was measured and works —
worst error 0.086 against a realised 0.78 — but the coefficient a caller reads is the per-sample
departure from Hardy–Weinberg, which the cohort fit already measures from two million positions
because it too is an average over them
([report](../../reports/ng_census_inbreeding_budget_2026-09-11.md)). So the runs estimator is a
second opinion rather than the prior, and it is not worth 35,000 lines.

**What is given up, stated plainly:** ng will no longer be able to estimate autozygosity from the
genomic *distribution* of heterozygosity. `parameter_prepass_generic.md` §6.3 calls that the
non-circular estimator and prefers it precisely because the homozygote excess is measured against a
population expectation the same fit produced. **The owner has accepted that trade** — the caller
does not need it — and the circularity now stands uncorrected and must stay stated in the output.

## The keep-list, which is the whole difficulty

Four things inside the doomed trees have live consumers on the census route and in calling. They are
rehomed first, so that deleting the trees afterwards is a deletion and not a rewrite.

| what | why it survives | where it goes |
|---|---|---|
| `generic/depth_bins.rs` | `DepthBinEdges`, `DepthBin`, `CENSUS_DEPTH_BIN_COUNT` — the census's own depth ladder, read by `joint::{census, census_file, fit, contamination}`, `run::gatherer`, `calling::parameters_file::bindings`, `estimate-contamination` | `parameter_estimation/depth_bins.rs` |
| `generic/calibration.rs` | `MintedReadErrors`, `minted_error_by_read_group` — what a library's own base qualities claimed, accumulated by the census writer and read by `calling::likelihood` | `parameter_estimation/calibration.rs` |
| `generic::DEFAULT_ERROR_RATE` | the rate calling falls back to | `parameter_estimation/mod.rs` |
| `ssr::{RepeatCount, Stratum, StratumKey}` | how calling names a repeat-tract stratum; `joint::census::Stratum` is a *different* type (plain numbers) and `run::census_fit` converts between them | `parameter_estimation/repeat_strata.rs` |

Everything else in `generic/`, `ssr/` and `fitting/` is reached only from inside those trees, from
their own tests, or from examples.

## Stages

Each stage compiles and its tests pass before the next begins.

- **0 — the budget back to two million.** `CensusSelection::SHIPPED.generic_target`, and the two
  command doc comments that quote it. ✅
- ✅ **1 — rehome the four survivors** and repoint every importer. No behaviour changes; the test
  suite must stay green on the same count.
- ✅ **2 — delete the trees.** `generic/`, `ssr/`, `fitting/`, `parameter_estimation/subsample.rs`
  (reached only from the two deleted paths), and `joint/census_runs.rs` — the census-to-runs-model
  bridge built on 2026-09-11, which has no purpose once the estimator it feeds is gone.
- ✅ **3 — delete the route's join into calling.** `RunParameters::from_prepass` and its tests; the
  `ParameterEstimationError` variants only the deleted fits raise.
- ✅ **4 — delete the runs-of-homozygosity *source*** from the census moments:
  `InbreedingSource::RunsOfHomozygosity`, `PanelInbreeding::FittedFromRuns`, the below-the-floor
  warning, and their tests. What is left is `Supplied` and `HomozygoteExcess`, which is the world
  after this change.
- ✅ **5 — delete the examples** that drive the removed code.
- ✅ **6 — the docs.** Every dangling intra-doc link (the pre-commit gate runs `cargo doc`), and a
  note at the head of each spec section describing a route that no longer exists.
- ✅ **7 — verify.** `cargo check`, `cargo test`, `cargo clippy`, `cargo doc`, and the oracle below.

## The oracle

**The parameters a cohort fits must not change.** `tmp/params_2M.toml` was written before any of
this, by the two-million binary over the four tomato psps in `tmp/psp_2M`. After the removal the
same command over the same psps must produce the same fitted numbers — the removal touches no code
the census route executes, so anything that moves is a mistake this catches.

## What it came to

**Done 2026-09-11.** All seven stages, and the oracle holds: `estimate-parameters` over the four
tomato psps in `tmp/psp_2M` writes a parameters file **byte-identical** to the one the
pre-removal binary wrote — md5 `bcc95d64d20d3edc7c6ee9c4ccc09933`, 96,373 bytes both times.

| | before | after |
|---|---:|---:|
| lines under `src/ng/parameter_estimation/` | 61,536 | 26,851 |
| library tests | 6,748 | 6,200 |

**34,685 lines and 548 tests went**, plus four examples and `RunParameters::from_prepass`. `cargo check --lib --tests`, `cargo test --lib`,
`cargo clippy --lib --tests` and `cargo doc --no-deps --lib` are all clean. Twelve intra-doc links
into the deleted modules were repaired, and five design documents carry a banner saying what in them
now describes code that does not exist.

**Two examples still do not compile, and neither is this change's**:
`ng_cohort_merge_real_cost` and `ng_candidate_selection_probe` are drifted against
`cohort_merge`'s API and contain no reference to anything removed here. PROJECT_STATUS records that
breakage from 2026-09-08.

## The follow-up, done the same day

**`estimate-parameters` now writes the coefficient it fitted, and the parameters file is the only
way one reaches a calling run.** Three rungs, resolved once and carried by the warrant on each row
(`DeclaredInbreeding::of_each_sample_over`): `supplied` where `--inbreeding` was given — which
overrides the fit, per the owner's ruling of 2026-08-27 — `fitted_here` where the cohort's own fit
measured each sample's homozygote excess, and `defaulted` at zero where neither. `--inbreeding`
became an option rather than a default of zero, because *"these plants are not inbred"* and *"use
what you measured"* are different instructions.

**Two things it will not do.** A sample whose fitted excess is exactly one has no coefficient to
become — the excess is a fraction in `[0, 1]` and a coefficient one in `[0, 1)` — so it is left out
and falls to the rung below rather than being coerced to something indistinguishable from one. And
a single-sample run gets `defaulted` even though its fit ran, because one genome's totals cannot
identify an excess: it comes back zero whatever the truth, and the fit already marks it so.

**Measured on real reads**, the four tomato accessions over 8 Mb: the fit writes 0.9128, 0.9836,
0.9822 and 0.9464, `fitted_here`, each over about 1.9 million covered positions, and calling with
them rather than with zero changes **4,497 genotypes of 191,752 (2.35%)** and cuts the heterozygous
call rate from 6.77% to 2.61%. The two psps that are the same plant sequenced twice get 0.9836 and
0.9822 — a difference of 0.0015, which is the consistency check this plumbing can be given.

**⚠ Those coefficients are a four-sample artefact and must not be read as tomato's biology.** The
63-accession fit put the panel at 0.23 to 0.90, median 0.78; at four samples the allele-frequency
curve is barely constrained and the homozygote excess absorbs what the curve cannot hold, which is
how it lands near the ceiling. **What is verified here is the plumbing, not the number** — and the
effect on calls at a coefficient this high is correspondingly larger than the panel would give.
Whether any of it is an *improvement* is still not established, for the reason §7 of the budget
report gives: the replicate pair's agreement rises from 95.04% to 98.07%, but "both heterozygous"
collapses from 1,230 pairs to 172, so most of that gain is a prior that has nearly forbidden
heterozygotes rather than one that has found the right ones.
