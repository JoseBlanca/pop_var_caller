# ng read likelihoods — handoff after Milestone B

*2026-08-24. Branch `ng-calling-likelihoods`, worktree `../pop_var_caller-calling-likelihoods`,
plan [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). Written for
whoever picks this up next.*

## Where the plan stands

**Milestones A and B are complete — steps A1, A2, B1, B2 of fifteen.** The next step is C1, the
contamination mixture.

| step | state |
|---|---|
| A1 evidence views | ✅ committed |
| A2 parameter views, floors, row contract | ✅ committed |
| B1 `m(a, g)`, the error spread | ✅ committed |
| B2 the closed form + production differential | ✅ committed |
| C1 the mixture, no `c = 0` branch | next |
| C2 `q(o)` — the batch allele frequency | after C1 |
| D1 partial observations | Milestone D |
| E1–E3 the stutter distribution | Milestone E |
| F1–F3, G1, H1–H2 the STR path | Milestones F to H |

Each step is one commit, with an implementation report under
`doc/devel/reports/implementations/ng_calling_likelihood_*`.

## What is built

`src/ng/calling/likelihood/`:

- **`mod.rs`** — the evidence views (`GenericObservation`, `GenericSampleEvidence`,
  `SsrSampleEvidence`), the parameter views (`ReadGroupCalibration`, `ContaminationView`), the
  floors, the buffers (`GenericEvidenceBuffer`, `SsrRowScratch<ModelScratch>`), and in its module
  documentation the row contract and the three parameter tiers.
- **`generic.rs`** — `fill_log_error_spreads` + `LogErrorSpreadTable`, and
  `genotype_log_likelihood_row`: spec §3.3's closed form, reconciled term-by-term against
  production's `standard_log_likelihood`.

Suite at B2: **4,268 passed / 0 failed / 14 ignored**, 81 of them in this module.
`cargo doc --no-deps` at `main`'s 23 unresolved links.

## The five things a successor must not rediscover

1. **The merge's allele numbering is not the candidate numbering.** `SupportedAllele.allele`
   indexes every distinct sequence the cohort showed; `AlleleId` indexes what selection kept, and a
   prune renumbers between them. The mapping is an argument
   (`GenericObservation::fill_from_supported_alleles`), never an assumption, and a dropped row's
   quality comes back rather than vanishing.

2. **The evidence view borrows the staging buffer, so the row's scratch is a separate type.** A row
   taking `&mut` the same object the evidence borrows cannot be called — `E0499`. That is why
   `GenericEvidenceBuffer` exists and why the generic row's own scratch is deferred to D1, which is
   the step that first has buffers to put in it.

3. **A partial read's witness is a set of runs with holes**, not one run. The projection restricted
   to it is a *gather*, so D1 needs a buffer sized by the widest witness; and the witness counts
   locus positions while the bases are what the read showed over them, so neither may index the
   other. Spec §5.3 said the singular until 2026-08-24.

4. **`SequencingBatches` is specified and not built.** C2 needs the batch a read group belongs to
   (`arch/parameter_prepass_joint_fit.md` §1.6). Take it as an argument, defaulting to one batch
   holding the cohort, which is that type's own stated default.

5. **The aggregation identity is not bitwise**, and neither is order independence. Both claims were
   in the specification and both are corrected, with the measurement: a relative 2 × 10⁻¹⁴ over 864
   combinations from 2 to 300 reads. **The model requirement is untouched** — no term may be a
   non-linear function of a per-read quality — and that is the half worth testing.

## What C1 and C2 have to do

**C1** — spec §3.6: `n_o · log[(1 − c)·own(o|g) + c·q(o)]`, evaluated in probability space with one
logarithm, `q(o)` taken as a parameter so this step is about the mixture and nothing else. **No
`c == 0` branch**: the two forms agree to a few ulp and the tolerance is a named test constant. The
row signature grows `contamination: &[ContaminationView]`, which A2 built and B2 deliberately did
not take.

**C2** — the frequency of the observation's own allele at this locus, over the samples in its
sequencing batch, **recomputed every iteration**. This is the item the owner moved into this plan
on 2026-08-24; the parameter pre-pass owes nothing. See (4) above.

## Two habits this run paid for repeatedly

**The review fan-out earns its cost, and the findings that mattered were not style.** Across four
steps it found: an allele numbering assumed rather than checked; a borrow shape that would have
stopped the next milestone at its first call site; a calibration laundering provenance; a term
worth 1,620 Phred deletable with every test green; and three specification claims that were false.

**Every number written about this work was measured, and several were wrong the first time.** A
"nine fields" that was eight, a factor of ten that was 3.5, a bound of 2 ulps that was 109 at
depth, "about 10⁻¹⁵ nats" that was fourteen times larger, an inequality pointing the wrong way. The
rule that catches these is asserting a number in a test rather than stating it in prose.

## Open items for the owner

Nothing blocks C1. Two things are recorded and unresolved:

- **The divisor buffer's owner** should be an eighth field of `CallingScratch`
  (`arch/calling_em_loop.md` §2) — per locus, not per sample, refilled once per locus and not once
  per pass, generic-path only. Recorded in `arch/read_likelihoods.md` §3; the loop plan has not
  been edited.
- **`SsrContamination`'s `OPEN:`** — the prior's seed builder fills one entry per *candidate* and
  the mixture wants a probability per observed *length*. Milestone H2 meets it; nothing before then
  needs it.
