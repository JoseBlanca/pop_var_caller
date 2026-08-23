# Code Review: ng_calling_prior_a1
**Date:** 2026-08-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A1 of the genotype-prior plan — the `genotype_prior/` folder and `types.rs`'s new diversity scalar
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step A1 of
  [`calling_prior.md`](../../ng/impl_plan/calling_prior.md), branch `ng-calling-prior`.
- **Reviewed against:** base commit `1742e3d6` plus `tmp/a1.patch`. Every sub-agent re-pointed its
  own worktree at that commit and applied the patch, and each confirmed
  `src/ng/calling/genotype_table.rs` and `src/ng/calling/genotype_prior/mod.rs` present before
  reviewing.
- **In-scope files:**
  - [src/ng/types.rs](../../../../src/ng/types.rs)
  - [src/ng/calling/genotype_prior/mod.rs](../../../../src/ng/calling/genotype_prior/mod.rs)
  - [src/ng/calling/genotype_prior/dirichlet_multinomial.rs](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs)
  - [src/ng/calling/genotype_prior/plug_in.rs](../../../../src/ng/calling/genotype_prior/plug_in.rs)
  - [src/ng/calling/genotype_prior/seed_spectrum.rs](../../../../src/ng/calling/genotype_prior/seed_spectrum.rs)
  - [src/ng/calling/genotype_prior/seed_ssr.rs](../../../../src/ng/calling/genotype_prior/seed_ssr.rs)
  - [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs)
- **Deliberately out of scope:** frozen production (`src/var_calling/`, `src/ssr/`,
  `src/genetics.rs` — read as reference only); the later plan steps that fill the four empty files;
  the three aggregate gates already red on `main` before this branch existed.
- **Categories dispatched, and why these five.** `reliability`, `errors` and `naming` are
  always-on. `defaults` because the diff lands a default-acting value on a scientific result.
  `module_structure` because the diff creates a folder and places a shared scalar.
  `idiomatic`, `refactor_safety` and `smells` were **not** dispatched: the diff declares one
  newtype, one constant and one enum variant, all following a shape the file repeats five times,
  and the remaining 60% of it is prose — the three categories would have re-covered `naming`'s
  and `module_structure`'s ground. `unsafe_concurrency` and `tooling` have no trigger here.

### 2. Verdict

**Approve-with-changes.** No Blocker. Four Major findings, all of them about what the tests fail
to pin rather than about the code's behaviour, and three of the four were demonstrated by a
mutation that survived the submitted suite. Two Minor findings correct claims the added prose
makes about other code — one of them a mechanism that never happened.

### 3. Execution status

Run in the dev container from the branch worktree's own `scripts/dev.sh`, and passed verbatim to
every sub-agent so none re-ran them:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 13.35s` |
| `cargo test --lib` | 0 | `test result: ok. 4005 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 639.43s` |

Not run, and why: `cargo clippy --all-targets`, `cargo test --all-targets` and
`cargo doc --no-deps --lib` are red on `main` in files this branch does not touch (18 clippy errors
across two benches and one example; a panic at `benches/psp_writer_perf.rs:386`; 17 unresolved
intra-doc links). `cargo audit` was not run — no dependency changed.

Findings labelled "Needs verification": **0.** Every finding below carries either a mutation that
was run or a file:line that was read.

**Mutation totals across the five agents: 27 run, 14 survived, 1 changed no behaviour.** The one
that changed nothing was `get()` returning `self.0.abs()`, which only clears the sign bit of
`-0.0`; it is reported because it is the boundary of what the type promises, not as a hazard.
Every other survivor was proved to change behaviour with a probe before being written up.

### 4. Open questions and assumptions

1. **Do later steps of this plan deliver the diversity as `Estimate<ExpectedHeterozygosity>` with
   a `Provenance`, or as a bare value?** The `errors` and `defaults` agents both raised the
   fallback's provenance; the answer decides whether M4's second half is A2's work or D2's.
   Affects **M3**, **M4**.
2. **May a name the architecture writes down be changed by the implementation?** Three findings
   propose renames of items `arch/calling_priors.md` §2.1 and the plan's file list spell out.
   Affects **Mi1**, **Mi7**, **Mi8**. Raised to the owner at Checkpoint A rather than decided here.

### 5. Top 3 priorities

1. **M1** — the new scalar is missing from the one test that pins the accept/reject boundary
   exactly, and two wrong constructors survive because of it.
2. **M3** — the fallback constant's *value* is unpinned, and the test that looks like it pins it
   passes at `0.1` and at `1.0`.
3. **M2** — a bare `f64` constant seeds three sibling newtypes silently, in the module the STR
   path imports from, under a doc comment warning against exactly that.

### 6. Findings

#### Major

**M1: src/ng/types.rs:1349 — `ExpectedHeterozygosity` is missing from the boundary sweep**
**Categories:** reliability, errors (convergent). **Confidence:** High.
`the_constrained_rates_accept_exactly_the_probabilities_and_round_trip` sweeps the whole `f64` line
for the other three constrained rates, asserting acceptance holds *exactly when* the value is a
finite probability and that an accepted value comes back bit for bit. The diff added the fourth
rate to four point tests and left this one at three. Two mutations survived the submitted suite as
a result: a constructor widened to `(-0.25..=1.25)` (accepts 1.1) and one quantising to six
decimals (turns `0.00062831853` into `0.000628`). `θ` reaches the prior as `α_alt`, so a drifted
bound or a rounding constructor is a silently wrong prior rather than an error.
**Fix:** one line in the existing arm, plus the count in its doc.

**M2: src/ng/types.rs:603 — the fallback is a bare `f64`, outside the newtype discipline the
section exists to enforce**
**Categories:** defaults, errors (convergent). **Confidence:** High.
The section header states the rule — five types and not one shared `Probability`, so an inbreeding
coefficient cannot be handed to something expecting an error rate. `DEFAULT_SPECIES_DIVERSITY_FALLBACK`
is declared outside it, and the `defaults` agent compiled the consequence: the constant seeds
`InbreedingF`, `ErrorRate` and `GenotypeFrequency` alike, because `1e-3` is a legal value of each.
The second half is worse than the first: the type's own doc three lines above warns that the STR
path must never take a SNP-scale diversity, and A1 then places the only diversity constant in the
shared vocabulary — the module the STR path imports from as a matter of course.
**Fix:** make it a value of the type it defaults, following `AlleleId::REFERENCE`'s precedent in
the same file.

**M3: src/ng/types.rs:1287 — the fallback's value is unpinned, and the test claims a guard it does
not provide**
**Categories:** reliability, errors, defaults (convergent). **Confidence:** High.
`the_species_diversity_fallback_is_a_constructible_heterozygosity` asserts
`try_new(CONST).unwrap().get() == CONST` — the constant on both sides — so it pins membership of
`[0, 1]` and nothing about the value. Its doc claims it would catch "a percentage, say, or a
per-kilobase rate"; both land *inside* `[0, 1]` and both were run as mutations and survived:
`1e-1` (the percentage reading), `1.0` (the per-kilobase reading), and `5e-4` and `0.5` besides.
This is the value a run with no fitted `θ` calls genotypes with, and it is not loud when wrong —
at `θ = 0.1` the hom-ref prior weight falls from 0.9985 to 0.85 and every such run tips toward
variant calls with nothing failing.
**Fix:** assert the value, and keep a second assertion that the value satisfies the type's own
predicate.

**M4: src/ng/types.rs:591–603 — the value is exported wider than the machinery that must report
it**
**Category:** defaults. **Confidence:** High.
The constant's doc states the contract — a run that lands on it must carry that fact into its
output — and A1 ships no shape that can carry it. `SeedRegime::FallbackDiversity` is scheduled for
`calling/genotype_prior/mod.rs` at A2, one module down, while the constant is `pub` from
`crate::ng::types`: a consumer outside `calling::genotype_prior` can import the value without ever
importing the enum that reports it. Against production, ng is currently weaker on both halves —
production makes its fallback overridable by passing it as a parameter and records
`DiversitySource::{Estimated, CliOverride, PriorFallback}` on all three branches. Two caveats the
agent checked rather than assumed: production's `cli_override` is a door with nothing behind it
(`pipeline.rs` passes `None`; no flag exists), and its `.source` field has no consumer outside its
own tests — so on provenance actually reaching output, the two are equal and both open.
**Fix, in two parts:** the doc obligation now, the structural one at A2 or D2 (see open question 1).

#### Minor

**Mi1: src/ng/types.rs:603 — `DEFAULT_SPECIES_DIVERSITY_FALLBACK` says one idea twice and names the
wrong quantity.** `DEFAULT_` and `_FALLBACK` both mean "used when nothing was supplied", and the
doc two lines above reads "not a default for the estimate" — the name asserts what its own doc
denies. `DIVERSITY` re-imports the ambiguity the type's doc spends a paragraph warning against.
The agent compiled and ran `SPECIES_FALLBACK_EXPECTED_HETEROZYGOSITY` as a replacement.

**Mi2: src/ng/types.rs:565–568 — "that substitution is the mistake production's STR path made"
names a mechanism production did not use.** Production's STR path hardcodes
`SFS_THETA = 0.01` (`src/ssr/cohort/freebayes_emit.rs:42`), freebayes' default, marked "Fixed, not
a per-run knob". That is a **population-scaled** θ for an Ewens prior — a different quantity in
different units from an expected heterozygosity in `[0, 1]`, and a different number from this
file's `1e-3`. Nothing was substituted from the SNP path. The spec states it correctly (§5); the
doc comment compressed it into a mechanism that never occurred.

**Mi3: src/ng/types.rs:599 — "a tomato panel is more diverse than a human one" points the opposite
way to this project's own fitted tomato number, and carries no size.** Spec §4.1 records tomato1's
fitted diversity as **6 in 10,000** — 0.0006, *below* the 1e-3 fallback. The spec hedges the claim
("Soft:"); the doc comment drops the hedge and states it flat. A reader who takes it at face value
overrides the fallback upward for tomato, which is the wrong direction for the one tomato cohort
this project has measured.

**Mi4: src/ng/types.rs:1266–1278 — the `NaN`/infinity assertions check only that *an* error came
back.** The loop uses `is_err()`, dropping the variant check that the type's own doc recommends
(`matches!(err, Err(DomainError::ErrorRate(r)) if r.is_nan())`). Returning
`DomainError::GenotypeFrequency` for an infinity passed all 48 tests, and a run then reports
"genotype frequency inf is not a finite probability in [0, 1]" — the wrong-fit misdirection the
separate variant exists to prevent, on exactly the path a bad fitted diversity arrives by.

**Mi5: src/ng/types.rs:689–697 — the reason given for a separate `DomainError` variant does not
distinguish it from `GenotypeFrequency`.** The doc argues "a heterozygosity averaged over sites and
the share of sites carrying one genotype are different quantities" — but ng's *observed*
heterozygosity is both, and is literally returned as a `GenotypeFrequency`
(`parameter_estimation::generic`'s `observed_heterozygosity`). The real distinction is whose two
chromosomes: this type draws them from the cohort, that one from one individual — which is why
inbreeding drives the second down and not the first.

**Mi6: src/ng/types.rs:570 — "Source:" names one of the two routes the architecture names.** Arch
§2.1 says "`JointFit::expected_heterozygosity`, **or the histogram route's mean**". As written, a
reader on the generic route concludes no fitted value reaches them.

**Mi7: seed_spectrum.rs / seed_ssr.rs — the two filenames are not parallel.** `seed_ssr` names the
locus class it serves; `seed_spectrum` names the input it reads, so a reader has to already know
that only the SNP/indel path has a spectrum. The crate's settled word for the other branch is
`generic` (`parameter_estimation/generic`, `locus_generation/generic`).

**Mi8: plug_in.rs — the module name reads as "plugin / extension point" before it reads as the
statistical term.** The concept the module owns has a name the geneticist reader already knows and
which the module's own first line uses: Hardy–Weinberg. The plug-in character survives in the type
the module will hold, `PlugInWrightPrior`.

**Mi9: src/ng/types.rs:1190 — the rewritten rationale entrenches an `InbreedingF` behaviour the
design says is wrong.** The doc now argues that an inbreeding coefficient of exactly one is "a real
answer, so a half-open check would reject valid data". Arch §2.1 says the opposite in as many
words, and says this exact assertion "moves to the rejection list". The code is right — the plan
defers the tightening to `calling_prerequisites.md` Milestone A — but whoever performs that
tightening meets a test whose comment argues the change is a bug.

**Mi10: src/ng/types.rs:649, :654, :1210 — three statements in the error type's own contract prose
are now stale.** The enum doc forward-references `Theta` as a future arrival —
`ExpectedHeterozygosity` *is* θ; "the three rate constructors" is now four; and the rejection
test's doc says "for all three" and omits the new type from the share-or-split list, which is the
one doc in the file the next contributor reads as precedent.

**Mi11: src/ng/calling/genotype_prior/mod.rs:27 — "softmax" is undefined machine-learning jargon**
in the sentence that explains what the module's output *is*. The spec's own wording for the same
step is "the caller multiplies that by what the reads say and normalises".

**Mi12: mod.rs:38–39, dirichlet_multinomial.rs — "Dirichlet-multinomial" and "marginalized" carry
the argument without ever being defined**, in the folder where they name the module and name the
shipping default. The definitions exist and are short (spec §1, §2.1), so this is a relocation.

**Mi13: plug_in.rs:7–9 — "largest at one sample and low depth" asserts a size the project has
measured.** Spec §2.2: on GIAB, each sample called on its own at 5×, 83.6% → 94.6% genotype
accuracy at true variants and 214 → 8 sites where a two-copy carrier was called heterozygous.

**Mi14: mod.rs:25, 28 — "Cheap arithmetic" and "the expensive call" state a cost with no measure.**
Arch §3.2 gives it: one `lgamma` per (allele, non-zero count) pair plus one `logsumexp` per
homozygous genotype, against one addition per allele for the other function.

#### Nits

`mod.rs:43` says "panel" where the same paragraph says "cohort" for the same samples. `mod.rs:39`
says "(plan step B)" where the file it points at splits the work into B1 and B2.
`seed_spectrum.rs:1` opens on "concentration", defined only on the parent module's rustdoc page.
`seed_ssr.rs` uses `ssr` and "STR" in one line without saying once that this is the crate's
convention. `src/ng/types.rs:472` counts the shared predicate's users at six and the crate's
probabilities at seven; both are now one higher. `src/ng/types.rs:576` accepts `-0.0` and returns
it with its sign bit intact, where `Phred` in the same file normalises deliberately and proves it —
the choice is untested and undocumented either way. `src/ng/types.rs:444` lost both halves of the
placement rationale ("Neither step exists yet", and the clause naming what the placement prevents),
leaving the paragraph reading as though two live consumers exist; there are none.
`src/ng/mod.rs:15` still says `calling` holds "so far the vocabulary its four sub-modules will
share", while the in-scope `calling/mod.rs` was updated to record that one sub-module now exists.

### 7. Out of scope observations

- `ExpectedHeterozygosity` admits up to `1.0`, while an expected heterozygosity at a **biallelic**
  site cannot exceed 0.5. Arch §2.1 specifies `[0, 1]`, so the code matches the design; whether the
  bound should be tighter is a question for the design, and it is not biallelic-only — a
  multi-allelic site's gene diversity can approach 1.
- `arch/module_layout.md`'s tree still lists `GenotypeTable`, `GenotypeIdx`, `AlleleId` and
  `Genotype` as contents of `calling/mod.rs`; all four have moved. Pre-existing, already recorded,
  and A1 does not make it worse.
- Production's `DiversitySource` field is documented as "recorded in the VCF header downstream" and
  has no consumer outside its own tests. Frozen production; noted because ng is about to build the
  same thing and should not inherit the gap.

### 8. Missing tests to add now

Grouped by what they guard. All four were checked by running the mutation they are written for.

1. **`the_constrained_rates_accept_exactly_the_probabilities_and_round_trip`** (extend, not add) —
   input class: the whole `f64` line, plus the dense `-2.0..3.0` arm either side of the bounds.
   Catches a moved bound and a clamping or rounding constructor.
2. **`the_species_fallback_is_one_difference_per_thousand_bases`** — input class: the constant's
   value. Catches the percentage slip (`0.1`), the per-kilobase slip (`1.0`) and any silent drift.
3. **`expected_heterozygosity_rejection_names_its_own_quantity`** — input class: the rendered
   `Display` of a rejection at `1.5`. Catches a message naming a different fit.
4. **`expected_heterozygosity_carries_negative_zero_verbatim`** — input class: `-0.0`, reachable
   from the fitted density, whose segregating mass is a `.max(0.0)`. Catches an unnoticed
   normalisation in the accessor, and records which of the two spellings the type promises.

### 9. What's good

- The `DomainError` variant is placed so the enum's variant order mirrors the newtype declaration
  order exactly — a reader scanning either list finds the same sequence (`errors`).
- Every added `NaN` assertion uses `is_err()`/`matches!` rather than `assert_eq!`, respecting the
  IEEE-equality trap the enum's own doc warns about (`errors`).
- The five files are exactly the five the architecture's tree names, and each stub's doc describes
  the content that tree assigns that filename (`module_structure`).
- Each stub states in bold which plan step fills it, so an empty module cannot quietly outlive its
  reason (`module_structure`).
- The type name matches the pre-pass field the value comes from, so the same concept has one
  spelling across the crate (`naming`).

### 10. Commands to re-verify

- `<worktree>/scripts/dev.sh cargo fmt --check`
- `<worktree>/scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `<worktree>/scripts/dev.sh cargo test --lib`
- New, and fast: `<worktree>/scripts/dev.sh cargo test --lib ng::types::tests`

Per-category files are left as an audit trail in
`tmp/review_2026-08-21_ng-calling-prior-a1/`.
