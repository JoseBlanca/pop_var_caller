# ng parameters file — D3: the fourth binding demotes rather than refusing

**Date:** 2026-08-31
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone D, step D3
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §2, §2.1, §6, §13 test 5
**Code:** `census_disagreement`, `demoted_to_no_better_than_supplied`, `to_run_parameters_for` and
`ParametersForThisRun`, in
[bindings.rs](../../../../src/ng/calling/parameters_file/bindings.rs)

---

## 1. What it does

Spec §6's fourth binding is the census — which store of evidence the numbers were fitted from — and
it is the one that does not refuse. Two censuses of one cohort differ in which loci were kept or at
what depth, not in what a plant's genome is, so the numbers remain interpretable: §2.1 keeps them
and demotes the file.

`to_run_parameters_for(reference, read_groups, census)` is the door a run calls, and it is the whole
of §6 in one place: `validate`, then the three refusals, then the census, then the projection.

## 2. Why the demotion happens to the file

§2.1's trap for the coder is that **demotion is per-file, not per-number** — there is no state in
which some of a file's numbers stay fitted. The shortest way to *be* that is to demote the file,
whose warrants are five public fields, and project it once. Demoting afterwards would mean reaching
into `RunParameters`, `StratumFits` and `Estimate` separately, which is five chances to leave one
behind.

**And the walk destructures without `..` at every level**, the idiom `CensusIdentity::of` already
uses over `RecordingTerms`: a warranted number added anywhere in the file — a new section, a new key
in a section, or a second warranted key on an existing row — **stops the demotion compiling**.
Measured both ways: a sixth `WarrantedValue` on `StatedConstants` and a second one on
`SubstitutionRateRow` each give `error[E0027]: pattern does not mention field`.

## 3. `weaker_of`, and why §13's fifth test cannot be met as written

**Owner's ruling of 2026-08-30**: the demotion is `Provenance::weaker_of(file's warrant, Supplied)`
and never an assignment. `Provenance` ranks `Supplied` **above** `Defaulted` — a number the run was
handed says nothing about this data, and a stated constant says less than nothing — so assigning
would *promote* every defaulted number into a claim that somebody chose a value nobody chose.

**⚑ So *every warrant `Supplied`* is not true of a demoted file and cannot be.** Spec §13's fifth
test says it, spec §2's own table contradicts it by ranking the four weakest-last with `Defaulted`
last, and the code follows §2. On the module's own fixture two numbers stay `Defaulted` after a
demotion and both must: read group 1's calibration multiplier of exactly 1.0, and the repeat-tract
outlier weight. **What is true is that no warrant survives stronger than `Supplied` and none was
promoted**, and that is what the test asserts. Recorded in `PROJECT_STATUS.md`; the spec is
untouched.

## 4. What the demotion cannot reach, and why that is a defect and not a nicety

Five numbers carry a `Warrant` and all five are demoted. The rest carry other vocabularies with no
*handed over* state: a slippage number says it came off the stratum's own fit, its period's curve,
or a blend; the prior seed says which rung it came from; a contamination fraction says which reads
it was fitted from.

**So a demoted file still says *this run's own* about numbers that are not this run's.**
`SeedRung::FittedCurve` reads "both moments came off **the run's own** fitted population curve", and
after a demotion the run that fitted it is a different run.

**This is open and it is the owner's.** `PROJECT_STATUS.md` records it, offers three ways out, and
recommends the one D3 did *not* take — refusing such a file like the other three bindings. D3 builds
what the plan and §2.1 describe; if the owner takes that recommendation, this method and the door
above it go. **The doc comment says so** rather than reading as settled, which is what an earlier
draft of it did.

## 5. Tests

Ten added. The ones that carry the step:

| test | what it holds |
|---|---|
| `the_door_demotes_and_not_only_the_method` | **the whole of what D3 composes** — a warrant on the far side of the door, including the bottom rung the base commit made carriable and the outlier weight that must not move |
| `the_demotion_changes_no_number_a_locus_is_scored_against` | ploidy, read groups, prior seed, contamination, batching, inbreeding, every calibration scale, every substitution rate, and the slippage fit and length spectrum at every stratum the file names |
| `no_warrant_survives_the_demotion_stronger_than_supplied` | every warranted number, against `weaker_of` itself, with the fixture asserted to hold something to demote **and** something that must not be promoted |
| `the_outlier_weight_the_project_guessed_is_not_promoted_to_one_somebody_chose` | the ruling's own case, both directions |
| `a_term_renamed_with_its_digest_unmoved_disagrees` | the half of the census comparison a moved digest cannot exercise |
| `a_census_naming_a_different_number_of_terms_disagrees` | both directions of a length mismatch |
| `the_door_runs_validate_first_and_then_the_three_refusals` | a file at odds with itself is named by `validate`, not blamed on the run |
| `a_demoted_file_still_validates_and_still_projects` | the regression on the base commit's `validate` change |

## 6. Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib ng::calling::parameters_file`: **170 passed, 0 failed** (160 before).
- `cargo test --lib`: **5,541 passed, 0 failed** (5,531 before).
- `cargo doc --no-deps`: **25 unresolved-link errors**, the baseline. *(The first draft of the
  rewritten module header added a twenty-sixth, `[`Self::from_toml`]` in a module-level doc where
  there is no self.)*

## 7. Mutation testing

Between them the two reviewers ran fourteen mutants and **four survived**. Three were real holes and
are fixed; all three are now killed, each re-run on the final tree:

| mutant | before | now |
|---|---|---|
| **the door notices the disagreement and does not demote** | survived | fails `the_door_demotes_and_not_only_the_method` |
| **a sixth warranted number added upstream, `fitted_here`** | survived | **does not compile**, and so does a second warranted key on a row |
| **the census comparison ignores a term's name** | survived | fails `a_term_renamed_with_its_digest_unmoved_disagrees` |
| the census comparison names the run's term, not the file's | survived | fails the same test |
| the demotion assigns `Supplied` instead of taking the weaker | 2 | 3 |
| the demotion also rewrites the sequencing batching | survived (the "same numbers" test did not compare it) | 4 |
| `census_disagreement` always `None` | 2 | 2 |
| each of the four other demotion sites, dropped | 1 each | 1 each |
| the door skips `validate`; the door skips the refusals | 1 each | 1 each |
| the different-lengths tail, dropped | 1 | 1 |
| **the outlier weight's demotion, dropped** | survived | **survives, and is equivalent** |

**The last one cannot be killed and should not be chased.** `validate` holds that key to `supplied`
or `defaulted`, and `weaker_of(·, Supplied)` is the identity on both — so the line changes nothing
observable. It stays because §2.1 has no per-number exemption and the destructured walk has to be
total.

## 8. What the review found

Two agents in isolated worktrees, correctness and design fidelity. They agreed on all three Majors.

- **The demotion was untested where it matters.** Every warrant assertion was made against the
  demotion *method* called directly, and the one test that went through the door looked at the
  agreeing case. So deleting the demotion from the door — the whole of what D3 composes — left all
  168 tests green, with the symptom being identical genotypes and a run that overstates every
  warrant it prints. That is precisely this step's own silent failure.
- **The guard against a sixth warranted number did not exist, and a comment said it did.** The test
  helper asserted a count of 8, which is the number of *rows the fixture happens to hold*; a
  reviewer added a sixth `WarrantedValue` to the shape and the suite stayed green. The walk
  destructures now.
- **Half of the census comparison was dead to the suite.** No test moved a term's *name*, so a
  build that renamed a term while its digest stood still would read as agreeing — and every
  `fitted_here` in a file from that build would survive into a run it was not fitted for.

Four comments said something the code does not do, and each is fixed: the count assertion above; a
citation calling an **open** `PROJECT_STATUS.md` question "the owner's" when its own recommendation
is the opposite of what D3 built; "both of those are warrants the demotion moves", where `defaulted`
is exactly the one it leaves alone; and "the one where that is visible" for a promotion that is
visible on at least two numbers.

Also taken: the "same numbers" test did not compare the sequencing batching, which is the population
a contaminating read is drawn from and so a number a locus is scored against; `ParametersForThisRun`
earns its place (folding the verdict into `RunParametersFromFile` would make `None` mean *no census
was compared* and *the census agreed* at once — the sentinel §5 exists to prevent) but its field is
`from_file` now, so a caller writes `from_file.parameters` rather than one word twice; and the
demotion is named `demoted_to_no_better_than_supplied`, because the shorter name asserts at every
call site the sentence §3 above exists to correct.

**And one hazard only the correctness pass found:** the demotion is public, and called *before*
`validate` it launders an illegal `fitted_here` outlier weight into a legal `supplied` one — so
`file.demoted_to_no_better_than_supplied().to_run_parameters()` accepts what
`file.to_run_parameters()` refuses. The door is safe because nothing is demoted until `validate` has
passed, which makes that ordering load-bearing for more than the message quality its comment gave as
the reason. Said there now.

## 9. What Milestones E and F are owed

Recorded in `PROJECT_STATUS.md`, none of it D3's to settle:

- **The door cannot be called by direct mode**, the mode this file is on the critical path of. It
  takes a `&CensusIdentity` and a direct-mode run has no census — no pre-pass, no psp. §6's fourth
  binding has no wording for a run with no census either. **This one needs the owner.**
- **§13 test 5's other half is not delivered**: the demotion is to be "visible in the run's report",
  and `RunParameterReport` holds three fields, none of them this. F2's.
- **What `[fitted_from].census` a demoted run writes back** (§7 makes writing unconditional). Write
  the file's and every re-run demotes again; write the run's own and the demotion is invisible for
  ever. F1's.
- **The verdict names a term and carries no values**, where D2's refusals were ruled to carry both.
  Whatever prints it will want them.
- **`validate` has no rule about the census block at all**, so a hand-reordered term list is
  silently demoted naming a term nothing differs on. A later `validate` rung.
