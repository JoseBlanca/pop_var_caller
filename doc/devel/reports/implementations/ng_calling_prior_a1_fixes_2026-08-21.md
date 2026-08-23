# ng genotype prior — A1 review fixes

*Fix-application report, 2026-08-21. Branch `ng-calling-prior`. Applies
[the A1 review](../reviews/ng_calling_prior_a1_2026-08-21.md) to step A1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md).*

## Summary

**All four Major findings and eleven of the fourteen Minors applied; three deferred to the
owner at Checkpoint A because they rename items the architecture and the plan write down.**
Nothing was disputed.

The six mutations the review recorded as surviving the submitted suite were **re-run against the
fixed code, and all six now fail it** (§4). The change that closes M2 was checked the same way: the
cross-type confusion the reviewer compiled is now a compile error.

## 1. What changed, finding by finding

### Major

**M1 — the boundary sweep now covers the new scalar.** One line added to the existing proptest arm
in `the_constrained_rates_accept_exactly_the_probabilities_and_round_trip`, and its doc says four
rates rather than three. This is the test that asserts acceptance holds *exactly when* a value is a
finite probability and that an accepted value comes back bit for bit.

**M2 — the fallback is now a value of the type it defaults.** `DEFAULT_SPECIES_DIVERSITY_FALLBACK:
f64` is gone; in its place `ExpectedHeterozygosity::SPECIES_FALLBACK: Self = Self(1e-3)`, following
`AlleleId::REFERENCE`'s precedent in the same file. The reviewer's demonstration — that the bare
constant seeds `InbreedingF`, `ErrorRate` and `GenotypeFrequency` alike — no longer compiles:

```
error[E0308]: mismatched types
   |  InbreedingF::try_new(ExpectedHeterozygosity::SPECIES_FALLBACK)
   |                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ expected `f64`, found `ExpectedHeterozygosity`
```

**This also resolves Mi1**, the naming agent's separate complaint that `DEFAULT_…_FALLBACK` says
one idea twice and names the value `DIVERSITY` where the type it seeds is a heterozygosity: the
type now carries the quantity and the constant carries only what makes it special. **It is a
deviation from `arch/calling_priors.md` §2.1**, which spells the item as a free `pub const … : f64`
at module level — raised at Checkpoint A so the architecture can be brought into line.

One consequence recorded in the code: `Self(1e-3)` does not pass through `try_new`, so the type
would have one value its own predicate never saw. The replacement test asserts both halves.

**M3 — the fallback's value is pinned.** `the_species_diversity_fallback_is_a_constructible_heterozygosity`
compared the constant with itself and so passed for every value in `[0, 1]`. It is replaced by
`the_species_fallback_is_one_difference_per_thousand_bases`, which asserts the value is `1e-3` and,
second, that `try_new(1e-3)` returns that same constant. Its doc now says what the old one only
claimed — naming the two slips that survived (`0.1`, the percentage reading; `1.0`, the
per-kilobase reading) and that both were run as mutations.

**M4 — the doc obligation applied; the structural half deferred.** `SPECIES_FALLBACK`'s doc now
states three things it did not: that the thing which will report a run landing on it is
`SeedRegime::FallbackDiversity`, and that **any code path reading the constant owes that report**;
that it must be overridable and **no door exists yet**, naming the calling-loop plan as where one
lands; and that nothing reads it today. The structural options the `defaults` agent proposed —
moving the fallback down beside `SeedRegime`, or making the only public door a constructor
returning `(value, regime)` — are **not applied**: both change where A2 puts things, and A2 is the
next step. Carried into A2's scope, and open question 1 of the review is the one that decides it.

### Minor — applied

| finding | what changed |
|---|---|
| Mi2 | The STR sentence no longer claims a substitution production never made. It now says production's STR path never measured repeat diversity at all: it hardcodes freebayes' default `SFS_THETA = 0.01`, and that number is a population-scaled mutation rate — a different quantity in different units from a heterozygosity in `[0, 1]`. |
| Mi3 | The tomato claim now carries its size and its direction: this project's own tomato panel fits **below** the fallback — 6 differences per 10,000 bases against the 10 per 10,000 the constant holds — while a diverse outcrosser would sit above. |
| Mi4 | The `NaN`/infinity loop asserts the **variant**, not merely that some error came back — `matches!(…, Err(DomainError::ExpectedHeterozygosity(_)))`. Done for all four rates in the loop, not only the new one: leaving the new type as the only one pinned would have been the odder half. |
| Mi5 | The `DomainError` variant's justification is rewritten to the distinction that actually holds — this type draws its two chromosomes from the **cohort**, the heterozygosity `GenotypeFrequency` carries draws them from **one individual**, and the two differ by a factor of `(1 − F)`. |
| Mi6 | "Source:" now names both pre-pass routes, and says what is true of each: the joint fit supplies the number; the histogram route supplies the ingredient — each sample's observed heterozygosity, of which `θ` is the mean of `Hobs / (1 − F)` — and **nothing computes that mean yet**. Checked before writing: `grep heterozygosity src/ng/parameter_estimation/generic/*.rs` finds `observed_heterozygosity() -> Option<GenotypeFrequency>` and no expected one. |
| Mi9 | The endpoints test's doc no longer argues that `InbreedingF::try_new(1.0)` succeeding is correct. It marks the assertion as *today's* behaviour, says why the prior needs `[0, 1)`, and names both the plan step that tightens it and the fitted producer that must be clamped first. |
| Mi10 | The enum doc no longer forward-references `Theta` as a future arrival; "three rate constructors" reads four; the rejection test's doc names four and lists `ExpectedHeterozygosity` among the types with their own variant. |
| Mi11 | "which is what a softmax consumes" is gone: "Only the differences between them matter: the loop adds what the reads say and rescales the row to sum to one, so a constant shared by every entry cancels." |
| Mi12 | The Dirichlet-multinomial is defined where it first does work, in both files that use the term, and "marginalizing" is defined as averaging over the unknown frequencies rather than fixing them at one estimate. |
| Mi13 | `plug_in.rs` now carries the measured size: on the GIAB trio, each sample called on its own at 5×, 83.6% → 94.6% genotype accuracy at true variants and 214 → 8 sites where a two-copy carrier was called heterozygous, with the emitted variant set byte-identical — and says plainly that this is one corner. |
| Mi14 | "Cheap arithmetic" and "the expensive call" are replaced by the costs the architecture states: one addition per allele against one `lgamma` per allele a genotype carries a copy of plus one `logsumexp` per homozygous genotype. |

### Nits — applied

The placement rationale in the section header is restored, adjusted to what A1 made true (step 7
has no module, step 8's holds no code). `src/ng/mod.rs`'s crate index records that `calling` now
holds the first of its four sub-modules. "panel" → "cohort" in `mod.rs` and `seed_ssr.rs`; "(plan
step B)" → "(plan steps B1–B2)"; `seed_spectrum.rs` defines "concentration" for a reader landing on
its own rustdoc page; `seed_ssr.rs` states the `ssr`/STR convention once; the shared predicate's
user count goes six → seven and the crate's probability count seven → eight. The `-0.0` question is
answered rather than left open: the type is a transparent wrapper, unlike `Phred`, and a new test
asserts it in bits.

### Deferred to Checkpoint A — three renames the design documents fix

Not applied, because applying them means editing `arch/calling_priors.md` and the plan's own file
list, which this loop does not do on its own authority:

- **Mi7** — `seed_spectrum.rs` → `seed_generic.rs`. The two filenames are not parallel: one names
  the locus class it serves, the other the input it reads.
- **Mi8** — `plug_in.rs` → `hardy_weinberg.rs`. In a Rust tree `plug_in` reads as an extension
  point first.
- The **M2/Mi1 rename is applied** rather than deferred, because it is inseparable from the type
  change M2 required — but it carries the same question, and it is on the Checkpoint A list.

### Not applied, with the reason

**A test pinning `SPECIES_FALLBACK` against production's `DEFAULT_DIVERSITY_PRIOR`** (proposed
under M3, and reachable — the constant is `pub` through `src/var_calling/mod.rs`). Two reasons.
`src/ng/calling/mod.rs` records that ng's single `use crate::var_calling::` sits in one greppable
place, and a second oracle site in `types.rs`'s tests would make that sentence false. And the
coupling the test would assert is not one the project wants: production is frozen, so the
assertion could only ever fail because *ng* moved its own constant deliberately. The doc claim was
weakened to match — it now says the value was taken from production with its reasoning, that the
two are not tied, and that ng may move this one.

## 2. Files touched

Measured with `git diff --numstat -- src/`, the index holding the reviewed state:

| file | + | − |
|---|---|---|
| `src/ng/types.rs` | 158 | 64 |
| `src/ng/calling/genotype_prior/mod.rs` | 18 | 9 |
| `src/ng/calling/genotype_prior/plug_in.rs` | 12 | 4 |
| `src/ng/calling/genotype_prior/dirichlet_multinomial.rs` | 3 | 2 |
| `src/ng/calling/genotype_prior/seed_spectrum.rs` | 3 | 2 |
| `src/ng/mod.rs` | 3 | 1 |
| `src/ng/calling/genotype_prior/seed_ssr.rs` | 1 | 1 |

## 3. Tests

Two added, one replaced, four extended.

| test | what it pins |
|---|---|
| `the_species_fallback_is_one_difference_per_thousand_bases` (replaces the constructibility test) | the value is `1e-3`, and the associated const agrees with the constructor it bypasses |
| `expected_heterozygosity_rejection_names_its_own_quantity` (new) | the rendered message names the diversity, not a neighbouring fit |
| `expected_heterozygosity_carries_negative_zero_verbatim` (new) | `-0.0` survives the round trip in bits — the opposite of `Phred`'s deliberate normalisation, recorded rather than left undecided |
| `the_constrained_rates_accept_exactly_the_probabilities_and_round_trip` (extended) | the exact accept/reject boundary and bit-for-bit round trip, for the fourth rate |
| `the_constrained_rates_reject_nan_and_the_infinities` (extended) | the **variant**, for all four rates, not merely that an error came back |
| `the_constrained_rates_accept_both_endpoints`, `each_constrained_rate_rejects_out_of_range_in_both_directions` | unchanged assertions, corrected rationale |

## 4. The mutations, re-run against the fixed code

Every survivor the review recorded, applied to the fixed tree one at a time and reverted, run with
`scripts/dev.sh cargo test --lib ng::types::tests` (baseline `50 passed; 0 failed`):

| mutation | before the fixes | after |
|---|---|---|
| `try_new` range widened to `(-0.25..=1.25)` | survived | **FAILED. 49 passed; 1 failed** |
| `try_new` quantises to six decimals | survived | **FAILED. 49 passed; 1 failed** |
| `SPECIES_FALLBACK` `1e-3` → `1e-1` | survived | **FAILED. 49 passed; 1 failed** |
| the `#[error]` message reworded to "genotype frequency {0} …" | survived | **FAILED. 49 passed; 1 failed** |
| non-finite inputs return `DomainError::GenotypeFrequency` | survived | **FAILED. 49 passed; 1 failed** |
| `get()` returns `self.0.abs()` | survived | **FAILED. 48 passed; 2 failed** |

**Six run, six killed, none surviving.** The tree was restored from a copy after each and the
baseline re-run green at the end.

Separately, the type-confusion probe the `defaults` agent compiled against the old code —
`InbreedingF::try_new(ExpectedHeterozygosity::SPECIES_FALLBACK)` — now fails to compile with
`error[E0308]: mismatched types … expected f64, found ExpectedHeterozygosity`.

## 5. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 3.00s` |
| `cargo test --lib ng::types::tests` | 0 | `test result: ok. 50 passed; 0 failed; 0 ignored; 0 measured; 3968 filtered out` |
| `cargo test --lib` | 0 | `test result: ok. 4007 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 658.71s` |

## 6. Follow-ups this leaves open

1. **The fallback's provenance has no carrier.** A2 mints `SeedRegime`; whether the fallback moves
   beside it, or a paired constructor becomes the only public door, is review open question 1.
2. **No override door.** Naming the calling-loop plan in the doc is not the same as having one.
3. **Two filename renames** await the owner (Mi7, Mi8), and the applied one (M2/Mi1) needs
   `arch/calling_priors.md` §2.1 updated to match.
