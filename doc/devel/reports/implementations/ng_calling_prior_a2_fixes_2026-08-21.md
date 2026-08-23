# ng genotype prior — A2 review fixes

*Fix-application report, 2026-08-21. Branch `ng-calling-prior`. Applies
[the A2 review](../reviews/ng_calling_prior_a2_2026-08-21.md) to step A2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md).*

## Summary

**All eight Majors and all thirteen Minors applied. Nothing disputed, nothing deferred.** The
file went from 525 lines to 925 (`git diff --numstat`: +863 / −2), of which the test module is
about half.

**The headline change: the seam's eight parameters became one checked bundle.** Three of the five
review agents identified the same defect independently and each compiled the same class of fix —
that the eight-argument method and the shape checks no implementation could be forced to call were
one fact seen from two sides. The trait is now
`fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF)`, and
`PriorRow::new` is the only way to build its argument.

**And a defect none of the three prototypes caught, found by re-running their own claim.** The
review's proposal was that a private-fielded bundle makes the checks unskippable. **It does not, in
the arrangement they compiled** — see §2. That is the one place this fix pass goes beyond what the
review asked for, and it is the part that makes the rest true.

## 1. The bundle, and why the reviewers' version did not work

The proposal was `PriorRow<'a>` with private fields and a checking `new`, declared in
`genotype_prior/mod.rs`. The claim was that no implementation could then reach the fields.

**A private field is visible to a module's descendants**, and
`genotype_prior::dirichlet_multinomial`, `::plug_in`, `::seed_spectrum` and `::seed_ssr` are all
descendants of `genotype_prior` — not siblings. Those four files are where every implementation of
`GenotypePriorModel` and every seed builder will live. Measured, before the fix: a probe added to
`dirichlet_multinomial.rs` built a `PriorRow` field by field, **compiled, and ran green**. Each
reviewer's prototype was exercised from `mod tests`, which is a descendant too, so none of them met
this.

The repair is one level of nesting. `Concentration`, `PriorRow` and `SpectrumSeed` now live in
`mod checked` inside `mod.rs` and are re-exported; the four sub-modules are siblings of `checked`
rather than descendants, so the fields are out of reach. Measured, after: the same probe now fails
to compile —

```
error[E0451]: fields `concentration`, `genotype_allele_counts`, `log_multinomial_coeffs`,
`homozygous_allele_for`, `per_allele_scratch` and `out` of struct
`calling::genotype_prior::checked::PriorRow` are private
```

The module's doc comment records the measurement, so the nesting cannot be flattened back as
tidying.

## 2. Major findings

| finding | what changed |
|---|---|
| **M1** — five clippy errors invisible to `cargo clippy --lib` | The `catch_unwind` table that carried four boxed closures is gone, replaced by six `#[should_panic]` tests — which removes the `type_complexity` and all four `redundant_locals` lints at source rather than silencing them. **`--lib --tests` is now one of this step's gates**, and `cargo clippy --lib --tests --all-features -- -D warnings` finishes clean. |
| **M2, M3, M4** — the eight-argument seam and the unenforceable checks | `PriorRow` bundles all six buffers; the checks are its constructor; the trait method takes three arguments; `assert_row_shapes` and `#[allow(clippy::too_many_arguments)]` are both deleted. Dyn-compatibility is kept and pinned by a test. |
| **M5** — "held in release" exercised by nothing | Confirmed and closed. Downgrading every `assert_eq!` in `PriorRow::new` to `debug_assert_eq!` leaves the debug run at **20 passed** and the release run at **FAILED. 11 passed; 6 failed**. `cargo test --release --lib ng::calling::genotype_prior` is now one of this step's gates, and the type's doc says which command pins which half. |
| **M6** — neither seam test could fail on a wrong genotype order | The seam test moved from biallelic to **triallelic**, whose expected row `[0, ln2, 0, ln2, ln2, 0]` is not a palindrome, so a reversed or permuted walk cannot pass. It also pins the two views the row's values never touch — the copy-count table and the homozygous lookup — against the layout the bundle documents. |
| **M7** — an empty `out` passed every check | `PriorRow::new` refuses it in release, with a test. Every locus has at least the all-reference genotype at any ploidy, so a zero-length row is a wiring bug for the same reason an empty concentration is. |
| **M8** — a mis-sized `out` made every message blame a correct array | Every length message now names **both** buffers it compared and where each number came from, and a test asserts that a six-genotype row over a three-genotype table produces a message naming `` `out` is sized for 6 genotypes `` and `` `log_multinomial_coeffs` holds 3 ``. |

## 3. Minor findings

`homozygous_allele_for` restored as the parameter name, matching the architecture, the plan and
the table's own field (Mi1). `Concentration::new`'s doc no longer claims to match production's
check — it states that it is **tighter** than production's `α > 0`, so `1e-13` passes there and
panics here, and why (Mi2). The `[CandidateAlleles]` link is qualified and `cargo doc` no longer
errors on it (Mi3). `assert_row_shapes` is gone, so Mi4, Mi5 and Mi14 dissolve with it — the
method is now the verb phrase `fill_genotype_log_priors`. `SeedRegime`'s two neutral variants each
state the `(1, θ)` shape they share and what separates them (Mi6), and `data_dominated` became
`census_sites_outweigh_regularizer` (Mi7). The module doc defines **seed** once, says what refines
it (nothing), and names the STR cohort EM's different meaning of the word so the collision is
visible (Mi8). The 2:1 ratio now carries its formula and two worked values — 1.998:1 at a human θ,
1.98:1 at ten times that (Mi9). `SpectrumSeed` has private fields, a checking `new` and three
by-value accessors, with the alternative total allowed to be exactly zero because a fully invariant
cohort is a real answer and the flooring belongs to the per-locus expansion (Mi10). The panic
payload is read as `String` **or** `&'static str` and says so when it is neither, instead of
defaulting to an empty message and blaming the assertion's wording (Mi11). Three value checks
gained tests — `+∞`, an over-long scratch, and an out-of-range homozygous id, the last as a
`debug_assert` since it is a value precondition (Mi12). A shape-coverage test runs ploidy 1, 2, 3, 4
and 8 across 1 to 6 alleles (Mi13).

Nits: the four dead `let view = view;` rebinds are gone with the closures; `Concentration` derives
`PartialEq` and takes `self` by value in both accessors; the floor renders as `1e-12` rather than
`0.000000000001`; the genotype-count multiply is `checked_mul`; "the GIAB trio" reads "GIAB's three
samples, each called on its own"; and every test-local binding is named for what it holds.

## 4. Tests

**Twenty, against nine before.** Six `#[should_panic]` cases replace the closure table, one per
mis-shape, so each failure names itself in `cargo test` output rather than through a loop variable.
New coverage: an empty row; an over-long scratch; a non-finite concentration entry; ploidy 1, 3, 4
and 8 and a single-allele locus; the full triallelic row and both table views; the length message
naming its yardstick; `SeedRegime::NeutralShape` actually constructed; and `SpectrumSeed`'s two
refusals plus its zero-alternative case.

## 5. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | `Finished dev profile … in 10.32s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `test result: ok. 20 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `running 17 tests … ok. 17 passed; 0 failed` |
| `cargo test --lib` | 0 | `test result: ok. 4027 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 647.25s` |

Seventeen and not twenty under `--release`: the three `#[cfg(debug_assertions)]` value-check tests
are compiled out there, which is the point — they pin the half of the split that is debug-only.

## 6. Deviations this leaves against the design documents

Four, all in `arch/calling_priors.md`, all recorded rather than absorbed silently:

1. **§2.1** spells the species-range fallback as a module-level `pub const … : f64`; it is
   `ExpectedHeterozygosity::SPECIES_FALLBACK` (from A1).
2. **§2.2** shows `Concentration` as a free type in the module; it is inside `mod checked`, for the
   measured reason in §1 above.
3. **§2.3** shows `SpectrumSeed` with public fields; they are private behind a checking constructor.
4. **§3.2** shows a six-parameter row function taking flat slices; it takes a `PriorRow` and an
   `InbreedingF`. **The contract is unchanged** — the same six buffers, all caller-owned, nothing
   allocated, the same flat views, no back-reference into the loop — which is what §7's decision
   is about and what the arch itself says is the deliverable, signatures being illustrative.

**Still open for the owner** (unchanged from Checkpoint A): whether `seed_spectrum.rs` and
`plug_in.rs` should be renamed, which touches the arch's file tree and the plan's scope section.
