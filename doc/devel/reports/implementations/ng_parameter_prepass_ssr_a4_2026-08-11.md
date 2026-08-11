# ng step 4, the STR path — A4: the three slippage rates

*Implementation report, 2026-08-11. Step A4 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), run under the
plan-driven-implementation skill. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §1.1, §3, §4.5 and
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.4. Includes the review
that followed it and the fixes applied — four agents, 24 mutations, 4 survivors.*

## What the step is

Three constrained rates and the model that holds them. A read at a repeat tract can **slip**,
showing a whole motif copy more or fewer than the allele it was drawn from, and three numbers
describe that: how often a read slips at all, which way it slips when it does, and how far. The
fourth number this path fits — the per-base substitution rate — is not slippage and is a division
rather than a search, so it is not here.

**Three types and not one shared probability.** All three are fractions in `[0, 1]`, so one type
would let the gain share be handed to something expecting the level and compile. The review made
this argument sharper than the first draft did: the two are not reliably far apart, and at tomato
dinucleotides of 12–15 repeats the level reaches 0.150 against a gain share of 0.17 — within
1.1-fold, where a transposition would look like nothing at all.

## Changes made

**[`src/ng/parameter_estimation/ssr/slippage.rs`](../../../../src/ng/parameter_estimation/ssr/slippage.rs)**
— `SlipRate`, `SlipGainShare`, `SlipStepDecay`, each `MismatchFraction`'s shape (private field,
checked `try_new`, `get`); `SlippageModel` holding the three, with `new` over the checked types and
`try_new` over three plain fractions; a `Display` that names each number's own denominator.

**[`src/ng/types.rs`](../../../../src/ng/types.rs)** — three `DomainError` variants, and
`checked_probability` widened from private to `pub(crate)`.

**Recorded deviation: the visibility change.** The predicate's own doc comment argues against a
second spelling of the same range test ("… is how one of them ends up written `0.0..1.0` and
rejecting a genotype frequency of exactly one"), so the three new rates reuse it rather than adding
a fourth. A reviewer checked the widening leaks nothing: the signature already takes
`fn(f64) -> DomainError` and `DomainError` is `pub`. The same reviewer found the argument was
understated — `SiteNoise::try_new` in `generic/` already spells the test by hand, so the drift is
not hypothetical. The doc comment now says so.

## Tests

Eight in the file. Both endpoints accepted and round-tripped; **both bounds** rejected on all three
types, under each one's own error variant; `NaN` and both infinities refused; the model carrying its
three numbers in their own roles; **whichever** column is bad reported, including the first; the
**leftmost** bad column reported when more than one is; the rendering asserted character for
character; and a barely-slipping stratum proven to render differently from one that never slips.

## What the review changed

Four agents in isolated worktrees — reliability, errors + defaults, naming + idiomatic, and the
quantitative-claim check. **24 mutations, 4 survivors, 0 changed-no-behaviour.** One Blocker, one
Major, five Minors, and four wrong claims.

**Blocker — the first column was never tested.** Replacing `SlipRate::try_new(slip_rate)?` with an
unchecked construction left all six original tests green, and under it
`try_new(NaN, 0.17, 0.065)` returned `Ok` carrying a `NaN` level. The code was right and nothing
held it there. Fixed by testing all three columns.

**Major — the rejection tests crossed one bound per type.** `SlipRate` was never offered a value
above one and `SlipGainShare` never a negative; both widenings survived. `types.rs` states this
standard for its own three older rates in as many words. Fixed.

**Minor — the rendering could not tell "barely" from "never".** At `{:.4}` the bottom of the
measured range, 0.00091, prints as `0.0009`, and anything under 0.00005 prints as `0.0000` — the
same text as a genuine zero, which `SlipRate`'s doc calls a real answer. Nearly half the loci this
path sees sit in strata at that level. Now five decimals, with a test that separates the three
cases.

**Minor — "first failure wins" was documented and unpinned.** A fixture with one bad column gives
the same answer under every ordering of the three checks. Now pinned with two bad columns.

**Minor — the gain share did not say what it was a share of.** Its two readings differ by a factor
of the slip rate — about fiftyfold at a level of 0.02 — and neither the type name nor the message
disambiguated. The error message now reads "gain share of slipped reads", and the `Display` line
names all three denominators: `0.02010 of reads slipping, 0.170 of those gaining, 0.065 of those
taking a further step`.

**Four wrong claims, every one an addition of mine rather than a quotation** — the same split as the
previous step, and as every milestone of the sibling plan. 22 of 25 design figures checked out. The
one that mattered: I wrote that 5,072 reads sat "a repeat short **of the reference**", when that
measurement is against each unit's **modal observed length** — the origin the design rejects for the
accumulator, and the distinction §4.1 exists to establish. Also corrected: "0.17 of it gaining" in a
doc comment where the code emits `0.170`; "the two differ by two orders of magnitude", which holds
only against the bottom stratum; and "reproduces the same inversion", where the modal origin returns
0.48 — losses still ahead by 1.1-fold, so the same *size* of collapse, not an inversion.

**Declined, with reasons:** renaming `SlipGainShare` to `GainShareOfSlips` (diverges from the
architecture's name; the message and the rendering now carry the denominator, which is where a
reader meets it); replacing `SlippageModel::try_new`'s three positional `f64`s (it is the one door a
table of numbers must come through, and the review confirmed every other call site is protected —
transposing `new`'s arguments is `error[E0308]`); a macro over the three identical newtypes
(`types.rs` spells out four such types without one).

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --all-features` | **3,399 passed, 0 failed, 10 ignored** across 12 binaries |
| `cargo test --lib parameter_estimation::ssr` | 22 passed |

Counted from the tree rather than recalled: `grep -c '#\[test\]'` on `slippage.rs` gives **8**, and
the suite moved 3,391 → 3,399 with this step.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-a4/` — four per-category files and the patch each agent
reviewed.
