# Code review — ng calling loop E2: the run's frozen parameters

**Scope:** the working-tree diff of step E2 of
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md), on top of `a2db2ae1` — a new
`src/ng/calling/run_parameters.rs`, one field and one lookup on `FrozenParameters`, and the call
churn that field caused.
**Date:** 2026-08-26. **Verdict: request changes** — **1 Blocker, 5 Majors, 5 Minors, 5 Nits**, and
**6 of 24 claims wrong**. All applied; see [the fix report](fixes_applied_2026-08-26_e2.md).

**One agent, in its own worktree**, carrying reliability, step 8a's re-derivation and a lighter
naming/errors/smells pass — proportionate to a diff that is one type, one constructor and one
lookup. **19 mutations run, 4 survived, 3 changed no behaviour.** Its findings could not be written
to the shared scratch directory from inside the worktree; they were returned in full and are
synthesised here.

---

## Blocker — a mechanism stated in six places and true in none

The module said a gap in the read-group ids "slides every later group's calibration onto its
neighbour", with the symptom "a wrong genotype rather than a crash" — in the module doc, the
constructor's `# Panics`, the helper's doc, the assert message, a test's doc, and the
implementation report.

**Nothing slides.** Both dense vectors are built over `0..count` by *keyed lookup*, so slot *i*
holds read group *i*'s own value or `defaulted()`. Measured, with ids `{0, 2}`: read group 0 keeps
its own calibration, slot 1 is `Defaulted`, and **read group 2 is dropped entirely** — the vector
is two long, so `calibration_of` panics the first time a read from that library is scored.

**Why it matters more than a wording fix:** a reader of any of the six sites concludes some other
path must be re-checked for slid indices, and does not learn the real cost — a dropped library
whose failure arrives at a locus, long after the pre-pass finished, naming the read group and not
the run. **The check stays**, for that second reason. All six sites are corrected.

## Majors

**M1 — the contamination map's keys were never checked, and the test that read as though it
checked them could not.** An estimate for a read group past the axis was dropped in silence, which
leaves a contaminated library uncorrected. And the walk that builds the views was pinned by a
fixture giving *every* read group an estimate, so "one view per read group" and "one view per
estimate" were the same number: measured, a per-estimate walk passed the whole suite while charging
read group 1's 3% to read group 0 and inventing a fraction for a read group the run does not have.
**Fixed**: the axis is checked, the fixture has a fourth read group with no entry, and the off-axis
case has its own test.

**M2 — `UNMEASURED_READ_GROUP` fabricates a provenance.** Its `source` says
`TheWholeSamplesReads`, which that type documents as *"fitted from every read of the sample and
copied onto this read group"* — a positive claim about a number never fitted, on the one field a
consumer reads to answer whether two libraries of one plant may be said to differ. The fraction and
the counts honour *absent is not a fitted zero*; `source` cannot, because the type has no
*not measured* variant. **Documented at the constant and pinned by a test**; making it
unrepresentable belongs to `ContaminationView`'s owner.

**M3 — only one direction of the rate/total pairing was tested.** Dropping the minted map from the
key union left the suite green and the read group silently absent from the run. The reverse is the
likelier direction — the accumulator runs over every read of every library, and a fit can decline
to model one. **Fixed** by a test.

**M4 — the diff added a broken intra-doc link**, which the crate denies. `cargo doc --no-deps` is
not in the branch's gate, which is why nobody saw it. **Fixed**, and the module now contributes
none of the build's remaining 28 pre-existing ones.

**M5 — the ploidy in the substitution lookup's key was unpinned.** Every fixture was diploid, so
hard-coding two survived the suite; on a haploid run whose map is keyed at ploidy 1 the difference
is `Some(0.001)` against `None` — and `None` is what the lookup's own doc calls an ordinary
absence, so a whole haploid contig would quietly lose its substitution rates. **Fixed** by a
haploid test, and by one that pins the read group in the key too.

## Minors, applied

The `.zip()` pairing a rate with its total has an **unreachable** empty arm — the check upstream
guarantees both — and its comment claimed that arm was the `Defaulted` route; an `expect` planted
there never fired across the suite. The inbreeding fixture was `[0.0, 0.9, 0.0]`, a **palindrome**,
so reversing the stored order left a bit-identical vector and the suite green; it is now
`[0.1, 0.9, 0.5]` with every slot asserted. `read_group_count_of` → `checked_read_group_count_of`,
because it is where two of the four rules are enforced and the old name said only that it counted.
Two stale counts in the docs ("seven things" for eight parameters; "three maps" for four). A stray
mid-file import.

## Deferred, with a recommendation

**The substitution lookup returns a bare `Option` where its sibling returns a typed error.**
`StratumFits::at` separates *unknown read group*, *no such stratum* and *group not in the fit* —
and that last, its comment says, "is not a quiet library". The two lookups fill the same scoring
context, so the assembled context will carry the weaker answer. **Mirror the error type when the
tract's context assembly is built**, which is the step that first holds both.

## Six wrong claims of twenty-four

The Blocker's mechanism (three sites in the code, one in the report); "which is a wrong genotype
rather than a crash"; "either half missing is the honest `Defaulted`" (unreachable); "seven things"
where the constructor takes eight; "the three maps … keyed by `ReadGroupId`" where there are four
and one is keyed otherwise; and a test assertion whose message claimed evidence its fixture could
not give.

**Everything numeric checked out**, including the calibration scale of `0.500000000282` and the
report's attribution of that miss to the accumulator's fixed-point sum — the reviewer re-derived it
from the quantisation bound and got the same 5.64 × 10⁻¹⁰ relative shift.

## Verification

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0; `cargo doc --no-deps --lib` — no error in this module.
- `cargo test --lib` — `4733 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `687 passed; 0 failed; 3 ignored`.
- The five release-held checks downgraded together: **6 failures**, every check reached.
