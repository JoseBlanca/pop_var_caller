# ng step 4, the STR path — C2: one read group's reads, counted into the offset buckets

*Implementation report, 2026-08-11. Step C2 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md) — one of the plan's
**own commit, do not bundle** steps, because its failure is silent. With the review that followed
and the fixes applied: two agents, 27 mutations, 4 behaviour-changing survivors. Design authority:
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.3 and
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §3, §4.1.*

## What the step is

`tally_of`: how far each of one read group's reads sits from the **reference** tract's length, in
whole motif copies, counted into nine offset buckets plus the guard.

**Complete witnesses only, and that is why the plan gives this step its own commit.** A read that
reached both borders of the tract measured its length; a read that reached only one saw *part* of
it, so what it shows is a lower bound. Scoring a lower bound as a length reads as a read that lost
repeats — a direct bias in the direction split, the parameter this path exists to protect and the
one that inverts on real data when an estimator gets it wrong. The partial witnesses are counted
rather than dropped in silence, because leaving reads out is a size a run has to be able to see.

## Recorded deviations from the architecture

1. **A tally, not a `LocusShape`.** Arch §2.3 sketches `shape_of(locus, read_group, cap) ->
   LocusShape`. The read cap is a *random subsample* (C3), and a draw needs the whole tally to draw
   from — while a `LocusShape` refuses to hold more than `MAX_LOCUS_READS`. So the architecture's
   `shape_of` is this function followed by C3's draw. The module doc says so.
2. **The period is read from the locus's own motif rather than passed in** — see below; this began
   as a parameter and the review took it out.

## What the review changed

**Blocker — nothing pinned that a bucket *accumulates*.** Every fixture put at most one observation
in any bucket, so replacing `+=` with `=` left all 87 tests green. It is not a hypothetical shape:
observations are keyed on `(bases, witness, read group)`, so an interruption or a substitution
inside the tract gives two entries at the *same* length and therefore the same bucket, and the two
saturating ends pool many offsets by design. The mutant loses reads unevenly across buckets — which
is the direction split again. Reproduced here: it now fails two tests.

**Major — the period was a parameter, so a caller could measure reads in one unit and file them
under a stratum named in another.** The doc argued the caller had already asked `stratum_of` for
it and the two must agree — but passing it is exactly what lets them disagree, and the result would
be a full, plausible tally measured against the wrong origin. `tally_of` now takes the period from
the locus's own motif, as `stratum_of` does, so they cannot differ; a locus that is not one repeat
tract has no motif and tallies as empty.

**Major — the plan's `reads_discarded_by_cap` rule was neither stated in the code nor tested.**
The plan says a locus whose depth cap fired is entered at the depth observed, because the
generator's reservoir is a uniform subsample; skipping such loci would be depth-dependent
selection, the bias this step exists to remove. Every fixture set the field to zero, so an
early-return mutant survived and emptied the tally. Stated in the doc and pinned by a test, which
I reproduced failing on the mutant.

**Major — no property test, where the sibling function in the same file has one.** A conservation
law over an arbitrary list of observations now holds: whatever the lengths, the period and the mix
of libraries, the depth is the sum over that group's complete observations and the partial count is
the sum over its partial ones. It kills the Blocker's mutant without anyone having to think of the
case.

**Minor — two boundary lengths.** A complete witness showing *no* tract bases (the tract fully
deleted, which the generator really can emit) was droppable with no test failing, and it is the
furthest-short allele there is — exactly what the low end bucket exists for. And an empty reference
tract, which `stratum_of` deliberately accepts and pins, had no matching pin here. Both covered.

**One wrong number of mine.** The draft commit message said eight tests where the step added seven.
Everything else held: the fixtures' byte lengths and the offsets they imply, `DEFAULT_SSR_MAX_READS_PER_LOCUS`
at a thousand, and the spec's own record of the direction split inverting (0.9× over all loci
against 3.4× over known-homozygous ones).

**Also applied:** `reads()` renamed `depth()`, which is what `LocusShape` calls the same quantity;
the `i64` conversions now `expect` rather than saturate, since a saturated length would make the
subtraction overflow silently; and the module doc no longer calls this file's first reduction its
only one.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation::ssr` | **93 passed** (80 before this step) |
| `cargo test --lib --bins --tests --all-features` | **3,470 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` on `locus_offsets.rs` gives **23**, of which 10
are C1's, so this step adds **13**; the suite moved 3,457 → 3,470.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets` fails
in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-c2/` — two per-category files (reliability; errors, naming
and the numbers check) and the reviewed patch.
