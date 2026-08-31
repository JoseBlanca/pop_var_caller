# ng — the run's parameters as a file: implementation plan

**Status:** draft, 2026-08-28. The build order for the **parameters file**: the TOML artefact that
carries every number calling runs on, each beside its warrant. Design is settled in
[`../spec/parameters_file.md`](../spec/parameters_file.md). This turns that design into build
order; it is **not** a place for new design — the three open items are §12 there, and one of them
(question 1) is measurement rather than a decision.

**Why this plan comes first.** Direct mode cannot run without it
([`../spec/run_streaming.md`](../spec/run_streaming.md) §2), and direct mode is the build order
that spec chose because it needs no file format. So this is the first thing on the critical path to
ng calling anything from alignment files.

---

## Scope

**In:** a new `src/ng/calling/parameters_file/` module — the file's Rust shape, its TOML writer,
its reader, the four bindings and their refusals, the compiled-in defaults, and the rule that every
run writes the file it used.

**Out (later plans, or already owned):**

- **How any number is estimated** — step 4's eight specs and the code under
  `src/ng/parameter_estimation/`, all built. This plan reads their results and writes them down.
- **The default slippage numbers themselves** — fitted from GIAB HG002, and the measurement does
  not exist ([`../spec/parameters_file.md`](../spec/parameters_file.md) §12 question 1). This plan
  builds the slot and the provenance that marks it defaulted; a research note fills it.
- **What the VCF header prints from the file** — the emission step's document, unwritten.
- **A command that writes the defaults without running a caller** — the command surface's, and
  cheap to add once the writer exists.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The round trip is the heart, so it is built before anything reads a real fit.** Milestones A–C
  prove that a `RunParameters` written and read back is the same object, on hand-built values.
  Everything after that is bindings, defaults and wiring, none of which can be trusted until the
  round trip is.
- **Reuse over rewrite.** `Provenance`, `Estimate<T>` and `RunParameters` exist and are not
  redesigned; the file is their serialised form. `toml` is already a dependency
  ([`src/psp/header.rs`](../../../../src/psp/header.rs)).
- **Verify against a real fit, not a fixture.** The north-star test is a `RunParameters` assembled
  from the joint fit on real tomato records, written and read back to an equal object — not a
  hand-built one, which cannot exercise the shapes the fit actually produces.
- **Isolate the silent failures.** Three steps here produce a quietly-wrong answer rather than a
  crash: float round-tripping, the absent/zero/defaulted distinctions, and the demotion rule. Each
  lands as its own commit with its oracle green before and after.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds.** All `cargo` through `./scripts/dev.sh` by absolute path (`CLAUDE.md`).

## Preconditions (already in place)

- **`RunParameters` and its assembly** — [`src/ng/calling/run_parameters.rs`](../../../../src/ng/calling/run_parameters.rs),
  with the four rules that turn the pre-pass's outputs into what calling reads.
- **`Provenance` and `Estimate<T>`** — [`src/ng/parameter_estimation/mod.rs`](../../../../src/ng/parameter_estimation/mod.rs).
  `Supplied` already exists and already sits below `Borrowed` in the warrant ladder; this plan adds
  no variant.
- **`RunParameterReport`** — [`src/ng/calling/run_report.rs`](../../../../src/ng/calling/run_report.rs),
  today reaching only contamination, batching and the inherited outlier weight, and called only by
  tests. Milestone F makes it a view over the file rather than a parallel structure.
- **The `toml` crate**, already a dependency for production's psp header.
- **A real fit to test against** — the joint fit runs on real tomato records
  (`examples/ng_joint_records_walk.rs` and the reports beside it), which is where Milestone C's
  oracle input comes from.

**There is no architecture document, and by the owner's decision (2026-08-28) one is not written
first.** The coder proposes the TOML tree and the key names in Milestone A, and they are **revised
later once there is a real file to look at** — reviewing names on paper is worse than reviewing them
against a file a fitted run actually produced.

**What that licence covers and what it does not.** The *spelling* is the coder's: key names, nesting,
where a table starts, whether a row is inline. The *content* is not — §3 lists what must be there and
§5 lists the states that must stay distinguishable, and dropping or merging one of those is a
stop-and-ask, not a naming choice. **Say what was chosen and why in A1's implementation report**, so
the later revision has something to argue with; a format nobody wrote down cannot be revised, only
re-derived.

**The revision has a trigger, not a date:** the first time a person reads a file this writer produced
and has to ask what a key means. Until then the provisional names stand.

---

## The steps

### Milestone A — the shape (types, no I/O)

✅ **A1. The file's Rust shape, and the provisional TOML tree.** One type per §3 section: identity
and binding, ploidy, the per-read-group calibration rows, the contamination table and its batching,
the per-sample inbreeding rows, the prior seed, the repeat-tract tables, and the compiled-in
constants. Serde derives, no reading or writing yet. **The key names and the shape of the tree are
the coder's proposal** (see the preconditions): choose them, and write down in the step's report what
was chosen and what was rejected, so the later revision has a starting point.
*Depends:* —. *Source:* [`parameters_file.md`](../spec/parameters_file.md) §3.

✅ **A2. A value and its warrant, as one serialised shape.** How `Estimate<T>` appears in the file:
the value, the warrant, and the observation count behind it. One shape reused by every numeric row,
so a reader cannot meet two spellings of the same idea.
*Depends:* A1. *Source:* §2.

✅ **A3. Absence is a missing key, never a sentinel.** The five states of §5 expressed in the types:
an absent contamination table, a measured-and-clean row, a defaulted calibration scale, a stratum
with no length spectrum, a missing `(stratum, slippage group)` row. `Option<T>` throughout; no
in-band zeros.
*Depends:* A1, A2. *Source:* §5.

> **Checkpoint A:** the types cover §3 and the five states of §5 are expressible and distinct.
> Pause for review.

### Milestone B — writing

✅ **B1. `RunParameters` → the file shape.** The projection, with no TOML yet: every field of
`RunParameters` mapped to its row, the read-group and sample axes carrying their names as well as
their indices.
*Depends:* A3. *Source:* §3.

✅ **B2. The file shape → TOML text.** Serialisation, with the layout choices §4 names: per-sample
rows as one inline table a line, the repeat-tract tables as arrays of arrays rather than arrays of
tables.
*Depends:* B1. *Source:* §4.

✅ **B3. The provenance comments.** Each defaulted value gets a comment beside it saying where the
default came from — for the repeat-tract slippage numbers, which alignments and at what depth
(§8). **This is why the format is TOML**, so it is built here rather than deferred.
*Depends:* B2. *Source:* §4, §8.

> **Checkpoint B:** a fitted `RunParameters` writes a file a person can read. Pause for review.

### Milestone C — reading, and the round trip that is the whole point

✅ **C1. TOML text → the file shape.** Parsing, with a malformed file failing at a line number.
*Depends:* B2. *Source:* §4, §9.

✅ **C2. The file shape → `RunParameters`, and the reader's `validate`.** The reverse of B1,
including the dense read-group axis over `0..n` that `RunParameters` requires.

**Landed as two commits, 2026-08-30**: the `validate` half (`6d67dc43`) and the projection. The
split was because the projection needed constructors that did not exist on types outside this
module — **three of them, not the two predicted**: `RunParameters`, whose only constructor takes
the *fit's raw inputs* rather than assembled values; `StratumFits`, whose `over` likewise takes the
fit's own outcome types; and `SequencingBatches`, whose `declared` takes sets of read groups plus
the run's `ReadGroups` table and derives the sample column from it, which a file has already been
through once. **`validate` had no caller in a run until the projection became one.**

**The dense `0..n` axis turned out to be three tables of five, not five.** The slippage-group
declaration and the substitution rate are sparse by the writer's own rule — a row exists only where
the run had something to say — and requiring a cover there refused the file a defaults run writes.

**C2 also owns refusing a file that parses and means nothing** — owner's decision, 2026-08-28,
because no step owned it and §9 promises "a malformed file fails at read with a line number".
`validate` runs after parsing and before the projection, and covers the constraints no shape can
state: a value outside its documented range (an inbreeding coefficient outside `[0, 1)`, a curve
weight outside `[0, 1]`, a substitution rate that is not a probability), a length spectrum whose
share count is even or below three or does not sum to one, an empty sample list, a contamination
table that is empty **or in which no row has a measurement** — the uncontaminated run written
longhand — and a measurement whose two evidence counts are both zero.

**Two of those already have a failing test waiting to be inverted**:
`the_shape_accepts_two_things_step_c2_must_refuse` pins that the last two are accepted today
(`parameters_file/mod.rs`), so landing the refusal flips an assertion rather than adding one
nobody remembered to write.
*Depends:* C1, B1. *Source:* §3, §6, §9.

✅ **C3. Float round-trip fidelity. Own commit, do not bundle.** Whether the `toml` crate emits
enough digits to recover every `f64` **has not been checked** (§4). Establish it; if it does not,
the fix is a serialiser that formats floats for round-trip, not a different format. **Oracle:** a
table of adversarial values — subnormals, values near the precision limit, the exact concentrations
a real fit produces — written and read back bit-identical.

**Landed 2026-08-30, and it found nothing.** Both formatters recover every double: twenty-two
adversarial values at four structural positions each, plus ten thousand pseudo-random bit patterns,
compared on `to_bits` rather than `==`. **Two of the brief's three adversarial categories are
covered and the third is not** — no value in the table comes from a real fit, which is C4's, and
the table carries four at a fit's *magnitudes* instead. **The spec's proposed fix turned out not to
fit either writer**: the one that writes the artefact already formats for round-trip by
construction, and the one that had genuinely not been checked writes a golden test file. **§4's
sentence saying this "has not been checked here" is now false and is the owner's to retire** —
recorded in `PROJECT_STATUS.md` rather than edited here. Mutating the writer's formatter to
`{value:.5?}` fails all four tests.
*Depends:* C1, C2. *Source:* §4, §13 test 1.

✅ **C4. The north-star round trip.** A `RunParameters` assembled from the joint fit on real tomato
records, written and read back, equal field for field: every float, every warrant, every count.
**This is goal 1 and the test the whole design rests on.**

**Landed 2026-08-30 on a run *shaped like* a fit's rather than on a real one**, by the owner's
ruling: no program in this tree produces a `RunParameters` from the joint fit — the joint walk
stops at a `JointFit`, every call to `assemble` is in a test module, and `assemble` needs the
minted-error totals the joint route does not fit. **⚑ The real-fit round trip is owed, and it is
owed a program rather than a test.** What is built instead goes through the fit's own doors —
`StratumFits::over` on `StratumOutcome`s and `assemble` on the raw per-read-group maps — with
three periods, 36 strata mixing fitted, derived and refused, two slippage groups, and the
substitution rate at the `(read group × stratum × ploidy)` grain.

**The review's Blocker was that the first draft of that fixture was not shaped like a fit's**:
within a period, the fitted and derived strata contradicted each other about whether the period had
a curve, and `LevelSource::Blend` — the fit's ordinary case — appeared nowhere in the file. Now the
first period has no curves and the other two have both, which is what decides which outcomes each
period's strata can have.

**It also measured what spec §9 prices and nothing had re-measured**: 157 bytes an inbreeding row
and 185 a substitution-rate row against §9's 146 each, which would put §9's 62 MB at 3,000 samples
nearer 79 MB — and both are floors, since a real cohort's read-group ids and bases-compared counts
are wider than this fixture's.
*Depends:* C3. *Source:* §1.2 goal 1, §13 test 1.

✅ **C5. The five states survive the round trip. Own commit, do not bundle.** One fixture per row of
§5, each built so that collapsing the two states it separates *changes an answer* — not merely so
they differ. An absent contamination table gives an uncontaminated run; a zero fraction with
non-zero counts gives a measured-and-clean read group; a defaulted scale of 1.0 does not read back
as fitted.

**Landed 2026-08-30, and the answer is not the same kind of thing for all five rows.** Two change
a number a locus is scored against — a stratum's length spectrum is the prior every tract of it is
seeded from, and a `(stratum × slippage group)` with no row sends the caller to the shipped stutter
model, which expects one read in twenty to come back a repeat short where a slip rate of zero
expects none. **Three change only what the run says about itself**, and that is the finding rather
than a weakness: a defaulted multiplier of 1.0 and a fitted one multiply by the same number, and
the three-term read likelihood at a fraction of zero equals the two-term one bit for bit. So the
warrant a call folds and the report an output prints are the only things that separate those
states — which is exactly what §5's own sentence says. Five mutations against the projection, one
a row, each failed its own row's test.
*Depends:* C4. *Source:* §5, §13 test 3.

> **Checkpoint C:** the round trip holds on a real fit and the five states are distinct through it.
> Pause for review. **Nothing after this milestone is worth reviewing until this one is green.**

### Milestone D — what the file is bound to

✅ **D1. The bindings, recorded.** The reference's content digest, the ordered sample list by name,
the read-group table, and the census recording terms the fit ran under — written by B1 and read by
C2.

**Landed 2026-08-30, and what was left was the two bindings that are not names.** The sample list
and the read-group table come off the run's own `ReadGroups` and were already written and read; the
other two are *derived*, and neither was derived anywhere — `of_run` took the reference digest as a
`&str`, its own documentation saying "nothing here can check that it is one", and took a
`CensusIdentity` nothing in the tree could build. **D2 cannot compare a file's bindings against a
run's until both sides come out of one function**, which is what `bindings.rs` is:
`CensusIdentity::of` over a census's `RecordingTerms` — twelve terms, named and ordered as
`RecordingTerms::first_disagreement` names and orders them — and one hex spelling.

**⚑ The review's Major was seven surviving mutants, every one a *part* of a value that is not a
scalar** — a block's contig, a stratum's period, the digest over all the kept loci. With the contig
undigested, two runs whose kept loci differ only in which chromosome a megabase sits on mint
byte-identical identities, which the fit calls unpoolable and D3 would not demote. Eight
part-by-part cases now, all 66 pairs for the order, and fourteen mutants each fail a test.

**The module's shared fixture changed twice, and both were the artefact teaching a reader something
false**: the reference digest was 16 hex characters where a reference digest spells as 32, and the
census carried one term where no run can produce fewer than twelve.
*Depends:* C4. *Source:* §3.1, §6.

✅ **D2. The three refusals.** A different reference, a sample list that does not match the run's,
a gap in the read-group ids. Each fails naming the field and the two values that differ, in the
shape the census's own refusal uses.

**Landed 2026-08-30, and the last clause of that sentence is the one that is wrong** — recorded in
`PROJECT_STATUS.md` and ruled by the owner on 2026-08-30: the census's own refusal is
`Freshness::{Rebuild, Refused}(&'static str)`, a field name with both compared values already
discarded, and §13 test 4 asks for the two values. **So this exceeds the census** rather than
matching it.

`refuse_if_not_this_runs_inputs` takes the same two arguments `of_run` writes from, so a file this
run wrote is a file this run accepts. **⚑ The review's Blocker was that it refused one anyway**:
the first draft joined the file's read-group table to the run's *by position*, where row order in
that table means nothing anywhere else in the module — so a file with two rows swapped validates,
projects, and is the same file. Joined on the read group's own number now.

**Three messages were rewritten because a reader could not act on them**, and one of the three was
found only by printing every message this code can produce and reading them: naming lanes by
`@RG ID` alone printed two identical lists beside *these differ* when only the numbering had moved.
`field` is a key path now — `fitted_from.read_groups[read_group = 1].library` — which is the
vocabulary `Meaningless` already uses.
*Depends:* D1. *Source:* §6, §13 test 4.

✅ **D3. The fourth binding demotes rather than refuses. Own commit, do not bundle.** A file fitted
from a different census of the same cohort is still interpretable, so every number in it is marked
`Supplied` and the run carries on. **Demotion is per-file, not per-number.** Its failure is silent:
the genotypes are identical and only the warrants move, so the test asserts both — same calls, every
warrant `Supplied`.

**Landed 2026-08-31**, with the tract ladder's bottom-rung warrant made carriable in its own commit
first, on the owner's ruling. `to_run_parameters_for` is the whole of §6 in one door: `validate`,
the three refusals, the census, the projection. **The demotion happens to the file** — five public
warrants, then one projection — which is the shortest way to be *per-file and not per-number*, and
the walk destructures without `..` so a sixth warranted number stops it compiling.

**⚑ "Every warrant `Supplied`" is the one clause of this step that cannot hold**, and the reason is
the owner's own ruling: the demotion is `weaker_of(file's warrant, Supplied)`, and `Provenance`
ranks `Supplied` above `Defaulted`, so assigning would *promote* every defaulted number into a
claim that somebody chose a value nobody chose. What the test asserts is **no warrant stronger than
`Supplied`, and none promoted**. Recorded in `PROJECT_STATUS.md`; §13 is the owner's to correct.

**⚑ And the door cannot be called by direct mode**, which has no census — the largest of three
questions this milestone hands back.
*Depends:* D2. *Source:* §2.1, §6, §13 test 5.

> **Checkpoint D:** a mismatched file is refused or demoted, and which one is never a surprise.
> Pause for review.

### Milestone E — the defaults, compiled in

✅ **E1. The defaults as named constants with their origin.** The base-quality calibration scale of
one, the repeat-tract outlier weight
([`DEFAULT_OUTLIER_WEIGHT`](../../../../src/ng/calling/likelihood/ssr.rs), 0.01, inherited and never
measured here), the flat concentration
([`STATED_FLAT_CONCENTRATION`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)), and
contamination's absence. Each marked `Defaulted` when used.

**Landed 2026-08-31, and the inventory is seven rows rather than these four.**
[`to_toml.rs`](../../../../src/ng/calling/parameters_file/to_toml.rs)'s own `origins` module —
written at B3 "so that step E1 has one list to reconcile against" — already carried the repeat-tract
substitution rate, the slippage numbers and the inbreeding coefficient. The last of those is the one
that mattered: the *fit* is forbidden a default for it, so a defaults run had nothing to write —
**the owner ruled on 2026-08-31 that the run takes zero**, which is what unblocked step E2. The
prior's seed is outside the table because what marks its fallback is
`SeedRegime::FallbackDiversity` and not a warrant, so no `validate` rung is possible for it.

**⚑ The step added a `validate` rung and had to take it out again**, and that is its most useful
finding. *A `defaulted` base-quality multiplier is the compiled-in constant and nothing else* — the
rule the outlier weight and the fallback concentration both carry — **refused a file this caller had
just written**. On those two, `defaulted` means the run took *that* number, so the warrant fixes the
value; on a multiplier it is copied from the **error rate** the multiplier was built from, and the
pre-pass's rate ladder has a `Defaulted` bottom rung of its own.
[`a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts`](../../../../src/ng/calling/parameters_file/from_run_parameters.rs)
stops it coming back.

**Ruled by the owner, 2026-08-31: that is the intended behaviour, and spec §5's third row is the
sentence to correct.** A library's real error rate is never its reported sequencing quality — the
quality scores describe base calling, and the reads also carry mismapping, chimeras and damage — so
a read group nothing could be fitted for is charged a stated rate rather than taken at its word, and
on a real library that is the conservative direction: at HG002's measured mean minted error of
2.9055 × 10⁻⁴ the multiplier is 3.44, 5.4 Phred less confident than the instrument claimed. The
0.001 is a placeholder and is to be fitted from GIAB, like §8's slippage numbers. Recorded in
`PROJECT_STATUS.md`; the spec is untouched.
*Depends:* C4. *Source:* §8.

✅ **E2. A run with no fit and no supplied file assembles `RunParameters` from the defaults.**

**Landed 2026-08-31.** `RunParameters::of_defaults` is the third door, beside `assemble` — which
takes the fit's raw outputs and derives the run's read-group axis from them, so a run with no fit
hands it empty maps and gets the wrong complaint — and `of_gathered_values`, whose nine arguments
are nine chances to leave a default out. It takes what a run with no fit has: its read groups, its
ploidy, and what it was told about inbreeding.

**The step turned on the owner's ruling of 2026-08-31**: a coefficient nobody stated is zero, and a
user may state one value for the whole run or a different value for any sample (`DeclaredInbreeding`,
joined by name). The *fit* is still forbidden a default for it; what separates the two is how far a
wrong constant travels — a fitted diversity divides by `1 − F` and carries the mistake into every
number the fit emits, where a defaulted coefficient reaches the calls and stops.

**The slippage group is declared and empty, which are two different things.** `StratumFits::at`
looks the read group up *before* the stratum, so declaring nothing would make every cell of every
tract answer `UnknownReadGroup` — *the run is not what it claims*, counted apart by
`TractScoringFits` — where a defaults run is exactly what it claims.

**⚑ The writer could not write such a run at all**, so two of its assertions were relaxed, narrowly:
the rates map covers every read group or is empty, and a missing rate is legal exactly where that
read group's calibration is `Defaulted`. **F1 inherits the third source** — a run scoring from a
supplied file has `Supplied` calibrations and no rate map, and `of_run` panics on it. Recorded in
`PROJECT_STATUS.md`.
*Depends:* E1. *Source:* §8.

☐ **E3. The slippage slot, and what a run does without it.** The per-(stratum × slippage group)
numbers have no compiled default until the GIAB measurement exists. **Trace what `StratumFits`
already does with a missing row** — the behaviour exists for partially-fitted runs and nobody has
looked — and make the gap visible in the run's report rather than silent. **Do not invent a
fallback**; if the traced behaviour scores a tract rather than refusing it, say so and stop for a
ruling.
*Depends:* E2. *Source:* §8, §12 question 1.

> **Checkpoint E:** a defaults run produces parameters, and says which of them were guessed.
> Pause for review.

### Milestone F — every run writes what it used

☐ **F1. One writer, three sources.** Supplied file, defaults, or fit — after assembly the run
cannot tell them apart, and writes the file beside its VCF unconditionally.
*Depends:* E2, B3. *Source:* §7.

☐ **F2. `RunParameterReport` becomes a view over the file.** It reaches only three things today and
its only callers are tests; the file is the superset, so the report reads from the same place rather
than being assembled separately.
*Depends:* F1. *Source:* §10 reuse map.

☐ **F3. A hand-written minimal file runs.** Someone writes the smallest file that calls a cohort —
no contamination, defaults where allowed — **from the spec alone, without reading this plan or the
code**. That is goal 3, and it is the only test here that cannot be automated.
*Depends:* F1. *Source:* §1.2 goal 3, §13 test 6.

> **Checkpoint F:** the file is a real input format, not an internal dump. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the five states of §5 are distinct in the types; a state that cannot be expressed is a compile error, not a test. **The key names are not verified here** — they are provisional by decision, and B2's first real output is what they get judged against |
| B | a fitted `RunParameters` writes a file that parses, with defaulted values carrying their origin as comments |
| C | **round-trip equality on a `RunParameters` from the real joint fit** (C4), plus adversarial float fidelity (C3) and one fixture per §5 row where collapsing two states changes an answer (C5) |
| D | each of the three refusals fires and names the differing field; the fourth binding gives identical genotypes with every warrant `Supplied` |
| E | a defaults run assembles and its report names every guessed number; the missing-slippage behaviour is traced and stated, not invented |
| F | a person writes a working file from the spec alone |

## Out of scope (next plans)

- **Direct mode itself** — [`run_driver_direct_mode.md`](run_driver_direct_mode.md), which consumes
  this file.
- **The GIAB slippage measurement** — a research note under
  [`../research/`](../research/), owed by §8.
- **`dump-parameters`** — the command surface's, trivial once F1 exists.
- **Whatever the VCF header prints from this** — the emission step's document, unwritten.
