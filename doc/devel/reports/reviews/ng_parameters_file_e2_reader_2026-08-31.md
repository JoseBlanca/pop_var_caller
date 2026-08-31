# Reading the parameters file a defaults run writes

I am a geneticist who runs variant callers. A colleague ran `ng` on a
tomato cohort with no parameter fit at all and gave me the parameters file
it wrote. I read it top to bottom with nothing else open. This is what I
could and could not get out of it.

**How the file was produced.** There is no checked-in copy of it, so I
built it: a throwaway `#[test]` beside
`a_defaults_run_writes_a_file_that_reads_back_as_the_same_run` in
`src/ng/calling/parameters_file/defaults.rs`, calling
`RunParameters::of_defaults` over the same three-lane / two-plant cohort
and printing `file.to_toml()`. I printed two versions: one where nobody
said anything about inbreeding (`DeclaredInbreeding::nothing_said()`),
which is the real "run with defaults", and one where `TS-1` was declared
at `F = 0.9`, to see what a declaration looks like. The probe is removed;
`cargo test --lib ng::calling::parameters_file` is back to 184 passed, 0
failed. The full produced file is beside this report as
`defaults_run_as_written.toml`; the parts I quote are pasted inline below.

Worktree: `/Users/jose/devel/pop_var_caller-e2-rev3` at `accd45a1` with
`tmp/e2_step.patch` applied. Nothing in the main repo was touched.

---

## The short answer

The file is well written wherever it is describing a fit. It is not
written for the case in front of me. Nine sections in, I still had not
been told the one fact that governs how I should read everything else:
**no fit ran**. I had to work that out by noticing that every warrant on
every line says `defaulted` and every table is empty.

Worse, in the two places where the file does speak to my situation, it
tells me things that are false:

- the section that carries my cohort's inbreeding coefficient says the
  caller *should not have been able to write that line*;
- the file's own editing rules say *do not mark an inbreeding coefficient
  `defaulted`* — and then every inbreeding row in the file is marked
  `defaulted`;
- the header block calls the whole thing `[fitted_from]`, "What these
  numbers were fitted from", above twelve digest lines. Nothing here was
  fitted from anything.

Two things I would have acted wrongly on: I would have believed my repeat
tracts were untouched by any stutter model (they were scored under
HipSTR's shipped constants, which the file never names), and I would have
filed a bug report about the inbreeding line rather than fixing my `F`.

---

## Blockers — I would act wrongly on these

### B1. The inbreeding row tells me the caller is broken, not that I have a knob to turn

This is what a defaults run writes for every sample, and it is the first
place I looked, because my cohort is a selfing tomato landrace and `F = 0`
is badly wrong for it:

```toml
[inbreeding]
# How inbred each sample is, as a fraction in [0, 1) — one row a sample, counted
# over the reference positions that sample's reads covered.
by_sample = [
    # no coefficient was fitted for this sample, and inbreeding has no default:
    # a run should not be able to write this line
    { sample = "TS-1", inbreeding_coefficient = { value = 0.0, warrant = "defaulted" } },
    # no coefficient was fitted for this sample, and inbreeding has no default:
    # a run should not be able to write this line
    { sample = "Ailsa Craig", inbreeding_coefficient = { value = 0.0, warrant = "defaulted" } },
]
```

**What I thought it meant.** Literally what it says: the program emitted a
line it believes it cannot emit. My first action would have been to stop
and report a bug, not to edit `0.0` to `0.9`. If I *had* edited it, I
would have done so believing I was working around a defect rather than
using a supported route.

**Why it's wrong.** The text is stale. `origins::INBREEDING_COEFFICIENT`
in `to_toml.rs:537` still carries the pre-ruling sentence, while
`of_defaults` now deliberately writes exactly this row for every sample.
The comment and the writer disagree about whether the row is legal.

**What would have told me.** Something in my language and about my plants,
e.g.:

> nobody said how inbred this sample is, so it is scored as an outcrosser
> (`0.0` = Hardy–Weinberg; the prior multiplies the heterozygote branch by
> `1 − F`). If your cohort selfs, say so — a landrace near `F = 0.9`
> scored at `0.0` has every homozygous stretch of its genome treated as a
> surprise. Set the value, set `warrant = "supplied"`, or declare it at
> run time.

Every clause of that already exists — in the Rust doc comment on
`DEFAULT_INBREEDING_COEFFICIENT`. None of it is in the file, which is the
only thing I have.

### B2. The file's own editing rules forbid the line the file just wrote

From the header, which I read first:

> `# On most of them there is no built-in number for it to be, so`
> `# writing one is a claim about this caller that no build makes — do not mark an`
> `# inbreeding coefficient or a substitution rate `defaulted`.`

Then `[inbreeding]` contains two rows marked `defaulted`. `validate()`
accepts them.

**What I thought it meant.** That the file I was handed was corrupt, or
had been hand-edited by my colleague. I would have asked him what he had
changed.

**What would have told me.** Either the rule needs the exception ("a run
that fitted nothing writes `defaulted` here and that is the one legitimate
case"), or the row needs a different warrant. As it stands, the one
sentence in the file that talks about my exact key says the key is wrong.

### B3. I cannot tell what stutter model my repeat tracts were scored under, and there is nothing to edit

My reads are PCR-amplified, so I expect more slippage at tracts than a
PCR-free library. Here is everything the file gives me:

```toml
slippage_by_stratum_and_group = []
...
length_spectrum_by_stratum = []
length_spectrum_by_period = []
...
substitution_rate_by_stratum = []
```

and, above them, this in the section prose:

> `# `slippage_by_stratum_and_group` is keyed by a stratum **and** a slippage`
> `# group, and a triple with no row means that group put no read in that stratum —`
> `# so one stratum can have a row for one group and none for another.`

**What I thought it meant.** That my reads never landed on a repeat tract,
so slippage never entered the calculation and nothing here concerned me. I
would have moved on.

**What is actually true.** Every cell falls back to
`StutterModel::hipstr_shipped` — one read in twenty a whole repeat short,
one in twenty a whole repeat long, plus one in a hundred each way for
part-repeat slips (`src/ng/alignment/stutter.rs:308`). Symmetric, taken
from another caller, fitted on no organism in particular. That is a
*stronger* claim about my data than "nothing happened", and it is the
exact number a PCR-library user needs to argue with: PCR stutter is
contraction-biased, so a symmetric 5%/5% is wrong for me in both magnitude
and shape.

The `[repeat_tracts]` prose is nearly a page long and never says the word
"HipSTR", never gives 0.05, and never distinguishes *no reads reached this
stratum* from *no stratum was ever fitted, so all of them fall back*. The
empty-array case is not the same as the missing-row case, and the prose
only covers the second.

Nor can I change it. I would have to reverse-engineer the row schema from
a fitted file I do not have, and the header tells me slippage rows carry
no warrant so there is nowhere to record that I did:

> `# The slippage numbers, the prior's two concentrations and the length spectrum`
> `# rows carry no warrant — they say where they came from another way, and there`
> `# is nowhere in them to record that you changed one. Note such an edit`
> `# elsewhere.`

"Note such an edit elsewhere" is not an instruction I can follow — elsewhere
where, and who reads it?

**What would have told me.** One line above the empty array:

> no stratum was fitted, so every repeat tract in this run was scored
> under HipSTR's shipped constants: 5 reads in 100 report a whole repeat
> short, 5 in 100 a whole repeat long, symmetric. PCR libraries usually
> slip more than this and slip short more often than long.

The module header of `defaults.rs` says this almost word for word. The
file does not.

Same shape, smaller stakes: `substitution_rate_by_stratum = []` means
every cell takes `DEFAULT_SSR_SUBSTITUTION_RATE` = 0.001 at the tract. The
file gives me an empty array and no number.

### B4. `[fitted_from]` claims a fit that did not happen

```toml
[fitted_from]
# What these numbers were fitted from. A run whose reference, samples or read
# groups do not match these is refused; one whose census does not match keeps
# the numbers and reports every one of them as supplied rather than fitted.
```

followed by a reference digest, the sample list, the read-group table, and
twelve census digest lines with names like `"the loci actually kept"`,
`"per-stratum locus counts"`, `"depth ladder edges"`, `"per-position depth
cap"`.

**What I thought it meant.** That a fit had run over those loci at those
depths, and that the numbers below came out of it. This is the second
section of the file and it set my reading of everything after it. It is
why it took me until `[stated_constants]` — the last section — to be
confident nothing was fitted at all.

**What would have told me.** The section is doing real work (it is what
lets a later run refuse a mismatched reference), so it should stay — but
under a heading that is true when the count of fitted numbers is zero:
"what this run was pointed at", with a line saying no numbers were fitted
from it.

---

## Majors — I could not act on these

### M1. Nothing in the file says "this run fitted nothing"

The single most important fact about the file I am holding is not stated
anywhere in it. It is only recoverable by reading all 245 lines and
noticing that every `warrant` says `defaulted`, every array is empty, and
`rung = "stated_heterozygosity"`. A one-line summary near the top —
"0 of 9 parameter groups were fitted from your data" — would have saved
the whole exercise, and would be exactly as useful on a *partial* fit,
which is the common case.

### M2. The contamination explanation disappears in the file that needs it most

The fitted example file explains, at the point of the section, the state I
am actually in:

> `#   - this whole section absent  -> nobody identified any contamination`
> `#   - a row with no `measurement` -> this lane could not be measured`
> `#   - a zero share with non-zero counts -> measured, and found clean`

My file has no `[contamination]` section, and therefore does not have that
text either. So the reader who most needs to be told what an absent
section means is the only reader who is not told. The only trace is the
generic header line "An absent key is not a zero", which does not tell me
which absence I am looking at or that contamination is now being modelled
as exactly zero.

I would have concluded either "contamination was checked and was fine" or
"contamination is not something this caller does". Both are wrong; the
truth is "you were not asked, so your reads are scored as clean".

### M3. Declaring `F` produces a row that breaks the file's own rule

In the run where I declared `F = 0.9` for `TS-1`:

```toml
{ sample = "TS-1", inbreeding_coefficient = { value = 0.9, warrant = "supplied", observations = { covered_positions = 0 } } },
```

The header says:

> `# **If you edit one, change its warrant to "supplied" and delete its `observations`**`

and

> `# An absent key is not a zero. ... a zero means it was measured and found to be zero.`

So the program writes, for a number I typed, an `observations` count of
zero — which by the file's own definition means *measured over zero
covered positions*, and which by the file's own editing rule I am
supposed to have deleted. Three sentences, no consistent reading. I could
not tell whether the row was fine or whether I had done something wrong.

### M4. The file records no version, no date, no command

`format_version = 1` is the format's version, not the caller's. There is
no run date, no `ng` version, no record of the command line or of the flag
that selected "run with defaults". For a file whose stated purpose is to
say what a run used so a later run can reproduce it, that is the first
thing I looked for and did not find. If two of these land in my project
directory I cannot tell which run wrote which, or which build's constants
they contain — and the header itself concedes the constants can move
("`DEFAULT_ERROR_RATE` is itself a placeholder", per the module docs).

### M5. `[ordinary_site_prior]` is honest about being human but gives me no way to act

```toml
reference_concentration = 1.0
alternative_concentration_total = 0.001
rung = "stated_heterozygosity"
```

The prose is good — it names `stated_heterozygosity` as "a stated
heterozygosity taken from human data ... which is the one that rests on
nothing this run measured". That is the clearest sentence in the file.

But: I am calling tomato. `0.001` is 10 differences per 10,000 bases, a
human figure; this project's own tomato panel sits at about 6 per 10,000
(`ExpectedHeterozygosity::SPECIES_FALLBACK` doc). The file gives me the
number in Dirichlet-concentration units, never in differences per
kilobase, so I cannot check it against anything I know about my species
without knowing that `alternative_concentration_total` *is* the expected
heterozygosity. And it carries no warrant, so per the header there is
"nowhere to record that you changed one".

I did check: the reader does take these two floats literally
(`to_run_parameters.rs:163`), so editing works. The file should say so.

---

## Minors — I had to re-read them

### m1. The header's editing rule mostly does not apply here

> `# **If you edit one, change its warrant to "supplied" and delete its `observations`**`

In a defaults file no value has `observations`, so half the instruction is
inert. Not wrong, just spent effort.

### m2. The multiplier paragraph in the header describes something my file does not contain

> `# On`
> `# `base_quality_calibration.by_read_group[...].error_probability_multiplier``
> `# there is a built-in number, 1.0, and the key still is not checked: `defaulted``
> `# there is copied from the error rate the multiplier was built from, and a run`
> `# whose rate itself was defaulted writes a `defaulted` multiplier that is not`
> `# 1.0.`

My three rows are all exactly `1.0`. I re-read this twice trying to work
out whether mine were the good case or the bad one. The paragraph is about
a situation (a fit that ran and could not measure one read group) that a
defaults file never contains.

### m3. `[sequencing_batches]` writes explicit tables under a flag that says nobody chose them

`batching_was_declared = false`, then a full `by_read_group` and
`by_sample` table assigning everything to batch 0. The prose does explain
the flag, and it explains it well. But the tables still read as a
determination somebody made, and on a first pass I took them as one.

### m4. `slippage_group_by_read_group` cannot be told from an operator's declaration

```toml
slippage_group_by_read_group = [
    { read_group = 0, slippage_group = 0 },
    { read_group = 1, slippage_group = 0 },
    { read_group = 2, slippage_group = 0 },
]
```

The prose says "the run declares it". I have three lanes across two
libraries and would very much like to know whether somebody decided they
slip alike or whether that is the fall-through. `[sequencing_batches]`
solves the same problem with `batching_was_declared`; this table has no
equivalent.

---

## What read well

Genuinely, and worth keeping:

- **The base-quality multiplier gloss.** "Above one says the instrument
  was optimistic and the reads are worse than it claimed; below one says
  they are better; one leaves the qualities exactly as reported. It is not
  a multiplier on the Phred score, which moves the other way." I knew what
  the number meant and which way to push it, immediately. This is the
  model the rest of the file should follow.
- **The fallback concentration row.** `fallback_length_spectrum_concentration
  = { value = 1.0, warrant = "defaulted" }` with "this run fitted no
  stratum on its own tracts, so there was no median to take and this is
  the caller's own constant" above it. It is the one place in the whole
  file that says *your run fitted nothing here, and here is the constant
  that stood in*. That sentence, repeated over the inbreeding rows and the
  empty slippage array, would fix B1 and B3.
- **The census "do not edit these to make a run match" paragraph.** It
  anticipates the exact thing an impatient user does, and says what it
  buys them.
- **`rung = "stated_heterozygosity"` and its explanation** — see M5. The
  honesty is there; only the units are missing.
- **The `-span..+span` warning on `shares_by_repeat_offset`.** "An array
  read as starting at zero shifts every length this stratum expects by one
  repeat." That is a real off-by-one someone would hit, named before they
  hit it.

---

## Part 3 — the module header of `src/ng/calling/parameters_file/defaults.rs`

The brief says this prose is meant to tell somebody in my position which
of the caller's numbers are guesses and how badly. Judged that way, it is
addressed to the wrong reader: it is a design record with a user-facing
table wedged in the middle. Below, only the things that hit me as a
geneticist.

### Sentences I could not act on

- **"a default taken at the tract, not written in the file"** (the
  substitution-rate row of the seven-row table). I do not know what "at
  the tract" means as a place where a default is taken. Three readings
  seemed possible.
- **"**So `defaulted` is not a warrant this file's substitution-rate rows
  can legitimately carry**, and nothing checks that they do not."** I am
  told a bad state exists and is unguarded. There is no action for me. And
  it sits three paragraphs from a table row that says the inbreeding
  default *is* legal — which is the pair I then have to hold in my head
  when I look at my own file and find `defaulted` on the inbreeding rows.
- **"§3.5 requires at least one row and a defaults run now has one for
  every sample."** I have no §3.5. This is the closing sentence of the
  paragraph I most needed, and it lands on a document I do not have.
- **"Recorded in `PROJECT_STATUS.md`."**, **"Both halves are recorded at
  Checkpoint B."** Same.
- **"what a run can look at is how its repeat-tract calls move when the
  number is edited"** (the outlier weight). That is a sensitivity analysis
  with no guidance on what a concerning movement looks like or what value
  to try instead. It is honest, but it is not something I can do on Friday
  afternoon.

### Terms used before they are defined

`warrant`, `stratum`, `slippage group`, `rung`, `census`, `the projection`,
`the door`, `minted` (as in "the instrument minted"), `provenance`. Some
are defined later; several are never defined here at all. `stratum` in
particular does load-bearing work in the table before the file's own
definition ("a class of tract, spelled as `period` and `reference_repeats`")
appears anywhere the reader has been sent.

The bare section numbers — §3.7, §3.8, §5's first row, §5's third row,
§8's third bullet, §12 question 1 — appear about fifteen times. Each is a
claim I am asked to accept on the authority of a document I was not given.

### Claims that read wrong to a geneticist

**The multiplier-is-conservative paragraph.**

> "A library's real error rate is never its reported sequencing quality:
> the quality scores describe base calling, and the reads also carry
> mismapping, chimeras and damage. So a read group the fit could not
> measure is charged a stated rate rather than taken at its word, and on
> any real library that pushes the reads the *conservative* way — on
> HG002's mean minted error of 2.9055 × 10⁻⁴ the multiplier is 3.44, 5.4
> Phred less confident than the instrument claimed."

The arithmetic checks (0.001 / 2.9055e-4 = 3.44; 10·log₁₀3.44 = 5.4 dB).
Two objections:

1. **Mismapping does not belong in a base-quality multiplier.** Mapping
   confidence is a per-read property carried by MAPQ, and this caller
   uses MAPQ elsewhere. Folding mismapping into a base-error scale either
   double-counts it or, worse, spreads a whole-read property across bases
   as if independent — which understates its effect at a site where a
   whole misplaced read supports one allele. The sentence should say
   *residual base-calling error the instrument's model does not capture*,
   and drop mismapping.
2. **"Conservative" is doing work a number should do.** 3.44 is
   conservative for a library at Q35. It is *anti*-conservative for a
   library whose true error rate is above 10⁻³ — an older chemistry, a
   degraded sample, a long-read run — where charging a flat 0.001 makes
   the reads look *better* than they are. The claim "on any real library
   that pushes the reads the conservative way" is stated without its
   range. Per this project's own standing rule about the range rather than
   the example, it needs the boundary: above what minted error does the
   flat rate stop being conservative?

**"one read in twenty comes back a whole repeat short and one in twenty a
whole repeat long".** Correct against
`StutterModel::hipstr_shipped` (0.05 / 0.05), and the contrast with
HipSTR's contraction-biased *fitted* values is exactly the right thing to
tell me. Two gaps: it omits the part-repeat shares (0.01 each way), and
this excellent paragraph is in the source file rather than in the TOML,
which is the only artefact I have (B3).

**"one chromosome's worth of belief spread flat".** I read "chromosome"
literally on the first pass and could not make it mean anything. It means
one prior observation in Dirichlet units. Say "as much weight as a single
observed chromosome", or just "one pseudo-count".

### The inbreeding paragraph — the fit may not default it while the run may

Quoted in full because this is where the argument breaks:

> "The **fit** may not default it: `parameter_estimation::generic::fallback`'s
> header states the rule — *"The inbreeding coefficient has one rung and it
> is not a default … it is the parameter that differs most between an
> outcrosser and a selfing landrace, and a cohort's diversity divides by
> `1 − F`, so a wrong constant would be amplified rather than absorbed"* —
> and that stands. **The run may**, by the owner's ruling ... The two are
> not in conflict because they are different acts: the fit declines to
> *infer* a coefficient from data that cannot carry one, where the run
> *declares* what it believes about its plants."

**The distinction does not hold for the case that produced my file.** My
colleague ran with `nothing_said()`. Nobody declared anything. The
paragraph calls the resulting `0.0` "what the run *declares* about its
plants in the absence of any statement" — but a declaration nobody made is
not a declaration; it is exactly the default the fit was forbidden to
take. The prose renames the act and calls the conflict resolved.

**And the stated reason for the prohibition applies to calling too.** The
argument against defaulting in the fit is that "a cohort's diversity
divides by `1 − F`, so a wrong constant would be amplified rather than
absorbed". The rebuttal offered is "a run that fitted nothing computes no
diversity". True — and irrelevant to the harm I actually take. In the
genotype prior the heterozygote branch is multiplied by `1 − F`. At the
true `F = 0.9` of my landrace, scoring at `F = 0` inflates that branch
tenfold. That is amplification of a wrong constant, in the calls, which is
the output I keep. The paragraph never mentions it. Its own doc comment on
`DEFAULT_INBREEDING_COEFFICIENT` does — "a landrace at `F = 0.9` scored at
zero is told every homozygous stretch of its genome is a surprise" — so
the module contains the counter-argument to its own rationale, two
screens apart, and does not reconcile them.

**"Zero is the value at which the genotype prior does nothing"** is the
sentence I would push back on hardest. `F = 0` is not the absence of an
assumption; it is Hardy–Weinberg, which is a specific and strong claim
about a population — random mating, no selfing, no structure. It is the
right *arithmetic* analogue of a multiplier of one, and the module says so
carefully. But the module's own multiplier paragraph already makes exactly
this point about the multiplier — "A multiplier of one declines to
recalibrate; it does not abstain from a claim ... which asserts the
instrument was right" — and then does not make it about `F`. On a selfing
crop, a plant sequenced at 3× and scored under Hardy–Weinberg will be
called heterozygous where it is homozygous, and that is the same class of
error the multiplier paragraph is careful about.

**What is missing entirely: what happens at one sample.** A single sample
has no cohort to estimate `F` from, so it will *always* land on this
default, and the genotype prior's `1 − F` is doing full-strength work
there. The paragraph discusses the fit's refusal and the run's declaration
and never says what a single-sample run gets. Given this project's own
rule that a method must state what it does at one sample, that is the
sentence I most wanted.

**What would have made the paragraph work.** Drop the "different acts"
framing, which does not survive `nothing_said()`, and replace it with the
honest version:

> A run that fitted nothing has to put *some* coefficient in the genotype
> prior, and zero is the one that leaves Hardy–Weinberg in place. That is a
> real assumption, not an abstention, and it is wrong in a known direction
> on any selfing or structured cohort: heterozygotes are over-called. The
> fit is still forbidden from taking it, because a fitted diversity divides
> by `1 − F` and would carry the error into every downstream number,
> whereas a defaulted call carries it only into the calls. The file marks
> every such coefficient `defaulted`, and an operator who knows the crop
> should say so.

---

## Where the fixes are

Every Blocker above is a comment string, not a redesign:

| finding | file to change |
|---|---|
| B1, B2 | `to_toml.rs:537` — `origins::INBREEDING_COEFFICIENT`, and the header's "do not mark an inbreeding coefficient ... `defaulted`" clause |
| B3 | `to_toml.rs` — a note on the empty `slippage_by_stratum_and_group` / `substitution_rate_by_stratum` arrays naming the fallback and its numbers |
| B4, M1 | `to_toml.rs` — `[fitted_from]`'s heading note, plus a "0 of N fitted" line near `format_version` |
| M2 | `to_toml.rs` — emit the three-state contamination note even when the section is absent |
| M3 | `defaults.rs` — `DeclaredInbreeding::of_each_sample` writes `observations: 0` on a `Supplied` row; that count should not be written |
