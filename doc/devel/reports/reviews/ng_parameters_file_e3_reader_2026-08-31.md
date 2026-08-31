# Reading the parameters file a no-fit run wrote — PCR-amplified tomato cohort

I run variant callers on tomato. My libraries are **PCR-amplified whole-genome
Illumina** (a handful of amplification cycles during library prep — not targeted
STR amplification), so I expect more stutter at repeat tracts than a PCR-free
library would give, and I expect it to run **short** far more often than long.

I called the cohort with `ng` with **no parameter fit at all**, and this is the
file the run wrote. I have read it top to bottom. Line numbers below are the
file's own, as produced (264 lines total).

I also diffed it against
`doc/devel/reports/reviews/ng_parameters_file_e2_defaults_run_as_written_2026-08-31.toml`,
the same file before this change, so I know which sentences are new.

---

## What I read — the repeat-tract section as produced

```toml
[repeat_tracts]
# Everything about repeat tracts. A **stratum** is a class of tract, and every
# row keyed by one spells it as `period` — how many bases one repeat unit is —
# and `reference_repeats` — how many copies of it the reference carries. A
# **slippage group** is a set of read groups whose reads are taken to slip
# alike; the run declares it, and `slippage_group_by_read_group` below is that
# declaration.
#
# `slippage_by_stratum_and_group` is keyed by a stratum **and** a slippage
# group, and a triple with no row means that group put no read in that stratum —
# so one stratum can have a row for one group and none for another. A stratum
# with no row in `length_spectrum_by_stratum` was never fitted on its own tracts
# and falls to its period's pooled one, or to the flat shape below. Neither
# absence is a zero.
#
# Three numbers a stratum: `share_of_reads_that_slip` — how often a read reports
# a tract length other than its allele's; `shorter_share` — of the reads that
# slip, the share showing a shorter tract; `fall_off` — how fast two-repeat
# slips fall off against one-repeat slips. `expected_slipped_reads` is
# fractional because it is how many reads the fitted share says slipped, not a
# count anybody labelled.
#
# `share_of_reads_that_slip_origin` says where the first of the three came from
# and `shorter_share_and_fall_off_origin` where the other two did: this
# stratum's own fit, its period's curve, or a blend of the two, with the curve
# itself written down so an interpolation can be told from a measurement. A row
# whose two shares were not fitted here at all has no
# `shorter_share_and_fall_off_origin` key. Each origin carries its own
# `expected_slipped_reads` where this stratum fitted a slip share of its own,
# and neither carries one where the number was taken whole from a curve — so a
# row showing the same count twice is not a duplicate, and a row showing none
# fitted nothing of its own.
#
# The curves under `shorter_share_and_fall_off_origin` also record
# `curve_fitted_on`, which says what that curve itself was fitted on: this
# period's own strata (`this_period`), or those same strata where there were too
# few to score the shape (`this_period_unscored`), or the other periods pooled
# (`other_periods`), or a stated constant where no period had anything to fit
# (`built_in_default`). The curve under `share_of_reads_that_slip_origin` is a
# different fit and has no such key.
#
# `fallback_length_spectrum_concentration` is what a tract falls back to where
# neither its own stratum nor its period was fitted: this many chromosomes'
# worth of belief, spread flat over whatever lengths the tract offers, so a
# larger number makes the prior harder for the reads to move. A run that fitted
# any stratum states the median of the concentrations those fits produced and
# marks it `fitted_here`; a run that fitted none states a built-in constant and
# marks it `defaulted`; a run handed this file by another run marks it
# `supplied`, and so should you if you change it. It carries no `observations`
# in any of those cases — a median over strata is not a measurement with a
# sample size.
# this run fitted no stratum on its own tracts, so there was no median to take
# and this is the caller's own constant
fallback_length_spectrum_concentration = { value = 1.0, warrant = "defaulted" }
slippage_group_by_read_group = [
    { read_group = 0, slippage_group = 0 },
    { read_group = 1, slippage_group = 0 },
    { read_group = 2, slippage_group = 0 },
]
# This table is empty, which is not the same as a missing row: no stratum was
# fitted at all, so every repeat tract in this run was scored under the stutter
# model this caller ships — 5 reads in 100 report a whole repeat short and 5 in
# 100 a whole repeat long, with 1 in 100 short and 1 in 100 long for a part
# repeat. Those are HipSTR's shipped constants, not a fit, and they are
# symmetric where real slippage is usually biased towards the shorter tract. A
# PCR library slips more than this. There is nothing in this file to edit: fit
# the run, or read the calls at repeat tracts as resting on somebody else's
# chemistry.
slippage_by_stratum_and_group = []
# `shares_by_repeat_offset` runs from -span to +span in whole repeat units from
# the **reference** tract length, so the middle entry is the reference length
# itself, the count is odd, and the shares sum to one. An array read as starting
# at zero shifts every length this stratum expects by one repeat.
length_spectrum_by_stratum = []
length_spectrum_by_period = []
# How often a base reads wrong inside a tract — per read group as well as per
# stratum, because that is a property of the chemistry, and per `ploidy` as
# well, because that is the set of genotypes the fit scored these tracts
# against. That third key is why a row repeats the number at the top of the
# file: a cohort called at two ploidies carries a row for each. Counted in bases
# compared, not reads: a read crossing a tract contributes as many bases as it
# crosses.
# This table is empty, which is not the same as a missing row: nothing was
# fitted for any read group at any stratum, so every tract's cells took the
# caller's stated 0.001 — about one base in a thousand read wrong inside a
# tract. Base quality inside tracts is usually worse than outside them, so on
# real reads that number is likely optimistic.
substitution_rate_by_stratum = []

[stated_constants]
# The numbers no fit produces, written out so that what a run inherited is
# visible and editable rather than buried in the binary.
# inherited from the existing caller and never measured here: the share of
# repeat-tract reads that came from nowhere the model can explain
repeat_tract_outlier_weight = { value = 0.01, warrant = "defaulted" }
```

Two comment blocks are new (lines 228–236 and 251–255). Everything else in this
section is unchanged from the previous copy.

---

## What read well — keep these sentences

Say this first because most of what follows is complaint, and the complaints are
about a note that is a large improvement on nothing.

1. **"This table is empty, which is not the same as a missing row"** (line 228).
   This is the sentence that fixed the reading. Before it, I would have read
   `slippage_by_stratum_and_group = []` against the paragraph at line 177 and
   concluded *my reads never landed on a repeat tract*, which for a tomato
   genome is nonsense but is exactly what the paragraph licenses. Keep this
   clause verbatim.

2. **"Those are HipSTR's shipped constants, not a fit"** (line 232). Naming the
   other caller is right. It tells me where to go and read for myself, and
   "not a fit" is the two words that stop me treating the number as a
   measurement of anything.

3. **"they are symmetric where real slippage is usually biased towards the
   shorter tract"** (lines 232–233). This is the single most useful sentence in
   the section for me, because contraction bias is the thing I know about my own
   chemistry and the thing that would make the calls wrong. Keep it.

4. **"read the calls at repeat tracts as resting on somebody else's chemistry"**
   (lines 235–236). Blunt and correct. That is the disclaimer I will paste into
   the methods section if I publish these calls.

5. Numbers written as **"5 reads in 100"** rather than 0.05. I can hold that
   against what I know. Do not convert these to decimals.

6. In the substitution note, **"about one base in a thousand read wrong inside a
   tract"** (lines 253–254) — the gloss on `0.001` is what makes the raw number
   mean something. Keep the gloss; it is the pattern the rest of the file should
   follow.

---

## 1. Can I tell what stutter model my tracts were scored under, and can I compare it to PCR?

**What model: yes.** Line 229 tells me plainly, and line 232 names its source.
This is a clear improvement — before this change the empty table said nothing at
all.

**In a form I can compare with PCR: no, not without arithmetic and not against
the file's own columns.** Three separate problems.

### 1a. The note quotes the model in a different vocabulary from the table it sits on — **Major**

> "Three numbers a stratum: `share_of_reads_that_slip` … `shorter_share` … `fall_off`" (lines 184–187)
>
> "5 reads in 100 report a whole repeat short and 5 in 100 a whole repeat long, with 1 in 100 short and 1 in 100 long for a part repeat" (lines 230–232)

Forty-four lines apart, the file teaches me one parameterisation and then
describes the defaults in a different one. **I cannot line the note up against
the table it is attached to**, and I cannot line it up against a fitted run's
rows later, which is the comparison I actually need — "here is what the defaults
claimed, here is what my data says".

What I would have wanted, in the file's own three words:

- `share_of_reads_that_slip` = **12 in 100** (0.05 + 0.05 + 0.01 + 0.01)
- `shorter_share` = **50 in 100** — dead even, which is the whole complaint
- `fall_off` — 95 in 100 of slips are a single step

I had to compute the first two myself from four numbers. **The 12-in-100 figure
is the one I compare against PCR**, and it is the one number the note does not
give me. My PCR+ WGS libraries are the mild case, but 12% of reads misreporting
tract length is at the low end for anything mononucleotide and roughly plausible
for a short tetranucleotide — I can only say that after doing the addition.

*What would have told me:* one extra sentence — "In this section's own three
numbers that is `share_of_reads_that_slip` = 0.12, `shorter_share` = 0.50,
`fall_off` = 0.95."

### 1b. "a whole repeat short" is ambiguous, and the two readings differ — **Minor**

I read "5 reads in 100 report a whole repeat short" as *5 in 100 come back
exactly one repeat unit shorter*. That is not what 0.05 is. In
`src/ng/alignment/stutter.rs` the 0.05 is the **direction share** — all
whole-repeat contractions of any size — and the one-step share of 0.95 then
splits it. P(exactly one repeat shorter) is 0.0475, or 4.75 in 100.

The size of the error is small and would not change what I do. But it is the
kind of error I cannot detect from the file, and the fix is one word: "a whole
repeat *or more* short", or "shorter by a whole number of repeats".

### 1c. "part repeat" is never defined anywhere in the file — **Major**

> "with 1 in 100 short and 1 in 100 long for a part repeat" (lines 231–232)

I grepped the whole 264-line file. **"part repeat" occurs exactly here and
nowhere else.** I do not know what it means. My guesses were (a) a tract whose
reference length is not a whole number of units, (b) a read that ends inside the
tract, (c) an indel inside the tract that is not a multiple of the period. It
turns out to be (c) — a sequencing indel or an interruption, sized in base pairs
— which I only learned by reading the Rust source, which I should not have to do.

Two of the four numbers in this note are therefore numbers I cannot interpret.

*What would have told me:* "…and 1 in 100 each way for a **part-repeat** change —
an insertion or deletion inside the tract that is not a whole number of units, so
a sequencing indel or an interruption rather than slippage."

### 1d. The note does not say the same four numbers are used for every period and every tract length — **Major**

> "every repeat tract in this run was scored under the stutter model this caller ships" (line 229)

I read "every repeat tract" as scope — *no tract escaped this*. What I did not
take from it, and what matters far more to me than PCR-versus-PCR-free, is that
**the same 5/5/1/1 is applied to a 25-base poly-A run and to a 5-copy
tetranucleotide alike**. Slippage rises steeply with period going down and with
tract length going up; a mononucleotide run of 20 A's in a PCR library stutters
several times as hard as 12 in 100. The whole point of the `stratum` key — the
file's own first paragraph defines a stratum as `period` × `reference_repeats` —
is that these numbers are supposed to differ across it, and here a single flat
pair is standing in for the entire surface.

**This is the finding that would change my behaviour**: I would drop
mononucleotide and dinucleotide calls entirely rather than merely caveat them,
and the file did not tell me to.

*What would have told me:* "One pair of numbers stands in for every stratum, so a
long mononucleotide run and a short tetranucleotide are scored identically — real
slippage differs several-fold across that range, and short-period long tracts are
where this is furthest wrong."

---

## 2. What would I do? Is there an action I can take?

The note names two options:

> "There is nothing in this file to edit: fit the run, or read the calls at repeat tracts as resting on somebody else's chemistry." (lines 234–236)

**Option two I can take** — and it is what I will do today. I will treat every
repeat-tract call as provisional, and I will not report STR allele frequencies
for this cohort.

**Option one I cannot take, because the file does not tell me how.** — **Major**

"Fit the run" is not an instruction. It does not tell me what flag to pass, what
extra input a fit needs (more coverage? a truth set? a second pass over the
BAMs?), how long it takes, or whether my 63 accessions at about three reads a
position have enough evidence to fit slippage at all. **That last is the question
I would actually ask**, and I suspect the answer at 3× is "not much", in which
case the file is pointing me at a door that is closed for my data — and should
say so.

*What would have told me:* the flag, plus one clause on what a fit needs.
Something like: "fit the run (`ng call --fit-parameters`, a second pass over the
same alignments; strata with too few reads still fall back, so a thin cohort will
fit some periods and not others)". Even naming the flag alone would move this
from Major to Minor.

**A smaller conflict in the same breath** — **Minor**. "There is nothing in this
file to edit" (line 234) sits 29 lines above:

> "The numbers no fit produces, written out so that what a run inherited is visible and **editable** rather than buried in the binary." (lines 259–260)
> `repeat_tract_outlier_weight = { value = 0.01, warrant = "defaulted" }` (line 263)

`repeat_tract_outlier_weight` is a repeat-tract number, it is editable, and the
file invites me to edit it. I stopped and re-read to work out whether the two
sentences disagreed. They do not — "nothing to edit" scopes to the slippage
table — but the scope is not on the page. "There is nothing **in this table** to
edit" fixes it.

---

## 3. Does the new note sit well with the long-standing missing-row paragraph?

**Reading in the order the file presents them: not quite. The correction arrives
50 lines after the misreading it corrects.** — **Minor**, but a cheap fix.

At line 177 I read:

> "`slippage_by_stratum_and_group` is keyed by a stratum **and** a slippage group, and a triple with no row means that group put no read in that stratum … Neither absence is a zero."

At that point I form a complete and wrong mental model: *absence here means my
reads did not cover it*. I then read 45 lines about origins, curve provenance and
fallback concentrations. By line 228, when the note arrives, I have moved on.

The note does the right thing — it names the distinction explicitly ("which is
not the same as a missing row") rather than assuming I remember — so I got there.
But I got there by being corrected, not by never forming the wrong model. The
paragraph at 177 never hints that an empty table is a separate case at all.

*What would have told me:* one clause at line 182, where the wrong model is
formed — "Neither absence is a zero. **An empty table is a third case and a
stronger claim; the note at the table itself says what.**" Then the note at 228 is
confirming something I am already expecting.

**They do not conflict.** They overlap only in the phrase "missing row", and the
new note uses it correctly to distinguish itself.

**One inconsistency I noticed while reading in order** — **Minor**. Four tables
in this section are empty. Two now carry a note; `length_spectrum_by_stratum` and
`length_spectrum_by_period` (lines 242–243) carry none. Having just been told
twice that an empty table is a stronger claim than a missing row, I paused at
those two and wondered whether their emptiness meant something *different*. It
does not — it is the same "nothing was fitted". The information is technically
present, in the `fallback_length_spectrum_concentration` note at lines 220–221
("this run fitted no stratum on its own tracts"), but that is attached to a
different key two lines above and I did not connect it. Either give those two a
one-line note as well, or say in the new note that the spectra are empty for the
same reason.

---

## 4. Is the substitution-rate note useful, or noise beside the slippage one?

**Useful, but it is currently costing the slippage note some of its force** —
**Minor**.

It is genuinely worth having: it turns `0.001` into "about one base in a thousand
read wrong inside a tract", and "base quality inside tracts is usually worse than
outside them, so on real reads that number is likely optimistic" is a fair and
correctly-hedged statement. As a geneticist I have no quarrel with it — Q30
inside a tract is optimistic but not wildly so, and the dominant tract error is
indel slippage, which the *other* note covers.

The problem is presentation. **The two notes open with the identical 13-word
clause** — "This table is empty, which is not the same as a missing row" — 23
lines apart. My eye read the second as a repeat of the first and nearly skipped
it. And the stakes are very different: the slippage note is telling me my STR
calls rest on another organism's chemistry, the substitution note is telling me
one base in a thousand is a slightly rosy estimate. Identical framing implies
equal weight.

*What I would do:* keep it, but vary the opening and shorten it. "Empty for the
same reason: nothing was fitted, so every tract's cells took the caller's stated
0.001 — about one base in a thousand read wrong inside a tract. Tract-internal
base quality is usually worse than that, so it is likely optimistic."

**One contradiction it walks into** — **Major**. The file's preamble says, at
lines 22–26:

> "A `defaulted` substitution rate is the one to avoid writing by hand: **nothing here defaults that number in this file**, so the warrant would be a claim about this caller that no build makes."

And then at line 252:

> "every tract's cells took **the caller's stated 0.001**"

Read in order, those say the caller has no default substitution rate, and then
that the caller has one and used it on every tract. I re-read both twice and
concluded one of them was a bug. It is not — the reconciliation is that the
default is applied *at the tract* and never written as a row here, so no row in
this file may legitimately carry a `defaulted` warrant — but the preamble's
phrasing "nothing here defaults that number" does not carry that meaning to a
reader. Note this is a new collision: this preamble sentence is itself new in
this change, and the note it collides with is new too.

*What would have told me:* preamble — "nothing **writes a row for** that number
in this file (the default is taken at the tract instead), so a `defaulted`
warrant on such a row would be a claim no build makes."

---

## 5. What I would dispute as a geneticist

### 5a. "A PCR library slips more than this" — **overstated and unqualified**, Major

Line 233–234, and this one is about me specifically.

It is directionally right and I am glad it is there. But as a flat claim it is
wrong at one end of my own data and it conflates two very different things:

- **Cycle count.** My reads are a **PCR+ WGS library** — a handful of
  amplification cycles during prep. Targeted STR amplification for genotyping
  runs 25–35 cycles, and stutter compounds per cycle. The two differ by close to
  an order of magnitude. Someone reading "PCR library" will map it to whichever
  they know. **I am the mild case, and the sentence does not let me tell.**
- **Period and length, again.** For a short tetranucleotide tract, 5 in 100
  contraction is roughly what a PCR+ WGS library gives — the shipped number is
  about right there, and "slips more than this" is simply not true. For a
  20-base mononucleotide run it is wrong by several-fold. **A single unqualified
  comparison across a range that spans an order of magnitude is not a claim I can
  use.**

*What would have told me:* "A PCR-amplified library slips more than this, and by
how much depends on the tract: a short tetranucleotide is close to these numbers,
a long mononucleotide run several times worse. Targeted STR amplification is
worse again than a PCR+ whole-genome prep."

### 5b. "real slippage is usually biased towards the shorter tract" — **right**, keep as is

No dispute. Contraction dominates expansion at essentially every period and
length, both in polymerase slippage during amplification and in the enzymology of
replication slippage generally. "Usually" is the correct hedge — there are
expansion-biased loci, but they are the exception.

### 5c. "Those are HipSTR's shipped constants" — **right, but missing the qualification that matters**, Minor

I checked `StutterModel::hipstr_shipped` in `src/ng/alignment/stutter.rs`: 0.05
/ 0.05 whole-repeat, 0.01 / 0.01 part-repeat, one-step share 0.95 both regimes.
The numbers are quoted accurately, and the doc comment there independently states
that HipSTR's *fitted* values are contraction-biased — so line 232's symmetry
complaint is properly sourced.

What is missing: **HipSTR does not normally run on these numbers.** Its default
behaviour is to learn a stutter model per locus by EM; these constants are what
it falls back to when learning is off or impossible. The note says "not a fit",
which is true of the numbers, but a reader can still come away thinking "well,
HipSTR uses these, and HipSTR is a well-regarded STR caller". **They are not what
HipSTR uses.** That makes borrowing them weaker than the sentence implies, and it
is worth one clause: "…HipSTR's shipped fallback constants, which HipSTR itself
normally replaces by fitting per locus."

### 5d. Not disputed but noted: "somebody else's chemistry" (line 236)

Precise and I like it. Worth noting that it is *also* somebody else's **organism**
— HipSTR's constants come from human data, and tomato tract composition is not
human tract composition. "Somebody else's chemistry, on no organism in
particular" would be stronger, and the phrase already exists in
`defaults.rs`'s module documentation.

---

## 6. Reading only this file, can I tell that nothing was fitted?

**For repeat tracts: yes, immediately, and this is the change's clear win.** Line
229, "no stratum was fitted at all", is unambiguous and sits directly above the
table it describes. Line 252 does the same for substitution rates. Both are
new. If I only cared about STRs I would be done in one glance.

**For the run as a whole: no.** — **Major**, and largely unchanged by this
edit. I still have to assemble the answer from six scattered places:

| where | what told me |
|---|---|
| lines 97–110 | three rows saying "not calibrated … no usable error rate could be fitted" |
| lines 139–150 | "nobody said how inbred this sample is" (new wording, and much better) |
| line 168 | `rung = "stated_heterozygosity"` — only meaningful if I recall line 165 calling it "the one that rests on nothing this run measured" |
| line 229 | "no stratum was fitted at all" (new) |
| line 252 | "nothing was fitted for any read group at any stratum" (new) |
| **no `[contamination]` section exists** | **nothing told me. I would have missed this entirely.** |

The contamination one is the worst of these. The preamble promises at lines 36–38:

> "An absent key is not a zero. A missing section, a missing row and a missing key each mean the thing was not measured … **The sections below say which is which where it matters.**"

But there is no `[contamination]` section, so no section says it. **The file
cannot tell me about a section it does not contain.** I would finish this file
without ever learning that `ng` can correct for cross-sample contamination and
that my run did not. For a 63-accession cohort where index hopping is a real
possibility, that is something I would want to know.

**And the file's first section argues the opposite case.** Line 40 onward:

> `[fitted_from]`
> "What these numbers were fitted from. A run whose reference, samples or read groups do not match these is refused…"
> `[fitted_from.census]`
> "Which store of evidence these numbers were fitted from…"

Twelve populated census digests, a reference digest, a full read-group table.
**This is the first thing I read after the preamble and it told me numbers were
fitted.** It took me until line 229 to be sure they were not. The section is
doing its real job — identity, so a file cannot be reused against the wrong data
— but its heading is a claim about provenance that this file cannot support.

*What would have told me, and what I would most like added:* **one line in the
preamble, before `format_version`.** Something like:

```
# **This run fitted nothing.** Every number below is a compiled-in default, a
# stated constant, or something you declared — no number here was measured from
# your reads. What each one costs is said beside it. Contamination is not
# corrected for at all, which is why there is no [contamination] section.
```

That is the difference between "I can find out in 30 seconds" and "I have to read
245 lines and know what to look for". The per-section notes this change adds are
the right work; they just need a header that says the same thing once.

---

## 7. `RepeatTractFitsUsed` — what a future output stage will print about my run

Read from `src/ng/calling/run_report.rs`, lines 249–310.

### Would its explanation help me? Yes, more than the file's does.

Two things it does that the file does not:

1. It states **why a run-level line is needed at all** — "Nothing in the
   parameters file distinguishes *no read group put a read in that stratum* from
   *no stratum was ever fitted, so every tract falls back*, because both are the
   same empty table. This is what says which." That is exactly the confusion I
   had, named exactly.

2. **"so on a PCR library, which slips more than a PCR-free one and slips short
   more often than long, it is wrong in both magnitude and shape."** This is more
   explicit than the file's "A PCR library slips more than this", and *shape* is
   the word that matters. **If one of these two texts is going to say this, this
   is the better sentence — consider moving its clarity into the file.**

3. `read_groups_with_no_slippage_group` — "**Empty on a run that fitted
   nothing**, which declares every read group into one group and simply has no
   strata — being told nothing about slippage and being unable to look it up are
   different failures." I would have read an empty list there as *all good*.
   Being told directly that it is empty for a boring reason on a no-fit run
   stopped me from mis-reading it. Keep.

### What it leaves out that I would want printed beside my calls

- **A count I can turn into a fraction of my calls.** `strata_with_slippage: 0`
  and `fitted_substitution_rates: 0` are counts of *strata* and *cells*. Zero
  answers my question. **Any non-zero number does not.** If a later run prints
  "40 strata carry slippage", I still cannot tell whether that covered 90% of my
  STR calls or 5% of them. What I want is a run-level line in my units: *of your
  N repeat-tract calls, M were scored under the shipped model.* The doc says the
  per-locus counts hold that — but if the output stage prints only this type, the
  number never reaches me.
- **Which strata were fitted.** A run that fitted tetranucleotides and not
  mononucleotides is a completely different situation from the reverse, and I
  need to know which of my calls to throw away. A bare count cannot say. A short
  list — "fitted: period 4 at 5–12 repeats; not fitted: all period 1 and 2" —
  would let me filter.
- **The numbers themselves.** The 5/5/1/1 shares live in the doc comment. The
  struct carries only counts. Please make sure the output actually prints what
  the fallback model *was*, not just that it was used — otherwise the run log and
  the parameters file disagree on how much they tell me.
- **Something that survives into the VCF.** A line in a run log does not travel
  with the calls. When I hand this cohort to a collaborator, they get a VCF. What
  I actually want is a per-record `INFO` flag or `FILTER` on repeat-tract records
  whose slippage was defaulted, so the caveat cannot be separated from the calls.
  I understand that is a different stage; it is the thing I would ask for next.

### Anything that reads wrong to me

- **"one read in twenty comes back a whole repeat short"** (line 273). Same
  ambiguity as the file — 0.05 is all whole-repeat contractions, not
  specifically one repeat. Also, **the file's "5 in 100" is the better
  rendering**: four numbers in twentieths and hundredths ("one in twenty … one in
  a hundred") are harder to add up than four numbers all per 100. Make them
  match, and make them per 100.
- **"it is wrong in both magnitude and shape"** — *shape* I accept without
  reservation. *Magnitude* is asserted without a size and without saying at what
  period, and as argued in 5a it is not true for short tetranucleotides in a PCR+
  WGS prep. **Overstated.** Either give it a range or scope it to where it
  holds.
- **"Zero means every tract's cells take the stated 0.001"** on
  `fitted_substitution_rates` — inherits the same collision with the parameters
  file preamble's "nothing here defaults that number in this file" that I
  describe in section 4. Fix the preamble and this reads fine.

---

## Findings, classified

| # | Finding | Class |
|---|---|---|
| 1d | Note does not say one pair of numbers covers every period and tract length; I would keep mononucleotide calls I should discard | **Blocker** |
| 1a | Shipped model quoted in a different vocabulary from the table's own three numbers; cannot compare to the table or to a later fitted run | Major |
| 1c | "part repeat" defined nowhere in the file; two of the four numbers uninterpretable | Major |
| 2 | "fit the run" names an action with no flag, no input requirement, and no word on whether a 3× cohort can do it | Major |
| 4b | Preamble "nothing here defaults that number" contradicts note "the caller's stated 0.001" | Major |
| 5a | "A PCR library slips more than this" — unqualified across period, tract length, and cycle count | Major |
| 6 | Cannot tell from the file that the *run* fitted nothing; no `[contamination]` section and nothing says one is missing; `[fitted_from]` reads as evidence of fitting | Major |
| 1b | "a whole repeat short" reads as exactly −1; the number is all contractions (4.75 vs 5 in 100) | Minor |
| 2b | "nothing in this file to edit" vs `[stated_constants]` "visible and editable" 29 lines below | Minor |
| 3a | Missing-row paragraph at 177 does not hint an empty table is a third case; correction arrives 50 lines later | Minor |
| 3b | Two of four empty tables in the section carry a note, two do not | Minor |
| 4a | Both new notes open with the identical 13-word clause; the second reads as a repeat and the stakes are not equal | Minor |
| 5c | "HipSTR's shipped constants" without saying HipSTR normally replaces them by fitting | Minor |
| 7 | `RepeatTractFitsUsed` counts strata and cells, not loci or calls; non-zero tells me nothing about coverage of my calls | Major (for the future output stage) |

---

## What I will actually do with this cohort

1. Report SNP and indel calls with the caveat that base qualities were not
   recalibrated and inbreeding was scored at zero, which for tomato accessions —
   many of which self — over-calls heterozygotes. The file told me both, clearly.
2. **Discard every mononucleotide and dinucleotide repeat-tract call.** The file
   did not tell me to; I concluded it from knowing that one flat symmetric pair
   cannot cover that range.
3. Caveat the remaining repeat-tract calls as scored under HipSTR's fallback
   constants, and not report STR allele frequencies.
4. Ask how to fit, because the file does not say.
