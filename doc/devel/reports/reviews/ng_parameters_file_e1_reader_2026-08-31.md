# Reading the parameters file as a user — review of the `defaulted`-multiplier step

Reviewer's position: geneticist who runs callers. I read
`src/ng/calling/parameters_file/testdata/every_shape_as_written.toml` top to
bottom with nothing else open, then ran the caller's refusals for real and read
the printed strings. Base commit `8877316f` plus `tmp/e1_step.patch`, in
`/Users/jose/devel/pop_var_caller-e1-rev3`.

**How the messages were obtained.** A throwaway `#[test]` in
`src/ng/calling/parameters_file/defaults.rs` that edits
`a_file_using_every_shape()` and prints `validate()`'s error. Every message
quoted below is copied from that run's stdout, not from the source. The probe
has been removed again; the worktree holds only the patch.

**One correction to the brief.** The command form given in the task —
`/Users/jose/devel/pop_var_caller/scripts/dev.sh bash -c '… cd $W …'` — does
**not** test the worktree. `scripts/dev.sh` bind-mounts only its own
`PROJECT_DIR`, so `/Users/jose/devel/pop_var_caller-e1-rev3` does not exist
inside the container, the `cd` fails silently, and cargo compiles the main repo
(the log line reads `Compiling pop_var_caller v0.1.0
(/Users/jose/devel/pop_var_caller)` and the run reports `0 tests`). The form
that works is the **worktree's own** copy of the script:

    /Users/jose/devel/pop_var_caller-e1-rev3/scripts/dev.sh bash -c \
      'export PATH=/usr/local/cargo/bin:$PATH
       cd /Users/jose/devel/pop_var_caller-e1-rev3 || exit 9
       cargo test --lib ng::calling::parameters_file::defaults -- --nocapture'

which prints `Compiling pop_var_caller v0.1.0
(/Users/jose/devel/pop_var_caller-e1-rev3)`.

---

## What read well

These are not padding; they are the parts I could act on without re-reading.

- **`[base_quality_calibration]`'s section comment.** "Above one says the
  instrument was optimistic and the reads are worse than it claimed; below one
  says they are better; one leaves the qualities exactly as reported. It is not
  a multiplier on the Phred score, which moves the other way." That last
  sentence is the single most useful line in the file. Multiplying an *error
  probability* and multiplying a *Q score* go opposite ways, and I would have
  got it backwards. Saying so pre-emptively is exactly right.

- **The `[contamination]` three-state comment.** "this whole section absent ->
  nobody identified any contamination / a row with no `measurement` -> this lane
  could not be measured / a zero share with non-zero counts -> measured, and
  found clean", plus "one row a lane, because two lanes of one library can
  differ: index hopping happens on a flowcell, not in a tube." I knew what each
  state meant and why the grain is a lane. I checked and the code keeps the
  three apart (`was_measured()` in `to_run_parameters.rs`), so the comment is
  true.

- **`[sequencing_batches]`'s explanation of why the grain changes.** "The grain
  is a lane everywhere else in this file precisely because two lanes of one
  plant can differ; here they cannot, because a contaminating read is drawn
  against one set of neighbours and a sample has one genotype to draw." I had
  the question ("why is this table by sample when everything else is by read
  group?") and the file answered it in the same paragraph.

- **The `shares_by_repeat_offset` warning**: "An array read as starting at zero
  shifts every length this stratum expects by one repeat." That is the mistake I
  would have made.

- **The new refusal for an edited `defaulted` multiplier** (message C below) is,
  taken on its own, the best of the three `defaulted` refusals: it says what the
  warrant claims, what number that claim implies, and what to type instead.

---

## The five questions

### Q1. Turning one library's calibration off

**First, a wording problem in the question itself, which is also a problem in
the file.** "Stop trusting a library's base qualities" and "turn its calibration
off" are *opposite* actions here. Turning the calibration off means setting the
multiplier to 1.0, which charges every read exactly the error the instrument
claimed — that is trusting the instrument **completely**. To stop trusting a
library you would raise the multiplier above 1, not turn calibration off. The
file's own row comment feeds this confusion:

> ```
> # no calibration: this read group's reported qualities are used exactly as
> # they came, because no usable error rate could be fitted for it
> { read_group = 1, error_probability_multiplier = { value = 1.0, warrant = "defaulted" } },
> ```

I read "no calibration" as *this lane is untrusted / excluded*. It means the
opposite. **Minor** — the section header two lines up does say "one leaves the
qualities exactly as reported", so a careful reader recovers. What would have
told me: `# not calibrated: this read group's reported qualities are taken at
face value, …`.

**Second, and this is the substantive finding.** The edit I am most likely to
make is **not** the one the new check catches.

Read group 0 is the fitted row:

```
{ read_group = 0, error_probability_multiplier = { value = 1.0324, warrant = "fitted_here", observations = { reads = 812344 } } },
```

- If I do it the way the file's header teaches — set `value = 1.0`, set
  `warrant = "supplied"`, delete `observations` — **the file is accepted**.
  Correct.
- If I change only the number and leave the rest of the line alone —
  `value = 1.0`, still `warrant = "fitted_here"`, still
  `observations = { reads = 812344 }` — **the file is also accepted, silently.**
  The run then reports my typed 1.0 as fitted from 812,344 reads.

So the *only* warrant slip this key polices is the `defaulted` one, which is the
rarer of the two: read group 0 is where a person actually wants to intervene,
and read group 1 is already at 1.0 and has nothing worth changing.

**Major.** I do not think this is a bug in the patch — the code comment says so
outright, and it is right that a fitted 1.0 is legitimate. The problem is what
the reader learns. Having been refused once for leaving a `defaulted` warrant
behind, I conclude *the caller checks my warrants*, and then make the same class
of mistake at a fitted row with no warning at all. What would have told me: one
sentence in the file's header saying which slips it can catch and which it
cannot — e.g. "*A warrant left behind on a number you changed is refused only
where the file's own numbers can refute it: a `defaulted` value that is not the
caller's constant, and a `fitted_here` fallback concentration with no fitted
stratum. A changed number still labelled `fitted_here` or `borrowed` cannot be
told from a real one, and nothing will stop you.*"

**Now the case the check does catch.** Change read group 1's number, leave
`warrant = "defaulted"`:

> `the parameters file cannot be used:
> base_quality_calibration.by_read_group[read_group = 1].error_probability_multiplier
> is 1.25, and its warrant is `defaulted`, which says no rate was fitted for this
> read group and its qualities were left as the instrument reported them — which
> is a multiplier of 1.0; a number you changed is one the run was handed, so
> change the warrant beside it to `supplied``

**Yes, I can act on this**, and quickly: the key path is searchable, "1.25" is
in my file, and "change the warrant beside it to `supplied`" is a literal
instruction. Two small snags:

- It does not tell me to delete `observations`, which the file's header does
  ("change its warrant to \"supplied\" **and delete its `observations`**"). Here
  that is harmless because a `defaulted` value cannot carry `observations` in
  the first place, but I did not know that while reading, so I went back to the
  header to check. **Minor.**
- "the warrant **beside** it" is right for this file (both are on one line, in
  an inline table), but the outlier-weight refusal uses the same word for a key
  that is also inline, and the file's header says "change **its** warrant".
  Consistent enough. No action.

### Q2. Is the quoted number spelled the way my file spells it?

**For the three `defaulted` refusals, yes — and the patch is what makes that
true for one of them.** `{:?}` on an `f64` prints `1.0`, `{}` prints `1`. The
patch changes `the_repeat_tract_numbers_are_what_they_claim` from `{value}` to
`{value:?}` and from `{STATED_FLAT_CONCENTRATION}` to
`{STATED_FLAT_CONCENTRATION:?}`, so that message now says:

> `repeat_tracts.fallback_length_spectrum_concentration is `defaulted` at 3.5,
> and the constant a run falls back to is 1.0; …`

Before the patch the constant printed as `1`, which is not a string that occurs
anywhere in this file. Good catch, and the new multiplier refusal was written
with `{:?}` from the start (it prints `2.0`, not `2`, when I set the value to
two).

**But the patch fixed one message in each of the two functions it touched and
left the neighbour beside it broken.** Both remaining ones I hit by accident
while probing:

> `base_quality_calibration.by_read_group[read_group = 0].error_probability_multiplier
> **is 0**, and a multiplier on an error probability is above zero — a zero
> charges every read of the library the error floor`

> `repeat_tracts.fallback_length_spectrum_concentration **is 0**, and a
> concentration is above zero`

My file spells that value `0.0`. Searching a TOML file for `0` matches almost
every line in it. This is the *same* defect the patch just fixed, two and
fourteen lines away respectively, in the two functions the patch edits.
**Major** — I could not find the offending line by searching for the number the
message gave me. Any negative value has the same problem: `-1.0` prints as `-1`.

One latent inconsistency, no live symptom: `stated_constants.repeat_tract_outlier_weight`
still uses `{}` for both its own value and `DEFAULT_OUTLIER_WEIGHT`. It happens
never to bite, because a valid outlier weight is strictly inside (0, 1) and so
never prints as a whole number. But now two of the three constants are formatted
one way and the third another. **Minor.**

### Q3. Does the new refusal agree with the file's header?

**No. The header now says something false, and the patch is what made it
false.** Lines 10–16 of the file I hold:

> ```
> # **Two keys do not take every warrant.**
> # `repeat_tracts.fallback_length_spectrum_concentration` is `fitted_here` only
> # where this file holds a fitted stratum spectrum for it to be the median of,
> # and `defaulted` only at the built-in constant;
> # `stated_constants.repeat_tract_outlier_weight` is `defaulted` only at the
> # caller's own constant. Both take `supplied` freely, which is what you write
> # when you change one. Anything else is refused, and says so.
> ```

I read this as an exhaustive list: **two** keys are fussy about warrants,
everything else takes any of the four. The patch adds a third —
`base_quality_calibration.by_read_group[…].error_probability_multiplier` is now
`defaulted` only at 1.0 — and does not touch the header. So the file tells me
the calibration multiplier takes every warrant, and the caller refuses it.

**Blocker.** Not because the check is wrong — it is right — but because the one
paragraph in the file whose whole job is to enumerate the exceptions omits the
exception this change creates, and it counts them out loud ("Two keys"), so a
reader cannot even suspect the list is stale. I would have edited the multiplier
under a `defaulted` warrant believing the file had told me that was allowed.

What would have told me: the header saying "**Three keys do not take every
warrant**" and a clause for the multiplier — "…
`base_quality_calibration.by_read_group[…].error_probability_multiplier` is
`defaulted` only at a multiplier of 1.0, which is what leaving the qualities
alone means." The golden file `every_shape_as_written.toml` and whatever writes
that header both need the edit.

### Q4. Do the three `defaulted` numbers refuse the same way?

**No.** Here they are side by side, verbatim, with the leading `the parameters
file cannot be used: ` stripped:

| key | message |
|---|---|
| multiplier (new) | `… is 1.25, and its warrant is \`defaulted\`, which says no rate was fitted for this read group and its qualities were left as the instrument reported them — which is a multiplier of 1.0; a number you changed is one the run was handed, so change the warrant beside it to \`supplied\`` |
| outlier weight | `… is 0.05, and its warrant is \`defaulted\`, which says this run inherited the compiled-in 0.01; a number you changed is one the run was handed, so change the warrant beside it to \`supplied\`` |
| fallback concentration | `… is \`defaulted\` at 3.5, and the constant a run falls back to is 1.0; a number somebody chose is \`supplied\`` |

The first two are one message with one substituted clause. **The third is a
different message in three separate ways:**

1. **Different sentence shape.** "is *X*, and its warrant is `defaulted`" versus
   "is `defaulted` at *X*". Same facts, transposed. I had to re-read the third
   to be sure it was telling me the same thing.
2. **Different person, and it stops being an instruction.** "a number **you
   changed** is one the run was handed, **so change the warrant beside it to
   `supplied`**" versus "a number **somebody chose** is `supplied`". The first
   two tell me what to type. The third states a general fact about the format
   and leaves me to infer the action. "Somebody" is also strange in a message
   addressed to the person who edited the file — the somebody is me.
3. **Different name for the same idea.** "the compiled-in 0.01" (outlier
   weight), "a multiplier of 1.0" (new one), "the constant a run falls back to"
   (concentration). Three vocabularies for *the number the binary holds*.

**Major** — I could act on all three eventually, but only after re-reading the
third, and only because I had the other two beside it. What would have told me:
one shape for all three, the outlier weight's, which is the one that both names
the claim and names the fix.

**This also contradicts a claim the patch itself makes.** The doc comment on
`a_defaulted_value_that_is_not_the_binarys_own_number_is_refused` in
`defaults.rs` says:

> "All three refusals say the same thing, so the reader meets one shape of
> message rather than three."

They do not say the same thing and there are two shapes, not one. The test that
sentence sits on cannot notice, because it only asserts that each message
contains the key, the constant and the word `supplied` — three substrings that
two different sentences can both satisfy. **Blocker for the claim**, since the
next person to touch this will trust it rather than print the messages.

**A fourth message worth knowing about**, which I hit while probing and which is
the best-written of the lot: leaving `observations` beside a `defaulted`
warrant gives

> `… .error_probability_multiplier.observations is written beside a \`defaulted\`
> warrant, and a stated constant has nothing behind it; delete the
> \`observations\` table, or — if you changed the number — set the warrant to
> \`supplied\`, which keeps its count`

It names both branches and the consequence of each. If the three above were
rewritten to this standard the section would be finished.

### Q5. Is there a `defaulted` number I can change with no objection?

**Yes — two, and both surprised me.**

```
{ sample = "TS-1", inbreeding_coefficient = { value = 0.9, warrant = "defaulted" } }
```
→ **accepted.**

```
{ read_group = 0, period = 2, reference_repeats = 6, ploidy = 2, rate = { value = 0.4, warrant = "defaulted" } }
```
→ **accepted.**

I can mark a sample's inbreeding coefficient `defaulted` at 0.9, or a
repeat-tract substitution rate `defaulted` at 0.4 — forty errors in a hundred
bases — and the run takes the file and reports both as numbers it inherited from
its own binary. There is no compiled-in inbreeding coefficient and no
compiled-in substitution rate, so `defaulted` on those keys is a claim about the
caller that is simply not true of any build.

**Is it a surprise?** Yes, in the specific sense that matters. Having been
refused on three keys for exactly this — a `defaulted` warrant over a number the
binary does not hold — I formed the rule "`defaulted` means the caller's own
constant, and the caller checks it." Two keys break that rule silently, and they
are the two where a wrong number does the most damage per line: an inbreeding
coefficient of 0.9 changes every genotype call for that sample, and a
substitution rate of 0.4 changes every repeat-tract call for that read group.

**Blocker.** Whichever way it is resolved — refuse `defaulted` on those two keys
because nothing defaults them, or say in the file's header that `defaulted` is
unchecked outside the three named keys — the current state teaches a rule and
then breaks it without saying so. What would have told me: the header paragraph
from Q3, extended: "*On every other key `defaulted` is not checked against
anything, because there is no constant for it to be; a `defaulted` inbreeding
coefficient or substitution rate is a claim nothing can refute, so do not write
one.*"

---

## `defaults.rs`'s module header

Read as the module says it should be read — "which of the caller's numbers are
guesses and how badly". Everything below is the prose at the top of the file,
not the tests.

### Blockers

**B1 — "A multiplier of one asserts nothing about the chemistry."**

> "**A multiplier of one asserts nothing about the chemistry.** … A run that took
> it is not guessing at a quantity; it is declining to change one."
> (also in `likelihood/mod.rs`: "this one asserts nothing about the chemistry at
> all")

This reads as wrong to me. Charging every read exactly the error probability its
Q-score claims is not the absence of an assertion about the chemistry; it is the
assertion that **the instrument's quality scores are correct**, which is
routinely false and is the entire reason base-quality recalibration exists as a
step. The second clause — "declining to change one" — is the accurate statement,
and it is *not* the same statement.

Why this is a Blocker and not a wording nit: the module's stated job is to tell
me which defaults are risky. It sorts the multiplier into a category ("asserts
nothing") that no other default is in, next to two it calls "inherited guesses",
and the effect is to tell me an uncalibrated read group carries no risk. It
carries a different risk, not none — and unlike the other two it is *per read
group*, so a cohort can be half-calibrated and half not, which is the case I
would most want flagged. What would have told me: "*A multiplier of one changes
no read's error probability, so a run that took it is scoring that library on
the instrument's own Q-scores. That is an assumption, not the absence of one —
it is just not an assumption this module can put a number on.*"

**B2 — the section that says the fifth number has no default, then gives its
default.**

> "# The fifth number, which has no default and is not one of these
>
> **The per-(stratum × slippage group) slippage numbers are to be fitted from the
> GIAB HG002 alignments and compiled in like the rest, and that measurement does
> not exist** … Until it does, a run with no slippage fit writes no slippage
> rows, and what a repeat tract is then scored under is decided one level down,
> at the tract: `inference::repeat_tract_parameters` gives such a cell
> `StutterModel::hipstr_shipped` and `Provenance::Defaulted` …"

The heading says the slippage numbers have no default. Two sentences later they
have one: HipSTR's shipped stutter model, marked `Defaulted`. I read the heading,
believed a run without a slippage fit would refuse or emit nothing, read on, and
found it quietly scores my tracts under somebody else's constants.

And the thing I most needed is not there. I call **tomato**. HipSTR's shipped
model was fitted on human data. Being told my repeat tracts fall back to it is
precisely the "how badly is this a guess" the module opened by promising, and
the module does not say what those numbers are, what organism they came from, or
which direction they are likely wrong in for a plant. What would have told me:
"*Its default is not a number in this module but a whole model: a tract with no
slippage row is scored under `StutterModel::hipstr_shipped`, HipSTR's own
constants, fitted on human sequencing. On any other organism that is a borrowed
guess of unknown size, and `inference::repeat_tract_parameters` counts how many
cells took it — check that count.*"

### Majors

**M3 — every citation is a section number with no document.**

The header cites, in order: "spec §8", "spec §8, owner's decision of
2026-08-28", "spec §7", "spec §8", "spec §3.8", "spec §5's first row", "spec
§8's third bullet, §12 question 1", "spec §5's third row", "spec §2.1". It never
once says which file that is. (`likelihood/mod.rs` names it —
`doc/devel/ng/spec/parameters_file.md` — but I only found that because I was
reading the patch.) Every one of these is a place where the header says "the
reason is elsewhere" and then does not say where. **Name the document once, at
the top.**

**M4 — the header promises a size and then says nobody has one.**

> "…which of them are honest"
> "Both would be different numbers on different data and nobody knows by how
> much"

The opening sentence sets up a grading; the body's answer for the two guesses is
"nobody knows by how much". That may be the honest state of affairs, but then
the module is not doing what its first line says. What would have told me — even
without a measurement — is the *direction and the observable*: which way each
number pushes calls if it is too high, and what in the output I could look at to
notice. For the outlier weight: does a too-low 0.01 make repeat-tract calls
over-confident or under-confident, and would I see it as excess homozygosity, as
inflated GQ, as something else? That is actionable; "nobody knows" is not.

**M5 — "the *simple* case of that model rather than the weak one".**

> "A run told nothing about contamination is uncontaminated, and the read
> likelihood then computes its plain formula — the *simple* case of that model
> rather than the weak one."

I cannot act on this sentence. "The weak one" is never defined, here or
anywhere in the file. Weak how — a weaker prior, a weaker likelihood, a
degenerate case? The whole point of the sentence is a contrast, and one side of
it is a word I have to guess at.

The same sentence's first clause is a second problem: "A run told nothing about
contamination **is uncontaminated**." The written parameters file is careful
here — "this whole section absent -> **nobody identified any** contamination" —
and the module states as a fact about the data what the file states as ignorance.
For a geneticist those are different: I will report a cohort differently if I
believe it was checked and clean than if I believe nobody looked. What would have
told me: "*A run told nothing about contamination is **scored as** uncontaminated
— the model has no fraction to mix in, so it computes its plain formula. That is
a modelling default, not a finding about your samples.*"

**M6 — "the tract ladder", used before it is defined.**

The table row reads "the tract ladder's fallback concentration", and the phrase
"the ladder's bottom rung" appears twice more. Neither "ladder" nor "rung"
occurs anywhere in the parameters file I was handed; the file calls this key
`fallback_length_spectrum_concentration` and describes it as "what a tract falls
back to where neither its own stratum nor its period was fitted". I could
eventually match them up, but only by matching the constant name. A term that
does load-bearing work in a table row needs one clause defining it.

### Minors

**m7 — "the same double".**

> "**A defaulted value and a fitted one are the same double.**"

"Double" is a programming word for a 64-bit float. In a paragraph aimed at
somebody deciding whether to trust a number, "the same number" says it without
the detour.

**m8 — the bookkeeping does not add up on a first read.**

The heading says "# The four, and the three different kinds of thing they are";
the table lists four rows, one of which is not a number ("**absence** — no
`[contamination]` section"); the `validate.rs` comment added by the same patch
calls the multiplier "the third of the three constants"; and then a whole section
introduces "The fifth number". Four, three, three, fifth. Each count is
defensible on its own terms and together they cost me a re-read to reconcile.

**m9 — "the repeat-tract outlier weight, one a run".**

The table's first column. I parsed "one a run" as a typo before working out it
means "one value per run" (as against the multiplier's "per read group" on the
line above, which is spelled out). Spell it out here too.

**m10 — "production's 0.01".**

> "The outlier weight is production's 0.01, which nothing in this project has
> measured"

"Production" is undefined for a reader in my position — production of what?
Interestingly the parameters file itself says it better, in one line: "inherited
from the existing caller and never measured here". Use the file's own phrasing.

**m11 — the heading quotes a spelling the file does not use.**

> `# What "marked \`Defaulted\` when used" is worth`

`Defaulted` with a capital D is the Rust variant name. The file I edit spells it
`defaulted`. The module elsewhere is careful about this — there is a whole
function, `the_word_for`, whose doc comment exists to stop exactly this — so the
heading is out of step with its own file's rule.

**m12 — "the library" for a key indexed by read group.**

Not new in this patch, but it sits in the function the patch edits and I met it
while probing: the zero-multiplier refusal says "a zero charges every read of
**the library** the error floor", while the key it names is
`by_read_group[read_group = 0]` and the new refusal beside it says "no rate was
fitted for this **read group**". A library and a read group are different things
in this file — `[contamination]` spends a paragraph on the fact that two lanes
of one library can differ — so using them interchangeably in two adjacent
messages about the same key is a small trap. The zero message is also strictly
wrong: it charges every read of that *read group*, not of the whole library.

---

## Summary of severities

| # | Finding | Severity |
|---|---|---|
| Q3 | File header says "**Two** keys do not take every warrant"; the patch makes it three and does not update the header | Blocker |
| Q4 | `defaults.rs` claims all three refusals "say the same thing… one shape of message"; there are two shapes, and the test cannot notice | Blocker |
| Q5 | `defaulted` is accepted at any value on `inbreeding_coefficient` and on repeat-tract `rate`, after being refused on three other keys | Blocker |
| B1 | "A multiplier of one asserts nothing about the chemistry" — it asserts the instrument's Q-scores are right | Blocker |
| B2 | "The fifth number, which has no default" then gives its default (`hipstr_shipped`), and never says it is human-fitted | Blocker |
| Q1 | Only the `defaulted` warrant slip is caught; the likelier slip — editing a `fitted_here` number and leaving warrant and `observations` — is accepted silently, and nothing says so | Major |
| Q2 | `is 0,` where the file says `0.0` — in both functions the patch touched, beside the message it fixed | Major |
| Q4 | Third `defaulted` refusal has a different shape, a different person, and stops being an instruction | Major |
| M3 | Nine "spec §N" citations, no document named | Major |
| M4 | Promises to say how badly each default is a guess; answers "nobody knows by how much", with no direction or observable | Major |
| M5 | "the *simple* case… rather than the weak one" — "the weak one" undefined; and "is uncontaminated" overstates "nobody looked" | Major |
| M6 | "the tract ladder" / "rung" used before defined, and absent from the file it describes | Major |
| Q1 | "no calibration" row comment reads as *untrusted* when it means *taken at face value* | Minor |
| Q2 | Outlier weight still formats with `{}` where the other two now use `{:?}` — latent, no live symptom | Minor |
| Q1 | New refusal omits "delete its `observations`", which the file's header teaches | Minor |
| m7 | "the same double" | Minor |
| m8 | four / three / three / fifth | Minor |
| m9 | "one a run" | Minor |
| m10 | "production's 0.01" — undefined; the file says it better | Minor |
| m11 | heading quotes `Defaulted`, the file spells `defaulted` | Minor |
| m12 | "the library" for a key indexed by read group, in adjacent messages | Minor |
