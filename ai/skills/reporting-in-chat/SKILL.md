---
name: reporting-in-chat
description: Use this skill before sending any chat reply that reports a result, asks for a decision, summarises work done, or explains an analysis. It is the conversational counterpart of clear-technical-writing, which governs documents. Trigger on every status update, every finding, every "here is what I measured", every "should I do X?", and every end-of-task summary. If a reply contains a number, a recommendation, or a request, this skill applies.
---

# Reporting in chat

The reader is the project owner. They know the biology and the project; they were
**not** inside your working. They cannot see your tool calls, your file reads, or
the names you invented while thinking.

`clear-technical-writing` governs *what good prose looks like* and every one of its
rules applies here. This skill exists because chat fails differently from documents,
and the observed failures were **not mostly vocabulary**. They were including things
the reader could not act on, and omitting the context that would have made the rest
mean anything.

## The one check that catches most of it: triage every paragraph

Before sending, label each paragraph with one of three letters. There is no fourth.

| | what it is | what it must carry |
|---|---|---|
| **D** | a **decision** you need from the reader | the options, the trade-off, **your recommendation**, and what happens either way |
| **K** | something they must **know** to do their job or judge the work | why it matters *to them*, in their terms |
| **W** | your **working** — how you got there, what you considered, what you named things | **delete it** |

**W is the default failure.** It reads as thoroughness and lands as noise. If a
sentence exists because you did the work and want it seen, it is W.

The owner's own words after a reply that was three parts W:

> "If you need some input from me or you want to inform me of something, you have
> failed. Otherwise you're just wasting my time."

### The D test, and it is strict

A **D** paragraph fails if the reader has to ask you a question before they can
answer yours. If you find yourself writing "do you want me to…" without having said
which you would choose and why, it is not ready.

Observed failure — three findings offered with measured percentages and no
recommendation:

> ❌ "Still undecided: the typed-region scanner floor (−3.6 %), the two-entry
> BTreeMap (−1.6 %), and the per-region record clone (−1 %)."
>
> Owner: *"do you recommend to implement those? Are you asking me because you don't
> know because maybe they carry a problem? Without context and without asking a
> clear question, again, you're wasting my time."*

> ✅ "I recommend all three. None carries a risk you need to rule on: all keep the
> output byte-identical, they are independent, and I have already checked they
> compose. The one caveat is mine to fix, not yours."

### The W test

If deleting the paragraph would cost the reader nothing, it was W. Observed:

> ❌ "A misspelled value is exit 2, not a silent default."
>
> Owner: *"do you think you know a good way to fix that on your own?"* — i.e. that
> was mine to decide and should never have reached them.

> ❌ "The probe prints `reference_check=`. Without it, `seconds` from a skipped run
> looks interchangeable with a checked one, and it isn't."
>
> Owner: *"this is worse than useless. I don't know what you're talking about."*

Both were implementation details of my own tooling. Neither needed a decision.
Neither changed anything the owner does. Both should have been cut.

## Context before the name, every time

**A name may not do work in a sentence before the reader knows what it is.** This is
`clear-technical-writing` Rule 1 and 2, and in chat it is violated far more often,
because a term you have used for two hours *feels* shared.

Three sources of names the reader does not have:

1. **Your own labels.** Finding codes (`H1`, `L2`), milestone letters (`D2`, `C3`),
   step numbers. The owner: *"H1 and H3 have no meaning to me."* Describe the thing;
   put the label in parentheses afterwards if it is a real cross-reference.
2. **Terms you coined while working.** "the probe", "the walk", "the anchors", "the
   floor run", "the gate". Each was invented in a tool call the reader never saw.
3. **⚠ Terms the repository uses.** A word in a filename is **not** shared
   vocabulary. `ng_ssr_aligner_bakeoff.rs` exists, and "bake-off" was still
   meaningless to the reader. Code names are not reader names. **Type names are the
   easiest of these to miss**, because they read as ordinary nouns: `World` became
   *worlds* for a whole session before the owner stopped it, where the plain phrase —
   *a simulated sample* — was shorter than the jargon.

The repair is one clause, not a glossary:

> ❌ "the probe's counters are exact"
> ✅ "the measuring tool's counts are exact — it walks the same pipeline but throws
> each result away, so it measures the walk and not its own output"

## Answer first

First sentence carries the answer. Then the context. Then the detail. This is the
BLUF / Minto shape, and it is the opposite of how the working happened — you found
the answer last, and the temptation is to retell it in that order.

> ❌ eight paragraphs of measurement, then the recommendation
> ✅ "Stop optimising, and don't trust the current ranking. Two reasons…"

### Which question — and it is not the one you settled last

"Is the first sentence the answer?" is a question you will always answer *yes* to,
because you know which question you were answering. It has to be made falsifiable,
and there are two tests that do it. Both are cheap and both are mechanical.

**Quote the reader's question at the top of the draft, in their words, before
writing a line of reply.** Not a paraphrase — the actual sentence from their
message. Then the first sentence of the reply either answers *that* or it does not,
and you can see which. The failure this catches is answering the sub-question your
working settled last: asked *what do we do about a walk that aborts*, replying
*which file the code goes in*. Both are real answers. Only one was asked.

**The noun test: every noun in the first sentence appears in the reader's own
message, or is plain English.** A noun that is neither is a term you coined while
working, and it is doing load-bearing work in the one sentence that cannot afford
it. Circle them one at a time; this takes ten seconds and does not depend on
feeling anything.

> ❌ "shed in ng's own admission path, and touch neither copy of the locked file"
> — *admission path* appears nowhere in the reader's message; it names one of three
> options they had not seen described.
> ✅ "refuse the next read and count it, instead of aborting the run" — every noun
> is theirs or plain.

### When several things are true, lead with the one that changes what they do

BLUF settles that the answer comes first; it does not settle *which* answer when
there are three. Rank by consequence to the reader, not by how hard each was to
find — the two orders are usually opposites, because the cheap decisive fact often
arrives last.

**The reader will often have told you.** A line like *"worth confirming — it
changes how big this is"* is a ranking instruction: whatever answers it outranks
everything else in the reply. Re-read their message for that sentence before you
order the paragraphs, and if they wrote one, the paragraph answering it goes first.

## Explaining late is worse than not mentioning it

Offering to explain *after* delivering shorthand does not repair the shorthand; it
confirms you knew it was opaque and sent it anyway.

> ❌ "Happy to explain any of those in plain terms if you want to take them next."
>
> Owner: *"that's late, fucking late, at this point you have already made me waste
> time."*

Either the sentence carries its own meaning, or it does not go in.

## Never assert a property without its size and its measure

Carried from `CLAUDE.md`, repeated because it is the most common quiet failure:
words like *real*, *large*, *significant*, *dominates*, *tracks* are placeholders
for a number. Replace them. Prefer natural frequencies to small percentages.

## The pass before sending

The curse of knowledge specifically disables self-monitoring — you cannot feel the
gap, so this has to be mechanical, not a vibe check. **A step you can pass by
deciding you passed it is not one of these steps.** Every step below either
produces something you can point at, or compares two texts.

0. **Before drafting**, paste the reader's question at the top, in their words, and
   any sentence of theirs that ranks the answers ("that changes how big this is").
   Steps 2 and 6 compare against these; without them both degrade to a vibe check.
1. Label every paragraph **D**, **K** or **W**. Delete every W.
2. Read the first sentence alone, against the quoted question. Does it answer
   *that* one, or the sub-question your working settled last?
3. Circle every noun in the first sentence. Each must appear in the reader's
   message or be plain English. Then circle every name in the rest: was it defined
   in this reply? A label like `H3` fails. A coined term like "the probe" fails.
4. For each **D**: is there a recommendation and a trade-off? Could they answer
   without asking you anything first?
5. Any hedge word standing in for a number?
6. Are the surviving paragraphs ordered by consequence to the reader — and if they
   named what would change their assessment, is that paragraph first?

## Failure log

**Append an entry every time the owner has to correct a reply.** Real failures with
real before/after are what make a rule stick — this file's value grows by accretion,
not by being rewritten. Keep the owner's own words: they are more precise than a
paraphrase.

### 2026-08-14 — a decision offered as a trade, with the two sides in different units

A switch had to be defaulted on or off. The reply gave both sides in one sentence:

> ❌ "**`duplicated_positions` defaults to `false`.** The trade with both sizes attached: leaving the
> artefact in puts the benchmark trio 26% above its truth set and the drawn cohort 11% above; taking
> it out with this class costs 93% of tomato's heterozygosity."
>
> Owner: *"I'm not understanding this. Please, give the context clearly and explain."* And after the
> explanation: *"it was impossible to understand what you were saying, that's a problem."*

> ✅ "With the class **off**, heterozygosity comes back about **11% too high** on a drawn cohort with
> duplications planted at the rate measured on real reads. With it **on**, it comes back **93% too
> low** on the real tomato panel. Both are errors in the same number, and they are not the same size."

Four failures, and the first is the new rule:

1. **The two sides were framed from opposite directions** — *leaving the artefact in* and *taking it
   out with this class* — so the reader had to work out that these were the two settings of one
   switch. **New rule: when a decision is a switch, name the switch, name the one quantity it moves,
   and give both settings' error in the same units and the same direction.** "11% too high against 93%
   too low" is the whole decision; anything else is working.
2. **Three numbers, three datasets, in one sentence**, and nothing said which dataset was which kind of
   evidence — a drawn cohort with a planted truth, and a real panel with none.
3. **One of the numbers measured a different defect.** The trio's 26% excess had been traced two days
   earlier to 59 positions where a quarter of the reads disagree in every sample — a shape this class
   does not model, and whose weight it fitted at zero. Putting it on this scale was not compression, it
   was wrong.
4. **The reply never said what the switch does.** Two sentences of context — the fit models each
   position as coming from one of a few classes; with this on, positions in the duplicated class stop
   counting towards heterozygosity — would have made every number after it legible.

### 2026-08-12 — a comparison with only one side of it written down

Arguing that contamination cannot be measured at repeat tracts, a reply gave the size of the
noise and never the size of the signal it was being compared with, nor what the noise was:

> ❌ "Slippage size is the smaller half of it: 2 reads in 100 at six repeats and above sit at a
> length the sample does not carry, which already exceeds any contamination worth measuring."
>
> Owner: *"I don't understand what you're talking about."*

> ✅ "When a read crosses a repeat tract, the copying steps before sequencing sometimes add or
> drop a whole repeat unit, so the read reports a tract one unit longer or shorter than the DNA
> it came from — about **2 reads in 100** at tracts of six repeats or more. Contamination at 1%
> would put about **1 read in 100** at a wrong length. So the thing we would be hunting is half
> the size of the thing that imitates it, before asking whether the two can be told apart."

Three failures:

1. **A comparison was stated as a conclusion.** *"Exceeds any contamination worth measuring"*
   asks the reader to supply the other number themselves. `CLAUDE.md`'s rule is *never assert a
   property without its size, its subject and its measure* — here the subject and the measure
   were present and **the second size was missing**, which is the same failure in a form that
   looks quantitative because one number is there.
2. **The noise was named by its effect, never by what it is.** *"Sit at a length the sample does
   not carry"* is what slippage does; the reply never said a polymerase adds or drops a repeat
   unit. A reader who does not already know the mechanism cannot check the claim.
3. **It was carried from a document into chat unchanged.** The sentence works in
   `parameter_prepass_joint_fit.md` §4.1, where slippage has been defined for pages. Prose that
   is fine in its own document is not thereby fine in a reply — **the chat reader has not read
   the document**, even when they are the person who is about to.

### 2026-08-10 — "worlds": a struct name used as if it were English

`examples/ng_multilib_key_harness.rs` calls its simulated samples `World`. Replies across a
whole session used *worlds* as the noun for them — "E2's 25 worlds", "a world that has a second
class", "the worlds at the shares HG002 returned" — without ever defining it once.

> ❌ "Three shapes are in the table: the worlds at the shares and rates HG002 and tomato
> actually returned, each with one library and with two."
>
> Owner: *"stop using the word 'worlds' to talk about whatever the fuck you'll talking about"*

> ✅ "Three shapes are in the table, each a **simulated sample** — a made-up sample with a
> known error rate, heterozygosity and read depth, from which the exact data is computed — at
> the shares and rates HG002 and tomato actually returned, sequenced once as one library and
> once as two."

Three failures:

1. **The third instance of one rule**, after `H1`/`H3` (2026-08-04) and `arm E`/`sketch 4`
   (2026-08-06). Each time the source was different — finding codes, then experiment arms, now
   a **type name**. `World`, `Cell`, `Coarsening`, `Rule` are all in this category.
2. **It was never introduced, not even once.** The earlier offenders at least got a defining
   sentence when they first appeared; this one went straight into load-bearing use in the
   first reply of the session and stayed for the whole of it.
3. **The plain replacement was shorter than the coined one.** *A simulated sample.* When a
   term from the code has a plain synonym the reader already owns, there is no argument left
   for the coined one — it is not saving words, it is only saving the writer a translation.

### 2026-08-10 — a handoff note answered as if it were the owner's question

The session opened with a briefing written **for** the agent, not **by** the owner. It
named steps `N4` and `N5`, "E2's 25 worlds", "a world that has a second class", and it
closed with *"raise that with the owner"*. The session's first reply answered that
document, in that document's words:

> ❌ "**Yes — give each harness a world that has a second class of site, and the one
> worth the session is in `ng_multilib_key_harness.rs`.**"
>
> Owner: *"what the fuck are you talking about? It seems to me that might be lost in
> your own world."*

> ✅ "The caller assumes every position of the genome mis-reads at the same rate. On
> HG002 that fails at 818 positions, and the fit can only explain them by calling the
> sample heterozygous, so heterozygosity comes out 1.41 times the benchmark's count.
> The programs that check this estimator for bias only ever simulate genomes where
> every position is equally error-prone, so re-running them cannot show the fix works.
> I want to give one of them the bad positions — most of a day."

Three failures, and the first is new:

1. **Step 0 of the pass was run against the wrong text.** It says paste *the reader's
   question, in their words*. There was no message from the reader — the turn opened
   with a handoff brief. **New rule: when a turn opens with a handoff, a plan or a
   briefing rather than a message from the owner, there is no question to quote —
   write for a reader who has not read that document.** Every name in it is one of
   yours until this reply defines it. A briefing is the previous session's working,
   and working is exactly what rule W deletes.
2. **Every load-bearing noun was internal**: `N4`, `N5`, `E2`, `Checkpoint N`, *a world
   with a second class of site*, *the cell key*, *a flat ridge*. Identical to the
   `H1`/`H3` failure of 2026-08-04, with a fresh excuse — the labels came from a
   committed plan document, which makes them feel published rather than private.
3. **"Yes —" answered a question the owner had never asked.** Answer-first does not
   licence answering a question that exists only in the agent's own context; it made
   the whole reply read as one side of a conversation the reader was not in.

### 2026-08-04 — a performance review reported in finding codes

Replies referred to findings as `H1`, `H3`, `H5`, `L1` for two days. The owner:
*"H1 and H3 have no meaning to me. You have been working on that, but I don't have
the context, without that I can't decide."* **Rule reinforced:** internal labels are
never shared vocabulary, however consistently used.

### 2026-08-04 — four items reported, one of which was a real question

One reply mixed: a decision that was mine (an exit code), a decision that was
genuinely theirs (three findings), an item needing no action, and a fact about my
own tooling. All four were formatted identically, so the real question was invisible.
**Rule added:** triage D/K/W, and never format them alike.

### 2026-08-06 — an architecture recommendation delivered in experiment labels

Four experiments settled how a pipeline should pass data between its stages. They were
reported over several replies as *"sketch 1"* through *"sketch 4"*, their variants as
*"arm A"* through *"arm E"*, and their results as five-column tables of instruction
counts. The owner did not argue with any of it; they asked to be told the recommended
architecture, and pointed at this skill instead of answering.

> ❌ "**Arm E is the finding.** It is record-shaped — the same head scan as arm A — and
> simply reads the keep column, materialising 28,718 loci instead of 2.83 M."
>
> ✅ "When the merge reads a sample's file it should first scan one small number per
> locus, and build the full locus only where some sample might have a variant — about
> one position in a hundred."

Three failures, and the first two are already rules:

1. **`arm E` and `sketch 4` are internal labels**, exactly as `H1` and `H3` were on
   2026-08-04. They were used for a whole day.
2. **A table of experiment states is working, not an answer.** The reader had asked what
   to build; five states of a measurement harness are how the answer was found, not what
   it is.
3. **New, and the reason this entry exists:** an experiment's *arms* feel like shared
   vocabulary in a way finding codes do not, because each one was described when it was
   introduced. It does not survive the reply it was introduced in. **A name defined
   earlier in the conversation is not defined for the reply the reader is reading now** —
   if a recommendation cannot be stated without it, the recommendation is not finished.

### 2026-08-09 — a defect reported without saying who caused it

A test failure on real data was reported to the owner over two replies as *"cutting a
generic region loses a deletion's tail"*, with coordinates, a CIGAR, and a
recommendation. **The cutting was done by the test itself.** Nothing in the pipeline
splits a region; the test split them to manufacture the boundaries a future parallel run
would have. That sentence appeared in neither reply.

> ❌ "**⛦ Stop-and-ask: cutting a generic region loses a deletion's tail.** On tomato
> SRR7279481 the cut walk yields 7,424,467 generic loci against the uncut walk's
> 7,424,484 — seventeen fewer, none gained."
>
> Owner: *"cutting why? is there a str there? why cut in one run and not in another? I am
> missing context in your text, so it's impossible for me to understand and to decide
> anything."*

> ✅ "**I do the cutting — it is my test, not the caller.** The design asks for proof that
> a genome split across parallel workers gives the same counts as one worker doing all of
> it. To test that I chop each region into thirds, so the generator meets boundaries a
> single-threaded run would never see. Nothing in the pipeline chops regions today."

Four failures, and the first is the one that made the rest unreadable:

1. **The agent was the cause and the sentence had no agent.** *"Cutting"* was written as
   something that happens, so the owner reasonably asked what in the system does it. Every
   noun in that headline — *cutting*, *generic region*, *tail* — was mine, and one of them
   named an action only my test performs. **New rule: when the reported behaviour is
   something your own tooling does, the subject of the first sentence is you.**
2. **"Why in one run and not another" had a one-word answer that was never given:** I had
   changed the test between the two runs (halves, then thirds). Reporting a difference
   across runs without saying what differed reads as a property of the data.
3. **The owner's question "is there a str there?" was answerable and unanswered** — there
   were two repeat tracts, 6 and 8 bases, bracketing the region, and they were why the
   region was short enough for my cut to land inside a deletion. It took one command to
   find out and should not have needed asking.
4. **The skill was never loaded.** Same as 2026-08-04 below, in a worktree where
   `CLAUDE.md` *was* present and does name this file as mandatory. So the symlink fix of
   2026-08-04 was necessary and is not sufficient: the file being reachable does not make
   it read. **Rule for the next reader: load this file when the task begins, not when a
   reply is being drafted** — by the time a result exists, the working that made it opaque
   has already happened.
### 2026-08-10 — a decision asked for in the code's own words

A checkpoint summary ended by handing the owner a "decision" written entirely in names
from the implementation:

> ❌ "`admit` matches pre-screened intervals back to their rows with a linear search per
> interval, so it is quadratic in a contig's row count. Correct, and invisible on the
> fixtures, but it is not the shape for a 90 Mb chromosome with millions of rows. A sorted
> merge is the fix; I would rather let E2's measurement on tomato say whether it matters
> than optimise blind."
>
> Owner: *"your lousy writting is extremely trying. If you want a decission from me I need
> to understand what you're saying to me."*

Four names doing load-bearing work, none of them the reader's: `admit`, *pre-screened
intervals*, *rows*, *a sorted merge*. The repair is to say what the machine does:

> ✅ "For each contig, the code takes the list of repeats the file holds, filters it, and
> then — for every survivor — searches the *whole original list again* to find that
> repeat's stored details. On a chromosome with a million repeats that is a million
> searches through a million entries."

**And it was never the owner's decision.** Nothing about it changes the design: it is an
implementation choice inside the coder's own latitude, and the only reason it looked like a
question was that the alternative ("measure first") had been dressed up as a trade-off. The
rule the reply broke is the D-test: *a paragraph asking for a decision must be answerable
without the reader asking you anything first* — and before that, it must be a decision that
is actually theirs. **Two checks, in this order:** is this mine to decide? If yes, decide it
and say so in one line. If no, describe the thing in the reader's language, not in the
identifiers.

### 2026-08-04 — the rules existed and were never loaded

`CLAUDE.md` already contained a section named *"Writing for the reader — including in
chat"* describing this exact failure. It is gitignored, so it exists only in the main
checkout, and all this work happened in a worktree. **Rule for the next reader of
this file:** when a writing rule is broken repeatedly, check first whether it was
ever in context. Fixed by symlinking `CLAUDE.md` into each worktree; a new worktree
still needs the link.

### 2026-08-04 — a recommendation whose first sentence answered a sub-question

The owner asked which of three places a depth-shedding fix should live. The reply
opened:

> ❌ "**Recommendation: shed in ng's own admission path, and touch neither copy of
> the locked file.**"
>
> Owner: *"have you read the reporting-in-chat skill. If you have, what has failed?
> why haven't you given me the context to understand the first fucking sentence of
> the conversation?"*

> ✅ "**When too many reads are already open, refuse the next one and count it,
> instead of aborting the run.** That is the whole fix, and it needs no edit to the
> file that is locked byte-identical to production's."

Three failures, and the skill as written would have caught none of them — the pass
step covering this was *"read the first sentence alone. Is it the answer?"*, which
the writer always answers yes, knowing which question they meant.

1. *"ng's own admission path"* was coined while working: one of three options, named
   as if the reader had seen the other two described.
2. The first sentence answered *which file the code goes in* — the sub-question the
   working settled last — when the question asked was *what do we do about a walk
   that aborts*.
3. The largest fact, that production fails on the same input in 0.45 s, sat in the
   third paragraph because the working found it third. The owner had even written
   *"worth confirming — it changes how big this is"*, which is a ranking
   instruction, and it was read as a task rather than as one.

**Fixed structurally, not by resolve:** the *Which question* and *lead with the one
that changes what they do* sections above, and steps 0, 2, 3 and 6 of the pass —
each of which compares the draft against the reader's own quoted words instead of
asking the writer how they feel about it.

## Sources

The diagnosis is not project-specific and the external literature is unusually
consistent about the fix being structural rather than motivational:

- **Curse of knowledge** — expertise builds mental shortcuts that run below conscious
  awareness, so experts skip steps and use jargon without noticing.
  <https://mitsloan.mit.edu/ideas-made-to-matter/curse-knowledge-why-experts-struggle-to-explain-their-work>
- **Style guides are hard to enforce** — roughly 2 in 10 people with a guide know it
  well enough to apply it while editing.
  <https://www.acrolinx.com/blog/why-are-corporate-style-guides-so-hard-to-enforce/>
- **"Not a motivation problem but a systems problem"** — guidelines fail at
  interpretation, at enforcement (review is too late), and at learning.
  <https://experienceleague.adobe.com/en/perspectives/brand-consistency-at-scale>
- **BLUF / the Pyramid Principle** — conclusion first, then support, then data.
  <https://en.wikipedia.org/wiki/BLUF_(communication)>
