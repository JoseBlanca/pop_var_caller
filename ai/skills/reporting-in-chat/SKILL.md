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
   meaningless to the reader. Code names are not reader names.

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
