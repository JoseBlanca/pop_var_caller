# Bad writing, caught by the owner — examples and criticisms

*Working document, started 2026-08-03. Not a skill yet. It collects sentences the owner sent
back, with what he said about them, so that the patterns behind them can eventually become rules
in [`clear-technical-writing`](../clear-technical-writing/SKILL.md) or
[`spec-authoring`](../spec-authoring/SKILL.md).*

**How to use it.** Add an entry every time a sentence gets sent back. Record the owner's words
verbatim — his criticism is the ground truth, and paraphrasing it loses the thing that made it
land. Then work out what the sentence was doing wrong and whether an existing skill should have
caught it. When a pattern has three or four entries it is probably ready to become a rule.

**Why it is worth keeping.** Every entry below was written *while following*
`clear-technical-writing`. So the interesting column is not "what was wrong" — it is **what the
skill does not currently say** that would have stopped it.

---

## 1. The word that names a category instead of its contents

**Written:** "`read/filtering.rs` holds **policy only** — no file reading, no conversion, no
loop."

**Owner:** *"bad writting. 'holds policy only' means absolutely nothing."*

**What was wrong.** "Policy" was reached for as shorthand and never unpacked. The file holds nine
filters, their thresholds, the name of each drop reason, and two functions that apply them. Every
one of those is concrete; "policy" is the bag they were put in so they would not have to be
listed. The reader cannot act on a bag.

It went unnoticed because the sentence *sounds* precise — it has a strong claim ("only") and
three negations after it. The negations were doing real work; the noun they attached to was not.

**Became:** "`read/filtering.rs` holds **the keep-or-drop rules and the thresholds they use, and
nothing else**: it opens no files, converts nothing, and drives no loop over reads."

**Test to apply.** Replace the word with the list of what it contains. If you cannot write the
list, you do not know what you meant. If you can, the list is the sentence.

**Other words in the same family, watch for them:** machinery, infrastructure, concerns,
plumbing, surface, semantics, story, shape.

**Skill gap.** `clear-technical-writing` Rule 4 covers *borrowed* jargon that decodes to nothing
("kernel"). It does not cover the plain English word that is *too general to be information*.
"Policy" is not jargon; it is a category label standing where the category's members belong.

---

## 2. The term the document turns on, never defined

**Written:** "**No change to what the cursor keeps.** It keeps decoded, filtered reads…"

**Owner:** *"in order to understand this a reader needs to know what's a 'filtered reads'. Is it a
RawRecord or a AlignedRead?"*

**What was wrong.** The whole document was about the boundary between an undecoded record and a
decoded read. It used "record" for one and "read" for the other, consistently, and **never once
said that was the convention.** Writing it felt fine because the distinction was clear in my
head; the reader has only the page.

**Became:** the distinction defined in §1 before the diagram that uses it, and the sentence
changed to name the type and where it lives — `VecDeque<AlignedRead>`, held above the whole
chain, never raw records.

**Test to apply.** If a distinction is load-bearing, define both sides *before* the first
sentence that leans on either. A distinction you are confident about is the one most likely to go
undefined.

**Skill gap.** `clear-technical-writing` Rule 3 says define every acronym and term on first use.
Neither "record" nor "read" is an acronym or a term of art — they are ordinary words — so the
rule as written did not fire. The rule should reach *any* word doing load-bearing work whose
meaning is narrower here than in general use.

---

## 3. The claim that overstates in order to sound rigorous

**Written:** "**No performance change sought.** If a measurement moves, something is wrong."

**Owner:** *"that's plain stupid and completely unnecessary. Getting a simpler and more
understandable architecture is a good thing, no need to improve performance. Regressing
significantly in performance would be bad and we would had to reconsider, but maybe with the new
architecture we even improve performance, that wouldn't be wrong at all."*

**What was wrong.** Two things, and the second is the worse one.

The claim is **literally false**: a faster result would not be wrong. It was written to sound
disciplined — *we are so careful that any deviation is a defect* — which is a pose, not a fact.

And it framed speed as a **goal the change had to hit**. It is a **constraint**, and an
asymmetric one: an improvement needs no defending, a significant slowdown means the design has to
be reconsidered. Stating it as a target invites someone to tune for it, which is not what the
change is for.

**Became:** "**Speed is a constraint here, not a goal.** … A *significant* slowdown would be a
reason to reconsider the design; an improvement is welcome and needs no defending. Neither is
what this is for. Either way, measure it." Plus a number in the verification section for it to be
checked against, marked as one run on one machine.

**Test to apply.** Read the claim as an absolute and ask whether the opposite outcome would
really be bad. If not, you have written a pose. Then ask: is this thing a *goal* or a
*constraint*? Constraints are usually asymmetric, and saying which direction is fine is more
useful than forbidding both.

**Skill gap.** `spec-authoring` has *"Facts, not arguments"* and *"never manufacture a
rationale"*. This is the sibling failure it does not name: **manufacturing a standard** — an
invented acceptance criterion, stricter than the truth, adopted because strictness reads as care.

---

## 4. The defensive lead — arguing with a critic who is not in the room

**Written:** "**It is not a matter of taste — it falls out of what each filter reads.** Six of the
nine read two integers off the raw aligned read. The other three read fields that only the
conversion produces."

**Owner:** *"bad writting, this is in the best of cases an empty sentence, but it is worse than
that because it is difficult to understand."*

**What was wrong.** Three faults compounding:

1. **"It is not a matter of taste"** pre-empts an objection nobody made. It argues rather than
   informs.
2. **"falls out of"** is a metaphor doing the work a verb should do. Falls out how? It gestures
   at "is determined by" without committing to it.
3. **"what each filter reads"** forward-references the very thing the next two sentences and the
   table are about to establish. At that point in the page the reader does not yet know filters
   read different things, so the clause explains nothing and asks to be held in memory.

The compound effect is what the owner named: not merely empty, but *harder to read than
nothing* — the reader stops to decode a sentence that had no content to deliver.

**Became:** the sentence deleted. The section now opens on the fact, and the heading ("Why the
division falls where it does") already told the reader what is being answered.

**Test to apply — the delete test.** Cut the sentence and read the passage. If it is not worse,
the sentence was not carrying anything. Apply it hardest to any sentence beginning *"It is not
just…"*, *"This is not merely…"*, *"Far from being…"* — the contrastive opening is where a
strawman usually hides.

**Skill gap.** `spec-authoring` warns against *inventing a strawman alternative in order to
defeat it* — but at the scale of a **decision record**. The same move at the scale of a
**sentence** is not covered, and it is far more frequent: a lead-in that defends the paragraph
instead of starting it.

---

## 5. The missing head noun — and the antecedent a deletion carried off with it

**Written:** "Six of the nine read two integers off the raw aligned read."

**Owner:** *"More bad writting. Six of the nine what?"*

**What was wrong.** The noun is absent. "Six of the nine" needs a head — six of the nine
*filters* — and the sentence never supplies one.

**Two things made it, and the second is the one worth learning from.**

*First*, the noun was carried by the sentence in front of it. The passage originally opened
"…what each **filter** reads", so "six of the nine" had something to bind to. That was the
defensive lead in entry 4. **Deleting it orphaned the sentence after it** — the cut was right and
the re-read was not thorough enough. Entry 4's own test says *cut it and read the passage*; the
cut happened and the reading did not.

*Second*, even with the noun restored the sentence leaned on a back-reference it had not earned.
"The nine filters" assumes the reader knows there are nine and which nine. This document
mentions the number twice before §2, and both times as another document's content — so a reader
arriving at §2 knows a count and nothing else. And "two integers" was coy: the table three lines
below names them, so there was no reason for the prose not to.

**Became:** "Step 1 runs nine filters; `read_filtering.md` §3 is where they are defined. **Six of
them read two values that are on the raw aligned read already** — the SAM flag and the mapping
quality. **The other three read fields that only the conversion produces.**"

**Two tests to apply.**

- **After deleting a sentence, read the next two.** A deletion can remove the antecedent for a
  pronoun, a bare number, or a partitive ("six of the nine", "the other three", "both", "the
  latter"). The delete test is not finished when the sentence is gone.
- **A bare count or partitive needs its noun *and* its set.** "Six of the nine X" fails if the
  reader has not been told there are nine X and where they are defined. Introduce the set in the
  same breath or point at where it lives.

**Skill gap.** `clear-technical-writing` Rule 7 asks for short active sentences and one idea per
sentence. It does not say that **a sentence must survive being read on its own** — which is the
property both faults here break, and the property most easily destroyed by an edit somewhere
else. Nothing in either skill covers *damage caused by a correct deletion*.

---

## 6. Explaining a design by ruling out the alternative, instead of saying what it does

**Written:** "**So 'all nine filters run on the raw aligned read' is not reachable.** It would
need either a second copy of the mismatch rule and the CIGAR scan written against noodles' types
— the failure this module guards against hardest — or the uppercase and the CIGAR conversion run
before the filter, which is the conversion under another name."

**Owner:** *"better reasoning. converting a rawalignedread to an alignedread has a cost, so we
filter what we can in the rawlignedreads, then we convert and finally we filter the
alignedreads."*

**What was wrong.** The paragraph closed the section by defeating an option, so the reader was
left holding a *negative*: they now know one arrangement will not work, and still have to infer
the arrangement that does. The owner's sentence carries the whole design in one line — a cost, and
the ordering that follows from it — and needs no alternative at all.

**Note what is *not* wrong with it.** The option really had been raised, in conversation, so this
is not the strawman of entry 4. The fault is the shape of the explanation, not its honesty: a
design explained by elimination is weaker than the same design explained by its reason, even when
every word of the elimination is true.

**Became:** "**Converting is not free, and that is what fixes the order.** Building an aligned
read copies the name, uppercases the sequence, rebuilds the CIGAR as ng's own operations and
works out the adaptor boundary. So: reject on what the raw aligned read already carries, convert
what survives, then reject on what the conversion produced. **A read dropped by the first six
never pays for a conversion.**" The constraint that follows — that the last three cannot move
earlier without a second copy of two rules — is kept, but demoted to the paragraph after, where
it is a consequence rather than the argument.

**Test to apply.** When a section ends by saying what cannot be done, ask what the *reason* for
the actual design is and lead with that. If the reason exists, the elimination is at best a
footnote; if you cannot state the reason without the elimination, you have not found it yet.

**Skill gap.** `clear-technical-writing` Rule 1 is *explain before you formalize* — plain English
before the formula. There is no companion rule for **explain positively before you eliminate**,
and it is the same instinct: the reader wants the thing, not the boundary around the thing.

---

## 7. Dressing a decision up as a necessity

**Written:** "**Both alternatives to the current placement are wasteful, and they fail in
opposite directions**, which is what makes the answer forced rather than preferred."

**Owner:** *"bad writting. You're tying to be sensationalist. We are not forced to do anything,
for all we care we could print out the bam file and do the anylisis manually. We *prefer* to have
a simpler architecture and maybe even an improved performance. So given our goals we decide
what's the best approach."*

**What was wrong.** "Forced" was reached for because it sounds stronger than "chosen". It is
false — nothing forces any of this — and the falseness is the smaller half.

**The real damage is that it deletes the goals from the page.** A decision presented as
inevitable carries no criteria, so a reader cannot tell what it was optimising for, and therefore
cannot tell whether it still applies when the situation changes. `spec-authoring` says a decision
recorded without its rationale rots. *"Forced"* is worse than a missing rationale: it tells the
reader there is nothing to reconsider.

It also quietly insults the alternatives. (a) and (b) would both produce correct output. They
lose on a cost, against goals this document states — which is a comparison a reader can check and
disagree with. "Forced" is not.

**Became:** the lead now says there are three placements, all producing the same reads, differing
in the work they spend. The conclusion says *"So we choose (c)… Given what this change is for — a
simpler shape, and no significant slowdown — (c) is the one that costs nothing to get. (a) and
(b) are not wrong; they would produce the same output, and if the goals were different they might
be the better trade."*

**Test to apply.** Whenever a design reads as the only possibility, ask: *what would have to be
true for one of the others to win?* If there is an answer, the design is a **choice** and the
sentence owes the reader the goals it was chosen against. If there genuinely is no answer, say
what makes it impossible — that is a fact, and it will be specific.

**Watch for:** forced, inevitable, the only way, has to be, no choice but, cannot be otherwise.
Each is worth one check before it survives a revision.

**Skill gap.** `spec-authoring` requires alternatives-considered *and* the reason the rejected
ones lost. It does not warn about the inverse move: keeping the alternatives on the page and then
overstating the verdict, so the section looks like it weighed something while telling the reader
the weighing was a formality. That reads as more rigorous than an honest "we preferred this, for
these goals", and it is much less useful.

---

## Patterns so far

| # | pattern | one-line test |
|---|---|---|
| 1 | a category word standing in for its contents | replace it with the list; if you cannot write the list, you do not know what you meant |
| 2 | a load-bearing distinction never defined | define both sides before the first sentence that leans on either |
| 3 | an invented standard, stricter than the truth | would the opposite outcome really be bad? is this a goal or a constraint? |
| 4 | a lead-in that argues instead of informs | delete it and re-read; watch every *"it is not just…"* |
| 5 | a count or partitive with no head noun, or no set to belong to | read the sentence alone: "six of the nine **what**, out of which nine?" |
| 6 | a design explained by ruling out the alternative | state the reason for the design; if you cannot, you have not found it |
| 7 | a choice presented as a necessity | what would have to be true for another option to win? if there is an answer, it is a choice — name the goals |

**A common root for 1–4, provisionally.** All four are sentences written for how they would
*sound* to someone judging the document, rather than for what a colleague needs from it.
`spec-authoring` already says the bar is not "could this survive review", because that puts an
adversary in the room. These are what that adversary does to individual sentences — and the skill
only names the damage at the level of sections and decisions.

That may be the rule worth extracting: **the adversary-in-the-room failure is a sentence-level
disease, not only a document-level one.**

**Entries 3 and 7 share a third root: overstatement that reads as rigour.** One invents a
standard stricter than the truth, the other promotes a preference to an inevitability. Both make
the document sound more disciplined and leave the reader with less — in 3 an acceptance criterion
nobody should hold to, in 7 a decision with its criteria removed. The tell in both is a word
doing emphasis rather than work: *wrong*, *forced*.

**Entry 5 is a different root again, and worth keeping separate.** Nothing about it is defensive — it
is a sentence that could not stand alone, half of it caused by fixing entry 4. That points at a
second candidate rule: **every revision pass should re-read the neighbours of anything it
changed**, because the failure a good edit introduces is invisible to the person who made it.
They still have the deleted sentence in their head.
