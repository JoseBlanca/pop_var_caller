# Spec review — `ng/spec/read_preparation.md` (step 2, generic path)

**Date:** 2026-07-25 · **Branch:** `ng-read-preparation` · **Status:** changes requested
**Target:** [doc/devel/ng/spec/read_preparation.md](../../ng/spec/read_preparation.md) (306 lines)
**Bar applied:** [`ai/skills/spec-authoring/SKILL.md`](../../../../ai/skills/spec-authoring/SKILL.md)
— the Revision pass and the two-reader test (does it save the *coder* a day; can the *manager*
tell cost / deps / whether it is still the right thing).

**Code read for grounding:** [src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs),
[left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs),
[delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs),
[src/ng/read/filtering.rs](../../../../src/ng/read/filtering.rs),
[src/ng/ref_seq.rs](../../../../src/ng/ref_seq.rs),
[src/ng/locus_generation/ssr.rs](../../../../src/ng/locus_generation/ssr.rs),
[src/pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs),
[indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs),
[src/pileup/per_sample/read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs),
[baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs),
[baq_stream.rs](../../../../src/pileup/per_sample/baq_stream.rs),
[src/bam/alignment_input.rs](../../../../src/bam/alignment_input.rs).
**Docs read:** [impl_plan/alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md),
[spec/alignment.md](../../ng/spec/alignment.md), [spec/ng_proposal.md](../../ng/spec/ng_proposal.md),
[arch/ng_step_interfaces.md](../../ng/arch/ng_step_interfaces.md).

Severity: **BLOCKER** (must resolve before a coder starts) · **MAJOR** (real defect or gap) ·
**MINOR** (wording / smaller gap). Every finding cites a file:line or a spec §.

---

## 1. Verdict

The consolidation is the right call and the central claim holds up under checking: read
preparation *is* a generic-path-only step, and the STR path really has no read preparation
(verified below). Length (306 lines) is comfortably inside the exemplars. The retirement of the
two path specs is clean, and nothing was silently dropped.

The defects cluster in one place: **§6, the interface** — which is also the one section a coder
would build directly from. Three of its four load-bearing statements are wrong against the code
it cites (where it runs, what the scratch holds, how the read arrives), and its signature cannot
express the error model §7 mandates — repeating, verbatim, the mistake step 1's spec made and its
implementation had to correct. Second cluster: the **reuse map** credits production with a
reusable function that is private and in a different file, while missing the `pub` one that does
the whole job.

And one manager-level gap: the spec never says what v1 *buys*. The repo has already measured
that lever and found it flat.

**Next-slice test.** A coder handed only this spec could not build `LeftAlignPreparer`
unambiguously: they would clone every read, allocate a reference buffer per read, put the work on
the serial walker, be unable to surface a fetch error, and be missing `qname`. All five are
fixable in §3/§6 without changing a decision.

---

## 2. What checks out (verified, not assumed)

| spec claim | verdict | evidence |
|---|---|---|
| Left-alignment is built: `AlignmentNormalizer` + three impls | ✅ | trait [alignment/mod.rs:618](../../../../src/ng/alignment/mod.rs#L618); `StructuredLeftAligner` [left_align_structured.rs:65](../../../../src/ng/alignment/left_align_structured.rs#L65), `FixpointLeftAligner` same file, `RepeatedLeftAligner` [left_align_repeated.rs:84](../../../../src/ng/alignment/left_align_repeated.rs#L84); default alias [mod.rs:655](../../../../src/ng/alignment/mod.rs#L655) |
| The re-align aligner (algorithm 2, affine) is NOT built and is gated | ✅ | only `Output = RepeatSpan` impls exist in `src/ng/alignment/`; plan [alignment_best_path.md:192-206](../../ng/impl_plan/alignment_best_path.md) — Milestone E, **☐ E1**, gate quoted, "*Depends:* D2, and the gate" |
| The STR path has no read preparation | ✅ | `SsrGenerator<R, MF, A: RepeatDelimiter>` [ssr.rs:1337](../../../../src/ng/locus_generation/ssr.rs#L1337) holds the delimiter and calls it per read (`classify_read` [ssr.rs:656](../../../../src/ng/locus_generation/ssr.rs#L656)); no prepared-read type anywhere on that path |
| `RawRefSeq` is the right (and existing) view for v1 | ✅ | [ref_seq.rs:180](../../../../src/ng/ref_seq.rs#L180) — its doc names "the left-alignment / mismatch-fraction path" as the consumer |
| G2/F1 are ng filters #9/#8, and #9 sees the raw CIGAR | ✅ | [filtering.rs:266-271](../../../../src/ng/read/filtering.rs#L266) (#9 before #8), [read_processor.rs:171-175](../../../../src/pileup/per_sample/read_processor.rs#L171) |
| Reassembly is out of scope, not deferred, on a recorded owner decision | ✅ | [ng_proposal.md:477-483](../../ng/spec/ng_proposal.md) ("Partly settled, 2026-07-23") |
| Cross-refs `alignment.md` §4.1 / §4.2 / §6 | ✅ | §4.1 affine best-path, §4.2 repeat-aware, §6 alignment normalization |
| `left_align_indels` rewrites only the CIGAR; the offset does not move | ✅ | called with `remove_deletions_at_ends = false` [indel_norm.rs:412](../../../../src/pileup/walker/indel_norm.rs#L412); ng's wrapper says so and does not move `reference_offset` [left_align_structured.rs:26-33](../../../../src/ng/alignment/left_align_structured.rs#L26) |

The one thing worth adding to that list, because it changes the size of the job: **ng's 1a *is*
production's `left_align_indels`, behind a wrapper** ([left_align_structured.rs:1-8, 50-56](../../../../src/ng/alignment/left_align_structured.rs)).
Left-alignment parity is therefore already banked — see **Mi4**.

---

## 3. Findings

### BLOCKER

**B1 — §6: `-> Option<PreparedRead>` cannot express the fatal error §7 mandates. Step 1's spec
made this exact mistake and the code had to overrule it.**

§7 says a failed reference fetch is "a broken run, **fatal**, surfaced as such, never folded into
a per-read `None`". But the only channel the signature has is `Option`, and the fetch it will call
is fallible: `RawRefSeq::fetch_raw_into(...) -> Result<(), RefSeqError>`
([ref_seq.rs:183-189](../../../../src/ng/ref_seq.rs#L183)). The two statements cannot both hold.

This is not hypothetical. Step 1's spec specified `Item = MappedRead`; the implementation changed
it and left a note saying why: "**The item is a `Result`** (unlike the spec's original
`Item = MappedRead` …)" — [filtering.rs:601-617](../../../../src/ng/read/filtering.rs#L601),
`type Item = Result<MappedRead, ReadFilterError>` [filtering.rs:791](../../../../src/ng/read/filtering.rs#L791),
with the fatal-fetch rule spelled out at [filtering.rs:246-252](../../../../src/ng/read/filtering.rs#L246).

**Fix:** `fn prepare_read(&self, read: MappedRead, scratch: &mut Self::Scratch) -> Result<Option<PreparedRead>, ReadPrepError>`
— `Ok(None)` = no usable observation (tallied), `Err` = the run is over. Same edit is owed in
[arch/ng_step_interfaces.md:261](../../ng/arch/ng_step_interfaces.md), which carries the same
signature.

---

**B2 — §6: "Where it is invoked" is inverted. Production deliberately moved this work *off* the
walker; the spec puts it back.**

§6: "**Where it is invoked — by the pileup, per read, as it walks a non-STR stretch.** Following
production, where `process_read` runs as the walker ingests reads."

Production does the opposite, and says why:

> These steps used to be split: G2/F3/F1 ran serially inside `AlignmentMergedReader::next`, only
> BAQ ran in parallel. Folding them into one function lets the *whole* per-read cost run on the
> worker threads, shrinking the serial floor to the coordinate merge + walker.
> — [read_processor.rs:17-24](../../../../src/pileup/per_sample/read_processor.rs#L17)

The call site is a rayon `map_init` with per-worker state
([baq_stream.rs:331-334](../../../../src/pileup/per_sample/baq_stream.rs#L331), `use rayon::prelude::*`
line 13; also [read_pipeline.rs:191](../../../../src/pileup/per_sample/read_pipeline.rs#L191)). The
walker only *consumes* `PreparedRead`s. §8 of the spec ("embarrassingly parallel") contradicts §6.

This is blocking rather than major because it is the step's placement — the coder builds to it and
the manager schedules from it, and following it re-creates the serial floor production spent a
refactor removing.

**Fix:** state that preparation runs in the **parallel per-read stage**, one `Scratch` per worker
(production's `map_init` shape), and that the walker consumes the output. Keep the "composes with
the gatherer, does not fuse into it" sentence — that part is right and is what keeps the bake-off
surface alive.

---

### MAJOR

**M1 — §3: the reproduced `PreparedRead` is not production's. It omits `qname` and re-types five
fields, while claiming "reused as-is … field for field … reproduced, not redefined".**

Production ([walker/mod.rs:236-295](../../../../src/pileup/walker/mod.rs#L236)):
`chrom_id: u32`, `alignment_start: u32`, `alignment_end: u32`, `mapq: u8`,
`adaptor_boundary: Option<u32>`, and **`qname: Arc<str>`** (line 263) — twelve fields; the spec
lists eleven, with `ContigId` / `Position` / `MapQual` in place of the raw integers.

Both halves cost something. The re-typing tells the coder ng's newtypes apply when the whole point
of reuse is that they do not (reuse means production's `u32`s cross the seam). And `qname` is real
per-read work the manager cannot see: `qname_to_arc` is one `Arc<str>` allocation per read with a
UTF-8-lossy fallback ([baq_engine.rs:434-440](../../../../src/pileup/per_sample/baq_engine.rs#L434)),
on the path §6 refuses a virtual call for. It is load-bearing downstream — the walker keys mate
pairing on it ([walker/mod.rs:225-227, 264-271](../../../../src/pileup/walker/mod.rs#L225)).

It fails loudly rather than silently (`PreparedRead` is deliberately not `#[non_exhaustive]` and
has no `Default` — [walker/mod.rs:228-235](../../../../src/pileup/walker/mod.rs#L228)), which is
why this is Major and not Blocker.

**Fix:** paste the struct verbatim from `walker/mod.rs`, or drop the block and link it. Add one
line on the `Arc<str>` per read.

---

**M2 — §9 reuse map: the named reusable function is private and in the wrong file; the one that is
`pub` and does the whole job is demoted to "model".**

- `mapped_to_prepared` is at [baq_engine.rs:410](../../../../src/pileup/per_sample/baq_engine.rs#L410),
  **not** `pileup/walker/mod.rs`, and is **private** (`fn`, no `pub`) — ng cannot call it. The row
  says "**reuse as-is**".
- `prepare_passthrough` is at [baq_engine.rs:405](../../../../src/pileup/per_sample/baq_engine.rs#L405),
  **not** `read_processor.rs`.
- The fact the spec misses: `prepare_passthrough` **is `pub`, in a `pub mod`**
  ([per_sample/mod.rs:20](../../../../src/pileup/per_sample/mod.rs#L20)), and it does *all* of v1's
  `PreparedRead` construction — `alignment_end`, `mate_role`, `qname`, `mq_log_err`,
  adaptor-boundary propagation, raw-qual copy. So v1 is: fetch the raw window → normalize →
  `prepare_passthrough(read, chrom_id)`. That is a materially smaller job than the spec implies,
  and it makes the parity oracle exact by construction for every field except the CIGAR.

**Fix:** correct both paths; mark `mapped_to_prepared` private/not-reusable; and *decide* the
open question the row hides — call production's `prepare_passthrough`, or re-build it in ng. The
precedent for calling it is already set: `StructuredLeftAligner` wraps production's `pub(crate)`
`left_align_indels` and records the debt ([left_align_structured.rs:50-56](../../../../src/ng/alignment/left_align_structured.rs#L50)).
Whichever way, it is a decision, not a "model".

---

**M3 — §5: the safety argument for moving F1 to step 1 rests on a debug assertion that checks a
different quantity.**

§5: "left-alignment provably preserves the mismatch count — a debug-assert in production's
`left_align_indels` guarantees it — so ng's order gives the identical verdict."

Three problems, all checkable:

1. The assertion is `#[cfg(debug_assertions)]` ([indel_norm.rs:413-421](../../../../src/pileup/walker/indel_norm.rs#L413))
   — compiled out of the build the project runs. The repo is blunt about exactly this elsewhere:
   "**a debug assertion compiles out of the release build this project runs**"
   ([alignment/mod.rs:430-434](../../../../src/ng/alignment/mod.rs#L430)).
2. It asserts over `count_mismatches` ([indel_norm.rs:432](../../../../src/pileup/walker/indel_norm.rs#L432)),
   which is **case-sensitive and unweighted**.
3. F1 thresholds a different number: `read_exceeds_mismatch_fraction`
   ([alignment_input.rs:1052-1100](../../../../src/bam/alignment_input.rs#L1052)) counts only
   mismatches with `q >= bq_floor`, uppercases the reference base, skips non-ATGC on either side,
   and divides by comparable bases.

So the assertion does not establish that F1's verdict is invariant. The argument that *does* work
is already written, in ng's own step-1 code: "left-alignment only shifts indels across equal
bases, so it does not change the match/mismatch tally"
([filtering.rs:239-242](../../../../src/ng/read/filtering.rs#L239)).

**Fix:** use that argument, present it as an argument rather than a guarantee, and name the check
that would settle it — assert the F1 verdict is unchanged before/after left-alignment across the
parity fixture.

---

**M4 — §7 vs §9: ng's fatal-fetch rule diverges from the parity oracle, and the spec does not say
so. The fixture will fail.**

Production treats a failed fetch as *skip F3/F1, pass the read through unchanged* — zero ref span,
unknown `ref_id`, repository miss, `pos == 0`, or a window past the contig end
([read_processor.rs:160-164](../../../../src/pileup/per_sample/read_processor.rs#L160) and
[:181](../../../../src/pileup/per_sample/read_processor.rs#L181);
`fetch_raw_slice` [:86-121](../../../../src/pileup/per_sample/read_processor.rs#L86)). ng chooses
fatal, deliberately and consistently with step 1
([filtering.rs:246-252](../../../../src/ng/read/filtering.rs#L246)).

That is a defensible ng decision — but it means byte-parity against `--no-baq` holds *only* on
reads whose window fetches. On the rest ng aborts the run where production emits an
un-left-aligned read. A coder writing the fixture from §9 alone would treat the first such read as
a bug.

**Fix:** record the divergence beside the oracle in §9 ("parity is claimed where the fetch
succeeds; where production skips, ng fails the run — deliberate, see §7"), and say the fixture must
either exclude those cases or assert ng's behaviour explicitly.

---

**M5 — §6: `Scratch` names the wrong occupant, and the right one is the reason it exists.**

§6: the scratch "holds the normalizer's working buffers now". `AlignmentNormalizer` has **no**
`Scratch` associated type, on purpose: "There is no `Scratch` associated type, because these
algorithms fill no matrix; a normalizer that wants a reusable buffer owns it as its own field"
([alignment/mod.rs:606-611](../../../../src/ng/alignment/mod.rs#L606)) — and `StructuredLeftAligner`
is a unit struct with no buffers at all ([left_align_structured.rs:65](../../../../src/ng/alignment/left_align_structured.rs#L65)).

What v1 actually needs a scratch *for* is the **reference window**: `prepare_read` takes `&self`,
and `fetch_raw_into` writes into a caller-owned `Vec<u8>`
([ref_seq.rs:183-189](../../../../src/ng/ref_seq.rs#L183)). Follow the spec as written and the
coder allocates a `Vec` per read on the billions-of-calls path. Production's shape is the model:
a per-worker `RawContigRefCache` built in `map_init`
([baq_stream.rs:328](../../../../src/pileup/per_sample/baq_stream.rs#L328)).

**Fix:** `type Scratch` = the reused reference-window buffer (plus the aligner matrices when the
re-align mode lands). Drop the normalizer-buffers claim.

---

**M6 — §6: `&MappedRead` forces a per-read clone of `seq`, `qual` and `cigar`, on the path the
same section refuses a virtual call for.**

Production takes the read **by value** precisely to avoid this: "Consumes the read by value so its
buffers move into the `PreparedRead`" ([read_processor.rs:205-208](../../../../src/pileup/per_sample/read_processor.rs#L205));
`pub fn prepare_passthrough(read: MappedRead, chrom_id: u32)` moves `read.cigar`, `read.seq`
straight into the output ([baq_engine.rs:405-431](../../../../src/pileup/per_sample/baq_engine.rs#L405)).
And ng's step 1 already hands back **owned** reads (`Item = Result<MappedRead, ReadFilterError>`,
[filtering.rs:791](../../../../src/ng/read/filtering.rs#L791)), so by-value costs the caller
nothing.

**Fix:** `read: MappedRead`. Same edit in [arch/ng_step_interfaces.md:261](../../ng/arch/ng_step_interfaces.md).

---

**M7 — the manager's question is unanswered: what does v1 buy? The repo has already measured that
lever, and it is flat.**

The proposal calls step 2 "the deepest axis", and records that the only arm still live is
trust-the-mapper vs **per-read re-align** ([ng_proposal.md:477-483](../../ng/spec/ng_proposal.md)).
That arm is the gated, unbuilt mode. Meanwhile the v1 lever has been measured:

> across 63,757 real indel-bearing reads the three produce **identical** output — so the *choice*
> of normalizer does not change calling, and normalization placement is not the lever behind the
> indel deficit … **Owner decision, 2026-07-24**
> — [alignment/mod.rs:646-651](../../../../src/ng/alignment/mod.rs#L646)

So v1 read preparation is **plumbing expected to change no call**; the step's entire experimental
value sits in the mode that is gated on an open question. That is a perfectly good thing to
build — but the spec does not say it, and §2 presents three modes as if the choice were live. This
is the failure the coder cannot see: read as a goal statement, the spec has quietly become "do
what production's F3 does".

**Fix:** one paragraph in §1 — v1 is the control arm, with the screen as the evidence, and the
step's open experiment is the re-align mode. It also answers the manager's cost question: the
price of §4's open question *is* the price of the interesting half of this step.

---

### MINOR

**Mi1 — §2: the "pass through" mode duplicates a fast path the ported code already has.**
`left_align_indels` early-returns on a no-indel CIGAR before touching the reference, allocating
nothing ([indel_norm.rs:409-411](../../../../src/pileup/walker/indel_norm.rs#L409)), and ng's
wrapper says it outright: "the caller does **not** need to pre-filter no-indel reads; this is
already essentially free on them" ([left_align_structured.rs:75-79](../../../../src/ng/alignment/left_align_structured.rs#L75)).
A coder reading "three modes, matched enum" writes a mode enum plus a per-read CIGAR scan that
duplicates production's own. **Fix:** keep pass-through as vocabulary; add the line that v1 need
not implement it as a branch.

**Mi2 — §5 does not say which reference window to fetch, and the two candidates differ.**
Production fetches exactly `[pos, pos + cigar_ref_span(cigar))`, raw
([read_processor.rs:180-181](../../../../src/pileup/per_sample/read_processor.rs#L180)); ng's
normalizer takes the whole stretch and slices from `Alignment::reference_offset`
([left_align_structured.rs:88-95](../../../../src/ng/alignment/left_align_structured.rs#L88)).
Also unsaid: `normalize` operates on `&mut Alignment` ([alignment/mod.rs:626](../../../../src/ng/alignment/mod.rs#L626)),
not on a `Vec<CigarOp>`, so the preparer must wrap the read's CIGAR in an `Alignment` and take it
back out. **Fix:** two sentences — fetch production's exact window, `reference_offset = 0`, wrap
and unwrap.

**Mi3 — §6 names the raw view but not the trap it exists for.** On a soft-masked (lowercase)
reference, `left_align_indels`' case-sensitive comparison shifts indels differently against
uppercased bytes — which is why production carries two reference sources
([read_processor.rs:26-34](../../../../src/pileup/per_sample/read_processor.rs#L26)). Silent
wrong-variant, one sentence to prevent.

**Mi4 — §9 overstates what the parity fixture proves, and understates what it needs.**
"One parity fixture, which also proves left-alignment in isolation" — left-alignment parity is
already banked: ng's 1a *is* production's function behind a wrapper, byte-parity asserted at
Milestone B1 ([left_align_structured.rs:1-8](../../../../src/ng/alignment/left_align_structured.rs#L1)).
What the fixture actually proves is the **window fetch and the field wiring**. And it needs one
setting stated: run production's `process_read(read, None, &mut raw_ref, &cfg)` in-process with
`max_read_mismatch_fraction: None` — otherwise production drops on F1 and ng does not, because ng
moved F1 to step 1, and the keep-sets diverge. The established pattern for the harness is the
`#[cfg(test)]` parity module ([delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs),
`scanner_parity.rs`), which also keeps shipping ng code free of production deps.

**Mi5 — §4/§11 miss a scheduling consequence the plan already recorded.** Milestone E is gated on
*this spec's* §4, and "If the gate is still open at Checkpoint D, **stop there**; this milestone
moves to read preparation's plan" ([alignment_best_path.md:192-197](../../ng/impl_plan/alignment_best_path.md)).
So the affine aligner's cost is now scheduled against this step — and there is no
`impl_plan/read_preparation.md` yet. One line in §11 under the open question.

**Mi6 — §3's "may want hoisting out of `pileup/walker/`" re-opens a settled question.** Production
is frozen (owner, 2026-07-16); the repo already records this exact debt for `CigarOp`, with its
resolution: "porting this module back is the natural moment to do it"
([alignment/mod.rs:66-73](../../../../src/ng/alignment/mod.rs#L66)). **Fix:** point at that
precedent instead of restating the wish.

---

## 4. Against the five questions

| | question | verdict |
|---|---|---|
| 1 | What are we building? | **Yes**, and the scope boundary (§1 locus-independence test, the "does not" list) is the strongest part of the doc. |
| 2 | How does it relate / what does it depend on? | **Partly.** Deps on step 1 and `alignment/` are precise. Missing: where it runs (**B2**), that Milestone E's cost lands here (**Mi5**), and that its consumer (the ng pileup) does not exist yet. |
| 3 | What will bite the coder? | **Weakest.** §6 is wrong on three counts (**B2, M5, M6**), §3 is missing a field (**M1**), and the real traps — the exact window, the `Alignment` wrap, the case-sensitivity, the F1-disabled parity config — are absent (**Mi2, Mi3, Mi4**). |
| 4 | Decided vs open? | **Yes** — §11 is honest, the gating is clear, BAQ's deferral has a home. One decision is hidden inside a reuse-map row (**M2**). |
| 5 | How do we know it works? | **Named but not runnable**, and it will fail as specified (**M4**), while claiming more than it proves (**Mi4**). |

**Manager's re-read of the goals (skill §Revision pass 6):** they have drifted to "do what
production's F3 does" — see **M7**. That is the one finding a coder-only review could not produce.

---

## 5. Suggested order of work

1. **B1, B2** — the two that change what gets built.
2. **M1, M5, M6** — the §3/§6 corrections; all three are mechanical once B2 fixes the placement.
3. **M2, M4, Mi4** — the reuse/oracle cluster; they settle together, and they shrink the slice.
4. **M7** — one paragraph, but it is the one the manager is reading for.
5. **M3, Mi1–Mi3, Mi5, Mi6** — wording and traps.

None of these needs the spec to grow. **M1, M2, Mi4** replace text with shorter, truer text; **M7**
adds a paragraph; the rest are single sentences. Net effect should be roughly length-neutral.
