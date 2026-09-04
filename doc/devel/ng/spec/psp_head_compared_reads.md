# ng — the psp record head carries the keep rule's denominator, and the locus kind's tag

*Status: settled, 2026-09-04 — the owner's decision, taken in the conversation preparing the
psp→calling performance work and routed here because this branch owns the head's encoder. One
new head field, no consumer in this plan's scope. Owning specs this amends:
[`psp_file_format.md`](psp_file_format.md) §4.3 (the head's layout),
[`psp_record_encoding.md`](psp_record_encoding.md) §2.3 (the record). The requirement it
discharges: [`run_streaming.md`](run_streaming.md) §3.3's ⚑. The rule it serves:
[`cohort_merge.md`](cohort_merge.md) §4.3. Implementation plan:
[`../impl_plan/psp_head_compared_reads.md`](../impl_plan/psp_head_compared_reads.md).*

---

## 1. What this is

**Two changes to the record head, one addition and one move.** Added:
`reads-compared-with-reference` — how many of the sample's reads were compared, whole, against
the reference over this record's locus; with it the head answers the cohort merge's keep rule
on its own at every depth, where today it answers it only at low depth. Moved: the
`locus-kind` tag — generic, repeat tract, or bundle — comes forward from the body (§3.1),
because two of the merge's pre-assembly decisions read it and a head that cannot say a
record's kind cannot serve them.

The head is the fixed prefix every stored observation opens with — position offset, reference
span, non-reference read count, body byte count, then the chain-id changes — so a reader can
judge a record and skip its body without decoding it
([`psp_file_format.md`](psp_file_format.md) §4.3;
[`RECORD_HEAD_FIELDS`, `src/ng/psp/record.rs:138-144`](../../../../src/ng/psp/record.rs)). This
spec adds a sixth entry to that list and changes nothing else about the format.

## 2. The problem — the rule needs two numbers and the head carries one

The keep rule decides which loci the cohort builds: a locus is built when **some single sample**
shows at least `max(floor, share × its compared reads)` non-reference reads — floor 2 reads,
share 2 in 100 (`MinAltReads::required_of` and `reached_by`,
[`src/ng/run/cohort_merge/mod.rs:500,526`](../../../../src/ng/run/cohort_merge/mod.rs);
[`cohort_merge.md`](cohort_merge.md) §4.3). Two numbers per sample per locus: the non-reference
count — **in the head** — and the compared-reads count the share is taken of — **today only in
the body**.

At three reads a position the floor is the whole rule and the head suffices. At three hundred
the share is what does the filtering — two non-reference reads in 300 is the sequencing error
rate, and the bar there is six — and it cannot be computed from the head.
[`run_streaming.md`](run_streaming.md) §3.3 flagged exactly this when the requirement was
written — *"a cheap read that returns only the non-reference count cannot answer it"* — and the
settled head never picked the flag up.

**The consumer is the next plan, not this one.** The cheap-numbers read — the cohort reader
that folds heads across samples and builds bodies only for the loci the cohort kept — is
explicitly out of [`run_driver_psp_mode.md`](../impl_plan/run_driver_psp_mode.md)'s scope and
belongs to its successor. The field goes in **now** because of when a format change is free
(§4), not because anything in this branch reads it.

## 3. The decision

**Field `reads-compared-with-reference`, a variable-length integer, between
`non-reference-reads` and `record-body-byte-count`.** The placement keeps the rule's two
numbers adjacent and nothing else constrains it — a soft choice.

**Semantics: exactly the second element of
`SampleLocusObservations::non_reference_and_compared_reads()`**
([`src/ng/locus_generation/mod.rs:264-277`](../../../../src/ng/locus_generation/mod.rs)): the
sum of `num_obs` over the record's `Complete`-witness observations, saturating `u32`. That is
the same subset the numerator is counted over — partial reads and
`reads_without_observation` are in neither half, so **numerator ≤ denominator on every
well-formed record**. A rule taking the share of *depth* instead would raise the bar with reads
that could never clear it, which is why depth was ruled out when the requirement was flagged
([`run_streaming.md`](run_streaming.md) §3.3).

**Derived by the encoder, never supplied** — the same contract as `non-reference-reads`, and
the deriving call already computes both and discards this one at the line that writes the head
([`src/ng/psp/record.rs:1719`](../../../../src/ng/psp/record.rs):
`let (non_reference_reads, _) = record.non_reference_and_compared_reads();`).

**Two checks come with it:**

- **head-only:** a head whose non-reference count exceeds its compared count is malformed. This
  is the one validity check a skipping reader can run without ever touching a body, and it is
  free.
- **head-against-body:** `decode_the_body_of` already re-derives the numerator from the built
  body and refuses a head that disagrees
  ([`src/ng/psp/record.rs:1992-2002`](../../../../src/ng/psp/record.rs)). The derivation
  returns both numbers in one pass, so the denominator joins the same comparison at no extra
  walk.

`RecordHead` gains the field
([`src/ng/psp/record.rs:81-106`](../../../../src/ng/psp/record.rs)), so every consumer of a
head — `records_where`, the walk, the future cheap-numbers read — gets it without decoding
anything.

### 3.1 The locus kind's tag moves from the body to the head

**Two of the merge's decisions read a record's kind before any evidence is assembled**, so a
summaries-first read needs it head-side:

- **the span verdict**: `max_cohort_locus_span` governs *generic* loci only — an STR
  observation's span is its reference tract, which may lawfully exceed the bound
  ([`cohort_merge.md`](cohort_merge.md) §3.1). A reader that cannot tell a tract record from a
  generic one fails every tract wider than 50 bases as an over-wide locus.
- **the never-mix assertion**: a cohort locus must not hold a generic and an STR member at
  once, checked per member at close
  ([`src/ng/run/cohort_merge/close.rs:647-653`](../../../../src/ng/run/cohort_merge/close.rs)).

Today the tag is the body's first-listed kind field
([`BODY_FIELDS`' `locus-kind`, `src/ng/psp/record.rs:220`](../../../../src/ng/psp/record.rs),
written by `put_kind` at `:1093`, read at `:1351`). **The tag moves to the head — one varint,
the same three values — and the body keeps what goes with it**: the `SsrDetail` (motif and
flanks) stays body-side, present exactly when the head's tag says repeat tract, because
nothing pre-assembly reads it. A move rather than a copy, so there is no second copy to check
against the first; the body decoder takes the kind from the head it is already handed
(`decode_the_body_of` receives the located head). Whether the detail's presence is gated on
the tag by a check or by construction is the coder's.

### Why now: the window a format change is free in

`RECORD_HEAD_FIELDS` is the file's fingerprint of its own head layout, and
`RecordLayout::from_manifest` refuses a file declaring anything else — the module's own doc
says what that means for timing: *"It costs nothing today because no psp exists; from
Milestone F it costs a version"*
([`src/ng/psp/record.rs:133-137`](../../../../src/ng/psp/record.rs)). The same shape as the
owner's 2026-09-03 ruling on the header fields ([`run_streaming.md`](run_streaming.md) §6):
no version bump, because the only refusable files are pre-ruling scratch files, and a refusal
naming the missing field is the accepted behaviour. **So this lands before
[`run_driver_psp_mode.md`](../impl_plan/run_driver_psp_mode.md) Milestone F**, after which the
same change costs a format version.

### The alternative that lost

**Over-admit from the numerator alone, and decode the passing sample's body for the
denominator.** Treat "non-reference ≥ floor" as a provisional pass; for the records that clear
it, build the body, read the denominator there, apply the exact rule. It is exact, and it costs
nothing in the file. It lost on two grounds: at three hundred reads a position with a
1-in-1,000 base error rate, error alone puts two or more non-reference reads at about 4
positions in 100 — several times what survives the full rule, each one a body decode the head
exists to avoid, at exactly the depth where a body is largest; and it makes the head unable to
answer the question [`psp_file_format.md`](psp_file_format.md) §4.3 defines it by — *"the two
questions a reader has before it decides whether it wants the record"*. It was the right
answer if the field had had to wait for a version bump; it does not.

## 4. What it costs

**The kind tag is a move, not an addition** — the varint leaves the body as it joins the head,
so its net cost is zero to first order (compression context may move it a little either way;
H3's measurement covers both changes at once). The denominator is the addition:

**Unmeasured — the plan measures it rather than guessing** (plan step H3). What is on the
record and is superseded by that measurement: the head costs 9.2 % of the file at three reads
a position and 5.8 % at 279, and the four scalar head fields compressed to 0.077 bytes a
record measured on their own ([`psp_file_format.md`](psp_file_format.md) §4.3) — figures taken
on a five-field head. The expectation the measurement checks: one more small varint — one raw
byte a record at three reads a position, two at 279 — compressing worse than the numerator,
because it tracks depth and so varies record to record where the numerator is almost always
zero.

## 5. Non-goals, and what this does not do

- **No consumer.** The cheap-numbers read, the fold over heads, the two-phase source — all the
  successor plan's ([`run_driver_psp_mode.md`](../impl_plan/run_driver_psp_mode.md), *Out of
  scope: psp-mode performance*).
- **No change to the keep rule**, its defaults, or where it lives.
- **No version machinery.** The land-before-Milestone-F window is the whole versioning story.

## 6. How we know it works

- **Round trip:** encode → decode equality on the existing suite, the new field and the moved
  tag included — every `LocusKind`, the tract detail surviving in the body
  (`every_locus_kind_round_trips`, [`record.rs:2366`](../../../../src/ng/psp/record.rs)).
- **Both checks provoked:** a head over-reading its denominator refused head-only; a head
  disagreeing with its body refused at build, in both directions.
- **The layout refusal:** a manifest declaring the old five-field head is refused with a
  message naming `reads-compared-with-reference` as missing.
- **The pinned bytes move deliberately:** `the_fixture_encodes_to_these_exact_bytes`
  ([`src/ng/psp/record.rs:2282`](../../../../src/ng/psp/record.rs)) is updated as part of the
  change, not discovered broken.

## 7. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the derivation | `non_reference_and_compared_reads` ([`locus_generation/mod.rs:264`](../../../../src/ng/locus_generation/mod.rs)) | called as-is; the discarded second element is written instead of dropped |
| head plumbing | `RECORD_HEAD_FIELDS`, `RecordHead`, `read_record_head`, `RecordLayout::from_manifest` ([`record.rs`](../../../../src/ng/psp/record.rs)) | one entry, one field, one read, one declaration — no new mechanism |
| the consistency check | `decode_the_body_of` ([`record.rs:1992`](../../../../src/ng/psp/record.rs)) | the denominator joins the existing head-vs-body comparison |
| the kind codec | `put_kind` / `read_locus_kind` ([`record.rs:1093,1351`](../../../../src/ng/psp/record.rs)) | moved to the head's encode/decode; the tract detail's body fields unchanged |
| the cost probes | [`examples/ng_psp_head_encoding.rs`](../../../../examples/ng_psp_head_encoding.rs), [`examples/ng_psp_skip_value.rs`](../../../../examples/ng_psp_skip_value.rs) | extended with the sixth field and rerun on the same corpora |

## 8. Open questions

None. The field's placement inside the head is marked soft (§3); everything else is either
settled above or measured by the plan.
