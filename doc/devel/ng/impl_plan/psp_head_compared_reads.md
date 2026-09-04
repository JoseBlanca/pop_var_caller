# psp record head — the keep rule's denominator: implementation plan

*Draft, 2026-09-04. Turns the settled design in
[`../spec/psp_head_compared_reads.md`](../spec/psp_head_compared_reads.md) into build order —
no new design here. It runs **inside the `ng-psp-mode` branch**, beside
[`run_driver_psp_mode.md`](run_driver_psp_mode.md), and its one hard ordering constraint is
against that plan: **everything here lands before that plan's Milestone F**, while a head
layout change still costs nothing (spec §3, "the window").*

## Scope

**In:** the `reads-compared-with-reference` head field end to end — writer, reader, the two
checks, the fixtures, the owning documents — and the re-measurement of the head's cost.

**Out:**

- **any consumer of the field.** The cheap-numbers read and the two-phase source are the
  psp-mode performance plan's ([`run_driver_psp_mode.md`](run_driver_psp_mode.md), *Out of
  scope: next plans*), which is being specced in the owner's other conversation.
- **format versioning.** Not needed inside the window; after Milestone F it would be, which is
  why the ordering constraint exists.

## Principles (how the order was chosen)

- **The format change is one commit, own commit** — its failure mode is silent (a wrong
  denominator is a wrong keep verdict in the *next* plan, not a panic here), so it lands alone
  with its round-trip and consistency oracles green, bisectable.
- **Documents follow the code they describe** in the same milestone, so no reader of the owning
  specs meets a five-field head that no longer exists.
- **Measure after, never guess** — the head-cost figures in the specs are superseded by this
  change; the plan reruns the probes that produced them rather than adjusting numbers by
  argument.

## Preconditions (already in place — verified 2026-09-04)

- The head's field list and its layout guard:
  [`RECORD_HEAD_FIELDS`, `src/ng/psp/record.rs:138-144`](../../../../src/ng/psp/record.rs);
  `RecordLayout::from_manifest` refuses a file declaring anything else (`record.rs:133-137`).
- The derivation, already computing the denominator and discarding it at the writing line:
  [`record.rs:1719`](../../../../src/ng/psp/record.rs), calling
  [`non_reference_and_compared_reads`, `locus_generation/mod.rs:264-277`](../../../../src/ng/locus_generation/mod.rs).
- The head-against-body check the denominator joins:
  [`decode_the_body_of`, `record.rs:1992-2002`](../../../../src/ng/psp/record.rs).
- The byte-pinning fixture that must move deliberately:
  `the_fixture_encodes_to_these_exact_bytes` ([`record.rs:2282`](../../../../src/ng/psp/record.rs)).
- The cost probes: [`examples/ng_psp_head_encoding.rs`](../../../../examples/ng_psp_head_encoding.rs)
  (head bytes, varint vs fixed, from a production `.psp`) and
  [`examples/ng_psp_skip_value.rs`](../../../../examples/ng_psp_skip_value.rs) (the skipping
  walk on a real `.ngpsp`).
- [`run_driver_psp_mode.md`](run_driver_psp_mode.md) **Milestone F not yet reached** — confirm
  before starting; if F has landed, stop and take the versioning question back to the spec.

## The steps

*One milestone, lettered **H** so its steps cannot be confused with
[`run_driver_psp_mode.md`](run_driver_psp_mode.md)'s A–G in the same conversation.*

### Milestone H — the head carries the denominator

**H1. ✅ The two head changes, end to end. Own commit, do not bundle.**
`RECORD_HEAD_FIELDS` gains `("reads-compared-with-reference", FieldEncoding::Varint)` between
`non-reference-reads` and `record-body-byte-count`, and the `locus-kind` tag moves in from the
body (spec §3.1) — `put_kind`/`read_locus_kind` relocate to the head codec, the `SsrDetail`
stays body-side, `BODY_FIELDS` loses the tag entry. `RecordHead` gains both fields;
`write_a_record` keeps both halves of the derivation and writes them; `read_record_head` reads
them and refuses a head whose numerator exceeds its denominator (the head-only check, a new
malformed-field case); `decode_the_body_of` compares the denominator against the body it
built, beside the numerator's existing comparison, and takes the kind from the head it is
handed. Fixtures updated with the change: the pinned bytes, the layout-refusal test gaining
the missing-field cases, `every_locus_kind_round_trips`, round trips over records with partial
witnesses (where the two counts differ) and with no complete observation (where both are
zero). Oracle: the whole psp suite green, the refusals provoked in both directions.
*Depends:* —. *Source:* spec §3, §3.1, §6.

**H2. ✅ The owning documents follow.** May share H1's commit; if split, H2 lands immediately
after. [`psp_file_format.md`](../spec/psp_file_format.md) §4.3: the field joins the head
diagram and the field-by-field list, and the measured-cost paragraph gains one line saying its
figures predate the sixth field and H3 re-takes them.
[`psp_record_encoding.md`](../spec/psp_record_encoding.md) §2.3: the head line updated.
[`run_streaming.md`](../spec/run_streaming.md) §3.3: the ⚑ paragraph marked **met**, pointing
at the spec. [`cohort_merge.md`](../spec/cohort_merge.md) §13: the position-summary bullet's
field list gains the denominator.
*Depends:* H1. *Source:* spec §1 (the owning specs), §2 (the ⚑).

**H3. ☐ The cost, measured.** Extend `ng_psp_head_encoding.rs` with the sixth field (its
corpus is production records, whose allele support gives the same approximation the numerator
already uses — state the approximation in the output as the probe's header already does for
`record-body-byte-count`) and rerun on the two corpora the specs quote (the tomato accession
at three reads a position, HG002); rerun `ng_psp_skip_value.rs` on a freshly written store.
Update every quoted head-cost figure (the 9.2 % / 5.8 %, the 0.077 bytes a record) where it
stands, each with the date it was re-taken. Report bytes a record before and after the field.
*Depends:* H1. *Source:* spec §4.

> **Checkpoint H: the head answers the keep rule at every depth, the file says so, the cost is
> a number. Pause for review.**

## Verification summary

| milestone | proven by |
|---|---|
| H — the field | round-trip suite green; head-only and head-vs-body refusals provoked both directions; layout refusal names the field; pinned-bytes fixture updated deliberately; head cost re-measured on both corpora, not adjusted |

## Out of scope (next plans)

- **The field's consumer** — the cheap-numbers read over folded heads, the two-phase
  `ObservationSource`, and the rest of psp-mode performance: the successor plan being specced
  in the owner's other conversation, slotted by
  [`run_driver_psp_mode.md`](run_driver_psp_mode.md)'s *Out of scope* list.
- **A format version story for post-F changes** — owed by whichever change first needs one.
