# ng psp store — H1: fixes applied

*2026-08-28. Applying [the H1 review](ng_psp_h1_2026-08-28.md) to step H1 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The three Blockers

### B1 — the cross-encoder arm now compares everything production carries, and cannot quietly stop

**Reproduced here first.** A store written with `placed_left` read one too high *and* every
region a base too long passed the whole 74,623-record hg002 run with a clean report. Both defects
are in the harness's own mapping, so the round-trip arm carries them on both sides and cannot see
them — which is the entire reason the second arm exists.

Three comparisons added: `placed_left`, the reference bases, and the record's **extent** against
the reference allele it was derived from. And the structural fix that makes the class impossible
to recur:

**Both comparison functions destructure with no `..`.** A field added to
`SampleLocusObservations` or `SequenceObservation` is now a compile error at both arms rather than
a field the oracle silently stops comparing. The four fields production has no counterpart for are
written `read_witness: _`, `read_group: _` — named, not swallowed. This is the discipline
`encode_record_body` already keeps in `src/ng/psp/record.rs`, for the same reason and in the same
words.

The doc comment that said the opposite now says what is true: four fields have no counterpart on
that arm, and everything else production carries is compared.

### B2 — every comparison in both arms is now held by a test

Two table-driven tests, one per arm, sixteen rows and thirteen rows: each perturbs one field of a
record that agrees and requires the check to refuse it **naming that field**. A comparison added
later without a row is a comparison nothing holds, and the table is the inventory as well as the
guard.

A third test asserts both arms *accept* a record that agrees, so a check that panicked
unconditionally could not be mistaken for one that discriminates.

⚠ **Writing them found a test that could not fail.** `an_allele_sequence_the_prototype_disagrees_with_is_refused`
changed allele **0**'s sequence — which, now that the reference bases are compared against that
same allele, is refused as a reference-bases difference before the per-allele loop is reached. The
test passed without the per-allele sequence comparison existing at all. It changes allele 1 now,
and says why in its doc.

### B3 — CI's clippy command is green, and the report's claim about it is corrected

The two lints move onto the prototype's **own items** — `#[allow(unused_assignments)]` on
`push_with_head`, `#[allow(clippy::too_many_arguments)]` on `encode_streaming` — where they cover
its own target as well as the included copy. A lint attribute changes no behaviour, which is the
ground the visibility widening was already accepted on. `cargo clippy --all-targets --all-features
-- -D warnings`, the command CI runs, is now clean; it had been failing since `b0e1a54a`.

What stays at the `mod` declaration is `dead_code` alone, with a comment saying it is this
harness's doing rather than the prototype's and naming the condition for removing it.

**And the prototype now says it is an oracle**, at the top of its own file, where someone about to
"fix" that dead store will see it.

## 2. The Majors

- **M1 — a fourth witness shape.** Every witness the corpus minted was reproducible by a witness
  stored as one prefix length and one suffix length, the hole included: its two runs sat on the two
  borders. `a_witness_for` now mints a fourth shape — a single run flush with **neither** border —
  `CorpusShape::interior_witnesses` counts it, and `assert_the_corpus_can_fail` requires one. Two
  tests hold the distinction: that the fourth shape touches neither border, and that the first
  three all touch one.
  Measured on the two corpora: **1,118 interior witnesses on hg002 and 2,303 on tomato.**
- **M2 — `Scales` is spelled out**, no `..default()`, so a field added to it has to be decided
  here rather than defaulted to the prototype's value. `Manifest { .., ..as_this_build_writes_it() }`
  is deliberately left as it is, and the reviewer distinguished the two: that constructor is by
  definition the right source for anything the harness does not override.
- **M4 — the header, the trailer and the block index are read back.** The header is compared
  field by field (destructured, no `..`), with the provenance split in two because `create`
  records the compression level into it — so the header in the file is deliberately not the one it
  was handed, and every parameter handed in must still be there. The trailer is compared to the
  marker `finish` was given. The index is checked **two ways**, and only the second has teeth:
  entering by ordinal must give a suffix of the whole walk, *and* each entry's coordinate must be
  the coordinate of that block's first record.
  ⚠ **The second assertion exists because the first did not catch the defect.** With only the
  suffix check, an index whose every entry named a position one base too far **passed** a
  20,000-record run — `records_from` searches the index on the coordinate, and entering by ordinal
  never reads it. With the coordinate check it fails at block 0.
  All three store-level defects are now caught: the contig length by the header check, the index
  coordinate by the new assertion, and the truncated trailer by `PspReader::open` refusing the file
  outright — which is the reader's own guard, not this harness's check, and is recorded as such.
- **M5 — the block count comes from `block_index().len()`**, with the writer's tally kept as a
  cross-check rather than as the source.
- **M6 — `open_production` names the file and says what to do**: a psp written before the current
  column set is the likeliest failure and arrived as a serde error with the whole TOML header
  inlined and no path in it.
- **M7 — the prototype's three format switches are named** at the call site.

## 3. The Minors taken

`--label` in the usage line and the module doc; `report` split into `assert_the_run_proves_something`
and `print_the_report`, so the two run-level claims cannot be lost behind a later `--quiet`;
`write_the_ng_store` returning a named struct; `compare` renamed to
`walk_the_three_streams_in_lockstep`; `at`/`which`/`g`/`w`/`p` replaced by `record_ordinal`,
`pushed` and `read_back`; `READ_GROUPS_IN_THE_CORPUS` shared between the modulus and the tally,
with the tally **panicking** rather than silently dropping a group past its end; the two coprime
moduli given their reason; `note_the_drift`'s parameters named for their streams; `--limit 0`
refused where it is parsed rather than surfacing as *the ng store holds records the source does
not*; the flag `expect` messages saying what the flag needs; the work directory and the genomic
grid printed in the report; `Bp` imported rather than fully qualified; the prototype-drift line
kept with a comment saying it equals the ng line by construction and what its differing would
mean; two untested branches of `assert_the_corpus_can_fail` given tests; the *residual derivation*
message rewritten in plain words.

**Not taken, and why:** the shared-mapping extraction (a `#[path]`-included
`examples/ng_psp_corpus/`) is a follow-up — the reviewer filed it as a note rather than a merge
condition, and it touches `ng_psp_head_encoding.rs`, which is outside this step. The renamed
`as_an_ng_record_with_synthesised_fields` closes the half of it that matters now: two functions in
one directory no longer share one name with deliberately different behaviour.

## 4. My own numbers, corrected

Three wrong claims, all of them about my own work, all found by the verification pass:

| claim | right value |
|---|---|
| "seven declarations became `pub(crate)`" | **eight** — and the sentence's own list had eight items |
| "it has two findings there … allowed at the `mod` declaration" | the attribute names **three** lints and suppresses **24** findings, 22 of them `dead_code` |
| the module doc's "two of them … are synthesised here" | **four** are — and this one inverts the point the report's §2.5 exists to make |

The third is the expensive kind: a reader trusting it would believe the corpus still leaves the
two counts at production's zero, which is the precise defect this step found and fixed.

**Forty-nine other figures were checked correct**, including every cell of the two-corpus table and
both reader-guard mechanisms word for word. The split is the same one this plan keeps producing:
figures quoted from a document are right, figures describing the author's own fixture are where to
look.

## 5. What the corpora say now

Both runs pass, at the shipped `GRID_BP = 10_000`, with every store-level check in place.

| | hg002 chr21 | tomato SRR7279481 |
|---|---:|---:|
| records compared | 74,623 | 7,687,686 |
| observations | 123,263 | 8,049,472 |
| blocks the ng store cut | 598 | 879 |
| chain ids compared | 150,720 | 555,022 |
| witnesses: complete / partial / with a hole / flush with neither border | 30,975 / 92,288 / 514 / 1,118 | 2,012,355 / 6,037,117 / 344 / 2,303 |
| worst summed-log-error against the source | 0.000122064 | 0.000122030 |

*The witness counts moved from the pre-fix run because the shapes went from three to four: one
observation in four is now `Complete` rather than one in three.*

## 6. Left open

- **`LocusKind` is `Generic` on every record of both corpora**, so no parity run has moved the
  `Ssr` arm of the kind codec. It is compared on the round-trip pair, and the module doc now says
  so rather than listing it under "on both pairs".
- **The lockstep walk itself has no unit test**, and neither does the skip arm that keeps the
  prototype's stream moving past a record with no ng record — on both corpora that count is 0, so
  the line never executes. Closing it needs a small production `.psp` committed as a fixture,
  which is a bigger change than this step.
- **The two corpora are not in the repository**, and their provenance is not recorded. An agent
  reviewing this had to survey 700 files to find a readable substitute, and two production `.psp`
  files still under `tmp/` are unreadable by this build.
- **The record limit is applied in two places** — the writing pass and the comparison pass — and
  they agree because both now use `take(limit)`. One definition would be better.
