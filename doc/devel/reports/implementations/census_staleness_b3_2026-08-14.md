# B3 — which pileup a census came from, and what to do when it is not this one

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step B3 — implementation plan 2.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§6.1, §7.13.
**Date:** 2026-08-14.

---

## 1. What landed

`PileupIdentity::of_header(header, records)` — the digest of a pileup's header beside its record
count — and `freshness(named, in_hand)`, which says what a run should do with a census file:

| | the pileup in hand |
|---|---|
| the census names the same one | **Fresh** — use it |
| it names a different one | **Rebuild**, naming what differs: the header, or the record count |
| it names one and there is none here | **Refused**, naming the pileup it was built from |
| it names none at all | Rebuild where a pileup is here, Refused where none is |

**Nothing reads a modification time**, which is the point: a modification time changes when a file
is copied and does not change when a file's contents are rewritten in place, so it answers a
question nobody asked.

## 2. What this step deliberately stops short of, and why

**Rebuilding is a verdict here, not an action.** Building a census from a pileup is the second
producer, which is milestone C's C1 — there is nothing to rebuild *with* yet.

**And ng has no pileup file.** `src/ng/locus_generation/pileup/` is the *locus generator* for
ordinary positions, not a file format; nothing in ng writes a pileup today. So the header bytes
this digests are supplied by the caller rather than read off a pileup, and which bytes exactly is
the pileup writer's business when it exists. What `of_header` promises is only what spec §6.1 asks
of it: two pileups with the same header and the same record count get the same identity, and no
others do.

That is the honest seam. The alternative — inventing a pileup header format here to digest — would
be designing the other half of the two-phase run inside a staleness check.

## 3. Tests

| test | what it pins |
|---|---|
| `a_census_naming_another_pileup_is_rebuilt_where_it_can_be_and_refused_where_it_cannot` | all four verdicts, with the two mismatches told apart by *which* value differs — a changed analysed-region set (the likeliest accident, and part of the header) against records added under the same header |
| `touching_a_census_changes_nothing_about_which_pileup_it_names` | the census file is rewritten byte for byte, re-opened, and names the same pileup with the same verdict |

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::parameter_estimation::joint::census_file` | `12 passed; 0 failed` (10 before) |
| `cargo test --lib` | `3,604 passed; 0 failed; 11 ignored` (3,602 before) |

**No oracle run.** This step adds a check nothing calls yet: `grep` finds `freshness` only in its
own module and its tests, so the walk and the fit cannot have moved.
