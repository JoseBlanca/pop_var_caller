# ng — the locus witness representation, Milestone A: review fixes applied

*Fix-application report, 2026-07-30. Review:
[ng_locus_witness_representation_a_2026-07-30.md](../reviews/ng_locus_witness_representation_a_2026-07-30.md).
Implementation:
[ng_locus_witness_representation_a_2026-07-30.md](ng_locus_witness_representation_a_2026-07-30.md).
Two commits: `0f43317` (the tests) and `89030aa` (the vocabulary).*

Every finding is accounted for below — applied, deferred with a home, or disputed with a
reason. Nothing was silently dropped.

---

## 1. Findings table

| ID | Title | Decision | Where |
|---|---|---|---|
| M1 | `witness_order` (pileup) has no test that can distinguish its two components | **Applied** | `0f43317` |
| M2 | the STR copy is untestable by the one test that names its order | **Applied** | `0f43317` |
| M3 | the fabrication deliverable is asserted nowhere | **Applied** | `0f43317` |
| Mi1 | a public error variant still names the removed type | **Applied** | `89030aa` |
| Mi2 | ~40 sites where a doc and the item it documents disagree | **Applied** | `89030aa` |
| Mi3 | `RecordWitness`'s doc states two false things; its field names a deleted variant | **Applied** | `89030aa` |
| Mi4 | ~40 prose sites still call a witness "coverage" | **Applied** | `89030aa` |
| Mi5 | `SequenceObservation`'s doc contradicts itself | **Applied** | `89030aa` |
| Mi6 | the generic dump's TSV header is pinned by no test | **Applied** | `0f43317` |
| Mi7 | "one read, one observation in its observation" | **Applied** | `89030aa` |
| Mi8 | `ssr.rs`'s "cell" assertion — the string A6 claimed it moved | **Applied** | `89030aa` |
| Mi9 | "observation" collides with `num_obs` in a differential message | **Applied** | `89030aa` |
| Mi10 | none of the three anchors discriminates on the renamed surface | **Applied** (the census invariant) | `0f43317` |
| Mi11 | two `witness_order` functions with identical bodies | **Deferred → B1** | plan |
| Mi12 | three `pub(super)` items with no consumer outside their file | **Disputed** | below |
| Mi13 | `super::super::` vs crate-absolute imports of the witness vocabulary | **Deferred → B1** | plan |
| Mi14 | the commit messages' accounting is wrong in four places | **Applied** (corrections recorded) | §3 |
| Mi15 | three copies of `witness_label`, already drifted | **Deferred → D4** | plan |
| Nits | articles, stale locals, retired test names, assertion messages | **Applied** | `89030aa` |
| Nit | comment wrap (46 lines pushed past the files' wrap) | **Deferred**, with a reason | §4 |

## 2. What was applied, and how it was verified

### The tests (`0f43317`)

The three Majors were one shape: **the suite could not see the surface the milestone
renamed.** Every test was written against a mutation and then run under it — a test that
passes is not evidence, a test that fails under the defect it names is.

| test | mutation | result under it |
|---|---|---|
| `witness_order_ranks_partials_by_offset_before_length` | the comparator's two `u16` components exchanged | 276 passed, **1 failed** — this one |
| `tally_orders_two_partials_of_one_sequence_by_offset_before_length` | same, on the STR copy | 276 passed, **1 failed** |
| the census's per-run invariant | `Partial`'s two fields exchanged where `witness_of` builds them | 271 passed, **6 failed**, naming `seed 0x5eed0001 case 3: locus 19: a partial run 2+0 is empty or reaches past its own locus of 5 positions` |
| the fabrication floor | `measure_fabrication` stops measuring | fails with "class 1 fired on 1484 loci but the fabrication deliverable is 0 reads / 0 bases" |
| the generic dump's header assertion | the `read_witness` column renamed to a placeholder | fails on the column line |

Before these, **every one of those mutations left the suite green.**

The census invariant is an **assertion, not a divergence class**, and it is placed in
`classify_locus` before the classification: a run outside its own locus is not a difference
from production, it is ng being wrong on its own terms, and it must be checked at every
locus of both census passes rather than in a fourth walk.

### The vocabulary (`89030aa`)

Named items, the public error variant, `RecordWitnessCounts`, the "coverage" prose sweep,
and five sentences the substitution broke. The commit message enumerates them. Two points
worth repeating here:

- **`RecordSpanExceedsCoverageRun` → `RecordSpanExceedsWitnessRun` was an open question**,
  not a mechanical call: one reviewer argued for deferring it to D4 because the `#[error]`
  text is user-visible and Milestone A forbids moving output bytes. Resolved by renaming
  now — the contract's oracle is the *generators' emitted data* (the dump TSV, the locus
  stream), no test asserts the message text, and A2 already moved a TSV column header under
  the same plan step. Recorded so the decision is visible rather than assumed.
- **`RecordWitness` → `RecordWitnessCounts` was already sanctioned**: the arch doc's
  boy-scout note asks for exactly it "when the file is next touched".

## 3. Corrections to the record (Mi14)

Four statements in the milestone's own commits and report are wrong. The commits are not
rewritten; the corrections live here and in the implementation report.

1. **A2 says "No expectation was edited."** One was, and it is the predicted one: the header
   literal in `ng_ssr_loci_dump`'s `render_emits_the_spec_9_header_and_tsv_rows`. A string
   inside `assert_eq!` is not an identifier and does not move when a rename tool runs — it
   was edited by hand. This matters because the plan's tripwire is "any test that changes
   expectations did more than rename": the sentence disarmed the tripwire instead of naming
   the one benign instance.
2. **The implementation report says "The one output change in the whole milestone is the STR
   dump's column header."** Two headers moved — the generic dump's as well. The plan
   authorises both ("the `read_coverage` column in **both** dump tools' TSV output"); what
   was missing is that the checkpoint oracle only covered the STR dump. Mi6's new assertion
   closes that.
3. **A6 says "five sentences were rewritten"** and names three; the other two are
   `generator.rs`'s `total_obs` doc and `mod.rs`'s "tally cell" → "tally **bucket**" — the
   latter a third vocabulary choice made silently. And "two article fixes" was a mechanical
   ~23 a/an flips.
4. **A6 claims `ssr.rs`'s "one row per (allele, read group) cell" moved.** It did not — it
   is an assertion message and A6 was comments-only. It was the last "cell" in the
   observation sense in the tree, and it is fixed in `89030aa`.

## 4. Deferred, each with a home

- **Mi11 — the two `witness_order` copies (`open_record.rs`, `ssr.rs`), byte-identical.**
  `open_record.rs`'s comment justifies withholding an `Ord` impl because it "would export
  *this file's* sorting convention to every other consumer"; the STR tally independently
  invented the same convention, so the order is the type's, not the file's. A reviewer
  verified that deriving `PartialOrd, Ord` on `ReadWitness` and reducing both call sites
  leaves 275 passing. **Home: B1**, which creates `witness.rs` and moves the type — the move
  should absorb the comparator rather than carry two copies across. Noted in the plan.
- **Mi13 — `super::super::` imports of the witness vocabulary** in three files against
  crate-absolute in two, where in `open_record.rs` the same spelling resolves to two
  different modules. **Home: B1**, for the same reason: after the move,
  `grep crate::ng::locus_generation::witness` should answer "who depends on this".
- **Mi15 — three copies of `witness_label` across the example dumps**, identical down to a
  shared seven-line comment and already drifted (`partial:interior` in files whose other arms
  use underscores). **Home: D4**, which rewrites the labels; the drift should be decided
  there, not inherited.
- **The comment re-wrap.** The substitution pushed ~46 comment lines past the files' wrap and
  `cargo fmt` does not reflow comments. An automated pass was written, run, and **reverted
  whole**: at a 96-column bound it re-wrapped 233 lines — most of them lines the milestone
  never touched — and left orphan words ("— and / the"). Doing it properly means either a
  width the files actually use or a hand pass over the touched blocks; either is churn that
  would bury the vocabulary diff. Best done with Milestone B, which rewrites the same
  comment blocks.

## 5. Disputed

- **Mi12 — demote `witness_of`, `ObservationKey` and `KeyedObservation` to private.** The
  finding is accurate: all three are `pub(super)`, nothing outside `open_record.rs` uses
  them, and demoting compiles clean. It is nonetheless declined, because
  [arch §2](../../ng/arch/locus_witness_representation.md) specifies `pub(super) fn
  witness_of` in the interface this plan is implementing and C2 changes its signature there.
  The visibility is a design decision, not rename residue. Recorded so it is not re-found.

## 6. Validation

Run in the container from **this worktree's** `scripts/dev.sh`, after each of the two
commits:

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --lib --bins --tests --examples --all-features`: **2,808 passed, 0 failed**
  (2,806 before the review; +2 unit tests — the census invariant and the fabrication floor
  are assertions inside existing tests).
- `cargo test --release --lib ng::locus_generation`: **277 passed, 0 failed, 1 ignored.**
- The STR dump on the tomato CRAM (`SRR7279503.p1.bench.cram`, `SL4.0ch01`): **byte-identical
  to A2's output**, i.e. the whole review-fix pass moved no emitted byte.
- The pre-existing `benches/psp_writer_perf.rs` failure is unchanged and unrelated.

## 7. One process note

The fix pass itself produced the run's clearest example of why "minimal diff discipline"
is a rule and not a preference: an automated comment re-wrap, run because a Nit asked for
one, touched 233 lines across five files — five times the finding's scope — and was
reverted by restoring the files and re-applying the vocabulary changes as a single scripted
pass. The scripted pass is in `tmp/vocab_pass.sh`, which is scratch and not committed.
