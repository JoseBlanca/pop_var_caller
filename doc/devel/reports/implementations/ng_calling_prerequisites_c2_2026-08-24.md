# ng — calling prerequisites, C2: partial observations survive collation

**2026-08-24**, branch `ng-calling-prerequisites`. Step C2 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §5.4 and §5.1.

**A read that ran out inside a locus now reaches the built locus carrying the stretch it
witnessed, and still contributes no allele.** Those are the two halves, and they pull in opposite
directions: the evidence has to survive, and it must never be mistaken for an allele.

---

## 1. What changed and why

The merge dropped every observation whose witness was not `Complete` before deriving alleles, so a
built cohort locus carried no partial reads at all and the censored likelihood term of
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §5 had nothing to read. §5.4
(corrected 2026-08-21) makes carrying them a requirement on this module; C1 added the row and left
it empty; this step fills it.

**The one real decision is the coordinate axis.** The mint measures a witness against the *record*
it belongs to, and a cohort locus can hold several of one sample's records and open before any of
them. A row carrying an unshifted run would name the wrong positions, and the consequence is not
noise: at the fixture in the first test, the unshifted run names the locus's first base, where the
read's `A` matches the reference's `A` — so a consumer restricting a candidate's projection to that
position would judge the partial *compatible* with the substitution, which projects to `ACTTA` and
carries `T` at the position the read actually saw. The shift is the difference between compatible
and not.

**That shift is the whole of "projected over that stretch rather than the whole locus span".** A
partial is never padded out to the span — that is what would make it read as a short allele
(§5.1) — it is only placed on the locus's axis so a consumer can restrict a candidate to it (§5.3).

## 2. Changes made

**[`src/ng/run/cohort_merge/build.rs`](../../../../src/ng/run/cohort_merge/build.rs)**

- `MemberPlacement::witnessed_across_locus` — a partial's witnessed runs, moved onto the locus by
  how far the member starts into it. It takes the runs rather than the sequence, so its one `None`
  means one thing: the shifted stretch is not representable. Refusing rather than clamping, because
  a clamp would shorten a witness into a claim about ground the read never saw.
- `partials_of_sample` — one row per `(record, sequence)` with a partial witness and at least one
  read, sorted ascending by `(witnessed stretch, read group, bases)`. Nothing is merged: the mint
  has already pooled every read that showed the same `(bases, witness, read group)`, and two of a
  sample's records are disjoint stretches.
- `assemble` calls it in place of the empty vector.
- Eight tests here, one in the locus generator (C1's), and the shared driver fixture gains a
  partial.

## 3. Deviations from the plan

**The plan says "the projection gains a witnessed-stretch variant instead of the panic (the panic
stays for the code path that must never see one)".** There is no new projection variant and the
panic is untouched. What a partial needs is not a projection of its *bases* — they are carried as
the mint recorded them — but of its *stretch*, and that is what
`MemberPlacement::witnessed_across_locus` is. `compose_into`'s `Complete` assertion and
`projectable_sequences` are both unchanged, so the path that must never see a partial still cannot.

**Partials are gathered by walking the records rather than through `alleles_of_sample`'s
callback.** The derivation's two branches differ in how *reads* are placed and a partial is placed
nowhere, so the walk is the same either way and threading a second emitter through would have
coupled two questions that do not interact.

## 4. What the reviews changed

Three agents, each in its own worktree: what must not have moved, every claim re-measured, and the
strength of the tests under mutation. **56 mutation runs between them, 12 survivors — and the three
agents' survivor lists agree on the same three.**

**The reversion test passes.** Feeding `assemble`'s `partials` field an empty vector while leaving
every test in place reddens four of them. One does not participate and it is worth naming: **the
test that a partial no read is behind yields no row asserted only that the vector was empty**, so it
would have passed on the tree before this step existed. It now carries a second sequence four reads
*are* behind and asserts one row, so it says what it means.

**Three of the row's own properties were asserted against values a broken build would also
produce.** Each was found independently by all three agents, each is fixed by a test that fails when
the property is broken, and each was re-run to confirm the failure:

- **`q_sum` was only ever `0.0`**, which is what the fixture helper hands out — so writing a
  constant zero into every row passed all 250 tests. The fixture now mints `-18.5`; the constant
  fails.
- **The `bases` third of the sort key was never exercised**, so cutting the key down to
  `(stretch, read group)` passed all 250 — on a key that is then not total, in the exact case the
  field's own doc names: a substitution witnessed partially is two rows over one stretch. A new
  four-row test varies all three components and fails when `bases` is dropped.
- **The order of the key's parts was not pinned either.** Reordering it to
  `(read group, stretch, bases)` also passed all 250. The same four-row test fails on that.

**The overflow refusal turned out to be undefended and its doc comment wrong.** Shifting a stretch
that would end past `u16::MAX` refuses rather than clamping, and the comment justified this as
unreachable — *"no locus this module builds can reach: the span is bounded by `MaxCohortLocusSpan`,
50 reference bases by default"*. **Both halves are false.** The span bound is the operator's to set
and holds any `NonZeroU32`, and the locus closer exempts repeat-tract loci from it outright. So the
reachable case is a satellite tract above 65,535 bases, and what happens there is a silent loss of
exactly the evidence this step exists to carry. Measured on a hand-built locus of 70,001 bases: a
partial at offset 65,532 keeps its row, one at 65,534 loses it, one at 66,000 loses it. The comment
now says that; a test pins all three; and both wrong answers now fail it — clamping the end
(which shortens a two-position witness to one) and truncating the offset (which moves a stretch
65,536 bases to the left of where the read was, and which the earlier boundary fixture could not
tell apart from the refusal, because both ends collapsed to an empty run).

**What to do instead of a silent loss is a decision about a failure mode, and it is left open**
rather than taken here — a panic aborts a run over legitimate input, and a count needs a field of
its own, since `reads_removed_as_evidence` means something narrower. Recorded beside the code and in
§6.

**One read at two of a sample's records becomes two rows, and nothing in them said so.** A read that
starts inside one record and runs out inside the next is partial at both, so it leaves one row per
record: two stretches, each claiming one read, for one molecule. §5.3 scores one term per
observation on the understanding that the reads behind it saw the same stretch. The rows carry no
chain id, so a consumer cannot fold them back. **Not fixed — folding them is a design change, and
this step's brief is to carry partials, not to redefine them** — but named in the type's doc and
pinned by a test, so the shape is a recorded fact rather than an accident. There is a gain buried
in it: before this step such a read reached *no* field of the table, not even
`reads_removed_as_evidence`, so the two rows are the first time it is visible at all.

**No driver comparison in this module could see the field, and the claim that they would was
mine.** C1's report said the determinism suite would start covering `partials` "with nothing new
written". It reads through the field — the comparison is on the whole `Debug` rendering — but every
record in every fixture behind it is minted `Complete`, so the field renders empty in every entry
and two drivers agree on it by both building nothing. That matters more after this step than
before: **`partials` is now the only part of a sample's support that can differ between drivers
without `alleles` or `supported` differing too**, because a sample all of whose observations are
partial contributes to nothing else. The shared three-sample fixture's deletion record now carries
one partial, and the width sweep asserts it arrives — 1 sample at 1 locus of the 59 — so the fixture
cannot drift back.

**Four prose corrections in this file.** Three doc comments still described the field as *"empty
until C2 fills it"* and its conventions as ones *"C2 has to establish"*; a test comment listed
adaptor masking as a cause of a holed witness, where the fold's own note says masking truncates from
one side and a hole comes from an interior `N` or a ref-skip; and the projection's doc said a
complete sequence "has the whole locus", where completeness is measured against the *record* and the
whole-locus width comes from the merge's padding rule rather than from the witness.

**One measured efficiency change, and one measured non-change.** Building the member's placement is
now deferred to the first partial in a record, so a record with none pays a scan of the witness
discriminant and nothing else. At 3,000 samples × 10 records × 32 reads with no partial anywhere —
the shape before this step, and the ordinary locus after it — placement was 28 % of this function's
353 µs, against 7,289 µs for the whole assembly. The intermediate `Vec` per row went with it: the
runs are canonical, so checking the last one's end is enough to make the shift plain addition over
an iterator. **The duplicated placement work is not otherwise worth restructuring** — with no
partial anywhere the whole function is 3–5 % of the build, and with half the sequences partial it is
29–43 %, which is the work the step exists to do.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | **253 passed, 0 failed** — the merge's own suite |
| `cargo test --lib` | **4,165 passed, 0 failed, 14 ignored**, 577.19 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 redundant-explicit-link-target warnings — the recorded baseline, unchanged |

The merge's suite held 246 tests before this step and holds 253 after: five new here, one renamed
and inverted from C1's marker, and the driver width sweep gains an assertion rather than a test.

**Every mutation quoted above was applied to this tree and run**, and each was undone with an edit
rather than a checkout, with `git diff HEAD`'s `+`/`-` lines read back afterwards.

## 6. Follow-ups

- **Nothing reads `partials` yet.** The consumer is the evidence view of
  [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md).
- **One read can back several rows, and the evidence view has to answer for it** — either the rows
  carry the read, or a read's stretches over several of a sample's records are folded into one row
  with a hole between them, which is a shape `WitnessedLocusPositions` already expresses. §4 has the
  measurement and the test that pins today's answer.
- **What a locus wider than `u16::MAX` should do with its partials is unresolved.** Today they are
  lost without a word. A panic aborts a run over legitimate input; a count needs a field of its own.
  It belongs with whoever takes the repeat path through the merge, since that is the only path that
  builds a locus that wide.
- **The repeat path's locus-existence amendment is still owed** and is not this plan's:
  §5.4.2 requires one line of the keep rule to change on the STR path, where a sample carrying an
  allele too long for a read to span shows no complete observation and is read as quiet. The rule
  is pinned as it stands by
  `locus_generation::tests::a_partial_observation_is_counted_in_neither_half_of_the_keep_rule`,
  added in C1.
