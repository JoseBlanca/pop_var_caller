# ng calling foundations — B1 + B2: the `calling/` folder and its first two shared types

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, on top of A2 (`f33c5198`).
Steps B1 and B2 of [`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md),
Milestone B.*

## 1. Plan, and why two steps share one loop

**B1 and B2 ran as one implement-review-fix-commit iteration, named here and in the commit.** B1's
whole content is a module file with a doc comment and a `pub mod calling;` line — a stage a review
would have nothing to say about, and a commit that could not carry a test. The
plan-driven-implementation skill allows tightly-coupled adjacent steps to share one iteration when
one of them would leave a stage near-empty, provided the merge is explicit rather than silent. B3
stays its own iteration.

- **B1** — `src/ng/calling/`, the one folder for steps 6 to 9, wired into `src/ng/mod.rs`.
- **B2** — `CandidateAlleles` and `ExpectedAlleleCopies` in `calling/mod.rs`.

Design authority: [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 and §Module
home; [`arch/module_layout.md`](../../ng/arch/module_layout.md) for the folder and principle 1b;
[`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §1.3 and §9 for what expected allele
copies are and why they travel in the output.

## 2. Assumptions and deviations

**`CandidateAlleles.alleles` is private, where the architecture sketches a public `Vec`.** The arch
writes `pub struct CandidateAlleles { pub alleles: Vec<Box<[u8]>>, pub kind: LocusKind }`. The plan
restates the contract the sketch cannot hold — "REFERENCE at index 0, always present" — and asks
for a test of it. A public `Vec` can be emptied, reordered, or seeded with the reference somewhere
other than the front, and nothing in the type would notice, so the field is private and the four
things a caller needs are methods: `new(reference, kind)` (the only constructor, seeding index 0),
`reference()`, `alleles()` (the flat view), `get(AlleleId) -> Option<&[u8]>`, and `push(allele) ->
AlleleId`. `kind` stays public: it is a fact about the locus and carries no invariant of this
type's.

This is the same trade A2 made for `Genotype`, and for the same reason — an invariant stated only
in prose is one a later step breaks. It is a deviation from the sketch's *shape*, not from its
contract; recorded here and flagged to the reviewers rather than absorbed silently.

**`get` is a checked lookup rather than an index, and that is a promise A1 made.** A1's review
recorded it as a follow-up for this step: `AlleleId`'s doc comment says an out-of-range id is
caught when the table is read, and until now nothing read a table. An id carries no locus, so an id
minted at one locus is a legal `u16` at the next; indexing would panic at best and, on a longer
table, silently return a different locus's allele.

**The final prune is not built.** The plan says the table is owned "because a discovery round
appends and the final prune shrinks it — a later plan's behaviour, this plan's shape". `push` is
built because without it the ownership argument has nothing behind it and no test can add an
allele. The prune is not, because *which* alleles to drop is a policy the calling loop owns
(`max_candidate_alleles`, `arch/calling_em_loop.md` §2.1), and a `retain` written here would have
to guess it. The loop plan adds the method it needs.

**`calling/mod.rs` does not yet declare `genotype_table`,** which the plan's B1 text asks for. That
file is built by C1 of this same plan, and declaring a module whose file does not exist does not
compile. The declaration lands with the file.

**Nothing checks that `ExpectedAlleleCopies` is the same length as its locus's allele table.** The
two are parallel by contract and neither type can see the other. Stated in both doc comments;
flagged to the reviewers as the assumption most worth challenging, since a length mismatch is a
silent wrong answer rather than a crash.

## 3. Changes made

Two files, **+260 / −1** (`git diff --stat`).

- **`src/ng/calling/mod.rs`** (new) — the module doc explaining why steps 6 to 9 are one folder,
  plus the two types and six tests.
- **`src/ng/mod.rs`** — `pub mod calling;` and one sentence in the crate module doc naming what has
  landed.

## 4. Tests added

Six, in `src/ng/calling/mod.rs`:

| test | what it pins |
|---|---|
| `the_reference_allele_is_id_zero_and_stays_there` | The plan's named invariant: a table cannot exist without its reference, and two `push`es do not move it. Also that `push` returns ids 1 and 2, in order. |
| `the_flat_view_is_in_id_order_with_the_reference_first` | The order a scorer walks is the order `get` resolves and the order the expected copies run in — checked by resolving every id the view offers. |
| `an_id_this_table_does_not_hold_resolves_to_nothing` | The A1 follow-up: an id from another locus resolves to `None`, where an implementation that indexed would panic or return the wrong allele. |
| `the_locus_kind_travels_with_the_allele_table` | `kind` survives construction on both the generic and the STR-bundle path, so neither becomes the only one that compiles. |
| `expected_copies_are_fractional_and_kept_as_given` | Copies come back bit for bit — a clamp or a round would show. |
| `expected_copies_cannot_be_built_parallel_to_nothing` | The empty-vector panic, with its message. |

## 5. Validation

Run in the dev container (`./scripts/dev.sh`, by its absolute path in this worktree), verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 3.84s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |

*(Those are the numbers as submitted for review, and the six tests below are the six it saw. The
review added six more and renamed three methods; the post-fix state is in
[the fix-application report](ng_calling_foundations_b12_fixes_2026-08-21.md).)*

The three aggregate gates that are red on `main` — `clippy --all-targets`, `test --all-targets` and
`doc --no-deps --lib`, all in files this branch does not touch — are unchanged by this step and are
recorded in [A1's review](../reviews/ng_calling_a1_2026-08-21.md) §7.

## 6. Trade-offs and follow-ups

- **The parallel-length contract between the two types is unenforced** — see §2.
- **The prune is the loop plan's** — see §2.
- **`genotype_table` is declared at C1** — see §2.
- **Nothing consumes either type yet.** `LocusInference` carries both, and that arrives at B3.
