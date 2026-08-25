# ng read likelihoods — A1: the evidence, and the two allele tables nobody had told apart

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step A1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone A, on
top of `b10dc612`.*

## 1. What it is

`src/ng/calling/likelihood/` — the module that will hold step 7 — and the three types that say
what one sample's evidence at one locus **is**, before anything scores it. No arithmetic: the
plan makes this step types only, and the scoring lands in Milestones B to H.

- **`GenericObservation`** — the merge's fold of every read that showed one allele from one read
  group: which allele, which read group, how many reads, and Σ `ln P(error)` over them. Four
  numbers, where the merge's own row carries eight.
- **`GenericSampleEvidence<'a>`** — those rows, plus the pooled quality of reads matching no
  candidate at all, plus the partial observations with their witnessed stretch intact.
- **`SsrSampleEvidence<'a>`** — the locus generator's own observations at a repeat tract, plus
  that tract's motif and flanks, split by witness.

## 2. The finding that changed the step: two allele tables, one numbering assumed

**A first version of this step offered `GenericObservation::of_supported(&SupportedAllele)` — one
merge row in, one observation out — and it was wrong in a way nothing would have reported.**

The cohort merge unifies every distinct sequence the whole cohort showed into one table, and its
rows index that table. Candidate *selection* — step 6, which does not exist yet — keeps some of
those alleles and drops the rest, and **dropping allele *k* renumbers every allele above it**;
`CandidateAlleles`' own documentation says so, and says the method that does it has to return the
remapping. `AlleleId` is documented as an index into the surviving table. So the two numberings
agree exactly until the first prune, and after it every observation above the dropped allele is
scored against the wrong sequence.

The first version's own doc comment argued that the conversion had to be written once because
"the mapping is where the two shapes can silently drift" — and then guarded the wrong drift, the
`usize` to `u16` width, while assuming away the one that changes which allele a read is scored
against.

**What shipped instead makes the mapping an argument.** `of_supported_allele(row, allele)` takes
the candidate id; `fill_from_supported_alleles(rows, candidate_of_merge_allele, out)` takes
selection's whole mapping, `&[Option<AlleleId>]` indexed by the merge's allele index. A row the
mapping drops is not silently skipped: **the fill returns the pooled `q_sum` of the rows it
dropped**, because those reads are evidence and the data likelihood is compared between loci
(spec §3.3). That return value is *part* of `unmatched_q_sum` and not the whole of it — selection
pools its own leftovers too, and the specification that says how is not written. The doc says so
rather than letting a caller mistake one for the other.

**Two design documents were factually wrong about this and are corrected in this commit.** Arch
§2.1's sketch carried the comment `allele: AlleleId, // a(o) — resolved by the merge's byte
unification`, which is the assumption in one line. And arch §2.2 plus the plan's own preconditions
said the generator's `complete_observations()` "stays the only unguarded access" — neither half of
which is true: the field it reads is `pub`, so it is a helpful name and not an enforcement, and
this step spells the same split again on a view that holds a bare slice.

## 3. The four numbers, and the four left behind

A merge row carries eight numbers: the allele index and the read group on `SupportedAllele`, and
six on `AlleleSupport` — read count, forward-strand count, summed log error, summed mapping
quality, summed squared mapping quality, and how many reads started left of the anchor. The view
keeps four.

**The four dropped are dropped because no likelihood reads them.** Strand bias, the
mapping-quality multi-mapper test and the read-position-bias term are site filters, and they run
*after* genotyping (spec §1.4). The type is what keeps them out of the row.

**One of the four kept is not the shape it looks.** `q_sum` is a sum of logarithms, not a
probability. That is exactly why the closed form of spec §3.3 charges an unexplained observation
`q_sum` and not `n · log ε̄` — the two are the same number, and only the first reproduces bit for
bit when the same reads are aggregated differently.

## 4. The read group is part of the identity, and the type says so

The formula puts a logarithm outside a sum over alleles. An observation's reads may therefore be
pooled into one term **only if every one of them would have got the same number** — and two reads
showing the same bases from two lanes have different error rates (spec §2.3).

The merge was changed for exactly this (prerequisites B1): its rows are one per
`(allele, read group)`, ascending, and `SampleSupport::pooled_support_for` exists as the named
escape hatch for questions that really are about the sample. Nothing here calls it, and nothing
here *can*: the fill takes one row at a time and never sees a `SampleSupport`.

**The ascending order is checked here as well as promised.** The merge sorts once per sample in
`AlleleTable::assemble`, and its own `the_rows_are_ordered_by_allele_then_read_group` pins that —
but that is a test in another module, and this view is built from a staging buffer that a
selection mapping could in principle reorder. `GenericSampleEvidence::new` holds a debug assertion
that the rows are strictly ascending on the pair, with two tests: one that reverses the alleles,
one that reverses only the read groups within one allele, which a check comparing the allele alone
would let through.

## 5. A partial observation is a different claim, and lives in a different field

A read that entered the locus and ran off its own end does not say what the sample carries; it says
the sample carries **at least** this. Scored as though it were complete it mis-scores as a *short*
allele, because its bases are a prefix of the truth (spec §5.1).

So `partials` is its own field, holding the merge's own `PartialObservation` rows — bases and
witnessed run intact, because spec §5.3's compatibility rule compares an allele's projection
*restricted to the witnessed run* and cannot do that once the run is gone. Arch §2.1 asks for one
type rather than two, and this is it: the view borrows what the merge owns.

**The run and the bases are on different axes and their lengths are not interchangeable** — they
differ by the net indel the read carried — so the fixture here is deliberately built with five
witnessed positions against seven bases. A fixture where the two are equal is exactly the one a
consumer indexing the bases by a locus offset would pass, and the merge's own documentation says
outright that equality does not license indexing one with the other.

On the repeat path the split is a filter rather than a field, because the generator emits one list.
Two shapes there were earned in review rather than planned:

- **The two filters yield `(position in the slice, observation)`.** The STR row's emission cache is
  keyed by that position — the cache is what makes a row cost `observations × candidates` instead
  of `observations × genotypes`, which spec §8 calls the design and not an optimisation. Without
  the position, the cheap way to build that cache is to re-walk the raw field, which is the thing
  the split exists to discourage. Yielding it makes the guarded route the useful one.
- **The partial side is an exhaustive `match` on `ReadWitness`, not `!= Complete`.** A third
  variant is then a compile error at the one place that decides what reaches the censored term,
  rather than joining it silently.

Both return `impl Iterator<…> + use<'a>`. Without `use<'a>`, edition 2024 captures the elided
`&self` and the iterator cannot outlive the view — which the row builder needs it to. A review
mutation proved it: dropping `use<'a>` turns a probe that returns the iterator into `E0597` while
every test still passes.

## 6. What the tests pin

Fourteen tests. What each would catch:

| test | the defect it fails on |
|---|---|
| `the_view_of_a_merge_row_keeps_the_four_numbers_the_formula_reads` | any of the four numbers read from the wrong source field |
| `a_dropped_allele_renumbers_the_ones_above_it_and_its_quality_comes_back` | the identity mapping assumed — merge allele 2 becomes candidate 1 here — and the dropped row's quality discarded |
| `keeping_every_allele_drops_no_quality` | the pooled quality accumulated on the wrong branch (it would return −30.75 rather than 0) |
| `the_scratch_buffer_holds_no_trace_of_the_previous_sample` | the fill extending instead of clearing |
| `a_mapping_that_does_not_cover_the_merges_table_is_a_caller_bug` | a mapping built at another locus, indexed past its end |
| `the_last_allele_the_mapping_covers_is_accepted` | that same bound written one off |
| `the_observation_stays_four_scalars_wide` | a field owning heap arriving unnoticed, which breaks `Copy` and spec §8's no-allocation contract |
| `the_constructor_keeps_the_pooled_leftover_and_both_slices` | the constructor dropping `unmatched_q_sum` — no genotype moves, every QUAL does |
| `a_sample_that_showed_nothing_carries_no_reads_and_no_pooled_leftover` | `empty()` seeded with anything but zeros |
| `rows_out_of_pair_order_are_a_caller_bug` | the order check missing |
| `rows_out_of_read_group_order_within_one_allele_are_a_caller_bug` | that check comparing the allele alone |
| `the_repeat_view_tells_complete_reads_from_reads_that_ran_out` | either filter inverted or returning everything, or the positions wrong |
| `two_observations_of_the_same_bases_split_on_the_witness_alone` | the split written on the bases — *is this shorter than the tract?* — which every other fixture here would still pass |
| `the_two_repeat_filters_partition_every_observation_and_agree_on_its_position` | the two overlapping, or a future witness variant dropped by both, or the positions disagreeing with the slice |

**The fixtures carry awkward values on purpose.** The merge row's four dropped fields are 7, 611,
37,442 and 3 rather than zero, so a mapping that picked one of them up by accident cannot pass —
a review mutation reading `num_fwd` instead of `num_reads` failed with 7 against 17. The partial's
witnessed run and bases are deliberately different lengths. The repeat fixture interleaves
complete–partial–complete, so an order-dependent filter shows.

## 7. What is not here, and one property to carry forward

`unmatched_q_sum` has **no complete producer in this crate**, and that is the plan's decision
rather than an omission: the pooled leftover belongs to candidate selection, which has no
specification yet. The row takes it as an input and tests supply it from fixtures.

**The whole step still reverts green**, and three reviewers each demonstrated it independently: the
module has no consumer outside itself, so deleting the directory and the one `pub mod` line leaves
`cargo build --lib --tests` and the rest of the suite clean. That is inherent to a types-only
scaffold whose stated proof is "type-level compilation as the sibling plans' import surface", and
it is not a defect on its own. What it obliges is that **B2 must build its evidence through
`fill_from_supported_alleles` rather than open-coding the narrowing** — otherwise the mapping
argument that this step exists to make explicit becomes a function nobody calls.

**One half of that was closable now and is closed.** Deleting the single line `pub mod likelihood;`
orphaned the whole module: the crate still built, clippy was still clean, and
`cargo test --lib ng::calling::likelihood` reported `0 passed; ok` — a green run naming a module
that was no longer compiled, which is the exact shape of false green this plan's siblings have been
failing on. `calling/mod.rs` now re-exports the three types, so the same deletion is
`error[E0432]: unresolved import`. The re-export's reason is different from the one `genotype_table`
is re-exported for, and its doc comment says which.

The parameter views, the floors and the row contract are step A2.

## 8. What the reviews changed

Four category agents in their own worktrees — conformance to the design documents, test adequacy
by mutation, API design and reuse, and one written for this step: the two invariants it could
silently destroy. Between them **49 mutations: 38 caught, 9 survivors — 6 distinct once the
overlap between agents is removed, each with a two-output proof — and 2 that changed no
behaviour**. All four also ran a whole-step revert check. Beyond the mapping finding of §2:

- **A survivor that shipping builds would have carried.** The first version narrowed with an
  `assert!` over an `as` cast, and `[profile.release]` sets `debug-assertions = false`: editing it
  to `debug_assert!` left all tests passing while a release build renamed allele 65,536 as the
  reference allele. Measured under `-C debug-assertions=off`: the unmutated code panicked, the
  mutated one returned `Ok(AlleleId(0))`. The cast is gone entirely now, so there is nothing left
  to downgrade.
- **A survivor in the constructor.** `new` discarding `unmatched_q_sum` and storing `0.0` passed
  every test. It is the one term that cancels in genotyping and does not cancel in the data
  likelihood, so nothing but a direct assertion would have noticed.
- **A claim the type could not perform.** The struct doc said the view borrows the merge's rows;
  `partials` does, and `supported` cannot — the merge's row is 48 bytes against this one's 24,
  measured. Hence the staging buffer and the fill, and a doc that names which half borrows what.
- **A false "guard".** Two reviewers, independently, found that the word overstated what a `pub`
  field can promise; the prose now says what the split actually buys.
- **A test that could not have failed**, and one no fixture could reach. The `assert_ne!` opening
  `the_allele_index_and_the_read_group_are_not_crossed` was a tautology given the two `assert_eq!`s
  under it; `SsrSampleEvidence::detail` was stored and never read back. The first is gone, the
  second now asserted.

## 9. Validation

In the dev container, on the committed tree:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib` — green.
- `cargo doc --no-deps` — 23 unresolved links, the same 23 that are on `main`. This change adds
  none.
