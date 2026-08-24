# ng read likelihoods — B1: where a wrong read's probability goes

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step B1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone B, on
top of `bed8bafd`.*

## 1. What it is

`src/ng/calling/likelihood/generic.rs` and the first piece of real arithmetic in this plan:
`m(a, g)`, how many things a read could have shown given that it is wrong.

**3.0** where the observation differs from **every** allele the genotype carries by a substitution
at exactly one position, **1.0** otherwise.

## 2. Why there is a divisor at all, and why it is not always three

A read the genotype cannot produce is wrong. But *wrong how?* If the individual carries `A` at this
base and the read shows `C`, the chance of that particular misread is not the chance of any
misread — there were three bases it could have gone to. That is the physical fact, and it is what
the parameter pre-pass's own noise model already assumes: three bases to go wrong into, one to come
back to.

**The second case is where the model runs out.** Insertions, deletions and multi-position
differences have no finite set of things a wrong read could have shown, so any divisor would be
invented. Leaving the mass unspread is conservative, and conservative here means favouring the
reference — the direction a caller should err in when it has nothing to say.

**"Every allele the genotype carries", and the word is load-bearing.** A read one substitution from
one of a heterozygote's alleles and an insertion away from the other has no clean three-way spread
either. That is what makes `m` a property of the `(allele, genotype)` pair rather than of the
allele pair alone, and it is the one place an `any` would look just as reasonable and be wrong. A
test builds exactly that genotype — reference plus a deletion, against an observation one
substitution from the reference — and pins 1.0, beside the same observation against the reference
homozygote, which is 3.0. Under an `any` rule the first would be 3.0 and every other test in the
file would still pass.

## 3. Why this is its own commit

The plan says so, and gives the reason: **a wrong divisor is `log 3` per wrong read, in the wrong
direction, and nothing crashes.** Asserted rather than quoted: `log 3` is 1.0986 nats, which is
4.77 on the Phred scale.

The direction matters as much as the size. Dividing by three makes a wrong read *less* probable, so
calling a read wrong costs more, and the divisor therefore **favours the heterozygote** — the
opposite direction from the multinomial coefficient this plan's next step drops. The three vendored
callers disagree about this and so do the two halves of production: GATK divides by three,
freebayes and production's SNP path divide by nothing, and production's own STR substitution term
divides by three.

## 4. Two decisions inside the implementer's latitude

**The pair table.** "One substitution apart" is a property of the allele pair, so it is answered
once per locus into an `A × A` table and looked up in the genotype loop rather than recomputed for
every genotype that carries the allele. At six candidates that is 36 booleans against 21 genotypes
× 6 alleles of lookups.

**One spelling of the index.** `divisor_at` exists because a table indexed genotype-major and one
indexed allele-major are read as each other silently, and the result is the wrong reads divided by
three. The layout matches `GenotypeTableView::genotype_allele_counts` row for row, so a reader
holding one row of counts holds the matching row of divisors at the same offset — one convention,
not two.

## 5. What the tests pin

Fourteen tests. The ones that guard something no other test does:

| test | the defect it fails on |
|---|---|
| `a_genotype_carrying_an_indel_refuses_the_spread_for_every_observation` | **`all` written as `any`** — the rule's central word, and the only test that can see it |
| `an_indel_is_never_one_substitution_however_short` | the length check dropped, which would make a two-base deletion whose bases are a prefix of the reference look like a substitution |
| `an_allele_the_genotype_carries_is_not_one_substitution_from_itself` | *exactly one* written as *at most one*, which would also give an identical allele pair the spread |
| `the_table_is_genotype_major_and_matches_the_counts_row_for_row` | a transposed fill — and its fixture is deliberately not square, three alleles against six genotypes, because a square one cannot see it |
| `the_rule_reads_which_alleles_are_carried_and_not_how_many_copies` | a rule that read copy counts rather than membership, at a tetraploid |
| `the_predicate_is_symmetric` | a pair table filled in one triangle only |
| `the_spread_is_worth_one_point_one_nats_per_wrongly_explained_read` | the size claim in the prose going stale |

**The fixtures name genotypes by what they carry, not by an index.** `genotype_carrying` finds the
row with the copy counts a test means, so a test reads *the heterozygote carrying the reference and
the deletion* rather than *genotype 4* — which would be a fact about the table's ordering rather
than about the rule.

## 6. What an independent oracle says

**A reviewer implemented `m(a, g)` a second time from the specification's words alone and swept
the two against each other over 1,758,811 `(locus, ploidy, genotype, allele)` cells across 19,188
loci at ploidy 1, 2, 3 and 4. Zero disagreements.**

The oracle is independent in the way that matters rather than a paraphrase: it decides "one
substitution apart" by **set membership** — building the explicit set of sequences one substitution
from an allele and asking whether the observation is in it — so *exactly one, not at most one* and
*never across lengths* are structural in it rather than conditions it checks. It enumerates the
genotype as a multiset of allele ids of its own rather than reading the shipped table's copy
counts, and it converts each shipped row back to a multiset before comparing, so a genotype-major
table read allele-major is a disagreement rather than an invisible reshuffle.

The grid covers equal-length pairs differing at 0, 1, 2 and 3 positions, insertions and deletions
of 1, 2 and 3 bases, duplicate alleles, and up to six alleles — 126 genotypes at a tetraploid.

## 7. What the reviews changed

**The genotype coordinate was a bare `usize`, and no assertion could have saved it.** `divisor_at`
typed the allele and left the genotype loose, so an allele index handed in where the genotype
belongs compiled — and because a shape always has at least as many genotypes as alleles, it was
*always in range* and returned a real but wrong divisor. That is this step's own named failure
shape, arriving through the function written to prevent it. It now takes `GenotypeIdx`, which is
what the sibling table takes for the same reason.

**Two assertions had no test, and one that did could not fail.** The buffer-length check survived
being weakened from "exactly" to "at least" — and that equality is what makes `divisor_at`'s bound
meaningful, since a longer table admits a genotype index that should have been out of range.
`divisor_at`'s own range check survived deletion entirely, returning a neighbouring genotype's
cell. Both now have a test.

**The test written for the length check was the one test that could not catch its removal.** All
three of its pairs were prefix or suffix relations, and the comparison `zip`s — which truncates to
the shorter sequence and therefore sees *zero* differing positions on exactly those pairs. The case
that bites is `ACGT` against `ATG`: truncation leaves one difference, so without the check it comes
back a substitution. The doc comment also argued from padding, which the code does not do.

**And a test that was a test of `f64::ln`.** The size claim computed `3.0_f64.ln()` from a literal
and never touched the table — it would have passed with the fill deleted. It now takes the
difference between a spread allele's divisor and an unspread one's, at the same locus and genotype,
so the number and the code that produces it fail together. **The plan's B2 differential has to
avoid the same shape**: if the `÷3` effect is written as a literal there, the differential passes
with this step deleted.

**A repeat tract can no longer be handed to this function.** Its substitution term is a different
rate on a different model, so a divisor table filled for one is a number with no meaning reaching
the wrong row builder.

**And the accessor could not enforce the one thing it existed for.** It took the allele count as a
loose argument, so nothing tied it to the stride the buffer was actually filled at — and a reviewer
measured what that costs: reading a three-allele diploid table at a stride of two returns a real
divisor from the wrong row on **six of twelve lookups, with nothing to panic about.** That is the
step's own named failure shape, arriving through the function written to prevent it, and no
arrangement of assertions closes it. It is now a `DivisorTable` that carries its own stride and is
built against the genotype view the fill used, so there is nothing left for a caller to supply
wrongly. It hands out a genotype's whole **row** as well as one cell, which is the shape the row
function wants — its inner loop already holds the matching row of copy counts.

The crate had argued this case against itself twice already, and B1 had not listened:
`CandidateAlleles::bases_of` returns an `Option` because indexing "would hand back a real but wrong
allele without complaint", and `GenotypeIdx` carries the same warning about a row meaning different
genotypes at different shapes.

## 8. What this step costs, measured

**About 15 million loci a genome** — the two anchors agree: HG002 chromosome 21 gives 251,792
generic loci over 46.7 Mb, extrapolating to 16.6 M; tomato `SL4.0ch01` gives 1,718,914 over 90.9
Mb, extrapolating to 14.8 M.

| alleles | ploidy | genotypes | fill | without the pair table |
|---|---|---|---|---|
| 2 | 2 | 3 | 44 ns | 50 ns |
| 6 | 2 | 21 | 626 ns | 1,478 ns |
| 6 | 4 | 126 | 2,350 ns | 10,360 ns |
| 16 | 4 | 3,876 | 288 µs | 615 µs |

The pair table earns 2.1× to 4.4× everywhere above two alleles, and is a wash at the biallelic
case. **Nothing here is quadratic in a way that bites**: at the worst shape the fill is about 3% of
the row function it feeds, and at six candidates and a diploid it is a tenth of a percent. The
quadratic that hurts is the genotype count itself, and that belongs to the genotype table.

**The per-call allocation is 3 to 7 ns of a 44 ns fill** — about a seventh, or a tenth of a second
per genome. It is contract-legal where it is (the no-allocation rule is about the per-*sample*
loop, and this is one level outside it) and it is not a performance problem. It should move into
the caller's scratch anyway, but as a rider on naming that scratch's owner rather than for its own
sake.

## 9. Two things recorded for the next step

**The divisor table's owner had no name anywhere.** It belongs in `CallingScratch` as an eighth
field, because it is per *locus* and that is the type allocated once per worker and reused per
locus — **not** in the row's own per-sample scratch, which would invite a refill per sample of a
quantity that does not vary by sample. It is refilled once per locus and not once per pass, and it
is generic-path only. Recorded in arch §3.

**Whether the table should store `m` or `log m` is B2's to settle, and it has to settle it.** The
row charges `log m` once per `(observation, genotype)`, so as it stands B2 calls `.ln()` in its
inner loop — measured at 1.392 ns a term against 0.553 ns if the table held the logarithm, about 26
seconds single-threaded over a high-depth sample's loci. It cannot be deferred past B2 because
storing the logarithm makes *divisor* the wrong word: nobody divides by 1.0986.

## 10. Validation

In the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,250 passed / 0 failed / 14 ignored**; the likelihood module holds 63,
  of which 22 are this step's.
- `cargo doc --no-deps` — 23 unresolved links, the same 23 that are on `main`.
