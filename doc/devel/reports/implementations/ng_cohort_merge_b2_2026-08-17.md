# ng cohort merge — B2: unification into one allele table

*Implementation report, 2026-08-17. Step B2 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §4.2 and the owner's ruling of 2026-08-17, which the
plan carries in full.*

## 1. Plan

Unify the members' projections into one table of distinct alleles, the reference among them
(spec §4.2, *Unification*). Two projections that come out as the same bytes are one allele,
wherever they came from.

**The ruling changes what a sample's projections are**, which is the whole of this step:

> Either we know the read covered the whole locus, and its allele is elongated with what it
> showed; or we know it did not cover it, and it is removed as evidence. Not being able to
> decide which is an error that must never happen.

So a sample's alleles are derived **per read**, not per record. B0 (committed) is what makes
that decidable: every read the generic mint folds now carries a chain id, the
reference-matching ones included, so a read's presence at a record is a fact the merge can
read.

## 2. What "per read" means, record by record

A cohort locus can span several of one sample's records — its SNP at 12 and its SNP at 14
inside another sample's deletion over 10–14. For each of that sample's reads, by chain id:

- **Named at every one of the sample's records**: it showed something at each, and those
  somethings are composed in coordinate order into one allele. The ground *between* two
  records — where this sample minted nothing, because none of its reads departed from the
  reference there — is filled from the locus's reference, which is what "this sample had
  nothing to say here" means.
- **Missing from any one of them**: removed as evidence. It may never have covered that
  position; a depth cap may have discarded it there, the cap acting per position and leaving
  no identities; or it may have been `Partial` there, having seen only part of that record.
  All three are the same fact — what the read showed there was never recorded — so all three
  take the same branch, and the code has one branch for them.

**Where a sample has one record inside the locus, no read is consulted** and each of that
record's sequences is projected on its own, which is B1 unchanged. That is the same answer
the rule above gives whenever ids are present, and it is the **only** answer available on
the STR path, which records no chain ids at all and needs none: an STR locus is one record,
so `ReadWitness` already says whether a read spanned it. It is also the ordinary generic
locus, where every sample has one record.

## 3. Changes made

All in [`build.rs`](../../../../src/ng/run/cohort_merge/build.rs):

- **`AlleleTable`** — `over(&ClosedLocus)` gathers the locus's reference (B1's
  `LocusReferenceBases::over`), seeds the table with it, and then adds what each sample's
  reads showed, in the run's sample order. `alleles()` is the distinct alleles, the
  reference first; `index_of(bases)` is the lookup B3 will attribute support through;
  `reference()` hands back the gathered reference so B3 need not re-gather it.
- **`AlleleIndex`** — the alleles and the byte-keyed lookup that unifies them, as one type,
  so that pushing an allele and recording where it went cannot become two steps a caller can
  do one of. Keyed rather than scanned because at the top of the committed cohort range a
  locus can hold one distinct allele per sample, where a scan is quadratic in that count;
  production keys the same table the same way. **Nothing iterates the map** — the table's
  order is the `alleles` vector's own — so the hasher cannot reach the output.
- **`alleles_of_sample`** — the rule of §2, handing each allele to a callback. The
  single-record branch projects; the multi-record branch groups the sample's `(read, record,
  sequence)` triples by chain id and composes.
- **`ReadAlleleScratch`** — one sorted `Vec` of those triples plus the compose buffer,
  refilled per sample rather than allocated per sample. **A sorted list rather than a map of
  lists**, which is the shape the pileup walk settled on for the same grouping-by-chain-id
  question and for the same reason (`resolve_mate_overlap_at_pos`): a map allocates a vector
  per read per locus.
- **`MemberPlacement::compose_into`** — one substitution: the reference from where the
  writing has reached up to this member's own region, then the sequence's bases in place of
  that region, answering how far the writing now reaches. `project_into` is now that plus
  the reference to the end of the locus, so the two spellings of the same arithmetic cannot
  drift apart, and the `Complete`-witness assertion lives in one place instead of two.

## 4. Judgement calls, recorded

Five choices the design documents do not make. None changes the ruling; each is at the code
with its reasoning.

- **A read that showed two different things at one record is removed** where the locus spans
  several of the sample's records, like a read that was not there: there is no one thing it
  showed, so nothing can be composed. It is a fragment whose mates overlap and disagree —
  the mint keeps both mates as observations under the one chain id it gave the fragment
  (`resolve_mate_overlap_at_pos`), which is a sequencing error rather than a haplotype.
- **At a sole record the same fragment keeps both**, and this is the one shape on which the
  two branches differ. It is not an optimisation of the other: with one record there is
  nothing to compose across and each mate's sequence already spans the locus, so both are
  complete evidence; composing across records needs the fragment as a unit, and an ambiguous
  one has to go. The review found the first draft claiming the branches agreed, and the
  claim was false on exactly this shape.
- **Two of one sample's records overlapping on the reference are refused** — structurally,
  once per sample, before any read is consulted. The generic mint cannot produce them: an
  event that falls inside an open record's footprint folds into that record and widens it
  rather than opening a second one (`find_overlapping`,
  `locus_generation/pileup/open_record.rs:2398`). The check has to be structural because the
  composition's own backstop is reached only by a read named at *both* records: where they
  carry different reads, every read fails the presence test and the sample's whole evidence
  disappears with nothing to say why.
- **An observation carrying reads and no chain id is refused where a locus spans several of
  a sample's records.** Dropping such reads instead would be silent — the locus would come
  back with fewer alleles and nothing to say why — and it is precisely the state the ruling
  says must never be reached. **This is a producer's guarantee, not a caller's**, so it owes
  the same follow-up B1's reference-width assertion owes: when observations are decoded from
  a psp file it becomes corrupt input and must become a `RunError` (arch §5). Recorded at
  the code.
- **A fragment whose two mates flank an unsequenced gap inside one locus is treated as
  having covered it.** Presence at every record is what the ruling makes decidable, and for
  a single read presence at two records implies coverage of the ground between them; for a
  read *pair* collapsed onto one id it does not, since the insert between the mates was
  never read. The locus would have to be narrow enough — 50 bases at the default bound —
  for both mates to have records inside it, which needs mates nearly adjacent. Nothing in
  the observations can distinguish it today; it is noted here rather than guessed at.

## 5. What the tests pin

27 tests added to `build.rs` (21 → 48), taking the module from 59 to 86 — 21 written with
the step and 6 more from [the review](../reviews/ng_cohort_merge_b2_2026-08-17.md), which is
where the two blind spots were found: no cross-record fixture carried an indel, and no locus
had two samples each holding several records. The ones that carry the step:

- **Unification.** Two samples showing the same change make one allele. One deletion at
  three placements inside one locus unifies — the plan's B2 fixture, with its explanation
  corrected: **projection is what normalises placement**, since dropping any one `C` of a
  `CCC` run gives the same string over the locus. What left-alignment upstream buys is one
  step earlier, and a second test shows it: two placements far enough apart close as *two*
  loci, each table holding one sample's half of the evidence, which unification can never
  repair.
- **The ruling's two branches, on one fixture.** A sample with records at 12 and 14 inside a
  deletion's locus; read 7 named at both composes `ACATC`, read 11 named at both (agreeing
  with the reference at 12) composes `ACGTC`, and read 9 — named at 12 only — is removed, so
  `ACATA` is **absent**. That absence is the assertion that separates the ruling from the
  superseded proposal, which would have credited read 9 with the reference at 14.
- **What B0 bought, asserted as an allele rather than as a presence.** `ACGTC` exists only
  because the record at 12 names the reads behind its *reference* sequence.
- **The other two removal reasons take the same branch**: a read a depth cap dropped at one
  record, and a read that was `Partial` there.
- **A chain id means nothing across samples** — a read 7 in another sample does not complete
  this sample's read 7.
- **The STR path's shape**: a single-record sample's alleles come out with every
  `chain_ids` empty.
- **Determinism**: the same fixture with the ids listed in the other order in every
  observation gives an identical table.
- **An allele is closed on the reference it consumed, not on its own length** — an insertion
  composing six bases over a five-base locus and a deletion composing four. Every other
  cross-record fixture is substitutions, where the two counts are equal at every step, which
  is what let the review's mutation of it pass 80 tests.
- **Two samples that each hold several records are derived independently, sharing a read
  id** — which is the ordinary case, the id space being per file. It is what makes the
  working buffer's per-sample reset observable at all.
- **The reads removed as evidence are counted** (`reads_removed_as_evidence`), one in the
  two-record fixture and none in a sample with a single record.
- **The refusals**: an observation with reads and no id across records; two of one sample's
  records overlapping, with a shared read and without; a sample's records supplied out of
  coordinate order.

## 6. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. (`--all-targets` is red on
  this branch and was red identically before this work: 49 pre-existing errors in
  `examples/`, `benches/` and other modules' test code.)
- `cargo test --lib ng::run::cohort_merge` — 86 passed, 0 failed.
- `cargo test --lib` — **3,709 passed, 0 failed, 11 ignored**. The parent commit's tree
  holds 3,693 tests in all (`cargo test --lib -- --list`, counted there), so this step adds
  27 and `build.rs` is the only file it changes.

## 7. What B2 does not do

- **No support.** Which reads back which allele, and with what moments, is B3's
  (`CohortObservation` / `SampleSupport`). `index_of` exists so that B3 can attribute by
  composing again rather than this step carrying an assignment for every observation.
- **The per-read moments do not split exactly**, and that is the question B3 opens: chain
  ids are per read, but `q_sum`, `mapq_sum`, `mapq_sum_sq`, `num_fwd` and `placed_left` are
  summed per observation. Where an observation's reads split across two cohort alleles the
  split is exact in read count and approximate in the moments. Production faces the same
  thing and subtracts constituent scalars (`per_group_merger.rs`,
  `project_compound_scalars`). Taken to the owner before B3 is coded.
- **The reads removed as evidence are not counted anywhere yet.** They are lost depth, and
  where that count belongs — beside `reads_without_observation` in `SampleSupport`, or
  nowhere — is B3's to settle.
</content>
</invoke>
