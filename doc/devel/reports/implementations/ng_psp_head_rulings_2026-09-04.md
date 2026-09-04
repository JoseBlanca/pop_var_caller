# Three of the owner's rulings on the record head and the read-group table

**Date:** 2026-09-04
**Branch:** `ng-psp-mode`
**Follows:** Milestone H of [psp_head_compared_reads.md](../../ng/impl_plan/psp_head_compared_reads.md), at its checkpoint

The owner reviewed Checkpoint H and ruled on three things. Two reverse decisions this branch had
already committed; the third is new work outside the plan, taken because he asked for it.

## 1. The locus kind comes back out of the head

**Ruled:** *"do we need that number in the front matter? I think it should be inside the records."*

**And the case for it is weaker than the spec argued — the kind is a function of the coordinate.**
A locus's kind is the kind of the *typed region* it falls in, and typed regions are computed from
the reference plus the repeat catalog before any read is looked at (`RegionKind::Generic`,
`SsrSegment`, `SsrBundle`, `Satellite`, mapping one-for-one onto `LocusKind` and to no locus at
all). Every psp records the segmentation inputs its typing used, and `PspVariantCaller::open`
already refuses a cohort whose file disagrees with the run's own segmentation, naming the sample
and the field. So a reader holding a coordinate can look the kind up. **That covers the two
decisions §3.1 cited** — the width bound, which governs generic loci only, and the never-mix
assertion — for the cheap head-only reader the successor plan wants as much as for a full one.

The tag is back at the end of the body, where it was. It stays *in the record* rather than
leaving it, because the body decoder needs it to know whether the tract's motif and flanks follow;
taking it out would make decoding a body depend on having the catalog to hand.

`RecordHead` keeps `reads_compared_with_reference` and loses `locus_kind`; `LocusKindTag` is gone
and `put_kind` / `read_locus_kind` are single functions again. Head fields 7 → 6, body 21 → 22,
declared total unchanged at 28.

**And the forward-compatibility consequence I raised at the checkpoint dissolves with it**: an
unknown kind tag is met in the body again, so a walk that declines a record carries on past it, and
`a_declined_records_body_is_never_decoded` rests on that once more. The owner's point stands
either way — a file from a later writer is the format version's problem, not this field's.

### What it cost, re-measured

Re-taken on the same two corpora with the kind in the body throughout, so the compared-read count
is measured alone:

| | tomato SRR7279481, 10.25 reads a position | HG002 chr21, 280.32 |
|---|---:|---:|
| head bytes a record, before the denominator | 6.195 | 8.661 |
| head bytes a record, now | 7.195 | 10.588 |
| compressed bytes a record, before | 4.944 | 15.831 |
| compressed bytes a record, now | 5.135 | 17.176 |
| **the field's share of the file** | **+3.86 %** | **+8.50 %** |

Against +3.91 % and +8.51 % when the kind was in the head as well — **so the kind's move was free
in both directions, which is what "a move, not an addition" predicted.** Every quoted figure in the
specs, the module docs and the bench moved to these.

## 2. One sample may not declare an `@RG ID` twice

**Ruled:** *"in practice it is more likely that it is a mistake, so that should be a hard error
reported to the user when the bam/cram files are opened."*

`build_read_groups` — the one door every entry point opens alignment files through — now refuses
it, before anything is walked, naming the sample, both files and the fix (`samtools addreplacerg`).
And `OpenPspCohort::open` refuses a stored file whose sample names one id twice, so a stored cohort
is held to what a walked cohort is held to. That reinstates spec §6.2's second clause, which this
branch had recommended dropping at Checkpoint E; the recommendation is withdrawn.

**⚠ Widened to the whole run the same day, on the owner's word — *"yes, we need to do that."***
It was first built scoped to one sample, because his two messages pointed at different scopes and
the narrower one was unambiguous. The wider rule is now what ships: no two read groups anywhere in
a run may share an id, whether they are one sample's or two, and a cohort of stored files is
refused for it exactly as a cohort of alignment files is.

**Across samples nothing merges, and it is refused anyway.** A lane's identity is the integer this
caller mints, so two same-named lanes of different plants would stay apart on their own. What the
rule buys is provenance: every report, every parameters file and every error message names a lane
by its id. Within a sample the collision is worse than untraceable — one lane's reads counted as
another's library, silently.

What was measured while the narrower rule was in place, and what it cost to widen:

- **A collision across samples merges nothing.** A read group's identity in a run is its sample
  together with its id, and the table keys on position, so two individuals whose files both say
  `ID:1` stay apart on their own.
- **It is what a cohort looks like when one pipeline aligned every sample.** Measured on this
  repository's own fixtures: the run-wide rule failed **88** tests, of which **76** were different
  samples sharing an id incidentally. The remaining 12 were the case the ruling is about — one
  sample across two files — and every one of those fixtures has been given distinct ids, which is
  what real per-lane files carry.
- On the tomato benchmark cohort the wider rule would never fire: all six accessions name their
  read group by run accession (`SRR7279481`, `SRR7279488`, …).

**The widening cost 69 fixtures**, every one of them a cohort whose samples all named their read
group `rg1`. The per-sample builders derive it from the plant now, `rg-{sample}`, and the psp
fixtures from theirs — which is what real data carries: the tomato benchmark cohort's six
accessions name their read groups by run accession, so the rule would never have fired on it.

**Collateral, and it is worth naming.** `reject_colliding_synthesized_libraries` is now unreachable
through `build_read_groups`: a synthesized library name is `sample:id:file-stem`, so two collide
only when all three agree, and the id check refuses that input first. The function stays — it is
the only thing between a synthesized name and a silent merge, and the id rule is a policy that
could be relaxed where this one could not — and it is now tested directly rather than through a
door that no longer opens.

## 3. The compared-read count is what a ratio may be taken of

**Ruled:** *"be careful because the reads excluded by a filter or by a cap shouldn't be in that
reads-compared-with-reference."*

**Checked, and it is already right — for a reason worth writing down rather than trusting.** The
field sums `num_obs` over the record's `Complete`-witness observations, and:

- a read a filter turned away never reaches the locus generator, so it is in no observation;
- a read the **depth cap** discarded is dropped *during the walk, before any record exists*
  (`pileup/parity.rs`: "the cap acts in the walk, before any record exists"), so it is in no
  observation either — the record counts those separately in `reads_discarded_by_cap`;
- a read that produced no observation is `reads_without_observation`, counted separately too;
- a read whose witness stopped inside the locus is in neither half, which is what keeps the
  numerator and the denominator over one subset.

The field's documentation now lists all four exclusions, and
`reads_a_filter_or_the_cap_turned_away_are_not_in_the_compared_count` pins the part a fixture can
reach: a record carrying 305 capped reads and 12 unobserved ones reports the same denominator as
the same record with both counters at zero.

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,162 passed**, 0 failed, 14 ignored (6,157 before these rulings).
- `cargo test --tests` — every integration target green.
- **Direct mode's VCF is unchanged**: six tomato accessions over 200 kb of SL4.0, **598 records,
  sha256 `5f0903cf…`** — the hash on record since Milestone D2.

## 4. The parameters file identifies by name, not by position — done

**Ruled:** *"The parameter file should not depend on the order to idenfity the RG or the samples,
it should use the sample name and RG id."*

It did depend on order, in two places at once: samples were compared name-against-name **at the
same position**, and every per-sample and per-read-group table was projected into the *file's* own
order. Both fail for one reason — a run's sample order and its read-group numbering are assigned
first-seen over the alignment files it was given, so the same cohort handed over in another order
lists its samples differently and numbers every lane differently, and a file fitted on that cohort
is still the right file.

Samples are matched as sets by name; read groups by **the sample and the `@RG ID` together** — the
pair and not the id alone, since an id is unique within a sample and nothing makes it unique
across them. A missing sample or lane is a hard fail naming it, in both directions. A matched
lane whose library differs is refused naming both libraries. `WhereTheFilesRowsBelong` says where
each of the file's rows goes: the file's own order where there is no run, the run's otherwise.

**The proof is where the values land, not that the file is accepted.** Accepting a re-ordered file
while still projecting into the file's order would be the silent mis-pairing the old refusal
existed to prevent — every plant handed its neighbour's inbreeding coefficient — and would pass a
test that only checked for `Ok`. Measured: replacing the run-order mapping with the identity leaves
241 tests green and fails exactly the one that asserts each plant keeps its own coefficient and
each lane its own evidence.

Two tests asserted the overturned behaviour and are inverted rather than deleted, each carrying why
it changed. One splits in two, because with the join on sample and id a row differing in either is
a missing lane rather than a per-field mismatch.

## Nothing further owed from the rulings

All four are built. What is next is the plan's own next milestone — the `call-from-psps`
subcommand — and it carries one open question of its own, recorded at Checkpoint E and unchanged:
**what a run over stored files says about each sample in its report.** Direct mode reports
per-sample walk tallies — regions handled, generator counts, read-filter drops — and a psp source
has none of them, because the walk happened in another process and what it counted is in that
file's provenance rather than its records. The caller returns the calling tallies alone. That is
the owner's to settle, not the implementer's.
