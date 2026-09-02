# ng STR observations — B3: routing parity, and the no-regression pin

*2026-09-02. Step B3 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §10](../../ng/spec/run_ssr_observations.md). Branch `ng-ssr-observations`.*

## Plan

B1 changed which ground a run sends down the repeat path. Two claims are owed: that the run's
partition is the one the measuring tool prints, and that where the routing did not move,
nothing else did either.

## Assumptions

**The parity test compares the run against the dump's own two lines rather than against the
dump's text.** `examples/ng_typed_region_dump.rs` prints a formatted TSV; a test that shelled
out to it and parsed the output would be testing the formatting. What matters is that the two
ask the catalog the same question, which is the call the test makes verbatim.

**The parity test also asserts the partition is exact**, which the dump cannot: it prints
whatever it is handed. Every base of the fixture contig is in exactly one region, in ascending
order, with no gap and no overlap — a property, not a comparison.

**The no-regression fixture puts every read on ground that is generic under both settings.** A
read over the six-base homopolymer would make the two runs *legitimately* differ — one calls
that ground, the other sets it aside — so the pin would be asserting something false. The
claim being pinned is narrower and is the one B1 makes: the change is which ground is generic,
never what happens on it.

## Changes made

Tests only; no source changed.

| test | what it holds |
|---|---|
| `the_runs_ground_partition_is_the_dumps_at_the_same_floors` | the run and the tool agree, at both sets of floors; and the partition is exact |
| `where_the_routing_did_not_move_the_vcf_is_byte_identical` | two settings of the routing switch, one cohort, identical VCFs |

Two fixture helpers came with them: `the_fixture_references_bases`, lifted out of B1's
reference builder so reads can be written against the same bases; and `a_read_showing`, because
the shared `read_named_with_length` writes a run of `A`, which on this reference is thirty
mismatches rather than one variant.

## What the fixtures had to get right, measured

**The reads are 30 bases because 30 is the shortest the filters keep**
(`DEFAULT_MIN_READ_LENGTH`). At 12 bases the run's own report said so — *"8 reads kept: 0, 8
dropped by the read filters — library rg1 of one: 8 too short"* — and the VCF was header-only.
A byte-identity test over two empty files would have passed and held nothing, which is why the
test asserts a record count above zero before it compares.

**Both runs write `calls.vcf`, in two directories.** The VCF header names the parameters file
beside it, so two runs writing `calling.vcf` and `catalog.vcf` differ on that line and on
nothing else — a difference in what the test asked for rather than in what it is testing. That
was the first failure the pin produced, and it was the pin working.

**The two runs really do route differently, and the test says by how much before it compares.**
At ng's floors the contig is 3 typed regions and 126 of 136 bases are called; at the catalog's
it is 5 regions and 120 bases. The run report shows both. Without that assertion the test could
pass by comparing two identical runs.

Each run writes 2 records over 8 reads.

## Validation

In the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --tests --examples --all-features --no-fail-fast` — **5,928 passed, 0
  failed, 14 ignored** in the library suite (5,926 before this step).
- **Mutation-tested**: restoring `routing_criteria` to `StrRepeatCriteria::default()` — the
  behaviour B1 removed — fails **both** of this step's tests.
- Unchanged by this step and still red: the three locus-dump tests and the psp writer bench,
  recorded in `PROJECT_STATUS.md`.

## Tradeoffs and follow-ups

- **The parity is against the tool's code path, not against an independent implementation of
  the classification.** `partition_resident_in` is that independent implementation — it finds
  repeats from the bases where the catalog reads them from a file — and it is already the
  catalog's own differential oracle. Pointing the run at it too would be a stronger claim than
  this one and is not what spec §10 asks for.
- **The no-regression pin is on a synthetic 136-base contig with one sample.** It holds the
  property; it says nothing about how much ground moves on a real genome. That is B4.
