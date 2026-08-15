# B1 — the census on disk: a header, a directory, and the sections

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step B1 — implementation plan 2.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§6.1, §6.2, §7.1; [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§2.2, §2.3.
**Date:** 2026-08-14.

---

## 1. What landed

A new module, `joint/census_file.rs`: `write_census`, `read_census`/`decode_census`, and
`decode_directory_of` — which says where every section sits without decoding one, and is what the
seeking reader of B2 will open with.

```text
  magic     8 bytes   "NGCENSUS"
  version   u16       1
  header              the sample's name, the twelve terms, and the pileup it was built from
  directory u32 n, then n × (section key, offset u64, length u64)
  sections            the bytes each directory entry points at
```

Every integer is little-endian and every offset is absolute, so a reader seeks once a section.
The writer encodes the sections first, then lays the directory out over them — the directory's own
size depends only on its keys, so it can be sized before the offsets are known and filled in
afterwards.

**The header is outside the sections because it is read before them.** The twelve recording terms
and the kept-loci digest are compared across every sample before anything large is decoded (spec
§5), so a cohort of a thousand that cannot be pooled says so after a few hundred bytes each.

**A separate module, not `census.rs`.** The architecture puts the types in `census.rs`, and they
stay there; this is the framing that puts them in a file, and `census.rs` is already 3,000 lines
about what the evidence *is*. `CensusError` stays in `census.rs`, where arch §2.3 puts it.

## 2. The one design question this step had to answer

**The seven selection values now travel as a table of `(field name, digest)` rather than as
values** — `SelectionTermsDigest`, and `RecordingTerms::selection` holds one.

The architecture asks for exactly this and gives the reason (§2.2): `SelectionTerms` holds
`StrRepeatCriteria` and `ScanParams` — two other modules' configuration, twenty-odd scalars
between them — and a codec over the values would have to track every field either type ever
grows. A field it failed to track would drop out of the comparison silently, which is the one
failure the whole check exists to prevent.

**What it costs, and it is deliberate: a census file can say whether it matches, not what it was
built at.** Nothing reads those seven as values — the only two uses in the tree are
`first_disagreement` and one test — so nothing is lost to a program; a person wanting a run's seed
reads it from the run.

**How a field is digested, and which way it fails.** The two configuration values are hashed
through their `Debug` rendering, so a field added anywhere inside them changes the digest without
anyone remembering to come back. The cost is that a changed `Debug` impl reads as changed
settings — which refuses a pooling rather than allowing a wrong one, and that is the direction to
fail in.

**The names are one list** (`SELECTION_FIELDS`), used by the writer, the reader and the refusal,
so a mismatch still names the same field it always did — and a file whose table is not this
build's own is refused rather than compared entry by entry against fields that mean something
else.

## 3. What the header carries about the pileup

`PileupIdentity` — a digest of the pileup's header and its record count, **never a modification
time** (spec §6.1). B1 stores it and B3 makes it act. `None` is a distinct state from a
zero-record identity, and §4's test says so.

## 4. Tests

Eight, in the new module.

| test | what it pins |
|---|---|
| `every_corner_state_survives_a_round_trip` | **the assertion the step rests on.** One sample carrying every corner spec §7.1 lists: a never-walked position, a walked-at-zero-depth one, one with reads and nothing non-reference, one at the depth cap with a 255-read allele, a multi-allelic one; and two strata of tracts with both end buckets saturated, a guard entry, a difference at read 299, and a locus the walk never reached |
| `a_stream_that_is_not_this_builds_census_is_refused` | empty input, a wrong magic, a wrong version, **and every one of the fixture's 775 truncations** — a length read off the end of the buffer has to be a refusal wherever it falls, not only at the top of the file |
| `the_directory_places_every_section_end_to_end_and_none_overlap` | the enumeration order is the key's own; sections abut; the last one ends at the end of the file |
| `writing_the_same_census_twice_gives_the_same_bytes` | what §7.12's byte-for-byte comparison between the two builders will rest on |
| `the_three_states_at_an_ordinary_position_come_back_apart` | the three states read off the *file*, not off the value that was written |
| `a_tracts_offsets_guard_and_difference_come_back_unchanged` | the saturating ends, the two per-stratum counts, the difference's read number and signed offset |
| `the_pileup_a_census_was_built_from_survives_and_its_absence_is_not_a_zero` | the two states of the identity |
| `a_census_with_one_half_empty_round_trips` | a census with no tracts, and one with no ordinary positions |

**Two mutations run by hand, both of the kind this step's failure mode is** — a plausible number
rather than a crash:

| mutation | result |
|---|---|
| the sparse entry's allele byte and count byte transposed on read | `4 passed; 4 failed` |
| a tract's censored-read count moved one field later in the section | `4 passed; 4 failed` |

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::parameter_estimation::joint::census_file` | `8 passed; 0 failed` |
| `cargo test --lib` | `3,600 passed; 0 failed; 11 ignored` (3,592 before) |
| the 88-second tomato oracle | byte-identical to A4's |
| the 74-second trio oracle | byte-identical to A4's |

**Both oracles are unchanged, and that is the check that matters for the terms' change of
shape**: the writer now digests the selection instead of carrying it, and the recording-terms
report both cohorts print is the same one it was.

**What is not measured yet: what a census file weighs on real reads.** Nothing writes one during
a walk — B4 is where the same cohort is fitted from memory and from files, and that is where the
size is read off a real sample rather than computed from the layout.
