# A1 — what a piece of the census is, and where it sits

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step A1 — implementation plan 2.
**Design authority:** [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.1a; [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§6.2.
**Date:** 2026-08-14.

---

## 1. What this step is for

**The fit never needs a whole sample's evidence at once, and this is the type that says so.** It
finishes the ordinary positions before it reads a tract — the tract half takes one number per
sample from the generic half and hands nothing back — and within the tract half it fits one band
of strata at a time. So the smallest piece anything ever asks for is *one read group's ordinary
positions* or *one read group's tracts for one stratum*, and that is what `SectionKey` names.

**Types only. Nothing holds a section yet**, which is the next step: `SampleCensusEvidence` still
carries its two public maps, and no code path constructs a `Section` outside this module's tests.

## 2. What landed

| type | what it is |
|---|---|
| `SectionKey` | `Generic(ReadGroupId)` or `Ssr(ReadGroupId, Stratum)` — and its `Ord` **is** the enumeration order |
| `Section` | the decoded contents, `Generic(GenericEvidence)` or `Ssr(SsrEvidence)`, with `answers(key)` |
| `ByteExtent` | an offset and a length, with `overlaps` |
| `Sections` | where a sample's sections are; today one state, `Resident`, built by a walk |

**`Sections::resident` refuses a section filed under a key of the other kind.** An
ordinary-position key holding tracts would index one half of the census by the other's
positions, and nothing downstream would notice.

**`Sections` is public and that is not what keeps a section from being retained.** What will do
that is the scoped access on the value that owns one (step A2): a section is lent for the length
of a call, and there is no field a caller could keep it in.

**Two names deliberately not introduced.** `Sections::Backed` and the reader it holds belong
with the file (milestone B); a second variant with nothing to read would be a state no code can
reach.

## 3. Tests

| test | what it pins |
|---|---|
| `sections_enumerate_ordinary_positions_first_then_by_read_group_and_stratum` | the enumeration order, including that period 1 at 8 repeats precedes period 2 at 6 — **the mutation it catches is the two variants swapping places**, which would read every sample's tracts before its ordinary positions |
| `a_section_key_says_whose_reads_and_which_stratum` | the two accessors, and that ordinary positions are in no stratum |
| `a_section_filed_under_the_wrong_kind_of_key_is_refused` | the constructor's panic |
| `a_resident_sample_lists_its_sections_and_hands_back_the_one_asked_for` | listing, lookup, and that a section never asked for is **absent** rather than empty |
| `byte_extents_meet_without_overlapping` | abutting is not overlapping; a zero-length extent overlaps nothing |

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors, and no warning in `census.rs` |
| `cargo test --lib ng::parameter_estimation::joint::census` | `40 passed; 0 failed` (35 before) |
