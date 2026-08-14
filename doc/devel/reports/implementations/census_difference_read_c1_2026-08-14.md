# C1 — a tract's mismatching bases are numbered by the locus's reads, not by one observation's

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step C1.
**Source:** the milestone-A review,
[census_rename_a2_2026-08-14.md](../reviews/census_rename_a2_2026-08-14.md) finding B1.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§3 and §7.3; [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.4.
**Date:** 2026-08-14.

---

## 1. What was wrong, and what it destroyed

At a repeat tract the census keeps every mismatching base as *(which read, how far into the
tract, which base)*. The read was numbered **from zero within one observation** — and the walk
folds reads carrying identical bases into a single observation with a count. So two reads that
differ from the reference **in different places** each arrived as their own observation of one
read, and both were written down as read 0.

**That is the one distinction the field exists to make.** The same substitution at the same
place on two reads is an allele — an interrupted repeat the population carries. Two
substitutions on one read is a bad read. A consumer grouping the entries by read saw the second
where the truth was the first.

The existing fixture handed the writer **one** observation with a count of two, which is
precisely the regime where per-observation and per-locus numbering agree.

## 2. The fix

The writer carries a per-(locus × read group) counter across the observation loop, and a
difference's read is that counter plus its offset within the observation.

**The counter advances over the reads that were compared** — the unslipped ones whose tract is
the reference's length, which are the only reads a mismatch can be read off at all — and not
over slipped, guarded or partially-covering reads. Every number spent on a read that cannot
carry a difference is headroom lost from a byte, which §4 is about.

## 3. Tests, and both fail against the code as it stood

| test | what it pins |
|---|---|
| `two_reads_carrying_different_interruptions_are_two_reads_and_not_one` | the review's own test: two reads, one interruption each at different offsets, come back as two reads |
| `a_locus_numbers_its_reads_from_zero_once_and_not_once_an_observation` | three observations at one tract — two reads sharing an interruption, one with its own, four clean — so the shared one is reads 0 and 1 and the lone one is read 2, not read 0 again |

**Measured, not assumed:** with the numbering put back the way it was, those two fail and the
other 32 census tests pass — `32 passed; 2 failed`. That is the whole of what the old suite
could see.

## 4. What the fix does not fix, and it is the deep end of the range

**A byte holds 255 and the per-locus read cap is 1,000.** Now that reads are numbered within the
locus rather than within an observation, a tract with more than 255 compared reads in one read
group writes every read from the 256th onwards as 255 — several reads reported as one, which is
the confusion this field exists to prevent, arriving at the deep end of the range this caller
commits to (`CLAUDE.md`: a few reads a position to several hundred).

It is reachable on the data already in use: the GIAB trio crosses 171,930 reads over 216 tracts
in three samples, about **265 crossing reads a tract a sample**.

**Measured on both oracles, after the fix.** The harness now reports the highest read a
difference sits on and how many sit at the byte's ceiling:

| | difference entries | highest read | at the ceiling |
|---|---:|---:|---:|
| HG002 | 583 | 255 | **75** |
| HG003 | 471 | 255 | **96** |
| HG004 | 1,129 | 255 | **294** |
| tomato, eight accessions | 16 to 303 | 3 to 35 | **0** |

**About one difference entry in five is at the ceiling on the trio** — 465 of 2,183 — and each
of those is a read the census can no longer tell from the other reads there. At 2.4 to 30.6
reads a position nothing comes near it: the highest read numbered is 35.

**This is not a regression the fix introduces so much as one it exposes.** The old numbering
saturated too — at 255 reads carrying *identical* bases in one observation — but it wasted the
byte differently, spending it again at every observation. What the fix does is make the
numbering mean what the field says it means, and at that point the byte is visibly too narrow
at the deep end.

### 4.1 Settled the same day: the field is two bytes — DECIDED 2026-08-14 (owner)

**`TractDifference::read` is a `u16`.** It reaches 65,535 where a locus enters at most 1,000
reads, so the numbering cannot collapse at any depth this record is written under.

**Measured on the trio after the widening**: the highest read a mismatch sits on is **350, 367
and 387** in the three samples, and **nothing is at the field's ceiling** — those are the reads
that were being written as 255. On tomato the highest is unchanged at 3 to 35.

**What it costs, measured on both oracles.** An entry goes from 8 bytes to 12, so a sample's
whole census grows by about a percent at the deep end and by nothing at the shallow one:

| | difference list | the whole census for that sample |
|---|---:|---:|
| HG004, 1,129 entries | 0.009 → **0.014 MB** | 0.432 → **0.437 MB** |
| HG002, 583 entries | 0.005 → **0.007 MB** | 0.286 → **0.288 MB** |
| tomato, deepest of eight | 0.002 → **0.003 MB** | 0.108 → **0.109 MB** |

**No fitted number moved** on either oracle; the only lines that differ are the sizes.

**One consequence for a document this plan may not edit.** Spec §6's measured difference-list
sizes — 0.2 MB at 2.4 reads a position and 2.2 MB at 30× across tomato's 462,701 tracts — were
taken at 8 bytes an entry. At 12 they become **0.3 MB and 3.3 MB a read group**. The arch
document's §1.4 and §3 have been corrected to say `u16` and why; the spec's two figures are for
the owner to fold in.

A test pins the property rather than the width: `three_hundred_reads_at_one_tract_are_three_hundred_reads`
gives one tract three hundred reads each carrying an interruption and asserts they come back
numbered 0 to 299.

## 5. Nothing reads the difference list yet, so no fitted number can move

The only consumers of `SsrEvidence::differences()` outside the census's own tests are the
harness's size report and the new diagnostic in §4. The substitution channel is fitted per
stratum from `bases_compared` and nothing else; the interrupted-repeat model that would read
this list is not built. So this step changes what the census records and moves nothing the fit
returns — which is why the plan puts it after milestone B's oracle rather than inside it.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::joint::census` | `35 passed; 0 failed` (32 before) |
| `cargo test --lib` | `3,584 passed; 0 failed; 11 ignored` (3,581 before) |
| the 88-second tomato oracle | no line differs but the new diagnostic |
| the 74-second trio oracle | no line differs but the new diagnostic |
