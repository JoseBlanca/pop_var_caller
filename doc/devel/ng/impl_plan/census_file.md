# ng census — implementation plan 2: the file, and reading it one section at a time

**Status:** draft, 2026-08-14. The build order for the one genuinely new artefact in this route: a
census file per sample beside its pileup, laid out so that the parameters fit can read one part of it
without decoding the rest. Design is settled in
[`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §6.1 and §6.2 and in
[`arch/parameter_prepass_joint_records.md`](../arch/parameter_prepass_joint_records.md) §1.1, §1.1a and
§2.2. **This plan turns that design into build order; it is not a place for new design.**

**It follows [plan 1](census_rename_and_encoding.md)**, which renames the types this plan serialises
and settles the encoding it writes. Starting this first would mean writing a format for types about to
be renamed and codes about to change meaning.

**Nothing in this plan is needed to compare the two routes.** A run that goes from alignments straight
to a fit holds everything in memory and always could; this plan exists for the run that walks each
sample once, months apart, and fits later — and for the cohort too large to hold.

---

## Scope

**In:**

- the byte layout of a census file, its directory, and the identity of the pileup it was built from;
- `Sections`, `SectionKey` and the scoped access that lends a section for the length of a call;
- `CohortCensusEvidence`, and `fit_jointly` taking it instead of a slice of whole record sets;
- the second producer — building a census from an existing pileup rather than from alignments;
- `LocusEvidence`, the per-locus item the walk emits and the fit iterates.

**Out (later plans, or not plans at all):**

- **how many samples the fit needs resident, and how it runs inside a memory ceiling** —
  [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11 questions 8 and 10.
  Both are measurements, and both are about the object this plan builds, so they come after it.
- **the duplicated-copy class in its coupled form**, and **the fourth shape** the benchmark trio
  showed — design work, owned by the fit specification.

---

## Principles (how the order was chosen)

- **The algorithmic heart before the plumbing.** The scoped access is the design's whole point — a
  section is lent and cannot be retained — so it is built and tested against a resident value, where
  no file exists to confuse a failure, before any bytes are written.
- **Simplest implementation first, as the oracle for the next.** A resident `SampleCensusEvidence` is
  the parity oracle for a file-backed one: same calls, same values, and §7.12's byte-for-byte
  comparison is what says the second is the first.
- **Types first, then implementation** (project rule).
- **Isolate a step whose failure is silent.** A codec that reads a field at the wrong offset produces
  a plausible number, not a crash. Every codec step lands with a round-trip over every corner state
  the records spec §7.1 lists.
- **Verify against ground truth.** The oracle is the direct run: the same cohort, fitted from memory
  and fitted from files, must give the same parameters.

---

## Preconditions (already in place)

- [Plan 1](census_rename_and_encoding.md) complete through checkpoint B — the types carry their final
  names and the encoding is settled.
- The pileup's header carries its reference, its analysed regions and its read filters, so a census
  file can name what it was built from.
- `examples/ng_joint_records_walk.rs` fits a cohort held in memory — the parity oracle for every step
  below.

---

## The steps

### Milestone A — sections, with nothing on disk

✅ **A1 — `SectionKey`, `Section`, `ByteExtent`, and `Sections::Resident`.** The types, and the
resident half only. `SectionKey`'s ordering is the enumeration order the contract promises.
*Depends:* —. *Source:* [arch records](../arch/parameter_prepass_joint_records.md) §1.1a.

✅ **A2 — `SampleCensusEvidence` holds its sections privately, and lends them scoped.** The two public
maps become a private field; access is `with_generic` and `with_strata`, which hand a closure borrows
and take them back. Nothing can retain a section.
*Depends:* A1. *Source:* [arch records](../arch/parameter_prepass_joint_records.md) §1.1, §2.2.

✅ **A3 — `CohortCensusEvidence`, and its scoped calls across every sample.** Built from a vector of
per-sample values, checking the twelve recording terms across all of them **before any section is
decoded**. The unit lent is a band of strata across every sample, because 68 of tomato's 141 strata
borrow from their neighbours.
*Depends:* A2. *Source:* [arch records](../arch/parameter_prepass_joint_records.md) §2.2;
[loci spec](../spec/parameter_prepass_joint_loci.md) §3.6.

✅ **A4 — `fit_jointly` takes `&mut CohortCensusEvidence`.** The estimator stops holding whole record
sets. **The fitted numbers must not move** on either cohort — this is a change of access, not of
arithmetic.
*Depends:* A3. *Source:* [arch fit](../arch/parameter_prepass_joint_fit.md) §2.1.

> **Checkpoint A:** the estimator reads by section, from memory, with the fitted numbers unchanged on
> the tomato cohort and the human trio. Pause for review.

### Milestone B — the file

✅ **B1 — the layout and its directory.** A header carrying the twelve recording terms, the kept-loci
digest and the pileup's identity; then a directory of `SectionKey → ByteExtent`; then the sections.
Terms and digest outside the sections, since they are compared before anything large is decoded.
**Own commit, do not bundle**, with a round-trip over every corner state.
*Depends:* A4. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.2.

✅ **B2 — `Sections::Backed`, and one read per section.** A call seeks once, fills a buffer once and
decodes from the slice; nothing is retained between calls, because there is no field to retain it in.
*Depends:* B1. *Source:* [arch records](../arch/parameter_prepass_joint_records.md) §1.1a.

☐ **B3 — the staleness check.** The census names the pileup by a digest of its header and its record
count, never modification time; a mismatch rebuilds where the pileup is reachable and fails naming the
field where it is not.
*Depends:* B1. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1.

☐ **B4 — assert the milestone.** The same cohort fitted from memory and from files gives the same
parameters, and §7.16's counting reader shows one section's read touches only that section's bytes.
*Depends:* B2, B3. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §7.15, §7.16.

> **Checkpoint B:** a census on disk, read a section at a time, with the fit's answers unchanged.
> Pause for review.

### Milestone C — the second producer

☐ **C1 — build a census from an existing pileup.** The same `CensusWriter`, driven from a pileup's
locus stream rather than from alignments, downstream of every filter and cap the pileup applied.
*Depends:* B4. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1;
[arch records](../arch/parameter_prepass_joint_records.md) §2.2.

☐ **C2 — the byte-for-byte agreement test, and it is the point of the milestone.** One sample's census
built during a walk and again from the pileup that walk wrote must be identical. **This is the test
that says whether the pileup holds everything a census needs**, and it fails on precisely the fields
that do not survive the round trip — so run it on a fixture carrying every corner state, including a
repeat tract, whose per-read length is the field most likely to be missing.
*Depends:* C1. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §7.12.

☐ **C3 — the subcommand.** A census built from an existing pileup, for the three cases that need it:
pileups older than this format, a census lost or built at knobs since changed, and a census wanted
larger than the one on disk.
*Depends:* C2. *Source:* [records spec](../spec/parameter_prepass_joint_records.md) §6.1.

> **Checkpoint C:** two producers, agreeing byte for byte. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the fitted parameters unchanged on both real cohorts, the change being access rather than arithmetic |
| B1 | a round trip over every corner state the records spec §7.1 lists |
| B2 | a counting reader: one section's read touches only that section's bytes (§7.16) |
| B3 | a changed pileup rebuilds; an absent one fails naming the field; a touched modification time changes nothing |
| B | memory-fitted and file-fitted parameters identical on both cohorts |
| C | one sample's census built two ways, byte for byte (§7.12) |

---

## Out of scope (next plans)

- **How many samples the fit needs at once, and how it runs inside a memory ceiling** — measurements
  against the object this plan builds
  ([`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11, questions 8 and 10).
- **The duplicated-copy class in its coupled form**, and **the fourth shape** — design work in the fit
  specification before either becomes build order.
