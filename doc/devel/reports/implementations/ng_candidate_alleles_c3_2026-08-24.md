# ng candidate alleles — C3: the leftover

*2026-08-24. Step C3 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
Milestone C's last step, on top of `1041e30e`. Design authority:
[`../../ng/spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §4.1, §5, §5.1, §8.*

---

## 1. What it builds

Every allele the admission rule or the cap removed keeps its reads' error mass in the arithmetic —
the SNP/indel genotype likelihood carries a term for reads matching no candidate. Nothing upstream
produces it, because nothing upstream drops anything, so selection owes it. `leftover_of` sums it
straight from the merge's own per-row `q_sum`, per sample.

**And the count that decides whether a sample is genotyped at all.**
`earned_reads_cut_by_the_cap` is this sample's reads on an allele **it** earned and the **cap** cut.
Keying it on the pool instead would no-call almost everybody: the rule drops sequencing error at
13,166 of 15,474 alternatives on the GIAB trio at 300×, and every sample carries a few error reads
at nearly every locus.

**Asking that needs no list of what the cap cut** — a correction a C2 reviewer supplied by writing
C3: a sample that cleared the rule for an allele is by construction a sample that made it a
candidate for the cap, so *dropped, and this sample reached the rule* is the whole test.

## 2. Changes made

- `mod.rs`: `leftover_of`; **`one_run_per_allele`**, the fold's own chunk-and-assert loop lifted
  into a shared helper so the two walks over a sample's rows cannot come to pool differently — a C2
  reviewer named that as the one thing the two steps must keep in step, since a sample with two
  libraries could otherwise clear the rule in one walk and not the other; and **spec §8's third
  assertion**, on a non-finite quality mass, which no earlier step could implement because none
  read `q_sum`.
- `generic.rs`: `select_generic` fills the leftover instead of zeroing it, and 15 tests.

## 3. What the review changed

Full account in [`../reviews/ng_candidate_alleles_c3_2026-08-24.md`](../reviews/ng_candidate_alleles_c3_2026-08-24.md).

**Three Blockers, all tests that could not fail.** The earned count could be the running pool total,
because every fixture gave the affected sample exactly one dropped allele. A leftover that skipped
row-less samples and padded the tail kept the right length and slid every later sample's value onto
its neighbour, unseen. And the denominator could be narrowed to the surviving alleles' reads —
systematically smaller than the right one, so it no-calls samples the rule never meant to touch.

**Two Majors.** Spec §8's third assertion was missing, and I had not noticed that C3 is where it
first becomes reachable. And the 400-sample test held *how many* samples the cap costs but not
*which*, at the cohort size where that claim is weakest.

**One refactor accident and two false claims of mine**: the lifted helper landed between
`summarise_alleles`'s doc comment and the function, so the fold lost its documentation and the
helper inherited a `# Panics` list of three it raises one of; the pool oracle claimed to pin an
addition order that its exact-binary-fraction masses cannot; and the helper's `# Panics` was
unconditional where the iterator is lazy.

## 4. Validation

- `cargo fmt --check` clean; both `clippy` gates clean;
- `cargo test --lib` **4,288 passed, 0 failed, 14 ignored** in 42.9 s, against 4,276 at `1041e30e`.

**Ten mutations, ten killed**, four of them survivors before the review's fixes.

## 5. Follow-ups — all at Checkpoint C

The range numbers are the substantive ones and are in §3 of the review report. The headline: **the
cap stops being a safety valve well before 400 samples**, where it binds at essentially every locus,
so spec §4.1's "rare: 23 of 53,935 tomato loci" is a fact about 63 accessions and does not carry.
