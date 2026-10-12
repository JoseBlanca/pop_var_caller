# Review: the allele window reaches past the recorded reads (commit 38b1fe56)

**Date:** 2026-10-12. **Branch:** `stutter-multiplier`. **Scope:** correctness of commit 38b1fe56
against [`str_allele_window.md`](../../ng/spec/str_allele_window.md) and its plan. Read only; nothing
was built or run. Line numbers are at 38b1fe56.

**Verdict.** The window change and calling's edge split are correct. One real defect, in the
counting fit (a prototype), which the wider window now exposes at the default. The rest is stale
documentation.

---

## Findings, most severe first

### 1. Medium — the counting fit loses the direction and the size of slips from alleles that sit on the end bucket

`src/parameter_estimation/joint/ssr_fit/composite_fit.rs:247-266`. The E-step reads a slip's
direction and size off `from = slot(allele)`. That slot is the allele's offset clamped to the
recorded range. At ±6 no allele class reached the end bucket. At ±9, the classes at −8 and −9 both
map to bucket 0, and +8 and +9 both map to bucket 16. Three consequences follow:

- **Every shorter slip from an allele at −8 or −9 lands in bucket 0, which equals `from`.** It is
  counted as slipped, but neither `shorter` nor `longer` is incremented (lines 248-252), and it
  adds no step (line 256). Longer slips from the same alleles are still counted, so the shorter
  share is biased towards longer. The long side is the mirror image.
- **A one-repeat longer slip from −9 lands in bucket 0 too**, so it is dropped from the direction
  count as well.
- **For an allele at −9, `distance` is measured from −8.** Every longer slip is therefore
  counted one step short, which biases the fall-off low.

**Failure scenario.** Take the commit's own drawn class: about 0.15 of alleles at −8, shorter share
0.8. Worked by hand (not run), the counted shorter share comes out near 0.77 against 0.80 drawn.
The test `alleles_past_a_narrow_window_are_counted_as_stutter_and_the_shipped_window_keeps_them`
checks only the level and the far shares, so it passes. The size of the error grows with the share
of alleles at ±8 or further, which is exactly the long homopolymer and dinucleotide classes this
change targets.

**Fix.** Split the slipped mass in a bucket between directions using the model's own terms rather
than the clamped slot. `read_probabilities` already computes the per-direction contributions; it
would have to return them separately. Compute the distance from the unclamped offset
`class − span`.

**Not affected:** the climb, which uses the full likelihood (item 3 below), and the library
multipliers' E-step (`library_stutter.rs:89-135`), which uses only "exact or slipped". That
question is answered correctly with a clamped slot.

### 2. Low — `DEFAULT_MAX_ROUNDS` and `QUADRATURE_POINTS` were tuned at 13 classes

`ssr_fit.rs:1104-1110` sets the limit at 40 rounds. That number was measured at 13 classes, where
walks took 6 to 36 rounds. `QUADRATURE_POINTS = 256` (`ssr_fit.rs:927`) was also measured at 13
classes, and spec §3 already says the error of the average over tract frequencies grows with the
class count. At 19 classes neither has been measured. The climb ran 11 times slower, not the 2.1
times the spec estimated, which suggests it now needs more rounds.

**Failure scenario:** walks at 19 classes run out of rounds at 40 and return unsettled answers.
Nothing in this commit would show that. **Ask:** when plan step 3 times the climb, report how each
walk ended (`ClimbEnding`) as well as its time.

### 3. Point 3 checked — `read_probabilities` handles an allele beyond the read span correctly

`ssr_fit.rs:117-144`. Take an allele at −9 with a read span of 8:

- The unslipped read is clamped into bucket −8.
- Shorter slips have no steps inside the range (`inside_steps = 0`), so the whole shorter tail,
  `fall_off^0 = 1`, goes to bucket −8.
- Longer slips have 17 steps inside the range, running from −8 to +8. The remaining tail,
  `fall_off^17`, saturates into +8.

The total is 1. The climb's per-allele read distributions (`ssr_fit.rs:4008`) and both fits' slot
functions iterate `-allele_span..=allele_span` and clamp to `read_span`, so they are correct. **No
defect.**

### 4. Point 2 checked — the edge split in `fill_candidate_weights` is correct

`src/calling/genotype_prior/seed_ssr.rs:310-343`. Every case asked about behaves as intended:

- **Off-by-one.** `offset_of` clamps to `[-span, span]`, so a candidate exactly at ±span and one
  past it both land on the edge.
- **Distinct lengths.** They are counted by repeat count. Equal counts give equal offsets, so
  removing duplicates by count is exact.
- **Division by zero.** It cannot happen: an edge's divisor is used only by a candidate on that
  edge, which makes the divisor at least 1.
- **Both sides past the window.** The two edges are counted independently.
- **The stated-flat rung.** It takes the `None` arm and is unchanged.
- **Total weight.** An edge's weight, summed over its distinct lengths, is still that class's
  weight, so `fill_ssr_seed`'s conditioning argument holds.
- **`fill_seed_share_per_candidate`.** It shares the same function, so the seed and the share
  cannot disagree.
- **The contamination consumer.** `repeat_tract_parameters.rs:494` normalises the shares and
  scatters them by byte length. It is unaffected beyond the intended change of values.
- **Discovery.** It runs as a pre-pass inside candidate selection, so the seed is built once over
  the final candidate set. No partial refill can leave a stale split.

I checked the test's expected numbers by hand: 0.025 and 0.04 in the new test, and 0.0303 in the
test that passes a wrong reference count. **No defect.**

### 5. Point 1 checked — nothing in `src/` is still sized for 13 classes or 91 genotypes

Every table is sized from `config.allele_span` or from the length of the spectrum:

- `ssr_fit.rs:1315`, `1740`, `2661` and `2882`;
- `library_stutter.rs:55` (asserts that the held spectrum has the configured length);
- `stratum_fits.rs:329`, which reads the span from the share count.

The parameters-file validator accepts any odd count of three or more (`validate.rs:1327-1349`), and
no check requires all strata in a file to have equal lengths. The census width
(`RECORDED_OFFSET_RANGE`) is untouched. The memory figures in `estimate_parameters.rs:135-136` and
`ssr_fit.rs:1008` and `1074` are right: 190 × 8 bytes × rows gives 0.71 and 15.3 GiB, and 96 kB a
tract at 63 samples. **No defect.**

### 6. Low — documentation still describing the old behaviour (point 4)

- **`ssr_fit.rs:55-58`.** The module doc still says "The records store `±4`; the lengths the fit
  may place allele mass on reach `±6`". Spec §3 flagged the same stale ±4 in `ALLELE_SPAN`'s
  comment; this second copy was missed.
- **`seed_ssr.rs:135-142` and `seed_ssr.rs:229-231`.** These say two spellings of one length
  "each receive that class's full weight". On an edge class, each now receives that length's share
  of the edge.
- **`seed_ssr.rs:174-175`.** This says "Three passes over the locus's candidate lengths". There
  are now two more, the edge counts, and each is quadratic in the candidates. It still allocates
  nothing.
- **`seed_ssr.rs:146-148` and `ssr_fit.rs:937-941` (`ALLELE_SPAN`'s comment).** These say alleles
  past ±9 produce the same reads as one at ±9. That holds for unslipped reads only. An allele at −9
  that slips two repeats longer reports −7; one at −15 reports −8. The fit models the −9 class as
  exactly −9, so alleles further out look like −9 alleles with fewer longer slips. The effect is
  small, because longer slips are about 2 reads in 100 here. The claim should be stated as
  approximate.
- **`examples/ng_ssr_strata_threads.rs:57-60`.** This says "Six is production's `ALLELE_SPAN`:
  thirteen allele classes and the ninety-one genotypes". The example draws at a span of 6, with a
  read span of 6 as well, but fits at the default span, now 9 (`..SsrFitConfig::default()`, line
  158). Its timings are no longer measured at the shape its doc describes.
- **`examples/ng_ssr_subset_against_every_sample.rs:156-157`.** The doc says "91 genotype pairs …
  at thirteen allele classes". The code computes the count from `config.allele_span` (line 228),
  so only the comment is stale.
