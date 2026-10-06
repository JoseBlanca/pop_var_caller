# Explaining a run's decisions at chosen loci (`--explain-loci`)

## 1. What it is for

A calling run decides many things at a locus between the evidence a psp stored and the record a VCF
holds: which candidate alleles to keep, which samples can be called, which genotype each sample
gets and how sure it is, what the site quality is after the artifact penalties, whether the
hidden-paralog filter removes the record, and whether the record is written at all. When a variant
is missing or wrong in the VCF, the reason is almost always one of those decisions, and none of
them is visible afterwards. `inspect-psp` shows what reached the psp; this shows what the run then
did with it.

It was asked for by a user tracing GIAB false negatives (2026-10-03). Of 27 indels missed on HG002,
HG003 and HG004 at 300×, they could trace 8 to the psp; the rest were lost after it.

## 2. Interface

```text
pop_var_caller call-from-psps       ... --explain-loci <regions.bed> <out.tsv>
pop_var_caller call-from-alignments ... --explain-loci <regions.bed> <out.tsv>
```

**It selects what is explained, never what is called.** The run calls the whole of its ground
exactly as it would without the option, and writes the same VCF byte for byte. Every locus the run
meets that overlaps a region of the BED gets its explanation written to `<out.tsv>`.

## 3. What it may cost

**Off, which is the default, it costs a normal run nothing it can measure.** That is a condition the
owner set on building it (2026-10-03). Concretely:

- every place that would explain something asks one question first — *is this locus in an explained
  region?* — and the answer is a constant `false` when the option is off;
- nothing is collected, allocated or formatted for a locus that is not explained;
- no decision is changed, reordered or recomputed to make it explainable. What is recorded is read
  off values the run computes anyway; where a value would otherwise be discarded before it can be
  read, it is read at the point it exists.

**On, it costs in proportion to the loci explained**, which a BED of a few hundred regions keeps
small. A BED covering a whole genome would explain every locus and is not what the option is for.

**Measured (2026-10-03)**, `call-from-psps` over four tomato accessions at about 3× and 20 regions,
ten runs each with the files cached: `main` before this work took a median of **3.30 s** (3.25 to
3.34), this build with the option off **3.28 s** (3.17 to 3.49), and with it on over 2 kb about
3.24 s. The difference is inside the run-to-run spread. The three VCFs are the same file.

## 4. The output

A TSV in **long form**: one row per fact, so that every kind of fact has the same columns and a
script filters on `step` rather than parsing a different layout per kind.

| column | what it holds |
|---|---|
| `contig`, `start`, `end` | the locus, 1-based and both ends included, as the run built it |
| `step` | which decision the row is about (below) |
| `sample` | the sample the row is about, or `.` for a fact about the whole locus |
| `subject` | what within the step: an allele's bases, a genotype, a partial read's bases, or `.` |
| `value` | the decision or the number |
| `detail` | supporting numbers as `key=value` pairs separated by `;`, or `.` |

The header carries `#explain_loci_columns_version=1`; it changes whenever a column or a step's
meaning does.

**Rows are in genome order of their locus, then in the order of the steps below.**

### 4.1 The steps

| `step` | one row per | `value` | `detail` |
|---|---|---|---|
| `locus` | locus | `generic`, `repeat_tract` or `repeat_bundle` | `alleles`, `covering_samples` |
| `allele` | allele the merge assembled | `kept` (as candidate *k*), `below_support`, `cut_by_allele_cap`, or `refused_not_periodic` | `cohort_reads`, `best_sample`, `best_sample_reads`, `best_sample_share` |
| `reads` | sample × allele with reads | complete reads showing it | `read_groups` |
| `partial` | sample × partial observation | reads | `witnessed`, `compatible_with` (candidate alleles it fits) |
| `callable` | sample | `callable` or `set_aside` | why, when set aside |
| `likelihood` | sample × genotype | the genotype's log-likelihood, natural log | . |
| `genotype` | sample | the called GT | `GQ`, `missing`, `reads_were_uninformative` |
| `site` | locus | QUAL as written | `uncorrected_qual`, `ABPEN`, `SPPEN`, `converged`, `passes`, `filter` |
| `paralog` | locus scored by the hidden-paralog filter | `kept`, `tagged` or `dropped` | `LR`, `posterior` |
| `outcome` | locus | `written`, `not_written`, `below_min_site_quality`, `strand_bias_at_or_above_cutoff`, `dropped_by_paralog_filter`, `nobody_to_call`, `tract_without_whole_repeats`, `bundle_set_aside`, `too_quiet`, `too_wide`, `over_depth_ceiling` | the reason in words |

**How partial reads were weighted** is the `partial` rows: a partial read compatible with a
genotype's alleles contributes the shares of the copies carrying them, one compatible with none of
them is charged as a sequencing error (`read_likelihoods.md` §5.3). The `compatible_with` list is the
whole of what decides which.

### 4.2 Loci the caller never saw

Three kinds of locus are dropped by the merge before calling, and each is explained with one `locus`
row and one `outcome` row:

- **too quiet** — no sample showed enough reads disagreeing with the reference (`cohort_merge.md`
  §3.3);
- **too wide** — wider than `--max-cohort-locus-span`;
- **over the depth ceiling** — some sample's read group is deeper than `--max-read-group-depth`.

A region of the BED in which no locus was built at all writes nothing: either no sample's psp holds
a record there, which `inspect-psp` shows, or the ground is outside the run's analysed regions.

## 5. Threads and order

A locus's rows come from three places at three times: the merge, while it assembles loci; calling,
which runs on several threads in `call-from-psps`'s default round driver; and the hidden-paralog
filter, which runs after every locus is called. **So the rows are collected in memory and written
once, sorted, when the run ends** — which is also when a locus's `outcome` is known, since a record
calling handed to the VCF may still be dropped by the filter. Holding them costs memory in
proportion to the loci explained and nothing otherwise.

**Calling's rows travel with the locus's own outcome** to the thread that folds outcomes back in
genome order, and are collected there, so no lock is taken on the calling threads.

## 6. Built in three steps

1. **Calling** — `locus` through `site`, and the `outcome` of every locus that reaches the caller.
2. **The hidden-paralog filter** — `paralog` rows, and the `dropped_by_paralog_filter` outcome.
3. **The merge's drops** — §4.2.

Each step is checked for the cost rule of §3: the VCF is byte-identical with the option on and off,
and a run with the option off takes the same time as before, on the tomato identity oracle's
cohort.
