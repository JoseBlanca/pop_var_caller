# A substitution haplotype is no longer re-read as a longer allele cut short (2026-10-30)

## What was wrong

pvc's calling has a rule for reads the mapper laid straight across the start of an insertion
(`reread_spellings_cut_short`, `src/run/records.rs`; `read_likelihoods.md` §5.3). Such a read
covers the whole locus but shows only the first bases of the longer allele, written as
substitutions. The rule turns every read whose spelling is exactly the start of a longer allele
into a partial read of that allele.

The rule also caught real haplotypes. At GIAB HG002 chr22:37788634 the truth is a heterozygous
`AACA` insertion and, 12 bases on, a homozygous `G>A`. The walk put both into one 13-base record:

| reads | spelling over the record | what it is |
|---|---|---|
| 16 | `GAACAAACAAACA` | the substitution alone |
| 18 | `GAACAAACAAACAAACA` | the insertion and the substitution |

The 16 reads spell exactly the first 13 bases of the 18. So the rule re-read them as the longer
allele cut short, the substitution-only allele was left with no reads, and the site was called
`1/1` with reads 0 and 18. DeepVariant called it `0/1` with 23 and 17. The reads themselves settle
it: the 16 continue 31 to 137 bases past the record, matching the reference.

The same pattern sat at chr22:17595058 (`GTTTTT` beside `GTTTTTTT`, truth `0/1` insertion plus `1/1`
substitution). There only the rule's read-count guard kept it off: 19 substitution reads against 18.

## The change

A fourth guard, `is_the_longer_allele_without_its_insertion`. A read cut short at the right border
spells the start of the longer allele over the reference, so from the first base where it departs
from the reference it shows the inserted copy, while the longer allele ends with reference bases
there. A haplotype carrying the longer allele's substitutions but not its insertion departs at those
substitutions, and the longer allele ends with exactly that stretch. So a spelling as long as the
reference is not re-read when the longer allele ends with it from its first departure on. The left
border is the mirror image. Spellings that carry an insertion of their own are untouched, and so
are deletion alleles, which have no insertion to take out.

**The sequences alone are not enough (found in review).** At a homozygote carrying both the
insertion and the substitution, reads cut short spell exactly what the substitution alone would,
so the sequence test passes for them too and would bring back the allele no haplotype carries,
the failure the rule was written to fix. The number of reads separates the two: reads cut short
are a minority of the carriers (34 against 107 at HG003 chr15:96140584), while one heterozygote's
two haplotypes bring about as many each (16 and 18 at chr22:37788634). So the spelling's reads stay
complete only when the sequence test passes and they are at least half as many as the longer
allele's. The decisive evidence is per read — a substitution read runs on along the reference past
the border, a read cut short has stopped — but the walk does not store how far a read ran past a
record, so that would need a psp format change. Still unresolved: a substitution haplotype with
fewer than half the longer allele's reads, and one whose substitution lies before the insertion
inside a run the insertion extends, which is a read cut short base for base.

The rule's existing tests all keep their answers. Three new tests:
`a_substitution_haplotype_the_longer_allele_carries_too_is_left_alone` (fails without the guard),
its left-border mirror, and `a_substitution_spelling_far_outnumbered_by_the_longer_allele_is_reread`
for the homozygote. The spec gains a paragraph under §5.3, after the rule's own measurements, and
its guard count reads four.

## Measured

HG002, chromosome 22, the colleague's benchmark inputs on rick (`~/giab_wgs_benchmark`), scored with
hap.py exactly as that benchmark scores. Both builds call from the same evidence, stored by main's
walk; only calling differs.

| depth | build | indels missed | false indels | indel F1 | SNPs missed | false SNPs |
|---|---|---|---|---|---|---|
| 30× | main 1b46e6f9 | 170 | 76 | 0.98106 | 327 | 339 |
| 30× | this change | 169 | 73 | 0.98137 | 327 | 339 |
| 5× | main 1b46e6f9 | 1,885 | 505 | 0.79579 | 7,143 | 1,659 |
| 5× | this change | 1,884 | 503 | 0.79603 | 7,144 | 1,660 |

At 30× the only error site that changes is chr22:37788634, now right. At 5× that site is fixed and
one homozygous SNP is lost, chr22:23873331, which deserves a note. There the truth is again an
insertion beside a homozygous SNP, with 2 reads of each haplotype. Main re-read the SNP's 2 reads as
cut short, called the insertion haplotype `1/1`, and so got the SNP right by accident. With the
guard both haplotypes keep their reads, and the reads favour `1/2` (log-likelihood −2.77 against
−5.30 for `0/1`), but the call is `0/1`, a genotype whose reference allele has no reads. The genotype
prior outweighs two-plus-two reads against a two-alternative heterozygote; that is a question for
the prior, not for this rule.

**The cohort corner.** 29 tomato accessions at about 3× (PRJNA454805, the first 30 CRAMs of the
rick benchmark set; two belong to one accession) over SL4.0ch01:1–10,000,000, both builds calling
from the same stored evidence: 20 of 207,606 records change. There is no truth set, so each changed
site was read against GATK HaplotypeCaller's joint calls for the same accessions (the suite's
`tomato_100` run), sample by sample, by hand. About 9 sample genotypes move towards GATK and about 4
away from it. The largest gain is SL4.0ch01:2906008, where GATK has 8 samples carrying the `T>C`
without the neighbouring insertion and main called 7 of them homozygous for the insertion; the fix
calls them `T>C`. The losses (1099337, 2744775, 3895439) share one shape: the substitution extends a
homopolymer the insertion also lengthens (`ATTTTG` → `ATTTTT` beside `ATTTTTT`), where at about 3×
single reads pass the half-the-reads test and reads cut short become a substitution allele GATK does
not call. That is the per-read limit stated above, at its worst at low depth.

Of the seven loci in the error set that showed the re-read pattern, two improve (37788634 fixed;
47934937 now `1/2`, as the truth has it) and one now carries both haplotypes but is still scored
wrong because one insertion falls in the neighbouring repeat-tract record (34479793). The other four
show the pattern only among reads that are not the deciding ones, and do not change.

## Tried and dropped

Widening the walk's record to the indel's slide along its carriers' own bases (the substitution lets
the insertion slide past the reference's stopping point) fixed the same site but, at 30×, merged a
GC-rich artifact haplotype with two true homozygous SNPs at chr22:50170399 and lost both, and it
cannot widen into repeat ground, so it does not reach 17595058. Branch `partial-prefix-walk-attempt`
keeps it.
