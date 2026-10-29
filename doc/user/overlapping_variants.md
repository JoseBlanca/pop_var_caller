# Overlapping variants

When two variants overlap, such as a SNP inside a deletion, a VCF can write them in two ways.
Callers differ in which way they use, and population-genetics tools read the two ways very
differently. One of them, separate records, can distort the results.

Let's imagine a reference `ATT` at positions 1 to 3 of `chrom1`, and three individuals:

| Individual | Haplotypes |
|---|---|
| ind1 | `A` / `A`: both copies carry a deletion of `TT` |
| ind2 | `ACT` / `ACT`: both copies carry a SNP, `T` to `C` at position 2 |
| ind3 | `A` / `ACT`: one copy of each |

Written **as one record**, whose alleles are the whole haplotypes. freebayes and pvc write it
this way:

```
#CHROM  POS  REF  ALT     ind1  ind2  ind3
chrom1  1    ATT  A,ACT   1/1   2/2   1/2
```

Written **as separate records**, one per event, each at its own position. GATK HaplotypeCaller
writes it this way. The SNP at position 2 lies inside the deletion, so its record needs the `*`
allele, which means "this position is deleted on this copy by a deletion that starts upstream":

```
#CHROM  POS  REF  ALT   ind1  ind2  ind3
chrom1  1    ATT  A     1/1   0/0   0/1
chrom1  2    T    C,*   2/2   1/1   1/2
```

Both say the same thing. The difference is what a downstream tool does with the second record.

## What goes wrong with the separate-records form

Most population-genetics tools expect biallelic SNPs and do not understand `*`. They usually
handle it in one of three ways, and each one distorts the results.

**`*` read as missing.** At position 2, ind1 becomes `./.` and ind3 becomes half or wholly
missing.

- The data go missing *because* an individual carries the deletion. Missing-data filters then
  remove sites or individuals according to the deletion's frequency, and imputation is biased.
- Among the copies left, `C` has frequency 1, so the site looks fixed, although `C` is on only
  half of the copies.
- ind3, a true heterozygote at this locus, is lost.

**`*` recoded as the reference allele.** ind1 becomes `T/T` and ind3 `T/C` at position 2. This
invents an allele, `T`, that no copy in the sample carries, and counts the deletion a second
time, as a SNP.

**`*` counted as another allele.**

- If `C` and `*` are both taken as "not reference", every individual becomes `1/1`. ind3, a
  heterozygote, becomes a false homozygote, and the site again looks fixed.
- If the site is kept with three alleles, the deletion is counted twice: once at position 1 and
  again as `*` at position 2, and once more for every other variable position it covers. That
  adds segregating sites, and it also creates perfect linkage disequilibrium between the two
  records.

Across a genome these sites are few, so genome-wide averages move little. But they cluster
where indels cluster: repeats, low-complexity sequence and highly polymorphic regions. **Window
scans are affected far more than averages**: windowed π and Tajima's D, F<sub>ST</sub> outlier
scans, selection scans. If a deletion's frequency differs between populations, so does the
missing data it causes, which can create a false signal of differentiation or hide a real one.

The common workaround, keeping only biallelic SNPs and dropping SNPs close to indels (for
example with `bcftools filter --SnpGap`), avoids these artifacts, but at a cost. It
systematically loses SNPs on deletion-carrying copies and under-samples diversity in exactly
those regions. And it throws away the sites with more than two alleles, which are the most
informative ones, because a site with more alleles tells more individuals apart.

## What pvc writes

- **Overlapping SNPs and indels go in one record**, whose alleles span all of them, as in the
  first example above.
- **An STR is one record**, whose alleles are the whole repeat sequences. Variants inside a repeat
  are never split into separate records.
- **Variants that do not overlap get separate records**, even when they are a base apart. Unlike
  freebayes, pvc does not merge them. This causes none of the problems above: neither event
  covers the other's position, so no `*` is needed.
- **One exception, at the edge of an STR.** When a deletion runs from an STR into the sequence
  beside it, pvc widens the STR record over those bases, and a SNP or indel record that overlaps
  them gets `*` for the copies the deletion removes. Removing this exception is in progress.

To see how many records with `*` a pvc VCF holds:

```
bcftools view -H -i 'ALT="*"' calls.vcf.gz | wc -l
```
