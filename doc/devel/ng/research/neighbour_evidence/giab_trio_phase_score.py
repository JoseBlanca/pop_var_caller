"""Scores the phase pvc writes for GIAB's HG002 against the phase its parents' genotypes give
(read phase spec §7.2), where the GIAB truth carries no phase of its own.

usage: uv run python giab_trio_phase_score.py <HG002 calls.vcf.gz> <child truth> <father truth> <mother truth>

The truth VCFs are GIAB v4.2.1's, cut to regions all three samples are confident in, so a site a
parent's truth does not list is homozygous reference in that parent. At a SNP where the child is
heterozygous and a parent homozygous, the child's allele from that parent is known, and the other
is from the other parent; a site where both parents are heterozygous is left out, and so is one
whose genotypes break Mendel's rule. A link is two consecutive phased records of the child in one
phase set, both such sites and both called with the truth's alleles; it is *wrong* when the copy
written first comes from the father at one and from the mother at the other.

Prints: true heterozygous SNPs, how many the trio resolves, how many the calls phase, the links
scored and the wrong ones, and the links the trio cannot score.
"""
import gzip
import sys


def snps(path):
    """{(contig, pos): (ref, alt, (allele, allele))} for single-base records."""
    out = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            alleles = [fields[3]] + fields[4].split(',')
            if any(len(allele) != 1 for allele in alleles):
                continue
            gt = fields[9].split(':')[0].replace('|', '/').split('/')
            if '.' in gt:
                continue
            out[(fields[0], int(fields[1]))] = (alleles, tuple(alleles[int(a)] for a in gt))
    return out


def paternal_base(child, father, mother):
    """The child's base from the father, or None when the trio does not say."""
    a, b = child
    if a == b:
        return None
    for parent, other in ((father, mother), (mother, father)):
        if parent[0] == parent[1] and parent[0] in (a, b):
            from_this = parent[0]
            from_other = b if from_this == a else a
            if from_other not in other:
                return None  # a Mendelian error
            return from_this if parent is father else from_other
    return None


def main():
    calls, child_path, father_path, mother_path = sys.argv[1:5]
    child, father, mother = snps(child_path), snps(father_path), snps(mother_path)
    paternal = {}
    true_hets = 0
    for key, (_, pair) in child.items():
        if pair[0] == pair[1]:
            continue
        true_hets += 1
        ref = child[key][0][0]
        f = father.get(key, (None, (ref, ref)))[1]
        m = mother.get(key, (None, (ref, ref)))[1]
        base = paternal_base(pair, f, m)
        if base is not None:
            paternal[key] = base
    phased = []  # (contig, pos, ps, copy-0 base, truth-het?)
    with gzip.open(calls, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            alleles = [fields[3]] + fields[4].split(',')
            keys = fields[8].split(':')
            values = fields[9].split(':')
            if 'PS' not in keys or '|' not in values[0] or any(len(a) != 1 for a in alleles):
                continue
            first, second = (alleles[int(x)] for x in values[0].split('|'))
            key = (fields[0], int(fields[1]))
            truth = child.get(key)
            right_genotype = truth is not None and sorted(truth[1]) == sorted([first, second])
            phased.append((key, values[keys.index('PS')], first, right_genotype))
    links = wrong = unscored = wrong_genotype_links = 0
    for (k0, ps0, c0, g0), (k1, ps1, c1, g1) in zip(phased, phased[1:]):
        if k0[0] != k1[0] or ps0 != ps1:
            continue
        if not (g0 and g1):
            wrong_genotype_links += 1
            continue
        if k0 not in paternal or k1 not in paternal:
            unscored += 1
            continue
        links += 1
        wrong += (c0 == paternal[k0]) != (c1 == paternal[k1])
    print(f'{calls}: {true_hets} true heterozygous SNPs, {len(paternal)} resolved by the trio; '
          f'{len(phased)} phased SNP columns; links scored {links}, wrong {wrong}; '
          f'{unscored} links the trio cannot score; {wrong_genotype_links} links touching a wrong '
          f'genotype')


if __name__ == '__main__':
    main()
