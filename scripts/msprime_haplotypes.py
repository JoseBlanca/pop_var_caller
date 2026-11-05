"""Haplotypes of a natural population for `ng_str_stutter_benchmark haplotypes=<file>`.

A coalescent simulation with recombination (msprime), so the cohort has what a real population
has and a mosaic of a few founders lacks: linkage that fades gradually with distance and many rare
haplotypes. The simulator lays its own events (STR segments, SNPs, indels) and borrows the
genealogy here: each event takes the next sites at or after its position.

usage: uv run --no-project --python 3.12 --with msprime==1.3.3 python scripts/msprime_haplotypes.py \
           <out file> --samples N --contigs C --contig-bp B [--inbreeding F] [--seed S]

**Selfing** is modelled the standard way (Nordborg and Donnelly 1997): a population selfing to an
inbreeding coefficient F coalesces as an outbred one of effective size Ne / (1 + F) whose
recombination rate is scaled by (1 − F). The autozygous runs themselves — the two copies of one
sample being one haplotype — are laid by the simulator (`inbreeding`, `run_kb`), not here.

**Site density**: the mutation rate is set so the sites fall on average `--site-spacing` bases
apart, far denser than the simulator's events, so each event finds a site close to it. These
mutations only carry the genealogy; the alleles the reads show are the simulator's.

Writes a header `#contig_bp=B<TAB>haplotypes=2N<TAB>...settings`, then one line a site:
`contig<TAB>position<TAB><one 0 or 1 a haplotype>`, contigs numbered from 0.
"""

import argparse

import msprime
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("out")
    parser.add_argument("--samples", type=int, required=True, help="diploid samples")
    parser.add_argument("--contigs", type=int, required=True)
    parser.add_argument("--contig-bp", type=int, required=True)
    parser.add_argument("--ne", type=float, default=10_000, help="effective size if outbred")
    parser.add_argument("--recombination-rate", type=float, default=1e-8, help="per base a generation")
    parser.add_argument("--inbreeding", type=float, default=0.0, help="F from selfing")
    parser.add_argument("--site-spacing", type=float, default=100.0, help="mean bases between sites")
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()

    ne = args.ne / (1.0 + args.inbreeding)
    recombination = args.recombination_rate * (1.0 - args.inbreeding)
    haplotypes = 2 * args.samples
    harmonic = sum(1.0 / i for i in range(1, haplotypes))
    mutation = 1.0 / (args.site_spacing * 4.0 * ne * harmonic)

    with open(args.out, "w") as out:
        out.write(
            f"#contig_bp={args.contig_bp}\thaplotypes={haplotypes}\tne={args.ne}\t"
            f"recombination_rate={args.recombination_rate}\tinbreeding={args.inbreeding}\t"
            f"site_spacing={args.site_spacing}\tseed={args.seed}\tmsprime={msprime.__version__}\n"
        )
        for contig in range(args.contigs):
            ancestry = msprime.sim_ancestry(
                samples=args.samples, ploidy=2, population_size=ne,
                sequence_length=args.contig_bp, recombination_rate=recombination,
                random_seed=args.seed + 2 * contig)
            sites = msprime.sim_mutations(
                ancestry, rate=mutation, model=msprime.BinaryMutationModel(),
                random_seed=args.seed + 2 * contig + 1)
            for variant in sites.variants():
                carriers = variant.genotypes.astype(np.uint8)
                if carriers.min() == carriers.max():
                    continue  # back-mutated to no variation
                out.write(f"{contig}\t{int(variant.site.position)}\t"
                          f"{(carriers + ord('0')).tobytes().decode()}\n")


if __name__ == "__main__":
    main()
