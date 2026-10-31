"""Inputs for genome-population-simulator (GPS): a reference segment, SNPs, and its configuration.

GPS is a colleague's simulator of diploid populations, kept in its own repository
(`~/genome_population_simulator`). The suite runs a release build of it from this checkout's `tmp/`
(GPS_BIN). It places SNPs at the frequencies it is given, draws each individual's genotypes under an
inbreeding coefficient F set per population, and writes the truth and paired Illumina reads.

**How GPS applies F, and why it suits this test.** F acts site by site and independently: each site
has the right shortage of heterozygotes, but there are no runs of homozygosity and every individual
of a population has the same expected F. pvc estimates a sample's F the same way, as its average
shortage of heterozygotes over positions (src/parameter_estimation/mod.rs), so the simulation holds
what pvc's estimator can see.

**Several F values in one cohort** are several GPS populations given the same SNPs at the same
frequencies (one table row per SNP and population), so F is the only thing that differs between
them.
"""

import os
import random

GPS_BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))), "tmp", "gps_build", "x86_64-unknown-linux-gnu", "release",
    "genome-population-simulator")

TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}


def group_id(f):
    """A population's name from its F: F000, F025, ..., F100."""
    return f"F{round(f * 100):03d}"


def sample_names(fs, per_group):
    return [f"{group_id(f)}_{i:02d}" for f in fs for i in range(1, per_group + 1)]


def group_of(sample):
    """The F a sample was simulated with, read from its name."""
    return int(sample.split("_")[0][1:]) / 100


def read_segment(fasta_lines):
    return "".join(l.strip() for l in fasta_lines if not l.startswith(">")).upper()


def write_snvs(path, contig, seq, fs, density, min_af, max_af, seed):
    """SNPs at about `density` a base, the same for every population.

    Positions are drawn without replacement over bases that are A, C, G or T. The alternative base
    is the transition two times in three and either transversion otherwise (a transition to
    transversion ratio of about 2, as in plants). Each SNP's ALT frequency is uniform between
    `min_af` and `max_af`: frequent alleles give every F group heterozygotes and both homozygotes to
    be called, where a neutral spectrum would leave most genotypes homozygous reference.
    Returns the number of SNPs.
    """
    rng = random.Random(seed)
    usable = [i for i, b in enumerate(seq) if b in "ACGT"]
    count = round(len(seq) * density)
    positions = sorted(rng.sample(usable, count))
    with open(path, "w") as fh:
        fh.write("contig\tposition\tref\talt\tpopulation_id\talt_af\n")
        for i in positions:
            ref = seq[i]
            if rng.random() < 2 / 3:
                alt = TRANSITIONS[ref]
            else:
                alt = rng.choice([b for b in "ACGT" if b not in (ref, TRANSITIONS[ref])])
            af = round(rng.uniform(min_af, max_af), 4)
            for f in fs:
                fh.write(f"{contig}\t{i + 1}\t{ref}\t{alt}\t{group_id(f)}\t{af}\n")
    return count


def _udi(rng, seen):
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(8))
        if s not in seen:
            seen.add(s)
            return s


def write_config(path, reference_name, contig, length, snvs_name, fs, per_group, coverage, seed):
    """GPS's general.toml: one population per F, PCR-free 150 bp pairs at `coverage`.

    PCR-free, so duplicates do not enter a test of the genotype prior; sequencing errors at GPS's
    default (enabled). `seed` makes the run repeatable.
    """
    rng = random.Random(seed)
    seed_hex = "".join(rng.choice("0123456789abcdef") for _ in range(64))
    names = sample_names(fs, per_group)
    lines = [
        'config_kind = "general"',
        'schema_version = "4"',
        f'seed = "{seed_hex}"',
        "",
        "[population_snvs]",
        'mode = "supplied"',
        f'input_path = "{snvs_name}"',
        "",
        "[reference]",
        'kind = "supplied"',
        'id = "tomato-segment"',
        f'fasta = "{reference_name}"',
        "",
        "[[reference.sequences]]",
        f'id = "{contig}"',
        'class = "chromosome"',
        f"length = {length}",
    ]
    for f in fs:
        ids = [n for n in names if n.startswith(group_id(f) + "_")]
        lines += ["", "[[initial_populations]]", f'population_id = "{group_id(f)}"',
                  f"individual_count = {len(ids)}",
                  "individual_ids = [" + ", ".join(f'"{i}"' for i in ids) + "]", "ploidy = 2"]
        if f == 0:
            lines.append('genotype_model = "hwe"')
        else:
            lines += ['genotype_model = "f_diploid"', f"inbreeding_coefficient = {f}"]
    lines += ["", "[vcf_output]", 'mode = "text"', "", "[technical_event]",
              'event_id = "reads"', "individual_ids = [" + ", ".join(f'"{i}"' for i in names) + "]",
              f"target_coverage_wgs = {coverage}", 'library_preparation = "pcr_free"']
    seen = set()
    for n in names:
        lines += ["", "[[technical_event.libraries]]", f'individual_id = "{n}"',
                  f'i7 = "{_udi(rng, seen)}"', f'i5 = "{_udi(rng, seen)}"']
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return names
