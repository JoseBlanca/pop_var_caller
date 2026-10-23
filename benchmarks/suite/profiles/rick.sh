# Machine profile for rick (x86-64 Linux, podman, 32 cores, 125 GB): sourced by run_suite.sh when it
# runs on rick. Thirty cores, not all 32, so the host's own load and the suite's sampler stay off the
# cores the callers are timed on; 80 GB as on kimura.
BENCH_THREADS=30
BENCH_MEM=80g
# The tomato cohort: 100 accessions copied from kimura's collection (tomato1's 62, then the next 38
# of project PRJNA454805 by sample name, as bench/speed.py orders them), chromosome 1 only. The
# reference is the machine's own copy of SL4.0.
BENCH_TOMATO_CRAMS=/home/joxi/tomato25/bench_crams
BENCH_TOMATO_REFERENCE=/home/joxi/tomato25/genome/S_lycopersicum_chromosomes.4.00.fa
# GIAB's 300x alignments of HG002, HG003 and HG004, already on this machine: the novoalign BAMs the
# suite would otherwise fetch from GIAB, as CRAMs of chromosomes 1 to 22 (every benchmark region is
# autosomal), on the data disk.
BENCH_GIAB_ALIGNMENTS='/media/giab/full_coverage/{sample}.chr1_22.300x.cram'
# Both data trees, mounted read-only into the dev container and the DeepVariant jobs.
BENCH_EXTRA_MOUNT=/home/joxi/tomato25:/media/giab
# Speed over all of chromosome 1 (90.9 Mb), as on kimura. Relative to the checkout.
BENCH_TOMATO_REGIONS=benchmarks/suite/data/tomato/SL4.0ch01.bed
BENCH_TOMATO_SIZES=1,2,5,10,20,50,100
