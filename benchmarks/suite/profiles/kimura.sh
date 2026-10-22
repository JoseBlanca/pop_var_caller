# Machine profile for kimura (x86-64 Linux, podman, 128 GB): sourced by run_suite.sh when it runs
# on kimura. The owner allows the suite at most 32 cores and 80 GB.
BENCH_THREADS=32
BENCH_MEM=80g
# The tomato collection: whole-genome CRAMs of about 2,100 accessions, one directory per project.
# The cohort is found in it by sample name (bench/speed.py); its reference beside it.
BENCH_TOMATO_CRAMS=/home/joxi/tomato25/snp_calling/crams
BENCH_TOMATO_REFERENCE=/home/joxi/tomato25/genome/S_lycopersicum_chromosomes.4.00.fa
# Both are under this directory, mounted read-only into the dev container and the DeepVariant jobs.
BENCH_EXTRA_MOUNT=/home/joxi/tomato25
# Speed is measured over all of chromosome 1 (90.9 Mb), not tomato1's 8 Mb: whole-genome CRAMs
# make larger ground possible, and 8 Mb is too little for start-up to stay small. Relative to the
# checkout.
BENCH_TOMATO_REGIONS=benchmarks/suite/data/tomato/SL4.0ch01.bed
# Sample counts for the speed curve: all real accessions on kimura, no copies.
BENCH_TOMATO_SIZES=1,2,5,10,20,50,100
