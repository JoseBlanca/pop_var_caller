#!/usr/bin/env bash
# **The gate for a behaviour-preserving change** — a refactor, a cleanup, a rename. It answers two
# questions about the pvc built from the working tree, in about six minutes on the owner's laptop:
#
#   1. Is every output byte-identical to the recorded baseline?
#   2. Did any step get slower, or use more memory, beyond the noise measured when the baseline
#      was recorded?
#
# Usage:
#
#   scripts/code_quality_oracle.sh              # check against scripts/code_quality_oracle.baseline
#   scripts/code_quality_oracle.sh --identity   # one pass, outputs only, no timing (about 2 min)
#   scripts/code_quality_oracle.sh --record     # run five times and write a new baseline
#
# Exit status: 0 identical and within noise; 1 an output differs from the baseline, or two passes
# of the same build disagree; 3 outputs identical but a step got slower or larger; 2 the run could
# not be made (missing input, failed build, a step that exited non-zero).
#
# **What runs.** One pass is thirteen steps, each a separate `pvc` process, at --threads 4:
#
#   catalog-tomato, catalog-grch38       repeat-catalog over the whole tomato and human references
#   tomato16-*                           16 tomato accessions over one 100 kb region: generate-psps,
#                                        estimate-parameters, then call-from-psps and
#                                        call-from-alignments with the fitted file, and
#                                        call-from-psps again with --refine-genotypes on
#   tomato63-call-from-alignments        all 63 accessions over two regions, --defaults
#   tomato1-call-from-alignments         one accession over all 80 regions (8 Mb), --defaults
#   hg002-*                              GIAB HG002 at about 290 reads a position over the first 50
#                                        benchmark regions (317 kb): generate-psps, a one-sample
#                                        estimate-parameters, call-from-psps and
#                                        call-from-alignments with that fit
#
# The tomato reads are those of benchmarks/tomato1 (about 11 reads a sample at a called record in
# these regions); the human reads are benchmarks/human_genome_bottle's. Both are git-ignored, so
# they are read from the main checkout (`CQ_DATA`, default the `benchmarks/` beside the git common
# directory); references come from `CQ_GENOMES` (default ~/genomes).
#
# **What is compared, and what is left out of it.** Checksums of: both catalogs; every psp; both
# parameters files; every VCF; the window coverage and per-sample histograms each calling step
# writes (NG_WINDOW_COVERAGE_FILE, the hidden-duplication filter's input, which never reaches the
# VCF). Left out, because they legitimately differ between two runs of one build:
#
#   VCF   the `##commandline=` and `##pvcLog=` header lines
#   psp   the header's `created =` and `command-line =` lines (the timestamp, and the command)
#   windows  row order (the rows are sorted first; the format does not promise an order)
#
# Every input is reached through a symbolic link under `CQ_WORK/in` and every path on a command
# line is relative, so the command lines — which a psp header and a VCF header record — are the
# same in every worktree. Two checks come free with the passes: each pass after the first must
# reproduce the first pass's checksums (a run-to-run difference is exit 1), and the two calling
# routes must agree with each other (call-from-psps and call-from-alignments with one parameters
# file write the same VCF, on both tomato16 and hg002).
#
# **How time and memory are judged.** Every step is measured with /usr/bin/time: wall seconds, peak
# resident memory, and on macOS the CPU instructions retired, which hardly move with machine load
# and so catch a slowdown that wall time hides in its noise. A check makes three passes and takes
# each step's fastest (minimum) figure; the baseline holds the minimum, median and maximum of its
# own five passes. A step is flagged when its minimum exceeds the baseline's minimum by more than
# its tolerance: 2%, or the baseline's own spread ((max − min) / min) where that is larger. The
# whole run's total wall time is judged the same way. The tolerance used is printed beside every
# step, so a reader can see which rule applied.
#
# The timing half is only meaningful on the machine and thread count it was recorded on, and on a
# quiet machine: the baseline names its host, and a check on another host compares outputs only.
# The 1-minute load average is printed before each pass; above 4 the wall times are suspect (the
# laptop has six fast cores).
#
# **Recording a new baseline** — needed after an owner-approved change that moves an output, or
# when the timing has shifted for a reason outside the code (a new OS, a new toolchain):
#
#   scripts/code_quality_oracle.sh --record
#   git diff scripts/code_quality_oracle.baseline    # which checksums moved
#
# and say in the commit which lines moved and why. Recording also keeps a copy of every output in
# `CQ_WORK/outputs/baseline`, which a later check diffs against to show *what* changed, not only
# that something did. A worktree that has never recorded has no such copy and reports checksums
# alone; record once on the base commit to get one.
#
# Settings (environment):
#   CQ_BIN      a pvc binary to use instead of building one (CARGO_TARGET_DIR=target-host, release)
#   CQ_REPS     passes to make (default 3 to check, 5 to record)
#   CQ_WORK     working directory (default tmp/code_quality_oracle)
#   CQ_BASELINE the baseline file (default scripts/code_quality_oracle.baseline)
#   CQ_DATA, CQ_GENOMES  where the reads and the references are (see above)
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)

mode=check
case "${1:-}" in
  "") ;;
  --record) mode=record ;;
  --identity) mode=identity ;;
  -h|--help) sed -n '2,${/^#/!q;p;}' "$0"; exit 0 ;;
  *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
esac

threads=4
work=${CQ_WORK:-$root/tmp/code_quality_oracle}
baseline=${CQ_BASELINE:-$root/scripts/code_quality_oracle.baseline}
genomes=${CQ_GENOMES:-$HOME/genomes}
if [ -z "${CQ_DATA:-}" ]; then
  common=$(git -C "$root" rev-parse --path-format=absolute --git-common-dir)
  CQ_DATA=$(cd "$common/.." && pwd)/benchmarks
fi
data=$CQ_DATA
case $mode in
  record) reps=${CQ_REPS:-5} ;;
  identity) reps=1 ;;
  *) reps=${CQ_REPS:-3} ;;
esac

die() { echo "code_quality_oracle: $*" >&2; exit 2; }

# The working directory is wiped between passes, so it must be one this script owns.
case "$work" in
  "$root"/tmp/*) ;;
  *) die "CQ_WORK must lie under $root/tmp, not $work" ;;
esac

if command -v md5sum >/dev/null 2>&1; then
  md5_stdin() { md5sum | cut -d' ' -f1; }
else
  md5_stdin() { md5 -q; }
fi

# ---------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------
tomato_fa=$genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa
human_fa=$genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
tomato_crams=$data/tomato1/crams
tomato_bed=$data/tomato1/regions.bed
hg002_cram=$data/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram
hg002_bed=$data/human_genome_bottle/HG002_bench_azar_merged_1000.bed

for needed in "$tomato_fa" "$tomato_fa.fai" "$human_fa" "$human_fa.fai" "$tomato_bed" \
              "$hg002_cram" "$hg002_cram.crai" "$hg002_bed" "$tomato_crams/SRR7279481.p1.bench.cram"; do
  [ -r "$needed" ] || die "missing input $needed (set CQ_DATA / CQ_GENOMES)"
done
tomato_count=$(ls "$tomato_crams"/*.bench.cram | wc -l | tr -d ' ')
[ "$tomato_count" = 63 ] || die "expected 63 tomato CRAMs in $tomato_crams, found $tomato_count"

in=$work/in
mkdir -p "$in/tomato" "$in/hg002" "$work/bin"
for f in "$tomato_fa" "$tomato_fa.fai" "$tomato_crams"/*.bench.cram "$tomato_crams"/*.bench.cram.crai; do
  ln -sf "$f" "$in/tomato/"
done
for f in "$human_fa" "$human_fa.fai" "$hg002_cram" "$hg002_cram.crai"; do
  ln -sf "$f" "$in/hg002/"
done
# The regions are copied, not linked, so a later edit to the benchmark's BED shows as a different
# slice (the baseline records each file's checksum) rather than as moved outputs.
head -n 1 "$tomato_bed" > "$in/tomato/regions1.bed"
head -n 2 "$tomato_bed" > "$in/tomato/regions2.bed"
cp "$tomato_bed" "$in/tomato/regions80.bed"
head -n 50 "$hg002_bed" > "$in/hg002/regions50.bed"

# ---------------------------------------------------------------------------------------------
# The binary
# ---------------------------------------------------------------------------------------------
if [ -n "${CQ_BIN:-}" ]; then
  bin=$CQ_BIN
else
  echo "=== building pvc (release, target-host) ==="
  CARGO_TARGET_DIR="$root/target-host" cargo build --release --bin pvc --quiet \
    --manifest-path "$root/Cargo.toml" || die "the build failed"
  bin=$root/target-host/release/pvc
fi
[ -x "$bin" ] || die "no binary at $bin"
# Called as `pvc` through a link, so argv[0] — which the VCF and psp headers record — is the same
# whichever binary this is.
ln -sf "$bin" "$work/bin/pvc"
export PATH="$work/bin:$PATH"

# ---------------------------------------------------------------------------------------------
# The steps of one pass
# ---------------------------------------------------------------------------------------------
run=$work/run
measure=$work/measure.tsv

if /usr/bin/time -l true >/dev/null 2>&1; then
  time_flag=-l
elif /usr/bin/time -v true >/dev/null 2>&1; then
  time_flag=-v
else
  die "no /usr/bin/time that reports memory (BSD -l or GNU -v)"
fi

# Wall seconds, peak resident MB and instructions retired (0 where the platform does not count
# them), from one /usr/bin/time report.
parse_time() {
  if [ "$time_flag" = -l ]; then
    awk '/ real / && !w { w = $1 }
         /maximum resident set size/ { r = $1 / 1048576 }
         /instructions retired/ { i = $1 }
         END { printf "%.3f\t%.1f\t%.0f\n", w, r, i + 0 }' "$1"
  else
    awk -F': ' '/Elapsed \(wall clock\)/ { n = split($2, p, ":"); w = 0
                                          for (k = 1; k <= n; k++) w = w * 60 + p[k] }
                /Maximum resident set size/ { r = $2 / 1024 }
                END { printf "%.3f\t%.1f\t0\n", w, r }' "$1"
  fi
}

# step NAME COMMAND... — run one pvc process from inside $run, measure it, record the figures.
step() {
  local name=$1; shift
  local log=$run/logs/$name
  if ! /usr/bin/time "$time_flag" -o "$log.time" "$@" > "$log.out" 2>&1; then
    echo "step $name failed; its output, last lines ($log.out):" >&2
    tail -n 15 "$log.out" >&2
    exit 2
  fi
  local figures
  figures=$(parse_time "$log.time")
  printf '%s\t%s\t%s\n' "$rep" "$name" "$figures" >> "$measure"
  printf '  %-34s %s\n' "$name" "$(echo "$figures" | awk -F'\t' '{printf "%7.2f s %8.1f MB", $1, $2}')"
}

T=../in/tomato
H=../in/hg002
tomato_ref="--reference $T/S_lycopersicum_chromosomes.4.00.fa --catalog catalog/S_lycopersicum_chromosomes.4.00.repeats.parquet"
human_ref="--reference $H/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna --catalog catalog/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.repeats.parquet"

# The 16 and the 63 accessions as --alignment arguments, in file-name order.
tomato16_alignments=""
for cram in $(cd "$in/tomato" && ls *.bench.cram | head -n 16); do
  tomato16_alignments="$tomato16_alignments --alignment $T/$cram"
done
tomato63_alignments=""
for cram in $(cd "$in/tomato" && ls *.bench.cram); do
  tomato63_alignments="$tomato63_alignments --alignment $T/$cram"
done

# A psp's header: the 4-byte magic, the TOML body's length as a little-endian u64, the body. The
# body's text, for reading keys out of it.
psp_header() {
  local length
  length=$(od -An -t u8 -j 4 -N 8 "$1" | tr -d ' ')
  head -c $((12 + length)) "$1" | tail -c +13
}

# The tomato16 psps as --psp arguments **in the order their alignments were given**, read from each
# header's `input-alignments`, so call-from-psps writes its sample columns in call-from-alignments'
# order and the two VCFs can be compared whole.
tomato16_psps() {
  local cram psp pairs=$run/logs/psp_alignments.txt
  : > "$pairs"
  for psp in tomato16/psps/*.psp; do
    echo "$(psp_header "$psp" | sed -n 's/^input-alignments = \["\(.*\)"\]$/\1/p') $psp" >> "$pairs"
  done
  for cram in $(cd "$in/tomato" && ls *.bench.cram | head -n 16); do
    awk -v c="$cram" '$1 == c { printf " --psp %s", $2 }' "$pairs"
  done
}

# shellcheck disable=SC2086  # the argument strings above are split on purpose
one_pass() {
  mkdir -p "$run/logs" "$run/catalog" "$run/tomato16" "$run/tomato63" "$run/tomato1" "$run/hg002"
  cd "$run"
  step catalog-tomato pvc repeat-catalog --reference $T/S_lycopersicum_chromosomes.4.00.fa \
    --output catalog/S_lycopersicum_chromosomes.4.00.repeats.parquet \
    --threads $threads --force --quiet --log logs/catalog-tomato.pvc.log
  step catalog-grch38 pvc repeat-catalog \
    --reference $H/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --output catalog/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.repeats.parquet \
    --threads $threads --force --quiet --log logs/catalog-grch38.pvc.log

  step tomato16-generate-psps pvc generate-psps $tomato_ref $tomato16_alignments \
    --regions $T/regions1.bed --output-dir tomato16/psps --force --quiet
  local psps
  psps=$(tomato16_psps)
  [ "$(echo $psps | wc -w | tr -d ' ')" = 32 ] || die "expected 16 tomato16 psps, found: $psps"
  step tomato16-estimate-parameters pvc estimate-parameters $tomato_ref $psps \
    --output tomato16/parameters.toml --threads $threads --force --quiet \
    --log logs/tomato16-estimate-parameters.pvc.log
  NG_WINDOW_COVERAGE_FILE=tomato16/from_psps.windows \
  step tomato16-call-from-psps pvc call-from-psps $tomato_ref $psps \
    --parameters tomato16/parameters.toml --threads $threads --output tomato16/from_psps.vcf \
    --force --quiet --log logs/tomato16-call-from-psps.pvc.log
  NG_WINDOW_COVERAGE_FILE=tomato16/from_alignments.windows \
  step tomato16-call-from-alignments pvc call-from-alignments $tomato_ref $tomato16_alignments \
    --regions $T/regions1.bed --parameters tomato16/parameters.toml --threads $threads \
    --output tomato16/from_alignments.vcf --force --quiet \
    --log logs/tomato16-call-from-alignments.pvc.log
  step tomato16-call-refined pvc call-from-psps $tomato_ref $psps \
    --parameters tomato16/parameters.toml --refine-genotypes on --threads $threads \
    --output tomato16/refined.vcf --force --quiet --log logs/tomato16-call-refined.pvc.log

  NG_WINDOW_COVERAGE_FILE=tomato63/from_alignments.windows \
  step tomato63-call-from-alignments pvc call-from-alignments $tomato_ref $tomato63_alignments \
    --regions $T/regions2.bed --defaults --threads $threads \
    --output tomato63/from_alignments.vcf --force --quiet \
    --log logs/tomato63-call-from-alignments.pvc.log
  step tomato1-call-from-alignments pvc call-from-alignments $tomato_ref \
    --alignment $T/SRR7279481.p1.bench.cram --regions $T/regions80.bed --defaults \
    --threads $threads --output tomato1/from_alignments.vcf --force --quiet \
    --log logs/tomato1-call-from-alignments.pvc.log

  step hg002-generate-psps pvc generate-psps $human_ref \
    --alignment $H/HG002_reads_selected_1000_rg.cram --regions $H/regions50.bed \
    --output-dir hg002/psps --force --quiet
  step hg002-estimate-parameters pvc estimate-parameters $human_ref --psp hg002/psps/HG002.psp \
    --output hg002/parameters.toml --threads $threads --force --quiet \
    --log logs/hg002-estimate-parameters.pvc.log
  NG_WINDOW_COVERAGE_FILE=hg002/from_psps.windows \
  step hg002-call-from-psps pvc call-from-psps $human_ref --psp hg002/psps/HG002.psp \
    --parameters hg002/parameters.toml --threads $threads --output hg002/from_psps.vcf \
    --force --quiet --log logs/hg002-call-from-psps.pvc.log
  NG_WINDOW_COVERAGE_FILE=hg002/from_alignments.windows \
  step hg002-call-from-alignments pvc call-from-alignments $human_ref \
    --alignment $H/HG002_reads_selected_1000_rg.cram --regions $H/regions50.bed \
    --parameters hg002/parameters.toml --threads $threads --output hg002/from_alignments.vcf \
    --force --quiet --log logs/hg002-call-from-alignments.pvc.log
  cd "$root"
}

# ---------------------------------------------------------------------------------------------
# Normalising and checksumming the outputs
# ---------------------------------------------------------------------------------------------

# normalise KIND FILE — the bytes that are compared, on stdout.
normalise() {
  case $1 in
    vcf) grep -v '^##commandline=\|^##pvcLog=' "$2" ;;
    psp)
      local length
      length=$(od -An -t u8 -j 4 -N 8 "$2" | tr -d ' ')
      head -c 12 "$2"
      psp_header "$2" | grep -v '^created = \|^command-line = '
      tail -c +$((13 + length)) "$2"
      ;;
    windows) LC_ALL=C sort "$2" ;;
    file) cat "$2" ;;
  esac
}

# Every compared output of a pass, as `KIND PATH` relative to $run.
artefacts() {
  echo "file catalog/S_lycopersicum_chromosomes.4.00.repeats.parquet"
  echo "file catalog/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.repeats.parquet"
  local psp
  for psp in $(cd "$run" && ls tomato16/psps/*.psp); do echo "psp $psp"; done
  echo "file tomato16/parameters.toml"
  for side in tomato16/from_psps tomato16/from_alignments tomato63/from_alignments \
              hg002/from_psps hg002/from_alignments; do
    echo "vcf $side.vcf"
    echo "windows $side.windows"
    echo "file $side.windows.histograms"
  done
  echo "vcf tomato16/refined.vcf"
  echo "vcf tomato1/from_alignments.vcf"
  echo "psp hg002/psps/HG002.psp"
  echo "file hg002/parameters.toml"
}

# digests OUT [KEEP_DIR] — one `digest NAME MD5` line per output, and the slice the run walked; with
# KEEP_DIR, a copy of each normalised output too.
digests() {
  local out=$1 keep=${2:-} kind path
  : > "$out"
  for bed in tomato/regions1.bed tomato/regions2.bed tomato/regions80.bed hg002/regions50.bed; do
    echo "digest in/$bed $(md5_stdin < "$in/$bed")" >> "$out"
  done
  artefacts | while read -r kind path; do
    [ -s "$run/$path" ] || { echo "digest $path MISSING" >> "$out"; continue; }
    if [ -n "$keep" ]; then
      mkdir -p "$keep/$(dirname "$path")"
      normalise "$kind" "$run/$path" > "$keep/$path"
      echo "digest $path $(md5_stdin < "$keep/$path")" >> "$out"
    else
      echo "digest $path $(normalise "$kind" "$run/$path" | md5_stdin)" >> "$out"
    fi
  done
}

# What changed in one output, read against the baseline's copy of it when this worktree has one.
explain_difference() {
  local path=$1 old=$work/outputs/baseline/$1 new=$work/outputs/current/$1
  if [ ! -f "$old" ]; then
    echo "      (no copy of the baseline's output in $work/outputs/baseline to diff against)"
    return
  fi
  case $path in
    *.vcf)
      local old_n new_n
      old_n=$(grep -vc '^#' "$old" || true); new_n=$(grep -vc '^#' "$new" || true)
      echo "      records: $old_n in the baseline, $new_n now;" \
           "$(diff "$old" "$new" | grep -c '^[<>]' || true) lines differ. First differences (< baseline, > now):"
      diff "$old" "$new" | grep '^[<>]' | head -n 6 | cut -c1-200 | sed 's/^/        /'
      ;;
    *.psp)
      if cmp -s <(psp_header "$old") <(psp_header "$new"); then
        echo "      the header is the same; the records after it differ"
      else
        echo "      the header differs:"
        diff <(psp_header "$old") <(psp_header "$new") | grep '^[<>]' | head -n 6 | cut -c1-200 \
          | sed 's/^/        /'
      fi
      ;;
    *.parquet)
      echo "      the catalog differs (sizes $(wc -c < "$old" | tr -d ' ') and $(wc -c < "$new" | tr -d ' ') bytes)"
      ;;
    *)
      diff "$old" "$new" | grep '^[<>]' | head -n 6 | cut -c1-200 | sed 's/^/        /'
      ;;
  esac
}

# ---------------------------------------------------------------------------------------------
# The passes
# ---------------------------------------------------------------------------------------------
: > "$measure"
host=$(hostname -s)
cpu=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)
status=0
echo "binary: $bin"
echo "inputs: $data, $genomes"
started=$(date +%s)

for rep in $(seq 1 "$reps"); do
  # Each pass starts empty, so a step that stopped writing an output cannot pass on a stale one.
  rm -rf "$run"
  load=$(sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}' || cut -d' ' -f1 /proc/loadavg)
  echo "=== pass $rep of $reps (1-minute load $load) ==="
  one_pass
  if [ "$rep" = 1 ]; then
    rm -rf "$work/outputs/current"
    digests "$work/digests.txt" "$work/outputs/current"
    # The two calling routes, handed one parameters file, must write the same VCF.
    for cohort in tomato16 hg002; do
      for what in vcf windows windows.histograms; do
        if ! cmp -s "$work/outputs/current/$cohort/from_psps.$what" \
                    "$work/outputs/current/$cohort/from_alignments.$what"; then
          echo "DIFFERS: $cohort's call-from-psps and call-from-alignments wrote different .$what"
          status=1
        fi
      done
    done
    if [ "$mode" != record ]; then
      [ -r "$baseline" ] || die "no baseline at $baseline; record one with --record"
      changed=0
      while read -r _ path sum; do
        expected=$(awk -v p="$path" '$1 == "digest" && $2 == p { print $3 }' "$baseline")
        if [ -z "$expected" ]; then
          echo "DIFFERS: $path is new (not in the baseline)"; changed=1
        elif [ "$expected" != "$sum" ]; then
          echo "DIFFERS: $path"; explain_difference "$path" || true; changed=1
        fi
      done < "$work/digests.txt"
      while read -r kind path _; do
        [ "$kind" = digest ] || continue
        grep -q "^digest $path " "$work/digests.txt" \
          || { echo "DIFFERS: $path is in the baseline but was not produced"; changed=1; }
      done < "$baseline"
      if [ "$changed" = 0 ]; then
        echo "outputs: all $(wc -l < "$work/digests.txt" | tr -d ' ') checksums match the baseline"
      else
        status=1
      fi
    fi
  else
    digests "$work/digests.again.txt"
    if ! cmp -s "$work/digests.txt" "$work/digests.again.txt"; then
      echo "DIFFERS between passes 1 and $rep of the same build (the output is not deterministic):"
      diff "$work/digests.txt" "$work/digests.again.txt" | grep '^>' | sed 's/^> digest /  /' || true
      status=1
    fi
  fi
done
elapsed=$(( $(date +%s) - started ))

# ---------------------------------------------------------------------------------------------
# Timing and memory
# ---------------------------------------------------------------------------------------------

# Per step (and TOTAL, the per-pass sum of wall time): min, median and max of wall, memory and
# instructions over the passes made.
summarise() {
  awk -F'\t' '
    { order[$2] = order[$2] ? order[$2] : ++steps; name[order[$2]] = $2
      n[$2]++; w[$2, n[$2]] = $3; r[$2, n[$2]] = $4; i[$2, n[$2]] = $5
      total[$1] += $3; reps[$1] = 1 }
    function stats(arr, s, k,   a, m, x, y, t) {
      m = n[s]
      for (x = 1; x <= m; x++) a[x] = arr[s, x]
      for (x = 2; x <= m; x++) { t = a[x]; for (y = x - 1; y >= 1 && a[y] > t; y--) a[y + 1] = a[y]; a[y + 1] = t }
      return sprintf("%s %s %s", a[1], (m % 2 ? a[(m + 1) / 2] : (a[m / 2] + a[m / 2 + 1]) / 2), a[m])
    }
    END {
      for (s = 1; s <= steps; s++) {
        k = name[s]
        printf "perf %s %s %s %s\n", k, stats(w, k), stats(r, k), stats(i, k)
      }
      m = 0; for (p in reps) { m++; tw["TOTAL", m] = total[p] }
      n["TOTAL"] = m
      printf "perf TOTAL %s 0 0 0 0 0 0\n", stats(tw, "TOTAL")
    }' "$measure"
}
summarise > "$work/perf.txt"

if [ "$mode" = record ]; then
  if [ "$status" != 0 ]; then
    echo "not recording: the passes disagreed with each other or the two calling routes differ" >&2
    exit 1
  fi
  {
    echo "# Baseline for scripts/code_quality_oracle.sh. Re-record with --record only when an output is meant"
    echo "# to change (say in the commit which lines moved and why) or the timing moved for a reason outside"
    echo "# the code. Recorded $(date -u +%Y-%m-%dT%H:%MZ) from $(git -C "$root" rev-parse --short HEAD)$(git -C "$root" diff --quiet HEAD -- src Cargo.toml Cargo.lock || echo ' (with uncommitted changes)'), $reps passes."
    echo "#"
    echo "# perf lines: step, then min median max of wall seconds, of peak resident MB, and of instructions"
    echo "# retired. They are compared only on the host and thread count below."
    echo "host $(echo "$host $cpu" | tr ' ' '_')"
    echo "threads $threads"
    cat "$work/digests.txt"
    cat "$work/perf.txt"
  } > "$baseline"
  rm -rf "$work/outputs/baseline"
  cp -R "$work/outputs/current" "$work/outputs/baseline"
  echo
  echo "recorded $baseline ($reps passes, ${elapsed}s); per step, min / median / max:"
  awk '$1 == "perf" { printf "  %-34s wall %s / %s / %s s   memory %s / %s / %s MB\n", $2, $3, $4, $5, $6, $7, $8 }' "$baseline"
  exit 0
fi

if [ "$mode" = identity ]; then
  echo "identity only: ${elapsed}s; timing not compared"
  exit "$status"
fi

base_host=$(awk '$1 == "host" { print $2 }' "$baseline")
base_threads=$(awk '$1 == "threads" { print $2 }' "$baseline")
this_host=$(echo "$host $cpu" | tr ' ' '_')
perf_status=0
echo
if [ "$base_host" != "$this_host" ] || [ "$base_threads" != "$threads" ]; then
  echo "timing NOT compared: the baseline was recorded on $base_host at $base_threads threads," \
       "this is $this_host at $threads"
else
  echo "timing and memory, fastest of $reps passes against the baseline's fastest" \
       "(tolerance: 2% or the baseline's own spread, whichever is larger):"
  printf '  %-30s %20s %22s %26s\n' step "wall s (tol)" "memory MB (tol)" "instructions G (tol)"
  awk -v report="$work/perf.verdict" '
    NR == FNR { if ($1 == "perf") { for (f = 3; f <= 11; f++) base[$2, f] = $f }; next }
    $1 == "perf" {
      flagged = ""
      line = sprintf("  %-30s", $2)
      for (m = 0; m < 3; m++) {
        f = 3 + 3 * m
        b = base[$2, f]; bmax = base[$2, f + 2]; now = $f
        if (b == "" ) { line = line sprintf(" %26s", "(new step)"); continue }
        if (b + 0 == 0) { line = line sprintf(" %26s", "-"); continue }
        tol = (bmax - b) / b; if (tol < 0.02) tol = 0.02
        change = now / b - 1
        mark = change > tol ? "  SLOWER" : (change < -tol ? "  better" : "")
        if (m == 1 && change > tol) mark = "  LARGER"
        if (change > tol) flagged = flagged " " (m == 0 ? "wall" : (m == 1 ? "memory" : "instructions"))
        scale = (m == 2 ? 1e9 : 1)
        line = line sprintf(" %8.2f %+6.1f%% (%4.1f%%)%s", now / scale, 100 * change, 100 * tol, mark)
      }
      print line
      if (flagged != "") print $2 ":" flagged > report
    }' "$baseline" "$work/perf.txt"
  if [ -s "$work/perf.verdict" ]; then perf_status=3; fi
  rm -f "$work/perf.verdict"
fi

echo
if [ "$status" != 0 ]; then
  echo "VERDICT: OUTPUT DIFFERS (exit 1) — ${elapsed}s"
  exit 1
elif [ "$perf_status" != 0 ]; then
  echo "VERDICT: outputs identical, but slower or larger beyond noise (exit 3) — ${elapsed}s"
  exit 3
fi
echo "VERDICT: identical, and within noise — ${elapsed}s"
exit 0
