#!/bin/bash
# **Run-level invariances on real reads**: the same cohort called in ways that must not change the
# calls, compared byte for byte. Spec doc/devel/ng/spec/run_streaming.md §12 names most of these
# and the unit tests hold them on one-contig fixtures; this runs them where real reads, a real
# catalog and a fitted parameters file give a defect somewhere to show up.
#
#   fit       estimate-parameters at 1 thread and at --threads, and with the psps reversed:
#             the same parameters file
#   psp       call-from-psps at 1 thread, at --threads, with the psps reversed, and with a
#             narrow locus-building region: the same VCF
#   direct    call-from-alignments at 1 thread, at --threads, and with the alignments reversed:
#             the same VCF; and with the fitted parameters, psp mode's VCF
#   regions   call-from-alignments over the whole BED, and once per BED line with the results
#             joined: the same records. Only with the hidden-dup filter off (--hidden-dup-fdr 0),
#             whose false-discovery rate is a genome-wide quantity and so legitimately depends on
#             what else was called.
#
# A reversed cohort writes its sample columns reversed, so every VCF is normalised before the
# comparison: its sample columns sorted by name, and the two header lines that record the command
# and the log path dropped. Nothing else is touched.
#
# Every check runs and reports; the script exits 1 at the end if any differed.
#
# Usage:
#
#   scripts/run_invariance_oracle.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <output-dir> <alignment.cram>...
#
# PVC_BIN picks the binary; without it, the newest of target{,-container}/{soak,release}/pvc.
# The soak build is preferred because it carries debug assertions and overflow checks, so every
# run here doubles as a run under those. THREADS (default 6) sets the many-thread side.
set -eu

if [ "$#" -lt 5 ]; then
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

reference=$1; catalog=$2; regions=$3; out=$4
shift 4
alignments=("$@")
threads=${THREADS:-6}

root=$(cd "$(dirname "$0")/.." && pwd)
bin=${PVC_BIN:-}
if [ -z "$bin" ]; then
  for candidate in "$root"/target{,-container}/{soak,release}/pvc; do
    if [ -x "$candidate" ] && { [ -z "$bin" ] || [ "$candidate" -nt "$bin" ]; }; then
      bin=$candidate
    fi
  done
fi
if [ -z "$bin" ]; then
  echo "no build of pvc; build one first (cargo build --profile soak --bin pvc)" >&2
  exit 1
fi
echo "binary: $bin"
mkdir -p "$out"

failures=0

# Sample columns sorted by name, `##commandline` and `##pvcLog` dropped. awk on macOS has no
# asort, hence the insertion sort over the header's names.
normalise_vcf() {
  awk -F'\t' -v OFS='\t' '
    /^##commandline=|^##pvcLog=/ { next }
    /^##/ { print; next }
    /^#CHROM/ {
      n = 0
      for (i = 10; i <= NF; i++) { n++; name[n] = $i; col[n] = i }
      for (i = 2; i <= n; i++) {
        for (j = i; j > 1 && name[j - 1] > name[j]; j--) {
          t = name[j]; name[j] = name[j - 1]; name[j - 1] = t
          t = col[j]; col[j] = col[j - 1]; col[j - 1] = t
        }
      }
    }
    {
      line = $1
      for (i = 2; i <= 9 && i <= NF; i++) line = line OFS $i
      for (k = 1; k <= n; k++) line = line OFS $(col[k])
      print line
    }' "$1"
}

# compare <label> <reference-file> <other-file> [vcf]
compare() {
  local label=$1 a=$2 b=$3 kind=${4:-}
  if [ "$kind" = vcf ]; then
    normalise_vcf "$a" > "$a.normalised"; a=$a.normalised
    normalise_vcf "$b" > "$b.normalised"; b=$b.normalised
  fi
  if cmp -s "$a" "$b"; then
    echo "SAME       $label ($(grep -vc '^#' "$a" || true) non-header lines)"
  else
    echo "DIFFERENT  $label"
    diff "$a" "$b" | head -20 | sed 's/^/    /'
    echo "    ($(diff "$a" "$b" | grep -c '^[<>]') differing lines)"
    failures=$((failures + 1))
  fi
}

# The parameters file a run scored with, from between its log's marker lines.
parameters_from_log() {
  sed -n '/^# ---- parameters file begins ----$/,/^# ---- parameters file ends ----$/p' "$1" \
    | sed '1d;$d'
}

reversed=()
for ((i = ${#alignments[@]} - 1; i >= 0; i--)); do reversed+=("${alignments[i]}"); done
as_args() { local flag=$1; shift; for f in "$@"; do printf '%s\n%s\n' "$flag" "$f"; done; }

echo "=== walking the alignments into psps ==="
psps=$out/psps
mkdir -p "$psps"
# shellcheck disable=SC2046
"$bin" generate-psps --reference "$reference" --catalog "$catalog" \
  $(as_args --alignment "${alignments[@]}") --regions "$regions" \
  --output-dir "$psps" --force > "$out/generate.log" 2>&1
psp_files=("$psps"/*.psp)
psp_reversed=()
for ((i = ${#psp_files[@]} - 1; i >= 0; i--)); do psp_reversed+=("${psp_files[i]}"); done

echo "=== fit ==="
fit() {
  local name=$1 nthreads=$2; shift 2
  "$bin" estimate-parameters --reference "$reference" --catalog "$catalog" \
    $(as_args --psp "$@") --threads "$nthreads" \
    --output "$out/$name.parameters.toml" --log "$out/$name.fit.log" --force \
    > "$out/$name.fit.stdout" 2>&1
}
fit fit_t1 1 "${psp_files[@]}"
fit fit_tn "$threads" "${psp_files[@]}"
fit fit_reversed "$threads" "${psp_reversed[@]}"
compare "fit: 1 thread vs $threads" "$out/fit_t1.parameters.toml" "$out/fit_tn.parameters.toml"
compare "fit: psps reversed" "$out/fit_tn.parameters.toml" "$out/fit_reversed.parameters.toml"
parameters=$out/fit_tn.parameters.toml

echo "=== psp mode ==="
call_psps() {
  local name=$1 nthreads=$2 extra=$3; shift 3
  # shellcheck disable=SC2086
  "$bin" call-from-psps --reference "$reference" --catalog "$catalog" \
    $(as_args --psp "$@") --parameters "$parameters" --threads "$nthreads" $extra \
    --output "$out/$name.vcf" --log "$out/$name.log" --force > "$out/$name.stdout" 2>&1
}
call_psps psp_t1 1 "" "${psp_files[@]}"
call_psps psp_tn "$threads" "" "${psp_files[@]}"
call_psps psp_reversed "$threads" "" "${psp_reversed[@]}"
call_psps psp_narrow "$threads" "--locus-builder-regions-len 7919" "${psp_files[@]}"
compare "psp: 1 thread vs $threads" "$out/psp_t1.vcf" "$out/psp_tn.vcf" vcf
compare "psp: psps reversed" "$out/psp_tn.vcf" "$out/psp_reversed.vcf" vcf
compare "psp: 7,919-base locus-building regions" "$out/psp_tn.vcf" "$out/psp_narrow.vcf" vcf

echo "=== direct mode ==="
call_direct() {
  local name=$1 nthreads=$2 extra=$3; shift 3
  # shellcheck disable=SC2086
  "$bin" call-from-alignments --reference "$reference" --catalog "$catalog" \
    $(as_args --alignment "$@") --threads "$nthreads" $extra \
    --output "$out/$name.vcf" --log "$out/$name.log" --force > "$out/$name.stdout" 2>&1
}
call_direct direct_t1 1 "--defaults --regions $regions" "${alignments[@]}"
call_direct direct_tn "$threads" "--defaults --regions $regions" "${alignments[@]}"
call_direct direct_reversed "$threads" "--defaults --regions $regions" "${reversed[@]}"
call_direct direct_fitted "$threads" "--parameters $parameters --regions $regions" \
  "${alignments[@]}"
compare "direct: 1 thread vs $threads" "$out/direct_t1.vcf" "$out/direct_tn.vcf" vcf
compare "direct: alignments reversed" "$out/direct_tn.vcf" "$out/direct_reversed.vcf" vcf
compare "direct vs psp mode, fitted parameters" "$out/psp_tn.vcf" "$out/direct_fitted.vcf" vcf

echo "=== regions ==="
call_direct regions_whole "$threads" "--defaults --hidden-dup-fdr 0 --regions $regions" \
  "${alignments[@]}"
mkdir -p "$out/regions_split"
line=0
while read -r contig start end _; do
  line=$((line + 1))
  printf '%s\t%s\t%s\n' "$contig" "$start" "$end" > "$out/regions_split/$line.bed"
  call_direct "regions_split/$line" "$threads" \
    "--defaults --hidden-dup-fdr 0 --regions $out/regions_split/$line.bed" "${alignments[@]}"
done < <(grep -v '^#\|^track\|^browser' "$regions")
# The joined file takes the header of the first piece; the records are each piece's, put in the
# order a whole run writes them — contigs as the header's `##contig` lines list them, then by
# position. Not BED order: a BED sorted by name puts chr10 before chr2, and a whole run does not.
{
  grep '^#' "$out/regions_split/1.vcf"
  for ((i = 1; i <= line; i++)); do grep -v '^#' "$out/regions_split/$i.vcf" || true; done |
    awk -F'\t' -v OFS='\t' -v header="$out/regions_split/1.vcf" '
      BEGIN {
        while ((getline row < header) > 0) {
          if (row ~ /^##contig=<ID=/) {
            id = row; sub(/^##contig=<ID=/, "", id); sub(/[,>].*$/, "", id)
            rank[id] = ++contigs
          }
        }
      }
      { print rank[$1], $2, $0 }' |
    sort -s -t$'\t' -k1,1n -k2,2n | cut -f3-
} > "$out/regions_joined.vcf"
compare "regions: whole BED vs one run per line, joined" \
  "$out/regions_whole.vcf" "$out/regions_joined.vcf" vcf

echo
if [ "$failures" -eq 0 ]; then
  echo "every invariance held"
else
  echo "$failures invariance(s) broken"
  exit 1
fi
