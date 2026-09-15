#!/bin/sh
# **The identity oracle for the ng promotion** (doc/devel/implementation_plans/
# promote_ng_to_production.md, milestone A): one cohort of real reads, called by whatever
# binary is built right now, reduced to a short list of checksums that a later build must
# reproduce exactly.
#
# It is not a new comparison. `ng_mode_equivalence_oracle.sh` already calls one cohort two
# ways and proves the two routes agree; this wraps it and writes down *what they agreed on*,
# so the question can be asked across builds instead of only within one run. The promotion
# changes no arithmetic anywhere, so every checksum below must survive every milestone —
# and the one line that cannot is named in `##commandline`, below.
#
# Usage:
#
#   scripts/promote_ng_oracle.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <output-dir> <alignment.cram>...
#
# Run it inside the dev container, because the binary it finds was built there:
#
#   DEV_EXTRA_MOUNT=<dir holding the crams and catalog> \
#     ./scripts/dev.sh scripts/promote_ng_oracle.sh ...
#
# **The PROMOTE_NG_* settings below do not cross into the container on their own**: `dev.sh`
# forwards only `CARGO_*`, `NG_*` and a few `RUST*` names. Set them inside the container command:
#
#   ./scripts/dev.sh sh -c 'PROMOTE_NG_BASELINE=none scripts/promote_ng_oracle.sh ...'
#
# **PROMOTE_NG_REGIONS caps how much genome is walked**, taking that many leading lines of
# the regions file (default 20). Identity is the question, not coverage: a byte that moves
# moves in the first region as readily as the hundredth, and a cheap oracle is one that gets
# run at every checkpoint rather than once. The count is written into the digest file, so two
# runs over different slices cannot be compared by accident.
#
# **The run is also fitted, and called with its fit** (added 2026-09-15,
# doc/devel/implementation_plans/portable_float.md, step C3). The two calling routes above run
# with `--defaults`, and default parameters almost never let one unit in the last binary place
# reach a call: before the maths went through `libm`, macOS and Linux wrote identical VCFs over
# 160 regions with default parameters, while their fitted parameters differed in 7 of 574 lines.
# So the fit's file, and the VCF called with it, are pinned too. The fit costs minutes — about 7 on
# the macOS host and 10 in the Linux container over the default 20 regions — against seconds for
# the rest; PROMOTE_NG_FIT=0 (exactly `0`) skips it, and the digest file says it was skipped.
#
# **The digests are compared with scripts/promote_ng_oracle.baseline**, the checksums recorded
# the last time the output was meant to change, and the script exits 1 if any differs or the
# baseline cannot be read. The comparison includes the slice walked — regions kept and alignment
# count — so a run over a different slice reports that rather than seven moved checksums. Set
# PROMOTE_NG_BASELINE to another file to compare against it, or to `none` to only record.
set -eu

if [ "$#" -lt 5 ]; then
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

reference=$1; catalog=$2; regions=$3; out=$4
shift 4

root=$(cd "$(dirname "$0")/.." && pwd)
regions_kept=${PROMOTE_NG_REGIONS:-20}

mkdir -p "$out"
head -n "$regions_kept" "$regions" > "$out/regions.bed"
kept=$(wc -l < "$out/regions.bed" | tr -d ' ')

# The wrapped oracle does the work: both routes, and its own comparison between them. A
# failure there is a failure here — the two routes disagreeing makes any checksum of either
# one meaningless — so its exit status is not swallowed.
"$root/scripts/ng_mode_equivalence_oracle.sh" \
  "$reference" "$catalog" "$out/regions.bed" "$out" "$@"

# **The fit, and the calls made with it**, over the psps the wrapped oracle stored. The binary is
# found the way the wrapped oracle finds it — the newer of the two release builds — so both halves
# pin one build; unlike there, a build that does not run on this machine (a macOS binary seen from
# the Linux container) is refused with a message rather than exec'd.
fit=${PROMOTE_NG_FIT:-1}
if [ "$fit" != 0 ]; then
  bin=""
  for candidate in "$root/target-container/release/pop_var_caller" \
                   "$root/target/release/pop_var_caller"; do
    if [ -x "$candidate" ] && { [ -z "$bin" ] || [ "$candidate" -nt "$bin" ]; }; then
      bin=$candidate
    fi
  done
  if ! "$bin" --version > /dev/null 2>&1; then
    echo "the newest release build, $bin, does not run here; touch the right one and rerun" >&2
    exit 1
  fi
  # The psps in the VCF's sample order, as the wrapped oracle names them.
  stored=""
  for sample in $(grep '^#CHROM' "$out/from_alignments.vcf" | cut -f10-); do
    stored="$stored --psp $out/psps/$sample.psp"
  done
  echo "=== the fit ==="
  # shellcheck disable=SC2086
  "$bin" estimate-parameters --reference "$reference" --catalog "$catalog" $stored \
    --output "$out/fitted.parameters.toml" --force > "$out/fit.log" 2>&1 \
    || { echo "the fit failed; see $out/fit.log" >&2; exit 1; }
  echo "=== psp mode, with the fit ==="
  # shellcheck disable=SC2086
  "$bin" call-from-psps --reference "$reference" --catalog "$catalog" $stored \
    --parameters "$out/fitted.parameters.toml" --threads 4 \
    --output "$out/from_fit.vcf" > "$out/from_fit.log" 2>&1 \
    || { echo "calling with the fit failed; see $out/from_fit.log" >&2; exit 1; }
  grep -v '^##commandline=' "$out/from_fit.vcf" > "$out/from_fit.comparable"
fi

# **What gets a checksum, and why these five.**
#
# `*.comparable` is the VCF with `##commandline` removed — every header line and every
# record of every call, and the one line the promotion is *expected* to change, since
# milestone E renames the binary that writes it.
#
# The parameters file is what the run fitted before calling: error rates, genotype
# frequencies, the inbreeding coefficient and its warrant. A change there moves every call
# that follows, so it is pinned beside the calls themselves.
#
# The window rows and the per-sample histograms are the hidden-duplication filter's input.
# They travel beside the records rather than in them, so a run that lost them writes exactly
# the VCF a correct run writes (spec window_coverage.md §10) — which is why the wrapped
# oracle checks them, and why they are pinned here too.
#
# **The psps are deliberately absent.** A psp's header carries the time it was written, so
# its file checksum differs between two runs of one binary and cannot serve as a baseline;
# skipping the header needs its byte length, which only the reader knows. Nothing is given
# up: a psp reaches this comparison through `call-from-psps`, so a change to what a psp
# carries shows in `from_psps.comparable`, and a change that does not reach the VCF is not a
# change to anything the caller decides.
{
  echo "# promote-ng identity oracle"
  echo "# reference: $(basename "$reference")"
  echo "# catalog: $(basename "$catalog")"
  echo "# regions kept: $kept (of $(wc -l < "$regions" | tr -d ' ') in $regions)"
  echo "# alignments: $#"
  for cram in "$@"; do echo "#   $(basename "$cram")"; done
  # The one descriptive line the baseline comparison reads: which slice was walked.
  echo "slice: $kept regions, $# alignments"
  # `|| true`: `grep -c` exits 1 on a count of zero, which `set -e` would turn into an abort
  # here rather than into the empty-output failure the wrapped oracle already reports.
  echo "# records: $(grep -vc '^#' "$out/from_alignments.vcf" || true)"
  artefacts="from_alignments.comparable from_psps.comparable from_alignments.parameters.toml
             from_alignments.windows.sorted from_alignments.windows.histograms"
  if [ "$fit" != 0 ]; then
    artefacts="$artefacts fitted.parameters.toml from_fit.comparable"
  else
    echo "# fit: skipped (PROMOTE_NG_FIT=0)"
  fi
  for artefact in $artefacts; do
    # `md5sum` on Linux, `md5 -r` on macOS; the container is Linux but the script is also
    # run on a host checkout, so neither is assumed.
    if command -v md5sum >/dev/null 2>&1; then
      sum=$(md5sum "$out/$artefact" | cut -d' ' -f1)
    else
      sum=$(md5 -q "$out/$artefact")
    fi
    echo "$sum  $artefact"
  done
} > "$out/digests.txt"

echo
echo "=== digests written to $out/digests.txt ==="
cat "$out/digests.txt"

# **Against the recorded baseline**, line by line: the slice walked, then every checksum. The `#`
# lines describe the run and are not compared. A skipped fit compares only the lines it wrote.
baseline=${PROMOTE_NG_BASELINE:-$root/scripts/promote_ng_oracle.baseline}
if [ "$baseline" != none ]; then
  echo
  if [ ! -r "$baseline" ]; then
    echo "no baseline at $baseline; set PROMOTE_NG_BASELINE=none to only record" >&2
    exit 1
  fi
  slice=$(grep '^slice: ' "$out/digests.txt")
  if ! grep -qxF "$slice" "$baseline"; then
    echo "a different slice from $baseline's ($(grep '^slice: ' "$baseline" || echo none)): $slice" >&2
    exit 1
  fi
  missing=0
  for line in $(grep -v -e '^#' -e '^slice: ' "$out/digests.txt" | tr ' ' '|'); do
    if ! grep -v '^#' "$baseline" | tr ' ' '|' | grep -qxF "$line"; then
      echo "DIFFERS from $baseline: $(echo "$line" | tr '|' ' ')"
      missing=1
    fi
  done
  if [ "$missing" = 1 ]; then
    exit 1
  fi
  echo "every checksum matches $baseline"
fi
