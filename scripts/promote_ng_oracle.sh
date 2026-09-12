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
# **PROMOTE_NG_REGIONS caps how much genome is walked**, taking that many leading lines of
# the regions file (default 20). Identity is the question, not coverage: a byte that moves
# moves in the first region as readily as the hundredth, and a cheap oracle is one that gets
# run at every checkpoint rather than once. The count is written into the digest file, so two
# runs over different slices cannot be compared by accident.
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
  # `|| true`: `grep -c` exits 1 on a count of zero, which `set -e` would turn into an abort
  # here rather than into the empty-output failure the wrapped oracle already reports.
  echo "# records: $(grep -vc '^#' "$out/from_alignments.vcf" || true)"
  for artefact in from_alignments.comparable \
                  from_psps.comparable \
                  from_alignments.parameters.toml \
                  from_alignments.windows.sorted \
                  from_alignments.windows.histograms; do
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
