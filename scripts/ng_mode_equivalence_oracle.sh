#!/bin/sh
# **The mode-equivalence oracle on real reads** (doc/devel/ng/spec/run_streaming.md §12.3): one
# cohort of alignment files, called two ways, compared byte for byte.
#
#   call-from-alignments            — walk the CRAMs and call
#   generate-psps + call-from-psps  — walk the CRAMs into stored files, then call those
#
# The unit test `pop_var_caller_exp::mode_equivalence` makes the same comparison on a fixture of
# one contig and two samples. This is where the claim has weight: real reads, a real catalog, and
# a cohort big enough that anything the psp fails to carry has somewhere to show up.
#
# **The window coverage is compared too, and it is not in the VCF.** Each sample's read depth and
# GC around a locus travels beside the record rather than in it, so a run whose windows were all
# absent — or all one cover behind — writes exactly the file a correct run writes
# (doc/devel/ng/spec/window_coverage.md §10). Both routes are asked to write theirs down, and the
# two files are compared row for row — including a check that some row carries a number, because
# two runs that measured nothing agree perfectly.
#
# Usage:
#
#   scripts/ng_mode_equivalence_oracle.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <output-dir> <alignment.cram>...
#
# It expects a release build of `pop_var_caller_exp` in `target-container/release` or `target/
# release` and takes whichever is newer, because a machine with no container runtime builds to
# the second (see CLAUDE.md).
set -eu

if [ "$#" -lt 5 ]; then
  # Every comment line from the shebang to the first line that is not one. **Not a fixed range**:
  # the last line of the header has moved twice as this file grew, and both times the usage the
  # message exists to print was the part that fell off the end.
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

reference=$1; catalog=$2; regions=$3; out=$4
shift 4

root=$(cd "$(dirname "$0")/.." && pwd)
bin=""
for candidate in "$root/target-container/release/pop_var_caller_exp" \
                 "$root/target/release/pop_var_caller_exp"; do
  if [ -x "$candidate" ] && { [ -z "$bin" ] || [ "$candidate" -nt "$bin" ]; }; then
    bin=$candidate
  fi
done
if [ -z "$bin" ]; then
  echo "no release build of pop_var_caller_exp; build one first" >&2
  exit 1
fi

psps=$out/psps
mkdir -p "$psps"

alignments=""
for cram in "$@"; do
  alignments="$alignments --alignment $cram"
done

echo "=== direct mode ==="
# shellcheck disable=SC2086
NG_WINDOW_COVERAGE_FILE="$out/from_alignments.windows" "$bin" call-from-alignments \
  --reference "$reference" --catalog "$catalog" $alignments \
  --regions "$regions" \
  --defaults --threads 4 --output "$out/run.vcf" > "$out/alignments.log" 2>&1
mv "$out/run.vcf" "$out/from_alignments.vcf"
mv "$out/run.parameters.toml" "$out/from_alignments.parameters.toml"

echo "=== the walk, stored ==="
# shellcheck disable=SC2086
"$bin" generate-psps \
  --reference "$reference" --catalog "$catalog" $alignments \
  --regions "$regions" \
  --output-dir "$psps" --force > "$out/generate.log" 2>&1

echo "=== psp mode ==="
# **The psps are named in the order the alignment files' samples were**, read out of the VCF
# direct mode just wrote, so the two files put their sample columns in one order. Whether ANY
# order gives the same calls is §12.6's question and not this one's.
stored=""
for sample in $(grep '^#CHROM' "$out/from_alignments.vcf" | cut -f10-); do
  stored="$stored --psp $psps/$sample.psp"
done
# shellcheck disable=SC2086
NG_WINDOW_COVERAGE_FILE="$out/from_psps.windows" "$bin" call-from-psps \
  --reference "$reference" --catalog "$catalog" $stored \
  --defaults --threads 4 --output "$out/run.vcf" > "$out/psps.log" 2>&1
mv "$out/run.vcf" "$out/from_psps.vcf"
mv "$out/run.parameters.toml" "$out/from_psps.parameters.toml"

echo "=== the comparison ==="
# **`##commandline` is the one line that cannot match and must not**: it records which command
# was typed, and the two routes are two different commands. Everything else — every other header
# line and every record — is compared byte for byte. This is spec §12.1's timestamp exemption in
# another place, not a weakening of §12.3; the unit test, where both routes run inside one
# process and record the same line, compares the files whole with nothing filtered out.
for side in from_alignments from_psps; do
  grep -v '^##commandline=' "$out/$side.vcf" > "$out/$side.comparable"
done
if cmp -s "$out/from_alignments.comparable" "$out/from_psps.comparable"; then
  echo "IDENTICAL apart from ##commandline: psp mode's VCF is direct mode's"
  printf 'records: '; grep -vc '^#' "$out/from_alignments.vcf"
  shasum -a 256 "$out/from_alignments.comparable" "$out/from_psps.comparable"
else
  echo "DIFFERENT"
  diff "$out/from_alignments.comparable" "$out/from_psps.comparable" | head -40
  exit 1
fi

if cmp -s "$out/from_alignments.parameters.toml" "$out/from_psps.parameters.toml"; then
  echo "and the parameters file beside each is identical too"
else
  echo "the parameters files differ"
  diff "$out/from_alignments.parameters.toml" "$out/from_psps.parameters.toml" | head -20
  exit 1
fi

# **The window coverage.** Rows carry a sample index, a contig, a position and the two fields'
# bit patterns. Both files are sorted first: today the records leave the merge on the calling
# thread in genome order and the two files are already identical unsorted, but nothing in either
# mode promises that, and a sort costs one pass over a file the run has already written.
#
# **Two guards, and the second is the one that matters.** An empty file is a failure — it would
# mean neither route recorded anything. And a file whose every row is an *absent* pair is a
# failure too: two runs that measured nothing produce byte-identical files, so without this the
# comparison passes loudest exactly when the measurement is off. An absent pair is two `NaN`s,
# whose `f32` bit patterns all begin `21474836` (0x7FC00000 and up), so a measured row is one
# whose depth field is not one of those.
for side in from_alignments from_psps; do
  if [ ! -f "$out/$side.windows" ]; then
    echo "$side wrote no window file: the binary is older than NG_WINDOW_COVERAGE_FILE"
    exit 1
  fi
  sort "$out/$side.windows" > "$out/$side.windows.sorted"
done
rows=$(wc -l < "$out/from_alignments.windows.sorted" | tr -d ' ')
if [ "$rows" -eq 0 ]; then
  echo "no window was recorded, so nothing about the measurement was compared"
  exit 1
fi
measured=$(awk -F'\t' '$5 != "-" && $5 + 0 == $5 && $5 < 2143289344' \
  "$out/from_alignments.windows.sorted" | wc -l | tr -d ' ')
if [ "$measured" -eq 0 ]; then
  echo "every recorded window is absent, so the two modes agree about having measured nothing"
  exit 1
fi
if cmp -s "$out/from_alignments.windows.sorted" "$out/from_psps.windows.sorted"; then
  echo "and the two modes' window coverage is identical, bit for bit: \
$rows rows, $measured of them a measurement"
else
  echo "the two modes' window coverage differs"
  diff "$out/from_alignments.windows.sorted" "$out/from_psps.windows.sorted" | head -20
  exit 1
fi
