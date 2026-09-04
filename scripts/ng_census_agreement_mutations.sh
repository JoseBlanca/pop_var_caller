#!/usr/bin/env bash
# Does the census-agreement test actually fail when the producer that reads a
# stored psp loses something?
#
# The test it exercises is `ng::run::census_from_psp::the_two_producers_agree`:
# one sample's census built while its reads are walked, and built again from the
# psp that walk wrote, compared byte for byte. A test that compares two things
# built the same way passes whether or not either is right, so this writes
# deliberate defects into the psp-driven producer and reports which the
# comparison catches.
#
# Each defect is applied on its own and the file is restored afterwards,
# including on a failure or an interrupt. Expect the first three to fail tests
# and the fourth to pass every one: a census records depth codes and allele
# counts and no per-read quality, so a change to a read's minted error cannot
# move a census byte. That is a fact about the format, not a hole in the test.
#
# Run it from anywhere: the repository root is derived from this file.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO/src/ng/run/census_from_psp.rs"
DEV="$REPO/scripts/dev.sh"
cp "$SRC" "$SRC.orig"
restore() { cp "$SRC.orig" "$SRC"; rm -f "$SRC.orig"; }
trap restore EXIT

run_case() {
  local name="$1" find="$2" replace="$3"
  cp "$SRC.orig" "$SRC"
  REPO="$REPO" FIND="$find" REPLACE="$replace" SRCF="$SRC" uv run --no-project python - <<'PY'
import os, pathlib
p = pathlib.Path(os.environ["SRCF"])
s = p.read_text()
find, replace = os.environ["FIND"], os.environ["REPLACE"]
assert find in s, "the mutation's anchor is not in the file"
p.write_text(s.replace(find, replace, 1))
PY
  echo "### mutation: $name"
  "$DEV" cargo test --lib the_two_producers_agree 2>&1 | grep -E "^test result:|FAILED|panicked at" | head -5
}

run_case "the producer skips repeat-tract loci" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            if !matches!(record.kind, crate::ng::locus_generation::LocusKind::Ssr(_)) {
                writer.add_locus(record);
            }
        }'

run_case "one read is lost at every locus" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            let mut record = record.clone();
            if !record.observations.is_empty() {
                record.observations.remove(0);
            }
            writer.add_locus(&record);
        }'

run_case "a read is credited to the wrong read group" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            let mut record = record.clone();
            for observation in &mut record.observations {
                observation.read_group = crate::ng::types::ReadGroupId(0);
            }
            writer.add_locus(&record);
        }'

# **Expected to be caught by nothing, and that is the finding.** A census holds a
# depth code per kept position per read group and the non-reference allele
# counts; no per-read quality reaches it. So the minted-error totals
# `RunParameters::assemble` wants cannot be read back out of a census file as the
# format stands.
run_case "one read's minted error arrives one step off" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            let mut record = record.clone();
            if let Some(observation) = record.observations.first_mut() {
                observation.q_sum = crate::ng::types::SummedLogError::from_steps(
                    observation.q_sum.steps() - 1,
                );
            }
            writer.add_locus(&record);
        }'

echo "### restoring and re-running clean"
cp "$SRC.orig" "$SRC"
"$DEV" cargo test --lib the_two_producers_agree 2>&1 | grep -E "^test result:" | head -2
