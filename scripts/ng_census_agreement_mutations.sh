#!/usr/bin/env bash
# Does the census-agreement test actually fail when the producer that reads a
# stored psp loses something?
#
# What it exercises is the test module `run::census_from_psp::the_two_producers_agree`:
# one sample's census built while its reads are walked, and built again from the
# psp that walk wrote, compared byte for byte. A test that compares two things
# built the same way passes whether or not either is right, so this writes
# deliberate defects into the psp-driven producer and reports which the
# comparison catches.
#
# Each defect is applied on its own and the file is restored afterwards,
# including on a failure or an interrupt. Expect all four to fail tests. The
# fourth used to be the exception and is not one any more: a census once held
# depth codes and allele counts and no per-read quality, so a change to a read's
# minted error could not move a census byte — but since the fit stage's Milestone
# C the census also carries each read group's minted read-error totals, and the
# fit needs them (ng_fit_stage_c_2026-09-05.md). Measured on 2026-09-10 at plan
# step E1: the fourth defect fails two of the module's three tests.
#
# **Two of three, and the third is not a survivor.** Only two of the module's
# tests compare the two censuses at all; the third, the_census_being_compared_
# holds_something, checks that the fixture's walk reached its kept loci, and a
# defect in one census's bytes is not what it is looking at.
#
# Run it from anywhere: the repository root is derived from this file.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO/src/run/census_from_psp.rs"
DEV="$REPO/scripts/dev.sh"
cp "$SRC" "$SRC.orig" || { echo "the file to mutate could not be backed up" >&2; exit 1; }
restore() { cp "$SRC.orig" "$SRC"; rm -f "$SRC.orig"; }
trap restore EXIT

# **Every step that writes the mutation is checked, and a failure stops the run.** Unchecked, a
# drifted anchor makes the python below raise, the tests then run against the *unmutated* file,
# and the harness prints `test result: ok` — reporting a defect as survived, which is the strongest
# wrong answer a mutation harness can give.
run_case() {
  local name="$1" find="$2" replace="$3"
  echo "### mutation: $name"
  cp "$SRC.orig" "$SRC" || { echo "  the file could not be put back before this case" >&2; exit 1; }
  if ! REPO="$REPO" FIND="$find" REPLACE="$replace" SRCF="$SRC" uv run --no-project python - <<'PY'
import os, pathlib
p = pathlib.Path(os.environ["SRCF"])
s = p.read_text()
find, replace = os.environ["FIND"], os.environ["REPLACE"]
assert find in s, "the mutation's anchor is not in the file"
p.write_text(s.replace(find, replace, 1))
PY
  then
    echo "  the mutation was not written — most likely its anchor has drifted out of $SRC" >&2
    exit 1
  fi
  "$DEV" cargo test --lib the_two_producers_agree 2>&1 | grep -E "^test result:|FAILED|panicked at" | head -5
}

run_case "the producer skips repeat-tract loci" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            if !matches!(record.kind, crate::locus_generation::LocusKind::Ssr(_)) {
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
                observation.read_group = crate::types::ReadGroupId(0);
            }
            writer.add_locus(&record);
        }'

# **Caught, and it was not always.** When this harness was written a census held a
# depth code per kept position per read group and the non-reference allele counts
# and no per-read quality, so nothing here could see a read's minted error move.
# The fit needed those totals, and the fit stage's Milestone C put them in the
# census — `CensusWriter::add_locus` folds every complete observation's `q_sum`
# into a per-read-group sum and `write_census` encodes it. So this defect now
# changes the bytes: measured 2026-09-10, both comparing tests fail, and on the
# fixture cohort of the day the two censuses first differed inside the read
# groups' minted totals.
run_case "one read's minted error arrives one step off" \
  'if let Some(record) = streamed.record.as_ref() {
            writer.add_locus(record);
        }' \
  'if let Some(record) = streamed.record.as_ref() {
            let mut record = record.clone();
            if let Some(observation) = record.observations.first_mut() {
                observation.q_sum = crate::types::SummedLogError::from_steps(
                    observation.q_sum.steps() - 1,
                );
            }
            writer.add_locus(&record);
        }'

echo "### restoring and re-running clean"
cp "$SRC.orig" "$SRC"
"$DEV" cargo test --lib the_two_producers_agree 2>&1 | grep -E "^test result:" | head -2
