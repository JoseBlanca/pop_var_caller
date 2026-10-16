"""Write the stuttering benchmark's parameters file with the **true** stutter in it.

    python3 scripts/ng_str_stutter_true_parameters.py <run-dir>

Reads `<run-dir>/ng_defaults.parameters.toml` (the `--defaults` run's parameters, which
ng_str_stutter_benchmark.sh copies out of that run's log) and replaces its two stutter tables
with the stutter the simulator drew the reads under, from `stutter_rates.tsv` and
`libraries.tsv`; everything else stays at the defaults. Writes
`<run-dir>/true_stutter.parameters.toml`, for `call-from-psps --parameters`.

Why (owner, 2026-10-09): the benchmark tests the STR calling, not the fit, and the fit is slow.
Called with the true stutter, ng shows its STR calls with a stutter model that is right; if that
agrees with the fitted run, later runs can skip the fit.

**How the rows are laid out, and why.**
- **One slippage group per read group**, because every library has its own multiplier on its
  level's rate (spec §3.3). The read groups' samples come from the file's own read-group table.
- **A row for every repeat count from 1 to 60 at every motif length.** The caller looks a row up
  by the *candidate allele's* repeat count, not the reference's
  (`src/calling/inference/repeat_tract_parameters.rs`), and a missing row silently falls back to
  the shipped stutter model; 60 is past any allele the simulator draws. The rate at a repeat
  count is the share of molecules that stutter at an allele of that length — the simulator's
  rate, times the library's multiplier, capped at 0.95 as the simulator caps it.
- `shorter_share` is the library level's; `fall_off`, the share of stutter that moves more than
  one unit, is 1 − 0.85 = 0.15 at every row. The simulator caps a move at four units, which the
  caller's geometric step does not; at 0.15 the mass beyond four units is 5 in 10,000.

The origins record nothing fitted (`expected_slipped_reads = 0.0`), so the run's report does not
present these numbers as measured.
"""

import os
import re
import sys

LONGEST_REPEATS = 60
FALL_OFF = 0.15
HIGHEST_STUTTER_SHARE = 0.95


def main():
    run = sys.argv[1]
    template = open(os.path.join(run, "ng_defaults.parameters.toml")).read()

    rates = {}
    with open(os.path.join(run, "stutter_rates.tsv")) as handle:
        handle.readline()
        for line in handle:
            level, period, repeats, rate = line.rstrip("\n").split("\t")
            rates[(level, int(period), int(repeats))] = float(rate)
    longest_tabled = max(repeats for _, _, repeats in rates)

    libraries = {}
    with open(os.path.join(run, "libraries.tsv")) as handle:
        handle.readline()
        for line in handle:
            sample, level, multiplier, shorter = line.rstrip("\n").split("\t")
            libraries[sample] = (level, float(multiplier), float(shorter))

    # The file's read groups, in its own order: `{ read_group = 0, declared_id = …, sample = … }`.
    read_groups = [(int(group), sample) for group, sample in re.findall(
        r'\{ read_group = (\d+), declared_id = "[^"]*", library = "[^"]*", sample = "([^"]*)" \}',
        template)]
    if sorted(sample for _, sample in read_groups) != sorted(libraries):
        sys.exit(f"read groups {read_groups} do not match the libraries {sorted(libraries)}")

    groups = "slippage_group_by_read_group = [\n" + "".join(
        f"    {{ read_group = {group}, slippage_group = {group} }},\n" for group, _ in read_groups) + "]"
    rows = []
    for group, sample in read_groups:
        level, multiplier, shorter = libraries[sample]
        for period in range(1, 7):
            for repeats in range(1, LONGEST_REPEATS + 1):
                if level == "none":
                    rate = 0.0
                else:
                    rate = rates[(level, period, min(repeats, longest_tabled))]
                share = min(rate * multiplier, HIGHEST_STUTTER_SHARE)
                rows.append(
                    f"    {{ period = {period}, reference_repeats = {repeats}, slippage_group = {group}, "
                    f"share_of_reads_that_slip = {share!r}, shorter_share = {shorter!r}, "
                    f"fall_off = {FALL_OFF!r}, share_of_reads_that_slip_origin = "
                    f'{{ smoothing = "this_stratum", expected_slipped_reads = 0.0 }}, '
                    f"shorter_share_and_fall_off_origin = {{ expected_slipped_reads = 0.0, "
                    f'shorter_share_smoothing = "this_stratum", fall_off_smoothing = "this_stratum" }} }},\n')
    table = "slippage_by_stratum_and_group = [\n" + "".join(rows) + "]"

    out, replaced = re.subn(r"slippage_group_by_read_group = \[\n(?:    .*\n)*\]", groups, template)
    if replaced != 1:
        sys.exit("no slippage_group_by_read_group table in the defaults file")
    out, replaced = re.subn(r"slippage_by_stratum_and_group = \[\]", table, out)
    if replaced != 1:
        sys.exit("the defaults file's slippage_by_stratum_and_group is not the empty table")
    with open(os.path.join(run, "true_stutter.parameters.toml"), "w") as handle:
        handle.write(out)
    print(f"{len(rows)} stutter rows over {len(read_groups)} read groups")


if __name__ == "__main__":
    main()
