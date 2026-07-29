# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Merge batched `ng_ssr_cohort_stutter` dumps into one, renumbering read groups.

A cohort of this size cannot be dumped in one process — ~50 open CRAM readers exhaust the dev
container — so it is run in batches. That makes merging non-trivial: the `read_group` column is an
identifier minted per run, so batch 2's group 0 is a *different* group from batch 1's group 0.
Concatenating the files would silently fuse them.

The stable identity is `(file, rg_id)`: the SAM specification makes `@RG ID` unique within its file,
and the file path distinguishes the rest. This assigns a global number per such pair and rewrites
every row's `read_group` to match.

Usage:
    combine_stutter_batches.py OUT.tsv BATCH.tsv [BATCH.tsv ...]
"""
import sys
from pathlib import Path

if len(sys.argv) < 3:
    sys.exit(__doc__)

out_path = Path(sys.argv[1])
batches = [Path(p) for p in sys.argv[2:]]

global_ids: dict[tuple[str, str], int] = {}
rg_rows: list[list[str]] = []
rg_columns: str | None = None
data_header: str | None = None
total_rows = 0

with open(out_path, "w") as out:
    # First pass writes nothing: the merged `#rg` table has to precede the rows, so the batches are
    # read twice rather than buffering millions of rows in memory.
    for batch in batches:
        with open(batch) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if line.startswith("#rg\t"):
                    fields = line.rstrip("\n").split("\t")[1:]
                    key = (fields[-1], fields[1])  # (file, rg_id)
                    if key not in global_ids:
                        new_id = len(global_ids)
                        global_ids[key] = new_id
                        rg_rows.append([str(new_id)] + fields[1:])
                elif line.startswith("#rg_columns"):
                    rg_columns = line.rstrip("\n")

    for row in rg_rows:
        out.write("#rg\t" + "\t".join(row) + "\n")
    if rg_columns:
        out.write(rg_columns + "\n")

    for batch in batches:
        local_to_global: dict[str, str] = {}
        with open(batch) as fh:
            for line in fh:
                if line.startswith("#rg\t"):
                    fields = line.rstrip("\n").split("\t")[1:]
                    key = (fields[-1], fields[1])
                    local_to_global[fields[0]] = str(global_ids[key])
                    continue
                if line.startswith("#"):
                    continue
                if line.startswith("sample\t"):
                    if data_header is None:
                        data_header = line
                        out.write(line)
                    continue
                # Column 2 is `read_group`; empty for the synthetic no_border / capped tallies.
                parts = line.split("\t")
                if len(parts) > 1 and parts[1]:
                    parts[1] = local_to_global[parts[1]]
                    line = "\t".join(parts)
                out.write(line)
                total_rows += 1

print(f"{len(batches)} batches -> {out_path}")
print(f"{len(global_ids)} read groups, {total_rows:,} rows")
