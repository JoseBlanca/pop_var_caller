"""Turn the repeat-tract fitting bench's rows into one table per measure.

The bench (examples/ng_str_fit_bench.rs, run over a grid by scripts/research/str_fit_grid.sh) writes
one long row per method, class and number. This reads every rows file in a directory and prints, as
Markdown, one table per measure of the research plan (doc/devel/ng/research/str_fit_research_plan.md
§2):

1. Recovery: how far each estimate is from the truth, over the seeds — bias and root mean square
   error, relative to the truth for the slippage level and the concentration. One line per cell,
   method and number, summarised over the classes (the median class); the per-class figures go to
   recovery_by_class.txt.
2. Error honesty: of the estimates that came with a standard error, the share within one and within
   two errors of the truth (about 68 and 95 in 100 if the errors are honest).
3. Held-out prediction: the log-likelihood a read and a locus of the held-out loci, as the method
   minus the truth's own on the same loci (zero or below; the closer to zero the better), with the
   Monte Carlo error.
4. Cost: seconds to fit every locus, seconds for the held-out fits, peak memory.

Usage:
    uv run --no-project --with numpy python scripts/research/str_fit_summarise.py DIR
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

CELL = ("period", "samples", "depth", "truth_variant")
SLIPPAGE = ("level", "shorter_share", "fall_off")
RELATIVE = {"level", "concentration"}
HEADLINE = ("level", "shorter_share", "fall_off", "between_locus_share", "share_0", "share_-1",
            "share_1", "share_-9", "share_9")


def number(text):
    return None if text in ("NA", "", None) else float(text)


def read_rows(directory):
    rows = []
    for path in sorted(Path(directory).glob("*.tsv")):
        if path.name.endswith(".truth.tsv") or path.name.startswith("recovery_by_class"):
            continue
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or "parameter" not in reader.fieldnames:
                continue
            for row in reader:
                rows.append(row)
    return rows


def cell_key(row):
    return tuple(row[column] for column in CELL)


def cell_label(key):
    period, samples, depth, truth = key
    return f"period {period}, {samples} samples, {float(depth):g} reads, {truth}"


def recovery(rows, by_class_path):
    # (cell, method, parameter, repeats) -> list of (truth, estimate, error)
    draws = defaultdict(list)
    for row in rows:
        parameter = row["parameter"]
        if not (parameter in SLIPPAGE or parameter in ("concentration", "between_locus_share")
                or parameter.startswith("share_")):
            continue
        truth = number(row["truth"])
        estimate = number(row["estimate"])
        draws[(cell_key(row), row["method"], parameter, row["repeats"])].append(
            (truth, estimate, number(row["standard_error"])))
    per_class = {}
    with open(by_class_path, "w") as out:
        out.write("cell\tmethod\tparameter\trepeats\tdraws\tanswered\ttruth\tbias\trmse\t"
                  "relative_bias\trelative_rmse\n")
        for (cell, method, parameter, repeats), values in sorted(draws.items()):
            answered = [(t, e) for t, e, _ in values if e is not None and t is not None]
            if answered:
                truth = np.array([t for t, _ in answered])
                estimate = np.array([e for _, e in answered])
                difference = estimate - truth
                bias = float(difference.mean())
                rmse = float(math.sqrt((difference ** 2).mean()))
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative = difference / truth
                relative_bias = float(relative.mean()) if np.all(truth > 0) else float("nan")
                relative_rmse = (float(math.sqrt((relative ** 2).mean()))
                                 if np.all(truth > 0) else float("nan"))
                truth_value = float(truth.mean())
            else:
                bias = rmse = relative_bias = relative_rmse = truth_value = float("nan")
            per_class[(cell, method, parameter, repeats)] = (
                len(values), len(answered), bias, rmse, relative_bias, relative_rmse)
            out.write(f"{cell_label(cell)}\t{method}\t{parameter}\t{repeats}\t{len(values)}\t"
                      f"{len(answered)}\t{truth_value:.6g}\t{bias:.6g}\t{rmse:.6g}\t"
                      f"{relative_bias:.6g}\t{relative_rmse:.6g}\n")

    print("## Recovery\n")
    print("Median over the classes of each class's bias and root mean square error over the seeds; "
          "relative to the truth for the level and the concentration. `answered` counts classes "
          "with an estimate in every seed, of all classes.\n")
    print("| cell | method | number | answered | median bias | median RMSE |")
    print("|---|---|---|---|---|---|")
    grouped = defaultdict(list)
    for (cell, method, parameter, repeats), value in per_class.items():
        grouped[(cell, method, parameter)].append(value)
    for (cell, method, parameter), values in sorted(grouped.items()):
        if parameter not in HEADLINE and parameter != "concentration":
            continue
        relative = parameter in RELATIVE
        biases = [v[4] if relative else v[2] for v in values if v[1] > 0]
        rmses = [v[5] if relative else v[3] for v in values if v[1] > 0]
        answered = sum(1 for v in values if v[1] == v[0] and v[0] > 0)
        suffix = " (relative)" if relative else ""
        if biases:
            print(f"| {cell_label(cell)} | {method} | {parameter}{suffix} | {answered} of "
                  f"{len(values)} | {np.nanmedian(biases):+.4f} | {np.nanmedian(rmses):.4f} |")
        else:
            print(f"| {cell_label(cell)} | {method} | {parameter}{suffix} | 0 of {len(values)} "
                  f"| NA | NA |")
    print()


def honesty(rows):
    counts = defaultdict(lambda: [0, 0, 0])
    for row in rows:
        error = number(row["standard_error"])
        truth = number(row["truth"])
        estimate = number(row["estimate"])
        parameter = row["parameter"]
        if parameter.startswith("heldout") or error is None or truth is None or estimate is None:
            continue
        group = parameter if not parameter.startswith("share_") else "shares"
        key = (cell_key(row), row["method"], group)
        counts[key][0] += 1
        distance = abs(estimate - truth)
        counts[key][1] += distance <= error
        counts[key][2] += distance <= 2 * error
    print("## Error honesty\n")
    print("Of the estimates that came with a standard error, the share whose truth lies within one "
          "and within two errors; honest errors give about 0.68 and 0.95.\n")
    if not counts:
        print("No method gave standard errors.\n")
        return
    print("| cell | method | number | estimates with an error | within 1 | within 2 |")
    print("|---|---|---|---|---|---|")
    for (cell, method, group), (total, one, two) in sorted(counts.items()):
        print(f"| {cell_label(cell)} | {method} | {group} | {total} | {one / total:.2f} | "
              f"{two / total:.2f} |")
    print()


def held_out(rows):
    by = defaultdict(dict)
    for row in rows:
        if row["repeats"] != "all" or not row["parameter"].startswith("heldout_"):
            continue
        key = (cell_key(row), row["method"], row["seed"])
        by[key][row["parameter"]] = (number(row["truth"]), number(row["estimate"]),
                                     number(row["standard_error"]))
    summary = defaultdict(list)
    for (cell, method, seed), values in by.items():
        if "heldout_ln_per_read" not in values:
            continue
        per_read = values["heldout_ln_per_read"]
        per_locus = values["heldout_ln_per_locus"]
        loci = values.get("heldout_loci", (None, float("nan"), None))[1]
        summary[(cell, method)].append((per_read[1] - per_read[0], per_read[2],
                                        per_locus[1] - per_locus[0], per_locus[2], loci))
    print("## Held-out prediction\n")
    print("The log-likelihood of the held-out loci's reads under the method's fit, minus under the "
          "truth, a read and a locus; mean over the seeds, with the Monte Carlo error of the "
          "method's score. Zero is as good as the truth.\n")
    if not summary:
        print("No held-out scores.\n")
        return
    print("| cell | method | seeds | loci a seed | per read | per locus | Monte Carlo error a locus |")
    print("|---|---|---|---|---|---|---|")
    for (cell, method), values in sorted(summary.items()):
        values = np.array(values, dtype=float)
        print(f"| {cell_label(cell)} | {method} | {len(values)} | {np.nanmean(values[:, 4]):.0f} | "
              f"{np.nanmean(values[:, 0]):+.5f} | {np.nanmean(values[:, 2]):+.4f} | "
              f"{np.nanmean(values[:, 3]):.4f} |")
    print()


def cost(rows):
    by = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row["repeats"] != "all":
            continue
        parameter = row["parameter"]
        if parameter in ("fit_seconds", "heldout_fit_seconds", "heldout_score_seconds",
                         "peak_rss_mib", "draw_seconds", "rows_with_reads"):
            value = number(row["estimate"])
            if value is not None:
                method = row["method"].removesuffix("_raw")
                by[(cell_key(row), method)][parameter].append(value)
    print("## Cost\n")
    print("Mean over the seeds. Peak memory is the process's, the drawn cohort included "
          "(`draw` rows: the draw alone, where a run without methods measured it).\n")
    print("| cell | method | fit, s | held-out fits, s | held-out scoring, s | peak memory, MiB |")
    print("|---|---|---|---|---|---|")
    for (cell, method), values in sorted(by.items()):
        if method in ("draw", "process", "truth"):
            continue

        def mean(name):
            return f"{np.mean(values[name]):.1f}" if values.get(name) else "NA"

        print(f"| {cell_label(cell)} | {method} | {mean('fit_seconds')} | "
              f"{mean('heldout_fit_seconds')} | {mean('heldout_score_seconds')} | "
              f"{mean('peak_rss_mib')} |")
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("directory")
    arguments = parser.parse_args()
    rows = read_rows(arguments.directory)
    if not rows:
        sys.exit(f"no rows files in {arguments.directory}")
    recovery(rows, Path(arguments.directory) / "recovery_by_class.txt")
    honesty(rows)
    held_out(rows)
    cost(rows)


if __name__ == "__main__":
    main()
