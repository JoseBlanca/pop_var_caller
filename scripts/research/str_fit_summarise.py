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

With `--compare B_raw,B_curves,D`, it prints instead the methods side by side: the held-out score
with its standard error over seeds and the paired difference D − B, the two held-out fallbacks at one
sample (the bench's `--both-fallbacks`), recovery and error honesty with shares grouped by size.

Usage:
    uv run --no-project --with numpy python scripts/research/str_fit_summarise.py DIR [DIR ...] [--compare ...]
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


def read_rows(directories):
    rows = []
    paths = sorted(path for directory in directories for path in Path(directory).glob("*.tsv"))
    for path in paths:
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
        if None in (per_read[0], per_read[1], per_locus[0], per_locus[1]):
            continue
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
                         "peak_rss_mib", "draw_seconds", "rows_with_reads",
                         "process_cpu_seconds"):
            value = number(row["estimate"])
            if value is not None:
                method = row["method"].removesuffix("_raw")
                by[(cell_key(row), method)][parameter].append(value)
    print("## Cost\n")
    print("Mean over the seeds. Peak memory is the process's, the drawn cohort included "
          "(`draw` rows: the draw alone, where a run without methods measured it).\n")
    print("| cell | method | fit, s | held-out fits, s | held-out scoring, s | process CPU, s | "
          "peak memory, MiB |")
    print("|---|---|---|---|---|---|---|")
    for (cell, method), values in sorted(by.items()):
        if method in ("draw", "process", "truth"):
            continue

        def mean(name):
            return f"{np.mean(values[name]):.1f}" if values.get(name) else "NA"

        print(f"| {cell_label(cell)} | {method} | {mean('fit_seconds')} | "
              f"{mean('heldout_fit_seconds')} | {mean('heldout_score_seconds')} | "
              f"{mean('process_cpu_seconds')} | {mean('peak_rss_mib')} |")
    print()


def mean_and_error(values):
    """The mean over seeds and its standard error (the seeds' spread over the root of their count)."""
    values = np.array([v for v in values if v is not None and not math.isnan(v)], dtype=float)
    if len(values) == 0:
        return float("nan"), float("nan"), 0
    error = float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else float("nan")
    return float(values.mean()), error, len(values)


def show(mean, error, digits=3):
    if math.isnan(mean):
        return "—"
    if math.isnan(error):
        return f"{mean:+.{digits}f}"
    return f"{mean:+.{digits}f} ± {error:.{digits}f}"


def held_out_by_seed(rows):
    """(cell, variant) -> {seed: (per-locus score minus the truth's, loci, loci on the fallback)}."""
    by = defaultdict(dict)
    for row in rows:
        if row["repeats"] != "all" or not row["parameter"].startswith("heldout_"):
            continue
        by[(cell_key(row), row["method"], row["seed"])][row["parameter"]] = (
            number(row["truth"]), number(row["estimate"]))
    out = defaultdict(dict)
    for (cell, method, seed), values in by.items():
        per_locus = values.get("heldout_ln_per_locus")
        if per_locus is None or None in per_locus:
            continue
        loci = values.get("heldout_loci", (None, float("nan")))[1]
        fallback = values.get("heldout_loci_on_fallback", (None, float("nan")))[1]
        out[(cell, method)][seed] = (per_locus[1] - per_locus[0], loci, fallback)
    return out


def comparison(rows, methods):
    """The head-to-head tables: held-out score, recovery and error honesty side by side."""
    scores = held_out_by_seed(rows)
    cells = sorted({cell for cell, _ in scores})
    print("## Held-out prediction, side by side\n")
    print("Log-likelihood a held-out locus under each fit minus under the truth (zero is as good as "
          "the truth; below zero is worse), at calling's fallback for a class without its own "
          "numbers. Mean over the seeds ± its standard error. The last columns are the paired "
          "difference, seed by seed, of D against B: above zero means D predicts better.\n")
    head = " | ".join(methods)
    pairs = [("D", "B_raw")] + [(m, "B_curves") for m in methods if not m.startswith("B")]
    pair_head = " | ".join(f"{m} − {b}" for m, b in pairs)
    print(f"| cell | seeds | loci a seed | {head} | {pair_head} |")
    print("|---|---|---|" + "---|" * (len(methods) + len(pairs)))
    for cell in cells:
        columns = []
        seeds_seen = set()
        loci = []
        for method in methods:
            values = scores.get((cell, method), {})
            seeds_seen |= set(values)
            loci += [v[1] for v in values.values()]
            columns.append(show(*mean_and_error([v[0] for v in values.values()])[:2]))
        differences = []
        for method, base in pairs:
            d = scores.get((cell, method), {})
            b = scores.get((cell, base), {})
            paired = [d[seed][0] - b[seed][0] for seed in d if seed in b]
            differences.append(show(*mean_and_error(paired)[:2]))
        print(f"| {cell_label(cell)} | {len(seeds_seen)} | {np.nanmean(loci):.0f} | "
              f"{' | '.join(columns)} | {' | '.join(differences)} |")
    print()

    one_sample = [cell for cell in cells if cell[1] == "1"]
    if one_sample:
        print("## Held-out prediction at one sample: the two fallbacks\n")
        print("Calling's fallback scores a class with no concentration of its own at a flat shape and "
              "the stated concentration; the other keeps the class's own shares at the stated "
              "concentration. D at one sample has no concentration in any class, so every class falls "
              "back; `D_working_F` is D scored at the concentration it worked at. The last column is "
              "the share of held-out loci whose class B left without a concentration. Mean over the "
              "seeds ± its standard error, a locus, against the truth; the paired column is seed by "
              "seed, above zero where D predicts better.\n")
        print("| cell | B_raw: calling / own | B_curves: calling / own | D: calling / own | "
              "D_working_F | D − B_curves, own shares | B classes on fallback |")
        print("|---|---|---|---|---|---|---|")
        for cell in one_sample:
            def pair(method):
                calling = scores.get((cell, method), {})
                own = scores.get((cell, f"{method}_own_shares"), {})
                return (show(*mean_and_error([v[0] for v in calling.values()])[:2]) + " / "
                        + show(*mean_and_error([v[0] for v in own.values()])[:2]))
            d = scores.get((cell, "D_own_shares"), {})
            b = scores.get((cell, "B_curves_own_shares"), {})
            paired = [d[seed][0] - b[seed][0] for seed in d if seed in b]
            calling = scores.get((cell, "B_curves"), {})
            on = [v[2] / v[1] for v in calling.values() if v[1]]
            working = scores.get((cell, "D_working_F"), {})
            print(f"| {cell_label(cell)} | {pair('B_raw')} | {pair('B_curves')} | {pair('D')} | "
                  f"{show(*mean_and_error([v[0] for v in working.values()])[:2])} | "
                  f"{show(*mean_and_error(paired)[:2])} | "
                  f"{np.mean(on) if on else float('nan'):.2f} |")
        print()

    # Recovery side by side: the median over classes of each class's RMSE over seeds.
    draws = defaultdict(list)
    for row in rows:
        parameter = row["parameter"]
        if parameter not in ("level", "shorter_share", "fall_off", "between_locus_share",
                             "share_0"):
            continue
        truth, estimate = number(row["truth"]), number(row["estimate"])
        if truth is None:
            continue
        draws[(cell_key(row), row["method"], parameter, row["repeats"])].append(
            None if estimate is None else
            ((estimate - truth) / truth if parameter == "level" else estimate - truth))
    numbers = ("level", "shorter_share", "fall_off", "between_locus_share", "share_0")
    print("## Recovery, side by side\n")
    print("Median over the classes of each class's root mean square error over the seeds: relative "
          "for the slip rate (`level`), absolute otherwise. A dash: the method gave no estimate "
          "(F at one sample is absent by design for B and D).\n")
    print("| cell | method | " + " | ".join(numbers) + " |")
    print("|---|---|" + "---|" * len(numbers))
    recovery_cells = sorted({key[0] for key in draws})
    for cell in recovery_cells:
        for method in methods:
            entries = []
            for parameter in numbers:
                per_class = []
                for (c, m, p, _), values in draws.items():
                    if c != cell or m != method or p != parameter:
                        continue
                    answered = [v for v in values if v is not None]
                    if answered:
                        per_class.append(math.sqrt(np.mean(np.square(answered))))
                entries.append(f"{np.median(per_class):.4f}" if per_class else "—")
            if all(entry == "—" for entry in entries):
                continue
            print(f"| {cell_label(cell)} | {method} | " + " | ".join(entries) + " |")
    print()

    # Honesty side by side, shares grouped by the size of their truth.
    counts = defaultdict(lambda: [0, 0, 0])
    for row in rows:
        error, truth, estimate = (number(row["standard_error"]), number(row["truth"]),
                                  number(row["estimate"]))
        parameter = row["parameter"]
        if parameter.startswith("heldout") or None in (error, truth, estimate):
            continue
        if parameter.startswith("share_"):
            group = ("shares ≥ 0.05" if truth >= 0.05 else
                     "shares 0.01–0.05" if truth >= 0.01 else "shares < 0.01")
        elif parameter in ("level", "shorter_share", "fall_off", "between_locus_share"):
            group = parameter
        else:
            continue
        key = (cell_key(row), row["method"], group)
        counts[key][0] += 1
        counts[key][1] += abs(estimate - truth) <= error
        counts[key][2] += abs(estimate - truth) <= 2 * error
    groups = ("level", "shorter_share", "fall_off", "between_locus_share", "shares ≥ 0.05",
              "shares 0.01–0.05", "shares < 0.01")
    print("## Error honesty, side by side\n")
    print("Share of estimates within one / two of their reported errors of the truth, pooled over "
          "classes and seeds (honest: about 0.68 / 0.95); in brackets, how many estimates had an "
          "error. Only the fits as they come carry errors (the curves replace them).\n")
    print("| cell | method | " + " | ".join(groups) + " |")
    print("|---|---|" + "---|" * len(groups))
    for cell in sorted({key[0] for key in counts}):
        for method in methods:
            entries = []
            for group in groups:
                total, one, two = counts.get((cell, method, group), (0, 0, 0))
                entries.append(f"{one / total:.2f} / {two / total:.2f} ({total})" if total else "—")
            if all(entry == "—" for entry in entries):
                continue
            print(f"| {cell_label(cell)} | {method} | " + " | ".join(entries) + " |")
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("directory", nargs="+", help="one or more directories of rows files")
    parser.add_argument("--compare", default=None,
                        help="comma-separated method variants to set side by side (e.g. "
                             "B_raw,B_curves,D,D_working_F); prints only the side-by-side tables")
    arguments = parser.parse_args()
    rows = read_rows(arguments.directory)
    if not rows:
        sys.exit(f"no rows files in {' '.join(arguments.directory)}")
    if arguments.compare:
        comparison(rows, arguments.compare.split(","))
        return
    recovery(rows, Path(arguments.directory[0]) / "recovery_by_class.txt")
    honesty(rows)
    held_out(rows)
    cost(rows)


if __name__ == "__main__":
    main()
