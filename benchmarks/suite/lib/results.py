"""The results table and the summary generated from it.

Every finished step's `.done` file carries its rows. The table for one pvc commit is the union of
that commit's pvc rows and the competitors' stored rows, so a pvc refresh never needs a competitor
rerun to produce a complete table. Each row keeps the date its step ran, which is how the table shows
that a competitor's figures are older than pvc's.
"""

import bisect
import csv
import json
import os
import time

COLUMNS = ["benchmark", "dataset", "samples", "depth", "region_bp", "variant_class", "caller",
           "caller_version", "pvc_commit", "stage", "metric", "value", "date", "machine",
           "host_load_mean", "source"]

# A step whose mean host load is above this ran while something else used the fast cores.
BUSY_LOAD = 7.0


def row(base, **fields):
    """A result row: the configuration in `base`, overridden or extended by `fields`."""
    out = dict(base)
    out.update(fields)
    return out


def time_rows(base, measurement, stage):
    """The time and memory rows of one measured stage."""
    if not measurement:
        return []
    return [
        row(base, stage=stage, metric="wall_seconds", value=measurement["wall_s"]),
        row(base, stage=stage, metric="cpu_seconds",
            value=round(measurement["user_s"] + measurement["sys_s"], 3)),
        row(base, stage=stage, metric="peak_rss_mb_largest_process",
            value=round(measurement["max_single_rss_bytes"] / 2**20, 1)),
        row(base, stage=stage, metric="peak_rss_mb_all_processes",
            value=round(measurement["peak_summed_rss_bytes"] / 2**20, 1)),
    ]


def _records(directory):
    for dirpath, dirnames, filenames in os.walk(directory):
        if ".done" in filenames:
            with open(os.path.join(dirpath, ".done")) as fh:
                yield json.load(fh)


class HostLoad:
    """The host's load average, sampled every 30 s by run_suite.sh on the host."""

    def __init__(self, path):
        self.times, self.loads = [], []
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) == 2:
                        self.times.append(float(parts[0]))
                        self.loads.append(float(parts[1]))

    def mean(self, start, end):
        if not self.times:
            return ""
        lo = bisect.bisect_left(self.times, start)
        hi = bisect.bisect_right(self.times, end)
        window = self.loads[lo:hi]
        if not window:
            # A step shorter than the sampling interval: take the nearest sample.
            i = min(max(lo, 0), len(self.loads) - 1)
            window = [self.loads[i]]
        return round(sum(window) / len(window), 2)


def build_table(ctx):
    """Write `results.tsv` for the current pvc commit; return its rows."""
    load = HostLoad(os.path.join(ctx.results_root, "host_load.tsv"))
    rows = []
    sources = [(os.path.join(ctx.results_root, "competitors"), "-"), (ctx.pvc_root(), ctx.commit)]
    for directory, commit in sources:
        for record in _records(directory):
            host_load = load.mean(record["started"], record["ended"]) if record.get("measurement") else ""
            for r in record.get("rows", []):
                full = {c: "" for c in COLUMNS}
                full.update(r)
                full["pvc_commit"] = commit if full.get("caller", "").startswith("pvc") else "-"
                full["date"] = record["date"]
                full["machine"] = record["machine"]
                full["host_load_mean"] = host_load if r.get("metric", "").endswith(
                    ("seconds", "_mb_largest_process", "_mb_all_processes")) else ""
                full["source"] = record["name"]
                rows.append(full)
    rows.sort(key=lambda r: (r["benchmark"], r["dataset"], _num(r["samples"]), _num(r["depth"]),
                             r["variant_class"], r["stage"], r["metric"], r["caller"]))
    os.makedirs(ctx.pvc_root(), exist_ok=True)
    path = os.path.join(ctx.pvc_root(), "results.tsv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _num(value):
    try:
        return float(str(value).rstrip("x"))
    except ValueError:
        return float("inf")


def _fmt(value):
    if isinstance(value, float) or (isinstance(value, str) and value.replace(".", "", 1).isdigit()
                                    and "." in value):
        v = float(value)
        if abs(v) < 1:
            return f"{v:.4f}"
        if abs(v) < 100:
            return f"{v:.1f}"
        return f"{v:,.0f}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def write_summary(ctx, rows, benchmarks):
    """Write `summary.md`: per benchmark, what it measured on, and its figures by caller."""
    lines = [f"# pvc benchmark figures — commit {ctx.commit}", "",
             f"Generated {time.strftime('%Y-%m-%d %H:%M')} from `results.tsv` beside this file.",
             "Every figure names its benchmark, its data, its depth and the caller version that",
             "produced it. Competitor figures carry the date they were measured; they are rerun",
             "only when their version or their input changes.", ""]
    machines = sorted({r["machine"] for r in rows if r["machine"]})
    for m in machines:
        lines.append(f"Machine: {m}")
    lines.append("")
    by_bench = {}
    for r in rows:
        by_bench.setdefault(r["benchmark"], []).append(r)
    order = [b for b in benchmarks if b in by_bench] + sorted(set(by_bench) - set(benchmarks))
    for bench in order:
        module = benchmarks.get(bench)
        lines.append(f"## {getattr(module, 'TITLE', bench)}")
        lines.append("")
        if module is not None and getattr(module, "DESCRIPTION", ""):
            lines.append(module.DESCRIPTION.strip())
            lines.append("")
        groups = {}
        for r in by_bench[bench]:
            key = (r["dataset"], r["samples"], r["region_bp"])
            groups.setdefault(key, []).append(r)
        for (dataset, samples, region_bp), group in sorted(groups.items(), key=lambda kv: (kv[0][0], _num(kv[0][1]))):
            region = f"{int(region_bp):,} bp" if str(region_bp).isdigit() else region_bp
            lines.append(f"### {dataset} — {samples} sample(s), {region}")
            lines.append("")
            callers = sorted({r["caller"] for r in group})
            versions = {}
            for r in group:
                versions.setdefault(r["caller"], set()).add(f"{r['caller_version']} ({r['date']})")
            for c in callers:
                lines.append(f"- {c}: {', '.join(sorted(versions[c]))}")
            lines.append("")
            lines.append("| depth | class | stage | metric | " + " | ".join(callers) + " |")
            lines.append("|---|---|---|---|" + "---:|" * len(callers))
            cells = {}
            for r in group:
                k = (r["depth"], r["variant_class"], r["stage"], r["metric"])
                value = _fmt(r["value"])
                if r["host_load_mean"] != "" and float(r["host_load_mean"]) > BUSY_LOAD:
                    value += " ⚠"
                cells.setdefault(k, {})[r["caller"]] = value
            for k in sorted(cells, key=lambda k: (_num(k[0]), k[1], k[2], k[3])):
                lines.append("| " + " | ".join(str(x) for x in k) + " | "
                             + " | ".join(cells[k].get(c, "") for c in callers) + " |")
            lines.append("")
    if any(r["host_load_mean"] != "" and float(r["host_load_mean"]) > BUSY_LOAD for r in rows):
        lines.append(f"⚠ marks a time or memory figure measured while the host's mean load was above "
                     f"{BUSY_LOAD:g}, so something else was using the fast cores; rerun that step "
                     "with `--redo-step <source>` before quoting it.")
        lines.append("")
    path = os.path.join(ctx.pvc_root(), "summary.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path
