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
           "host_load_start", "host_load_end", "other_cpu_cores", "source"]

# A time figure whose step ran while processes other than the suite used more than this many cores
# on average was not measured on a quiet machine, and is marked.
BUSY_OTHER_CPU = 1.0

# The read floor (lib/readfloor.py) set against the pvc step that reads the alignments:
# (floor stage, pvc caller, pvc stage).
READ_FLOOR_PAIRS = [
    ("read floor, one process per sample", "pvc (fitted)", "generate-psps"),
    ("read floor, every file at the run's threads", "pvc (defaults)", "call-from-alignments"),
]
TIME_METRICS = ("seconds", "_mb_largest_process", "_mb_all_processes")


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
    """The host's load, sampled every 5 s by run_suite.sh on the host: the one-minute load average,
    and the cores' worth of CPU used by processes other than the suite's virtual machine."""

    def __init__(self, path):
        self.times, self.loads, self.other = [], [], []
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) >= 2:
                        self.times.append(float(parts[0]))
                        self.loads.append(float(parts[1]))
                        self.other.append(float(parts[2]) if len(parts) > 2 else None)

    def _nearest(self, t):
        i = bisect.bisect_left(self.times, t)
        if i == len(self.times) or (i > 0 and t - self.times[i - 1] < self.times[i] - t):
            i -= 1
        return i

    def at(self, t):
        """The load average at the sample nearest `t`."""
        return round(self.loads[self._nearest(t)], 2) if self.times else ""

    def other_mean(self, start, end):
        """Mean cores used by other processes between `start` and `end`."""
        if not self.times:
            return ""
        lo, hi = bisect.bisect_left(self.times, start), bisect.bisect_right(self.times, end)
        window = [v for v in self.other[lo:hi] if v is not None]
        if not window:
            v = self.other[self._nearest(start)]
            window = [v] if v is not None else []
        return round(sum(window) / len(window), 2) if window else ""


def busy(row):
    """Was this time or memory figure measured while other processes used the fast cores?"""
    return row["other_cpu_cores"] not in ("", None) and float(row["other_cpu_cores"]) > BUSY_OTHER_CPU


def build_table(ctx):
    """Write `results.tsv` for the current pvc commit; return its rows."""
    load = HostLoad(os.path.join(ctx.results_root, "host_load.tsv"))
    rows = []
    sources = [(os.path.join(ctx.results_root, "competitors"), "-"), (ctx.pvc_root(), ctx.commit)]
    for directory, commit in sources:
        for record in _records(directory):
            measured = bool(record.get("measurement"))
            load_start = load.at(record["started"]) if measured else ""
            load_end = load.at(record["ended"]) if measured else ""
            other = load.other_mean(record["started"], record["ended"]) if measured else ""
            for r in record.get("rows", []):
                full = {c: "" for c in COLUMNS}
                full.update(r)
                full["pvc_commit"] = commit if full.get("caller", "").startswith("pvc") else "-"
                full["date"] = record["date"]
                full["machine"] = record["machine"]
                if r.get("metric", "").endswith(TIME_METRICS):
                    full["host_load_start"], full["host_load_end"] = load_start, load_end
                    full["other_cpu_cores"] = other
                full["source"] = record["name"]
                rows.append(full)
    rows.extend(_read_floor_ratios(rows))
    rows.sort(key=lambda r: (r["benchmark"], r["dataset"], _num(r["samples"]), _num(r["depth"]),
                             r["variant_class"], r["stage"], r["metric"], r["caller"]))
    os.makedirs(ctx.pvc_root(), exist_ok=True)
    path = os.path.join(ctx.pvc_root(), "results.tsv")
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _read_floor_ratios(rows):
    """pvc's wall time over the read floor's, for each configuration that has both."""
    walls = {}
    for r in rows:
        if r["metric"] == "wall_seconds":
            walls[(r["benchmark"], r["dataset"], str(r["samples"]), str(r["depth"]),
                   r["caller"], r["stage"])] = r
    out = []
    for floor_stage, pvc_caller, pvc_stage in READ_FLOOR_PAIRS:
        for key, floor in walls.items():
            if key[4] != "samtools (read floor)" or key[5] != floor_stage:
                continue
            pvc = walls.get(key[:4] + (pvc_caller, pvc_stage))
            if pvc is None or not float(floor["value"]):
                continue
            ratio = dict(pvc)
            ratio.update(stage=f"{pvc_stage} ÷ {floor_stage}", metric="wall_ratio",
                         value=round(float(pvc["value"]) / float(floor["value"]), 2),
                         host_load_start="", host_load_end="", other_cpu_cores="",
                         source=f"{pvc['source']} ÷ {floor['source']}")
            out.append(ratio)
    return out


def _read_floor_section(rows):
    """The summary's side-by-side table of the read floor and pvc's reading step."""
    walls = {(r["benchmark"], r["dataset"], str(r["samples"]), str(r["depth"]), r["caller"],
              r["stage"]): r for r in rows if r["metric"] == "wall_seconds"}
    lines = []
    for floor_stage, pvc_caller, pvc_stage in READ_FLOOR_PAIRS:
        for key, floor in sorted(walls.items(), key=lambda kv: (kv[0][0], kv[0][1], _num(kv[0][2]),
                                                                 _num(kv[0][3]))):
            if key[4] != "samtools (read floor)" or key[5] != floor_stage:
                continue
            pvc = walls.get(key[:4] + (pvc_caller, pvc_stage))
            if pvc is None:
                continue
            ratio = float(pvc["value"]) / float(floor["value"]) if float(floor["value"]) else 0
            quiet = "no ⚠" if busy(pvc) or busy(floor) else "yes"
            lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3]} | {pvc_caller} {pvc_stage} | "
                         f"{float(pvc['value']):,.1f} | {float(floor['value']):,.1f} | {ratio:.2f} | "
                         f"{quiet} |")
    if not lines:
        return []
    return ["## Read floor — pvc's time against the cost of only reading its input", "",
            "The read floor is `samtools view -c -M -L <regions>` over the same alignment files:",
            "every record decompressed and decoded, nothing written. No caller can read its input",
            "faster, so pvc's time over the floor is the share that is pvc's own work. pvc's",
            "`generate-psps` is set against the floor run the way it reads (one process per sample,",
            "as many at once as the run's cores, one decoding thread each for BAM and two for CRAM);",
            "`call-from-alignments` against the floor summed over the files, each read with the",
            "run's thread count. Every file was read once before each timed step, so every tool",
            "read from memory. A CRAM floor decodes exactly the fields pvc decodes",
            "(`required_fields=0x17ff`: all but the auxiliary tags other than the read group).",
            "Quiet: whether processes other than the suite used under one core on average during",
            "both runs.", "",
            "| benchmark | data | samples | depth | pvc step | pvc wall s | floor wall s | ratio | quiet |",
            "|---|---|---:|---|---|---:|---:|---:|---|"] + lines + [""]


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
                if busy(r):
                    value += " ⚠"
                cells.setdefault(k, {})[r["caller"]] = value
            for k in sorted(cells, key=lambda k: (_num(k[0]), k[1], k[2], k[3])):
                lines.append("| " + " | ".join(str(x) for x in k) + " | "
                             + " | ".join(cells[k].get(c, "") for c in callers) + " |")
            lines.append("")
    lines.extend(_read_floor_section(rows))
    if any(busy(r) for r in rows):
        lines.append(f"⚠ marks a time or memory figure measured while processes other than the suite "
                     f"used more than {BUSY_OTHER_CPU:g} core on average, so the machine was not quiet; "
                     "rerun that step with `--redo-step <source>` before quoting it.")
        lines.append("")
    path = os.path.join(ctx.pvc_root(), "summary.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path
