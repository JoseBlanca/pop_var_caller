"""The penalised fit's reference fixture: record it, re-run its cells, compare.

The research bench (`examples/ng_str_fit_bench.rs`) draws one motif length's repeat classes from a
known model and fits them. With `--fixture-out` it fits one penalised method on every locus and
writes everything the fit returned as one line of JSON: what it was given (`input`), each class's
numbers and standard errors at full precision, the curves and their strengths, and the round
counts. The fixture is those lines for every cell and seed of
`doc/devel/ng/research/str_fit/compare_D_prior_2026-10-14.md` (§2.1 and §2.2), recorded from the
prototype before the code moved (plan `doc/devel/ng/impl_plan/str_penalised_fit.md`, step C0).
Later steps re-run the cells and compare.

Commands (run from the repository root, through uv; macOS or Linux):

    uv run --no-project python scripts/research/str_fit_fixture.py list
    uv run --no-project python scripts/research/str_fit_fixture.py run --out tmp/fixture_runs \
        [--only p1_n1_d3_pcr p2_n60_d10_pcr_free_s1] [--jobs 4] [--stutter-prior-reads 0]
    uv run --no-project python scripts/research/str_fit_fixture.py collect tmp/fixture_runs \
        src/parameter_estimation/joint/ssr_fit/bench_fixtures/d_tract_prior_reference.jsonl.gz
    uv run --no-project python scripts/research/str_fit_fixture.py compare FIXTURE NEW \
        [--relative 1e-6] [--absolute 1e-12] [--require tolerance|exact]
    uv run --no-project python scripts/research/str_fit_fixture.py selftest

**`run`** writes one `<cell>_s<seed>.json` a run into the directory and appends each run's wall
time and peak memory (from `/usr/bin/time`: `-l` on macOS, `-v` on Linux) to `runs.tsv` there.
`--only` takes whole names: a cell (`p1_n1_d3_pcr`, every seed) or a run (`p1_n1_d3_pcr_s2`).
`--stutter-prior-reads` hands every run a stutter prior of that many reads a class in place of
production's one; the runs keep their method's name, so `0` compares the fit without the prior
with a fixture recorded before it.

- **A run whose output is already there is skipped**, so an interrupted directory can be resumed.
  **An output older than the bench binary is refused**: `run` stops before starting anything and
  names those outputs. Comparing a run made by an older build would judge the wrong code, and
  skipping it with a warning would let that comparison happen anyway. Run into a new directory.
- Cells below 2,000 samples run up to `--jobs` at a time (at most 4).
- **2,000-sample runs run alone, after the others, under a memory guard.** Each run starts in its
  own process group. Every 2 s the guard reads the memory available (`vm_stat` on macOS,
  `MemAvailable` in `/proc/meminfo` on Linux) and appends the reading to `memory.log`. Below
  `--memory-floor-gib` (8 GiB) it kills the whole group (`/usr/bin/time` and the bench together),
  waits until no process of the group is left, and starts no further 2,000-sample run.

Build the bench first:

    cargo build --release --features bench-fixtures --example ng_str_fit_bench

**`collect`** gathers a run directory into a fixture file, one line a run, sorted by truth, motif
length, samples, depth and seed. A target ending in `.gz` is gzipped with no file name and a zero
time in its header, so collecting the same runs again gives the same bytes.

**`compare`** reads the fixture and the new runs (a directory of `.json` files, or a `.jsonl` or
`.jsonl.gz` file) and matches them by cell, seed and method. Two runs with the same cell, seed and
method in one set are an error. For each run it first compares `input`, what the fit was given
(the sample-locus rows with reads, and each class's loci with reads). **A difference there is
reported as "the input moved"** and fails the run, because the draw changed, not the fit. Then it
gives two verdicts:

- **exact**: every number has the same bits (`-0.0` is not `0.0`), and every count is equal and of
  the same type (an integer count written as a float is a count difference).
- **within tolerance**: every count equal, and every number with
  `|new − fixture| ≤ relative × max(|new|, |fixture|) + absolute`. `last_change`, the
  convergence diagnostic, is a difference of two nearly equal numbers and moves far more than any
  estimate when the order of a sum changes; **it is reported but not judged here**. The exact
  verdict still compares it.

Counts are the round counts, the readings of F off pairs, the class-rounds counted, where each
curve was eased, whether each class was fitted, whether a number is present, and the shape of the
line. Times (`seconds`) and the printed configuration (`settings`) are not compared. It prints every
number that is not bit-identical (up to `--show` a run), largest relative difference first.

The exit status is 0 when every run compared passes `--require` (`tolerance` by default, or
`exact`); 1 when a run fails, its input moved, or the new set holds a run the fixture does not; and
2 when a fixture run is missing from the new runs and `--partial` was not given, or a set holds a
duplicate run.

**What the fixture does not exercise.** Every one of its 153 runs fits every class and returns
curves. The branches for an unfitted class (`"fitted": false`) and for a run with no curves
(`"curves": null`) are written by the bench and compared here, but no recorded run takes them.

**`selftest`** checks the memory and peak-memory parsers on captured macOS and Linux output, and
the comparison of single values (bits, types, tolerance).
"""

import argparse
import gzip
import json
import math
import os
import platform
import re
import signal
import struct
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
METHOD = "D_tract_prior"
LOCI_AT_MOST = "100"
MEMORY_FLOOR_GIB = 8.0
MEMORY_POLL_SECONDS = 2.0
# A number compared but not judged by the tolerance verdict (see the module's documentation).
DIAGNOSTIC_NUMBERS = {".last_change"}


def cells():
    """Every run of the compare report's grid: (run name, cell name, bench arguments, samples)."""
    grid = []

    def add(period, samples, depth, truth, seeds, between_locus_share=None):
        variant = truth if between_locus_share is None else f"{truth}_F{between_locus_share}"
        cell = f"p{period}_n{samples}_d{depth}_{variant}"
        for seed in seeds:
            args = ["--period", str(period), "--samples", str(samples), "--depth", str(depth),
                    "--seed", str(seed), "--truth", truth, "--loci-at-most", LOCI_AT_MOST,
                    "--methods", METHOD]
            if between_locus_share is not None:
                args += ["--between-locus-share", str(between_locus_share)]
            grid.append((f"{cell}_s{seed}", cell, args, samples))

    # §2.1, the default truth: 6 seeds at one sample, 3 elsewhere.
    for period in (1, 2):
        for samples, depths, seeds in ((1, (3, 10, 30, 300), range(1, 7)),
                                       (10, (3, 10, 30, 300), range(1, 4)),
                                       (60, (3, 10, 30, 300), range(1, 4)),
                                       (500, (3, 10), range(1, 4)),
                                       (2000, (3,), range(1, 4))):
            for depth in depths:
                add(period, samples, depth, "pcr", seeds)
    # §2.2: PCR-free on both motif lengths, F = 0.4 on two-base runs; one seed at 2,000 samples.
    other = (((1, 3), range(1, 4)), ((1, 300), range(1, 4)), ((60, 10), range(1, 4)),
             ((500, 3), range(1, 4)), ((2000, 3), range(1, 2)))
    for (samples, depth), seeds in other:
        for period in (1, 2):
            add(period, samples, depth, "pcr_free", seeds)
        add(2, samples, depth, "pcr", seeds, between_locus_share=0.4)
    return grid


def bench_binary():
    """The newer of the two builds (CLAUDE.md, "Container vs. host")."""
    found = [REPO / d / "release/examples/ng_str_fit_bench" for d in ("target", "target-container")]
    found = [p for p in found if p.exists()]
    if not found:
        sys.exit("no bench binary; cargo build --release --features bench-fixtures "
                 "--example ng_str_fit_bench")
    return max(found, key=lambda p: p.stat().st_mtime)


# ---- Memory: what is available now, and a run's peak, on macOS and on Linux -------------------

def parse_vm_stat(text):
    """GiB available from macOS `vm_stat`: free, inactive, speculative and purgeable pages."""
    pages, page = {}, 16384
    match = re.search(r"page size of (\d+) bytes", text)
    if match:
        page = int(match.group(1))
    for line in text.splitlines()[1:]:
        key, _, value = line.partition(":")
        try:
            pages[key.strip()] = int(value.strip().rstrip("."))
        except ValueError:
            pass
    free = sum(pages.get(k, 0) for k in ("Pages free", "Pages inactive", "Pages speculative",
                                         "Pages purgeable"))
    return free * page / 2**30


def parse_meminfo(text):
    """GiB available from Linux `/proc/meminfo`: its `MemAvailable` line, in kB."""
    match = re.search(r"^MemAvailable:\s+(\d+)\s*kB", text, re.MULTILINE)
    if not match:
        raise ValueError("no MemAvailable line in /proc/meminfo")
    return int(match.group(1)) / 2**20


def available_gib():
    if platform.system() == "Linux":
        return parse_meminfo(Path("/proc/meminfo").read_text())
    return parse_vm_stat(subprocess.run(["vm_stat"], capture_output=True, text=True).stdout)


def time_flag():
    """`/usr/bin/time`'s flag for the full report: BSD `-l` on macOS, GNU `-v` on Linux."""
    return "-v" if platform.system() == "Linux" else "-l"


def parse_peak_mib(text):
    """A run's peak resident memory in MiB from `/usr/bin/time`'s report, or None.

    BSD `time -l` (macOS) writes bytes: `  49496064  maximum resident set size`. GNU `time -v`
    (Linux) writes kB: `Maximum resident set size (kbytes): 48336`.
    """
    for line in text.splitlines():
        gnu = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", line)
        if gnu:
            return int(gnu.group(1)) / 2**10
        bsd = re.match(r"\s*(\d+)\s+maximum resident set size", line)
        if bsd:
            return int(bsd.group(1)) / 2**20
    return None


def group_alive(group):
    try:
        os.killpg(group, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def kill_group(group, wait_seconds=30.0):
    """Kill every process of a group, then wait until none is left. True when it emptied."""
    try:
        os.killpg(group, signal.SIGKILL)
    except ProcessLookupError:
        return True
    deadline = time.time() + wait_seconds
    while time.time() < deadline:
        if not group_alive(group):
            return True
        time.sleep(0.1)
    return False


# ---- run ---------------------------------------------------------------------------------------

def run_one(bench, out, name, args, log_lock, memory_floor=None):
    """One run of the bench. With `memory_floor`, a guard watches the memory available.

    Returns (message, killed): `killed` is True when the guard stopped the run.
    """
    target = out / f"{name}.json"
    if target.exists() and target.stat().st_size > 0:
        return f"skip {name} (already there)", False
    timing = out / f"{name}.time"
    started = time.time()
    killed = threading.Event()
    lowest = [math.inf]
    with open(out / f"{name}.log", "w") as log:
        # Its own session, so its own process group: the guard kills `time` and the bench together.
        process = subprocess.Popen(["/usr/bin/time", time_flag(), "-o", str(timing), str(bench),
                                    *args, "--fixture-out", str(target)],
                                   stderr=log, start_new_session=True)
        stop = threading.Event()

        def watch():
            with open(out / "memory.log", "a") as memory:
                while not stop.is_set():
                    now = available_gib()
                    lowest[0] = min(lowest[0], now)
                    memory.write(f"{time.strftime('%H:%M:%S')}\t{name}\t{now:.1f}\n")
                    memory.flush()
                    if now < memory_floor:
                        killed.set()
                        emptied = kill_group(process.pid)
                        memory.write(f"KILLED {name} at {now:.1f} GiB available; process group "
                                     f"{'gone' if emptied else 'STILL ALIVE after 30 s'}\n")
                        memory.flush()
                        return
                    stop.wait(MEMORY_POLL_SECONDS)

        watcher = None
        if memory_floor is not None:
            watcher = threading.Thread(target=watch, daemon=True)
            watcher.start()
        try:
            code = process.wait()
        except BaseException:
            # Interrupted (Ctrl-C reaches only this script, not the run's own session): take the
            # run down with it rather than leave it running.
            kill_group(process.pid)
            raise
        stop.set()
        if watcher:
            watcher.join()
        # The bench may outlive `time` by a moment; leave nothing of the group behind.
        if group_alive(process.pid):
            kill_group(process.pid)
    wall = time.time() - started
    peak = parse_peak_mib(timing.read_text()) if timing.exists() else None
    peak_text = "NA" if peak is None else f"{peak:.1f}"
    if killed.is_set():
        status = "KILLED(memory)"
        target.unlink(missing_ok=True)
    else:
        status = "ok" if code == 0 and target.exists() else f"FAILED({code})"
    with log_lock, open(out / "runs.tsv", "a") as runs:
        lowest_text = f"{lowest[0]:.1f}" if memory_floor is not None else "NA"
        runs.write(f"{name}\t{status}\t{wall:.1f}\t{peak_text}\t{lowest_text}\n")
    return f"{status} {name} in {wall:.1f} s, peak {peak_text} MiB", killed.is_set()


def chosen_runs(only):
    grid = cells()
    if not only:
        return grid
    known = {run for run, _, _, _ in grid} | {cell for _, cell, _, _ in grid}
    unknown = [o for o in only if o not in known]
    if unknown:
        sys.exit(f"--only: no cell or run named {unknown}; `list` prints them")
    wanted = set(only)
    return [c for c in grid if c[0] in wanted or c[1] in wanted]


def command_run(options):
    out = Path(options.out)
    bench = Path(options.bench) if options.bench else bench_binary()
    chosen = chosen_runs(options.only)
    if options.stutter_prior_reads is not None:
        extra = ["--stutter-prior-reads", str(options.stutter_prior_reads)]
        chosen = [(run, cell, [*args, *extra], samples) for run, cell, args, samples in chosen]
    built = bench.stat().st_mtime
    stale = [run for run, _, _, _ in chosen
             if (out / f"{run}.json").exists() and (out / f"{run}.json").stat().st_mtime < built]
    if stale:
        sys.exit(f"{len(stale)} outputs in {out} are older than the bench binary {bench}, e.g. "
                 f"{stale[:3]}; they were made by an older build. Run into a new directory.")
    out.mkdir(parents=True, exist_ok=True)
    small = [c for c in chosen if c[3] < 2000]
    big = [c for c in chosen if c[3] >= 2000]
    jobs = max(1, min(options.jobs, 4))
    lock = threading.Lock()
    print(f"bench {bench}; {len(small)} runs below 2,000 samples ({jobs} at a time), "
          f"{len(big)} at 2,000 (alone, killed below {options.memory_floor_gib:g} GiB available)",
          flush=True)
    with ThreadPoolExecutor(jobs) as pool:
        for message, _ in pool.map(lambda c: run_one(bench, out, c[0], c[2], lock), small):
            print(message, flush=True)
    for index, (name, _, args, _) in enumerate(big):
        message, killed = run_one(bench, out, name, args, lock,
                                  memory_floor=options.memory_floor_gib)
        print(message, flush=True)
        if killed:
            left = [n for n, _, _, _ in big[index + 1:]]
            sys.exit(f"memory fell below {options.memory_floor_gib:g} GiB during {name}; it was "
                     f"killed, and the {len(left)} 2,000-sample runs after it were not started")


# ---- collect -----------------------------------------------------------------------------------

def order_key(line):
    run = line["run"]
    return (run["truth_variant"], run["period"], run["samples"], run["depth"], run["seed"])


def run_key(line):
    return (line["run"]["cell"], line["run"]["seed"], line["method"])


def read_runs(path):
    """Every run in a directory of `.json` files or a `.jsonl` / `.jsonl.gz` file, by its key."""
    path = Path(path)
    if path.is_dir():
        lines = [json.loads(file.read_text()) for file in sorted(path.glob("*.json"))]
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            lines = [json.loads(text) for text in handle if text.strip()]
    runs, duplicates = {}, []
    for line in lines:
        key = run_key(line)
        if key in runs:
            duplicates.append(key)
        runs[key] = line
    if duplicates:
        print(f"error: {path} holds {len(duplicates)} runs twice (same cell, seed and method), "
              f"e.g. {[label(k) for k in duplicates[:3]]}", file=sys.stderr)
        sys.exit(2)
    return runs


def command_collect(options):
    runs = read_runs(options.runs)
    expected = {run for run, _, _, _ in cells()}
    have = {f"{cell}_s{seed}" for cell, seed, _ in runs}
    missing = sorted(expected - have)
    if missing and not options.partial:
        sys.exit(f"{len(missing)} runs missing, e.g. {missing[:5]}; --partial to collect anyway")
    target = Path(options.fixture)
    text = "".join(json.dumps(line, sort_keys=True, allow_nan=False) + "\n"
                   for line in sorted(runs.values(), key=order_key)).encode()
    partial = target.with_name(target.name + ".partial")
    with open(partial, "wb") as raw:
        if target.suffix == ".gz":
            # No file name and a zero time in the header: the same runs give the same bytes.
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as packed:
                packed.write(text)
        else:
            raw.write(text)
    os.replace(partial, target)
    print(f"{len(runs)} runs into {target} ({target.stat().st_size / 1e6:.2f} MB; "
          f"{len(text) / 1e6:.2f} MB unpacked)")


# ---- compare -----------------------------------------------------------------------------------

SKIPPED = {"seconds", "settings", "input"}


def same_bits(a, b):
    """The exact verdict on one value: same type, and for a float the same 64 bits."""
    if type(a) is not type(b):
        return False
    if isinstance(a, float):
        return struct.pack("<d", a) == struct.pack("<d", b)
    return a == b


def is_number(value):
    return isinstance(value, float)


def walk(fixture, new, path, counts, numbers):
    """Pair up every leaf of two runs: numbers (both floats) and counts (everything else).

    A count is an integer, a boolean, a string (`"NaN"`), `null`, or the shape of the line. A
    float paired with anything but a float is a count difference: an integer written as a float,
    or a number that appeared or disappeared.
    """
    if isinstance(fixture, dict) and isinstance(new, dict):
        for key in sorted(set(fixture) | set(new)):
            if key in SKIPPED and not path:
                continue
            if key not in fixture or key not in new:
                counts.append((f"{path}.{key}", fixture.get(key, "<missing>"),
                               new.get(key, "<missing>")))
            else:
                walk(fixture[key], new[key], f"{path}.{key}", counts, numbers)
    elif isinstance(fixture, list) and isinstance(new, list):
        if len(fixture) != len(new):
            counts.append((f"{path}.length", len(fixture), len(new)))
        for index, (a, b) in enumerate(zip(fixture, new)):
            walk(a, b, f"{path}[{index}]", counts, numbers)
    elif is_number(fixture) and is_number(new):
        numbers.append((path, fixture, new))
    elif not same_bits(fixture, new):
        counts.append((path, fixture, new))


def relative_difference(a, b):
    if a == b:
        return 0.0
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else math.inf


def within_tolerance(a, b, relative, absolute):
    """`|a − b| ≤ relative × max(|a|, |b|) + absolute`; equal values (also ±0) always pass."""
    return a == b or abs(a - b) <= relative * max(abs(a), abs(b)) + absolute


def input_differences(fixture, new):
    """The input a run was given, compared: a list of differences, or None when not recorded."""
    if "input" not in fixture or "input" not in new:
        return None
    counts, numbers = [], []
    walk(fixture["input"], new["input"], ".input", counts, numbers)
    return counts + [(p, a, b) for p, a, b in numbers if not same_bits(a, b)]


def label(run):
    return f"{run[0]} seed {run[1]} {run[2]}"


def command_compare(options):
    fixture = read_runs(options.fixture)
    new = read_runs(options.new)
    missing = sorted(set(fixture) - set(new))
    extra = sorted(set(new) - set(fixture))
    failed = moved = unrecorded = 0
    exact_runs = tolerance_runs = 0
    largest = (0.0, None)
    for key in sorted(set(fixture) & set(new)):
        moved_input = input_differences(fixture[key], new[key])
        if moved_input is None:
            unrecorded += 1
        elif moved_input:
            moved += 1
            print(f"{label(key)}: THE INPUT MOVED — the draw changed, not only the fit; "
                  f"{len(moved_input)} input values differ")
            for path, a, b in moved_input[: options.show]:
                print(f"    input {path}: fixture {a!r}, new {b!r}")
        counts, numbers = [], []
        walk(fixture[key], new[key], "", counts, numbers)
        differing = [(relative_difference(a, b), p, a, b) for p, a, b in numbers
                     if not same_bits(a, b)]
        differing.sort(key=lambda d: -d[0])
        over = [d for d in differing if d[1] not in DIAGNOSTIC_NUMBERS
                and not within_tolerance(d[2], d[3], options.relative, options.absolute)]
        exact = not moved_input and not counts and not differing
        within = not moved_input and not counts and not over
        exact_runs += exact
        tolerance_runs += within
        if differing and differing[0][0] > largest[0]:
            largest = (differing[0][0], (label(key), differing[0][1]))
        passed = exact if options.require == "exact" else within
        failed += not passed
        verdict = ("exact" if exact else "within tolerance" if within else
                   "INPUT MOVED" if moved_input else "OUTSIDE tolerance")
        print(f"{label(key)}: {verdict}; {len(numbers)} numbers, {len(differing)} not "
              f"bit-identical, {len(over)} outside the tolerance; counts differing: {len(counts)}")
        for path, a, b in counts[: options.show]:
            print(f"    count {path}: fixture {a!r}, new {b!r}")
        for rel, path, a, b in differing[: options.show]:
            note = " (diagnostic, not judged by the tolerance)" if path in DIAGNOSTIC_NUMBERS else ""
            print(f"    {path}: fixture {a!r}, new {b!r}, relative {rel:.3g}{note}")
    compared = len(set(fixture) & set(new))
    print(f"\n{compared} runs compared: {exact_runs} exact, {tolerance_runs} within "
          f"{options.relative:g} relative + {options.absolute:g} with every count equal; "
          f"{failed} fail --require {options.require}")
    if moved:
        print(f"{moved} runs were given a different input: the draw moved, not only the fit")
    if unrecorded:
        print(f"{unrecorded} runs have no input recorded on one side; their input was not checked")
    if largest[1]:
        print(f"largest relative difference {largest[0]:.3g} at {largest[1][0]} {largest[1][1]}")
    if extra:
        print(f"{len(extra)} new runs not in the fixture, e.g. {[label(k) for k in extra[:3]]}")
    if missing:
        print(f"{len(missing)} fixture runs not in the new runs, e.g. "
              f"{[label(k) for k in missing[:3]]}")
    if failed or extra:
        sys.exit(1)
    if missing and not options.partial:
        sys.exit(2)


# ---- selftest ----------------------------------------------------------------------------------

VM_STAT_SAMPLE = """\
Mach Virtual Memory Statistics: (page size of 16384 bytes)
Pages free:                               10000.
Pages active:                            500000.
Pages inactive:                           20000.
Pages speculative:                         3000.
Pages throttled:                              0.
Pages wired down:                        200000.
Pages purgeable:                           1000.
"""

MEMINFO_SAMPLE = """\
MemTotal:       131895132 kB
MemFree:         2034816 kB
MemAvailable:   97534212 kB
Buffers:          412300 kB
Cached:         89123456 kB
"""

BSD_TIME_SAMPLE = """\
        1.36 real         3.28 user         2.58 sys
            49496064  maximum resident set size
                   0  average shared memory size
            27656744  peak memory footprint
"""

GNU_TIME_SAMPLE = """\
\tCommand being timed: "ng_str_fit_bench --period 1"
\tUser time (seconds): 3.28
\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:01.36
\tAverage resident set size (kbytes): 0
\tMaximum resident set size (kbytes): 48336
\tExit status: 0
"""


def command_selftest(_options):
    checks = [
        ("vm_stat: free+inactive+speculative+purgeable pages × 16 KiB",
         math.isclose(parse_vm_stat(VM_STAT_SAMPLE), 34000 * 16384 / 2**30)),
        ("/proc/meminfo: MemAvailable kB", math.isclose(parse_meminfo(MEMINFO_SAMPLE),
                                                        97534212 / 2**20)),
        ("time -l (bytes)", math.isclose(parse_peak_mib(BSD_TIME_SAMPLE), 49496064 / 2**20)),
        ("time -v (kB)", math.isclose(parse_peak_mib(GNU_TIME_SAMPLE), 48336 / 2**10)),
        ("no peak line", parse_peak_mib("0.1 real") is None),
        ("-0.0 is not 0.0 bit for bit", not same_bits(-0.0, 0.0)),
        ("-0.0 and 0.0 within tolerance", within_tolerance(-0.0, 0.0, 1e-6, 1e-12)),
        ("int count vs the same float is not the same", not same_bits(133, 133.0)),
        ("true is not 1", not same_bits(True, 1)),
        ("same float, same bits", same_bits(0.1 + 0.2, 0.30000000000000004)),
        ("0.0 vs 1e-17 within tolerance (absolute floor)", within_tolerance(0.0, 1e-17, 1e-6,
                                                                            1e-12)),
        ("0.0 vs 1e-11 outside tolerance", not within_tolerance(0.0, 1e-11, 1e-6, 1e-12)),
        ("5e-7 relative within", within_tolerance(1.0, 1.0 + 5e-7, 1e-6, 1e-12)),
        ("2e-6 relative outside", not within_tolerance(1.0, 1.0 + 2e-6, 1e-6, 1e-12)),
    ]
    counts, numbers = [], []
    walk({"rounds": 133, "x": 1.0, "e": None}, {"rounds": 133.0, "x": None, "e": 1.0}, "",
         counts, numbers)
    checks.append(("walk sends int/float and float/null pairs to counts",
                   len(counts) == 3 and not numbers))
    failed = [name for name, ok in checks if not ok]
    for name, ok in checks:
        print(f"{'ok  ' if ok else 'FAIL'} {name}")
    if failed:
        sys.exit(1)
    print(f"{len(checks)} checks passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="print every run of the grid")
    run = commands.add_parser("run", help="run the grid's cells with the bench")
    run.add_argument("--out", required=True)
    run.add_argument("--only", nargs="*",
                     help="only these cells (every seed) or runs, by whole name")
    run.add_argument("--jobs", type=int, default=4)
    run.add_argument("--bench", help="the bench binary, in place of the newest build")
    run.add_argument("--memory-floor-gib", type=float, default=MEMORY_FLOOR_GIB,
                     help="kill a 2,000-sample run when available memory falls below this")
    run.add_argument("--stutter-prior-reads", type=float,
                     help="the stutter prior's weight in reads a class, in place of production's; "
                          "the runs keep the method's name, so 0 compares the fit without the "
                          "prior with a fixture recorded before it")
    collect = commands.add_parser("collect", help="gather a run directory into a fixture file")
    collect.add_argument("runs")
    collect.add_argument("fixture")
    collect.add_argument("--partial", action="store_true")
    compare = commands.add_parser("compare", help="compare new runs with the fixture")
    compare.add_argument("fixture")
    compare.add_argument("new")
    compare.add_argument("--relative", type=float, default=1e-6)
    compare.add_argument("--absolute", type=float, default=1e-12)
    compare.add_argument("--require", choices=("tolerance", "exact"), default="tolerance")
    compare.add_argument("--show", type=int, default=5)
    compare.add_argument("--partial", action="store_true",
                         help="fixture runs missing from the new runs are not an error")
    commands.add_parser("selftest", help="check the parsers and the value comparison")
    options = parser.parse_args()
    if options.command == "list":
        for name, _, args, _ in cells():
            print(name, " ".join(args))
    elif options.command == "run":
        command_run(options)
    elif options.command == "collect":
        command_collect(options)
    elif options.command == "selftest":
        command_selftest(options)
    else:
        command_compare(options)


if __name__ == "__main__":
    main()
