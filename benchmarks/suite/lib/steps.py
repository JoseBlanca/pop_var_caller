"""The step runner: what makes a run stoppable and resumable.

A run is a list of small steps, each one caller's one stage on one configuration, or one scoring
pass. A step owns one directory. It is finished when that directory holds a `.done` file, written
last, carrying the step's measurement, its result rows and its fingerprint. So:

- a resumed run skips every step with a `.done`, and reads its result from there;
- a directory without `.done` is a step that was cut short; it is deleted and the step starts over;
- `benchmark_results/STOP`, made by `run_suite.sh --stop`, is checked between steps, and the run
  ends cleanly before starting the next one.

**The fingerprint** is what the step's output depends on: the tool's version, the command, and the
inputs' paths, sizes and modification times. What a mismatch does depends on the step's kind:

| kind         | lives under            | on a fingerprint mismatch                                   |
|--------------|------------------------|-------------------------------------------------------------|
| `input`      | `inputs/`              | stop and say what changed (they are made once)              |
| `competitor` | `competitors/`         | stop and say what changed, unless `--rerun-competitors`     |
| `pvc`        | `pvc/<commit>/`        | rerun                                                       |
| `score`      | beside what it scores  | rerun                                                       |

A competitor is rerun only when asked for, because a competitor rerun costs hours and its calls do
not change unless its version or its input does.
"""

import hashlib
import json
import os
import shutil
import time

DONE = ".done"
SMALL_FILE_BYTES = 64 * 1024 * 1024


class SuiteError(Exception):
    pass


class StopRequested(Exception):
    """The run ends cleanly between steps: `--stop` was asked for, or `--stop-at` has passed."""


def file_identity(path):
    """What a fingerprint records about one input file: a checksum if small, else size and time."""
    path = os.path.realpath(path)
    if os.path.isdir(path):
        return {"path": path, "dir": sorted(os.listdir(path))}
    st = os.stat(path)
    if st.st_size <= SMALL_FILE_BYTES:
        digest = hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                digest.update(block)
        return {"path": path, "sha256": digest.hexdigest()}
    return {"path": path, "size": st.st_size, "mtime": int(st.st_mtime)}


class Step:
    """One unit of work.

    `name` is unique within a run and is also the step's directory below its kind's root.
    `action(step_dir)` does the work inside `step_dir` and returns a dict with any of
    `measurement` (from measure.run_commands), `rows` (result rows, see results.py) and `meta`.
    `attrs` carries what the selection options filter on: benchmark, depth, caller.
    """

    def __init__(self, name, kind, action, fingerprint, attrs=None, title=None):
        assert kind in ("input", "competitor", "pvc", "score")
        self.name = name
        self.kind = kind
        self.action = action
        self.fingerprint = fingerprint
        self.attrs = attrs or {}
        self.title = title or name


class Runner:
    def __init__(self, ctx):
        self.ctx = ctx
        self.done_count = 0
        self.skipped_count = 0

    def step_dir(self, step):
        return self.ctx.step_dir(step)

    def load_done(self, step):
        path = os.path.join(self.step_dir(step), DONE)
        if not os.path.exists(path):
            return None
        with open(path) as fh:
            return json.load(fh)

    def is_done(self, step):
        return self.load_done(step) is not None

    def run(self, step):
        """Run `step` unless it is already done with a matching fingerprint; return its record."""
        directory = self.step_dir(step)
        fingerprint = _canonical(step.fingerprint() if callable(step.fingerprint) else step.fingerprint)
        record = self.load_done(step)
        if record is not None:
            if (record.get("fingerprint") == fingerprint and not self._forced(step)
                    and not self._retime(step, record)):
                self.skipped_count += 1
                return record
            if not self._forced(step) and not self._retime(step, record):
                if step.kind == "input" or step.kind == "competitor":
                    changed = _difference(record.get("fingerprint"), fingerprint)
                    hint = ("rerun with --rerun-competitors to remake it"
                            if step.kind == "competitor" else "delete the directory to remake it")
                    raise SuiteError(
                        f"{step.name}: its stored output was made from different inputs "
                        f"({changed}); {hint}: {directory}")
        self._check_stop()
        if step.kind in ("pvc", "competitor") and self.ctx.wait_for_quiet:
            self._wait_for_quiet()
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        started = time.time()
        self.ctx.log(f"[{time.strftime('%H:%M:%S')}] {step.kind:10s} {step.name}")
        try:
            result = step.action(directory) or {}
        except BaseException:
            self.ctx.log(f"    failed or stopped; {directory} will be redone next time")
            raise
        ended = time.time()
        record = {
            "name": step.name,
            "kind": step.kind,
            "attrs": step.attrs,
            "fingerprint": fingerprint,
            "started": started,
            "ended": ended,
            "date": time.strftime("%Y-%m-%d", time.localtime(started)),
            "machine": self.ctx.machine,
            "measurement": result.get("measurement"),
            "rows": result.get("rows", []),
            "meta": result.get("meta", {}),
        }
        tmp = os.path.join(directory, DONE + ".tmp")
        with open(tmp, "w") as fh:
            json.dump(record, fh, indent=1, default=str)
        os.replace(tmp, os.path.join(directory, DONE))
        self.done_count += 1
        if result.get("measurement"):
            m = result["measurement"]
            self.ctx.log(f"    {m['wall_s']:.1f} s, peak {m['peak_summed_rss_bytes'] / 2**20:.0f} MB")
        return record

    def _check_stop(self):
        if os.path.exists(os.path.join(self.ctx.results_root, "STOP")):
            raise StopRequested("stopped as asked (run_suite.sh --stop)")
        if self.ctx.stop_at and time.time() >= self.ctx.stop_at:
            raise StopRequested("stopped: the --stop-at time has passed")

    def _wait_for_quiet(self):
        """Hold a timed step until the machine is quiet.

        Quiet: for the last `quiet_minutes`, processes other than the suite used under
        `quiet_cores` on average in every 5-s sample run_suite.sh recorded. Checked every 30 s;
        `--stop` and `--stop-at` still end the run while it waits. Without the host's samples (the
        suite started some other way) it does not wait.
        """
        from .results import HostLoad

        path = os.path.join(self.ctx.results_root, "host_load.tsv")
        logged = 0.0
        while True:
            load = HostLoad(path)
            now = time.time()
            if not load.times or now - load.times[-1] > 60:
                return
            window = [o for t, o in zip(load.times, load.other)
                      if t >= now - 60 * self.ctx.quiet_minutes and o is not None]
            enough = len(window) >= max(3, int(self.ctx.quiet_minutes * 60 / 5 * 0.8))
            if enough and max(window) < self.ctx.quiet_cores:
                if logged:
                    self.ctx.log(f"    quiet now; waited {(now - logged) / 60:.0f} min")
                return
            if not logged or now - logged > 600:
                if not enough:
                    why = "the load has not been sampled for long enough yet"
                else:
                    why = f"other processes are using {window[-1]:.1f} cores"
                self.ctx.log(f"[{time.strftime('%H:%M:%S')}] waiting for a quiet machine ({why}; "
                             f"need under {self.ctx.quiet_cores:g} for {self.ctx.quiet_minutes:g} min)")
                logged = logged or now
            time.sleep(30)
            self._check_stop()

    def _retime(self, step, record):
        """With --redo-busy, a pvc or competitor step whose time was taken on a busy machine."""
        if not self.ctx.redo_busy or step.kind not in ("pvc", "competitor"):
            return False
        if not record.get("measurement"):
            return False
        from .results import BUSY_OTHER_CPU, HostLoad

        load = HostLoad(os.path.join(self.ctx.results_root, "host_load.tsv"))
        other = load.other_mean(record["started"], record["ended"])
        return other != "" and other > BUSY_OTHER_CPU

    def _forced(self, step):
        if step.name in self.ctx.redo_steps:
            return True
        if step.kind == "competitor" and self.ctx.rerun_competitors:
            return True
        if step.kind == "score" and self.ctx.rescore:
            return True
        return False


def _canonical(value):
    return json.loads(json.dumps(value, sort_keys=True))


def _difference(old, new, prefix=""):
    """A short description of the first place two fingerprints differ."""
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new)):
            if old.get(key) != new.get(key):
                return _difference(old.get(key), new.get(key), f"{prefix}{key}.")
    if isinstance(old, list) and isinstance(new, list) and len(old) == len(new):
        for i, (a, b) in enumerate(zip(old, new)):
            if a != b:
                return _difference(a, b, f"{prefix}{i}.")
    return f"{prefix.rstrip('.') or 'fingerprint'}: {old!r} -> {new!r}"
