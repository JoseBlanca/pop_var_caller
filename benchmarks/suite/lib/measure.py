"""Run commands and measure their wall time, CPU time and peak memory.

The dev image has no GNU `time`, so this does the measuring itself, from two sources:

- **When a process ends**, the operating system reports its CPU time exactly (`os.wait4`), for the
  process and every descendant it waited for, so a shell pipeline is measured whole.
- **The largest single process's peak memory** is each process's own high-water mark (`VmHWM` in
  `/proc/<pid>/status`), sampled while it runs. Not `wait4`'s peak, which on Linux also counts the
  memory of the process that started it, as it was when it started it: with the suite's Python at
  817 MB, every command, a samtools count included, reported at least 817 MB.
- **While processes run**, a sampler adds up the resident memory of every live process the step
  started, five times a second. That sum is what a stage of several processes at once really
  costs; no single process's peak shows it.

Several commands can run at once (`max_concurrent`), which is how the one-sample stages are run:
one process per sample, as many at a time as the suite has cores.
"""

import os
import signal
import subprocess
import threading
import time

import psutil

SAMPLE_INTERVAL_S = 0.2


class CommandFailed(Exception):
    pass


def _high_water(pid):
    """A process's peak resident memory so far, in bytes (Linux `VmHWM`), or None."""
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        return None
    return None


class _MemorySampler(threading.Thread):
    """Adds up the resident memory of the live processes and all their descendants."""

    def __init__(self):
        super().__init__(daemon=True)
        self.pids = set()
        self.lock = threading.Lock()
        self.peak_summed_bytes = 0
        self.high_water = {}
        self.stopping = threading.Event()

    def add(self, pid):
        with self.lock:
            self.pids.add(pid)

    def remove(self, pid):
        with self.lock:
            self.pids.discard(pid)

    def run(self):
        while not self.stopping.is_set():
            with self.lock:
                pids = list(self.pids)
            total = 0
            for pid in pids:
                try:
                    root = psutil.Process(pid)
                    procs = [root] + root.children(recursive=True)
                except psutil.Error:
                    continue
                for proc in procs:
                    try:
                        total += proc.memory_info().rss
                    except psutil.Error:
                        pass
                    hwm = _high_water(proc.pid)
                    if hwm is not None and hwm > self.high_water.get(proc.pid, 0):
                        self.high_water[proc.pid] = hwm
            self.peak_summed_bytes = max(self.peak_summed_bytes, total)
            self.stopping.wait(SAMPLE_INTERVAL_S)


def run_commands(commands, log_dir, max_concurrent=1, cwd=None, env=None):
    """Run shell commands, at most `max_concurrent` at a time, and measure them together.

    `commands` is a list of (label, shell command). Each command's standard output and standard
    error go to `<log_dir>/<label>.stdout` and `.stderr`, and the command line itself to
    `<label>.cmd`; a command that writes data on standard output redirects it itself.

    Returns a dict: wall seconds from the first start to the last finish, user and system CPU
    seconds summed over every command, the largest single process's peak resident memory, and the
    peak of the summed resident memory of everything alive at once. Raises CommandFailed, after
    stopping the rest, if any command exits non-zero. Ctrl-C or SIGTERM stops every command.
    """
    os.makedirs(log_dir, exist_ok=True)
    pending = list(commands)
    running = {}
    finished = []
    sampler = _MemorySampler()
    sampler.start()
    start = time.time()
    failure = None
    try:
        while pending or running:
            while pending and len(running) < max_concurrent and failure is None:
                label, cmd = pending.pop(0)
                with open(os.path.join(log_dir, f"{label}.cmd"), "w") as fh:
                    fh.write(cmd + "\n")
                out = open(os.path.join(log_dir, f"{label}.stdout"), "wb")
                err = open(os.path.join(log_dir, f"{label}.stderr"), "wb")
                proc = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd], stdout=out,
                                        stderr=err, cwd=cwd, env=env, start_new_session=True)
                out.close()
                err.close()
                running[proc.pid] = (label, proc, time.time())
                sampler.add(proc.pid)
            if failure is not None and not running:
                break
            if failure is not None and pending:
                pending.clear()
            pid, status, usage = os.wait4(-1, 0)
            if pid not in running:
                continue
            label, proc, started = running.pop(pid)
            sampler.remove(pid)
            proc.returncode = os.waitstatus_to_exitcode(status)
            finished.append({
                "label": label,
                "exit": proc.returncode,
                "wall_s": time.time() - started,
                "user_s": usage.ru_utime,
                "sys_s": usage.ru_stime,
                # Linux reports ru_maxrss in kilobytes.
                "max_rss_bytes": usage.ru_maxrss * 1024,
            })
            if proc.returncode != 0 and failure is None:
                failure = f"{label} exited with {proc.returncode}; see {log_dir}/{label}.stderr"
                _stop_all(running)
    except BaseException:
        _stop_all(running)
        raise
    finally:
        sampler.stopping.set()
        sampler.join()
    if failure is not None:
        raise CommandFailed(failure)
    return {
        "wall_s": round(time.time() - start, 3),
        "user_s": round(sum(f["user_s"] for f in finished), 3),
        "sys_s": round(sum(f["sys_s"] for f in finished), 3),
        # A process shorter than one sampling interval (0.2 s) has no high-water sample; its peak
        # is then too small to matter for any figure quoted.
        "max_single_rss_bytes": max(sampler.high_water.values(), default=0),
        "peak_summed_rss_bytes": max(sampler.peak_summed_bytes,
                                     max(sampler.high_water.values(), default=0)),
        "processes": len(finished),
        "max_concurrent": max_concurrent,
        "commands": finished,
    }


def _stop_all(running):
    for _, proc, _ in running.values():
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.time() + 10
    for _, proc, _ in list(running.values()):
        try:
            proc.wait(timeout=max(0.1, deadline - time.time()))
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
    running.clear()


def merge(*measurements):
    """One measurement for several stages run one after another: times add, memory peaks max."""
    measurements = [m for m in measurements if m]
    return {
        "wall_s": round(sum(m["wall_s"] for m in measurements), 3),
        "user_s": round(sum(m["user_s"] for m in measurements), 3),
        "sys_s": round(sum(m["sys_s"] for m in measurements), 3),
        "max_single_rss_bytes": max((m["max_single_rss_bytes"] for m in measurements), default=0),
        "peak_summed_rss_bytes": max((m["peak_summed_rss_bytes"] for m in measurements), default=0),
    }


def check_alignments(paths):
    """Refuse a damaged alignment file before any tool reads it.

    `samtools quickcheck` reads each file's header and checks that it ends with its end-of-file
    marker. Some tools read a truncated CRAM silently up to the break (pvc and GATK did, with a
    tomato CRAM missing its last third), others abort (freebayes), so a damaged input would give
    each caller different data without anyone being told.
    """
    files = [p for p in paths if p.endswith((".bam", ".cram"))]
    if not files:
        return
    out = subprocess.run(["samtools", "quickcheck", "-v", *files], capture_output=True, text=True)
    if out.returncode != 0:
        bad = sorted(set(l.strip() for l in out.stdout.splitlines() if l.strip()))
        raise CommandFailed("damaged alignment file(s), refused before reading: " + ", ".join(bad)
                            + (f" ({out.stderr.strip()})" if out.stderr.strip() else ""))


def warm(paths):
    """Read each file once, untimed, so that the timed step after it reads from memory.

    Every timed step that reads alignment files warms them first (with their indexes and the
    reference), so every tool starts from the same cache state: the files already in memory.
    """
    import glob

    files = []
    for p in paths:
        for f in [p] + glob.glob(p + ".*i") + glob.glob(p + ".fai"):
            if os.path.isfile(f) and f not in files:
                files.append(f)
    for f in files:
        with open(f, "rb") as fh:
            while fh.read(1 << 24):
                pass
