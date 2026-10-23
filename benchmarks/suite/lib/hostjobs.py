"""Steps that run in another tool's official image, through the host job runner (host_jobs.sh).

A step writes a request into `<results>/host_jobs/` and waits for the runner's answer. Where the
runner cannot run the image (this Mac; any machine that is not x86-64 with podman or Docker), the
answer says why and the step is skipped — listed at the end of the run, not counted as a failure.
"""

import json
import os
import time
import uuid

from .measure import CommandFailed
from .steps import SkipStep

DEEPVARIANT_VERSION = "1.10.0"
GLNEXUS_VERSION = "1.4.1"


def _dir(ctx):
    return os.path.join(ctx.results_root, "host_jobs")


def capabilities(ctx):
    """(available, reason) as the host job runner found the machine."""
    path = os.path.join(_dir(ctx), "capabilities.json")
    if not os.path.exists(path):
        return False, "the host job runner is not running (start the suite with run_suite.sh)"
    with open(path) as fh:
        cap = json.load(fh)
    return bool(cap.get("available")), cap.get("reason", "")


def require(ctx):
    available, reason = capabilities(ctx)
    if not available:
        raise SkipStep(f"DeepVariant and GLnexus cannot run here: {reason}")


def skip_with(ctx, callset):
    """A step that scores a DeepVariant call set is skipped where DeepVariant was."""
    if callset.arm == "deepvariant":
        require(ctx)


def run(ctx, image, command, workdir, label, threads=None):
    """Run `command` in `image` ("deepvariant" or "glnexus") with `workdir` as its working
    directory; return the suite's measurement of it, or None when the image could not measure it."""
    require(ctx)
    jobs = _dir(ctx)
    os.makedirs(jobs, exist_ok=True)
    os.makedirs(os.path.join(workdir, "logs"), exist_ok=True)
    job = uuid.uuid4().hex[:12]
    result = os.path.join(workdir, f"{label}.measurement.json")
    request = {"image": image, "workdir": workdir, "result": result,
               "threads": threads or ctx.threads, "label": label, "command": command}
    with open(os.path.join(workdir, "logs", f"{label}.cmd"), "w") as fh:
        fh.write(f"[{image} image, through the host job runner] {command}\n")
    tmp = os.path.join(jobs, f"{job}.request.tmp")
    with open(tmp, "w") as fh:
        json.dump(request, fh)
    os.replace(tmp, os.path.join(jobs, f"{job}.request"))
    answer_path = os.path.join(jobs, f"{job}.answer")
    while not os.path.exists(answer_path):
        time.sleep(2)
    with open(answer_path) as fh:
        answer = json.load(fh)
    os.remove(answer_path)
    if answer.get("unavailable"):
        raise SkipStep(answer.get("error", "unavailable"))
    if not answer.get("ok"):
        log = answer.get("runner_log", "")
        tail = ""
        if log and os.path.exists(log):
            with open(log) as fh:
                tail = fh.read()[-1500:]
        raise CommandFailed(f"{label} in the {image} image failed (exit {answer.get('exit')}): "
                            f"{answer.get('error', '')} {tail}")
    if os.path.exists(result):
        with open(result) as fh:
            return json.load(fh).get("measurement")
    return None
