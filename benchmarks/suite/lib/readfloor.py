"""The read floor: the cost of only reading and decoding the alignment files.

No caller can be faster than reading its input, so set beside pvc's step that reads the alignments
it shows how much of pvc's time is pvc's own work. The command, per alignment file:

    samtools view -c -M -L regions.bed -@ <extra threads> file

`-c` decompresses and decodes every record but formats and writes nothing (`-o /dev/null` would
format every record as SAM text, work pvc never does). `-L` with `-M` reads exactly the regions pvc
calls, each read once even where regions overlap.

Two ways, each matching one pvc step (results.READ_FLOOR_PAIRS):

- **one process per sample**, as many at once as the run's cores, against `generate-psps`, which is
  run that way. Each process reads with the threads `generate-psps` uses for one file: one for a
  BAM, two for a CRAM (pvc decodes CRAM containers ahead on a second thread), so `-@ 0` and `-@ 1`;
  samtools' `-@` counts threads beyond the first.
- **every file in turn, at the run's thread count** (`-@ threads-1`), summed, against
  `call-from-alignments`.

**CRAM:** samtools 1.16 decodes every field of a CRAM record for a count — measured on the HG002
5 Mb CRAM: 15.4 s by default, 15.4 s forcing every field, 1.3 s asking for flags only. pvc decodes
every field except the auxiliary tags other than the read group, so a CRAM floor asks for exactly
those (`required_fields=0x17ff`). BAM has no such shortcut: every record is decoded whole.

Run once per dataset and machine, like a competitor, and reused while samtools' version and the
input are unchanged.
"""

import os
import shlex

from .measure import run_commands, warm
from .results import row, time_rows
from .steps import Step, file_identity

q = shlex.quote
# FLAG, QNAME, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL and the read group: everything
# pvc reads from a record (src/read/aligned_read.rs, container.rs), and no other auxiliary tag.
PVC_CRAM_FIELDS = "0x17ff"


def _command(path, reference, regions, extra_threads, out):
    cram = path.endswith(".cram")
    fields = f"--input-fmt-option required_fields={PVC_CRAM_FIELDS} " if cram else ""
    return (f"samtools view -c -M -L {q(regions)} -T {q(reference)} {fields}"
            f"-@ {extra_threads} {q(path)} > {q(out)}")


def floor_steps(ctx, prefix, reference, alignments_by_sample, regions, base):
    """The two read-floor steps for one configuration."""
    files = [(s, f) for s, v in alignments_by_sample.items()
             for f in (v if isinstance(v, list) else [v])]
    row_base = dict(base, caller="samtools (read floor)", caller_version=ctx.version("samtools"))

    def fingerprint(mode):
        return lambda: {"tool": ctx.version("samtools"), "mode": mode, "threads": ctx.threads,
                        "files": [file_identity(f) for _, f in files],
                        "regions": file_identity(regions),
                        "reference": file_identity(reference + ".fai"),
                        "cram_fields": PVC_CRAM_FIELDS}

    def per_sample(d):
        warm([f for _, f in files] + [reference])
        cmds = [(f"floor.{s}.{i}",
                 _command(f, reference, regions, 1 if f.endswith(".cram") else 0,
                          os.path.join(d, f"{s}.{i}.count")))
                for i, (s, f) in enumerate(files)]
        m = run_commands(cmds, os.path.join(d, "logs"), max_concurrent=ctx.threads)
        return {"measurement": m,
                "rows": time_rows(row_base, m, "read floor, one process per sample")
                + [row(row_base, stage="read floor, one process per sample",
                       metric="records_read", value=_records(d))]}

    def at_threads(d):
        warm([f for _, f in files] + [reference])
        cmds = [(f"floor.{s}.{i}",
                 _command(f, reference, regions, max(0, ctx.threads - 1),
                          os.path.join(d, f"{s}.{i}.count")))
                for i, (s, f) in enumerate(files)]
        m = run_commands(cmds, os.path.join(d, "logs"), max_concurrent=1)
        return {"measurement": m,
                "rows": time_rows(row_base, m, "read floor, every file at the run's threads")
                + [row(row_base, stage="read floor, every file at the run's threads",
                       metric="records_read", value=_records(d))]}

    attrs = {"caller": "samtools", "depth": base.get("depth")}
    return [Step(f"{prefix}/read_floor_per_sample", "competitor", per_sample,
                 fingerprint("per sample"), attrs=attrs),
            Step(f"{prefix}/read_floor_at_threads", "competitor", at_threads,
                 fingerprint("at threads"), attrs=attrs)]


def _records(d):
    total = 0
    for f in os.listdir(d):
        if f.endswith(".count"):
            with open(os.path.join(d, f)) as fh:
                total += int(fh.read().strip() or 0)
    return total
