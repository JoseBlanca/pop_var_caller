#!/usr/bin/env bash
# Survey where tracts start to stutter, across many libraries — run this ON RICK.
#
# ng routes a locus to the STR path when it is **likely to stutter**, not merely when it contains a
# repeat (doc/devel/ng/spec/parameter_prepass_ssr.md §5.1), and the per-period copy floors are where
# that line is drawn. Whether a tract stutters is a property of the **library** — PCR amplification
# stutters more than a PCR-free preparation — and every number ng has today comes from a single
# library per species. This runs the measurement across the archive so the floors rest on the axis
# that actually drives them.
#
# Usage:
#   ./ng_str_library_survey.sh OUT.tsv REF CRAM [CRAM ...]
#
# Example on rick — **rebuild first**, in the container, or the run refuses to start:
#   ./scripts/dev.sh cargo build --release --example ng_str_stutter_by_library
#   JOBS=32 ./scripts/ng_str_library_survey.sh ~/tmp/stutter_by_library.tsv \
#       /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa \
#       /media/tomato25_bams/crams/*/*.cram
#
# The defaults below are the ones this survey wants: one chromosome, copy floors one step
# under ng's at the three periods that can move, and ng's own 15 bp bundle radius. Nothing
# needs to be set beyond JOBS.
#
# ## It is built for an archive that is being written to while it runs
#
# The walk takes hours, so a file that was complete when the run started can be appended to before
# the walk reaches it, and one that was mid-write at the start can be finished long before the run
# ends. **A gate at the front cannot express either.** So:
#
#   * **one walk per file**, not one per batch — a file's failure costs that file and nothing else,
#     and the region-typing pass it repeats is a fraction of the read fetching it does not;
#   * **each file is checked immediately before its walk and again immediately after**, and a file
#     whose size or mtime moved in between has its result thrown away. This is the check that
#     matters: a half-written CRAM does not reliably *error*, it can simply run out of records part
#     way through the region, and that library would come back with thinner strata and a copy floor
#     read off them with nothing saying why;
#   * **files that were not ready are retried in later passes**, so a library still being written at
#     the start is picked up once it settles rather than lost;
#   * **results are kept in a work directory as they are produced**, so re-running the same command
#     resumes instead of starting over — which is also how you add libraries to a finished survey.
#
# Knobs, all env-overridable:
#   CONTIGS=SL4.0ch01   the walk. ~90 Mb gives ~200k loci, which settles these curves; the cost
#                       scales with the number of libraries, so a wider walk buys precision you do
#                       not need. Set to "" to walk the whole genome.
#   MIN_COPIES=         the copy floors to type from — one number for every period, or a table of
#                       six (`6,2,4,3,3,3`). Empty means ng's defaults, [6,4,4,3,3,3].
#
#                       **Do not set this to a low uniform value, and an earlier version of this
#                       script defaulted it to 2.** The floor is read twice — once by `prefilter`,
#                       which runs *before* bundling, and once by `classify`
#                       (`segment_criteria.rs:601`, `:985`) — so lowering it changes what counts as
#                       a neighbouring repeat, not only what is admitted as a locus. Two copies of a
#                       mononucleotide is any `AA`, which occurs every few bases, so every real
#                       tract acquires a neighbour inside the bundle threshold and the walk emits
#                       `SsrBundle` where this survey counts only `SsrSegment`. Measured over a
#                       2 Mb slice of tomato SL4.0ch01:
#
#                         floors          STR loci   bundles   bundle bp
#                         [6,4,4,3,3,3]      6,237     1,943      74,289   (ng's defaults)
#                         uniform 4          7,434    11,720     675,372
#                         uniform 3            848     7,950   1,623,636
#                         uniform 2              0         1   1,177,849   (one bundle, 1.2 Mb)
#                         [2,4,4,3,3,3]          0        18   1,599,157
#                         [6,2,4,3,3,3]      1,618    13,913   1,017,906
#
#                       **The loss is not confined to the period that moved.** Lowering period 2
#                       alone takes period-1 tracts at 6 copies from 2,678 loci to 225 — 92% gone at
#                       a period whose own floor never changed. So two settings do not extend one
#                       curve; they measure different loci. See the finding in
#                       doc/devel/ng/spec/parameter_prepass_ssr.md §5.1.
#
#                       **The default is one step down, `5,3,3,3,3,3`**, which is what the archive
#                       run wants. One step still costs loci — about half of every shared stratum —
#                       but the criterion that places the floors for periods 2-6 survives it: the
#                       guard share at a given stratum agrees between the swept and the default walk
#                       (dinucleotides at 4 repeats: 0.345 default against 0.344 swept). It is the
#                       off-reference share, which is all mononucleotides have, that comes back low.
#   JOBS=8              files walked at once. **The knob that decides whether this finishes.** One
#                       file over 90 Mb takes ~10 minutes, so 2,475 of them sequentially is 17 days
#                       and at 32-way is half a day. Each walk is independent — that is what the
#                       per-file design buys — so set this to about the core count, backing off if
#                       the archive's disks are the bottleneck rather than the CPU.
#   REGIONS=            a BED to walk instead of whole contigs. **The other way to make a large
#                       survey finish**: these curves need thousands of loci per stratum, not the
#                       ~200k that a whole chromosome gives, so a 10 Mb slice is ~9x faster and
#                       still ample. Overrides CONTIGS.
#   PASSES=3            how many times to come back for files that were not ready.
#   WAIT=900            seconds between passes. Long enough that a file being written has plausibly
#                       finished; the survey is hours long, so this costs nothing.
#   WORK=OUT.work       where per-file results accumulate. Delete it to force a full re-run.

set -euo pipefail

# ---------------------------------------------------------------------------
# **What counts as a finished file, used in all three places that ask**
# ---------------------------------------------------------------------------
#
# A result is finished when it holds at least one stratum row — not when it merely exists. The
# distinction is the whole point of this survey's last failure: 1,300 walks wrote a well-formed
# table with no rows in it, and every later pass skipped them because a file was there.
#
# **Three callers, and each one was a hole.** The worker skips a file it considers done; the pass
# loop counts one as walked; and the merge folds one into the output. Checking only the worker
# leaves a partially-finished run producing a table of phantom libraries — a `#rg` line and zero
# loci for every file not yet reached.
#
# **awk rather than grep**, because `grep -q -v` with several `-e` patterns does not mean the same
# thing in every grep on every machine (ugrep gets it wrong), and a portability bug here reads as
# "no results" rather than as an error.
has_rows() {
    [[ -s "$1" ]] || return 1
    awk -F'\t' '
        /^#/ { next }
        $1 == "read_group" || $1 == "library_key" { next }
        { found = 1; exit }
        END { exit !found }
    ' "$1"
}

# ---------------------------------------------------------------------------
# The per-file worker, re-entered through `--one`
# ---------------------------------------------------------------------------
#
# Walking one file is its own invocation so the parent can drive many at once with `xargs -P`. The
# worker reports back through `$WORK/<key>.status` rather than stdout, because several of them are
# writing at the same time and interleaved lines would be unreadable and unparseable.
if [[ "${1:-}" == "--one" ]]; then
    cram=$2
    key=$(basename "$cram" | tr -c 'A-Za-z0-9._-' '_')
    result="$SURVEY_WORK/$key.tsv"
    status="$SURVEY_WORK/$key.status"
    stamp_of() { stat -c '%s %Y' "$1" 2>/dev/null || stat -f '%z %m' "$1" 2>/dev/null || echo gone; }
    index_of() {
        local candidate
        for candidate in "$cram.crai" "$cram.bai" "${cram%.*}.crai" "${cram%.*}.bai"; do
            [[ -f "$candidate" ]] && { echo "$candidate"; return 0; }
        done
        return 1
    }
    # A result that exists but measures nothing is not done — drop it and walk again. This is what
    # replaces the previous run's empty results in place, so a re-run cleans up after itself.
    has_rows "$result" && exit 0
    rm -f "$result"
    if [[ ! -f "$cram" ]]; then echo "missing" > "$status"; exit 0; fi
    if ! index=$(index_of); then echo "no index yet" > "$status"; exit 0; fi
    if [[ "$cram" -nt "$index" ]]; then
        echo "index older than the file — still growing, or not re-indexed" > "$status"
        exit 0
    fi
    if [[ "$SURVEY_QUICKCHECK" == "1" ]] && ! samtools quickcheck "$cram" 2>/dev/null; then
        echo "no EOF block — still being written" > "$status"
        exit 0
    fi
    before=$(stamp_of "$cram")
    read -r -a survey_args <<< "$SURVEY_ARGS"
    if ! "$SURVEY_BIN" "${survey_args[@]}" "$SURVEY_REF" "$cram" \
            > "$result.partial" 2> "$SURVEY_WORK/$key.err"; then
        rm -f "$result.partial"
        echo "walk failed: $(tail -n 1 "$SURVEY_WORK/$key.err" 2>/dev/null | tr '\t' ' ')" > "$status"
        exit 0
    fi
    # **The check that catches the silent case.** A file appended to during its own walk yields a
    # result that looks like a thin library rather than an error, so the stamp — not the exit code —
    # decides whether the result is trustworthy.
    if [[ "$before" != "$(stamp_of "$cram")" ]]; then
        rm -f "$result.partial"
        echo "changed during its own walk — result discarded" > "$status"
        exit 0
    fi
    # **And the same check against a walk that exited 0 having measured nothing.** The binary refuses
    # an empty table itself, so this only fires for a binary that predates that refusal — which is
    # exactly the case that produced 1,300 empty results, so it is worth keeping rather than trusting
    # everyone to have rebuilt.
    if ! has_rows "$result.partial"; then
        rm -f "$result.partial"
        echo "walked successfully and measured nothing — no stratum rows (stale binary?)" > "$status"
        exit 0
    fi
    mv "$result.partial" "$result"
    rm -f "$status"
    exit 0
fi

OUT=${1:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
REF=${2:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
shift 2
(($# > 0)) || { echo "no CRAMs given" >&2; exit 1; }

# `${CONTIGS-...}` without the colon, so an explicitly empty CONTIGS means "the whole genome"
# rather than falling back to the default — which is what the usage above promises.
CONTIGS=${CONTIGS-SL4.0ch01}
REGIONS=${REGIONS:-}
# One step below ng's `[6,4,4,3,3,3]` at the three periods that can move — see the knob's note
# above for why it is one step and not a low uniform value. Set to "" for ng's defaults exactly.
MIN_COPIES=${MIN_COPIES-5,3,3,3,3,3}
PASSES=${PASSES:-3}
WAIT=${WAIT:-900}
WORK=${WORK:-$OUT.work}
JOBS=${JOBS:-8}

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE="$REPO_ROOT/examples/ng_str_stutter_by_library.rs"

# **Two target directories, and looking in only one of them is how a survey runs a stale binary.**
# `scripts/dev.sh` builds with `CARGO_TARGET_DIR=$PROJECT_DIR/target-container`, so a rebuild done
# the way CLAUDE.md prescribes lands *beside* the host tree rather than in it. A script hard-coded to
# `target/` then keeps running whatever was built on the host months ago — which is what happened on
# rick, where the binary's `#rg` line still had 12 fields. Take the newer of the two and say which.
BIN=""
for candidate in "$REPO_ROOT/target-container/release/examples/ng_str_stutter_by_library" \
                 "$REPO_ROOT/target/release/examples/ng_str_stutter_by_library"; do
    [[ -x "$candidate" ]] || continue
    if [[ -z "$BIN" || "$candidate" -nt "$BIN" ]]; then BIN=$candidate; fi
done
[[ -n "$BIN" ]] || {
    echo "not built: no ng_str_stutter_by_library in target/ or target-container/" >&2
    echo "  ./scripts/dev.sh cargo build --release --example ng_str_stutter_by_library" >&2
    exit 1
}

# **A binary older than its source is the failure this survey already paid for**, and it is silent:
# an out-of-date walk writes a well-formed table with the wrong columns in it. Refuse rather than
# warn — a warning on stderr at the head of a run measured in hours is a warning nobody sees.
if [[ -f "$SOURCE" && "$SOURCE" -nt "$BIN" ]]; then
    echo "$BIN is older than $SOURCE — rebuild before surveying:" >&2
    echo "  ./scripts/dev.sh cargo build --release --example ng_str_stutter_by_library" >&2
    exit 1
fi
echo "using $BIN" >&2
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }

mkdir -p "$WORK"

# ---------------------------------------------------------------------------
# The work directory belongs to one set of walk parameters
# ---------------------------------------------------------------------------
#
# **Resuming with a different walk would silently mix incomparable results.** The per-file results
# are keyed by filename alone, so a run over `SL4.0ch01` and a later one over a 10 Mb BED would land
# side by side in the same directory and merge into one table whose libraries were measured over
# different loci — a survey that looks whole and compares nothing. Since narrowing the walk is the
# obvious thing to reach for when a run is taking too long, this is a mistake worth making
# impossible rather than documenting.
PARAMS="contigs=${CONTIGS:-} regions=${REGIONS:-} min_copies=$MIN_COPIES"
STAMP="$WORK/.survey-params"
existing=$(find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 2>/dev/null | wc -l | tr -d ' ')
if [[ -f "$STAMP" ]]; then
    recorded=$(cat "$STAMP")
    if [[ "$recorded" != "$PARAMS" ]]; then
        echo "this work directory was built with different walk parameters:" >&2
        echo "  it holds: $recorded" >&2
        echo "  you asked: $PARAMS" >&2
        echo "Results from two different walks are not comparable. Use a different OUT path, or" >&2
        echo "delete $WORK to start over." >&2
        exit 1
    fi
elif ((existing > 0)); then
    # Written by a version of this script that did not stamp, so the parameters are unknown.
    echo "$WORK holds $existing result(s) but no record of how they were walked." >&2
    echo "If they came from the same walk you are asking for now, re-run with" >&2
    echo "  SURVEY_ADOPT_EXISTING=1" >&2
    echo "Otherwise delete $WORK. Mixing two walks makes a table that compares nothing." >&2
    [[ "${SURVEY_ADOPT_EXISTING:-0}" == "1" ]] || exit 1
    echo "$PARAMS" > "$STAMP"
else
    echo "$PARAMS" > "$STAMP"
fi

walk_args=()
if [[ -n "$REGIONS" ]]; then
    [[ -f "$REGIONS" ]] || { echo "regions BED not found: $REGIONS" >&2; exit 1; }
    walk_args=(--regions "$REGIONS")
elif [[ -n "$CONTIGS" ]]; then
    walk_args=(--contigs "$CONTIGS")
fi
if [[ -n "$MIN_COPIES" ]]; then
    walk_args+=(--min-copies "$MIN_COPIES")
fi

command -v samtools >/dev/null && HAVE_QUICKCHECK=1 || HAVE_QUICKCHECK=0
((HAVE_QUICKCHECK)) || echo "warning: samtools not on PATH — no EOF check, only indexes and stamps" >&2

key_of() { basename "$1" | tr -c 'A-Za-z0-9._-' '_'; }

export SURVEY_WORK="$WORK" SURVEY_REF="$REF" SURVEY_BIN="$BIN"
export SURVEY_QUICKCHECK="$HAVE_QUICKCHECK"
export SURVEY_ARGS="${walk_args[*]}"
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")

declare -a pending=("$@")
total=${#pending[@]}
echo "surveying $total file(s) over ${REGIONS:-${CONTIGS:-the whole genome}}, copy floors ${MIN_COPIES:-the ng defaults}, $JOBS at a time" >&2
echo "  work directory: $WORK (delete it to force a full re-run)" >&2

declare -a deferred_files=() deferred_reasons=()

for ((pass = 1; pass <= PASSES; pass++)); do
    ((${#pending[@]} > 0)) || break
    ((pass == 1)) || {
        echo "pass $pass: waiting ${WAIT}s for ${#pending[@]} unfinished file(s) to settle" >&2
        sleep "$WAIT"
    }
    # Clear last pass's verdicts so a file that has since become readable is not still carrying one.
    for cram in "${pending[@]}"; do rm -f "$WORK/$(key_of "$cram").status"; done

    echo "pass $pass: ${#pending[@]} file(s) to walk" >&2
    # **`xargs -P` rather than a hand-rolled job pool**, because it is the portable way to keep N
    # slots full: a `wait`-on-every-N loop stalls the whole group on its slowest member, and these
    # files differ several-fold in size.
    #
    # **The file is appended as the last argument rather than substituted with `-I`.** `-I` implies
    # one argument per run, so passing `-n 1` beside it made GNU xargs warn that the two are
    # mutually exclusive — harmless, since `-I` won and one file per worker is what was wanted, but
    # `-I` also has a history of interacting badly with `-P` across implementations. Appending needs
    # neither flag and the worker's signature (`--one <cram>`) already puts the file last.
    printf '%s\0' "${pending[@]}" \
        | xargs -0 -n 1 -P "$JOBS" "$SELF" --one \
        || echo "  (some workers reported failures; their reasons are in the work directory)" >&2

    declare -a still_pending=()
    deferred_files=()
    deferred_reasons=()
    done_count=0
    for cram in "${pending[@]}"; do
        key=$(key_of "$cram")
        if has_rows "$WORK/$key.tsv"; then
            done_count=$((done_count + 1))
            continue
        fi
        still_pending+=("$cram")
        deferred_files+=("$cram")
        deferred_reasons+=("$(cat "$WORK/$key.status" 2>/dev/null || echo 'no result and no reason recorded')")
    done
    echo "  pass $pass done: $done_count walked, ${#still_pending[@]} outstanding" >&2
    pending=("${still_pending[@]:-}")
    # `${arr[@]:-}` yields one empty element on an empty array under `set -u`; drop it.
    ((${#pending[@]} == 1)) && [[ -z "${pending[0]}" ]] && pending=()
done

# **Only files holding stratum rows are results.** A run stopped part way leaves the rest of the
# work directory as it was — including any empty result an earlier run wrote — and folding those in
# would put a `#rg` line and zero loci into the table for every library never actually walked.
# `if` and not `has_rows "$f" && echo "$f"`: under `set -e` an AND-list whose last command fails is
# fatal, so the first rejected file would kill the loop — and silently, since it runs in a pipeline.
list_results() {
    find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 | sort | while read -r f; do
        if has_rows "$f"; then echo "$f"; fi
    done
}
surveyed=$(list_results | wc -l | tr -d ' ')
((surveyed > 0)) || { echo "no file could be surveyed; nothing written" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Every result must have been walked under the same settings
# ---------------------------------------------------------------------------
#
# **The parameter stamp cannot cover this and it is the case that bites.** The stamp records what
# the *driver* was asked for, so it catches a resume with a different `--contigs` or `--min-copies`.
# It cannot see the copy floors and the bundle radius when they come from the binary's own
# defaults — so **rebuilding part way through a survey silently changes the walk**, and the results
# from before and after merge into one table that compares nothing. That is not hypothetical: ng's
# bundle radius moved from 20 bp to 15 on 2026-08-07 while a survey of the tomato archive was
# running, and only clearing the work directory kept the two apart.
#
# So each result carries the settings it was produced under (`#config`), and they must all agree.
configs=$(list_results | while read -r f; do awk -F'\t' '$1 == "#config" { print; exit }' "$f"; done | sort -u)
if (($(echo "$configs" | grep -c . ) > 1)); then
    echo "the results in $WORK were not all walked under the same settings:" >&2
    echo "$configs" | sed 's/^/  /' >&2
    echo "This happens when the binary is rebuilt part way through a survey — the copy floors and" >&2
    echo "the bundle radius come from its defaults, so a rebuild changes the walk with nothing to" >&2
    echo "say so. Delete $WORK and start over, or move the older results aside." >&2
    exit 1
fi

# Merge. **The numeric read_group is minted per invocation**, so rows are re-keyed onto
# `(file, rg_id)` — the stable identity, since the SAM specification makes `@RG ID` unique within
# its file — before the per-file results are joined.
echo "merging $surveyed result(s)" >&2
{
    # **The output says what it does not contain.** A survey missing a third of the archive because
    # those CRAMs were mid-write looks exactly like a survey of an archive that size, so the
    # exclusions travel with the numbers rather than in a terminal scrollback.
    # `walked=` is the BED when one was given, since `REGIONS` overrides `CONTIGS` — an earlier
    # version printed `contigs=` either way, so a BED-restricted survey claimed a whole chromosome.
    printf '#survey\twalked=%s\tmin_copies=%s\tgiven=%s\tsurveyed=%s\tmissing=%s\tpasses=%s\n' \
        "${REGIONS:-${CONTIGS:-ALL}}" "${MIN_COPIES:-default}" "$total" "$surveyed" \
        "${#deferred_files[@]}" "$PASSES"
    for ((i = 0; i < ${#deferred_files[@]}; i++)); do
        printf '#missing\t%s\t%s\n' "${deferred_files[$i]}" "${deferred_reasons[$i]}"
    done
    # Three passes so the file reads top to bottom — every `#rg`, then every `#floor`, then the
    # rows — rather than interleaving one library's three blocks with the next library's.
    # The results are one small file per library, so re-reading them twice more costs nothing.
    results=$(list_results)

    printf '#rg_columns\tlibrary_key\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\texperiment_origin\tplatform\tloci\treads\tmean_tract_bases_per_read\tfile\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); $2 = path[n] "::" $3; print }
        ' "$f"
    done

    # **What region typing gave each file, beside what its reads said.** A library whose strata are
    # thin because its tracts bundled is a different finding from one that is merely shallow, and
    # only this block separates them. Keyed by file, since typing is a property of the reference and
    # the walk rather than of a read group.
    printf '#typing_columns\tfile\tspans\tssr_loci\tssr_bundles\tssr_bundle_bp\tsatellites\trepeat_bp_with_no_locus\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' -v src="$f" '
            $1 == "#rg" { n = split($13, path, "/"); file = path[n]; next }
            $1 == "#typing" { $1 = "#typing"; print $1, file, $2, $3, $4, $5, $6, $7 }
        ' "$f"
    done

    printf '#floor_columns\tlibrary_key\tperiod\timplied_floor\tcriterion\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); stable[$2] = path[n] "::" $3; next }
            $1 == "#floor" { $2 = stable[$2]; print }
        ' "$f"
    done

    printf 'library_key\tperiod\trepeats\tloci\treads\toff_ref_reads\toff_ref_share\tnot_whole_reads\tguard_share\tend_bucket_reads\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); stable[$2] = path[n] "::" $3; next }
            /^#/ { next }
            $1 == "library_key" || $1 == "read_group" { next }
            { $1 = stable[$1]; print }
        ' "$f"
    done
} > "$OUT"

libraries=$(grep -c '^#rg	' "$OUT" || true)
# Exclude the column-header line, which is not a comment and is not a stratum either.
rows=$(grep -v '^#' "$OUT" | grep -vc '^library_key	' || true)
echo "wrote $OUT — $libraries library/libraries, $rows stratum rows" >&2
if ((${#deferred_files[@]} > 0)); then
    echo "${#deferred_files[@]} file(s) never became readable and are listed as #missing in the output:" >&2
    for ((i = 0; i < ${#deferred_files[@]}; i++)); do
        echo "  $(basename "${deferred_files[$i]}") — ${deferred_reasons[$i]}" >&2
    done
    echo "Re-run the same command later to pick them up; finished files are not re-walked." >&2
fi
echo >&2
echo "The blocks in it:" >&2
echo "  #survey/#missing  what was walked and what was not, with the reason" >&2
echo "  #rg               one line per read group: which library, and how much data it gave" >&2
echo "  #floor            the copy floor each period's own data implies, per library — the answer" >&2
echo "  rows              the per-(library, period, repeat count) curves behind the floors" >&2
echo >&2
echo "Before comparing libraries, join to rick_sample_manifest.sh on the \`@RG\` id: read length" >&2
echo "is a confound this tool cannot see, and mixing 100 bp with 150 bp libraries makes stutter" >&2
echo "appear to start at different tract lengths for purely geometric reasons." >&2
