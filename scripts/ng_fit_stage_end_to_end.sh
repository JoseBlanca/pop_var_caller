#!/usr/bin/env bash
#
# **Psp mode's three commands and its repair, run end to end on real reads** — plan step E1 of
# doc/devel/ng/impl_plan/psp_census_pair.md.
#
#   generate-psps        alignments        ->  <sample>.psp, census sealed inside it
#   estimate-parameters  psps              ->  cohort.parameters.toml
#   call-from-psps       psps + that file  ->  the VCF
#   regenerate-census    a copied psp      ->  its trailer rebuilt, and the file whole again
#
# and then the question the calling stage exists to answer: **what do the fitted numbers change?**
# The same cohort is called twice, once with --defaults and once with the file the fit wrote,
# and the two VCFs are compared record for record and genotype for genotype.
#
#   scripts/ng_fit_stage_end_to_end.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <out-dir> <alignment.cram>...
#
# Run it through the dev container — the release binaries it looks for are built there, and on
# macOS they will not run on the host:
#
#   ./scripts/dev.sh scripts/ng_fit_stage_end_to_end.sh ...
#
# **There is no census file anywhere in this run**, and no step that writes one: since
# psp_census_pair.md §3 a sample is one file, and the census the fit reads is in that file's
# trailer. What replaced the old walk-versus-rebuild `cmp` over two census files is stronger and
# is the last step here — a walked psp copied, the copy's census dropped, `regenerate-census` run
# on it, and the copy identical to the walked file **byte for byte, whole**: header, blocks,
# index, trailer, footer. The same oracle on fixtures is
# `a_regenerated_psp_is_the_walked_one_byte_for_byte`.
set -uo pipefail

if (( $# < 5 )); then
    # Every leading comment line, so the range cannot drift from the block it prints.
    awk 'NR == 1 { next } /^#/ { print; next } { exit }' "$0"
    exit 2
fi
reference=$1; catalog=$2; regions=$3; out=$4
shift 4

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# The newer of the two release trees: a container build points CARGO_TARGET_DIR at
# target-container, a host build does not, and a machine with no container runtime has only the
# second (CLAUDE.md).
newest_release_build_of() {
    local what="$1" found="" candidate
    for candidate in "$root/target-container/release/$what" "$root/target/release/$what"; do
        if [[ -x "$candidate" ]] && { [[ -z "$found" ]] || [[ "$candidate" -nt "$found" ]]; }; then
            found=$candidate
        fi
    done
    printf '%s' "$found"
}

bin=$(newest_release_build_of pop_var_caller)
if [[ -z "$bin" ]]; then
    echo "no release build of pop_var_caller; build one with" >&2
    echo "  ./scripts/dev.sh cargo build --release --bin pop_var_caller" >&2
    exit 1
fi
# **The last step's oracle needs a psp owed a rebuild**, and nothing that ships empties a psp's
# trailer — see the block at the head of this file.
drop_census=$(newest_release_build_of examples/ng_psp_drop_census)
if [[ -z "$drop_census" ]]; then
    echo "no release build of the ng_psp_drop_census example, which the last step needs; build it with" >&2
    echo "  ./scripts/dev.sh cargo build --release --example ng_psp_drop_census" >&2
    exit 1
fi

rm -rf "$out"
mkdir -p "$out/psps"

alignments=()
for cram in "$@"; do
    alignments+=(--alignment "$cram")
done

say() { printf '\n=== %s ===\n' "$1"; }

say "1. generate-psps"
"$bin" generate-psps \
    --reference "$reference" --catalog "$catalog" "${alignments[@]}" \
    --regions "$regions" --output-dir "$out/psps" > "$out/generate-psps.log" 2>&1 || {
    echo "generate-psps failed:" >&2; cat "$out/generate-psps.log" >&2; exit 1; }
tail -1 "$out/generate-psps.log"

say "the walk left one file a sample, with a census in it"
# **The walk's own report line**, which says both halves of that heading — "N samples: X bytes of
# psp, of which Y bytes are census". The psps on disk are counted separately and the two counts
# have to agree, so a walk that dropped a sample is caught here rather than becoming a smaller
# cohort three steps below. A census total of zero would mean every psp was sealed with an empty
# trailer.
#
# **Each parse insists on exactly one match**, and that is not caution for its own sake: the
# per-sample lines say "of which N bytes are **its** census" and are printed twice a sample, so a
# pattern one word looser would return twelve numbers for six samples — which is the shape this
# check was found in at Checkpoint A. A count other than one means the wording drifted and this
# check has stopped checking.
psps=0
for psp in "$out/psps"/*.psp; do
    [[ -e "$psp" ]] || break
    psps=$(( psps + 1 ))
done
if (( psps == 0 )); then
    echo "the walk left no psp; nothing below is worth reading" >&2
    exit 1
fi

one_number_from() {
    local what="$1" pattern="$2" found
    found=$(sed -n "$pattern" "$out/generate-psps.log")
    if [[ $(printf '%s\n' "$found" | grep -c '^[0-9][0-9]*$') != 1 ]]; then
        echo "generate-psps' report line has changed shape; the $what check no longer reads it" >&2
        return 1
    fi
    printf '%s' "$found"
}

walked=$(one_number_from "sample count" 's/^\([0-9][0-9]*\) samples: .* bytes of psp,.*/\1/p') || exit 1
inside=$(one_number_from "census total" 's/^[0-9][0-9]* samples: .* of which \([0-9][0-9]*\) bytes are census.*/\1/p') || exit 1
if (( walked != psps )); then
    echo "the walk says it wrote $walked samples and there are $psps psps on disk" >&2
    exit 1
fi
if (( inside == 0 )); then
    echo "the walk sealed every psp with an empty trailer" >&2
    exit 1
fi
echo "  $psps psps, one a sample, holding $inside bytes of census between them"
grep -m 1 "bytes of psp" "$out/generate-psps.log" || true

say "2. estimate-parameters"
# **From the psps themselves** (plan step C2): the census it reads is the one in each file's
# trailer, and the ground and the repeat criteria come from the headers, so none of the walk's
# settings is retyped here — there is no --census, and no flag saying what counts as a repeat.
"$bin" estimate-parameters \
    --reference "$reference" --catalog "$catalog" \
    --psp "$out/psps" --output "$out/cohort.parameters.toml" \
    > "$out/estimate-parameters.log" 2>&1 || {
    echo "estimate-parameters failed:" >&2; cat "$out/estimate-parameters.log" >&2; exit 1; }
tail -1 "$out/estimate-parameters.log"
# **The file, not only the exit status.** Without this the first thing to notice an empty
# parameters file is the second call below, which would report itself as the failure.
if [[ ! -s "$out/cohort.parameters.toml" ]]; then
    echo "the fit returned 0 and wrote no parameters file worth reading" >&2
    exit 1
fi

say "3. call-from-psps, twice"
# **Each pass is made to say which numbers it scored with**, off its own report line — "numbers
# behind the calls: N of 7 groups the file says were fitted". The defaults pass has to say 0 and
# the fitted pass more than 0. Without it, a fitted pass that quietly fell back to the defaults
# would write the same VCF twice and the comparison below would print zeros and exit 0 — which
# reads exactly like "the fitted numbers changed nothing".
#
# A cohort too small to fit any group would stop here, and that is right for this script: it
# exists to answer what fitting changes, and a run with nothing fitted cannot answer it.
for how in defaults fitted; do
    if [[ $how == defaults ]]; then
        numbers=(--defaults)
        expected="exactly 0"
    else
        numbers=(--parameters "$out/cohort.parameters.toml")
        expected="more than 0"
    fi
    "$bin" call-from-psps \
        --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        "${numbers[@]}" --threads 4 --output "$out/$how.vcf" \
        > "$out/call-$how.log" 2>&1 || {
        echo "call-from-psps --$how failed:" >&2; cat "$out/call-$how.log" >&2; exit 1; }
    fitted_groups=$(sed -n 's/^numbers behind the calls: \([0-9][0-9]*\) of [0-9][0-9]* groups.*/\1/p' \
        "$out/call-$how.log")
    if [[ $(printf '%s\n' "$fitted_groups" | grep -c '^[0-9][0-9]*$') != 1 ]]; then
        echo "call-from-psps' report line has changed shape; the $how pass cannot be told from the other" >&2
        exit 1
    fi
    if { [[ $how == defaults ]] && (( fitted_groups != 0 )); } ||
       { [[ $how == fitted ]] && (( fitted_groups == 0 )); }; then
        echo "the $how pass scored with $fitted_groups fitted groups of numbers, and it should be $expected" >&2
        exit 1
    fi
    echo "  $how: $(grep -cv '^#' "$out/$how.vcf") records, $fitted_groups groups of numbers fitted"
done

say "what the fitted numbers change"
# **The one thing asserted here is that both calls produced records.** Not that the fitted numbers
# changed something: at one sample, or on ground where the fit lands close to the defaults, calling
# twice and getting the same VCF is a legitimate answer, and a script that demanded a difference
# would fail on a correct run. What is never legitimate is two empty VCFs, which without this check
# print four lines of zeros and exit 0 — and so does a fitted call that quietly used the defaults.
awk '
    function gts(line,   n, i, f, out) {
        n = split(line, f, "\t")
        out = ""
        for (i = 10; i <= n; i++) {
            split(f[i], g, ":")
            out = out (i > 10 ? "\t" : "") g[1]
        }
        return out
    }
    FNR == NR {
        if ($0 ~ /^#/) next
        split($0, f, "\t"); key = f[1] ":" f[2] ":" f[4] ":" f[5]
        left[key] = gts($0); leftn++
        next
    }
    {
        if ($0 ~ /^#/) next
        split($0, f, "\t"); key = f[1] ":" f[2] ":" f[4] ":" f[5]
        rightn++
        if (!(key in left)) { only_right++; next }
        seen[key] = 1
        r = gts($0)
        if (r != left[key]) {
            records_differing++
            nl = split(left[key], a, "\t"); split(r, b, "\t")
            for (i = 1; i <= nl; i++) if (a[i] != b[i]) genotypes_differing++
        }
        total_genotypes += split(r, b, "\t")
    }
    END {
        for (k in left) if (!(k in seen)) only_left++
        printf "  records: %d with the defaults, %d with the fitted numbers\n", leftn, rightn
        printf "  records only one of them called: %d defaults-only, %d fitted-only\n",
               only_left + 0, only_right + 0
        printf "  of the %d records both called, %d differ in at least one genotype\n",
               leftn - (only_left + 0), records_differing + 0
        printf "  genotypes: %d differ out of %d compared\n",
               genotypes_differing + 0, total_genotypes + 0
        if (leftn + 0 == 0 || rightn + 0 == 0) {
            printf "at least one of the two calling passes emitted no record at all\n" > "/dev/stderr"
            exit 1
        }
    }
' "$out/defaults.vcf" "$out/fitted.vcf" || exit 1

say "what the parameters file says it fitted"
grep -m 1 -A 4 "groups of numbers in this file" "$out/cohort.parameters.toml" || true

say "4. regenerate-census, and the psp it rebuilds is the walked one byte for byte"
# **The whole-file oracle on real reads** (plan step D4). Each walked psp is copied, the copy's
# census is dropped, `regenerate-census` rebuilds it from the copy's own records, and the copy has
# to be the walked file again — header, blocks, index, trailer and footer, every byte.
#
# **Dropping the census first is what makes the comparison mean anything.** A psp whose census is
# the one this run would write is skipped and its records are never read (psp_census_pair.md §8),
# so a copy handed straight to the command comes back untouched and `cmp` passes without a rebuild
# having happened. The two comparisons below are therefore a pair: **before** the rebuild each
# copy must differ from its original, and **after** it each must match.
#
# The originals are left alone — the command is pointed at the copies' directory — so nothing
# above this line is disturbed by it.
copies="$out/copies"
mkdir -p "$copies"
copied=0
for psp in "$out/psps"/*.psp; do
    cp "$psp" "$copies/$(basename "$psp")" || {
        echo "copying $(basename "$psp") failed" >&2; exit 1; }
    copied=$(( copied + 1 ))
done
if (( copied != psps )); then
    echo "$copied of the $psps psps were copied; the comparison below would be against a file that is not there" >&2
    exit 1
fi

# **Both comparisons address the same pair, through one function**, so that they cannot come to
# compare different things — a second loop naming its own two paths is one edit away from
# comparing a copy with itself and passing.
walked_and_its_copy_differ() {
    ! cmp -s "$out/psps/$1" "$copies/$1"
}

"$drop_census" "$copies"/*.psp > "$out/drop-census.log" 2>&1 || {
    echo "dropping the copies' censuses failed:" >&2; cat "$out/drop-census.log" >&2; exit 1; }
armed=0
for psp in "$out/psps"/*.psp; do
    if ! walked_and_its_copy_differ "$(basename "$psp")"; then
        echo "  $(basename "$psp"): the copy is already identical, so the rebuild below proves nothing" >&2
        armed=$(( armed + 1 ))
    fi
done
if (( armed != 0 )); then
    echo "$armed of the $psps copies were not armed; nothing below is worth reading" >&2
    exit 1
fi

"$bin" regenerate-census \
    --reference "$reference" --catalog "$catalog" --psp "$copies" \
    > "$out/regenerate-census.log" 2>&1 || {
    echo "regenerate-census failed:" >&2; cat "$out/regenerate-census.log" >&2; exit 1; }
# **Named, not positional.** This log holds both streams: every sample's line goes to stderr as
# the sample is decided, and the report then goes to stdout with its summary ahead of its own
# per-sample lines. So the summary sits in the middle, and `head -1` and `tail -1` would each
# show a sample instead — which is what the two commands above use, because their summary is
# their last line.
grep -m 1 "^regenerated " "$out/regenerate-census.log" || {
    echo "regenerate-census returned 0 and printed no summary line" >&2; exit 1; }

differing=0
for psp in "$out/psps"/*.psp; do
    if walked_and_its_copy_differ "$(basename "$psp")"; then
        echo "  $(basename "$psp"): DIFFERS from the walked file after being rebuilt"
        differing=$(( differing + 1 ))
    else
        echo "  $(basename "$psp"): identical to the walked file, whole"
    fi
done
if (( differing != 0 )); then
    echo "  $differing of the $psps rebuilt psps are not the walked file"
    exit 1
fi
