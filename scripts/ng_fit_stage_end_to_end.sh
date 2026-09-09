#!/usr/bin/env bash
#
# **The four commands of psp mode, run end to end on real reads** — plan steps D1 and D2 of
# doc/devel/ng/impl_plan/parameter_prepass_runs.md.
#
#   generate-psps        alignments  ->  <sample>.psp, with its census inside
#   generate-census      psps        ->  <sample>.census, from the stored files
#   estimate-parameters  censuses    ->  cohort.parameters.toml
#   call-from-psps       psps + that file  ->  the VCF
#
# and then the question the last stage exists to answer: **what do the fitted numbers change?**
# The same cohort is called twice, once with --defaults and once with the file the fit wrote,
# and the two VCFs are compared record for record and genotype for genotype.
#
#   scripts/ng_fit_stage_end_to_end.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <out-dir> <alignment.cram>...
#
# Run it through the dev container — the release binary it looks for is built there, and on
# macOS it will not run on the host:
#
#   ./scripts/dev.sh scripts/ng_fit_stage_end_to_end.sh ...
#
# **It no longer checks the two routes against each other.** Since psp_census_pair.md §3 the walk
# seals its census into the psp's trailer instead of writing a file, and no shipped subcommand
# writes a psp's trailer back out — a script could carve it from the footer's offset and length,
# but a script has no business seeking into a psp. The comparison has two homes instead: on
# fixtures, the test `the_two_producers_agree_on_a_cohort_with_a_repeat_tract`; on real reads,
# scripts/ng_census_route_cost.sh, whose harness writes both routes' census bytes out for exactly
# that. **A stronger check reaches this script at plan step E1** — D4's oracle, a psp copied,
# regenerated with regenerate-census, and identical to the original byte for byte, whole.
set -uo pipefail

if (( $# < 5 )); then
    # Every leading comment line, so the range cannot drift from the block it prints.
    awk 'NR == 1 { next } /^#/ { print; next } { exit }' "$0"
    exit 2
fi
reference=$1; catalog=$2; regions=$3; out=$4
shift 4

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bin=""
for candidate in "$root/target-container/release/pop_var_caller_exp" \
                 "$root/target/release/pop_var_caller_exp"; do
    if [[ -x "$candidate" ]] && { [[ -z "$bin" ]] || [[ "$candidate" -nt "$bin" ]]; }; then
        bin=$candidate
    fi
done
if [[ -z "$bin" ]]; then
    echo "no release build of pop_var_caller_exp; build one first" >&2
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

say "2. generate-census, from the stored psps"
# **Into the psps' own directory**, because estimate-parameters opens each census against the psp
# it names and refuses one that has none — measured: pointing it at a directory of its own failed
# with "SRS3394606's census names a psp and there is none at .../rebuilt/SRS3394606.psp". The walk
# itself no longer leaves a file here, so nothing is overwritten.
"$bin" generate-census \
    --reference "$reference" --catalog "$catalog" \
    --psp "$out/psps" --output-dir "$out/psps" > "$out/generate-census.log" 2>&1 || {
    echo "generate-census failed:" >&2; cat "$out/generate-census.log" >&2; exit 1; }
tail -1 "$out/generate-census.log"

say "every psp got a census rebuilt beside it"
# **Per psp and by name, not two totals.** The next step opens each census against the psp beside
# it and refuses one that has none, so what matters is the pairing rather than the count — and a
# count would pass on six psps and six censuses whose stems did not match, which is the shape the
# first attempt at this repair actually produced.
psps=0
missing=0
for psp in "$out/psps"/*.psp; do
    psps=$(( psps + 1 ))
    if [[ ! -e "${psp%.psp}.census" ]]; then
        echo "  no census beside $(basename "$psp")" >&2
        missing=1
    fi
done
if (( psps == 0 )); then
    echo "the walk left no psp; nothing below is worth reading" >&2
    exit 1
fi
if (( missing != 0 )); then
    echo "nothing below is worth reading" >&2
    exit 1
fi
# **And the walk built a census, not an empty one.** The run's own total, off its report line —
# "N samples: X bytes of psp, of which Y bytes are census". It is the walk's own count rather than
# a reopen of the files, so it says a census that big was built; whether those bytes reached the
# psps is what E1's whole-file oracle will answer. A total of zero, or a line this cannot parse,
# stops the run: the second means the wording drifted and this check has stopped checking.
inside=$(sed -n 's/.*of which \([0-9][0-9]*\) bytes are census.*/\1/p' "$out/generate-psps.log")
if [[ -z "$inside" ]]; then
    echo "generate-psps' report line has changed shape; this check no longer reads it" >&2
    exit 1
fi
if (( inside == 0 )); then
    echo "the walk sealed every psp with an empty trailer" >&2
    exit 1
fi
echo "  $psps psps, each with a census beside it, and $inside bytes of census inside them"
grep -m 1 "bytes of psp" "$out/generate-psps.log" || true

say "3. estimate-parameters"
# **From the censuses step 2 rebuilt**, which are beside the psps they were built from — the
# walk no longer leaves any there. Plan step C2 is what moves this command onto the psps
# themselves, and step 2 goes with it.
"$bin" estimate-parameters \
    --reference "$reference" --catalog "$catalog" \
    --census "$out/psps" --output "$out/cohort.parameters.toml" \
    > "$out/estimate-parameters.log" 2>&1 || {
    echo "estimate-parameters failed:" >&2; cat "$out/estimate-parameters.log" >&2; exit 1; }
tail -1 "$out/estimate-parameters.log"

say "4. call-from-psps, twice"
for how in defaults fitted; do
    if [[ $how == defaults ]]; then
        numbers=(--defaults)
    else
        numbers=(--parameters "$out/cohort.parameters.toml")
    fi
    "$bin" call-from-psps \
        --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        "${numbers[@]}" --threads 4 --output "$out/$how.vcf" \
        > "$out/call-$how.log" 2>&1 || {
        echo "call-from-psps --$how failed:" >&2; cat "$out/call-$how.log" >&2; exit 1; }
    echo "  $how: $(grep -cv '^#' "$out/$how.vcf") records"
done

say "what the fitted numbers change"
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
    }
' "$out/defaults.vcf" "$out/fitted.vcf"

say "what the parameters file says it fitted"
grep -m 1 -A 4 "groups of numbers in this file" "$out/cohort.parameters.toml" || true
