#!/usr/bin/env bash
#
# **The four commands of psp mode, run end to end on real reads** — plan steps D1 and D2 of
# doc/devel/ng/impl_plan/parameter_prepass_runs.md.
#
#   generate-psps        alignments  ->  <sample>.psp and <sample>.census
#   generate-census      psps        ->  <sample>.census, again, from the stored files
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
# It also checks, for free, that the two routes to a census still agree on real reads: the
# censuses generate-psps wrote during the walk against the ones generate-census built from the
# stored psps afterwards.
set -uo pipefail

if (( $# < 5 )); then
    sed -n '2,24p' "$0"
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
mkdir -p "$out/psps" "$out/rebuilt"

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
"$bin" generate-census \
    --reference "$reference" --catalog "$catalog" \
    --psp "$out/psps" --output-dir "$out/rebuilt" > "$out/generate-census.log" 2>&1 || {
    echo "generate-census failed:" >&2; cat "$out/generate-census.log" >&2; exit 1; }
tail -1 "$out/generate-census.log"

say "do the two routes still agree on real reads?"
same=1
for walked in "$out/psps"/*.census; do
    rebuilt="$out/rebuilt/$(basename "$walked")"
    if cmp -s "$walked" "$rebuilt"; then
        echo "  $(basename "$walked"): identical"
    else
        echo "  $(basename "$walked"): DIFFERENT"
        same=0
    fi
done
(( same == 1 )) || { echo "the two producers disagree; nothing below is worth reading" >&2; exit 1; }

say "3. estimate-parameters"
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
