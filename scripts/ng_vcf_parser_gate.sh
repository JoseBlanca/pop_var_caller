#!/usr/bin/env bash
#
# Milestone C's gate: does an outside parser accept ng's VCF?
#
# The conventions ng inherits from production's repeat-tract writer were never pushed through
# an external parser, and the interleaving of SNP and repeat-tract records in one file is new.
# This asks bcftools — the parser everything downstream is built on — four questions:
#
#   1. does the file parse at all, plain and bgzf, with no warnings?
#   2. do the two encodings carry identical content?
#   3. does the `STR` flag separate the two kinds of record in one filter expression?
#   4. can the file be indexed, and does a region query return both records of the one legal
#      position tie?
#
# Question 4 is the one that cannot be answered by reading: an ordering fault is invisible until
# something indexes the file.
#
# Run it inside the dev container:
#   ./scripts/dev.sh ./scripts/ng_vcf_parser_gate.sh
set -euo pipefail

WORK="${1:-tmp/vcf_gate}"
mkdir -p "$WORK"

echo "== building and writing the fixture =="
cargo run --release --quiet --example ng_vcf_fixture -- "$WORK"

PLAIN="$WORK/fixture.vcf"
BGZF="$WORK/fixture.vcf.gz"
failures=0

check() {
    local what="$1"
    shift
    if "$@"; then
        echo "  ok   — $what"
    else
        echo "  FAIL — $what"
        failures=$((failures + 1))
    fi
}

echo
echo "== 1. it parses, and says nothing while doing so =="
for file in "$PLAIN" "$BGZF"; do
    err="$WORK/$(basename "$file").err"
    bcftools view "$file" > "$file.parsed" 2> "$err" || {
        echo "  FAIL — bcftools view $file exited non-zero"; failures=$((failures + 1)); }
    if [[ -s "$err" ]]; then
        echo "  FAIL — bcftools wrote to stderr for $file:"
        sed 's/^/         /' "$err"
        failures=$((failures + 1))
    else
        echo "  ok   — $file parses with no warnings"
    fi
done

echo
echo "== 2. the two encodings carry the same content =="
# bcftools stamps its own command line into the header, naming the input file; that one line is
# expected to differ and nothing else may.
grep -v '^##bcftools_viewCommand' "$PLAIN.parsed" > "$WORK/plain.body"
grep -v '^##bcftools_viewCommand' "$BGZF.parsed" > "$WORK/bgzf.body"
check "plain and bgzf agree" diff -q "$WORK/plain.body" "$WORK/bgzf.body"

echo
echo "== 3. the STR flag separates the two kinds of record =="
tracts=$(bcftools view -H -i 'STR=1' "$PLAIN" | wc -l | tr -d ' ')
others=$(bcftools view -H -e 'STR=1' "$PLAIN" | wc -l | tr -d ' ')
total=$(bcftools view -H "$PLAIN" | wc -l | tr -d ' ')
echo "  $tracts repeat-tract records, $others SNP/indel records, $total in all"
check "the two selections partition the file" test "$((tracts + others))" -eq "$total"
check "both kinds are present" test "$tracts" -gt 0 -a "$others" -gt 0
# Every tract record must carry the motif and period beside the flag.
missing=$(bcftools query -i 'STR=1' -f '%CHROM\t%POS\t%INFO/RU\t%INFO/PERIOD\n' "$PLAIN" \
    | grep -c '\.' || true)
check "every STR record carries RU and PERIOD" test "$missing" -eq 0

echo
echo "== 4. it indexes, and the one legal position tie survives a region query =="
bcftools index -f -t "$BGZF"
tie=$(bcftools view -H -r chr1:100-100 "$BGZF" | wc -l | tr -d ' ')
echo "  $tie records at chr1:100 — a SNP and the repeat tract padded onto it"
check "both records of the tie come back" test "$tie" -eq 2

echo
if [[ "$failures" -eq 0 ]]; then
    echo "PASS — every check"
else
    echo "FAIL — $failures check(s)"
    exit 1
fi
