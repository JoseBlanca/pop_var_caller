#!/usr/bin/env bash
#
# What each route to a census costs — wall time and peak resident memory.
#
#   scripts/ng_census_route_cost.sh <reference.fa> <catalog.parquet> <regions.bed> <cram-or-dir>
#
# A census can be built during the walk over the reads, which is what generate-psps
# does, or afterwards from the stored psp, which is what generate-census does. The two
# produce the same file byte for byte, so what separates them is what they cost. This
# runs each route in a process of its own — peak resident memory of two routes in one
# process is the larger of them and says nothing about either — and then checks that the
# censuses the two wrote really are identical, because a timing comparison between two
# different outputs would mean nothing.
#
# NG_SAMPLES and NG_REGIONS pass through to the harness (how many alignment files and
# how many BED intervals). Everything is written under the repository's own tmp/.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if (( $# != 4 )); then
    echo "usage: ng_census_route_cost.sh <reference.fa> <catalog.parquet> <regions.bed> <cram-or-dir>" >&2
    exit 2
fi
REFERENCE="$1"; CATALOG="$2"; BED="$3"; CRAMS="$4"

# The release binary, from whichever target tree holds the newer one: the container build
# points CARGO_TARGET_DIR at target-container, a host build does not, and a machine with no
# container runtime only has the second.
BIN=""
for candidate in "$REPO/target-container/release/examples/ng_census_route_cost" \
                 "$REPO/target/release/examples/ng_census_route_cost"; do
    if [[ -x "$candidate" ]] && { [[ -z "$BIN" ]] || [[ "$candidate" -nt "$BIN" ]]; }; then
        BIN="$candidate"
    fi
done
if [[ -z "$BIN" ]]; then
    echo "no ng_census_route_cost binary; build it with:" >&2
    echo "  $REPO/scripts/dev.sh cargo build --release --example ng_census_route_cost" >&2
    exit 1
fi

OUT="$REPO/tmp/ng_census_route_cost"
rm -rf "$OUT"
mkdir -p "$OUT"

for route in during-the-walk after-the-walk; do
    NG_WORK="$OUT/$route" "$REPO/scripts/peak_rss.sh" "$OUT/$route.measure" \
        "$BIN" "$route" "$REFERENCE" "$CATALOG" "$BED" "$CRAMS" > "$OUT/$route.log" 2>&1
    status=$?
    if (( status != 0 )); then
        echo "the $route route failed; its log:" >&2
        cat "$OUT/$route.log" >&2
        exit "$status"
    fi
done

echo
echo "route              wall_s   peak_rss_mb"
for route in during-the-walk after-the-walk; do
    read -r rss wall < "$OUT/$route.measure"
    # The harness's own clock over the work, which excludes the setup both routes share.
    working=$(awk -F'seconds=' '/^route=/ {print $2}' "$OUT/$route.log")
    printf "%-18s %6s   %11s   (wrapper wall %ss)\n" "$route" "$working" "$rss" "$wall"
done

echo
echo "do the two routes write the same censuses?"
same=1
for census in "$OUT/during-the-walk"/*.census; do
    [[ -e "$census" ]] || { echo "  the during-the-walk route wrote no census"; same=0; break; }
    other="$OUT/after-the-walk/$(basename "$census")"
    if [[ ! -e "$other" ]]; then
        echo "  $(basename "$census"): the after-the-walk route wrote none"
        same=0
    elif cmp -s "$census" "$other"; then
        echo "  $(basename "$census"): identical"
    else
        echo "  $(basename "$census"): DIFFERENT — the timings above compare two different things"
        same=0
    fi
done
(( same == 1 )) || exit 1
