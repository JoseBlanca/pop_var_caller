#!/usr/bin/env bash
# Fetch the tools the dev image lacks into this checkout's tmp/tools/, each pinned by version and
# checksum, each skipped if already there:
#
# - RTG Tools 3.13 (vcfeval). Run as `java -jar RTG.jar` inside the container: the `rtg` launcher
#   refuses the container's arm64 processor. The no-JRE build, since the container has its own Java.
# - Beagle 4.1 (27Jan18.7e1), refine_sim's imputation competitor. 4.1, not 5: Beagle 5 takes hard
#   genotypes only, and 4.1 is the last to impute from genotype likelihoods (`gl=`).
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$PROJECT_DIR/tmp/tools"
mkdir -p "$DEST"

sha256() {
    if command -v sha256sum > /dev/null; then sha256sum "$1" | awk '{print $1}'
    else shasum -a 256 "$1" | awk '{print $1}'
    fi
}

fetch() {  # fetch <url> <file> <sha256>
    echo "fetching $(basename "$2") into $(dirname "$2")"
    curl -fsSL -o "$2.part" "$1"
    local got
    got="$(sha256 "$2.part")"
    if [[ "$got" != "$3" ]]; then
        echo "$(basename "$2") download has checksum $got, expected $3; refusing it ($2.part)" >&2
        exit 1
    fi
    mv "$2.part" "$2"
}

RTG_VERSION=3.13
if [[ ! -f "$DEST/rtg-tools-$RTG_VERSION/RTG.jar" ]]; then
    zip="$DEST/rtg-tools-$RTG_VERSION-nojre.zip"
    fetch "https://github.com/RealTimeGenomics/rtg-tools/releases/download/$RTG_VERSION/rtg-tools-$RTG_VERSION-nojre.zip" \
        "$zip" bc5c6badb07d7e20d1c5c557bd6d571a022bbd9f58fa1e3840bcff9431a18f96
    (cd "$DEST" && unzip -q "$zip")
    rm -f "$zip"
fi

BEAGLE_JAR=beagle.27Jan18.7e1.jar
if [[ ! -f "$DEST/beagle/$BEAGLE_JAR" ]]; then
    mkdir -p "$DEST/beagle"
    fetch "https://faculty.washington.edu/browning/beagle/$BEAGLE_JAR" "$DEST/beagle/$BEAGLE_JAR" \
        6c94610b278fc108c3e80b1134226911be1fc92b7d378ba648ac3eb97c5a3207
fi
