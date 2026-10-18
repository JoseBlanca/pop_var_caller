#!/usr/bin/env bash
# Fetch the scoring tool the dev image lacks: RTG Tools 3.13 (vcfeval), into this checkout's
# tmp/tools/. Pinned by version and checksum. Does nothing if it is already there.
#
# RTG is run as `java -jar RTG.jar` inside the container: the `rtg` launcher refuses the
# container's arm64 processor. The no-JRE build is used because the container has its own Java.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
VERSION=3.13
URL="https://github.com/RealTimeGenomics/rtg-tools/releases/download/$VERSION/rtg-tools-$VERSION-nojre.zip"
SHA256="bc5c6badb07d7e20d1c5c557bd6d571a022bbd9f58fa1e3840bcff9431a18f96"
DEST="$PROJECT_DIR/tmp/tools"

[[ -f "$DEST/rtg-tools-$VERSION/RTG.jar" ]] && exit 0

mkdir -p "$DEST"
zip="$DEST/rtg-tools-$VERSION-nojre.zip"
echo "fetching RTG Tools $VERSION into $DEST"
curl -fsSL -o "$zip" "$URL"
got="$(shasum -a 256 "$zip" | awk '{print $1}')"
if [[ "$got" != "$SHA256" ]]; then
    echo "RTG Tools download has checksum $got, expected $SHA256; refusing it ($zip)" >&2
    exit 1
fi
(cd "$DEST" && unzip -q "$zip")
rm -f "$zip"
