#!/usr/bin/env bash
# Download the 3GPP AMR-NB test sequences (TS 26.074) and reference C code (TS 26.073) into the
# gitignored .cache/3gpp/, verify the pinned SHA-256 and unpack them. Idempotent: a verified
# download is not repeated.
# 3GPP material is "all rights reserved": a test oracle only, never commit it.
set -euo pipefail
cd "$(dirname "$0")/.."

URL=https://www.3gpp.org/ftp/Specs/archive/26_series/26.074/26074-j00.zip
SHA=e4d8771a807023738d7648e374913f0337c1c92e688bd5f4b572d44c5231e881
DIR=.cache/3gpp
ZIP=$DIR/26074-j00.zip

mkdir -p "$DIR"
if ! echo "$SHA  $ZIP" | sha256sum -c --status 2>/dev/null; then
  echo "vectors: downloading $URL"
  curl -fsSL --retry 3 -o "$ZIP.part" "$URL"
  mv "$ZIP.part" "$ZIP"
  echo "$SHA  $ZIP" | sha256sum -c --status || { echo "vectors: SHA-256 MISMATCH for $ZIP" >&2; exit 1; }
fi

# Outer zip holds one zip per test set; each set holds one zip per sub-suite.
if [ ! -f "$DIR/.unpacked" ]; then
  unzip -oq "$ZIP" -d "$DIR"
  for set in NODTX DTX DTX2; do
    mkdir -p "$DIR/${set,,}"
    unzip -oq "$DIR/TS_AMR_h01_${set}.zip" -d "$DIR/${set,,}"
    for inner in "$DIR/${set,,}"/*.zip; do unzip -oq "$inner" -d "$DIR/${set,,}"; done
  done
  touch "$DIR/.unpacked"
fi

# The TS 26.073 reference C code: the unit-test oracle (crates/ref26073-sys), never ported from.
REF_URL=https://www.3gpp.org/ftp/Specs/archive/26_series/26.073/26073-j00.zip
REF_SHA=8e37c063f1aba04bd226d71595b09fa65804329d267b09f347b2eb94459046fa
REF=$DIR/ref26073
REF_ZIP=$REF/26073-j00.zip
mkdir -p "$REF"
if ! echo "$REF_SHA  $REF_ZIP" | sha256sum -c --status 2>/dev/null; then
  echo "vectors: downloading $REF_URL"
  curl -fsSL --retry 3 -o "$REF_ZIP.part" "$REF_URL"
  mv "$REF_ZIP.part" "$REF_ZIP"
  echo "$REF_SHA  $REF_ZIP" | sha256sum -c --status || { echo "vectors: SHA-256 MISMATCH for $REF_ZIP" >&2; exit 1; }
fi
if [ ! -f "$REF/.unpacked" ]; then
  unzip -oq "$REF_ZIP" -d "$REF"
  unzip -oq "$REF/26073-i00_ANSI_C_source_code.zip" -d "$REF"
  touch "$REF/.unpacked"
fi

echo "vectors: ok ($DIR: nodtx, dtx, dtx2, ref26073 C code; sha256 verified)"
