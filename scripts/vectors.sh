#!/usr/bin/env bash
# Download the 3GPP test sequences and reference C code into the gitignored .cache/3gpp/, verify
# the pinned SHA-256 and unpack them. Idempotent: a verified download is not repeated.
#   AMR-NB: TS 26.074 vectors (nodtx, dtx, dtx2), TS 26.073 reference C (crates/ref26073-sys)
#   AMR-WB: TS 26.174 vectors (ts26174/T, ts26174/S), TS 26.173 reference C (crates/ref26173-sys)
# The WB pins are the latest versions (j00, Rel-19) as first downloaded on 2026-11-14 (3GPP is
# the only source, so trust on first use). 3GPP material is "all rights reserved": a test oracle
# only, never commit it, never port from it.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/fetch.sh

A=https://www.3gpp.org/ftp/Specs/archive/26_series
DIR=.cache/3gpp

fetch $A/26.074/26074-j00.zip e4d8771a807023738d7648e374913f0337c1c92e688bd5f4b572d44c5231e881 "$DIR/26074-j00.zip"
# Outer zip holds one zip per test set; each set holds one zip per sub-suite.
if [ ! -f "$DIR/.unpacked" ]; then
  unzip -oq "$DIR/26074-j00.zip" -d "$DIR"
  for set in NODTX DTX DTX2; do
    mkdir -p "$DIR/${set,,}"
    unzip -oq "$DIR/TS_AMR_h01_${set}.zip" -d "$DIR/${set,,}"
    for inner in "$DIR/${set,,}"/*.zip; do unzip -oq "$inner" -d "$DIR/${set,,}"; done
  done
  touch "$DIR/.unpacked"
fi

REF=$DIR/ref26073
fetch $A/26.073/26073-j00.zip 8e37c063f1aba04bd226d71595b09fa65804329d267b09f347b2eb94459046fa "$REF/26073-j00.zip"
if [ ! -f "$REF/.unpacked" ]; then
  unzip -oq "$REF/26073-j00.zip" -d "$REF"
  unzip -oq "$REF/26073-i00_ANSI_C_source_code.zip" -d "$REF"
  touch "$REF/.unpacked"
fi

WB=$DIR/ts26174
fetch $A/26.174/26174-j00.zip bb2eb62db39263bdf047b41d4a3693ad3785b559c29e116043a2cb4f082d908d "$WB/26174-j00.zip"
if [ ! -f "$WB/.unpacked" ]; then
  unzip -oq "$WB/26174-j00.zip" -d "$WB"
  unzip -oq "$WB/T.zip" -d "$WB"
  unzip -oq "$WB/S.zip" -d "$WB/S"
  touch "$WB/.unpacked"
fi

REFWB=$DIR/ref26173
fetch $A/26.173/26173-j00.zip 738251c4fba4aee33f7af856752d3081a3f54b20d1b977c2e04096ce0bb9ce4c "$REFWB/26173-j00.zip"
if [ ! -f "$REFWB/.unpacked" ]; then
  unzip -oq "$REFWB/26173-j00.zip" -d "$REFWB"
  unzip -oq "$REFWB/26173-i00_ANSI-C_source_code.zip" -d "$REFWB"
  touch "$REFWB/.unpacked"
fi

echo "vectors: ok ($DIR: nodtx, dtx, dtx2, ref26073; ts26174, ref26173; sha256 verified)"
