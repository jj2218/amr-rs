#!/usr/bin/env bash
# Fetch opencore-amr 0.1.6 and vo-amrwbenc 0.1.3 (both Apache-2.0, the opencore-amr project)
# into the gitignored .cache/opencore/: the C oracles for differential tests and the benchmark
# baselines (vo-amrwbenc: the AMR-WB encoder). SHA-256s pinned and cross-checked against Debian's
# opencore-amr_0.1.6-1.dsc and vo-amrwbenc_0.1.3-2.dsc. Idempotent.
set -euo pipefail
cd "$(dirname "$0")/.."

URL=https://downloads.sourceforge.net/project/opencore-amr/opencore-amr/opencore-amr-0.1.6.tar.gz
SHA=483eb4061088e2b34b358e47540b5d495a96cd468e361050fae615b1809dc4a1
DIR=.cache/opencore
TGZ=$DIR/opencore-amr-0.1.6.tar.gz

source scripts/fetch.sh
fetch "$URL" "$SHA" "$TGZ"
if [ ! -f "$DIR/.unpacked" ]; then
  rm -rf "$DIR/opencore-amr-0.1.6"
  tar xzf "$TGZ" -C "$DIR"
  touch "$DIR/.unpacked"
fi
VO_URL=https://downloads.sourceforge.net/project/opencore-amr/vo-amrwbenc/vo-amrwbenc-0.1.3.tar.gz
VO_SHA=5652b391e0f0e296417b841b02987d3fd33e6c0af342c69542cbb016a71d9d4e
VO_TGZ=$DIR/vo-amrwbenc-0.1.3.tar.gz
fetch "$VO_URL" "$VO_SHA" "$VO_TGZ"
if [ ! -f "$DIR/.vo-unpacked" ]; then
  rm -rf "$DIR/vo-amrwbenc-0.1.3"
  tar xzf "$VO_TGZ" -C "$DIR"
  touch "$DIR/.vo-unpacked"
fi
echo "opencore: ok ($DIR: opencore-amr-0.1.6, vo-amrwbenc-0.1.3; sha256 verified)"
