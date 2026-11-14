#!/usr/bin/env bash
# Fetch opencore-amr 0.1.6 (Apache-2.0) into the gitignored .cache/opencore/: the C oracle for
# differential tests and the benchmark baseline. SHA-256 pinned and cross-checked against
# Debian's opencore-amr_0.1.6-1.dsc. Idempotent.
set -euo pipefail
cd "$(dirname "$0")/.."

URL=https://downloads.sourceforge.net/project/opencore-amr/opencore-amr/opencore-amr-0.1.6.tar.gz
SHA=483eb4061088e2b34b358e47540b5d495a96cd468e361050fae615b1809dc4a1
DIR=.cache/opencore
TGZ=$DIR/opencore-amr-0.1.6.tar.gz

mkdir -p "$DIR"
if ! echo "$SHA  $TGZ" | sha256sum -c --status 2>/dev/null; then
  echo "opencore: downloading $URL"
  curl -fsSL --retry 3 -o "$TGZ.part" "$URL"
  mv "$TGZ.part" "$TGZ"
  echo "$SHA  $TGZ" | sha256sum -c --status || { echo "opencore: SHA-256 MISMATCH for $TGZ" >&2; exit 1; }
fi
if [ ! -f "$DIR/.unpacked" ]; then
  rm -rf "$DIR/opencore-amr-0.1.6"
  tar xzf "$TGZ" -C "$DIR"
  touch "$DIR/.unpacked"
fi
echo "opencore: ok ($DIR/opencore-amr-0.1.6; sha256 verified)"
