#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC_DIR="$ROOT_DIR/docs/diagrams/src"
OUT_DIR="$ROOT_DIR/docs/diagrams/generated"
CACHE_ROOT="${XDG_CACHE_HOME:-$HOME/.cache}/chem-spectra-app"
PNG_DIAGRAMS=(
  "architecture-layers"
  "request-flow"
  "system-overview"
)

export NPM_CONFIG_CACHE="${NPM_CONFIG_CACHE:-$CACHE_ROOT/npm}"
export NPM_CONFIG_AUDIT=false
export NPM_CONFIG_FUND=false
export NPM_CONFIG_UPDATE_NOTIFIER=false
export PUPPETEER_CACHE_DIR="${PUPPETEER_CACHE_DIR:-$CACHE_ROOT/puppeteer}"

mkdir -p "$OUT_DIR"

for file in "$SRC_DIR"/*.mmd; do
  name="$(basename "${file%.mmd}")"
  npx --yes --package @mermaid-js/mermaid-cli mmdc -i "$file" -o "$OUT_DIR/$name.svg"
done

for name in "${PNG_DIAGRAMS[@]}"; do
  npx --yes --package @mermaid-js/mermaid-cli mmdc -i "$SRC_DIR/$name.mmd" -o "$OUT_DIR/$name.png"
done
