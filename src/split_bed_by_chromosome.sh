#!/usr/bin/env bash

set -euo pipefail

BED_INPUT="data/processed/gnomad_lookup.bed"
OUT_DIR="data/processed/bed_by_chr"
CHROMS=(chr1 chr2 chr3 chr6 chr7 chr11 chr14 chr19)

mkdir -p "$OUT_DIR"

for chrom in "${CHROMS[@]}"; do
  grep "^${chrom}" "$BED_INPUT" > "${OUT_DIR}/${chrom}.bed" || true
done
