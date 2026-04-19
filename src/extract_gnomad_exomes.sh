#!/usr/bin/env bash

set -euo pipefail

OUT_DIR="data/processed/gnomad_by_chr"
BED_DIR="data/processed/bed_by_chr"
OUT_FILE="data/processed/gnomad_subset_with_header.tsv"
CHROMS=(1 2 3 6 7 11 14 19)

mkdir -p "$OUT_DIR"

for chrom in "${CHROMS[@]}"; do
  echo "Processing chr${chrom}..."
  bcftools query \
    -R "${BED_DIR}/chr${chrom}.bed" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\n' \
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/exomes/gnomad.exomes.v4.1.1.sites.chr${chrom}.vcf.bgz" \
    -o "${OUT_DIR}/chr${chrom}.tsv"
done

printf 'CHROM\tPOS\tREF\tALT\tAF\tAC\tAN\n' > "$OUT_FILE"
cat "${OUT_DIR}"/*.tsv >> "$OUT_FILE"
echo "Saved ${OUT_FILE}"
