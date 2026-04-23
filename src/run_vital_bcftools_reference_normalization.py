from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
REFERENCE_DIR = BASE_DIR / "data" / "external" / "reference"

MATCHED_IN = DATA_DIR / "arrhythmia_gnomad_matched.csv"
REFERENCE_FASTA = REFERENCE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REFERENCE_FASTA_GZ = REFERENCE_DIR / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

TEMP_VCF = DATA_DIR / "vital_bcftools_input.vcf"
NORM_VCF = DATA_DIR / "vital_bcftools_normalized.vcf"
DECOMP_VCF = DATA_DIR / "vital_bcftools_decomposed.vcf"

NORMALIZED_OUT = DATA_DIR / "vital_bcftools_reference_normalization.csv"
DECOMP_OUT = DATA_DIR / "vital_bcftools_decomposed_components.csv"
SUMMARY_OUT = DATA_DIR / "vital_bcftools_reference_normalization_summary.json"

WSL_PREFIX = "/mnt"


def bash_quote(text: str) -> str:
    return "'" + text.replace("'", "'\"'\"'") + "'"


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    parts = [part for part in resolved.parts[1:] if part not in {"\\", "/"}]
    return f"{WSL_PREFIX}/{drive}/" + "/".join(parts)


def run_wsl(script: str) -> None:
    completed = subprocess.run(
        ["wsl", "bash", "-lc", script],
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"WSL command failed ({completed.returncode}).\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )


def load_variants() -> pd.DataFrame:
    df = pd.read_csv(MATCHED_IN, dtype=str)
    df["pos"] = pd.to_numeric(df["pos"], errors="raise").astype(int)
    if df["variant_key"].duplicated().any():
        dupes = df.loc[df["variant_key"].duplicated(), "variant_key"].tolist()[:5]
        raise ValueError(f"Duplicate variant_key values in {MATCHED_IN.name}: {dupes}")
    return df


def write_input_vcf(df: pd.DataFrame) -> None:
    TEMP_VCF.parent.mkdir(parents=True, exist_ok=True)
    chroms = sorted(df["chrom"].astype(str).unique(), key=lambda value: (len(value), value))
    with TEMP_VCF.open("w", encoding="utf-8", newline="\n") as out:
        out.write("##fileformat=VCFv4.2\n")
        for chrom in chroms:
            out.write(f"##contig=<ID={chrom}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in df.itertuples(index=False):
            out.write(
                f"{row.chrom}\t{row.pos}\t{row.variant_key}\t{row.ref}\t{row.alt}\t.\t.\t.\n"
            )


def parse_vcf(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, pos, variant_key, ref, alt, qual, flt, info, *rest = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "variant_key": variant_key,
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                }
            )
    return rows


def main() -> None:
    if not REFERENCE_FASTA.exists():
        raise SystemExit(f"Missing reference FASTA: {REFERENCE_FASTA}")
    if not REFERENCE_FASTA_GZ.exists():
        raise SystemExit(f"Missing reference FASTA archive: {REFERENCE_FASTA_GZ}")
    if not REFERENCE_FASTA.with_suffix(REFERENCE_FASTA.suffix + ".fai").exists():
        raise SystemExit(
            f"Missing FASTA index: {REFERENCE_FASTA.with_suffix(REFERENCE_FASTA.suffix + '.fai')}"
        )

    df = load_variants()
    write_input_vcf(df)

    ref_wsl = to_wsl_path(REFERENCE_FASTA)
    in_wsl = to_wsl_path(TEMP_VCF)
    norm_wsl = to_wsl_path(NORM_VCF)
    decomp_wsl = to_wsl_path(DECOMP_VCF)

    run_wsl(
        "set -euo pipefail; "
        f"bcftools norm -f {bash_quote(ref_wsl)} -Ov -o {bash_quote(norm_wsl)} {bash_quote(in_wsl)}"
    )
    run_wsl(
        "set -euo pipefail; "
        f"bcftools norm -f {bash_quote(ref_wsl)} -m -both -Ov -o {bash_quote(decomp_wsl)} {bash_quote(in_wsl)}"
    )

    original = df[["variant_key", "chrom", "pos", "ref", "alt"]].copy()
    original = original.rename(
        columns={
            "chrom": "original_chrom",
            "pos": "original_pos",
            "ref": "original_ref",
            "alt": "original_alt",
        }
    )

    normalized = pd.DataFrame(parse_vcf(NORM_VCF))
    normalized = normalized.rename(
        columns={
            "chrom": "normalized_chrom",
            "pos": "normalized_pos",
            "ref": "normalized_ref",
            "alt": "normalized_alt",
        }
    )
    normalized = original.merge(normalized, on="variant_key", how="left", validate="one_to_one")
    normalized["representation_changed"] = (
        normalized["original_chrom"].astype(str) != normalized["normalized_chrom"].astype(str)
    ) | (
        normalized["original_pos"].astype(int) != normalized["normalized_pos"].astype(int)
    ) | (
        normalized["original_ref"].astype(str) != normalized["normalized_ref"].astype(str)
    ) | (
        normalized["original_alt"].astype(str) != normalized["normalized_alt"].astype(str)
    )

    components = pd.DataFrame(parse_vcf(DECOMP_VCF))
    if components.empty:
        components["component_index"] = pd.Series(dtype=int)
    else:
        components["component_index"] = components.groupby("variant_key").cumcount() + 1
    components = components.rename(
        columns={
            "chrom": "component_chrom",
            "pos": "component_pos",
            "ref": "component_ref",
            "alt": "component_alt",
        }
    )

    NORMALIZED_OUT.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(NORMALIZED_OUT, index=False, lineterminator="\n")
    components.to_csv(DECOMP_OUT, index=False, lineterminator="\n")

    summary = {
        "reference_fasta": str(REFERENCE_FASTA),
        "reference_archive": str(REFERENCE_FASTA_GZ),
        "total_variants": int(len(df)),
        "reference_normalized_variants": int(normalized["representation_changed"].sum()),
        "decomposed_rows": int(len(components)),
        "decomposed_variants": int(components["variant_key"].nunique()) if not components.empty else 0,
    }
    SUMMARY_OUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
