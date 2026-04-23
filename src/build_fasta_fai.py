from __future__ import annotations

import argparse
from pathlib import Path


def build_fai(fasta_path: Path) -> Path:
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
    with fasta_path.open("rb") as fasta, fai_path.open("w", encoding="ascii", newline="\n") as out:
        name: str | None = None
        seq_len = 0
        seq_offset = 0
        line_bases: int | None = None
        line_width: int | None = None

        while True:
            line = fasta.readline()
            if not line:
                if name is not None:
                    out.write(f"{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_width}\n")
                break

            if line.startswith(b">"):
                if name is not None:
                    out.write(f"{name}\t{seq_len}\t{seq_offset}\t{line_bases}\t{line_width}\n")

                header = line[1:].decode("ascii", errors="ignore").strip()
                name = header.split()[0]
                seq_len = 0
                seq_offset = fasta.tell()
                line_bases = None
                line_width = None
                continue

            stripped = line.rstrip(b"\r\n")
            if not stripped:
                continue
            if line_bases is None:
                line_bases = len(stripped)
                line_width = len(line)
            seq_len += len(stripped)

    return fai_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a FASTA .fai index without samtools.")
    parser.add_argument("fasta", type=Path, help="Path to FASTA file")
    args = parser.parse_args()

    fai_path = build_fai(args.fasta)
    print(f"Wrote {fai_path}")


if __name__ == "__main__":
    main()
