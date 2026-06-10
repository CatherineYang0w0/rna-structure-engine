from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


MISSING_TOKEN = "-999"


@dataclass(frozen=True)
class ConvertedProfile:
    name: str
    fasta: Path
    shape: Path
    vienna_shape: Path
    table: Path
    length: int


def _parse_float(value: str) -> float | None:
    text = value.strip()
    if not text or text.lower() == "nan":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def _clean_reactivity(value: float | None) -> float | None:
    if value is None:
        return None
    return max(0.0, value)


def _write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def convert_shapemapper_profile(
    profile: Path,
    outdir: Path,
    name: str,
    reactivity_column: str = "Norm_profile",
) -> ConvertedProfile:
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    with profile.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{profile} has no header")
        if reactivity_column not in reader.fieldnames:
            raise ValueError(f"{profile} is missing column {reactivity_column!r}")
        required = {"Nucleotide", "Sequence"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{profile} is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            rows.append(row)

    sequence = "".join(row["Sequence"].strip().upper().replace("T", "U") for row in rows)
    if not sequence:
        raise ValueError(f"{profile} produced an empty sequence")

    fasta = outdir / f"{name}.fa"
    shape = outdir / f"{name}.shape"
    vienna_shape = outdir / f"{name}.vienna.shape"
    table = outdir / f"{name}.converted.tsv"

    _write_fasta(fasta, name, sequence)
    with shape.open("w") as shape_handle, vienna_shape.open("w") as vienna_handle, table.open("w", newline="") as table_handle:
        writer = csv.writer(table_handle, delimiter="\t")
        writer.writerow(["position", "base", "raw_reactivity", "shape_reactivity"])
        for expected_pos, row in enumerate(rows, start=1):
            position = int(row["Nucleotide"])
            if position != expected_pos:
                raise ValueError(f"{profile} has non-contiguous positions at row {expected_pos}: {position}")
            base = row["Sequence"].strip().upper().replace("T", "U")
            raw = _parse_float(row[reactivity_column])
            clean = _clean_reactivity(raw)
            shape_value = MISSING_TOKEN if clean is None else f"{clean:.6f}"
            raw_value = "nan" if raw is None else f"{raw:.6f}"
            shape_handle.write(f"{position}\t{shape_value}\n")
            if clean is not None:
                vienna_handle.write(f"{position}\t{base}\t{clean:.6f}\n")
            writer.writerow([position, base, raw_value, shape_value])

    return ConvertedProfile(name=name, fasta=fasta, shape=shape, vienna_shape=vienna_shape, table=table, length=len(sequence))


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a ShapeMapper profile to FASTA and RNAstructure SHAPE files.")
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--reactivity-column", default="Norm_profile")
    args = parser.parse_args()

    converted = convert_shapemapper_profile(args.profile, args.outdir, args.name, args.reactivity_column)
    print(f"Wrote {converted.fasta}")
    print(f"Wrote {converted.shape}")
    print(f"Wrote {converted.vienna_shape}")
    print(f"Wrote {converted.table}")


if __name__ == "__main__":
    main()
