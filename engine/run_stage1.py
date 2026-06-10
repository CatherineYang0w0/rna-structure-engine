from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import struct
import subprocess
import zlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from engine.adapters.shapemapper_profile import ConvertedProfile, convert_shapemapper_profile
from engine.metrics import (
    compare_pair_sets,
    ct_pairs,
    dotbracket_pairs,
    parse_probability_plot_text,
    parse_vienna_dp_ps,
    read_dbn_structure,
    read_fasta_sequence,
    rolling_mean,
    shannon_entropy,
    write_entropy,
)
from engine.run_rnastructure import ct_to_dotbracket, draw_structure, require_tool, run_command


ROOT = Path(__file__).resolve().parents[1]
SCHEMES = ("rnastr-deigan", "rnafold-deigan", "rnafold-zarringhalam", "rnafold-washietl")


@dataclass(frozen=True)
class SchemeResult:
    scheme: str
    outdir: Path
    mfe_dbn: Path
    mea_dbn: Path | None
    entropy_tsv: Path
    energy_tsv: Path


def _write_dbn(path: Path, name: str, sequence: str, structure: str) -> None:
    path.write_text(f">{name}\n{sequence}\n{structure}\n")


def _read_shape_values(path: Path) -> list[tuple[int, str]]:
    values: list[tuple[int, str]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                values.append((int(fields[0]), fields[1]))
    return values


def _write_energy_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = ["scheme", "structure", "energy_type", "energy_kcal_mol", "note"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _run_rnastructure(converted: ConvertedProfile, outdir: Path, log_path: Path) -> SchemeResult:
    outdir.mkdir(parents=True, exist_ok=True)
    fold = require_tool("Fold")
    partition = require_tool("partition")
    probability_plot = require_tool("ProbabilityPlot")
    max_expect = require_tool("MaxExpect")
    efn2 = require_tool("efn2")

    scheme = "rnastr-deigan"
    ct = outdir / "mfe.ct"
    dbn = outdir / "mfe.dbn"
    pfs = outdir / "partition.pfs"
    probability_txt = outdir / "pairing_probability.txt"
    probability_svg = outdir / "pairing_probability.svg"
    entropy_tsv = outdir / "shannon_entropy.tsv"
    mea_ct = outdir / "mea.ct"
    mea_dbn = outdir / "mea.dbn"
    mfe_energy = outdir / "mfe_efn2.txt"
    mea_energy = outdir / "mea_efn2.txt"
    energy_tsv = outdir / "energy_summary.tsv"
    structure_svg = outdir / "mfe_structure.svg"
    structure_png = outdir / "mfe_structure.png"

    slope = "1.8"
    intercept = "-0.6"
    run_command(
        [
            fold,
            str(converted.fasta),
            str(ct),
            "--SHAPE",
            str(converted.shape),
            "--SHAPEslope",
            slope,
            "--SHAPEintercept",
            intercept,
        ],
        log_path,
    )
    ct_to_dotbracket(ct, dbn)
    run_command(
        [
            partition,
            str(converted.fasta),
            str(pfs),
            "--SHAPE",
            str(converted.shape),
            "--SHAPEslope",
            slope,
            "--SHAPEintercept",
            intercept,
            "--quiet",
        ],
        log_path,
    )
    run_command([probability_plot, str(pfs), str(probability_txt), "--text"], log_path)
    run_command([probability_plot, str(pfs), str(probability_svg), "--svg"], log_path)
    run_command([max_expect, str(pfs), str(mea_ct), "--structures", "1"], log_path)
    ct_to_dotbracket(mea_ct, mea_dbn)
    run_command([efn2, str(ct), str(mfe_energy), "--SHAPE", str(converted.shape), "--SHAPEslope", slope, "--SHAPEintercept", intercept, "--quiet"], log_path)
    run_command([efn2, str(mea_ct), str(mea_energy), "--SHAPE", str(converted.shape), "--SHAPEslope", slope, "--SHAPEintercept", intercept, "--quiet"], log_path)

    sequence = read_fasta_sequence(converted.fasta)
    entropy = shannon_entropy(len(sequence), parse_probability_plot_text(probability_txt))
    write_entropy(entropy_tsv, entropy, rolling_mean(entropy, 55))
    draw_structure(ct, structure_svg, structure_png, log_path)
    _write_energy_tsv(
        energy_tsv,
        [
            {
                "scheme": scheme,
                "structure": "mfe",
                "energy_type": "efn2",
                "energy_kcal_mol": _first_energy(mfe_energy),
                "note": "RNAstructure efn2 with SHAPE pseudoenergies",
            },
            {
                "scheme": scheme,
                "structure": "mea",
                "energy_type": "efn2",
                "energy_kcal_mol": _first_energy(mea_energy),
                "note": "RNAstructure efn2 with SHAPE pseudoenergies",
            },
        ],
    )
    return SchemeResult(scheme, outdir, dbn, mea_dbn, entropy_tsv, energy_tsv)


def _first_energy(path: Path) -> str:
    text = path.read_text()
    match = re.search(r"(-?\d+(?:\.\d+)?)", text)
    return match.group(1) if match else ""


def _run_capture(cmd: list[str], log_path: Path, cwd: Path, stdin_text: str | None = None) -> subprocess.CompletedProcess[str]:
    with log_path.open("a") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
    proc = subprocess.run(cmd, input=stdin_text, cwd=cwd, text=True, capture_output=True, check=True)
    with log_path.open("a") as log:
        if proc.stdout:
            log.write(proc.stdout)
        if proc.stderr:
            log.write(proc.stderr)
    return proc


def _rnafold_shape_args(
    scheme: str,
    converted: ConvertedProfile,
    outdir: Path,
    log_path: Path,
    sequence_text: str,
    washietl_sample_size: int,
) -> list[str]:
    if scheme == "rnafold-deigan":
        return [f"--shape={converted.shape}", '--shapeMethod=Dm1.8b-0.6']
    if scheme == "rnafold-zarringhalam":
        return [f"--shape={converted.shape}", "--shapeMethod=Z"]
    if scheme == "rnafold-washietl":
        rnapvmin = require_tool("RNApvmin")
        perturbation = outdir / "washietl_perturbation_vector.shape"
        cmd = [rnapvmin]
        if washietl_sample_size > 0:
            cmd.append(f"--sampleSize={washietl_sample_size}")
        cmd.append(str(converted.vienna_shape))
        proc = _run_capture(cmd, log_path, outdir, stdin_text=sequence_text)
        perturbation.write_text(proc.stdout)
        return [f"--shape={perturbation}", "--shapeMethod=W"]
    raise ValueError(f"Unsupported RNAfold scheme: {scheme}")


def _parse_rnafold_stdout(stdout: str) -> dict[str, str]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    result: dict[str, str] = {}
    structure_lines = [line for line in lines if line and line[0] in ".([{<"]
    if structure_lines:
        result["mfe_structure"] = structure_lines[0].split()[0]
        result["mfe_energy"] = _extract_energy(structure_lines[0], "paren")
    if len(structure_lines) > 1:
        result["ensemble_energy"] = _extract_energy(structure_lines[1], "bracket")
    for line in structure_lines:
        if "MEA=" in line:
            result["mea_structure"] = line.split()[0]
            result["mea_score"] = line.split("MEA=", 1)[1].split("}", 1)[0].strip()
    for line in lines:
        if line.startswith("frequency of mfe structure"):
            result["ensemble_note"] = line
    return result


def _extract_energy(line: str, kind: str) -> str:
    if kind == "paren":
        pattern = r"\(\s*(-?\d+(?:\.\d+)?)\)"
    elif kind == "bracket":
        pattern = r"\[\s*(-?\d+(?:\.\d+)?)\]"
    else:
        pattern = r"(-?\d+(?:\.\d+)?)"
    match = re.search(pattern, line)
    return match.group(1) if match else ""


def _run_rnafold_scheme(
    scheme: str,
    converted: ConvertedProfile,
    outdir: Path,
    log_path: Path,
    washietl_sample_size: int,
) -> SchemeResult:
    outdir.mkdir(parents=True, exist_ok=True)
    rnafold = require_tool("RNAfold")
    sequence = read_fasta_sequence(converted.fasta)
    header = f"{scheme}_{converted.name}"
    sequence_text = f">{header}\n{sequence}\n"
    shape_args = _rnafold_shape_args(scheme, converted, outdir, log_path, sequence_text, washietl_sample_size)
    stdout_path = outdir / "rnafold.stdout.txt"
    cmd = [rnafold, "-p", "--MEA", "--noPS", *shape_args]
    proc = _run_capture(cmd, log_path, outdir, stdin_text=sequence_text)
    stdout_path.write_text(proc.stdout)

    parsed = _parse_rnafold_stdout(proc.stdout)
    mfe_dbn = outdir / "mfe.dbn"
    mea_dbn = outdir / "mea.dbn"
    _write_dbn(mfe_dbn, f"{header}_mfe", sequence, parsed["mfe_structure"])
    _write_dbn(mea_dbn, f"{header}_mea", sequence, parsed.get("mea_structure", parsed["mfe_structure"]))

    dp_ps_candidates = sorted(outdir.glob("*_dp.ps"))
    if not dp_ps_candidates:
        raise FileNotFoundError(f"RNAfold did not produce a *_dp.ps file in {outdir}")
    dp_ps = dp_ps_candidates[0]
    probability_txt = outdir / "pairing_probability.tsv"
    _write_vienna_probability_tsv(dp_ps, probability_txt)
    entropy_tsv = outdir / "shannon_entropy.tsv"
    entropy = shannon_entropy(len(sequence), parse_vienna_dp_ps(dp_ps))
    write_entropy(entropy_tsv, entropy, rolling_mean(entropy, 55))

    energy_tsv = outdir / "energy_summary.tsv"
    _write_energy_tsv(
        energy_tsv,
        [
            {
                "scheme": scheme,
                "structure": "mfe",
                "energy_type": "RNAfold_MFE",
                "energy_kcal_mol": parsed.get("mfe_energy", ""),
                "note": "ViennaRNA RNAfold stdout",
            },
            {
                "scheme": scheme,
                "structure": "ensemble",
                "energy_type": "RNAfold_ensemble",
                "energy_kcal_mol": parsed.get("ensemble_energy", ""),
                "note": parsed.get("ensemble_note", "ViennaRNA RNAfold partition ensemble energy"),
            },
            {
                "scheme": scheme,
                "structure": "mea",
                "energy_type": "RNAfold_MEA_score",
                "energy_kcal_mol": parsed.get("mea_score", ""),
                "note": "MEA score, not an energy",
            },
        ],
    )
    return SchemeResult(scheme, outdir, mfe_dbn, mea_dbn, entropy_tsv, energy_tsv)


def _write_vienna_probability_tsv(dp_ps: Path, output: Path) -> None:
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["i", "j", "probability"])
        seen: set[tuple[int, int]] = set()
        for i, entries in parse_vienna_dp_ps(dp_ps).items():
            for j, probability in entries:
                pair = tuple(sorted((i, j)))
                if pair in seen:
                    continue
                seen.add(pair)
                writer.writerow([pair[0], pair[1], f"{probability:.9f}"])


def _write_applicability(path: Path, washietl_numeric_warning: bool) -> None:
    rows = []
    for scheme in SCHEMES:
        for algorithm in ("MFE/Fold", "partition", "MEA", "efn2_or_ensemble_energy", "Shannon_entropy", "ProbKnot"):
            status = "done"
            reason = ""
            if scheme == "rnafold-washietl" and washietl_numeric_warning and algorithm in {
                "partition",
                "MEA",
                "efn2_or_ensemble_energy",
                "Shannon_entropy",
            }:
                status = "done_with_numeric_warning"
                reason = "RNApvmin/RNAfold Washietl ensemble output contains NaN with the bounded sample-size setting; MFE structure is still emitted."
            if algorithm == "ProbKnot" and scheme == "rnastr-deigan":
                status = "deferred_to_stage_3"
                reason = "ProbKnot is an RNAstructure partition consumer and is planned for the pseudoknot phase."
            elif algorithm == "ProbKnot":
                status = "N/A"
                reason = "RNAfold cannot represent pseudoknots."
            rows.append({"scheme": scheme, "algorithm": algorithm, "status": status, "reason": reason})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scheme", "algorithm", "status", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _existing_result(scheme: str, outdir: Path) -> SchemeResult | None:
    mfe_dbn = outdir / "mfe.dbn"
    mea_dbn = outdir / "mea.dbn"
    entropy_tsv = outdir / "shannon_entropy.tsv"
    energy_tsv = outdir / "energy_summary.tsv"
    if all(path.exists() and path.stat().st_size > 0 for path in (mfe_dbn, mea_dbn, entropy_tsv, energy_tsv)):
        return SchemeResult(scheme, outdir, mfe_dbn, mea_dbn, entropy_tsv, energy_tsv)
    return None


def _compare(results: list[SchemeResult], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    pair_sets: dict[tuple[str, str], set[tuple[int, int]]] = {}
    for result in results:
        mfe_pairs = dotbracket_pairs(read_dbn_structure(result.mfe_dbn))
        if mfe_pairs:
            pair_sets[(result.scheme, "mfe")] = mfe_pairs
        if result.mea_dbn:
            mea_pairs = dotbracket_pairs(read_dbn_structure(result.mea_dbn))
            if mea_pairs:
                pair_sets[(result.scheme, "mea")] = mea_pairs

    rows: list[dict[str, str]] = []
    keys = sorted(pair_sets)
    for index, left in enumerate(keys):
        for right in keys[index + 1 :]:
            metrics = compare_pair_sets(pair_sets[left], pair_sets[right])
            rows.append(
                {
                    "left_scheme": left[0],
                    "left_structure": left[1],
                    "right_scheme": right[0],
                    "right_structure": right[1],
                    **{key: f"{value:.6f}" if isinstance(value, float) else str(value) for key, value in metrics.items()},
                }
            )
    tsv = outdir / "pairwise_structure_consistency.tsv"
    with tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_scheme",
                "left_structure",
                "right_scheme",
                "right_structure",
                "pairs_a",
                "pairs_b",
                "shared_pairs",
                "jaccard",
                "a_overlap",
                "b_overlap",
                "base_pair_distance",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    _write_consistency_matrix(outdir / "pairwise_mfe_jaccard.svg", outdir / "pairwise_mfe_jaccard.png", results)


def _has_numeric_warning(path: Path) -> bool:
    if not path.exists():
        return False
    return "nan" in path.read_text().lower()


def _write_consistency_matrix(svg: Path, png: Path, results: list[SchemeResult]) -> None:
    schemes = [result.scheme for result in results]
    pairs = {
        result.scheme: dotbracket_pairs(read_dbn_structure(result.mfe_dbn))
        for result in results
        if dotbracket_pairs(read_dbn_structure(result.mfe_dbn))
    }
    schemes = [scheme for scheme in schemes if scheme in pairs]
    matrix = [[compare_pair_sets(pairs[a], pairs[b])["jaccard"] for b in schemes] for a in schemes]
    cell = 74
    left = 210
    top = 50
    width = left + cell * len(schemes) + 30
    height = top + cell * len(schemes) + 120
    rects = []
    for r, row in enumerate(matrix):
        for c, value in enumerate(row):
            red = round(255 * (1 - float(value)))
            green = round(255 * float(value))
            color = f"#{red:02x}{green:02x}70"
            x = left + c * cell
            y = top + r * cell
            rects.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{color}" stroke="#ffffff"/>')
            rects.append(f'<text x="{x + cell / 2:.1f}" y="{y + cell / 2 + 5:.1f}" text-anchor="middle" font-family="Arial" font-size="13">{float(value):.2f}</text>')
    labels = []
    for i, scheme in enumerate(schemes):
        labels.append(f'<text x="{left - 8}" y="{top + i * cell + cell / 2 + 5:.1f}" text-anchor="end" font-family="Arial" font-size="12">{scheme}</text>')
        labels.append(f'<text x="{left + i * cell + cell / 2:.1f}" y="{top + len(schemes) * cell + 18}" text-anchor="middle" font-family="Arial" font-size="12" transform="rotate(35 {left + i * cell + cell / 2:.1f} {top + len(schemes) * cell + 18})">{scheme}</text>')
    svg.write_text(
        f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{left}" y="28" font-family="Arial" font-size="18">MFE base-pair Jaccard consistency</text>
{''.join(rects)}
{''.join(labels)}
</svg>
"""
    )
    _write_matrix_png(png, matrix)


def _write_matrix_png(path: Path, matrix: list[list[float]]) -> None:
    cell = 90
    width = cell * len(matrix)
    height = cell * len(matrix)
    data = bytearray([255, 255, 255] * width * height)
    for r, row in enumerate(matrix):
        for c, value in enumerate(row):
            red = round(255 * (1 - float(value)))
            green = round(255 * float(value))
            for y in range(r * cell, (r + 1) * cell):
                for x in range(c * cell, (c + 1) * cell):
                    offset = (y * width + x) * 3
                    data[offset : offset + 3] = bytes((red, green, 112))
    rows = [b"\x00" + data[y * width * 3 : (y + 1) * width * 3] for y in range(height)]

    def chunk(kind: bytes, payload: bytes) -> bytes:
        return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)

    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(b"".join(rows), level=9))
        + chunk(b"IEND", b"")
    )


def run_stage1(case: str, test: str, reagent: str, washietl_sample_size: int) -> None:
    if case != "xist" or test != "test-1":
        raise ValueError("Stage 1 currently supports CASE=xist and TEST=test-1")
    raw_profile = ROOT / "inputs" / test / case / "raw" / f"Ecoli_cellfree_16S_{reagent}_profile.txt"
    if not raw_profile.exists():
        raise FileNotFoundError(f"Missing input profile: {raw_profile}")

    root_out = ROOT / "outputs" / test / case
    stage_out = root_out / "stage1"
    prepared_dir = stage_out / "prepared"
    schemes_dir = stage_out / "schemes"
    comparisons_dir = stage_out / "comparisons"
    stage_out.mkdir(parents=True, exist_ok=True)
    log_path = stage_out / "run.log"
    if log_path.exists():
        log_path.unlink()

    converted = convert_shapemapper_profile(raw_profile, prepared_dir, f"ecoli_16s_cellfree_{reagent.lower()}")
    results: list[SchemeResult] = []
    rnastr_dir = schemes_dir / "rnastr-deigan"
    results.append(_existing_result("rnastr-deigan", rnastr_dir) or _run_rnastructure(converted, rnastr_dir, log_path))
    for scheme in ("rnafold-deigan", "rnafold-zarringhalam", "rnafold-washietl"):
        scheme_dir = schemes_dir / scheme
        results.append(
            _existing_result(scheme, scheme_dir)
            or _run_rnafold_scheme(scheme, converted, scheme_dir, log_path, washietl_sample_size)
        )
    _compare(results, comparisons_dir)
    washietl_warning = _has_numeric_warning((schemes_dir / "rnafold-washietl" / "rnafold.stdout.txt"))
    _write_applicability(stage_out / "algorithm_applicability.tsv", washietl_warning)
    (stage_out / "manifest.json").write_text(
        json.dumps(
            {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "stage": 1,
                "case": case,
                "test": test,
                "actual_dataset": "Busan et al. 2019 E. coli 16S rRNA cell-free SHAPE-MaP",
                "reagent": reagent,
                "schemes": SCHEMES,
                "washietl_sample_size": washietl_sample_size,
                "washietl_numeric_warning": washietl_warning,
                "input_profile": str(raw_profile.relative_to(ROOT)),
                "output_root": str(stage_out.relative_to(ROOT)),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"Wrote stage-1 outputs to {stage_out}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run stage-1 horizontal scheme expansion.")
    parser.add_argument("--case", default="xist")
    parser.add_argument("--test", default="test-1")
    parser.add_argument("--reagent", default="1M7")
    parser.add_argument("--washietl-sample-size", type=int, default=5)
    args = parser.parse_args()
    run_stage1(args.case, args.test, args.reagent, args.washietl_sample_size)


if __name__ == "__main__":
    main()
