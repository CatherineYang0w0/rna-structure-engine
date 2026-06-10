from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from engine.adapters.shapemapper_profile import convert_shapemapper_profile
from engine.metrics import parse_probability_plot_text, read_fasta_sequence, rolling_mean, shannon_entropy, write_entropy
from engine.plotting import plot_reactivity_entropy


ROOT = Path(__file__).resolve().parents[1]


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"Required RNAstructure executable not found in PATH: {name}")
    return path


def run_command(cmd: list[str], log_path: Path) -> None:
    with log_path.open("a") as log:
        log.write("$ " + " ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT, text=True)


def ct_to_dotbracket(ct: Path, output: Path) -> None:
    bases: list[str] = []
    pairs: list[int] = []
    with ct.open() as handle:
        header = handle.readline().strip()
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            bases.append(fields[1].replace("T", "U"))
            pairs.append(int(fields[4]))
    dots = ["."] * len(bases)
    for index, paired in enumerate(pairs, start=1):
        if paired > index:
            dots[index - 1] = "("
            dots[paired - 1] = ")"
    output.write_text(f">{header}\n{''.join(bases)}\n{''.join(dots)}\n")


def draw_structure(ct: Path, svg: Path, png: Path, log_path: Path) -> None:
    draw = require_tool("draw")
    run_command([draw, str(ct), str(svg), "--svg", "--number", "1"], log_path)
    if shutil.which("convert"):
        run_command(["convert", str(svg), str(png)], log_path)
    elif shutil.which("magick"):
        run_command(["magick", str(svg), str(png)], log_path)
    else:
        with log_path.open("a") as log:
            log.write("WARNING: no SVG-to-PNG converter found; structure PNG was not generated.\n")


def run_stage0(case: str, test: str, reagent: str, scheme: str) -> None:
    if scheme != "rnastr-deigan":
        raise ValueError("Stage 0 only supports SCHEME=rnastr-deigan")
    if case != "xist" or test != "test-1":
        raise ValueError("Stage 0 currently supports CASE=xist and TEST=test-1")

    raw_profile = ROOT / "inputs" / test / case / "raw" / f"Ecoli_cellfree_16S_{reagent}_profile.txt"
    if not raw_profile.exists():
        raise FileNotFoundError(f"Missing input profile: {raw_profile}")

    outdir = ROOT / "outputs" / test / case
    prepared_dir = outdir / "prepared"
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "run.log"
    if log_path.exists():
        log_path.unlink()

    name = f"ecoli_16s_cellfree_{reagent.lower()}_{scheme}"
    converted = convert_shapemapper_profile(raw_profile, prepared_dir, name)

    fold = require_tool("Fold")
    partition = require_tool("partition")
    probability_plot = require_tool("ProbabilityPlot")

    ct = outdir / f"{name}.ct"
    dbn = outdir / f"{name}.dbn"
    pfs = outdir / f"{name}.pfs"
    probability_txt = outdir / f"{name}.pairing_probability.txt"
    probability_svg = outdir / f"{name}.pairing_probability.svg"
    entropy_tsv = outdir / f"{name}.shannon_entropy.tsv"
    profile_png = outdir / f"{name}.reactivity_entropy.png"
    profile_svg = outdir / f"{name}.reactivity_entropy.svg"
    structure_svg = outdir / f"{name}.mfe_structure.svg"
    structure_png = outdir / f"{name}.mfe_structure.png"
    manifest = outdir / "manifest.json"

    # Deigan SHAPE pseudoenergy parameters requested for this engine:
    # DeltaG_SHAPE = m ln(reactivity + 1) + b, with m=1.8 and b=-0.6.
    shape_slope = "1.8"
    shape_intercept = "-0.6"

    run_command(
        [
            fold,
            str(converted.fasta),
            str(ct),
            "--SHAPE",
            str(converted.shape),
            "--SHAPEslope",
            shape_slope,
            "--SHAPEintercept",
            shape_intercept,
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
            shape_slope,
            "--SHAPEintercept",
            shape_intercept,
            "--quiet",
        ],
        log_path,
    )
    run_command([probability_plot, str(pfs), str(probability_txt), "--text"], log_path)
    run_command([probability_plot, str(pfs), str(probability_svg), "--svg"], log_path)

    sequence = read_fasta_sequence(converted.fasta)
    pairs = parse_probability_plot_text(probability_txt)
    entropy = shannon_entropy(len(sequence), pairs)
    write_entropy(entropy_tsv, entropy, rolling_mean(entropy, 55))
    plot_reactivity_entropy(
        converted.table,
        entropy_tsv,
        profile_png,
        profile_svg,
        f"E. coli 16S rRNA cell-free {reagent} SHAPE-MaP ({scheme})",
    )
    draw_structure(ct, structure_svg, structure_png, log_path)

    outputs = {
        "fasta": str(converted.fasta.relative_to(ROOT)),
        "shape": str(converted.shape.relative_to(ROOT)),
        "converted_profile": str(converted.table.relative_to(ROOT)),
        "ct": str(ct.relative_to(ROOT)),
        "dot_bracket": str(dbn.relative_to(ROOT)),
        "partition": str(pfs.relative_to(ROOT)),
        "pairing_probability": str(probability_txt.relative_to(ROOT)),
        "entropy": str(entropy_tsv.relative_to(ROOT)),
        "reactivity_entropy_png": str(profile_png.relative_to(ROOT)),
        "reactivity_entropy_svg": str(profile_svg.relative_to(ROOT)),
        "mfe_structure_svg": str(structure_svg.relative_to(ROOT)),
    }
    if structure_png.exists():
        outputs["mfe_structure_png"] = str(structure_png.relative_to(ROOT))

    manifest.write_text(
        json.dumps(
            {
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "stage": 0,
                "case": case,
                "test": test,
                "actual_dataset": "Busan et al. 2019 E. coli 16S rRNA cell-free SHAPE-MaP",
                "reagent": reagent,
                "scheme": scheme,
                "shape_slope": float(shape_slope),
                "shape_intercept": float(shape_intercept),
                "input_profile": str(raw_profile.relative_to(ROOT)),
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"Wrote outputs to {outdir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run stage-0 RNAstructure Deigan pipeline.")
    parser.add_argument("--case", default="xist")
    parser.add_argument("--test", default="test-1")
    parser.add_argument("--reagent", default="1M7")
    parser.add_argument("--scheme", default="rnastr-deigan")
    args = parser.parse_args()
    run_stage0(args.case, args.test, args.reagent, args.scheme)


if __name__ == "__main__":
    main()
