from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_fasta_sequence(path: Path) -> str:
    sequence: list[str] = []
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith(">"):
                sequence.append(text.upper())
    return "".join(sequence)


def parse_probability_plot_text(path: Path) -> dict[int, list[tuple[int, float]]]:
    pairs: dict[int, list[tuple[int, float]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith(("#", ";")):
                continue
            fields = text.split()
            if len(fields) < 3:
                continue
            try:
                i = int(fields[0])
                j = int(fields[1])
                minus_log10_probability = float(fields[2])
            except ValueError:
                continue
            probability = 10 ** (-minus_log10_probability)
            if probability <= 0 or not math.isfinite(probability):
                continue
            pairs[i].append((j, probability))
            pairs[j].append((i, probability))
    return pairs


def parse_vienna_dp_ps(path: Path) -> dict[int, list[tuple[int, float]]]:
    pairs: dict[int, list[tuple[int, float]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 4 or fields[3] != "ubox":
                continue
            try:
                i = int(fields[0])
                j = int(fields[1])
                sqrt_probability = float(fields[2])
            except ValueError:
                continue
            probability = sqrt_probability * sqrt_probability
            if probability <= 0 or not math.isfinite(probability):
                continue
            pairs[i].append((j, probability))
            pairs[j].append((i, probability))
    return pairs


def dotbracket_pairs(dotbracket: str) -> set[tuple[int, int]]:
    stacks: dict[str, list[int]] = {"(": [], "[": [], "{": [], "<": []}
    closers = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs: set[tuple[int, int]] = set()
    for index, char in enumerate(dotbracket.strip(), start=1):
        if char in stacks:
            stacks[char].append(index)
        elif char in closers:
            opener = closers[char]
            if stacks[opener]:
                left = stacks[opener].pop()
                pairs.add((left, index))
    return pairs


def read_dbn_structure(path: Path) -> str:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith(">")]
    if len(lines) < 2:
        raise ValueError(f"Could not read dot-bracket structure from {path}")
    return lines[1]


def ct_pairs(path: Path) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    with path.open() as handle:
        next(handle, None)
        for line in handle:
            fields = line.split()
            if len(fields) < 5:
                continue
            i = int(fields[0])
            j = int(fields[4])
            if j > i:
                pairs.add((i, j))
    return pairs


def compare_pair_sets(a: set[tuple[int, int]], b: set[tuple[int, int]]) -> dict[str, float | int]:
    shared = len(a & b)
    union = len(a | b)
    return {
        "pairs_a": len(a),
        "pairs_b": len(b),
        "shared_pairs": shared,
        "jaccard": shared / union if union else 1.0,
        "a_overlap": shared / len(a) if a else 1.0,
        "b_overlap": shared / len(b) if b else 1.0,
        "base_pair_distance": len(a ^ b),
    }


def shannon_entropy(length: int, pairs: dict[int, list[tuple[int, float]]]) -> list[float]:
    entropies: list[float] = []
    for position in range(1, length + 1):
        probabilities = [probability for _, probability in pairs.get(position, []) if probability > 0]
        paired_total = sum(probabilities)
        if paired_total > 1:
            probabilities = [probability / paired_total for probability in probabilities]
            paired_total = 1.0
        unpaired = max(0.0, 1.0 - paired_total)
        terms = probabilities + ([unpaired] if unpaired > 0 else [])
        entropies.append(-sum(probability * math.log(probability) for probability in terms if probability > 0))
    return entropies


def rolling_mean(values: list[float], window: int) -> list[float]:
    if window <= 1:
        return values[:]
    radius = window // 2
    smoothed: list[float] = []
    for index in range(len(values)):
        start = max(0, index - radius)
        end = min(len(values), index + radius + 1)
        smoothed.append(sum(values[start:end]) / (end - start))
    return smoothed


def write_entropy(path: Path, entropies: list[float], smoothed: list[float]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["position", "shannon_entropy", "shannon_entropy_smoothed"])
        for position, (entropy, smooth) in enumerate(zip(entropies, smoothed), start=1):
            writer.writerow([position, f"{entropy:.6f}", f"{smooth:.6f}"])


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute per-position Shannon entropy from RNAstructure pairing probabilities.")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--probability-text", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--window", type=int, default=55)
    args = parser.parse_args()

    sequence = read_fasta_sequence(args.fasta)
    pairs = parse_probability_plot_text(args.probability_text)
    entropies = shannon_entropy(len(sequence), pairs)
    write_entropy(args.output, entropies, rolling_mean(entropies, args.window))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
