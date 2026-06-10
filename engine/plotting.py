from __future__ import annotations

import argparse
import csv
import html
import struct
import zlib
from pathlib import Path


Color = tuple[int, int, int]


def _read_converted_profile(path: Path) -> tuple[list[int], list[float | None]]:
    positions: list[int] = []
    values: list[float | None] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            positions.append(int(row["position"]))
            token = row["shape_reactivity"]
            values.append(None if token == "-999" else float(token))
    return positions, values


def _read_entropy(path: Path) -> tuple[list[int], list[float], list[float]]:
    positions: list[int] = []
    raw: list[float] = []
    smoothed: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            positions.append(int(row["position"]))
            raw.append(float(row["shannon_entropy"]))
            smoothed.append(float(row["shannon_entropy_smoothed"]))
    return positions, raw, smoothed


def _scaled_points(
    positions: list[int],
    values: list[float | None],
    left: int,
    top: int,
    width: int,
    height: int,
) -> list[tuple[float, float]]:
    numeric = [value for value in values if value is not None]
    ymax = max(1.0, max(numeric) if numeric else 1.0)
    xmax = max(1, max(positions) if positions else 1)
    points: list[tuple[float, float]] = []
    for position, value in zip(positions, values):
        if value is None:
            continue
        x = left + ((position - 1) / max(1, xmax - 1)) * width
        y = top + height - (value / ymax) * height
        points.append((x, y))
    return points


def _polyline(points: list[tuple[float, float]]) -> str:
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _write_svg(
    path: Path,
    title: str,
    reactivity_positions: list[int],
    reactivities: list[float | None],
    entropy_positions: list[int],
    entropy_raw: list[float],
    entropy_smoothed: list[float],
) -> None:
    width, height = 1300, 650
    left, right = 78, 28
    panel_h = 230
    top1, top2 = 62, 360
    plot_w = width - left - right
    reactivity_points = _scaled_points(reactivity_positions, reactivities, left, top1, plot_w, panel_h)
    entropy_raw_points = _scaled_points(entropy_positions, entropy_raw, left, top2, plot_w, panel_h)
    entropy_smooth_points = _scaled_points(entropy_positions, entropy_smoothed, left, top2, plot_w, panel_h)
    escaped_title = html.escape(title)
    path.write_text(
        f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{left}" y="32" font-family="Arial, sans-serif" font-size="20" fill="#202124">{escaped_title}</text>
  <text x="18" y="{top1 + 150}" transform="rotate(-90 18 {top1 + 150})" font-family="Arial, sans-serif" font-size="13" fill="#333">SHAPE reactivity</text>
  <text x="18" y="{top2 + 150}" transform="rotate(-90 18 {top2 + 150})" font-family="Arial, sans-serif" font-size="13" fill="#333">Shannon entropy</text>
  <text x="{left + plot_w / 2 - 48:.2f}" y="632" font-family="Arial, sans-serif" font-size="13" fill="#333">16S rRNA position</text>
  <line x1="{left}" y1="{top1 + panel_h}" x2="{left + plot_w}" y2="{top1 + panel_h}" stroke="#555"/>
  <line x1="{left}" y1="{top1}" x2="{left}" y2="{top1 + panel_h}" stroke="#555"/>
  <line x1="{left}" y1="{top2 + panel_h}" x2="{left + plot_w}" y2="{top2 + panel_h}" stroke="#555"/>
  <line x1="{left}" y1="{top2}" x2="{left}" y2="{top2 + panel_h}" stroke="#555"/>
  <polyline points="{_polyline(reactivity_points)}" fill="none" stroke="#2166ac" stroke-width="1.1"/>
  <polyline points="{_polyline(entropy_raw_points)}" fill="none" stroke="#bdbdbd" stroke-width="0.7"/>
  <polyline points="{_polyline(entropy_smooth_points)}" fill="none" stroke="#b2182b" stroke-width="1.5"/>
  <text x="{left + plot_w - 185}" y="{top2 + 20}" font-family="Arial, sans-serif" font-size="12" fill="#b2182b">55-nt rolling mean</text>
</svg>
"""
    )


def _blank_rgb(width: int, height: int, color: Color = (255, 255, 255)) -> bytearray:
    data = bytearray()
    data.extend(color * width * height)
    return data


def _set_pixel(data: bytearray, width: int, height: int, x: int, y: int, color: Color) -> None:
    if 0 <= x < width and 0 <= y < height:
        offset = (y * width + x) * 3
        data[offset : offset + 3] = bytes(color)


def _draw_line(data: bytearray, width: int, height: int, p1: tuple[float, float], p2: tuple[float, float], color: Color) -> None:
    x0, y0 = round(p1[0]), round(p1[1])
    x1, y1 = round(p2[0]), round(p2[1])
    dx = abs(x1 - x0)
    dy = -abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    error = dx + dy
    while True:
        _set_pixel(data, width, height, x0, y0, color)
        if x0 == x1 and y0 == y1:
            break
        doubled = 2 * error
        if doubled >= dy:
            error += dy
            x0 += sx
        if doubled <= dx:
            error += dx
            y0 += sy


def _draw_polyline(data: bytearray, width: int, height: int, points: list[tuple[float, float]], color: Color) -> None:
    for p1, p2 in zip(points, points[1:]):
        _draw_line(data, width, height, p1, p2, color)


def _write_png(
    path: Path,
    reactivity_positions: list[int],
    reactivities: list[float | None],
    entropy_positions: list[int],
    entropy_raw: list[float],
    entropy_smoothed: list[float],
) -> None:
    width, height = 1300, 650
    left, right = 78, 28
    panel_h = 230
    top1, top2 = 62, 360
    plot_w = width - left - right
    data = _blank_rgb(width, height)

    axes = (85, 85, 85)
    for x1, y1, x2, y2 in [
        (left, top1 + panel_h, left + plot_w, top1 + panel_h),
        (left, top1, left, top1 + panel_h),
        (left, top2 + panel_h, left + plot_w, top2 + panel_h),
        (left, top2, left, top2 + panel_h),
    ]:
        _draw_line(data, width, height, (x1, y1), (x2, y2), axes)

    _draw_polyline(
        data,
        width,
        height,
        _scaled_points(reactivity_positions, reactivities, left, top1, plot_w, panel_h),
        (33, 102, 172),
    )
    _draw_polyline(
        data,
        width,
        height,
        _scaled_points(entropy_positions, entropy_raw, left, top2, plot_w, panel_h),
        (189, 189, 189),
    )
    _draw_polyline(
        data,
        width,
        height,
        _scaled_points(entropy_positions, entropy_smoothed, left, top2, plot_w, panel_h),
        (178, 24, 43),
    )

    rows = [b"\x00" + data[y * width * 3 : (y + 1) * width * 3] for y in range(height)]
    compressed = zlib.compress(b"".join(rows), level=9)

    def chunk(kind: bytes, payload: bytes) -> bytes:
        return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)

    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", compressed)
        + chunk(b"IEND", b"")
    )


def plot_reactivity_entropy(
    converted_profile: Path,
    entropy_tsv: Path,
    png: Path,
    svg: Path,
    title: str,
) -> None:
    reactivity_positions, reactivities = _read_converted_profile(converted_profile)
    entropy_positions, entropy_raw, entropy_smoothed = _read_entropy(entropy_tsv)
    png.parent.mkdir(parents=True, exist_ok=True)
    svg.parent.mkdir(parents=True, exist_ok=True)
    _write_svg(svg, title, reactivity_positions, reactivities, entropy_positions, entropy_raw, entropy_smoothed)
    _write_png(png, reactivity_positions, reactivities, entropy_positions, entropy_raw, entropy_smoothed)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot SHAPE reactivity and Shannon entropy profiles.")
    parser.add_argument("--converted-profile", required=True, type=Path)
    parser.add_argument("--entropy", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    plot_reactivity_entropy(args.converted_profile, args.entropy, args.png, args.svg, args.title)
    print(f"Wrote {args.png}")
    print(f"Wrote {args.svg}")


if __name__ == "__main__":
    main()
