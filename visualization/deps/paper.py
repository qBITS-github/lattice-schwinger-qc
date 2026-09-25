"""Read the paper's vector figures and map their strokes into data coordinates.

The SVG files are the figures published with arXiv:2607.02894. A solid colored
stroke in an overlay is this repository's curve, a black dashed stroke is the
paper's curve, and a dot is one of the paper's quantum-computing markers.
"""

import re
from pathlib import Path

import numpy as np

SVG_DIR = Path(__file__).resolve().parent / "svg"

# Matplotlib tab10 colors used in the paper, keyed to the quench field.
PAPER_HEX = {
    "1f77b4": 0.5,
    "ff7f0e": 1.0,
    "2ca02c": 1.5,
    "d62728": 2.0,
    "9467bd": 2.5,
    "8c564b": 3.0,
    "7f7f7f": None,
    "000000": None,
}
# Black dashes drawn on top of our solid curve. A match reads as one stroke.
PAPER_LS = (0, (3.2, 2.0))

X12 = [0, 2, 4, 6, 8, 10, 12]
X03 = [0, 1, 2, 3]


def load_svg(name):
    """Return the SVG text and the viewBox height."""
    text = (SVG_DIR / name).read_text()
    height = float(re.search(r'viewBox="0 0 [\d.]+ ([\d.]+)"', text).group(1))
    return text, height


def parse_path(d):
    """Polyline through M/L/H/V commands. Implicit pairs after M are lineto."""
    tokens = re.findall(r"[MLHVZmlhvz]|[-+]?\d*\.?\d+(?:e[-+]?\d+)?", d)
    pts = []
    cmd = None
    i = 0
    x = y = 0.0
    while i < len(tokens):
        tok = tokens[i]
        if re.fullmatch(r"[MLHVZmlhvz]", tok):
            cmd = tok
            i += 1
            continue
        if cmd in ("M", "L", "m", "l"):
            dx, dy = float(tokens[i]), float(tokens[i + 1])
            if cmd in ("m", "l") and pts:
                x, y = x + dx, y + dy
            else:
                x, y = dx, dy
            pts.append((x, y))
            if cmd == "M":
                cmd = "L"
            elif cmd == "m":
                cmd = "l"
            i += 2
        elif cmd in ("H", "h"):
            dx = float(tokens[i])
            x = x + dx if cmd == "h" else dx
            pts.append((x, y))
            i += 1
        elif cmd in ("V", "v"):
            dy = float(tokens[i])
            y = y + dy if cmd == "v" else dy
            pts.append((x, y))
            i += 1
        else:
            i += 1
    if len(pts) < 2:
        return None
    return np.array(pts, dtype=float)


def polylines(text):
    """Data polylines in y-up coordinates. Marker glyphs are skipped."""
    out = []
    for match in re.finditer(r"<path([^>]*)>", text):
        tag = match.group(1)
        dm = re.search(r'd="([^"]+)"', tag)
        sm = re.search(r'stroke="#([0-9a-fA-F]{6})"', tag)
        tm = re.search(r'transform="matrix\(1,0,0,-1,([^,]+),([^)]+)\)"', tag)
        if not (dm and sm and tm):
            continue
        tx, _ty = float(tm.group(1)), float(tm.group(2))
        if abs(tx) > 1e-6:
            continue
        pts = parse_path(dm.group(1))
        if pts is None or len(pts) < 4:
            continue
        out.append(
            {
                "hex": sm.group(1).lower(),
                "dash": "stroke-dasharray" in tag,
                "pts": pts,
            }
        )
    return out


def markers(text, height):
    """Marker centers as (hex, x, y_up)."""
    found = []
    for match in re.finditer(r"<path([^>]*)>", text):
        tag = match.group(1)
        tm = re.search(r'transform="matrix\(1,0,0,-1,([^,]+),([^)]+)\)"', tag)
        sm = re.search(r'stroke="#([0-9a-fA-F]{6})"', tag)
        if not (tm and sm):
            continue
        tx, ty = float(tm.group(1)), float(tm.group(2))
        if abs(tx) < 1e-3 and abs(ty - height) < 0.05:
            continue
        found.append((sm.group(1).lower(), tx, height - ty))
    return found


def axis_boxes(text):
    """Panel rectangles (x0, y0, x1, y1) in y-up coords, top-to-bottom then left-to-right."""
    spans = []
    for match in re.finditer(r'd="M([0-9.]+) ([0-9.]+)H([0-9.]+)"', text):
        x1, y, x2 = map(float, match.groups())
        if abs(x2 - x1) > 40:
            spans.append((round(min(x1, x2), 1), round(max(x1, x2), 1), round(y, 1)))
    groups = {}
    for x0, x1, y in spans:
        groups.setdefault((x0, x1), set()).add(y)
    boxes = []
    for (x0, x1), ys in groups.items():
        ys = sorted(ys)
        for a, b in zip(ys, ys[1:]):
            if b - a > 50:
                boxes.append((x0, a, x1, b))
    boxes.sort(key=lambda b: (-b[3], b[0]))
    return boxes


def tick_maps(text, box, x_values, y_values):
    """Map SVG pixels to data using the outward tick marks on one panel."""
    x0, y0, x1, y1 = box
    xt, yt = [], []
    for match in re.finditer(r'd="M([0-9.]+) ([0-9.]+)V([0-9.]+)"', text):
        x, a, b = map(float, match.groups())
        lo, hi = min(a, b), max(a, b)
        if 2 < hi - lo < 12 and abs(hi - y0) < 1.5 and x0 - 1 <= x <= x1 + 1:
            xt.append(x)
    for match in re.finditer(r'd="M([0-9.]+) ([0-9.]+)H([0-9.]+)"', text):
        a, y, b = map(float, match.groups())
        lo, hi = min(a, b), max(a, b)
        if 2 < hi - lo < 12 and abs(hi - x0) < 1.5 and y0 - 1 <= y <= y1 + 1:
            yt.append(y)
    xt = np.array(sorted(set(np.round(xt, 2))))
    yt = np.array(sorted(set(np.round(yt, 2))))
    if len(xt) != len(x_values) or len(yt) != len(y_values):
        raise RuntimeError(
            f"tick count x {len(xt)}!={len(x_values)} y {len(yt)}!={len(y_values)} "
            f"box={box} xt={xt} yt={yt}"
        )

    def axis_map(pix, vals):
        p0, p1 = float(pix[0]), float(pix[-1])
        v0, v1 = float(vals[0]), float(vals[-1])

        def convert(p):
            return v0 + (np.asarray(p, dtype=float) - p0) * (v1 - v0) / (p1 - p0)

        return convert

    return axis_map(xt, x_values), axis_map(yt, y_values)


def to_data(pts, fx, fy):
    """Apply the tick maps to a polyline."""
    return np.column_stack([fx(pts[:, 0]), fy(pts[:, 1])])


def in_box(pts, box, pad=8):
    """True when most of a polyline lies inside a panel."""
    x0, y0, x1, y1 = box
    x, y = pts[:, 0], pts[:, 1]
    return (
        (x > x0 - pad).mean() > 0.8
        and (x < x1 + pad).mean() > 0.8
        and (y > y0 - pad).mean() > 0.8
        and (y < y1 + pad).mean() > 0.8
    )


def paper_curves(text, boxes, tick_spec, dashed_only=True):
    """Map each polyline into the panel that contains it."""
    lines = polylines(text)
    placed = []
    for i, (box, (xv, yv)) in enumerate(zip(boxes, tick_spec)):
        fx, fy = tick_maps(text, box, xv, yv)
        for line in lines:
            if dashed_only and not line["dash"]:
                continue
            if line["pts"].shape[0] < 4:
                continue
            if not in_box(line["pts"], box):
                continue
            width = box[2] - box[0]
            if np.ptp(line["pts"][:, 0]) < 0.35 * width:
                continue
            placed.append((i, line["hex"], to_data(line["pts"], fx, fy)))
    return placed


def paper_markers(text, height, boxes, tick_spec):
    """Marker positions in data coordinates, tagged with their panel."""
    raw = markers(text, height)
    placed = []
    for i, (box, (xv, yv)) in enumerate(zip(boxes, tick_spec)):
        fx, fy = tick_maps(text, box, xv, yv)
        x0, y0, x1, y1 = box
        for hex_, x, y in raw:
            if x0 <= x <= x1 and y0 <= y <= y1:
                placed.append((i, hex_, float(fx(x)), float(fy(y))))
    return placed
