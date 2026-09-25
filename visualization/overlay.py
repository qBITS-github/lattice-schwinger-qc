"""Draw this repository's figures on the paper's curves.

A solid colored line is the exact diagonalization in ``data/ED/``. Markers in
those colors are this repository's VQE, VQD, or Trotter results. A black
dashed line is the curve in the paper's figure. A gray dot is one of the
paper's quantum-computing markers. The vector figures live in
``visualization/deps/svg``.

Run from anywhere: ``python visualization/overlay.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "deps"))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from ED import (
    load_charges,
    load_condensate_spectrum,
    load_conserved,
    load_critical_scaling,
    load_decay_rate,
    load_field_energy,
    load_vacuum_fidelity,
)
from paper import (
    PAPER_HEX,
    PAPER_LS,
    X03,
    X12,
    SVG_DIR,
    axis_boxes,
    in_box,
    load_svg,
    markers,
    paper_curves,
    paper_markers,
    polylines,
    tick_maps,
    to_data,
)
from style import COLORS, PAPER_EPS_MARKER, apply_style, figure_png, save_figure

# Field colors shared with the standalone plots. Indexed by the paper's ε values.
EPS_COLOR = {
    0.5: COLORS[0],
    1.0: COLORS[1],
    1.5: COLORS[2],
    2.0: COLORS[3],
    2.5: COLORS[4],
    3.0: COLORS[5],
}


def _report(name, paper_xy, ours_x, ours_y):
    y = np.interp(paper_xy[:, 0], ours_x, ours_y)
    err = np.max(np.abs(y - paper_xy[:, 1]))
    print(f"  {name}: max |ED − paper| = {err:.4g}")


def _style_handles(with_qc=True, with_ours=False):
    handles = [Line2D([0], [0], color=COLORS[0], lw=2.2, label="ED")]
    if with_ours:
        handles.append(Line2D([0], [0], color=COLORS[0], marker="o", ms=5, lw=0, label="QC"))
    handles.append(Line2D([0], [0], color="black", lw=1.3, ls=PAPER_LS, label="paper"))
    if with_qc:
        handles.append(Line2D([0], [0], color="0.25", marker="o", ms=4.5, lw=0, label="paper QC"))
    return handles


def _marker_for(eps):
    for key, color in EPS_COLOR.items():
        if np.isclose(eps, key):
            return color, PAPER_EPS_MARKER[key]
    return None, None


def _mark_ours(ax, times, values, eps, step):
    """This repository's quantum samples, in the color of that field."""
    if values is None:
        return
    color, marker = _marker_for(eps)
    if color is None:
        return
    idx = np.arange(0, len(times), step)
    ax.plot(
        np.asarray(times)[idx],
        np.asarray(values)[idx],
        ls="none",
        marker=marker,
        ms=4.6,
        color=color,
        zorder=5,
    )


def _qc_series(measured, key, eps, panel=None):
    if measured is None:
        return None
    hit = np.where(np.isclose(measured["eps"], eps))[0]
    if hit.size == 0:
        return None
    series = measured[key][int(hit[0])]
    if panel is None or np.ndim(series) == 1:
        return series
    return series[:, panel]


def _save(fig, number):
    save_figure(fig, figure_png(number, "overlay"))


def fig1(spec, quantum=None):
    """Fig. 1: condensate and spectrum, this repository over the paper."""
    text, height = load_svg("Fig1.svg")
    boxes = axis_boxes(text)
    ticks = [
        (X03, [-0.45, -0.40, -0.35, -0.30, -0.25]),
        (X03, [-5, -2.5, 0, 2.5, 5, 7.5]),
    ]
    curves = paper_curves(text, boxes, ticks, dashed_only=False)
    marks = paper_markers(text, height, boxes, ticks)
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), constrained_layout=True)
    eps = spec["eps"]
    for panel, hex_, xy in curves:
        if panel == 0 and hex_ == "000000":
            axes[0].plot(eps, spec["gamma"], color=COLORS[0], lw=2.2, zorder=3)
            axes[0].plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.25, zorder=4)
            _report("fig1 condensate", xy, eps, spec["gamma"])
        elif panel == 1 and hex_ == "1f77b4":
            axes[1].plot(eps, spec["E0"], color=COLORS[0], lw=2.2, zorder=3)
            axes[1].plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.25, zorder=4)
            _report("fig1 E0", xy, eps, spec["E0"])
        elif panel == 1 and hex_ == "ff7f0e":
            axes[1].plot(eps, spec["E1"], color=COLORS[1], lw=2.2, zorder=3)
            axes[1].plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.25, zorder=4)
            _report("fig1 E1", xy, eps, spec["E1"])
    for panel, _hex, x, y in marks:
        axes[panel].plot(x, y, ".", color="0.45", ms=3.5, alpha=0.75, zorder=2)
    vqe = None if quantum is None else quantum.get("vqe")
    vqd = None if quantum is None else quantum.get("vqd")
    if vqe is not None:
        idx = np.unique(np.round(np.linspace(0, len(vqe["eps"]) - 1, 31)).astype(int))
        axes[0].plot(vqe["eps"][idx], vqe["gamma"][idx], ls="none", marker="o", ms=4.5, color="#d62728", zorder=5)
        axes[1].plot(vqe["eps"][idx], vqe["E0"][idx], ls="none", marker="o", ms=4.5, color=COLORS[0], zorder=5)
    if vqd is not None:
        idx = np.unique(np.round(np.linspace(0, len(vqd["eps"]) - 1, 31)).astype(int))
        axes[1].plot(vqd["eps"][idx], vqd["E1"][idx], ls="none", marker="s", ms=4.2, color=COLORS[1], zorder=5)
    for i, ax in enumerate(axes):
        fx, fy = tick_maps(text, boxes[i], *ticks[i])
        ax.set_xlim(fx(boxes[i][0]), fx(boxes[i][2]))
        ax.set_ylim(fy(boxes[i][1]), fy(boxes[i][3]))
        ax.set_xlabel(r"$\varepsilon$")
    axes[0].set_title("(a) chiral condensate")
    axes[1].set_title("(b) spectrum")
    axes[0].set_ylabel(r"$\langle\bar{\psi}\psi\rangle\cdot a$")
    axes[1].set_ylabel(r"$E\cdot a$")
    axes[1].legend(
        handles=[
            Line2D([0], [0], color=COLORS[0], lw=2, label=r"$E_0$"),
            Line2D([0], [0], color=COLORS[1], lw=2, label=r"$E_1$"),
        ],
        loc="upper left",
    )
    fig.suptitle(r"Fig. 1  over the paper, $N=8$")
    handles = _style_handles(with_ours=vqe is not None or vqd is not None)
    if vqe is not None:
        handles.append(Line2D([0], [0], color="#d62728", marker="o", ls="none", ms=5, label="VQE"))
    if vqd is not None:
        handles.append(Line2D([0], [0], color=COLORS[1], marker="s", ls="none", ms=5, label="VQD"))
    fig.legend(handles=handles, loc="outside lower center", ncol=6, frameon=False)
    _save(fig, 1)


def fig2(scaling):
    """Fig. 2: critical field versus 1/N, ED over the paper."""
    text, height = load_svg("Fig2.svg")
    boxes = axis_boxes(text)
    ticks = [([0, 0.025, 0.05, 0.075, 0.10, 0.125], [0.50, 0.55, 0.60, 0.65, 0.70])]
    fx, fy = tick_maps(text, boxes[0], *ticks[0])
    inv_n = 1.0 / scaling["N"]
    slope, intercept = np.polyfit(inv_n, scaling["eps_c"], 1)
    apply_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    for line in polylines(text):
        if line["hex"] == "808080" and np.ptp(line["pts"][:, 0]) > 50:
            xy = to_data(line["pts"], fx, fy)
            ax.plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.4, zorder=4)
    paper_pts = []
    for _hex, x, y in markers(text, height):
        if boxes[0][0] <= x <= boxes[0][2] and boxes[0][1] <= y <= boxes[0][3]:
            paper_pts.append((float(fx(x)), float(fy(y))))
    paper_pts = np.array(sorted(paper_pts))
    ax.plot(inv_n, scaling["eps_c"], "o", ms=8, color=COLORS[0], zorder=3)
    x_line = np.array([0.0, float(fx(boxes[0][2]))])
    ax.plot(x_line, slope * x_line + intercept, color=COLORS[0], lw=2.0, zorder=3)
    ax.plot(paper_pts[:, 0], paper_pts[:, 1], "s", ms=7, mfc="none", mec="black", mew=1.4, zorder=4)
    ax.set_xlim(fx(boxes[0][0]), fx(boxes[0][2]))
    ax.set_ylim(fy(boxes[0][1]), fy(boxes[0][3]))
    ax.set_xlabel(r"$1/N$")
    ax.set_ylabel(r"$\varepsilon_c(N)$")
    ax.set_title(r"Fig. 2  ED over the paper")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=COLORS[0], lw=1.6, label="ED"),
            Line2D([0], [0], marker="s", color="black", mfc="none", lw=1.4, ls=PAPER_LS, label="paper"),
        ],
        loc="upper left",
    )
    print(f"  fig2 ED intercept {intercept:.4f}")
    _save(fig, 2)


def _time_grid(text, svg, spec, ours, ylabel, number, title, panel_titles, qc=True, measured=None, qc_key=None, step=4):
    """ours[panel] is a list of (eps, times, values)."""
    boxes = axis_boxes(text)
    curves = paper_curves(text, boxes, spec)
    marks = paper_markers(text, load_svg(svg)[1], boxes, spec) if qc else []
    apply_style()
    n = len(boxes)
    ncols = 2 if n <= 4 else 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.3 * ncols, 2.55 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    flat = axes.ravel()
    seen = set()
    for panel, hex_, xy in curves:
        eps = PAPER_HEX.get(hex_)
        color = EPS_COLOR.get(eps, "0.2") if eps is not None else "0.35"
        if panel < len(ours) and eps is not None:
            for e, t, y in ours[panel]:
                if np.isclose(e, eps) and (panel, eps) not in seen:
                    flat[panel].plot(t, y, color=color, lw=2.2, zorder=3)
                    _report(f"fig{number} p{panel} eps={eps:g}", xy, t, y)
                    seen.add((panel, eps))
                    if qc_key is not None:
                        _mark_ours(
                            flat[panel], t, _qc_series(measured, qc_key, eps, panel), eps, step
                        )
        flat[panel].plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.2, zorder=4)
    for panel, _hex, x, y in marks:
        flat[panel].plot(x, y, ".", color="0.4", ms=3.2, alpha=0.7, zorder=2)
    for i, box in enumerate(boxes):
        fx, fy = tick_maps(text, box, *spec[i])
        flat[i].set_xlim(fx(box[0]), fx(box[2]))
        flat[i].set_ylim(fy(box[1]), fy(box[3]))
        flat[i].set_title(panel_titles[i])
        if i % ncols == 0:
            flat[i].set_ylabel(ylabel)
        if i >= (nrows - 1) * ncols:
            flat[i].set_xlabel(r"$t/a$")
    for j in range(n, len(flat)):
        flat[j].axis("off")
    fig.suptitle(title)
    eps_shown = sorted({e for _panel, hex_, _xy in curves if (e := PAPER_HEX.get(hex_))})
    color_handles = [
        Line2D([0], [0], color=EPS_COLOR[e], lw=2.2, label=rf"$\varepsilon={e:g}$") for e in eps_shown
    ]
    style = _style_handles(with_qc=qc, with_ours=measured is not None)
    fig.legend(
        handles=color_handles + style,
        loc="outside lower center",
        ncol=min(6, len(color_handles) + len(style)),
        frameon=False,
    )
    _save(fig, number)


def fig3(conserved, quantum=None):
    """Fig. 3: total charge and energy drift, each panel with its own label."""
    text, height = load_svg("Fig3.svg")
    boxes = axis_boxes(text)
    spec = [(X12, [-0.04, -0.02, 0, 0.02, 0.04])] * 2
    curves = paper_curves(text, boxes, spec)
    marks = paper_markers(text, height, boxes, spec)
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), sharex=True, constrained_layout=True)
    times = conserved["t_over_a"]
    eps = conserved["eps"]
    ours = [conserved["Q_N"], conserved["delta_E"]]
    measured = None if quantum is None else quantum.get("measured")
    qc_keys = ("Q_N", "delta_E")
    ylabels = [r"$Q_N(t)$", r"$\Delta E(t)\cdot a$"]
    for panel, hex_, xy in curves:
        field = PAPER_HEX.get(hex_)
        color = EPS_COLOR.get(field, "0.45") if field is not None else "0.45"
        if field is not None and np.any(np.isclose(eps, field)):
            k = int(np.where(np.isclose(eps, field))[0][0])
            axes[panel].plot(times, ours[panel][k], color=color, lw=2.2, zorder=3)
            _mark_ours(axes[panel], times, _qc_series(measured, qc_keys[panel], field), field, 6)
            _report(f"fig3 p{panel} eps={field:g}", xy, times, ours[panel][k])
        axes[panel].plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.2, zorder=4)
    for panel, _hex, x, y in marks:
        axes[panel].plot(x, y, ".", color="0.4", ms=3.2, alpha=0.7, zorder=2)
    for i, box in enumerate(boxes):
        fx, fy = tick_maps(text, box, *spec[i])
        axes[i].set_xlim(fx(box[0]), fx(box[2]))
        axes[i].set_ylim(fy(box[1]), fy(box[3]))
        axes[i].set_ylabel(ylabels[i])
        axes[i].set_xlabel(r"$t/a$")
    fig.suptitle(r"Fig. 3  over the paper, $N=8$")
    fig.legend(
        handles=[Line2D([0], [0], color=EPS_COLOR[float(e)], lw=2.2, label=rf"$\varepsilon={e:g}$") for e in eps]
        + _style_handles(with_ours=measured is not None),
        loc="outside lower center",
        ncol=6,
        frameon=False,
    )
    _save(fig, 3)


def fig6(fidelity, quantum=None):
    """Fig. 6: vacuum fidelity. Legends sit above the axes, colors then line styles."""
    text, height = load_svg("Fig6.svg")
    boxes = axis_boxes(text)
    spec = [(X12, [0, 0.2, 0.4, 0.6, 0.8, 1.0])]
    curves = paper_curves(text, boxes, spec)
    marks = paper_markers(text, height, boxes, spec)
    times = fidelity["t_over_a"]
    eps = fidelity["eps"]
    measured = None if quantum is None else quantum.get("measured")
    apply_style()
    fig, ax = plt.subplots(figsize=(9.8, 4.15))
    fig.subplots_adjust(top=0.64, bottom=0.13, left=0.07, right=0.985)
    seen = set()
    for _panel, hex_, xy in curves:
        field = PAPER_HEX.get(hex_)
        if field is None or not np.any(np.isclose(eps, field)):
            continue
        color = EPS_COLOR[field]
        k = int(np.where(np.isclose(eps, field))[0][0])
        if field not in seen:
            ax.plot(times, fidelity["P_vac"][k], color=color, lw=2.2, zorder=3)
            _mark_ours(ax, times, _qc_series(measured, "P_vac", field), field, 4)
            seen.add(field)
        ax.plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.2, zorder=4)
        _report(f"fig6 eps={field:g}", xy, times, fidelity["P_vac"][k])
    for _panel, _hex, x, y in marks:
        ax.plot(x, y, ".", color="0.25", ms=4.5, zorder=2)
    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(0.0, 1.020)
    ax.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"$t/a$")
    ax.set_ylabel(r"$P_{\mathrm{vac}}(t)$")
    eps_handles = [
        Line2D([0], [0], color=EPS_COLOR[float(e)], lw=2.4, label=rf"$\varepsilon={float(e):g}$") for e in eps
    ]
    style = [
        Line2D([0], [0], color="0.2", lw=2.2, label="ED"),
    ]
    if measured is not None:
        style.append(Line2D([0], [0], color="0.2", marker="o", ls="none", ms=5, label="QC"))
    style.extend(
        [
            Line2D([0], [0], color="black", lw=1.3, ls=PAPER_LS, label="paper"),
            Line2D([0], [0], color="0.25", marker=".", ms=6, lw=0, label="paper QC"),
        ]
    )
    color_legend = ax.legend(
        handles=eps_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.28),
        ncol=6,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.4,
        columnspacing=1.15,
    )
    ax.add_artist(color_legend)
    ax.legend(
        handles=style,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.06),
        ncol=4,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.6,
        columnspacing=1.4,
    )
    fig.suptitle(r"Fig. 6  over the paper, $N=8$", y=0.985)
    _save(fig, 6)


def fig7(decay, quantum=None):
    """Fig. 7: early-time decay rate, this repository over the paper."""
    text, _height = load_svg("Fig7.svg")
    boxes = axis_boxes(text)
    spec = [([0.5, 1.0, 1.5, 2.0, 2.5, 3.0], [0, 0.2, 0.4, 0.6, 0.8, 1.0])]
    fx, fy = tick_maps(text, boxes[0], *spec[0])
    apply_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
    for line in polylines(text):
        if not in_box(line["pts"], boxes[0]):
            continue
        if np.ptp(line["pts"][:, 0]) < 0.35 * (boxes[0][2] - boxes[0][0]):
            continue
        xy = to_data(line["pts"], fx, fy)
        if line["dash"] and line["hex"] == "1f77b4":
            ax.plot(xy[:, 0], xy[:, 1], color="black", ls=PAPER_LS, lw=1.3, zorder=4)
            _report("fig7", xy, decay["eps"], decay["gamma_eff"])
        else:
            ax.plot(xy[:, 0], xy[:, 1], "o", color="0.25", ms=5, zorder=2)
    ax.plot(decay["eps"], decay["gamma_eff"], "o-", color=COLORS[0], lw=2.2, ms=8, zorder=3)
    measured = None if quantum is None else quantum.get("measured")
    if measured is not None:
        from paperfigs import _fit_decay

        gamma_qc, _err = _fit_decay(measured["times"], measured["P_vac"], measured["a"])
        ax.plot(measured["eps"], gamma_qc, ls="none", marker="s", ms=6, color=COLORS[1], zorder=5)
    ax.set_xlim(fx(boxes[0][0]), fx(boxes[0][2]))
    ax.set_ylim(fy(boxes[0][1]), fy(boxes[0][3]))
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$\gamma_{\mathrm{eff}}\cdot a$")
    ax.set_title(r"Fig. 7  over the paper, $N=8$")
    handles = [
        Line2D([0], [0], color=COLORS[0], marker="o", lw=1.6, label="ED"),
    ]
    if measured is not None:
        handles.append(Line2D([0], [0], color=COLORS[1], marker="s", lw=0, ms=6, label="QC"))
    handles.extend(
        [
            Line2D([0], [0], color="black", lw=1.3, ls=PAPER_LS, label="paper"),
            Line2D([0], [0], color="0.4", marker=".", lw=0, ms=6, label="paper QC"),
        ]
    )
    ax.legend(
        handles=handles,
        loc="upper left",
    )
    _save(fig, 7)


def draw(numbers=None):
    """Write overlay.png for the requested paper figures."""
    numbers = list(range(1, 9) if numbers is None else numbers)
    if not (SVG_DIR / "Fig1.svg").is_file():
        raise FileNotFoundError(f"Paper figures are missing from {SVG_DIR}")
    from paperfigs import load_quantum

    quantum = load_quantum()
    measured = quantum["measured"]
    spec = load_condensate_spectrum() if 1 in numbers else None
    scaling = load_critical_scaling() if 2 in numbers else None
    conserved = load_conserved() if 3 in numbers else None
    charges = load_charges() if 4 in numbers or 8 in numbers else None
    field = load_field_energy() if 5 in numbers else None
    fidelity = load_vacuum_fidelity() if 6 in numbers else None
    decay = load_decay_rate() if 7 in numbers else None
    missing = [
        name
        for name, series in (
            ("Fig. 1", spec if 1 in numbers else True),
            ("Fig. 2", scaling if 2 in numbers else True),
            ("Fig. 3", conserved if 3 in numbers else True),
            ("Figs. 4 and 8", charges if 4 in numbers or 8 in numbers else True),
            ("Fig. 5", field if 5 in numbers else True),
            ("Fig. 6", fidelity if 6 in numbers else True),
            ("Fig. 7", decay if 7 in numbers else True),
        )
        if series is None
    ]
    if missing:
        raise FileNotFoundError("Missing exact-diagonalization data for " + ", ".join(missing))

    if 1 in numbers:
        fig1(spec, quantum)
    if 2 in numbers:
        fig2(scaling)
    if 3 in numbers:
        fig3(conserved, quantum)

    if charges is not None:
        t = charges["t_over_a"]
        eps4 = charges["eps"]
    if 4 in numbers:
        spec4 = [(X12, [-0.5, -0.25, 0, 0.25, 0.5])] * 4
        ours_q = [[(float(e), t, charges["Q_i"][k, :, i]) for k, e in enumerate(eps4)] for i in range(4)]
        _time_grid(
            load_svg("Fig4.svg")[0],
            "Fig4.svg",
            spec4,
            ours_q,
            r"$Q_i(t)$",
            4,
            r"Fig. 4  over the paper, $N=8$",
            [rf"$i={i}$" for i in range(1, 5)],
            measured=measured,
            qc_key="Q_i",
        )

    if 5 in numbers:
        spec5 = [
            (X12, [0.94, 0.96, 0.98, 1.00]),
            (X12, [3.4, 3.5, 3.6]),
            (X12, [7.00, 7.25, 7.50, 7.75, 8.00]),
            (X12, [11.5, 12.0, 12.5, 13.0, 13.5, 14.0]),
        ]
        ours_he = [[(float(e), field["t_over_a"], field["H_E"][k])] for k, e in enumerate(field["eps"])]
        _time_grid(
            load_svg("Fig5.svg")[0],
            "Fig5.svg",
            spec5,
            ours_he,
            r"$H_E(t)\cdot a$",
            5,
            r"Fig. 5  over the paper, $N=8$",
            [rf"$\varepsilon={e:g}$" for e in field["eps"]],
            measured=measured,
            qc_key="H_E",
        )
    if 6 in numbers:
        fig6(fidelity, quantum)
    if 7 in numbers:
        fig7(decay, quantum)

    if 8 in numbers:
        spec8 = [([0, 5, 10], [-0.5, 0, 0.5])] * 8
        ours_n = [[(float(e), t, charges["q_n"][k, :, i]) for k, e in enumerate(eps4)] for i in range(8)]
        _time_grid(
            load_svg("Fig8.svg")[0],
            "Fig8.svg",
            spec8,
            ours_n,
            r"$q_n(t)$",
            8,
            r"Fig. 8  over the paper, $N=8$",
            [rf"$n={n}$" for n in range(1, 9)],
            measured=measured,
            qc_key="q_n",
        )


def main():
    """Write overlay.png for Figs. 1–8."""
    draw()


if __name__ == "__main__":
    main()
