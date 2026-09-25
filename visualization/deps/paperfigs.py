"""Draw Figs. 1–8 in the paper's layout from this repository's datasets.

A dashed curve is exact diagonalization. A marker is VQE, VQD, or
second-order Trotter, using the paper's colors and shapes. The files are
``figures/fig*/reproduction.png``.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

VIS = Path(__file__).resolve().parents[1]
if str(VIS) not in sys.path:
    sys.path.insert(0, str(VIS))

from style import (
    PAPER_DASH,
    PAPER_EPS_COLOR,
    PAPER_EPS_MARKER,
    drop_legacy_pngs,
    paper_style,
    reproduction_png,
    save_figure,
)

_QUANTUM = None

# Axis frames read off the published figures.
_XLIM_TIME = (-0.60, 12.60)
_FIG1_XLIM = (-0.05, 3.05)
_FIG4_YLIM = (-0.606, 0.606)
_FIG8_YLIM = (-0.628, 0.628)
_HE_WINDOWS = (
    ((-0.60, 12.60), (0.930, 1.004), (0.94, 0.96, 0.98, 1.00)),
    ((-0.60, 12.60), (3.326, 3.659), (3.4, 3.5, 3.6)),
    ((-0.60, 12.60), (6.840, 8.094), (7.00, 7.25, 7.50, 7.75, 8.00)),
    ((-0.60, 12.60), (11.395, 14.313), (11.5, 12.0, 12.5, 13.0, 13.5, 14.0)),
)


def load_quantum():
    """VQE, VQD, and Trotter series, or None for any dataset that is absent."""
    global _QUANTUM
    if _QUANTUM is not None:
        return _QUANTUM
    from Trotter import load_quench, measure
    from VQD import load_excited_scan
    from VQE import load_ground_states

    quench = load_quench()
    measured = None
    if quench is not None:
        print("measuring Trotter observables")
        measured = measure(quench)
    _QUANTUM = {
        "vqe": load_ground_states(),
        "vqd": load_excited_scan(),
        "measured": measured,
    }
    return _QUANTUM


def _style(eps):
    for key, color in PAPER_EPS_COLOR.items():
        if np.isclose(eps, key):
            return color, PAPER_EPS_MARKER[key]
    raise KeyError(eps)


def _row(eps_values, eps):
    if eps_values is None:
        return None
    hit = np.where(np.isclose(eps_values, eps))[0]
    return None if hit.size == 0 else int(hit[0])


def _marks(n, step):
    return np.arange(0, n, step)


def _field_marks(n):
    """About 31 markers across a field scan, matching Fig. 1."""
    return np.unique(np.round(np.linspace(0, n - 1, 31)).astype(int))


def _crossings(eps, gamma, jump=0.05):
    """External fields where the condensate jumps between plateaus."""
    groups = []
    for i in np.where(np.diff(gamma) > jump)[0]:
        if groups and i <= groups[-1][-1] + 1:
            groups[-1].append(int(i))
        else:
            groups.append([int(i)])
    return [0.5 * (eps[g[0]] + eps[g[-1] + 1]) for g in groups]


def _legend(ax, handles, ncol, y=1.02):
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.4,
        columnspacing=1.1,
    )


def _quench_handles(eps_values):
    handles = [
        Line2D([0], [0], color="black", ls=PAPER_DASH, lw=1.4, label="Exact"),
        Line2D([0], [0], color="black", marker="o", ls="none", ms=5, label="QC"),
    ]
    for eps in eps_values:
        color, marker = _style(eps)
        handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker=marker,
                ls=PAPER_DASH,
                ms=5,
                label=rf"$\varepsilon={eps:g}$",
            )
        )
    return handles


def _exact_qc(ax, times, exact, qc, eps, step):
    color, marker = _style(eps)
    if exact is not None:
        ax.plot(times, exact, color=color, ls=PAPER_DASH, lw=1.35, zorder=2)
    if qc is not None:
        idx = _marks(len(times), step)
        ax.plot(
            times[idx],
            qc[idx],
            ls="none",
            marker=marker,
            ms=5,
            color=color,
            markeredgewidth=0.45,
            zorder=3,
        )


def _panel_letter(ax, letter):
    ax.text(
        0.04,
        0.95,
        f"({letter})",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontweight="bold",
    )


def _finish(fig, number):
    drop_legacy_pngs(number)
    save_figure(fig, reproduction_png(number))


def _fit_decay(times, p_vac, a, t_min=0.0, t_max=1.0):
    mask = (times / a >= t_min) & (times / a <= t_max)
    t_fit = times[mask]
    gamma = np.empty(p_vac.shape[0])
    err = np.empty(p_vac.shape[0])
    for k in range(p_vac.shape[0]):
        window = p_vac[k, mask]
        coeff, cov = np.polyfit(t_fit, -np.log(window), 1, cov=True)
        gamma[k] = float(coeff[0]) * a
        err[k] = float(np.sqrt(cov[0, 0]) * a)
    return gamma, err


def reproduce_fig1(ed, quantum):
    """Fig. 1: condensate and spectrum, with VQE and VQD markers when saved."""
    if ed is None:
        print("skip Fig. 1 reproduction: exact condensate is missing")
        return
    vqe = quantum["vqe"]
    vqd = quantum["vqd"]
    if vqe is None:
        print("skip Fig. 1 VQE markers: data/VQE/ground_states_n8.npz is missing")
    if vqd is None:
        print("skip Fig. 1 VQD markers: data/VQD/e1_scan_n8.npz is missing")
    paper_style()
    eps = ed["eps"]
    crossings = _crossings(eps, ed["gamma"])
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.15))
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.08, right=0.98, wspace=0.28)
    axes[0].plot(eps, ed["gamma"], color="black", ls=PAPER_DASH, lw=1.5, zorder=2)
    if vqe is not None:
        idx = _field_marks(len(vqe["eps"]))
        axes[0].plot(
            vqe["eps"][idx],
            vqe["gamma"][idx],
            ls="none",
            marker="o",
            ms=5,
            color="#ff0000",
            zorder=3,
        )
    axes[1].plot(eps, ed["E0"], color="#1f77b4", ls=PAPER_DASH, lw=1.5, zorder=2)
    axes[1].plot(eps, ed["E1"], color="#ff7f0e", ls=PAPER_DASH, lw=1.5, zorder=2)
    if vqe is not None:
        idx = _field_marks(len(vqe["eps"]))
        axes[1].plot(
            vqe["eps"][idx],
            vqe["E0"][idx],
            ls="none",
            marker="o",
            ms=5,
            color="#1f77b4",
            zorder=3,
        )
    if vqd is not None:
        idx = _field_marks(len(vqd["eps"]))
        axes[1].plot(
            vqd["eps"][idx],
            vqd["E1"][idx],
            ls="none",
            marker="s",
            ms=4.5,
            color="#ff7f0e",
            zorder=3,
        )
    labels = (r"$\varepsilon_1$", r"$\varepsilon_2$")
    # (dx, y) in each panel's data coordinates, just beside the dashed line.
    notes = (
        ((0.08, -0.40), (-0.62, -0.27)),
        ((0.08, -2.55), (-0.78, 0.35)),
    )
    for ax, text_xy in zip(axes, notes):
        for xpos, (dx, y), name in zip(crossings, text_xy, labels):
            ax.axvline(xpos, color="0.45", ls=PAPER_DASH, lw=1.0, zorder=1)
            ax.text(xpos + dx, y, rf"{name}$\approx{xpos:.2f}$", fontsize=10, color="0.15")
        ax.set_xlim(*_FIG1_XLIM)
        ax.set_xlabel(r"$\varepsilon$")
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    axes[0].set_ylim(-0.458, -0.211)
    axes[0].yaxis.set_major_locator(plt.MultipleLocator(0.05))
    axes[0].set_ylabel(r"$\langle\bar{\psi}\psi\rangle\cdot a$")
    axes[0].text(0.04, 0.95, "(a)", transform=axes[0].transAxes, va="top", fontweight="bold")
    axes[1].set_ylim(-5.29, 9.05)
    axes[1].yaxis.set_major_locator(plt.MultipleLocator(2.5))
    axes[1].set_ylabel(r"$E\cdot a$")
    axes[1].text(0.04, 0.95, "(b)", transform=axes[1].transAxes, va="top", fontweight="bold")
    _legend(
        axes[0],
        [
            Line2D([0], [0], color="black", ls=PAPER_DASH, lw=1.5, label="ED"),
            Line2D([0], [0], color="#ff0000", marker="o", ls="none", ms=5, label="VQE"),
        ],
        ncol=2,
    )
    _legend(
        axes[1],
        [
            Line2D([0], [0], color="#1f77b4", ls=PAPER_DASH, lw=1.5, label=r"ED $E_0$"),
            Line2D([0], [0], color="#1f77b4", marker="o", ls="none", ms=5, label=r"VQE $E_0$"),
            Line2D([0], [0], color="#ff7f0e", ls=PAPER_DASH, lw=1.5, label=r"ED $E_1$"),
            Line2D([0], [0], color="#ff7f0e", marker="s", ls="none", ms=4.5, label=r"VQD $E_1$"),
        ],
        ncol=2,
        y=1.02,
    )
    _finish(fig, 1)


def reproduce_fig2(scaling):
    """Fig. 2: critical field versus 1/N and the linear extrapolation."""
    if scaling is None:
        print("skip Fig. 2 reproduction: exact critical-field data is missing")
        return
    paper_style()
    inv_n = 1.0 / scaling["N"]
    slope, intercept = np.polyfit(inv_n, scaling["eps_c"], 1)
    print(f"eps_c(N->inf) = {intercept:.6f}  (slope d(eps_c)/d(1/N) = {slope:.6f})")
    fig, ax = plt.subplots(figsize=(5.0, 3.9))
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.96)
    x_line = np.array([0.0, 0.135])
    ax.plot(x_line, slope * x_line + intercept, color="0.45", ls=PAPER_DASH, lw=1.6, zorder=2)
    ax.plot(inv_n, scaling["eps_c"], "o", color="black", ms=6.5, zorder=3)
    ax.text(0.008, intercept + 0.008, rf"$\varepsilon_c(\infty)\approx{intercept:.3f}$", fontsize=11)
    ax.set_xlim(-0.004, 0.135)
    ax.set_ylim(0.46, 0.71)
    ax.set_xticks([0.0, 0.025, 0.05, 0.075, 0.10, 0.125])
    ax.set_yticks([0.50, 0.55, 0.60, 0.65, 0.70])
    ax.set_xlabel(r"$1/N$")
    ax.set_ylabel(r"$\varepsilon_c(N)$")
    ax.legend(
        handles=[
            Line2D([0], [0], color="black", marker="o", ls="none", ms=6, label="ED data"),
            Line2D([0], [0], color="0.45", ls=PAPER_DASH, lw=1.6, label="Linear fit"),
        ],
        loc="upper left",
        frameon=False,
    )
    _finish(fig, 2)


def reproduce_fig3(conserved, quantum):
    """Fig. 3: total charge and energy drift. ε = 0 is the stationary vacuum."""
    if conserved is None:
        print("skip Fig. 3 reproduction: exact quench data is missing")
        return
    measured = quantum["measured"]
    paper_style()
    times = conserved["t_over_a"]
    fields = (0.0,) + tuple(float(e) for e in conserved["eps"])
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.08, right=0.98, wspace=0.12)
    series = (conserved["Q_N"], conserved["delta_E"])
    qc_keys = ("Q_N", "delta_E")
    for panel, (exact, key) in enumerate(zip(series, qc_keys)):
        _exact_qc(axes[panel], times, np.zeros_like(times), None, 0.0, 6)
        for k, eps in enumerate(conserved["eps"]):
            qc = None
            if measured is not None:
                row = _row(measured["eps"], eps)
                qc = None if row is None else measured[key][row]
            _exact_qc(axes[panel], times, exact[k], qc, float(eps), 6)
        _panel_letter(axes[panel], "ab"[panel])
        axes[panel].set_xlabel(r"$t/a$")
        axes[panel].set_xlim(0.0, 12.0)
        axes[panel].set_ylim(-0.05, 0.05)
        axes[panel].set_xticks([0, 2, 4, 6, 8, 10, 12])
        axes[panel].set_yticks([-0.04, -0.02, 0.0, 0.02, 0.04])
    axes[0].set_ylabel(r"$Q_N(t)$")
    axes[1].set_ylabel(r"$\Delta E(t)\cdot a$")
    fig.legend(
        handles=_quench_handles(fields),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=7,
        frameon=False,
        borderaxespad=0.0,
    )
    _finish(fig, 3)


def _charge_figure(times, eps, exact, qc, ylabel, titles, number, ylim, yticks, xticks, xlim):
    paper_style()
    n = exact.shape[-1]
    ncols = 2 if n <= 4 else 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(2.35 * ncols + 0.4, 2.15 * nrows + 0.7),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    fig.subplots_adjust(top=0.86, bottom=0.08, left=0.07, right=0.985, hspace=0.38, wspace=0.12)
    letters = "abcdefgh"
    for i, ax in enumerate(axes.ravel()):
        for k, eps_k in enumerate(eps):
            qc_row = None if qc is None else qc[k, :, i]
            _exact_qc(ax, times, exact[k, :, i], qc_row, float(eps_k), 4)
        ax.set_title(titles[i])
        _panel_letter(ax, letters[i])
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        if i % ncols == 0:
            ax.set_ylabel(ylabel)
        if i >= (nrows - 1) * ncols:
            ax.set_xlabel(r"$t/a$")
    fig.legend(
        handles=_quench_handles(tuple(float(e) for e in eps)),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.93),
        ncol=6,
        frameon=False,
        borderaxespad=0.0,
    )
    _finish(fig, number)


def reproduce_fig4(charges, quantum):
    """Fig. 4: spatial-point charge, one panel per spatial point."""
    if charges is None:
        print("skip Fig. 4 reproduction: exact charge data is missing")
        return
    qc = _qc_charge(charges, quantum, "Q_i")
    _charge_figure(
        charges["t_over_a"],
        charges["eps"],
        charges["Q_i"],
        qc,
        r"$Q_i(t)$",
        [rf"$i={i}$" for i in range(1, 5)],
        4,
        _FIG4_YLIM,
        [-0.5, -0.25, 0.0, 0.25, 0.5],
        [0, 2, 4, 6, 8, 10, 12],
        _XLIM_TIME,
    )


def reproduce_fig8(charges, quantum):
    """Fig. 8: site charge, one panel per lattice site."""
    if charges is None:
        print("skip Fig. 8 reproduction: exact charge data is missing")
        return
    qc = _qc_charge(charges, quantum, "q_n")
    _charge_figure(
        charges["t_over_a"],
        charges["eps"],
        charges["q_n"],
        qc,
        r"$q_n(t)$",
        [rf"$n={n}$" for n in range(1, 9)],
        8,
        _FIG8_YLIM,
        [-0.5, 0.0, 0.5],
        [0, 5, 10],
        _XLIM_TIME,
    )


def _qc_charge(charges, quantum, key):
    measured = quantum["measured"]
    if measured is None:
        return None
    rows = [_row(measured["eps"], eps) for eps in charges["eps"]]
    if any(row is None for row in rows):
        return None
    return np.stack([measured[key][row] for row in rows], axis=0)


def reproduce_fig5(field, quantum):
    """Fig. 5: electric-field energy, one panel per quench field."""
    if field is None:
        print("skip Fig. 5 reproduction: exact field-energy data is missing")
        return
    measured = quantum["measured"]
    paper_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.3))
    fig.subplots_adjust(top=0.80, bottom=0.07, left=0.08, right=0.98, hspace=0.42, wspace=0.22)
    letters = "abcd"
    for k, ax in enumerate(axes.ravel()):
        eps = float(field["eps"][k])
        qc = None
        if measured is not None:
            row = _row(measured["eps"], eps)
            qc = None if row is None else measured["H_E"][row]
        _exact_qc(ax, field["t_over_a"], field["H_E"][k], qc, eps, 4)
        xlim, ylim, yticks = _HE_WINDOWS[k]
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
        ax.set_yticks(yticks)
        ax.set_title(rf"$\varepsilon={eps:g}$")
        ax.set_xlabel(r"$t/a$")
        ax.set_ylabel(r"$H_E(t)\cdot a$")
        ax.annotate(
            f"({letters[k]})",
            xy=(0, 1),
            xytext=(-36, 3),
            textcoords="offset points",
            xycoords="axes fraction",
            fontweight="bold",
            annotation_clip=False,
        )
    fig.legend(
        handles=_quench_handles(tuple(float(e) for e in field["eps"])),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.93),
        ncol=6,
        frameon=False,
        borderaxespad=0.0,
    )
    _finish(fig, 5)


def reproduce_fig6(fidelity, quantum):
    """Fig. 6: vacuum fidelity, one curve per quench field."""
    if fidelity is None:
        print("skip Fig. 6 reproduction: exact fidelity data is missing")
        return
    measured = quantum["measured"]
    paper_style()
    fig, ax = plt.subplots(figsize=(9.6, 4.05))
    fig.subplots_adjust(top=0.72, bottom=0.14, left=0.08, right=0.985)
    for k, eps in enumerate(fidelity["eps"]):
        qc = None
        if measured is not None:
            row = _row(measured["eps"], eps)
            qc = None if row is None else measured["P_vac"][row]
        _exact_qc(ax, fidelity["t_over_a"], fidelity["P_vac"][k], qc, float(eps), 4)
    ax.set_xlim(0.0, 12.0)
    ax.set_ylim(0.0, 1.02)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"$t/a$")
    ax.set_ylabel(r"$P_{\mathrm{vac}}(t)$")
    fields = [0.5, 1.5, 2.5, 1.0, 2.0, 3.0]
    handles = [
        Line2D([0], [0], color="black", ls=PAPER_DASH, lw=1.4, label="Exact"),
        Line2D([0], [0], color="black", marker="o", ls="none", ms=5, label="QC"),
    ]
    for eps in fields:
        color, marker = _style(eps)
        handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                marker=marker,
                ls=PAPER_DASH,
                ms=5,
                label=rf"$\varepsilon={eps:g}$",
            )
        )
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=5,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.6,
        columnspacing=1.15,
    )
    _finish(fig, 6)


def reproduce_fig7(decay, fidelity, quantum):
    """Fig. 7: early-time decay rate, exact and Trotter, with fit uncertainties."""
    if decay is None or fidelity is None:
        print("skip Fig. 7 reproduction: exact decay data is missing")
        return
    a = fidelity["a"]
    _gamma, err_ed = _fit_decay(fidelity["t_over_a"] * a, fidelity["P_vac"], a, decay["t_min"], decay["t_max"])
    # The saved rate is the plotted value. The error bar is the fit uncertainty.
    gamma_ed = decay["gamma_eff"]
    paper_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.96)
    ax.errorbar(
        decay["eps"],
        gamma_ed,
        yerr=err_ed,
        color="#1f77b4",
        marker="o",
        ms=6,
        ls=PAPER_DASH,
        lw=1.3,
        capsize=3,
        zorder=3,
    )
    measured = quantum["measured"]
    if measured is not None:
        gamma_qc, err_qc = _fit_decay(measured["times"], measured["P_vac"], measured["a"])
        ax.errorbar(
            measured["eps"],
            gamma_qc,
            yerr=err_qc,
            color="#ff7f0e",
            marker="s",
            ms=5.5,
            ls="-",
            lw=1.3,
            capsize=3,
            zorder=4,
        )
    else:
        print("skip Fig. 7 QC markers: data/Trotter/quench_n8.npz is missing")
    ax.set_xlim(0.374, 3.125)
    ax.set_ylim(-0.017, 1.007)
    ax.set_xticks([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$\gamma_{\mathrm{eff}}\cdot a$")
    ax.legend(
        handles=[
            Line2D([0], [0], color="#1f77b4", marker="o", ls=PAPER_DASH, ms=6, label="Exact"),
            Line2D([0], [0], color="#ff7f0e", marker="s", ls="-", ms=5.5, label="QC"),
        ],
        loc="upper left",
        frameon=False,
    )
    _finish(fig, 7)


def write(numbers=None):
    """Write ``reproduction.png`` for the requested paper figures."""
    from ED import (
        load_charges,
        load_condensate_spectrum,
        load_conserved,
        load_critical_scaling,
        load_decay_rate,
        load_field_energy,
        load_vacuum_fidelity,
    )

    numbers = list(range(1, 9) if numbers is None else numbers)
    quantum = load_quantum()
    if 1 in numbers:
        reproduce_fig1(load_condensate_spectrum(), quantum)
    if 2 in numbers:
        reproduce_fig2(load_critical_scaling())
    if 3 in numbers:
        reproduce_fig3(load_conserved(), quantum)
    charges = load_charges() if (4 in numbers or 8 in numbers) else None
    if 4 in numbers:
        reproduce_fig4(charges, quantum)
    if 5 in numbers:
        reproduce_fig5(load_field_energy(), quantum)
    fidelity = load_vacuum_fidelity() if (6 in numbers or 7 in numbers) else None
    if 6 in numbers:
        reproduce_fig6(fidelity, quantum)
    if 7 in numbers:
        reproduce_fig7(load_decay_rate(), fidelity, quantum)
    if 8 in numbers:
        reproduce_fig8(charges, quantum)
