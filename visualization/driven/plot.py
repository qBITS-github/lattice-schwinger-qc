"""Render the eight figures that support the current driven-response question.

Run from any directory with no arguments. All inputs are completed experiments
under data/driven; this script never propagates a state.
Each figure has a fixed, documented question and selection. The output manifest
records its exact input hashes, including the original simulation metadata.
"""
import argparse
from hashlib import sha256
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/schwinger-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from matplotlib.colors import BoundaryNorm, LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/driven/results"
# Separate style: changing driven figures must not alter baseline renderers.
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
FLOOR = 1e-8
RUN_SUFFIX = ""

INPUTS = {}


def read_bytes(path):
    """Record every actual input, rather than identifying a plot by HEAD alone."""
    content = path.read_bytes()
    INPUTS[str(path.relative_to(ROOT))] = sha256(content).hexdigest()
    return content


def read_json(path):
    return json.loads(read_bytes(path))


def rows(run):
    directory = ROOT / "data/driven" / (run+RUN_SUFFIX)
    read_json(directory / "metadata.json")
    data = read_json(directory / "summary.json")
    # These figures compare the measured unit-coupling protocol. Fail rather
    # than silently relabel a rerun at different physical parameters.
    for row in data:
        if any(not np.isclose(row["model"][key], 1.) for key in ("a", "m", "g")):
            raise ValueError("Figure recipes require a=m=g=1; define new comparisons for other units")
        for t in row.get("trajectories", {}).values():
            # Early saved studies predate the boolean flag. Check their actual
            # diagnostics too, instead of interpreting a missing flag as a pass.
            valid = (t["max_norm_error"] < 1e-9 and t["max_charge_squared"] < 1e-12
                     and t["max_work_residual"] < max(1e-8, .01*abs(t["W"]))
                     and t["post_Y_drift"] < max(1e-8, .01*t["Y"])
                     and t["post_W_drift"] < max(1e-8, .01*abs(t["W"])))
            if t.get("uncertainty") is not None:
                valid &= all(t["uncertainty"][q] < max(1e-8, .01*abs(t[q])) for q in ("Y", "W"))
            if not valid or t.get("numerically_valid") is False:
                raise ValueError(f"Failed numerical checks in {run}, row {row['index']}")
    return data


def point(run, identifier):
    """Read a saved point only after checking its numeric-file checksum."""
    path = ROOT / "data/driven" / (run+RUN_SUFFIX) / "points" / (identifier + ".npz")
    info = read_json(path.with_suffix(".json"))
    read_bytes(path)
    if INPUTS[str(path.relative_to(ROOT))] != info["npz_sha256"]:
        raise ValueError(f"Corrupt saved point: {path}")
    with np.load(path, allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in saved.files}
    return arrays, info["result"]


def select(data, **pulse):
    return [r for r in data if all(np.isclose(r["pulse"][k], v, rtol=0, atol=1e-12)
                                  for k, v in pulse.items())]


def one(data, **pulse):
    matches = select(data, **pulse)
    if len(matches) != 1:
        raise ValueError(f"Expected one saved row for {pulse}; found {len(matches)}")
    return matches[0]


def spectrum():
    data = rows("spectrum_n8_refined_final")
    saved = [(r, point("spectrum_n8_refined_final", r["point"])[0]) for r in data]
    scale = max(float(d["weights"].max()) for _, d in saved)
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.8))
    for r, d in saved:
        ax[0].scatter(np.full(len(d["gaps"]), r["epsilon"]), d["gaps"], s=2, color="0.82", rasterized=True)
        bright = d["weights"] > FLOOR
        neutral = ax[0].scatter(np.full(bright.sum(), r["epsilon"]), d["gaps"][bright],
            s=4+22*d["weights"][bright]/scale, c=d["weights"][bright], cmap="viridis",
            norm=LogNorm(FLOOR, scale), rasterized=True)
        charge = ax[1].scatter(np.full(len(d["unrestricted_energies"]), r["epsilon"]),
            d["unrestricted_energies"]-r["E0"], c=d["unrestricted_charges"],
            cmap="coolwarm", norm=BoundaryNorm(np.arange(-4.5, 5.5), 256), s=2, rasterized=True)
    zero = next(d for r, d in saved if r["epsilon"] == 0)
    bright = zero["weights"] > FLOOR
    ax[2].vlines(zero["gaps"][bright], FLOOR, zero["weights"][bright], color=COLORS[0])
    ax[2].set_yscale("log")
    ax[2].set_ylim(FLOOR, 1.)
    ax[0].set(xlabel=r"Background $\epsilon_b$", ylabel=r"$(E_n-E_0)/g$", ylim=(0, 8),
              title="Neutral spectrum")
    ax[1].set(xlabel=r"Background $\epsilon_b$", ylabel=r"$(E-E_{0,Q=0})/g$",
              ylim=(-2, 8), title="Charge-resolved spectrum")
    ax[2].set(xlabel=r"$(E_n-E_0)/g$", ylabel=r"$|\langle n|V|0\rangle|^2/g^2$", xlim=(0, 8),
              title=r"Transitions at $\epsilon_b=0$")
    fig.colorbar(neutral, ax=ax[0], label=r"$|\langle n|V|0\rangle|^2/g^2$")
    fig.colorbar(charge, ax=ax[1], ticks=range(-4,5,2), label=r"Charge $Q$")
    return fig, {"rows": "all 79 static backgrounds", "display": "gaps up to 8g; colored transition weights >1e-8, common logarithmic scale"}


def weak_response():
    scan, checks = rows("case_b_zero_n8_final"), rows("weak_checks_n8_v2")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.8))
    slopes = {}
    for color, duration in zip(COLORS, (32., 64.)):
        group = sorted(select(scan, T_w=duration, epsilon_w=.01), key=lambda r: r["pulse"]["omega"])
        x = [r["pulse"]["omega"] for r in group]
        ax[0].plot(x, [r["Y"]/.01**2 for r in group], "o-", ms=3, color=color, label=rf"$gT_w={duration:g}$: TDSE")
        ax[0].plot(x, [r["trajectories"]["both"]["Y_LR"]/.01**2 for r in group], "--", color=color,
                   label=rf"$gT_w={duration:g}$: LR")
        ladder = sorted(select(checks, T_w=duration, omega=2.658), key=lambda r: r["pulse"]["epsilon_w"])
        a = np.array([r["pulse"]["epsilon_w"] for r in ladder])
        y = np.array([r["Y"] for r in ladder])
        slopes[str(duration)] = float(np.polyfit(np.log(a), np.log(y), 1)[0])
        ax[1].loglog(a, y, "o", color=color, label=rf"$gT_w={duration:g}$, slope {slopes[str(duration)]:.3f}")
        ax[1].loglog(a, y[0]*(a/a[0])**2, "--", color=color)
        for amplitude, style in zip((.01, .02, .04), ("-", "--", ":")):
            group = sorted(select(scan, T_w=duration, epsilon_w=amplitude), key=lambda r: r["pulse"]["omega"])
            # NaN leaves an actual gap instead of joining across unresolved data.
            departure = [r["Y"]/r["trajectories"]["both"]["Y_LR"]-1
                         if min(r["Y"], r["trajectories"]["both"]["Y_LR"]) > 1e-8 else np.nan for r in group]
            ax[2].plot([r["pulse"]["omega"] for r in group], departure, style, color=color,
                label=rf"$gT_w={duration:g},\ \epsilon_w={amplitude:g}$")
    ax[0].set(xlabel=r"Probe frequency $\omega/g$", ylabel=r"$P_{\mathrm{exc}}/\epsilon_w^2$",
              title=r"Weak probe, $\epsilon_w=0.01$")
    ax[1].set(xlabel=r"Probe amplitude $\epsilon_w$", ylabel=r"$P_{\mathrm{exc}}$", title=r"Amplitude scaling, $\omega/g=2.658$")
    ax[1].set_xticks([.0025, .005, .01], labels=[r"$2.5\times10^{-3}$", r"$5\times10^{-3}$", r"$10^{-2}$"])
    ax[1].xaxis.set_minor_locator(NullLocator())
    ax[2].axhline(0, color="k", lw=.5)
    ax[2].set(xlabel=r"Probe frequency $\omega/g$", ylabel=r"$P_{\mathrm{exc}}/P_{\mathrm{LR}}-1$",
              yscale="symlog", title="Departure from linear response")
    ax[2].set_yscale("symlog", linthresh=1e-3)
    return fig, {"resonance": "all frequencies at weak amplitude 0.01", "scaling": "all three amplitudes at omega/g=2.658, both durations",
                 "departures": "all amplitudes; both measured and predicted Y > 1e-8", "amplitude_slopes": slopes}


def case_a():
    data = rows("case_a_n8_final")
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.5))
    for i, background in enumerate((.6, 1.)):
        for j, amplitude in enumerate((.01, .02, .04)):
            group = sorted(select(data, epsilon_s=background, epsilon_w=amplitude), key=lambda r: r["pulse"]["tau"])
            x = [r["pulse"]["tau"] for r in group]
            color = COLORS[j]
            for kind, style in (("both", "-"), ("weak", "--")):
                y = np.array([r["trajectories"][kind]["Y"] for r in group])
                ax[0, i].plot(x, np.where(y > FLOOR, y, np.nan), style, color=color,
                    label=rf"$\epsilon_w={amplitude:g}$" if kind == "both" else None)
                unresolved = y <= FLOOR
                if unresolved.any():
                    ax[0, i].scatter(np.array(x)[unresolved], np.full(unresolved.sum(), FLOOR),
                                     marker="v", facecolors="none", edgecolors=color, s=18)
            if j == 0:
                strong = group[0]["trajectories"]["strong"]["Y"]
                if strong > FLOOR:
                    ax[0, i].axhline(strong, ls=":", color="0.25", label=r"$P_s$")
                else:
                    ax[0, i].axhline(FLOOR, ls=":", color="0.5", label=r"$P_s<10^{-8}$ (unresolved)")
            for target, quantity, divisor in ((ax[1, 0], "D_Y", 1), (ax[1, 1], "D_W", 8)):
                target.errorbar(x, [r[quantity]/divisor for r in group],
                    yerr=[r[quantity+"_uncertainty"]/divisor for r in group],
                    fmt="o-" if i == 0 else "s--", ms=3, color=color,
                    label=rf"$\epsilon_s={background:g},\ \epsilon_w={amplitude:g}$" if quantity == "D_W" else None)
        ax[0, i].set(yscale="log", ylabel=r"Excitation probability $P_{\mathrm{exc}}$", title=rf"$\epsilon_s={background:g}$")
    ax[1, 0].set(ylabel=r"$\Delta P_{\mathrm{exc}}$", title="Excitation nonadditivity")
    ax[1, 1].set(ylabel=r"$\Delta W/(Ng)$", title="Energy nonadditivity")
    for target in ax[1]:
        target.set_yscale("log")
        target.axhline(0, color="k", lw=.5)
    for target in ax.flat:
        target.set(xscale="log", xlabel=r"Pulse width $g\tau$")
    return fig, {"rows": "all 150 main-scan points", "ratios": "omitted; absolute excitation and signed differences are primary"}


def controls():
    data = rows("controls_n8_final")
    fig, ax = plt.subplots(2, 3, figsize=(7.2, 5.6))
    for i, background in enumerate((.6, 1.)):
        group = sorted(select(data, epsilon_s=background, epsilon_w=.01), key=lambda r: r["pulse"]["tau"])
        x = [r["pulse"]["tau"] for r in group]
        ax[i, 0].plot(x, [r["Y"] for r in group], "o-", color=COLORS[0], label="combined")
        slow = np.array([r["controls"]["same_peak_slow"]["Y"] for r in group])
        if np.all(slow <= FLOOR):
            ax[i, 0].plot(x, np.full(len(x), FLOOR), "v", mfc="none", color=COLORS[1],
                           label=r"slow: $P_{\mathrm{exc}}<10^{-8}$")
        else:
            ax[i, 0].plot(x, np.where(slow > FLOOR, slow, np.nan), "s--", color=COLORS[1], label="slow, same peak")
        ax[i, 1].plot(x, [r["controls"]["fixed_weak_budget"]["Y"] for r in group], "o-", color=COLORS[0], label=r"fixed $\int\epsilon_w(t)^2dt$")
        ax[i, 1].plot(x, [r["Y"] for r in group], "s--", color=COLORS[1], label=r"fixed $\epsilon_w=0.01$")
        r = one(group, tau=.5)
        ax[i, 2].plot([-10, 0, 10], [r["controls"]["delay_-10.0"]["Y"], r["Y"], r["controls"]["delay_10.0"]["Y"]],
                      "o-", color=COLORS[0], label="positive weak pulse")
        ax[i, 2].plot(0, r["controls"]["negative_weak"]["Y"], "x", ms=9, color=COLORS[1], label="negative weak pulse")
        for target in ax[i]:
            target.set(ylabel=rf"$P_{{\mathrm{{exc}}}}$ ($\epsilon_s={background:g}$)")
        for target in ax[i, :2]:
            target.set(xlabel=r"Pulse width $g\tau$", xscale="log", yscale="log")
            target.set_xticks([.2, .5, 1, 2], labels=["0.2", "0.5", "1", "2"])
            target.xaxis.set_minor_locator(NullLocator())
        ax[i, 2].set_xlabel(r"Delay $gt_d$ ($g\tau=0.5$)")
    for target, title in zip(ax[0], ("Equal peak field", "Equal weak-field budget", "Timing and polarity")):
        target.set_title(title)
    return fig, {"weak_amplitude": .01, "widths": "all saved at that amplitude", "delay_and_sign_width": .5,
                 "budget": "input integral, not absorbed energy"}


def response():
    data = rows("controls_n8_final")
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.5))
    population_sums = {}
    static = next(r for r in rows("spectrum_n8_refined_final") if r["epsilon"] == 0.)
    energies = point("spectrum_n8_refined_final", static["point"])[0]["energies"]
    gaps = (energies-energies[0])[1:]
    for i, background in enumerate((.6, 1.)):
        group = sorted(select(data, epsilon_s=background, epsilon_w=.01), key=lambda r: r["pulse"]["tau"])
        x = [r["pulse"]["tau"] for r in group]
        ys = np.array([r["trajectories"]["strong"]["Y"] for r in group])
        ax[0, 0].plot(x, np.array([r["Y"] for r in group])-ys, "o-", color=COLORS[i], label=rf"$\epsilon_s={background:g}$: TDSE")
        ax[0, 0].plot(x, np.array([r["Y_pump_prediction"] for r in group])-ys, "x--", color=COLORS[i], label=rf"$\epsilon_s={background:g}$: expansion")
        for key, style in (("Y_even", "o-"), ("Y_odd", "s--")):
            ax[0, 1].plot(x, [r[key] for r in group], style, color=COLORS[i], label=rf"$\epsilon_s={background:g}$: {key[2:].lower()}")
        r = one(group, tau=.5)
        # All 69 excited states enter. This decomposition sums exactly to D_Y;
        # unlike conditional populations, it attributes the signed excess.
        p = [np.array(r["trajectories"][k]["populations"])[1:] for k in ("both", "strong", "weak")]
        excess = p[0]-p[1]-p[2]
        if len(excess) != 69 or not np.isclose(excess.sum(), r["D_Y"], atol=1e-12, rtol=1e-6):
            raise ValueError("Incomplete or inconsistent state-resolved excitation")
        population_sums[str(background)] = float(excess.sum())
        target = ax[1, i]
        colors = [COLORS[0] if y >= 0 else COLORS[1] for y in excess]
        target.vlines(gaps, 0, excess, colors=colors, lw=1)
        target.scatter(gaps, excess, color=colors, s=5)
        target.axhline(0, color="k", lw=.5)
        target.axhspan(-1e-8, 1e-8, color="gray", alpha=.15)
        target.set(xlabel=r"Excitation energy $(E_n-E_0)/g$", ylabel=r"$\Delta p_n$",
                   xlim=(0, gaps.max()+.5), title=rf"$\epsilon_s={background:g},\ g\tau=0.5$")
        target.set_yscale("symlog", linthresh=1e-8)
        target.margins(y=.12)
    ax[0, 0].set(xlabel=r"Pulse width $g\tau$", xscale="log", yscale="log", ylabel=r"$\delta P=P_{s+w}-P_s$", title="Response about the strong pulse")
    ax[0, 1].set(xlabel=r"Pulse width $g\tau$", xscale="log", ylabel=r"$\delta P_{\mathrm{even}},\;\delta P_{\mathrm{odd}}$", title="Even and odd response")
    ax[0, 1].set_yscale("symlog", linthresh=1e-9)
    for target in ax[0]:
        target.set_xticks([.2, .5, 1, 2], labels=["0.2", "0.5", "1", "2"])
        target.xaxis.set_minor_locator(NullLocator())
    return fig, {"weak_amplitude": .01, "population_width": .5, "state_count": 69, "state_axis": "full excitation energy relative to the zero-field neutral ground",
                 "state_excess_sums": population_sums, "shading": "1e-8 reference floor; not a per-state confidence interval"}


def background_limits():
    data = rows("plateau_n8_final")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.8))
    for i, background in enumerate((.6, 1.)):
        selected = [one(data, epsilon_b=background, omega=2., ramp_duration=t) for t in (10., 20.)]
        ax[0].semilogy([10, 20], [1-r["ramp_up_ground_overlap"] for r in selected], "o-", color=COLORS[i], label=rf"$\epsilon_b={background:g}$")
        for color, duration in zip(COLORS, (10., 20.)):
            group = sorted(select(data, epsilon_b=background, ramp_duration=duration), key=lambda r: r["pulse"]["omega"])
            ax[i+1].errorbar([r["pulse"]["omega"] for r in group], [r["D_Y"] for r in group],
                yerr=[r["D_Y_uncertainty"] for r in group], fmt="o-", color=color, label=rf"$gT_r={duration:g}$")
        ax[i+1].axhline(0, color="k", lw=.5)
        ax[i+1].set(xlabel=r"Probe frequency $\omega/g$", ylabel=r"$\Delta P_{\mathrm{exc}}$", title=rf"Plateau $\epsilon_b={background:g}$")
    ax[0].set(xlabel=r"Ramp duration $gT_r$", ylabel=r"Preparation infidelity $1-F$", title="Background preparation")
    return fig, {"backgrounds": [.6, 1.], "ramps": [10, 20], "frequencies": "all six coarse probe frequencies",
                 "limit": "These are ramped states, not assumed stationary ground states; connecting lines guide the eye."}


def size_dependence():
    data = rows("size_case_a_final")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.8))
    central, edge = [], []
    for color, n in zip(COLORS, (8, 12, 16)):
        group = sorted([r for r in data if r["model"]["N"] == n], key=lambda r: r["pulse"]["tau"])
        x = [r["pulse"]["tau"] for r in group]
        for target, key, divisor in ((ax[0], "D_Y", 1), (ax[1], "D_W", n)):
            target.errorbar(x, [r[key]/divisor for r in group], yerr=[r[key+"_uncertainty"]/divisor for r in group],
                            fmt="o-", color=color, label=f"N={n}")
        r = one(group, tau=.5)
        v = [np.array(r["trajectories"][k]["occupation_average"]) for k in ("both", "strong", "weak")]
        occupation = v[0]-v[1]-v[2]
        central.append(float(occupation[n//2-2:n//2+2].mean()))
        edge.append(float(np.r_[occupation[:2], occupation[-2:]].mean()))
    ax[0].set(xlabel=r"Pulse width $g\tau$", ylabel=r"$\Delta P_{\mathrm{exc}}$", title="Excitation nonadditivity")
    ax[1].set(xlabel=r"Pulse width $g\tau$", ylabel=r"$\Delta W/(Ng)$", title="Energy nonadditivity")
    ax[2].plot([8, 12, 16], central, "o-", color=COLORS[0], label="central four sites")
    ax[2].plot([8, 12, 16], edge, "s--", color=COLORS[1], label="outer four sites")
    ax[2].set(xlabel=r"Number of sites $N$", ylabel=r"$\Delta\overline{o}_R$", title=r"Spatial response, $g\tau=0.5$")
    for target in ax[:2]:
        target.axhline(0, color="k", lw=.5)
    return fig, {"rows": "all nine fixed-protocol size points", "strong_peak": 1., "weak_amplitude": .01,
                 "spatial_width": .5, "spatial_definition": "vacuum-subtracted occupation(both) - occupation(strong) - occupation(weak)",
                 "limit": "fixed lattice spacing; only the strong=1 protocol was tested across sizes"}


def numerical_checks():
    run = "controls_n8_final"
    r = one(rows(run), epsilon_s=.6, epsilon_w=.01, tau=.5)
    d, _ = point(run, r["trajectories"]["both"]["point"])
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.5))
    ax[0, 0].semilogy(d["times"], np.maximum(abs(d["energy"]-d["energy"][0]-d["source_work"]), 1e-17), color=COLORS[0])
    ax[0, 0].set(xlabel=r"Time $gt$", ylabel=r"$|E(t)-E(t_i)-W_{\mathrm{src}}(t)|/g$", title="Energy-work consistency")
    ax[0, 1].semilogy(d["times"], np.maximum(abs(d["norm"]-1), 1e-17), color=COLORS[0])
    ax[0, 1].set(xlabel=r"Time $gt$", ylabel=r"$|\langle\psi|\psi\rangle-1|$",
                  title="Norm conservation")
    h, differences, sums = [], [], []
    for q in rows("independent_case_a_final"):
        adaptive = [q["trajectories"][k]["Y"] for k in ("both", "strong", "weak")]
        midpoint = [point("independent_case_a_final", q["trajectories"][k]["midpoint_point"])[1]["Y"] for k in ("both", "strong", "weak")]
        h.append(q["coordinates"]["analysis.midpoint_step"])
        differences.append(abs(midpoint[0]-midpoint[1]-midpoint[2]-q["D_Y"]))
        sums.append(float(np.abs(np.array(midpoint)-adaptive).sum()))
    ax[1, 0].loglog(h, differences, "o-", label=r"$|\Delta P_{\mathrm{mid}}-\Delta P_{\mathrm{ad}}|$")
    ax[1, 0].loglog(h, sums, "s--", label=r"$\sum_{j=s,w,s+w}|P_j^{\mathrm{mid}}-P_j^{\mathrm{ad}}|$")
    thresholds, tail_differences = [], []
    for tail in r["trajectories"]["both"]["tail_checks"]:
        threshold = tail["threshold"]
        values = []
        for kind in ("both", "strong", "weak"):
            check = next(t for t in r["trajectories"][kind]["tail_checks"] if t["threshold"] == threshold)
            values.append(point(run, check["point"])[1]["Y"])
        thresholds.append(threshold)
        tail_differences.append(abs(values[0]-values[1]-values[2]-r["D_Y"]))
    ax[1, 1].loglog(thresholds, tail_differences, "o-", label="tail-refinement difference")
    for target in ax[1]:
        target.axhline(r["D_Y"], color="k", ls=":", label=r"measured $\Delta P_{\mathrm{exc}}$")
    ax[1, 0].set(xlabel=r"Midpoint step $gh$", ylabel="Absolute probability difference", title="Independent propagation")
    ax[1, 1].set(xlabel=r"Tail threshold $\eta$", ylabel=r"$|\Delta P_{\mathrm{exc}}^{(\eta)}-\Delta P_{\mathrm{exc}}|$", title="Pulse-tail convergence")
    return fig, {"strong_peak": .6, "weak_amplitude": .01, "width": .5,
                 "D_Y": r["D_Y"], "reported_D_Y_uncertainty": r["D_Y_uncertainty"],
                 "limit": "Observed differences are numerical checks, not rigorous error bounds; the main scan retains its own refinements."}


PLOTS = (
    (spectrum, "01_accessible_states", "Which neutral excitations can the field couple to?"),
    (weak_response, "02_weak_probe", "Does the weak probe recover the accessible resonance and quadratic scaling, and where does that approximation fail?"),
    (case_a, "03_two_pulse_result", "How much excitation and energy exceed the separate-pulse comparison?"),
    (controls, "04_matched_controls", "Can peak height or input budget explain the change, and do timing and sign matter?"),
    (response, "05_response_and_states", "Does the strong-background response predict the change, and which states contribute to its signed excess?"),
    (background_limits, "06_background_preparation", "When do ramped backgrounds enhance or suppress the probe, and was their ground state prepared?"),
    (size_dependence, "07_size_dependence", "Does the selected effect persist at N=12,16, and is there a central spatial contribution?"),
    (numerical_checks, "08_numerical_checks", "Is the representative excess larger than independent-method and pulse-tail differences?"),
)


def main():
    global RUN_SUFFIX
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-suffix", default="", help="read a complete regenerated set with this suffix")
    args = parser.parse_args()
    if args.run_suffix and not re.fullmatch(r"[A-Za-z0-9_-]+", args.run_suffix):
        parser.error("Invalid run suffix")
    RUN_SUFFIX = args.run_suffix
    if any(p.is_symlink() for p in (OUT, *OUT.parents)):
        raise ValueError("Figure output must not contain symlinks")
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "legend.fontsize": 6, "legend.frameon": False, "axes.linewidth": .7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.1,
        "lines.markersize": 3.5, "errorbar.capsize": 2, "axes.grid": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
        "axes.formatter.use_mathtext": True,
    })
    source = Path(__file__).read_bytes()
    manifest = {"study": "driven", "plotter_sha256": sha256(source).hexdigest(), "run_suffix": RUN_SUFFIX,
                "notation": {"P_exc": "saved Y", "Delta_P_exc": "saved D_Y", "Delta_W": "saved D_W"},
                "resolution_floor": FLOOR,
                "captions_sha256": sha256((OUT/"captions.md").read_bytes()).hexdigest(), "plots": {}}
    for render, name, question in PLOTS:
        INPUTS.clear()
        fig, selection = render()
        for target in fig.axes:
            if target.get_legend_handles_labels()[0]:
                target.legend()
        if name == "03_two_pulse_result":
            # Six shared curve labels belong outside the data, not over the
            # small excess being compared in the lower panels.
            handles, labels = fig.axes[3].get_legend_handles_labels()
            fig.axes[3].get_legend().remove()
            fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=6)
            for target in fig.axes[:2]:
                target.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=.95)
        # Only data axes receive panel labels; colorbar axes are separate.
        data_axes = [a for a in fig.axes if a.get_label() != "<colorbar>"]
        for letter, target in zip("abcdefghijklmnopqrstuvwxyz", data_axes):
            target.set_title(f"({letter})  {target.get_title()}", loc="left", pad=7)
            target.set_title("", loc="center")
        fig.tight_layout(pad=1.1, w_pad=1.1, h_pad=1.8,
                         rect=(0, .075 if name == "03_two_pulse_result" else 0, 1, 1))
        outputs = {}
        for extension in ("png", "pdf"):
            path = OUT / f"{name}.{extension}"
            fig.savefig(path, dpi=300, bbox_inches="tight", metadata={
                "Title": "Driven Schwinger study: "+name,
                **({"CreationDate": None, "ModDate": None} if extension == "pdf" else {}),
            })
            outputs[path.name] = sha256(path.read_bytes()).hexdigest()
        plt.close(fig)
        manifest["plots"][name] = {"question": question, "selection": selection,
                                   "inputs_sha256": dict(INPUTS), "outputs_sha256": outputs}
        print(name, flush=True)
    (OUT / "plot_source.txt").write_bytes(source)
    (OUT / "plot_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n")
    print(f"Saved {len(PLOTS)} figures in {OUT}")


if __name__ == "__main__":
    main()
