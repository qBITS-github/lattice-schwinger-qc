"""Run the retained experiments: prepare -> propagate -> observe -> save.

scripts/driven/run.py supplies a recipe and its exact scan-row selection.
This module expands that grid, runs matched trajectories, and records the
response predictions and numerical comparisons used by the eight figures.
"""
from copy import deepcopy
from dataclasses import replace
from itertools import product
import json

import numpy as np

from .io import StudyStore
from .model import build_sector_operators, sector_ground_state, transition_spectrum, lanczos_spectrum
from .observables import make_observer, summarize
from .pulses import Constant, Sauter, Sine, Sum, pulse_components
from .response import linear_response, pump_expansion, squared_drive
from .evolution import propagate, sampling_grid, step_cap


def scan_points(config):
    """Dot-separated keys select explicit nested JSON fields, with stable indices."""
    scan = config.get("scan", {})
    for index, values in enumerate(product(*scan.values())):
        point = deepcopy(config)
        coordinates = dict(zip(scan, values))
        for key, value in coordinates.items():
            target = point
            parts = key.split(".")
            for part in parts[:-1]: target = target[part]
            target[parts[-1]] = value
        yield index, point, coordinates


def ground(store, ops, model, eps):
    specification = {"kind": "ground_state", "model": model, "epsilon": eps,
                     "method": "dense_eigh_below_1025_else_eigsh", "tol": 1e-12}
    def compute():
        energy, state = sector_ground_state(ops, eps)
        residual = np.linalg.norm(ops.hamiltonian(eps)@state-energy*state)
        return {"energy": np.array(energy), "state": state}, {"residual": float(residual)}
    arrays, info, identifier = store.point(specification, compute)
    return float(arrays["energy"]), arrays["state"], identifier


def run_trajectory(store, ops, config, vacuum, pulse, interval, settings=None):
    """Shared source-work and spatial readout; state history is optional.

    Only endpoint/final states are mandatory. Large-N runs stream observables
    through the adaptive observer instead of retaining the state trajectory.
    """
    settings = dict(config["integrator"] if settings is None else settings)
    method = settings.pop("method", "DOP853")
    settings["max_step"] = settings.get("max_step", (.05 if method=="DOP853" else .02)/ops.g)
    readout = config["readout"]
    start, end = interval
    post = readout.get("post_pulse_duration", 10/ops.g)
    events = [*pulse.boundaries, end]
    times = sampling_grid(start, end+post, readout.get("sampling_interval", .05/ops.g), events)
    specification = {"kind": "trajectory", "model": config["model"], "pulse": pulse.metadata(),
                     "interval": interval, "readout": readout, "initial_state": config["initial_state"],
                     "background": 0., "method": method, "settings": settings,
                     "seed": config.get("seed"), "store_states": config.get("analysis", {}).get("store_states", False)}
    def compute():
        observer = make_observer(ops, vacuum, pulse)
        endpoint_state = []
        ramp_states = []
        def observe(t, psi):
            if abs(t-end) < 1e-10: endpoint_state.append(psi.copy())
            if config["pulse"]["kind"] == "plateau_sine" and abs(t-config["pulse"]["ramp_duration"]) < 1e-10:
                ramp_states.append(psi.copy())
            return observer(t, psi)
        Eref = float(np.vdot(vacuum, ops.H0@vacuum).real)
        trajectory = propagate(ops, pulse, vacuum, times, method=method, energy_shift=Eref,
                               observer=observe, store_states=specification["store_states"], **settings)
        arrays = {"times": times, "integration_times": trajectory.integration_times,
                  "source_work": trajectory.work, "endpoint_state": endpoint_state[0],
                  "final_state": trajectory.final_state, **trajectory.observations}
        if trajectory.states is not None: arrays["states"] = trajectory.states
        if ramp_states: arrays["ramp_state"] = ramp_states[0]
        info = summarize(trajectory, end, ops.N, post)
        info["integrator"] = trajectory.settings
        info["tail_treatment"] = pulse.metadata()
        info["preparation"] = {"method": config["initial_state"], "Q": ops.Q, "epsilon": 0.}
        return arrays, info
    return store.point(specification, compute)


def run_spectrum_point(store, cfg):
    ops = build_sector_operators(**cfg["model"])
    background = cfg.get("analysis", {}).get("background", 0.)
    E0, vacuum, _ = ground(store, ops, cfg["model"], background)
    steps = cfg.get("analysis", {}).get("lanczos_steps", 128)
    specification = {"kind": "spectrum", "model": cfg["model"], "background": background,
                     "lanczos_steps": steps, "unrestricted": cfg["model"]["N"] <= 8}
    def compute():
        spectrum = (transition_spectrum(ops, background) if len(ops.basis) <= 1024 else
                    lanczos_spectrum(ops, background, steps=steps, ground=(E0, vacuum)))
        arrays = {"gaps": spectrum["gaps"], "weights": spectrum["weights"],
                  "basis": ops.basis, "ground_state": vacuum, "ground_energy": np.array(E0)}
        if "energies" in spectrum: arrays["energies"] = spectrum["energies"]
        if ops.N <= 8:
            energies, charges = [], []
            for Q in range(-ops.N//2, ops.N//2+1):
                sector = build_sector_operators(**{**cfg["model"], "Q": Q})
                levels = np.linalg.eigvalsh(sector.hamiltonian(background).toarray())
                energies.extend(levels); charges.extend([Q]*len(levels))
            arrays["unrestricted_energies"] = np.array(energies)
            arrays["unrestricted_charges"] = np.array(charges)
        return arrays, {"variance": spectrum["variance"], "captured_weight": spectrum["captured_weight"],
                        "method": "eigenpairs" if "energies" in spectrum else "Lanczos quadrature measure",
                        "steps": spectrum.get("steps"), "sector_dimension": len(ops.basis)}
    arrays, info, identifier = store.point(specification, compute)
    bright = arrays["weights"] > 1e-8
    return {"point": identifier, "epsilon": background, "E0": E0,
            "bright_gap": float(arrays["gaps"][bright][0]) if bright.any() else None, **info}


def physical_error(base, refined):
    return {name: max(abs(base[name]-refined[name]), 1e-10) for name in ("Y", "W")}


def with_tail(pulse, threshold):
    if isinstance(pulse, Sauter): return replace(pulse, threshold=threshold)
    if isinstance(pulse, Sum): return Sum(tuple(with_tail(p, threshold) for p in pulse.components))
    return pulse


def run_pulse_point(store, cfg):
    """Four matched controls, finite-bandwidth response, optional pump analysis."""
    ops = build_sector_operators(**cfg["model"])
    background = 0.
    E0, vacuum, _ = ground(store, ops, cfg["model"], background)
    strong, weak = pulse_components(cfg["pulse"])
    both = Sum((strong, weak)); interval = both.interval
    pulses = {
        "none": Constant(0., *interval), "strong": strong, "weak": weak, "both": both}
    analysis = cfg.get("analysis", {})
    spectrum = (transition_spectrum(ops, background) if len(ops.basis) <= 1024 else
                lanczos_spectrum(ops, background, steps=analysis.get("lanczos_steps", 128), ground=(E0, vacuum)))
    row = {"dimension": len(ops.basis), "E0": E0, "trajectories": {}}
    for name, pulse in pulses.items():
        arrays, info, identifier = run_trajectory(store, ops, cfg, vacuum, pulse, interval)
        prediction = linear_response(spectrum, pulse)
        record = {"point": identifier, **info, "Y_LR": prediction["Y"], "W_LR": prediction["W"]}
        if len(ops.basis) <= 1024:
            populations = abs(spectrum["vectors"].conj().T@arrays["endpoint_state"])**2
            record["populations"] = populations.tolist()
            record["unresolved_population"] = float(np.vdot(arrays["endpoint_state"], arrays["endpoint_state"]).real-populations.sum())
        if analysis.get("convergence", False):
            refined_settings = {**cfg["integrator"], "max_step": info["integrator"]["max_step"]/2,
                                "resolution_scale": .5,
                                "rtol": cfg["integrator"].get("rtol", 1e-10)/10,
                                "atol": cfg["integrator"].get("atol", 1e-12)/10}
            _, refined, refined_id = run_trajectory(store, ops, cfg, vacuum, pulse, interval, refined_settings)
            record["uncertainty"] = physical_error(info, refined)
            record["refined_point"] = refined_id
        else: record["uncertainty"] = None
        if analysis.get("tail_checks") and cfg["pulse"]["kind"] == "two_sauter":
            record["tail_checks"] = []
            for threshold in analysis["tail_checks"]:
                extended = with_tail(pulse, threshold)
                extended_interval = with_tail(both, threshold).interval
                if name == "none": extended = Constant(0., *extended_interval)
                _, tail, tail_id = run_trajectory(store,ops,cfg,vacuum,extended,extended_interval)
                differences = {x: abs(info[x]-tail[x]) for x in ("Y","W")}
                record["tail_checks"].append({"threshold":threshold,"point":tail_id,**differences})
            if record["uncertainty"] is not None:
                for quantity in ("Y","W"):
                    record["uncertainty"][quantity] += max(x[quantity] for x in record["tail_checks"])
        if analysis.get("independent_midpoint", False):
            midpoint_settings = {**cfg["integrator"],"method":"midpoint","max_step":analysis.get("midpoint_step",.005/ops.g)}
            _,midpoint,midpoint_id = run_trajectory(store,ops,cfg,vacuum,pulse,interval,midpoint_settings)
            record["midpoint_point"] = midpoint_id
            if record["uncertainty"] is not None:
                for quantity in ("Y","W"):
                    record["uncertainty"][quantity] += abs(info[quantity]-midpoint[quantity])
        record["numerically_valid"] = bool(info["max_norm_error"]<1e-9 and info["max_charge_squared"]<1e-12
            and info["max_work_residual"]<max(1e-8,.01*abs(info["W"]))
            and info["post_Y_drift"]<max(1e-8,.01*info["Y"])
            and info["post_W_drift"]<max(1e-8,.01*abs(info["W"])))
        if record["uncertainty"] is not None:
            record["numerically_valid"] &= (record["uncertainty"]["Y"]<max(1e-8,.01*info["Y"])
                and record["uncertainty"]["W"]<max(1e-8,.01*abs(info["W"])))
        row["trajectories"][name] = record
    # Nonadditivity is an absolute, signed difference. No small-denominator
    # enhancement ratio is needed to establish or compare these effects.
    s, w, b = [row["trajectories"][name] for name in ("strong", "weak", "both")]
    row.update({"Y": b["Y"], "W": b["W"], "D_Y": b["Y"]-s["Y"]-w["Y"],
                "D_W": b["W"]-s["W"]-w["W"], "effect_resolved": False})
    if all(item["uncertainty"] is not None and item["numerically_valid"] for item in (s, w, b)):
        for quantity in ("Y", "W"):
            row["D_"+quantity+"_uncertainty"] = sum(item["uncertainty"][quantity] for item in (s,w,b))
        row["effect_resolved"] = abs(row["D_Y"]) > max(1e-8, 5*row["D_Y_uncertainty"])
    if analysis.get("spectral_convergence", False) and len(ops.basis) > 1024:
        doubled = lanczos_spectrum(ops, background, steps=2*analysis.get("lanczos_steps", 128), ground=(E0,vacuum))
        prediction = linear_response(doubled, weak)
        row["spectral_Y_difference"] = abs(prediction["Y"]-linear_response(spectrum, weak)["Y"])
        row["spectral_variance"] = spectrum["variance"]
    if cfg["pulse"]["kind"] == "plateau_sine":
        plateau_ground = sector_ground_state(ops, cfg["pulse"]["epsilon_b"])[1]
        both_arrays, _, _ = run_trajectory(store, ops, cfg, vacuum, both, interval)
        row["ramp_up_ground_overlap"] = float(abs(np.vdot(plateau_ground, both_arrays["ramp_state"]))**2)
    if analysis.get("pump_response", False):
        amplitude = cfg["pulse"]["epsilon_w"]
        if not isinstance(weak, (Sauter, Sine)) or amplitude == 0:
            raise ValueError("Pump response requires nonzero Sauter/sine weak amplitude")
        unit = replace(weak, amplitude=1.)
        specification = {"kind": "pump_expansion", "model": cfg["model"], "strong": strong.metadata(),
                         "unit_weak": unit.metadata(), "interval": interval, "integrator": cfg["integrator"]}
        def compute():
            result = pump_expansion(ops, vacuum, strong, unit, interval,
                max_step=min(cfg["integrator"]["max_step"], step_cap(both,ops.g)),
                rtol=cfg["integrator"].get("rtol",1e-10), atol=cfg["integrator"].get("atol",1e-12))
            return result, {"definition": "psi0+lambda*psi1+lambda**2*psi2; includes psi2 interference"}
        expansion, _, identifier = store.point(specification, compute)
        row["pump_expansion_point"] = identifier
        row["pump_coefficients"] = expansion["coefficients"].tolist()
        row["Y_pump_prediction"] = float(np.polynomial.polynomial.polyval(amplitude, expansion["coefficients"]))
    row["weak_squared_drive"] = squared_drive(weak)
    if analysis.get("controls", False):
        controls = {}
        if isinstance(weak, Sauter) and isinstance(strong, Sauter):
            slow = Sauter(strong.amplitude+weak.amplitude, strong.width,
                          threshold=strong.threshold, taper_widths=strong.taper_widths)
            controls["same_peak_slow"] = slow
            controls["negative_weak"] = Sum((strong, replace(weak, amplitude=-weak.amplitude)))
            for delay in analysis.get("delays", [-strong.width, strong.width]):
                controls[f"delay_{delay}"] = Sum((strong, replace(weak, center=float(delay))))
            reference_tau = analysis.get("budget_reference_tau", 1.)
            # Compute the actual compact-waveform budget, including the taper.
            target = squared_drive(replace(weak, width=reference_tau))
            scaled = weak.amplitude*np.sqrt(target/max(squared_drive(weak), 1e-300))
            if abs(scaled) <= analysis.get("weak_amplitude_ceiling", .04):
                controls["fixed_weak_budget"] = Sum((strong, replace(weak, amplitude=scaled)))
            row["fixed_budget_amplitude"] = scaled
            row["fixed_budget_omitted"] = "fixed_weak_budget" not in controls
            row["fixed_budget_target"] = target
        row["controls"] = {}
        for name, pulse in controls.items():
            control_interval = (min(interval[0], pulse.interval[0]),max(interval[1],pulse.interval[1]))
            _, info, identifier = run_trajectory(store, ops, cfg, vacuum, pulse, control_interval)
            row["controls"][name] = {"point": identifier, "Y": info["Y"], "W": info["W"]}
        if "negative_weak" in row["controls"]:
            plus = row["trajectories"]["both"]["Y"]
            minus = row["controls"]["negative_weak"]["Y"]
            row["Y_even"] = (plus+minus)/2-row["trajectories"]["strong"]["Y"]
            row["Y_odd"] = (plus-minus)/2
    return row


def run(config_path, kind, run_id, indices=None, limit=None, resume=False):
    """Execute the exact rows of one retained recipe; each row is saved at once.

    The JSON scan order is stable (first coordinate changes slowest). A limit
    bounds trial work without changing the physical parameters of those rows.
    Existing data require identical code/configuration and explicit resume.
    """
    config = json.loads(config_path.read_text())
    wanted = None if indices is None else set(indices)
    points = [(i,c,x) for i,c,x in scan_points(config) if wanted is None or i in wanted]
    if limit is not None:
        if limit < 1: raise ValueError("limit must be positive")
        points = points[:limit]
    if not points: raise ValueError("Selection contains no points")
    store = StudyStore(run_id, config, kind, resume,
                       {"indices": [i for i,_,_ in points]})
    for index, point, coordinates in points:
        row = run_spectrum_point(store, point) if kind == "spectrum" else run_pulse_point(store, point)
        row.update({"index": index, "coordinates": coordinates, "model": point["model"], "pulse": point["pulse"]})
        store.table([row])
        print(f"{run_id} [{index}]: {coordinates} Y={row.get('Y')} D_Y={row.get('D_Y')}", flush=True)
    return store.data
