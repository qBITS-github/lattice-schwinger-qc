"""Actual compact waveforms, their derivatives, supports and normalization.

All times use the repository's numerical units. A config with times in 1/g
must supply the numerical values divided by g when g is changed.
"""
from dataclasses import asdict, dataclass
import math
import numpy as np


def ramp(u):
    if np.isscalar(u):
        u=min(1.,max(0.,u))
        return u**3*(10-15*u+6*u*u)
    u = np.clip(u, 0, 1)
    return u**3*(10-15*u+6*u*u)


def ramp_dot(u):
    if np.isscalar(u): return 30*u*u*(1-u)**2 if 0 < u < 1 else 0.
    return np.where((u > 0) & (u < 1), 30*u*u*(1-u)**2, 0.)


@dataclass(frozen=True)
class Constant:
    value: float
    start: float
    end: float

    @property
    def interval(self): return self.start, self.end
    @property
    def features(self): return [self.end-self.start]
    @property
    def boundaries(self): return [self.start, self.end]
    def epsilon(self, t): return self.value if np.isscalar(t) else np.zeros_like(np.asarray(t), dtype=float)+self.value
    def epsilon_dot(self, t): return 0. if np.isscalar(t) else np.zeros_like(np.asarray(t), dtype=float)
    def metadata(self): return {"kind": "constant", **asdict(self)}


@dataclass(frozen=True)
class Sauter:
    amplitude: float
    width: float
    center: float = 0.
    threshold: float = 1e-4
    taper_widths: float = 1.

    def __post_init__(self):
        if self.width <= 0 or not 0 < self.threshold < 1 or self.taper_widths <= 0:
            raise ValueError("Require width,taper_widths > 0 and 0 < tail threshold < 1")

    @property
    def L(self): return self.width*np.arccosh(1/np.sqrt(self.threshold))
    @property
    def taper_duration(self): return self.width*self.taper_widths
    @property
    def interval(self):
        radius = self.L+self.taper_duration
        return self.center-radius, self.center+radius
    @property
    def boundaries(self):
        return [self.interval[0], self.center-self.L, self.center, self.center+self.L, self.interval[1]]
    @property
    def features(self): return [self.width, self.taper_duration]

    def _parts(self, t):
        x = np.asarray(t)-self.center
        # exp form avoids cosh overflow far outside support.
        e = np.exp(-2*np.abs(x)/self.width)
        shape = 4*e/(1+e)**2
        u = (np.abs(x)-self.L)/self.taper_duration
        return x, shape, u

    def epsilon(self, t):
        if np.isscalar(t):
            L=self.width*math.acosh(1/math.sqrt(self.threshold))
            x=t-self.center
            if abs(x)>=L+self.taper_duration: return 0.
            e=math.exp(-2*abs(x)/self.width)
            return self.amplitude*4*e/(1+e)**2*(1-ramp((abs(x)-L)/self.taper_duration))
        x, shape, u = self._parts(t)
        inside = (np.asarray(t) > self.interval[0]) & (np.asarray(t) < self.interval[1])
        return np.where(inside, self.amplitude*shape*(1-ramp(u)), 0.)

    def epsilon_dot(self, t):
        if np.isscalar(t):
            L=self.width*math.acosh(1/math.sqrt(self.threshold))
            x=t-self.center
            if abs(x)>=L+self.taper_duration: return 0.
            e=math.exp(-2*abs(x)/self.width)
            u=(abs(x)-L)/self.taper_duration
            return self.amplitude*4*e/(1+e)**2*(-2*math.tanh(x/self.width)/self.width*(1-ramp(u))
                                               -ramp_dot(u)*math.copysign(1.,x)/self.taper_duration)
        x, shape, u = self._parts(t)
        inside = (np.asarray(t) > self.interval[0]) & (np.asarray(t) < self.interval[1])
        value = self.amplitude*shape*(-2*np.tanh(x/self.width)/self.width*(1-ramp(u))
                                     - ramp_dot(u)*np.sign(x)/self.taper_duration)
        return np.where(inside, value, 0.)

    def metadata(self):
        return {"kind": "sauter", **asdict(self), "support": self.interval,
                "taper_duration": self.taper_duration, "normalization": "component_peak"}


@dataclass(frozen=True)
class Sine:
    amplitude: float
    omega: float
    duration: float
    phase: float = 0.
    start: float = 0.

    def __post_init__(self):
        if self.duration <= 0 or self.omega <= 0:
            raise ValueError("Sine duration and angular frequency must be positive")

    @property
    def interval(self): return self.start, self.start+self.duration
    @property
    def boundaries(self): return list(self.interval)
    @property
    def features(self): return [self.duration, 1/self.omega]
    def epsilon(self, t):
        if np.isscalar(t):
            x=t-self.start
            return (self.amplitude*math.sin(self.omega*x+self.phase)*math.sin(math.pi*x/self.duration)**2
                    if 0<x<self.duration else 0.)
        x = np.asarray(t)-self.start
        return np.where((x > 0) & (x < self.duration),
                        self.amplitude*np.sin(self.omega*x+self.phase)*np.sin(np.pi*x/self.duration)**2, 0.)
    def epsilon_dot(self, t):
        if np.isscalar(t):
            x=t-self.start
            return (self.amplitude*(self.omega*math.cos(self.omega*x+self.phase)*math.sin(math.pi*x/self.duration)**2
                    +math.sin(self.omega*x+self.phase)*math.pi/self.duration*math.sin(2*math.pi*x/self.duration))
                    if 0<x<self.duration else 0.)
        x = np.asarray(t)-self.start
        value = self.amplitude*(self.omega*np.cos(self.omega*x+self.phase)*np.sin(np.pi*x/self.duration)**2
                   + np.sin(self.omega*x+self.phase)*np.pi/self.duration*np.sin(2*np.pi*x/self.duration))
        return np.where((x > 0) & (x < self.duration), value, 0.)
    def metadata(self): return {"kind": "windowed_sine", **asdict(self), "normalization": "carrier_coefficient"}


@dataclass(frozen=True)
class Plateau:
    amplitude: float
    ramp_duration: float
    hold_duration: float

    def __post_init__(self):
        if self.ramp_duration <= 0 or self.hold_duration <= 0:
            raise ValueError("Ramp and hold durations must be positive")

    @property
    def interval(self): return 0., 2*self.ramp_duration+self.hold_duration
    @property
    def boundaries(self): return [0., self.ramp_duration, self.ramp_duration+self.hold_duration, self.interval[1]]
    @property
    def features(self): return [self.ramp_duration, self.hold_duration]
    def epsilon(self, t):
        if np.isscalar(t):
            return self.amplitude*ramp(t/self.ramp_duration)*ramp((self.interval[1]-t)/self.ramp_duration)
        t = np.asarray(t)
        return self.amplitude*ramp(t/self.ramp_duration)*ramp((self.interval[1]-t)/self.ramp_duration)
    def epsilon_dot(self, t):
        if np.isscalar(t):
            u=t/self.ramp_duration;v=(self.interval[1]-t)/self.ramp_duration
            return self.amplitude/self.ramp_duration*(ramp_dot(u)*ramp(v)-ramp(u)*ramp_dot(v))
        t = np.asarray(t); u = t/self.ramp_duration; v = (self.interval[1]-t)/self.ramp_duration
        return self.amplitude/self.ramp_duration*(ramp_dot(u)*ramp(v)-ramp(u)*ramp_dot(v))
    def metadata(self): return {"kind": "quintic_plateau", **asdict(self)}


@dataclass(frozen=True)
class Sum:
    components: tuple

    @property
    def interval(self): return min(p.interval[0] for p in self.components), max(p.interval[1] for p in self.components)
    @property
    def features(self): return [x for p in self.components for x in p.features]
    @property
    def boundaries(self): return sorted(set(x for p in self.components for x in p.boundaries))
    def epsilon(self, t): return sum(p.epsilon(t) for p in self.components)
    def epsilon_dot(self, t): return sum(p.epsilon_dot(t) for p in self.components)
    def metadata(self): return {"kind": "sum", "components": [p.metadata() for p in self.components], "support": self.interval}


def pulse_components(p):
    """Return strong and weak pulses; zero amplitudes retain matched support."""
    kind = p["kind"]
    if kind == "two_sauter":
        common = {"threshold": p.get("tail_relative_threshold", 1e-4),
                  "taper_widths": p.get("taper_duration_in_component_widths", 1.)}
        return (Sauter(p["epsilon_s"], p["T_s"], **common),
                Sauter(p["epsilon_w"], p["tau"], p.get("t_d", 0.), **common))
    if kind == "windowed_sine":
        weak = Sine(p["epsilon_w"], p["omega"], p["T_w"], p.get("phi", 0.))
        return Constant(0., *weak.interval), weak
    if kind == "plateau_sine":
        strong = Plateau(p["epsilon_b"], p["ramp_duration"], p["T_w"]+2*p.get("buffer", 2.))
        weak = Sine(p["epsilon_w"], p["omega"], p["T_w"], p.get("phi", 0.),
                    p["ramp_duration"]+p.get("buffer", 2.))
        return strong, weak
    raise ValueError(f"Unknown pulse kind: {kind}")
