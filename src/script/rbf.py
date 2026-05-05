from math import log
import numpy as np
from dataclasses import dataclass, field
from typing import Callable, Optional, override, Any


@dataclass 
class RBFKernel:
    name : str = ""
    extra_param : dict[str, float | str] | None = None
    phi : Callable[..., Any] | None = None
    dphi : Callable[..., Any] | None = None
    ddphi : Callable[..., Any] | None = None

    def __call__(self, r, p = None):
        return np.vectorize(self.phi)(r, p)

    def d(self, r, p=None):
        return np.vectorize(self.dphi)(r, p)

    def dd(self, r, p=None):
        return np.vectorize(self.ddphi)(r, p)



class _RBFRegistry:
    def __init__(self) -> None:
        self._kernels : dict[str, RBFKernel] = {}

    def register(self, kernel: RBFKernel):
        self._kernels[kernel.name] = kernel 

    def get(self, name: str) -> RBFKernel:
        return self._kernels[name]
    
    def __getattr__(self, name: str, /) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        return self.get(name)

    def names(self):
        return list(self._kernels.keys())


RBF = _RBFRegistry()

def _gen_extra_param(name : str, min, max, default, step) -> dict[str, float | str] :
    return {
            "name" : name,
            "min" : min,
            "max" : max,
            "default" : default,
            "step" : step
            }


RBF.register(RBFKernel(
    name    = "pow3",
    phi     = lambda r, _ : r**3,
    dphi    = lambda r, _ : 3*r**2,
    ddphi   = lambda r, _ : 6*r
    ))
RBF.register(RBFKernel(
    name    = "thinplate",
    phi     = lambda r, _ : 0 if r == 0 else r**2*np.log(r),
    dphi    = lambda r, _ : 0 if r == 0 else 2*r*np.log(r) + r,
    ddphi   = lambda r, _ : 2 * np.log(r) + 3
    ))
RBF.register(RBFKernel(
    name    = "multiquadric",
    extra_param= _gen_extra_param("epsilon", 0.1, 5.0, 0.5, 0.01),
    phi     = lambda r, e : np.sqrt(1 + e**2*r**2),
    dphi    = lambda r, e : e**2*r / np.sqrt(1 + e**2*r**2),
    ddphi   = lambda r, e : e**2 / (1 + e**2*r**2)**1.5,
    ))
RBF.register(RBFKernel(
    name    = "inverse-multiquadric",
    extra_param= _gen_extra_param("epsilon", 0.1, 5.0, 0.5, 0.01),
    phi     = lambda r, e : 1.0 / np.sqrt(1 + e**2*r**2),
    dphi    = lambda r, e : -(e**2)*r / (1 + e**2*r**2)**1.5,
    ddphi   = lambda r, e : e**2*(2*e**2*r**2 - 1) / (1 + e**2*r**2)**2.5
    ))


RBF.register(RBFKernel(
    name    = "gaussian",
    extra_param= _gen_extra_param("sigma", 0.001, 5.0, 0.5, 0.001),
    phi     = lambda r, s : np.exp(-r**2 / (2*s**2)),
    dphi    = lambda r, s : (-r / (s**2)) * np.exp(-r**2 / (2*s**2)),
    ddphi   = lambda r, s : (r**2/s**4 - 1/s**2) * np.exp(-r**2 / (2*s**2))
    ))


RBF.register(RBFKernel(
    name    = "wendlandC2",
    extra_param= _gen_extra_param("support radius", 0.0, 10.0, 2.0, 0.05),
    phi     = lambda r, p : 0.0 if (r/p) >= 1 else (1-(r/p))**4 * (4*(r/p)+1),
    dphi    = lambda r, p : 0.0 if (r/p) >= 1 else -20 * (r/p) * (1-(r/p))**3 / p ,
    ddphi   = lambda r, p : 0.0 if (r/p) >= 1 else 20 * (1-(r/p))**2 * (1 - 4*(r/p)) / p**2 
    ))
