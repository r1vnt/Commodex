"""Commodex thermodynamics implementation and public interfaces."""
from __future__ import annotations

from bisect import bisect_right
from collections.abc import Callable, Mapping
from collections.abc import Mapping
from dataclasses import dataclass, field
from dataclasses import dataclass, field, replace
from decimal import Decimal
from decimal import Decimal, InvalidOperation
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from types import MappingProxyType
from typing import Any
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.exceptions import CommodityUnitError


# ---------------------------------------------------------------------------
# Real-Gas Compressibility (Z-Factor) & Actual-to-Standard Volume Engine
# (Dranchuk–Abou-Kassem EoS + Optional AGA-8 Detail / CoolProp HEOS / GERG-2008)
# ---------------------------------------------------------------------------


class SolverError(CommodityUnitError):
    """Base exception for numerical equation-of-state and property solver failures."""


from core.units.exceptions import ConversionConvergenceError


class SolverConvergenceError(SolverError, ConversionConvergenceError):
    """Raised when an iterative numerical solver fails to converge or hits a singular state."""

    def __init__(self, message: str, *, result: SolverResult | None = None) -> None:
        super().__init__(message)
        self.result = result


from core.units.exceptions import InvalidPhysicalStateError


class SolverDomainError(SolverError, InvalidPhysicalStateError):
    """Raised when physical state, composition, or correlation inputs lie outside valid bounds."""





class SolverBackendError(SolverError, ArithmeticError):
    """Raised when an external thermodynamic/EOS backend fails or returns a non-physical value."""





class SolverBackendUnavailableError(SolverBackendError, ImportError):
    """Raised when an explicitly requested optional EOS backend package is not installed."""





@dataclass(frozen=True, slots=True)
class SolverResult:
    """Structured diagnostic result returned by hardened numerical solvers."""

    value: float
    converged: bool
    iterations: int
    residual: float
    method: str
    backend_version: str | None
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise TypeError("SolverResult.value must be a real float.")
        val = float(self.value)
        if not isinstance(self.converged, bool):
            raise TypeError("SolverResult.converged must be a bool.")
        if self.converged and not math.isfinite(val):
            raise ValueError("Converged SolverResult.value must be finite.")
        if isinstance(self.iterations, bool) or not isinstance(self.iterations, int) or self.iterations < 0:
            raise ValueError("SolverResult.iterations must be a non-negative integer.")
        if isinstance(self.residual, bool) or not isinstance(self.residual, (int, float)):
            raise TypeError("SolverResult.residual must be a real float.")
        res = float(self.residual)
        if self.converged and (not math.isfinite(res) or res < 0.0):
            raise ValueError("Converged SolverResult.residual must be finite and non-negative.")
        if not isinstance(self.method, str) or not self.method.strip():
            raise ValueError("SolverResult.method must be a non-empty string.")
        if self.backend_version is not None and (
            not isinstance(self.backend_version, str) or not self.backend_version.strip()
        ):
            raise ValueError("SolverResult.backend_version must be a non-empty string or None.")
        if isinstance(self.warnings, (str, bytes)):
            raise TypeError("SolverResult.warnings must be a tuple/iterable of warning strings, not a bare string.")
        cleaned_warnings = tuple(str(w) for w in self.warnings)
        object.__setattr__(self, "value", val)
        object.__setattr__(self, "residual", res)
        object.__setattr__(self, "method", self.method.strip())
        object.__setattr__(
            self,
            "backend_version",
            self.backend_version.strip() if self.backend_version is not None else None,
        )
        object.__setattr__(self, "warnings", cleaned_warnings)

    def __float__(self) -> float:
        return self.value





def _solve_z_aga8_detail(
    p_kpa: float,
    t_k: float,
    molar_composition: Mapping[str, float],
    prior_warnings: tuple[str, ...] = (),
) -> SolverResult:
    """Strict 21-component AGA-8 Detail via pyaga8; never substitutes another EOS."""
    from core.units.conversions import _AGA8_COMPONENTS
    from core.units.conversions import _checked_backend_z
    from core.units.conversions import _detect_package_version
    from core.units.conversions import _validated_composition
    composition = _validated_composition(molar_composition)
    try:
        import pyaga8
    except ImportError as error:
        raise SolverBackendUnavailableError(
            "method='aga8' requires pyaga8; choose 'heos' or 'gerg2008' explicitly for CoolProp."
        ) from error

    backend_ver = _detect_package_version(pyaga8, "pyaga8")
    # pyaga8's public API takes a Composition record, pressure in kPa and
    # temperature in K. A positional list and the old p/t/density spellings
    # are not accepted by the published 0.1.18 bindings.
    component_fields = (
        "methane", "nitrogen", "carbon_dioxide", "ethane", "propane",
        "isobutane", "n_butane", "isopentane", "n_pentane", "hexane",
        "heptane", "octane", "nonane", "decane", "hydrogen", "oxygen",
        "carbon_monoxide", "water", "hydrogen_sulfide", "helium", "argon",
    )
    try:
        mixture = pyaga8.Composition()
        for component, field_name in zip(_AGA8_COMPONENTS, component_fields):
            setattr(mixture, field_name, composition.get(component, 0.0))
        aga = pyaga8.Detail()
        aga.set_composition(mixture)
        aga.pressure, aga.temperature = p_kpa, t_k
        aga.calc_density()
        aga.calc_properties()
        raw_z = aga.z
    except SolverError:
        raise
    except Exception as error:
        raise SolverBackendError(f"pyaga8 AGA-8 Detail evaluation failed: {error}") from error

    z_val = _checked_backend_z(raw_z, "pyaga8 AGA-8 Detail")
    raw_iters = getattr(aga, "iterations", 0)
    iters = int(raw_iters) if isinstance(raw_iters, int) and not isinstance(raw_iters, bool) and raw_iters >= 0 else 0
    raw_res = getattr(aga, "residual", 0.0)
    res_val = abs(float(raw_res)) if isinstance(raw_res, (int, float)) and not isinstance(raw_res, bool) and math.isfinite(raw_res) else 0.0
    return SolverResult(
        value=z_val,
        converged=True,
        iterations=iters,
        residual=res_val,
        method="aga8",
        backend_version=backend_ver,
        warnings=prior_warnings,
    )





def _solve_z_coolprop(
    p_kpa: float,
    t_k: float,
    molar_composition: Mapping[str, float],
    backend_scheme: str,
    prior_warnings: tuple[str, ...] = (),
) -> SolverResult:
    """Evaluate CoolProp with explicit HEOS:: or GERG2008:: backend identity."""
    from core.units.conversions import _checked_backend_z
    from core.units.conversions import _detect_package_version
    from core.units.conversions import _validated_composition
    if backend_scheme not in {"HEOS", "GERG2008"}:
        raise ValueError(f"Unsupported CoolProp backend scheme: {backend_scheme!r}")
    # Preserve CoolProp's exact model distinction: HEOS:: is never reported as GERG-2008.
    method_label = "HEOS" if backend_scheme == "HEOS" else "GERG-2008"
    norm_comp = _validated_composition(molar_composition)
    try:
        import CoolProp as cp_root
        from CoolProp.CoolProp import PropsSI
    except ImportError as error:
        raise SolverBackendUnavailableError(
            f"CoolProp backend {backend_scheme!r} requires the optional CoolProp package."
        ) from error

    cp_sub = getattr(cp_root, "CoolProp", None)
    backend_ver = _detect_package_version(cp_root, "CoolProp")
    if backend_ver is None and cp_sub is not None:
        get_param = getattr(cp_sub, "get_global_param_string", None)
        if callable(get_param):
            try:
                ver_str = get_param("version")
                if isinstance(ver_str, str) and ver_str.strip():
                    backend_ver = ver_str.strip()
            except Exception:
                backend_ver = None

    cp_map = {
        "CH4": "Methane", "N2": "Nitrogen", "CO2": "CarbonDioxide",
        "C2H6": "Ethane", "C3H8": "Propane", "iC4H10": "IsoButane",
        "nC4H10": "n-Butane", "iC5H12": "Isopentane", "nC5H12": "n-Pentane",
        "C6H14": "n-Hexane", "C7H16": "n-Heptane", "C8H18": "n-Octane",
        "C9H20": "n-Nonane", "C10H22": "n-Decane", "H2": "Hydrogen",
        "O2": "Oxygen", "CO": "CarbonMonoxide", "H2O": "Water",
        "H2S": "HydrogenSulfide", "He": "Helium", "Ar": "Argon",
    }
    mixture = "&".join(f"{cp_map[name]}[{fraction:.17g}]" for name, fraction in norm_comp.items())
    fluid_spec = f"{backend_scheme}::{mixture}"
    try:
        raw_z = PropsSI("Z", "P", p_kpa * 1000.0, "T", t_k, fluid_spec)
    except SolverError:
        raise
    except Exception as error:
        raise SolverBackendError(f"CoolProp {backend_scheme} evaluation failed: {error}") from error

    z_val = _checked_backend_z(raw_z, f"CoolProp {backend_scheme}")
    return SolverResult(
        value=z_val,
        converged=True,
        iterations=0,
        residual=0.0,
        method=method_label,
        backend_version=backend_ver,
        warnings=prior_warnings,
    )





def _solve_z_dak(
    pressure: pint.Quantity,
    temperature: pint.Quantity,
    relative_density: float,
    co2_mol_frac: float,
    h2s_mol_frac: float,
    molar_composition: Mapping[str, float] | None,
    *,
    tol: float,
    max_iter: int,
    raise_on_nonconvergence: bool,
    prior_warnings: tuple[str, ...] = (),
) -> SolverResult:
    """Hardened Dranchuk–Abou-Kassem (1975) + Sutton (1985) + Wichert–Aziz (1972) solver."""
    from core.units.conversions import _DAK_A
    from core.units.conversions import _DAK_UNSUPPORTED_NON_HYDROCARBONS
    from core.units.registry import ureg
    if molar_composition is not None:
        unsupported_non_hc = sorted(set(molar_composition) & _DAK_UNSUPPORTED_NON_HYDROCARBONS)
        if unsupported_non_hc:
            raise SolverDomainError(
                f"DAK/Sutton is a hydrocarbon approximation; use a composition EOS for components {unsupported_non_hc}."
            )

    if isinstance(relative_density, bool) or not isinstance(relative_density, (int, float, Decimal)):
        raise TypeError("relative_density must be a real numeric scalar.")
    sg = float(relative_density)
    if not math.isfinite(sg) or not (0.50 <= sg <= 1.50):
        raise SolverDomainError(f"relative_density out of valid hydrocarbon range [0.50, 1.50]: {relative_density}")

    if (
        isinstance(co2_mol_frac, bool)
        or isinstance(h2s_mol_frac, bool)
        or not isinstance(co2_mol_frac, (int, float, Decimal))
        or not isinstance(h2s_mol_frac, (int, float, Decimal))
    ):
        raise TypeError("co2_mol_frac and h2s_mol_frac must be real numeric scalars.")
    co2 = float(co2_mol_frac)
    h2s = float(h2s_mol_frac)
    if (
        not math.isfinite(co2)
        or not math.isfinite(h2s)
        or not (0.0 <= co2 < 1.0)
        or not (0.0 <= h2s < 1.0)
        or (co2 + h2s) >= 1.0
    ):
        raise SolverDomainError("CO2 and H2S mole fractions must be finite and in [0, 1) with sum < 1.")

    p_psia = float(pressure.to(ureg.psia).magnitude)
    t_rankine = float(temperature.to(ureg.degR).magnitude)

    # 1. Sutton (1985) pseudocritical pressure (psia) and temperature (°R)
    p_pc = 756.8 - 131.0 * sg - 3.6 * (sg ** 2)
    t_pc = 169.2 + 349.5 * sg - 74.0 * (sg ** 2)

    # 2. Wichert–Aziz (1972) sour-gas correction
    dak_warnings: list[str] = list(prior_warnings)
    a_sour = co2 + h2s
    if a_sour > 0.0:
        if a_sour > 0.25:
            dak_warnings.append(
                f"Sour-gas acid mole fraction (CO2 + H2S = {a_sour:.4f}) exceeds primary Wichert-Aziz calibration range (0.25)."
            )
        b_sour = h2s
        eps = 120.0 * (a_sour ** 0.9 - a_sour ** 1.6) + 15.0 * (b_sour ** 0.5 - b_sour ** 4)
        t_pc_adj = t_pc - eps
        denom = t_pc + b_sour * (1.0 - b_sour) * eps
        if denom <= 0.0:
            raise SolverDomainError("Invalid pseudocritical state denominator after sour-gas correction.")
        p_pc_adj = (p_pc * t_pc_adj) / denom
        p_pc, t_pc = p_pc_adj, t_pc_adj

    if not math.isfinite(p_pc) or not math.isfinite(t_pc) or p_pc <= 0.0 or t_pc <= 0.0:
        raise SolverDomainError("Invalid pseudocritical state after sour-gas correction.")

    # 3. Pseudoreduced pressure and temperature domain checks
    p_pr = p_psia / p_pc
    t_pr = t_rankine / t_pc
    if not math.isfinite(t_pr) or not (0.70 <= t_pr <= 5.0):
        raise SolverDomainError(
            f"Pseudoreduced temperature T_pr={t_pr:.4f} is outside valid DAK applicability [0.70, 5.00]."
        )
    if not math.isfinite(p_pr) or p_pr < 0.0 or p_pr > 30.0:
        raise SolverDomainError(
            f"Pseudoreduced pressure P_pr={p_pr:.4f} exceeds valid DAK applicability [0.0, 30.0]."
        )
    if t_pr < 1.0 and p_pr >= 1.0:
        raise SolverDomainError(
            f"Subcritical high-pressure state (T_pr={t_pr:.4f} < 1.0, P_pr={p_pr:.4f} >= 1.0) lies in the liquid/two-phase region outside DAK gas applicability."
        )
    if 0.0 < p_pr < 0.20:
        dak_warnings.append(
            f"Pseudoreduced pressure P_pr={p_pr:.4f} is below the primary DAK (1975) fitted floor of 0.20; low-pressure asymptotic regime applied."
        )
    if t_pr > 3.0:
        dak_warnings.append(
            f"Pseudoreduced temperature T_pr={t_pr:.4f} exceeds the primary DAK (1975) fitted ceiling of 3.00."
        )

    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = _DAK_A
    c1 = a1 + a2 / t_pr + a3 / (t_pr ** 3) + a4 / (t_pr ** 4) + a5 / (t_pr ** 5)
    c2 = a6 + a7 / t_pr + a8 / (t_pr ** 2)
    c3 = a9 * (a7 / t_pr + a8 / (t_pr ** 2))
    c4 = a10 / (t_pr ** 3)
    rhs = 0.27 * p_pr / t_pr
    conv_threshold = tol * max(1.0, abs(rhs))

    def eval_dak(r_val: float) -> tuple[float, float]:
        r2 = r_val * r_val
        exp_t = math.exp(-a11 * r2)
        f = (
            r_val
            + c1 * r2
            + c2 * (r_val ** 3)
            - c3 * (r_val ** 6)
            + c4 * (r_val ** 3) * (1.0 + a11 * r2) * exp_t
            - rhs
        )
        df = (
            1.0
            + 2.0 * c1 * r_val
            + 3.0 * c2 * r2
            - 6.0 * c3 * (r_val ** 5)
            + c4 * r2 * (3.0 + a11 * r2 * (3.0 - 2.0 * a11 * r2)) * exp_t
        )
        return f, df

    # 4. Newton–Raphson solver for reduced density (rho_r)
    rho_r = rhs
    f_val, df_val = eval_dak(rho_r)
    if not math.isfinite(f_val):
        failed = SolverResult(
            value=float("nan"),
            converged=False,
            iterations=0,
            residual=float("inf"),
            method="DAK",
            backend_version=None,
            warnings=(*dak_warnings, "Initial DAK residual evaluated to a non-finite value."),
        )
        if raise_on_nonconvergence:
            raise SolverConvergenceError("DAK solver encountered a non-finite initial residual.", result=failed)
        return failed

    if abs(f_val) <= conv_threshold:
        z_init = 1.0 if rhs == 0.0 else rhs / rho_r
        if not math.isfinite(z_init) or z_init <= 0.0:
            raise SolverConvergenceError(f"DAK solver produced non-physical Z-factor: {z_init!r}.")
        return SolverResult(
            value=z_init,
            converged=True,
            iterations=0,
            residual=abs(f_val),
            method="DAK",
            backend_version=None,
            warnings=tuple(dak_warnings),
        )

    clamped = False
    for iteration in range(1, max_iter + 1):
        if not math.isfinite(df_val) or abs(df_val) < 1e-15:
            est_z = rhs / rho_r if rho_r > 0.0 and math.isfinite(rhs / rho_r) else float("nan")
            failed = SolverResult(
                value=est_z,
                converged=False,
                iterations=iteration,
                residual=abs(f_val),
                method="DAK",
                backend_version=None,
                warnings=(*dak_warnings, "DAK solver encountered a singular or non-finite derivative."),
            )
            if raise_on_nonconvergence:
                raise SolverConvergenceError("DAK solver has a singular/nonfinite derivative.", result=failed)
            return failed

        step = f_val / df_val
        candidate_rho = rho_r - step
        if candidate_rho <= 1e-12:
            clamped = True
            rho_r = 1e-12
        else:
            rho_r = candidate_rho

        if not math.isfinite(rho_r):
            failed = SolverResult(
                value=float("nan"),
                converged=False,
                iterations=iteration,
                residual=abs(f_val),
                method="DAK",
                backend_version=None,
                warnings=(*dak_warnings, "DAK solver produced a non-finite reduced density iterate."),
            )
            if raise_on_nonconvergence:
                raise SolverConvergenceError("DAK solver produced nonfinite reduced density.", result=failed)
            return failed

        try:
            f_val, df_val = eval_dak(rho_r)
        except OverflowError:
            f_val, df_val = float("inf"), float("nan")

        if not math.isfinite(f_val):
            failed = SolverResult(
                value=float("nan"),
                converged=False,
                iterations=iteration,
                residual=float("inf"),
                method="DAK",
                backend_version=None,
                warnings=(*dak_warnings, "DAK residual overflowed during Newton iteration."),
            )
            if raise_on_nonconvergence:
                raise SolverConvergenceError("DAK solver residual overflowed during iteration.", result=failed)
            return failed

        if abs(f_val) <= conv_threshold:
            z_val = rhs / rho_r
            if not math.isfinite(z_val) or z_val <= 0.0:
                raise SolverConvergenceError(f"DAK solver converged to non-physical Z-factor: {z_val!r}.")
            if clamped:
                dak_warnings.append("DAK Newton step clamped at positive reduced-density floor (1e-12).")
            return SolverResult(
                value=z_val,
                converged=True,
                iterations=iteration,
                residual=abs(f_val),
                method="DAK",
                backend_version=None,
                warnings=tuple(dak_warnings),
            )

    if clamped:
        dak_warnings.append("DAK Newton step clamped at positive reduced-density floor (1e-12).")
    dak_warnings.append(f"DAK solver failed to converge within {max_iter} iterations (residual={abs(f_val):.3e}).")
    est_z = rhs / rho_r if rho_r > 0.0 and math.isfinite(rhs / rho_r) else float("nan")
    failed = SolverResult(
        value=est_z,
        converged=False,
        iterations=max_iter,
        residual=abs(f_val),
        method="DAK",
        backend_version=None,
        warnings=tuple(dak_warnings),
    )
    if raise_on_nonconvergence:
        raise SolverConvergenceError(
            f"DAK solver failed to establish residual convergence in {max_iter} iterations.",
            result=failed,
        )
    return failed





def solve_z_factor(
    pressure: pint.Quantity,
    temperature: pint.Quantity,
    relative_density: float = 0.60,
    *,
    method: str = "auto",
    molar_composition: Mapping[str, float] | None = None,
    co2_mol_frac: float | None = None,
    h2s_mol_frac: float | None = None,
    tol: float = 1e-10,
    max_iter: int = 50,
    raise_on_nonconvergence: bool = True,
) -> SolverResult:
    """Solve real-gas compressibility factor (Z) and return structured SolverResult diagnostics.

    Distinguishes strictly between:
    - Convergence: ``SolverResult(converged=True, ...)``
    - Non-convergence: raises ``SolverConvergenceError`` (or returns ``SolverResult(converged=False, ...)``
      when ``raise_on_nonconvergence=False``)
    - Out-of-domain physical/composition inputs: raises ``SolverDomainError``
    - External EOS backend failure / missing backend: raises ``SolverBackendError`` /
      ``SolverBackendUnavailableError``

    CoolProp ``HEOS::`` and ``GERG2008::`` are treated and reported as distinct backends
    (``method="HEOS"`` vs ``method="GERG-2008"``); ``HEOS::`` is never labelled as GERG-2008.
    """
    from core.units.conversions import _absolute_reference_temperature
    from core.units.conversions import _validated_composition
    from core.units.conversions import convert_pressure
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(pressure, pint.Quantity) or not isinstance(temperature, pint.Quantity):
        raise TypeError("pressure and temperature must be pint.Quantity instances.")
    if not isinstance(raise_on_nonconvergence, bool):
        raise TypeError("raise_on_nonconvergence must be a bool.")
    validate_registry_quantity(pressure)
    validate_registry_quantity(temperature)

    if not isinstance(method, str):
        raise TypeError("method must be a string.")
    raw_method = method.strip().lower()
    method_key = "gerg2008" if raw_method == "gerg-2008" else raw_method
    if method_key not in {"ideal", "auto", "dak", "aga8", "heos", "gerg2008"}:
        raise ValueError(f"method must be ideal, auto, dak, aga8, heos or gerg2008, got {method!r}")

    if (
        isinstance(tol, bool)
        or not isinstance(tol, (int, float))
        or not math.isfinite(tol)
        or tol <= 0
        or isinstance(max_iter, bool)
        or not isinstance(max_iter, int)
        or max_iter <= 0
    ):
        raise ValueError("tol must be finite/positive and max_iter a positive integer.")

    if molar_composition is not None:
        molar_composition = _validated_composition(molar_composition)
        for name, supplied in (("CO2", co2_mol_frac), ("H2S", h2s_mol_frac)):
            if supplied is not None and not math.isclose(
                float(supplied), molar_composition.get(name, 0.0), rel_tol=0.0, abs_tol=1e-12
            ):
                raise SolverDomainError(f"Explicit {name} fraction conflicts with complete composition.")
        co2_mol_frac = molar_composition.get("CO2", 0.0)
        h2s_mol_frac = molar_composition.get("H2S", 0.0)
    co2_resolved = 0.0 if co2_mol_frac is None else co2_mol_frac
    h2s_resolved = 0.0 if h2s_mol_frac is None else h2s_mol_frac

    if method_key in {"heos", "gerg2008"} and molar_composition is None:
        raise SolverDomainError("An explicit molar_composition is required for CoolProp HEOS and GERG-2008 models.")
    if method_key == "aga8" and not molar_composition:
        raise SolverDomainError("method='aga8' requires a non-empty `molar_composition` dictionary.")

    try:
        temperature = _absolute_reference_temperature(temperature, "Gas EOS temperature")
        pressure = convert_pressure(pressure, ureg.pascal)
    except (TypeError, pint.DimensionalityError):
        raise
    except ValueError as error:
        raise SolverDomainError(str(error)) from error

    p_kpa = float(pressure.to(ureg.kilopascal).magnitude)
    t_k = float(temperature.to(ureg.kelvin).magnitude)
    if not math.isfinite(p_kpa) or not math.isfinite(t_k) or p_kpa < 0.0 or t_k <= 0.0:
        raise SolverDomainError(f"Invalid thermodynamic state: pressure={p_kpa} kPa, temperature={t_k} K.")

    if method_key == "ideal" or p_kpa == 0.0:
        canonical_zero_method = {
            "ideal": "ideal",
            "auto": "ideal",
            "dak": "DAK",
            "aga8": "aga8",
            "heos": "HEOS",
            "gerg2008": "GERG-2008",
        }[method_key]
        return SolverResult(
            value=1.0,
            converged=True,
            iterations=0,
            residual=0.0,
            method=canonical_zero_method,
            backend_version=None,
            warnings=(),
        )

    if method_key == "heos":
        assert molar_composition is not None
        return _solve_z_coolprop(p_kpa, t_k, molar_composition, "HEOS")
    if method_key == "gerg2008":
        assert molar_composition is not None
        return _solve_z_coolprop(p_kpa, t_k, molar_composition, "GERG2008")
    if method_key == "aga8":
        assert molar_composition is not None
        return _solve_z_aga8_detail(p_kpa, t_k, molar_composition)

    fallback_warnings: list[str] = []
    if method_key == "auto" and molar_composition:
        try:
            return _solve_z_aga8_detail(p_kpa, t_k, molar_composition)
        except SolverBackendUnavailableError:
            pass
        # Optional CoolProp fallback uses HEOS:: and MUST be identified as HEOS, never strict GERG-2008.
        try:
            heos_msg = (
                "AGA-8 backend (pyaga8) unavailable; auto selected CoolProp HEOS (HEOS::), "
                "not AGA-8 or strict GERG-2008."
            )
            res_heos = _solve_z_coolprop(p_kpa, t_k, molar_composition, "HEOS", prior_warnings=(heos_msg,))
            warnings.warn(heos_msg, RuntimeWarning, stacklevel=2)
            return res_heos
        except SolverBackendUnavailableError:
            dak_msg = "AGA-8 and CoolProp HEOS backends unavailable; auto selected DAK, not AGA-8/GERG-2008."
            warnings.warn(dak_msg, RuntimeWarning, stacklevel=2)
            fallback_warnings.append(dak_msg)

    return _solve_z_dak(
        pressure,
        temperature,
        relative_density,
        co2_resolved,
        h2s_resolved,
        molar_composition,
        tol=tol,
        max_iter=max_iter,
        raise_on_nonconvergence=raise_on_nonconvergence,
        prior_warnings=tuple(fallback_warnings),
    )





def calculate_z_factor(
    pressure: pint.Quantity,
    temperature: pint.Quantity,
    relative_density: float = 0.60,
    *,
    method: str = "auto",
    molar_composition: Mapping[str, float] | None = None,
    co2_mol_frac: float | None = None,
    h2s_mol_frac: float | None = None,
    tol: float = 1e-10,
    max_iter: int = 50,
    return_result: bool = False,
) -> float | SolverResult:
    """Computes real-gas compressibility factor (Z), optionally returning full SolverResult diagnostics.

    Methods
    -------
    - "ideal"   : Exact ideal-gas Z = 1.0.
    - "dak"     : Built-in Dranchuk–Abou-Kassem (1975) 11-constant EoS with Sutton (1985)
                  pseudocritical properties and Wichert–Aziz (1972) sour-gas corrections.
    - "aga8"    : Strict 21-component AGA-8 Detail via optional pyaga8.
    - "heos"    : Explicit CoolProp Helmholtz Equation of State (HEOS::); distinct from GERG-2008.
    - "gerg2008": Explicit CoolProp GERG-2008 (GERG2008:: / ISO 20765-2); never falls back to HEOS.
    - "auto"    : With `molar_composition`, prefers AGA-8 Detail (`pyaga8`), then CoolProp `HEOS`
                  (explicitly identified as `method="HEOS"` with a warning, never GERG-2008),
                  then built-in `DAK`.
    """
    if not isinstance(return_result, bool):
        raise TypeError("return_result must be a bool.")
    solved = solve_z_factor(
        pressure,
        temperature,
        relative_density,
        method=method,
        molar_composition=molar_composition,
        co2_mol_frac=co2_mol_frac,
        h2s_mol_frac=h2s_mol_frac,
        tol=tol,
        max_iter=max_iter,
        raise_on_nonconvergence=True,
    )
    return solved if return_result else solved.value


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasPropertyBackend(_safety_QuantitySnapshot):
    method_id: str
    version: str
    source: str
    temperature_range: tuple[pint.Quantity, pint.Quantity]
    pressure_range: tuple[pint.Quantity, pint.Quantity]
    components: frozenset[str]
    gas_constant: pint.Quantity
    z_calculator: Callable[[GasComposition, GasReferenceConditions], Decimal] | None = None
    liquid_density_calculator: Callable[[GasComposition, LiquidState], pint.Quantity] | None = None
    dew_point_calculator: Callable[[GasComposition, pint.Quantity, str], pint.Quantity] | None = None
    ideal: bool = False

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_pressure
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.composition import GAS_MOLAR_MASSES
        from core.units.registry import ureg
        for name in ("method_id", "version", "source"):
            object.__setattr__(self, name, _converters_petroleum_text(getattr(self, name), name))
        temperatures = tuple(_converters_petroleum_temperature(value) for value in self.temperature_range)
        pressures = tuple(_converters_petroleum_pressure(value) for value in self.pressure_range)
        if len(temperatures) != 2 or len(pressures) != 2 or temperatures[0] >= temperatures[1] or pressures[0] >= pressures[1]:
            raise ValueError("Gas backend temperature/pressure bounds must be increasing pairs.")
        object.__setattr__(self, "temperature_range", temperatures)
        object.__setattr__(self, "pressure_range", pressures)
        components = frozenset(self.components)
        if not components or components - GAS_MOLAR_MASSES.keys():
            raise ValueError("Declare a non-empty supported-component set.")
        object.__setattr__(self, "components", components)
        gas_constant = _converters_petroleum_positive_quantity(self.gas_constant, "joule/mol/kelvin", "EOS gas constant")
        object.__setattr__(self, "gas_constant", ureg.Quantity(gas_constant, "joule/mol/kelvin"))
        for name in ("z_calculator", "liquid_density_calculator", "dew_point_calculator"):
            if getattr(self, name) is not None and not callable(getattr(self, name)):
                raise TypeError(f"{name} must be callable or None.")
        if not isinstance(self.ideal, bool):
            raise TypeError("ideal must be boolean.")

    def validate(self, composition: GasComposition, temperature: pint.Quantity, pressure: pint.Quantity) -> None:
        from core.units.composition import GasComposition
        if not isinstance(composition, GasComposition):
            raise TypeError("Supply a GasComposition.")
        unsupported = composition.mole_fractions.keys() - self.components
        if unsupported:
            raise ValueError(f"Backend {self.method_id!r} does not support {sorted(unsupported)}.")
        if not self.temperature_range[0] <= temperature <= self.temperature_range[1] or not self.pressure_range[0] <= pressure <= self.pressure_range[1]:
            raise ValueError("Gas state is outside backend applicability.")

    def z_factor(self, composition: GasComposition, reference: GasReferenceConditions) -> Decimal:
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        self.validate(composition, reference.temperature, reference.absolute_pressure)
        if reference.water_mole_fraction is not None and abs(reference.water_mole_fraction - composition.water_mole_fraction) > Decimal("1e-12"):
            raise ValueError("Metering water fraction conflicts with the complete gas composition.")
        if self.z_calculator is None:
            raise ValueError("Backend does not provide gas compressibility.")
        z = _converters_petroleum_number(self.z_calculator(composition, reference), "Z")
        if z <= 0 or self.ideal and z != 1:
            raise ValueError("Backend returned invalid compressibility.")
        if reference.compressibility_factor is not None and abs(reference.compressibility_factor - z) > max(abs(z), abs(reference.compressibility_factor)) * Decimal("1e-8"):
            raise ValueError("Supplied metering Z conflicts with the selected backend.")
        return z

    def liquid_density(self, composition: GasComposition, state: LiquidState) -> pint.Quantity:
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.registry import ureg
        self.validate(composition, state.temperature, state.absolute_pressure)
        if self.liquid_density_calculator is None:
            raise ValueError("Backend does not provide liquid density; supply a measured cargo density.")
        return ureg.Quantity(_converters_petroleum_positive_quantity(self.liquid_density_calculator(composition, state), "kg/meter**3", "LNG density"), "kg/meter**3")





def _converters_gas_coolprop_mixture(composition: GasComposition, backend: str) -> str:
    from core.units.assays.energy.lng import _converters_gas_CP_NAMES
    return backend + "::" + "&".join(f"{_converters_gas_CP_NAMES[name]}[{float(value):.17g}]" for name, value in composition.mole_fractions.items())





def builtin_gas_backend(method: str, *, temperature_range: tuple[pint.Quantity, pint.Quantity],
                        pressure_range: tuple[pint.Quantity, pint.Quantity]) -> GasPropertyBackend:
    """Explicit EOS choice and operating bounds; never auto-fallback or alias HEOS to GERG."""
    from core.units.composition import GAS_MOLAR_MASSES
    from core.units.conversions import UNIVERSAL_GAS_CONSTANT_J_MOL_K
    from core.units.registry import ureg
    from core.units.thermodynamics import SolverBackendError
    from core.units.thermodynamics import SolverBackendUnavailableError
    from core.units.thermodynamics import SolverDomainError
    from core.units.thermodynamics import calculate_z_factor
    raw_method = method.strip().lower()
    method = "gerg2008" if raw_method == "gerg-2008" else raw_method
    if method not in {"ideal", "dak", "aga8", "gerg2008", "heos"}:
        raise ValueError("Choose ideal, dak, aga8, gerg2008, heos, or provide an explicit GasPropertyBackend.")

    def z_factor(composition, state):
        if method == "ideal":
            return Decimal(1)
        if method == "dak" and any(name in composition.mole_fractions for name in ("H2", "O2", "CO", "H2O", "He", "Ar")):
            raise SolverDomainError("DAK/Sutton is a hydrocarbon approximation; use a composition EOS for this gas.")
        z = Decimal(str(calculate_z_factor(state.absolute_pressure, state.temperature,
                                           relative_density=float(composition.ideal_relative_density), method=method,
                                           molar_composition=composition.mole_fractions,
                                           co2_mol_frac=float(composition.mole_fractions.get("CO2", 0)),
                                           h2s_mol_frac=float(composition.mole_fractions.get("H2S", 0)))))
        if method in {"heos", "gerg2008"}:
            try:
                from CoolProp.CoolProp import PhaseSI
            except ImportError as error:
                raise SolverBackendUnavailableError(f"CoolProp backend {method!r} requires CoolProp.") from error
            try:
                phase = PhaseSI("P", float(state.absolute_pressure.magnitude), "T", float(state.temperature.magnitude),
                                _converters_gas_coolprop_mixture(composition, "GERG2008" if method == "gerg2008" else "HEOS"))
            except Exception as error:
                raise SolverBackendError(f"CoolProp PhaseSI failed for {method!r}: {error}") from error
            if phase not in {"gas", "supercritical_gas", "supercritical"}:
                raise SolverDomainError(f"Requested gas state is not a single gas phase: {phase}.")
        return z

    def liquid_density(composition, state):
        from core.units.registry import ureg
        try:
            from CoolProp.CoolProp import PhaseSI, PropsSI
        except ImportError as error:
            raise SolverBackendUnavailableError(f"CoolProp backend {method!r} requires CoolProp.") from error
        backend = "GERG2008" if method == "gerg2008" else "HEOS"
        mixture = _converters_gas_coolprop_mixture(composition, backend)
        pressure, temperature = float(state.absolute_pressure.magnitude), float(state.temperature.magnitude)
        try:
            phase = PhaseSI("P", pressure, "T", temperature, mixture)
        except Exception as error:
            raise SolverBackendError(f"CoolProp PhaseSI failed for {backend}: {error}") from error
        if phase != "liquid":
            raise SolverDomainError(f"Requested LNG state is not a single liquid phase: {phase}.")
        try:
            rho = PropsSI("Dmass", "P", pressure, "T", temperature, mixture)
        except Exception as error:
            raise SolverBackendError(f"CoolProp PropsSI liquid density failed for {backend}: {error}") from error
        return ureg.Quantity(rho, "kg/meter**3")

    def dew_point(composition, pressure, kind):
        from core.units.registry import ureg
        if kind != "hydrocarbon":
            raise ValueError("A mixture dew flash cannot identify a water dew point; supply a dedicated water-equilibrium backend.")
        if composition.water_mole_fraction:
            raise ValueError("Hydrocarbon dew point requires a declared dry composition or a dedicated multi-condensate backend.")
        if any(composition.mole_fractions.get(name, 0) > 0 for name in ("CO2", "H2S")):
            raise ValueError("Acid-gas dew points require a component-phase-specific backend; a mixture flash is not identified as hydrocarbon condensation.")
        try:
            from CoolProp.CoolProp import PropsSI
        except ImportError as error:
            raise SolverBackendUnavailableError(f"CoolProp backend {method!r} requires CoolProp.") from error
        backend = "GERG2008" if method == "gerg2008" else "HEOS"
        try:
            t_dew = PropsSI("T", "P", float(pressure.magnitude), "Q", 1, _converters_gas_coolprop_mixture(composition, backend))
        except Exception as error:
            raise SolverBackendError(f"CoolProp PropsSI dew-point flash failed for {backend}: {error}") from error
        return ureg.Quantity(t_dew, "K")

    r = {"aga8": Decimal("8.31451"), "gerg2008": Decimal("8.314472")}.get(method, Decimal(str(UNIVERSAL_GAS_CONSTANT_J_MOL_K)))
    backend_desc = {
        "heos": "Explicit CoolProp HEOS (Helmholtz Equation of State; distinct from strict GERG-2008); declared engineering bounds; no automatic EOS substitution",
        "gerg2008": "Explicit CoolProp GERG-2008 (ISO 20765-2); declared engineering bounds; no automatic EOS substitution",
    }.get(method, "Explicit " + method + "; declared engineering bounds; no automatic EOS substitution")
    return GasPropertyBackend(method, "runtime backend", backend_desc,
                              temperature_range, pressure_range, frozenset(GAS_MOLAR_MASSES), ureg.Quantity(r, "joule/mol/kelvin"),
                              z_factor, liquid_density if method in {"heos", "gerg2008"} else None,
                              dew_point if method in {"heos", "gerg2008"} else None, method == "ideal")





def gas_density(composition: GasComposition, reference: GasReferenceConditions, *, backend: GasPropertyBackend) -> pint.Quantity:
    """rho = P*M/(Z*R*T), using declared EOS R and pinned conventional molecular masses."""
    z = backend.z_factor(composition, reference)
    return _converters_gas_density_from_z(composition, reference, backend, z)





def _converters_gas_density_from_z(composition, reference, backend, z):
    from core.units.registry import ureg
    mass = composition.molar_mass.to("kg/mol").magnitude
    density = reference.absolute_pressure.magnitude * mass / (z * backend.gas_constant.magnitude * reference.temperature.magnitude)
    return ureg.Quantity(density, "kg/meter**3")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasProperties(_safety_QuantitySnapshot):
    density: pint.Quantity
    molar_mass: pint.Quantity
    relative_density: Decimal
    relative_density_basis: str
    hhv: pint.Quantity
    lhv: pint.Quantity
    gross_wobbe: pint.Quantity
    net_wobbe: pint.Quantity
    metering_reference: GasReferenceConditions
    combustion_temperature: pint.Quantity
    backend_id: str
    calorific_source: str





def gas_properties(composition: GasComposition, reference: GasReferenceConditions, *, backend: GasPropertyBackend,
                    calorific_values: GasCalorificValues | None = None, air_compressibility: Decimal | None = None,
                    measured_relative_density: Decimal | None = None) -> GasProperties:
    from core.units.assays.energy.petroleum import _converters_petroleum_number
    values = calorific_values if calorific_values is not None else composition.calorific_values()
    if reference.combustion_reference_temperature is not None and abs(reference.combustion_reference_temperature.to("K").magnitude - values.combustion_temperature.magnitude) > Decimal("1e-8"):
        raise ValueError("Combustion reference temperatures differ; normalize the calorific data explicitly.")
    z = backend.z_factor(composition, reference)
    if air_compressibility is not None and measured_relative_density is not None:
        raise ValueError("Supply air compressibility or measured relative density, not both.")
    if measured_relative_density is not None:
        relative = _converters_petroleum_number(measured_relative_density, "measured relative density")
        relative_basis = "measured at the declared metering state"
    else:
        if air_compressibility is None:
            if not backend.ideal:
                raise ValueError("Real-gas Wobbe requires explicit air compressibility or measured relative density.")
            air_compressibility = Decimal(1)
        air_z = _converters_petroleum_number(air_compressibility, "air Z")
        if air_z <= 0:
            raise ValueError("Air Z must be positive.")
        relative = composition.ideal_relative_density * air_z / z
        relative_basis = "molar-mass ratio with declared air/gas Z"
    if relative <= 0:
        raise ValueError("Relative density must be positive.")
    density = _converters_gas_density_from_z(composition, reference, backend, z)
    hhv = (values.hhv * density).to("MJ/meter**3")
    lhv = (values.lhv * density).to("MJ/meter**3")
    calculated_reference = replace(reference, compressibility_factor=z, water_mole_fraction=composition.water_mole_fraction)
    return GasProperties(density, composition.molar_mass, relative, relative_basis, hhv, lhv,
                         hhv / relative.sqrt(), lhv / relative.sqrt(), calculated_reference, values.combustion_temperature,
                         backend.method_id, values.source)





def gas_water_content(composition: GasComposition, reference: GasReferenceConditions, *, backend: GasPropertyBackend,
                       target_unit: str = "mg/meter**3") -> pint.Quantity:
    from core.units.composition import GAS_MOLAR_MASSES
    water_mass_fraction = composition.water_mole_fraction * GAS_MOLAR_MASSES["H2O"] / composition.molar_mass.magnitude
    return (gas_density(composition, reference, backend=backend) * water_mass_fraction).to(target_unit)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasDewPoint(_safety_QuantitySnapshot):
    kind: str
    temperature: pint.Quantity
    absolute_pressure: pint.Quantity
    method_id: str
    source: str





def gas_dew_point(composition: GasComposition, pressure: pint.Quantity, *, kind: str, backend: GasPropertyBackend) -> GasDewPoint:
    """Phase-specific equilibrium result from an explicit backend; no water/HC substitution."""
    from core.units.assays.energy.petroleum import _converters_petroleum_pressure
    from core.units.assays.energy.petroleum import _converters_petroleum_temperature
    pressure = _converters_petroleum_pressure(pressure)
    if kind not in {"water", "hydrocarbon"}:
        raise ValueError("Dew-point kind must be water or hydrocarbon.")
    if kind == "water" and composition.water_mole_fraction == 0:
        raise ValueError("Dry gas has no finite water dew point in this model.")
    if kind == "hydrocarbon" and not any(name in composition.mole_fractions for name in ("CH4", "C2H6", "C3H8", "iC4H10", "nC4H10", "iC5H12", "nC5H12", "C6H14", "C7H16", "C8H18", "C9H20", "C10H22")):
        raise ValueError("A gas without hydrocarbons has no hydrocarbon dew point.")
    if composition.mole_fractions.keys() - backend.components or not backend.pressure_range[0] <= pressure <= backend.pressure_range[1]:
        raise ValueError("Dew-point input is outside backend applicability.")
    if backend.dew_point_calculator is None:
        raise ValueError("Selected backend does not provide the requested phase-equilibrium calculation.")
    temperature = _converters_petroleum_temperature(backend.dew_point_calculator(composition, pressure, kind))
    backend.validate(composition, temperature, pressure)
    return GasDewPoint(kind, temperature, pressure, backend.method_id, backend.source)





@dataclass(frozen=True, slots=True)
class SaturatedGas:
    composition: GasComposition
    reference: GasReferenceConditions
    water_model_id: str


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class WaterSaturationModel(_safety_QuantitySnapshot):
    method_id: str
    source: str
    temperature_range: tuple[pint.Quantity, pint.Quantity]
    pressure_range: tuple[pint.Quantity, pint.Quantity]
    saturation_calculator: Callable[[GasComposition, GasReferenceConditions], Decimal]
    dew_point_calculator: Callable[[GasComposition, pint.Quantity], pint.Quantity]

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_pressure
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        object.__setattr__(self, "method_id", _converters_petroleum_text(self.method_id, "water model ID"))
        object.__setattr__(self, "source", _converters_petroleum_text(self.source, "water model source"))
        temperatures = tuple(_converters_petroleum_temperature(value) for value in self.temperature_range)
        pressures = tuple(_converters_petroleum_pressure(value) for value in self.pressure_range)
        if len(temperatures) != 2 or len(pressures) != 2 or temperatures[0] >= temperatures[1] or pressures[0] >= pressures[1]:
            raise ValueError("Declare increasing water-model bounds.")
        object.__setattr__(self, "temperature_range", temperatures)
        object.__setattr__(self, "pressure_range", pressures)
        if not callable(self.saturation_calculator) or not callable(self.dew_point_calculator):
            raise TypeError("Water model requires saturation and dew-point calculators.")

    def saturate(self, dry_composition: GasComposition, reference: GasReferenceConditions) -> SaturatedGas:
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        if dry_composition.water_mole_fraction != 0:
            raise ValueError("Supply an explicitly dry recipe; call composition.dry() before saturation.")
        if not self.temperature_range[0] <= reference.temperature <= self.temperature_range[1] or not self.pressure_range[0] <= reference.absolute_pressure <= self.pressure_range[1]:
            raise ValueError("Saturation state is outside water-model applicability.")
        water = _converters_petroleum_fraction(self.saturation_calculator(dry_composition, reference), "saturated water fraction", allow_one=False)
        composition = dry_composition.with_water(water)
        # A dry-gas Z is not a wet-gas Z. Require fresh EOS evaluation afterwards.
        reference = replace(reference, water_mole_fraction=water, compressibility_factor=None)
        return SaturatedGas(composition, reference, self.method_id)

    def dew_point(self, composition: GasComposition, pressure: pint.Quantity) -> GasDewPoint:
        from core.units.assays.energy.petroleum import _converters_petroleum_pressure
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        pressure = _converters_petroleum_pressure(pressure)
        if composition.water_mole_fraction == 0:
            raise ValueError("Dry gas has no finite water dew point.")
        if not self.pressure_range[0] <= pressure <= self.pressure_range[1]:
            raise ValueError("Pressure is outside water-model applicability.")
        temperature = _converters_petroleum_temperature(self.dew_point_calculator(composition, pressure))
        if not self.temperature_range[0] <= temperature <= self.temperature_range[1]:
            raise ValueError("Water dew point is outside model applicability; no extrapolation.")
        return GasDewPoint("water", temperature, pressure, self.method_id, self.source)





def ideal_water_saturation_model() -> WaterSaturationModel:
    """Low-pressure liquid-water Antoine/ideal-partial-pressure estimate, 273.15–303 K.

    This does not model ice, fugacity enhancement, hydrates or pipeline-pressure
    gas/water equilibria. Use a dedicated backend for those conditions.
    """
    from core.units.registry import ureg
    a, b, c = Decimal("5.40221"), Decimal("1838.675"), Decimal("-31.737")

    def saturation(composition, reference):
        pressure_bar = Decimal(10) ** (a - b / (reference.temperature.magnitude + c))
        return pressure_bar * 100000 / reference.absolute_pressure.magnitude

    def dew(composition, pressure):
        from core.units.registry import ureg
        partial_bar = composition.water_mole_fraction * pressure.magnitude / 100000
        temperature = b / (a - partial_bar.log10()) - c
        return ureg.Quantity(temperature, "K")

    return WaterSaturationModel("ideal-water-Antoine-273.15-303K", "NIST WebBook water / Bridgeman-Aldrich coefficients; ideal partial pressure; conservative low-pressure engineering domain",
                                (ureg.Quantity(Decimal("273.15"), "K"), ureg.Quantity(303, "K")),
                                (ureg.Quantity(10000, "Pa"), ureg.Quantity(200000, "Pa")), saturation, dew)





def _converters__interpolation_validate_real_scalar(value: Any, label: str) -> Decimal:
    """Validate a finite numeric scalar (int, float, or Decimal), rejecting bool and strings."""
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be an int, float, or Decimal, got {type(value).__name__}.")
    try:
        dec = value if isinstance(value, Decimal) else Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"{label} must be finite, got {value!r}.") from exc
    if not dec.is_finite():
        raise ValueError(f"{label} must be finite, got {value!r}.")
    return dec





def _converters__interpolation_match_scalar_type(result_float: float | Decimal, *sources: Any) -> float | Decimal:
    """Return padded IEEE 754 Decimal if any source is Decimal, else float."""
    from core.units.precision import pad_ieee754_decimal

    use_decimal = False
    for src in sources:
        mag = getattr(src, "magnitude", src)
        if isinstance(mag, Decimal):
            use_decimal = True
            break
    if use_decimal:
        dec = result_float if isinstance(result_float, Decimal) else Decimal(str(result_float))
        return pad_ieee754_decimal(dec)
    return float(result_float)





class MonotonicCubicTable:
    """Shape-preserving strictly monotonic cubic Hermite (PCHIP) interpolator with exact inverse."""

    __slots__ = ("xs", "ys", "hs", "slopes")

    def __init__(self, pairs: list[list[float]] | tuple[tuple[float, float], ...]) -> None:
        if len(pairs) < 2:
            raise ValueError("MonotonicCubicTable requires at least two knot pairs.")
        xs = tuple(float(p[0]) for p in pairs)
        ys = tuple(float(p[1]) for p in pairs)
        n = len(xs)
        for val in xs + ys:
            if not math.isfinite(val):
                raise ValueError("All table coordinates must be finite.")
        hs: list[float] = []
        deltas: list[float] = []
        for i in range(n - 1):
            dx = xs[i + 1] - xs[i]
            dy = ys[i + 1] - ys[i]
            if dx <= 0.0 or dy <= 0.0:
                raise ValueError(
                    f"Table knots must be strictly increasing at index {i}: "
                    f"({xs[i]}, {ys[i]}) -> ({xs[i + 1]}, {ys[i + 1]})."
                )
            hs.append(dx)
            deltas.append(dy / dx)

        if n == 2:
            slopes = [deltas[0], deltas[0]]
        else:
            slopes = [0.0] * n
            # Endpoint slopes (non-centered shape-preserving three-point PCHIP)
            d0 = ((2.0 * hs[0] + hs[1]) * deltas[0] - hs[0] * deltas[1]) / (hs[0] + hs[1])
            if d0 <= 0.0:
                d0 = 0.0
            elif d0 > 3.0 * deltas[0]:
                d0 = 3.0 * deltas[0]
            slopes[0] = d0

            dn = ((2.0 * hs[-1] + hs[-2]) * deltas[-1] - hs[-1] * deltas[-2]) / (hs[-1] + hs[-2])
            if dn <= 0.0:
                dn = 0.0
            elif dn > 3.0 * deltas[-1]:
                dn = 3.0 * deltas[-1]
            slopes[-1] = dn

            # Interior slopes (Brodlie / Fritsch-Butland weighted harmonic mean)
            for k in range(1, n - 1):
                w1 = 2.0 * hs[k] + hs[k - 1]
                w2 = hs[k] + 2.0 * hs[k - 1]
                slopes[k] = (w1 + w2) / ((w1 / deltas[k - 1]) + (w2 / deltas[k]))

        self.xs: tuple[float, ...] = xs
        self.ys: tuple[float, ...] = ys
        self.hs: tuple[float, ...] = tuple(hs)
        self.slopes: tuple[float, ...] = tuple(slopes)

    @property
    def min_x(self) -> float:
        return self.xs[0]

    @property
    def max_x(self) -> float:
        return self.xs[-1]

    @property
    def min_y(self) -> float:
        return self.ys[0]

    @property
    def max_y(self) -> float:
        return self.ys[-1]

    @staticmethod
    def _hermite_val_and_deriv(
        t: float,
        yk: float,
        yk1: float,
        hk_dk: float,
        hk_dk1: float,
    ) -> tuple[float, float]:
        t2 = t * t
        t3 = t2 * t
        val = (
            (2.0 * t3 - 3.0 * t2 + 1.0) * yk
            + (t3 - 2.0 * t2 + t) * hk_dk
            + (-2.0 * t3 + 3.0 * t2) * yk1
            + (t3 - t2) * hk_dk1
        )
        deriv = (
            (6.0 * t2 - 6.0 * t) * (yk - yk1)
            + (3.0 * t2 - 4.0 * t + 1.0) * hk_dk
            + (3.0 * t2 - 2.0 * t) * hk_dk1
        )
        return val, deriv

    def evaluate(self, x: float) -> float:
        """Evaluate monotonic cubic Hermite interpolant at x in [min_x, max_x]."""
        if x < self.xs[0] or x > self.xs[-1]:
            raise ValueError(f"x={x} out of table bounds [{self.xs[0]}, {self.xs[-1]}].")
        if x == self.xs[0]:
            return self.ys[0]
        if x == self.xs[-1]:
            return self.ys[-1]
        k = bisect_right(self.xs, x) - 1
        if k >= len(self.hs):
            k = len(self.hs) - 1
        if x == self.xs[k]:
            return self.ys[k]
        hk = self.hs[k]
        t = (x - self.xs[k]) / hk
        val, _ = self._hermite_val_and_deriv(
            t, self.ys[k], self.ys[k + 1], hk * self.slopes[k], hk * self.slopes[k + 1]
        )
        return val

    def invert(self, y: float) -> float:
        """Invert monotonic cubic Hermite interpolant at y in [min_y, max_y]."""
        if y < self.ys[0] or y > self.ys[-1]:
            raise ValueError(f"y={y} out of table bounds [{self.ys[0]}, {self.ys[-1]}].")
        if y == self.ys[0]:
            return self.xs[0]
        if y == self.ys[-1]:
            return self.xs[-1]
        k = bisect_right(self.ys, y) - 1
        if k >= len(self.hs):
            k = len(self.hs) - 1
        if y == self.ys[k]:
            return self.xs[k]

        yk = self.ys[k]
        yk1 = self.ys[k + 1]
        hk = self.hs[k]
        hk_dk = hk * self.slopes[k]
        hk_dk1 = hk * self.slopes[k + 1]

        lo, hi = 0.0, 1.0
        t = (y - yk) / (yk1 - yk)
        tol = 1e-15 * max(1.0, abs(y))

        for _ in range(24):
            val, deriv = self._hermite_val_and_deriv(t, yk, yk1, hk_dk, hk_dk1)
            err = val - y
            if abs(err) <= tol or (hi - lo) <= 1e-15:
                break
            if err < 0.0:
                lo = t
            else:
                hi = t
            if deriv > 0.0:
                cand = t - (err / deriv)
                t = cand if lo < cand < hi else 0.5 * (lo + hi)
            else:
                t = 0.5 * (lo + hi)

        return self.xs[k] + t * hk


__all__ = ('GasDewPoint', 'GasProperties', 'GasPropertyBackend', 'ISO6976Table', 'SaturatedGas', 'SolverBackendError', 'SolverBackendUnavailableError', 'SolverConvergenceError', 'SolverDomainError', 'SolverError', 'SolverResult', 'ThermodynamicBasis', 'WaterSaturationModel', 'builtin_gas_backend', 'calculate_z_factor', 'convert_combustion_heat', 'gas_density', 'gas_dew_point', 'gas_properties', 'gas_water_content', 'ideal_water_saturation_model', 'solve_z_factor', 'water_latent_heat')
_PUBLIC_ALIASES = {'water_latent_heat': ('heat_content', 'water_latent_heat'), 'convert_combustion_heat': ('heat_content', 'convert_combustion_heat'), 'ISO6976Table': ('heat_content', 'ISO6976Table'), 'ThermodynamicBasis': ('heat_content', 'ThermodynamicBasis')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)

# Static declarations for lazy public interfaces.
from typing import TYPE_CHECKING, Any, Callable
if TYPE_CHECKING:
    from core.units.composition import GasCalorificValues
    from core.units.composition import GasComposition
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.heat_content import ISO6976Table
    from core.units.assays.energy.petroleum import LiquidState
    from core.units.heat_content import ThermodynamicBasis
    from core.units.heat_content import convert_combustion_heat
    from core.units.heat_content import water_latent_heat
