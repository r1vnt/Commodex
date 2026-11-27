"""Commodex uncertainty implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, fields, is_dataclass
from decimal import Context, Decimal, ROUND_HALF_EVEN
from enum import Enum
from functools import lru_cache
from hashlib import sha256
from types import MappingProxyType
from typing import Any, Final
import json
import math
import pint




# ---------------------------------------------------------------------------
# First-Order Taylor Error Propagation (Independent & Correlated Covariance)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class UncertaintyPropagationResult:
    """Detailed first-order Taylor uncertainty breakdown with covariance terms."""

    value: Decimal
    standard_uncertainty: Decimal
    variance: Decimal
    relative_uncertainty: Decimal | None
    sensitivities: Mapping[str, Decimal]
    variance_contributions: Mapping[str, Decimal]
    covariance_contribution: Decimal





def _audit_lookup_pair(
    matrix: Mapping[tuple[str, str], Any] | Sequence[Sequence[Any]] | None,
    names: Sequence[str],
    i: int,
    j: int,
    *,
    label: str,
) -> Decimal | None:
    from core.units.provenance import _audit_dec
    if matrix is None:
        return None
    if isinstance(matrix, Mapping):
        ni, nj = names[i], names[j]
        has_ij = (ni, nj) in matrix
        has_ji = (nj, ni) in matrix
        if has_ij and has_ji:
            vij = _audit_dec(matrix[(ni, nj)], f"{label}[{ni}, {nj}]")
            vji = _audit_dec(matrix[(nj, ni)], f"{label}[{nj}, {ni}]")
            if abs(vij - vji) > Decimal("1e-15") * max(Decimal(1), abs(vij), abs(vji)):
                raise ValueError(f"Asymmetric {label} entries for ({ni!r}, {nj!r}): {vij} != {vji}.")
            return vij
        if has_ij:
            return _audit_dec(matrix[(ni, nj)], f"{label}[{ni}, {nj}]")
        if has_ji:
            return _audit_dec(matrix[(nj, ni)], f"{label}[{nj}, {ni}]")
        return None
    if len(matrix) != len(names) or any(len(row) != len(names) for row in matrix):
        raise ValueError(f"{label} matrix dimensions must match {len(names)}x{len(names)}.")
    vij = _audit_dec(matrix[i][j], f"{label}[{i}][{j}]")
    vji = _audit_dec(matrix[j][i], f"{label}[{j}][{i}]")
    if abs(vij - vji) > Decimal("1e-15") * max(Decimal(1), abs(vij), abs(vji)):
        raise ValueError(f"Asymmetric {label} matrix at ({i}, {j}): {vij} != {vji}.")
    return vij





def _validate_correlation_psd(matrix):
    """LDL decomposition of the whole dimensionless covariance shape.

    Checking only J*Cov*J.T misses an indefinite matrix when the chosen
    sensitivity vector happens to point in a positive-variance direction.
    Normalisation prevents tiny but valid standard uncertainties being treated
    as absolute zero merely because another input has a much larger scale.
    """
    from decimal import localcontext
    from core.units.provenance import _audit_HIGH_PREC_CTX
    with localcontext(_audit_HIGH_PREC_CTX):
        size = len(matrix)
        lower = [[Decimal(0) for _ in range(size)] for _ in range(size)]
        diagonal = [Decimal(0) for _ in range(size)]
        tolerance = Decimal("1e-24")
        for i in range(size):
            for j in range(i + 1):
                residual = matrix[i][j] - sum(
                    (lower[i][k] * lower[j][k] * diagonal[k] for k in range(j)), Decimal(0)
                )
                if i == j:
                    if residual < -tolerance:
                        raise ValueError("Covariance matrix must be positive semi-definite.")
                    diagonal[i] = max(Decimal(0), residual)
                    lower[i][i] = Decimal(1)
                elif diagonal[j] > tolerance:
                    lower[i][j] = residual / diagonal[j]
                elif abs(residual) > tolerance:
                    raise ValueError("Covariance matrix must be positive semi-definite.")


def propagate_uncertainty(
    sensitivities: Mapping[str, Decimal | float | int],
    standard_uncertainties: Mapping[str, Decimal | float | int],
    *,
    value: Decimal | float | int = Decimal(0),
    covariance: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
    correlations: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
) -> UncertaintyPropagationResult:
    """Propagate standard uncertainties via first-order Taylor expansion with optional covariance.

    Computes ``sigma_y^2 = sum_i (df/dx_i)^2 * sigma_i^2 + 2 * sum_{i<j} (df/dx_i)(df/dx_j) Cov(x_i, x_j)``.
    Supply at most one of ``covariance`` or ``correlations``.
    """
    from core.units.precision import pad_ieee754_decimal
    from core.units.provenance import _audit_HIGH_PREC_CTX
    from core.units.provenance import _audit_dec
    if covariance is not None and correlations is not None:
        raise ValueError("Supply either covariance or correlations, not both.")

    unknown = set(standard_uncertainties) - set(sensitivities)
    if unknown:
        raise ValueError(f"Uncertainties supplied for parameters without sensitivities: {sorted(unknown)}.")

    if isinstance(covariance, Mapping):
        known_names = set(sensitivities)
        bad_cov_keys = {k for pair in covariance for k in pair if k not in known_names}
        if bad_cov_keys:
            raise ValueError(f"Covariance references unknown parameters: {sorted(bad_cov_keys)}.")
    if isinstance(correlations, Mapping):
        known_names = set(sensitivities)
        bad_corr_keys = {k for pair in correlations for k in pair if k not in known_names}
        if bad_corr_keys:
            raise ValueError(f"Correlations reference unknown parameters: {sorted(bad_corr_keys)}.")

    names = tuple(sensitivities.keys())
    jac: dict[str, Decimal] = {}
    sigmas: dict[str, Decimal] = {}
    for name in names:
        jac[name] = _audit_dec(sensitivities[name], f"sensitivity[{name}]")
        raw_sigma = standard_uncertainties.get(name, Decimal(0))
        sig = _audit_dec(raw_sigma, f"standard_uncertainty[{name}]")
        if sig < Decimal(0):
            raise ValueError(f"Standard uncertainty for {name!r} must be non-negative, got {raw_sigma!r}.")
        sigmas[name] = sig

    # Validate diagonal overrides if a full covariance/correlation matrix was passed
    for idx, name in enumerate(names):
        if covariance is not None:
            diag_cov = _audit_lookup_pair(covariance, names, idx, idx, label="covariance")
            if diag_cov is not None:
                if diag_cov < Decimal(0):
                    raise ValueError(f"Diagonal variance for {name!r} must be non-negative, got {diag_cov}.")
                expected_var = _audit_HIGH_PREC_CTX.multiply(sigmas[name], sigmas[name])
                if name in standard_uncertainties and abs(diag_cov - expected_var) > Decimal("1e-12") * max(Decimal(1), expected_var):
                    raise ValueError(
                        f"Diagonal covariance for {name!r} ({diag_cov}) conflicts with sigma^2 ({expected_var})."
                    )
                sigmas[name] = _audit_HIGH_PREC_CTX.sqrt(diag_cov)
        if correlations is not None:
            diag_corr = _audit_lookup_pair(correlations, names, idx, idx, label="correlations")
            if diag_corr is not None and abs(diag_corr - Decimal(1)) > Decimal("1e-12"):
                raise ValueError(f"Diagonal correlation for {name!r} must be 1.0, got {diag_corr}.")

    diag_contribs: dict[str, Decimal] = {}
    diag_sum = Decimal(0)
    for name in names:
        term = _audit_HIGH_PREC_CTX.multiply(
            _audit_HIGH_PREC_CTX.multiply(jac[name], jac[name]),
            _audit_HIGH_PREC_CTX.multiply(sigmas[name], sigmas[name]),
        )
        diag_contribs[name] = pad_ieee754_decimal(term)
        diag_sum = _audit_HIGH_PREC_CTX.add(diag_sum, term)

    correlation_matrix = [[Decimal(0) for _ in names] for _ in names]
    for i, name in enumerate(names):
        correlation_matrix[i][i] = Decimal(0) if sigmas[name].is_zero() else Decimal(1)
    cov_sum = Decimal(0)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            ni, nj = names[i], names[j]
            si, sj = sigmas[ni], sigmas[nj]
            max_cov = _audit_HIGH_PREC_CTX.multiply(si, sj)
            cov_ij = Decimal(0)
            if covariance is not None:
                found_cov = _audit_lookup_pair(covariance, names, i, j, label="covariance")
                if found_cov is not None:
                    if abs(found_cov) > max_cov * Decimal("1.000000000001") + Decimal("1e-30"):
                        raise ValueError(
                            f"|Cov({ni}, {nj})| = {abs(found_cov)} exceeds Cauchy-Schwarz bound sigma_i*sigma_j = {max_cov}."
                        )
                    cov_ij = found_cov
            elif correlations is not None:
                found_rho = _audit_lookup_pair(correlations, names, i, j, label="correlations")
                if found_rho is not None:
                    if abs(found_rho) > Decimal(1):
                        raise ValueError(
                            f"Correlation rho({ni}, {nj}) must lie in [-1, 1], got {found_rho}."
                        )
                    cov_ij = _audit_HIGH_PREC_CTX.multiply(found_rho, max_cov)

            if max_cov.is_zero():
                if not cov_ij.is_zero():
                    raise ValueError("Covariance involving zero standard uncertainty must be zero.")
            else:
                correlation_matrix[i][j] = correlation_matrix[j][i] = _audit_HIGH_PREC_CTX.divide(cov_ij, max_cov)
            if not cov_ij.is_zero():
                cross = _audit_HIGH_PREC_CTX.multiply(
                    Decimal(2),
                    _audit_HIGH_PREC_CTX.multiply(_audit_HIGH_PREC_CTX.multiply(jac[ni], jac[nj]), cov_ij),
                )
                cov_sum = _audit_HIGH_PREC_CTX.add(cov_sum, cross)

    if covariance is not None or correlations is not None:
        _validate_correlation_psd(correlation_matrix)
    total_var = _audit_HIGH_PREC_CTX.add(diag_sum, cov_sum)
    if total_var < Decimal(0):
        # Allow tiny rounding cancellation near perfectly anti-correlated singularities, reject non-PSD matrices
        scale_ref = max(Decimal("1e-30"), abs(diag_sum) * Decimal("1e-12"))
        if abs(total_var) <= scale_ref:
            total_var = Decimal(0)
        else:
            raise ValueError(
                f"Implied covariance matrix is not positive semi-definite (propagated variance = {total_var} < 0)."
            )

    std_u = _audit_HIGH_PREC_CTX.sqrt(total_var)
    val_dec = _audit_dec(value, "value")
    rel_u = (
        pad_ieee754_decimal(_audit_HIGH_PREC_CTX.divide(std_u, abs(val_dec)))
        if not val_dec.is_zero()
        else None
    )
    return UncertaintyPropagationResult(
        value=pad_ieee754_decimal(val_dec),
        standard_uncertainty=pad_ieee754_decimal(std_u),
        variance=pad_ieee754_decimal(total_var),
        relative_uncertainty=rel_u,
        sensitivities=MappingProxyType({k: pad_ieee754_decimal(v) for k, v in jac.items()}),
        variance_contributions=MappingProxyType(diag_contribs),
        covariance_contribution=pad_ieee754_decimal(cov_sum),
    )





def propagate_callable_uncertainty(
    func: Callable[[Mapping[str, Decimal]], Decimal | float | pint.Quantity],
    inputs: Mapping[str, Decimal | float | int],
    standard_uncertainties: Mapping[str, Decimal | float | int],
    *,
    covariance: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
    correlations: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
    step_fraction: Decimal = Decimal("1e-7"),
) -> UncertaintyPropagationResult:
    """Propagate parameter uncertainties through ``func(params)`` via central finite differences."""
    from core.units.provenance import _audit_HIGH_PREC_CTX
    from core.units.provenance import _audit_dec
    base_params = {k: _audit_dec(v, f"inputs[{k}]") for k, v in inputs.items()}
    unknown = set(standard_uncertainties) - set(base_params)
    if unknown:
        raise ValueError(f"standard_uncertainties keys missing from inputs: {sorted(unknown)}.")

    def _eval_scalar(p: Mapping[str, Decimal]) -> Decimal:
        from core.units.provenance import _audit_dec
        raw = func(p)
        if isinstance(raw, pint.Quantity):
            return _audit_dec(raw.magnitude, "func output")
        return _audit_dec(raw, "func output")

    y0 = _eval_scalar(base_params)
    sensitivities: dict[str, Decimal] = {}
    for name, x0 in base_params.items():
        if name not in standard_uncertainties and covariance is None and correlations is None:
            sensitivities[name] = Decimal(0)
            continue
        sigma = _audit_dec(standard_uncertainties.get(name, Decimal(0)), f"sigma[{name}]")
        step = max(abs(x0) * step_fraction, sigma * Decimal("1e-4"), Decimal("1e-10"))
        plus_p = dict(base_params)
        minus_p = dict(base_params)
        plus_p[name] = _audit_HIGH_PREC_CTX.add(x0, step)
        try:
            y_plus = _eval_scalar(plus_p)
        except Exception:
            y_plus = y0
            plus_p[name] = x0
        minus_cand = _audit_HIGH_PREC_CTX.subtract(x0, step)
        minus_p[name] = minus_cand
        try:
            y_minus = _eval_scalar(minus_p)
            denom = _audit_HIGH_PREC_CTX.subtract(plus_p[name], minus_p[name])
        except Exception:
            # One-sided forward difference at domain boundary (e.g. x0 == 0)
            y_minus = y0
            denom = _audit_HIGH_PREC_CTX.subtract(plus_p[name], x0)
        if denom.is_zero():
            sensitivities[name] = Decimal(0)
        else:
            sensitivities[name] = _audit_HIGH_PREC_CTX.divide(_audit_HIGH_PREC_CTX.subtract(y_plus, y_minus), denom)

    return propagate_uncertainty(
        sensitivities,
        standard_uncertainties,
        value=y0,
        covariance=covariance,
        correlations=correlations,
    )





def uncertain_quantity(nominal, standard_uncertainty, unit, *, tag=None):
    from core.units.registry import Q_
    from core.units.registry import parse_unit
    if isinstance(nominal, bool) or isinstance(standard_uncertainty, bool):
        raise TypeError("Uncertainty values must be numbers, not booleans.")
    value, uncertainty = float(nominal), float(standard_uncertainty)
    if not math.isfinite(value) or not math.isfinite(uncertainty) or uncertainty < 0:
        raise ValueError("Supply finite nominal value and nonnegative standard uncertainty.")
    try:
        from uncertainties import ufloat
    except ImportError as error:
        raise ImportError("Install commodex[uncertainty] for uncertain quantities.") from error
    return Q_(ufloat(value, uncertainty, tag=tag), parse_unit(unit))


__all__ = ('UncertaintyPropagationResult', 'propagate_callable_uncertainty', 'propagate_uncertainty', 'uncertain_quantity')
_PUBLIC_ALIASES = {'Q_': ('registry', 'Q_'), 'parse_unit': ('registry', 'parse_unit')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
