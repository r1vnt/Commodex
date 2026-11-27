"""Commodex precision implementation and public interfaces."""
from __future__ import annotations

from contextlib import contextmanager
from decimal import Context, Decimal, ROUND_HALF_EVEN
from decimal import localcontext
from typing import Any, Final
import math
import pint
import sys

@contextmanager
def decimal_precision(digits=28):
    """Temporarily select Decimal precision without changing global settings."""
    if type(digits) is not int or digits < 1:
        raise ValueError("Decimal precision must be a positive integer.")
    with localcontext() as context:
        context.prec = digits
        yield context

IEEE754_FLOAT64_MANTISSA_BITS: Final[int] = sys.float_info.mant_dig  # 53
IEEE754_FLOAT64_MIN_SIGNIFICANT_DIGITS: Final[int] = sys.float_info.dig  # 15
IEEE754_FLOAT64_DECIMAL_DIGITS: Final[int] = 17  # DBL_DECIMAL_DIG
IEEE754_FLOAT64_EPSILON: Final[Decimal] = Decimal(2) ** -52

def pad_ieee754_decimal(value: Decimal | float | int | str, *, sig_digits: int = IEEE754_FLOAT64_DECIMAL_DIGITS,) -> Decimal:
    """Convert ``value`` to ``Decimal`` with at least ``sig_digits`` significant digits.

    - Inputs with fewer than ``sig_digits`` significant digits are padded with
      trailing zeros to ``sig_digits`` significant digits (for example,
      ``100.00 -> Decimal('100.00000000000000')`` and
      ``0 -> Decimal('0.0000000000000000')``).
    - Inputs that already carry ``>= sig_digits`` significant digits retain all
      of their existing digits without truncation.
    """
    if isinstance(value, bool) or not isinstance(value, (Decimal, int, float, str)):
        raise TypeError(f"Expected Decimal, float, int, or str for IEEE 754 decimal formatting, got {type(value).__name__}.")
    if isinstance(sig_digits, bool) or not isinstance(sig_digits, int) or sig_digits < 1:
        raise ValueError(f"sig_digits must be a positive integer, got {sig_digits!r}.")
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"IEEE 754 decimal value must be finite, got {value!r}.")
        dec = Decimal(str(value))
    elif isinstance(value, Decimal):
        if not value.is_finite():
            raise ValueError(f"IEEE 754 decimal value must be finite, got {value!r}.")
        dec = value
    elif isinstance(value, int):
        dec = Decimal(value)
    else:
        dec = Decimal(value.strip())
        if not dec.is_finite():
            raise ValueError(f"IEEE 754 decimal value must be finite, got {value!r}.")

    if dec.is_zero():
        sign = "-" if dec.is_signed() else ""
        return Decimal(f"{sign}0." + ("0" * max(1, sig_digits - 1)))

    tup = dec.as_tuple()
    exponent = tup.exponent
    assert isinstance(exponent, int)  # Non-finite Decimals were rejected above.
    adj = dec.adjusted()
    # Number of significant digits currently represented when written out
    target_exp = adj - (sig_digits - 1)
    if exponent <= target_exp:
        # Already has >= sig_digits significant digits (and fractional precision)
        return dec

    # For large integers in [10**(sig_digits-1), 10**24] that already have >= sig_digits
    # integer digits, preserve full integer representation plus at least '.0' if exponent > 0.
    if target_exp > 0:
        if exponent <= 0:
            return dec
        if adj <= 24:
            target_exp = 0

    needed_prec = max(38, abs(adj) + abs(target_exp) + sig_digits + 8)
    ctx = Context(prec=needed_prec, rounding=ROUND_HALF_EVEN)
    quantum = Decimal(1).scaleb(target_exp, context=ctx)
    return ctx.quantize(dec, quantum)

def to_ieee754_float(value: Decimal | float | int) -> float:
    """Return the exact IEEE 754 binary64 floating-point value for ``value``."""
    if isinstance(value, bool) or not isinstance(value, (Decimal, int, float)):
        raise TypeError(f"Expected real numeric scalar, got {type(value).__name__}.")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Numeric value must be finite, got {value!r}.")
    return result

def preserve_ieee754_scalar(decimal_value: Decimal, *reference_inputs: Any, sig_digits: int = IEEE754_FLOAT64_DECIMAL_DIGITS,) -> Decimal | float:
    """Return padded ``Decimal`` if any reference input uses ``Decimal``, else exact ``float``."""
    use_decimal = False
    for item in reference_inputs:
        mag = item.magnitude if isinstance(item, pint.Quantity) else item
        if isinstance(mag, Decimal):
            use_decimal = True
            break
    if use_decimal:
        return pad_ieee754_decimal(decimal_value, sig_digits=sig_digits)
    return to_ieee754_float(decimal_value)

def ensure_quantity_ieee754(quantity: pint.Quantity, *, sig_digits: int = IEEE754_FLOAT64_DECIMAL_DIGITS,) -> pint.Quantity:
    """Ensure a Pint Quantity retains IEEE 754 binary64 precision and Decimal trailing zeros."""
    mag = quantity.magnitude
    if isinstance(mag, Decimal):
        padded = pad_ieee754_decimal(mag, sig_digits=sig_digits)
        if padded.as_tuple() == mag.as_tuple():
            return quantity
        return quantity._REGISTRY.Quantity(padded, quantity.units)
    if isinstance(mag, (int, float)) and not isinstance(mag, bool):
        return quantity._REGISTRY.Quantity(to_ieee754_float(mag), quantity.units)
    return quantity

def format_ieee754(value: Decimal | float | int | pint.Quantity, *, sig_digits: int = IEEE754_FLOAT64_DECIMAL_DIGITS,) -> str:
    """Format a scalar or ``pint.Quantity`` with at least ``sig_digits`` significant digits and trailing zeros."""
    if isinstance(value, pint.Quantity):
        mag_str = format_ieee754(value.magnitude, sig_digits=sig_digits)
        return f"{mag_str} {value.units}".strip()
    padded = pad_ieee754_decimal(value, sig_digits=sig_digits)
    # Render without scientific notation for ordinary commodity magnitudes and zero
    adj = padded.adjusted()
    if padded.is_zero() or -12 <= adj <= 24:
        return format(padded, "f")
    return str(padded)

__all__ = ('IEEE754_FLOAT64_DECIMAL_DIGITS', 'IEEE754_FLOAT64_EPSILON', 'IEEE754_FLOAT64_MANTISSA_BITS', 'IEEE754_FLOAT64_MIN_SIGNIFICANT_DIGITS', 'decimal_precision', 'ensure_quantity_ieee754', 'format_ieee754', 'pad_ieee754_decimal', 'preserve_ieee754_scalar', 'to_ieee754_float')
_PUBLIC_ALIASES: dict[str, tuple[str, str]] = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
