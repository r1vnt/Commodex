"""Executable verification of CrudeAssay unit and price conversions."""

from decimal import Decimal

import pint

from core.units.conversions import CrudeAssay
from core.units.registry import format_quantity, parse_quantity, ureg


def verify_crude_assay() -> None:
    brent = CrudeAssay.brent()
    assert abs(brent.specific_gravity - 0.8333) < 0.0001
    assert abs(brent.bbl_per_mt - 7.55) < 0.01

    cargo_volume = parse_quantity("600,000 bbl")
    fob_price_bbl = parse_quantity("$75.50/bbl")
    with ureg.context(brent.as_context()):
        cargo_mass = cargo_volume.to("mt")
        fob_price_mt = fob_price_bbl.to("USD / tonne")

    expected_mass = 600_000 / brent.bbl_per_mt
    expected_price = 75.50 * brent.bbl_per_mt
    assert abs(float(cargo_mass.magnitude) - expected_mass) < 0.1
    assert abs(float(fob_price_mt.magnitude) - expected_price) < 0.001
    assert format_quantity(cargo_mass, precision=2) == "79415.42 mt"
    assert format_quantity(fob_price_mt, precision=2) == "$570.42/mt"

    try:
        cargo_volume.to("mt")
    except pint.DimensionalityError:
        pass
    else:
        raise AssertionError("Volume-to-mass conversion unexpectedly succeeded outside its assay context.")

    assert isinstance(cargo_mass.magnitude, Decimal)


if __name__ == "__main__":
    verify_crude_assay()
    print("Crude assay verification passed.")
