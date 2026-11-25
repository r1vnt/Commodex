"""Shared calorific values, analytical bases and thermal efficiency.

Measured HHV/GCV and LHV/NCV pairs take precedence. Missing values may be
estimated from fuel-bound hydrogen and moisture using latent heat at 25 C
(2.442 MJ/kg water, 8.936 kg water/kg hydrogen), or an explicitly supplied
NCV/GCV ratio. This latent-water model is not a substitute for an ISO 1928
or ISO 18125 laboratory calculation (constant-volume/pressure corrections).

AR/AD moisture is percent of the wet sample; ash is always percent of dry
matter; hydrogen excludes hydrogen bound in moisture and is on ``basis``.
Reference temperature/absolute pressure describe the supplied density, not
an equation of state. Tagged standard gas volumes require matching density.
See https://www.fao.org/4/j4504e/j4504e08.htm for moisture/NCV basis relations.

Examples (all quantities use the application's singleton registry)::

    from core.units.conversions import BiomassAssay
    from core.units.registry import Q_

    fuel = HeatContentAssay.from_assay(BiomassAssay.enplus_a1_pellets())
    energy = convert_heat_content(Q_(100, "metric_ton"), "MWh", assay=fuel)
    dry_ncv = fuel.calorific_value("ncv", basis="DB")
    output = thermal_to_output(energy, efficiency=Q_(40, "percent"))

For any assay without calorific fields, declare measured values explicitly::

    fuel = HeatContentAssay.from_assay(
        crop, ncv=Q_(18, "MJ/kg"), moisture_pct=12,
        density=Q_(650, "kg/m^3"),
    )
    energy = conversions.convert(Q_(2, "kg"), "MJ", assay=crop, heat_content=fuel)

``convert(..., heat_content=...)`` selects this thermal specification, including
its declared density. Other commodity processing/yield modes remain explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from itertools import permutations
from typing import Any, Final

import pint

from .registry import parse_unit, ureg
from .conversions import _as_decimal, _scale_factor

__all__: Final[tuple[str, ...]] = (
    "HeatContentAssay", "heat_content_from_assay", "gcv_to_ncv", "ncv_to_gcv",
    "as_received_to_dry", "dry_to_as_received", "create_heat_content_context",
    "convert_heat_content", "efficiency_from_heat_rate", "heat_rate_from_efficiency",
    "thermal_to_output", "output_to_thermal",
)


def _number(value: Any, label: str, *, positive: bool = False) -> Decimal:
    result = _as_decimal(value)
    if not result.is_finite() or (result <= 0 if positive else result < 0):
        raise ValueError(f"{label} must be finite and {'positive' if positive else 'non-negative'}.")
    return result


def _fraction(value: Any, label: str) -> Decimal:
    result = _number(value, label) / 100
    if result >= 1:
        raise ValueError(f"{label} must be in [0, 100).")
    return result


def _basis(value: str) -> str:
    result = value.upper().strip()
    result = {"DRY": "DB", "D": "DB", "ADB": "AD"}.get(result, result)
    if result not in {"AR", "AD", "DB", "DAF"}:
        raise ValueError("Analytical basis must be AR, AD, DB or DAF.")
    return result


def _kind(value: str) -> str:
    result = value.lower().strip()
    if result in {"hhv", "gcv", "gar", "gross"}:
        return "gcv"
    if result in {"lhv", "ncv", "nar", "net"}:
        return "ncv"
    raise ValueError("Heating-value basis must be GCV/HHV or NCV/LHV.")


def _quantity(value: pint.Quantity, label: str) -> pint.Quantity:
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{label} must be a Pint quantity with units.")
    _number(value.magnitude, label, positive=True)
    # Use Decimal internally without converting a caller's float arithmetic
    # into an assertion of additional measurement precision.
    return ureg.Quantity(_as_decimal(value.magnitude), value.units)


def _is_volume_dimension(dimensionality: Any) -> bool:
    return dimensionality == (ureg.meter ** 3).dimensionality or (
        len(dimensionality) == 1 and "gas_volume" in next(iter(dimensionality))
        and next(iter(dimensionality.values())) == 1
    )


def _input_type(result: pint.Quantity, source: pint.Quantity) -> pint.Quantity:
    if isinstance(source.magnitude, Decimal):
        return result
    return ureg.Quantity(float(result.magnitude), result.units)


@dataclass(frozen=True, slots=True)
class HeatContentAssay:
    """Calorific data usable alongside any commodity assay.

    ``gcv``/``ncv`` are energy per mass, volume or mole quantities. Supply
    density and/or molar_mass to bridge those bases. ``standard_density`` can
    additionally bridge a registry-tagged standard gas volume. All densities
    describe this declared analytical basis and reference state.

    Hydrogen-carrier adapters report delivered H2 energy after cracking,
    matching HydrogenAssay, rather than the carrier's combustion heat.
    """

    name: str = "custom_heat_content"
    gcv: pint.Quantity | None = None
    ncv: pint.Quantity | None = None
    basis: str = "AR"
    moisture_pct: Decimal | float | None = None
    as_received_moisture_pct: Decimal | float | None = None
    air_dried_moisture_pct: Decimal | float | None = None
    hydrogen_pct: Decimal | float | None = None
    ash_db_pct: Decimal | float | None = None
    density: pint.Quantity | None = None
    standard_density: pint.Quantity | None = None
    molar_mass: pint.Quantity | None = None
    reference_temperature: pint.Quantity | None = None
    reference_pressure: pint.Quantity | None = None
    ncv_to_gcv_ratio: Decimal | float | None = None
    latent_heat_mj_per_kg_water: Decimal | float = Decimal("2.442")
    water_per_kg_hydrogen: Decimal | float = Decimal("8.936")

    def __post_init__(self) -> None:
        object.__setattr__(self, "basis", _basis(self.basis))
        if self.gcv is None and self.ncv is None:
            raise ValueError("Supply at least one measured GCV/HHV or NCV/LHV quantity.")
        for field in ("moisture_pct", "as_received_moisture_pct", "air_dried_moisture_pct", "ash_db_pct"):
            value = getattr(self, field)
            if value is not None:
                _fraction(value, field)
        if self.hydrogen_pct is not None and _number(self.hydrogen_pct, "hydrogen_pct") > 100:
            raise ValueError("hydrogen_pct must be in [0, 100].")
        for field in ("latent_heat_mj_per_kg_water", "water_per_kg_hydrogen"):
            _number(getattr(self, field), field, positive=True)
        if self.ncv_to_gcv_ratio is not None:
            if _number(self.ncv_to_gcv_ratio, "ncv_to_gcv_ratio", positive=True) > 1:
                raise ValueError("ncv_to_gcv_ratio must be in (0, 1].")
        if self.basis in {"DB", "DAF"} and self.moisture_pct not in (None, 0):
            raise ValueError("Dry-basis assays cannot declare nonzero sample moisture.")
        declared_moisture = self.as_received_moisture_pct if self.basis == "AR" else self.air_dried_moisture_pct if self.basis == "AD" else None
        if self.moisture_pct is not None and declared_moisture is not None:
            if _as_decimal(self.moisture_pct) != _as_decimal(declared_moisture):
                raise ValueError("Sample moisture must agree with the moisture declared for its analytical basis.")
        dimensions = set()
        for field in ("density", "standard_density"):
            value = getattr(self, field)
            if value is None:
                continue
            value = _quantity(value, field)
            volume_dim = ureg.kg.dimensionality / value.dimensionality
            if not _is_volume_dimension(volume_dim):
                raise ValueError(f"{field} must have mass per volume units.")
            key = str(volume_dim)
            if key in dimensions:
                raise ValueError("Supply only one density for each volume dimensionality.")
            dimensions.add(key)
            object.__setattr__(self, field, value)
        if self.molar_mass is not None:
            object.__setattr__(self, "molar_mass", _quantity(self.molar_mass, "molar_mass").to("kg/mol"))
        if self.reference_temperature is not None:
            _quantity(self.reference_temperature.to("kelvin"), "reference_temperature")
        if self.reference_pressure is not None:
            pressure = self.reference_pressure
            if any(str(unit).lower().endswith(("psig", "barg", "kpag")) for unit in pressure.units._units):
                raise ValueError("Reference pressure must be absolute, not gauge.")
            _quantity(pressure.to("pascal"), "reference_pressure")
        for field in ("gcv", "ncv"):
            value = getattr(self, field)
            if value is not None:
                object.__setattr__(self, field, self._to_mass_specific(_quantity(value, field)))
        if self.gcv is not None and self.ncv is not None and self.ncv > self.gcv:
            raise ValueError("NCV/LHV cannot exceed GCV/HHV.")

    def _densities(self) -> tuple[pint.Quantity, ...]:
        return tuple(q for q in (self.density, self.standard_density) if q is not None)

    def _to_mass_specific(self, value: pint.Quantity) -> pint.Quantity:
        candidates = [value]
        candidates.extend(value / density for density in self._densities())
        if self.molar_mass is not None:
            candidates.append(value / self.molar_mass)
        for candidate in candidates:
            if candidate.is_compatible_with("MJ/kg"):
                return candidate.to("MJ/kg")
        # Volumetric-only gas and molar-only specifications do not imply a
        # mass density. They still support energy parity on their native basis.
        denominator_dim = ureg.MJ.dimensionality / value.dimensionality
        if _is_volume_dimension(denominator_dim) or denominator_dim == ureg.mol.dimensionality:
            return value
        raise ValueError("Calorific units need energy per mass, or matching density/molar_mass metadata.")

    def _moisture(self, basis: str) -> Decimal:
        if basis in {"DB", "DAF"}:
            return Decimal(0)
        value = self.moisture_pct if basis == self.basis else None
        if basis == "AR":
            value = self.as_received_moisture_pct if self.as_received_moisture_pct is not None else value
        if basis == "AD":
            value = self.air_dried_moisture_pct if self.air_dried_moisture_pct is not None else value
        if value is None:
            raise ValueError(f"Moisture is required for {basis} analytical basis.")
        return _fraction(value, "moisture_pct")

    def _organic_fraction(self, basis: str) -> Decimal:
        if basis == "DAF":
            return Decimal(1)
        if self.ash_db_pct is None:
            raise ValueError("ash_db_pct is required for DAF conversion.")
        return (1 - self._moisture(basis)) * (1 - _fraction(self.ash_db_pct, "ash_db_pct"))

    def _latent_loss(self) -> Decimal:
        if self.hydrogen_pct is None:
            raise ValueError("Missing opposite heating value: supply hydrogen_pct or an explicit NCV/GCV ratio.")
        water = self._moisture(self.basis) + _as_decimal(self.water_per_kg_hydrogen) * _as_decimal(self.hydrogen_pct) / 100
        return _as_decimal(self.latent_heat_mj_per_kg_water) * water

    def calorific_value(self, kind: str = "ncv", *, basis: str | None = None) -> pint.Quantity:
        """Returns MJ/kg (or native volume/molar units); supplied pairs remain exact."""
        kind = _kind(kind)
        value = getattr(self, kind)
        if value is None:
            other = self.gcv if kind == "ncv" else self.ncv
            if self.ncv_to_gcv_ratio is not None:
                ratio = _as_decimal(self.ncv_to_gcv_ratio)
                value = other * (ratio if kind == "ncv" else 1 / ratio)
            else:
                if not other.is_compatible_with("MJ/kg"):
                    raise ValueError("Hydrogen/moisture estimation requires mass-specific calorific data.")
                loss = ureg.Quantity(self._latent_loss(), "MJ/kg")
                value = other - loss if kind == "ncv" else other + loss
            _number(value.magnitude, "Derived calorific value", positive=True)
        target = self.basis if basis is None else _basis(basis)
        if target == self.basis:
            return value
        if not value.is_compatible_with("MJ/kg"):
            raise ValueError("Analytical basis changes require mass-specific calorific data.")
        # Remove sample-water evaporation loss before scaling dry matter;
        # add the target moisture loss afterwards. This also preserves a
        # supplied measured pair without replacing it with an approximation.
        source_m = self._moisture(self.basis)
        target_m = self._moisture(target)
        if "DAF" in {self.basis, target}:
            factor = self._organic_fraction(target) / self._organic_fraction(self.basis)
        else:
            factor = (1 - target_m) / (1 - source_m)
        latent = _as_decimal(self.latent_heat_mj_per_kg_water)
        magnitude = value.magnitude
        if kind == "ncv":
            magnitude += latent * source_m
        magnitude *= factor
        if kind == "ncv":
            magnitude -= latent * target_m
        _number(magnitude, "Converted calorific value", positive=True)
        return ureg.Quantity(magnitude, "MJ/kg")

    def energy_density(self, unit: str | pint.Unit = "MJ/kg", *, kind: str = "ncv", basis: str | None = None) -> pint.Quantity:
        """Converts mass, molar or reference-volume energy density explicitly.

        Density metadata applies only to the declared basis; basis-normalized
        mass values need a new density before volumetric conversion.
        """
        target = parse_unit(unit)
        value = self.calorific_value(kind, basis=basis)
        if value.is_compatible_with(target):
            return value.to(target)
        if basis is not None and _basis(basis) != self.basis:
            raise ValueError("Volumetric/molar conversion after a basis change requires new density/molar metadata.")
        candidates = [value * density for density in self._densities()]
        if self.molar_mass is not None:
            candidates.append(value * self.molar_mass)
        for candidate in candidates:
            if candidate.is_compatible_with(target):
                return candidate.to(target)
        raise ValueError("Target energy-density units require matching density or molar_mass metadata.")

    def as_context(self, *, basis: str = "lhv", name: str | None = None) -> pint.Context:
        return create_heat_content_context(self, basis=basis, name=name)

    @classmethod
    def from_assay(cls, assay: Any, **overrides: Any) -> HeatContentAssay:
        return heat_content_from_assay(assay, **overrides)


def heat_content_from_assay(assay: Any, **overrides: Any) -> HeatContentAssay:
    """Adapts fuel assays, or adds explicit calorific data to any other assay.

    No heat is inferred for metals, ores, crops or nuclear fuel. Their callers
    supply gcv/ncv quantities. Existing carrier recovery and measured/estimated
    pairs are preserved. Overrides are validated by HeatContentAssay.
    """
    from .conversions import (
        GasAssay, CoalAssay, CleanFuelAssay, BiofuelAssay, NGLAssay,
        BiomassAssay, PetcokeAssay, HydrogenAssay,
    )
    if isinstance(assay, HeatContentAssay):
        from dataclasses import replace
        return replace(assay, **overrides)
    data: dict[str, Any] = {"name": getattr(assay, "name", type(assay).__name__)}
    if isinstance(assay, GasAssay):
        data.update(gcv=ureg.Quantity(assay.gcv_btu_per_scf, "btu_it/scf"),
                    ncv=ureg.Quantity(assay.ncv_btu_per_scf, "btu_it/scf"))
        # GasAssay defines scf at 60 F / 14.73 psia. Its conversion methods
        # carry this reference correction; derive a matching mass density.
        if assay.relative_density is not None:
            from .conversions import DRY_AIR_DENSITY_15C_KG_M3, FT3_PER_M3, _DENSITY_RATIO_ISO15
            rho = assay.relative_density * DRY_AIR_DENSITY_15C_KG_M3 / (FT3_PER_M3 * _DENSITY_RATIO_ISO15)
            data["standard_density"] = ureg.Quantity(rho, "kg/scf")
        data.update(reference_temperature=ureg.Quantity(60, "degF"), reference_pressure=ureg.Quantity(14.73, "psia"))
    elif isinstance(assay, CoalAssay):
        data.update(gcv=ureg.Quantity(assay.gj_per_mt_gar, "MJ/kg"), ncv=ureg.Quantity(assay.gj_per_mt_nar, "MJ/kg"),
                    moisture_pct=assay.total_moisture_pct, air_dried_moisture_pct=assay.inherent_moisture_pct,
                    hydrogen_pct=assay.hydrogen_ar_pct)
        if assay.ash_pct is not None:
            data["ash_db_pct"] = _as_decimal(assay.ash_pct) / (1 - _as_decimal(assay.total_moisture_pct) / 100)
    elif isinstance(assay, (CleanFuelAssay, BiofuelAssay, NGLAssay, BiomassAssay, PetcokeAssay, HydrogenAssay)):
        mass = ureg.Quantity(Decimal(1), "kg")
        data.update(gcv=assay.mass_to_energy(mass, basis="hhv") / mass,
                    ncv=assay.mass_to_energy(mass, basis="lhv") / mass)
        if isinstance(assay, (BiomassAssay, PetcokeAssay)):
            data.update(moisture_pct=assay.moisture_pct, ash_db_pct=assay.ash_db_pct,
                        density=ureg.Quantity(assay.bulk_density_kg_per_m3, "kg/m^3"))
        if isinstance(assay, BiomassAssay):
            data["hydrogen_pct"] = _as_decimal(assay.hydrogen_db_pct) * (1 - _as_decimal(assay.moisture_pct) / 100)
        if isinstance(assay, CleanFuelAssay):
            if assay.liquid_density_kg_per_m3 is not None:
                data["density"] = ureg.Quantity(assay.liquid_density_kg_per_m3, "kg/m^3")
            if assay.density_kg_per_nm3 is not None:
                data["standard_density"] = ureg.Quantity(assay.density_kg_per_nm3, "kg/Nm3")
        if isinstance(assay, (BiofuelAssay, NGLAssay, HydrogenAssay)):
            data["density"] = ureg.Quantity(assay.density_kg_per_m3, "kg/m^3")
        if isinstance(assay, NGLAssay) and assay.molecular_weight_g_per_mol is not None:
            data["molar_mass"] = ureg.Quantity(assay.molecular_weight_g_per_mol, "g/mol")
        if isinstance(assay, HydrogenAssay):
            data.update(reference_temperature=ureg.Quantity(assay.reference_temp_c, "degC"),
                        reference_pressure=ureg.Quantity(assay.reference_pressure_bar, "bar"))
    data.update(overrides)
    return HeatContentAssay(**data)


def gcv_to_ncv(gcv: pint.Quantity, **metadata: Any) -> pint.Quantity:
    """Returns supplied NCV or an explicit hydrogen/moisture/ratio estimate."""
    return _input_type(HeatContentAssay(gcv=gcv, **metadata).energy_density(gcv.units, kind="ncv"), gcv)


def ncv_to_gcv(ncv: pint.Quantity, **metadata: Any) -> pint.Quantity:
    return _input_type(HeatContentAssay(ncv=ncv, **metadata).energy_density(ncv.units, kind="gcv"), ncv)


def as_received_to_dry(value: pint.Quantity, moisture_pct: Any, *, kind: str = "gcv") -> pint.Quantity:
    data = {_kind(kind): value}
    return _input_type(HeatContentAssay(**data, moisture_pct=moisture_pct).energy_density(value.units, kind=kind, basis="DB"), value)


def dry_to_as_received(value: pint.Quantity, moisture_pct: Any, *, kind: str = "gcv") -> pint.Quantity:
    # Reverse the wet-to-dry relation, including NCV sample-water loss.
    m = _fraction(moisture_pct, "moisture_pct")
    dry = _quantity(value, "calorific value").to("MJ/kg")
    result = dry * (1 - m)
    if _kind(kind) == "ncv":
        result -= ureg.Quantity(Decimal("2.442") * m, "MJ/kg")
    _number(result.magnitude, "Converted calorific value", positive=True)
    return _input_type(result.to(value.units), value)


def create_heat_content_context(assay: HeatContentAssay, *, basis: str = "lhv", name: str | None = None) -> pint.Context:
    """Bridges energy, mass, molar and explicitly declared volume states."""
    value = assay.calorific_value(basis)
    if not value.is_compatible_with("MJ/kg"):
        denominator = ureg.MJ / value.units
        nodes = [(denominator, value.magnitude), (ureg.MJ, Decimal(1))]
        heat = None
    else:
        heat = value.to("MJ/kg").magnitude
        nodes = [(ureg.kg, heat), (ureg.MJ, Decimal(1))]
    for density in assay._densities() if heat is not None else ():
        # Retain the density's declared volume unit (including standard tags).
        volume = ureg.kg / density.units
        rho = density.to(ureg.kg / volume).magnitude
        nodes.append((volume, heat * rho))
    if assay.molar_mass is not None and heat is not None:
        nodes.append((ureg.mol, heat * assay.molar_mass.to("kg/mol").magnitude))
    context = pint.Context(name or f"heat_content_{assay.name}_{_kind(basis)}")

    def transform(source: pint.Unit, destination: pint.Unit, factor: Decimal):
        def apply(_registry: pint.UnitRegistry, quantity: pint.Quantity) -> pint.Quantity:
            q = quantity.to(source)
            return ureg.Quantity(q.magnitude * _scale_factor(q.magnitude, factor), destination)
        return apply

    for (source, source_heat), (destination, destination_heat) in permutations(nodes, 2):
        context.add_transformation(source.dimensionality, destination.dimensionality,
                                   transform(source, destination, source_heat / destination_heat))
    return context


def convert_heat_content(quantity: pint.Quantity, target_unit: str | pint.Unit, *, assay: Any,
                         basis: str = "lhv", analytical_basis: str | None = None, **metadata: Any) -> pint.Quantity:
    """Converts quantities, rates, quotes or a specific calorific value.

    Commodity assays are adapted automatically; nonfuel assays need explicit
    gcv/ncv metadata. ``basis`` means HHV/LHV; ``analytical_basis`` means AR/DB.
    Basis normalization of calorific densities is explicit, never a hidden
    conversion of two identically dimensioned mass units.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError("Expected a Pint quantity.")
    if not _as_decimal(quantity.magnitude).is_finite():
        raise ValueError("Heat-content conversion input must be finite.")
    spec = heat_content_from_assay(assay, **metadata)
    target = parse_unit(target_unit)
    # A calorific-density input represents the selected heating value on the
    # declared analytical basis. Preserve the supplied value rather than
    # silently substituting the assay's benchmark.
    density_dims = [ureg.MJ.dimensionality / unit.dimensionality for unit, _ in _heat_nodes(spec)]
    if quantity.dimensionality in density_dims:
        from dataclasses import replace
        supplied = replace(spec, gcv=quantity if _kind(basis) == "gcv" else None,
                           ncv=quantity if _kind(basis) == "ncv" else None)
        result = supplied.energy_density(target, kind=basis, basis=analytical_basis)
        return _input_type(result, quantity)
    if analytical_basis is not None and _basis(analytical_basis) != spec.basis:
        raise ValueError("Use calorific-density quantities for analytical basis changes.")
    # Execute the thermal context here rather than calling the master router:
    # conversions.convert delegates energy crossings to this function.
    from .conversions import _convert_contextual_quantity
    context = create_heat_content_context(spec, basis=basis)
    with ureg.context(context):
        result = _convert_contextual_quantity(quantity, target, [context])
    return _input_type(result, quantity)


def _heat_nodes(assay: HeatContentAssay) -> list[tuple[pint.Unit, None]]:
    nodes = [(ureg.kg, None)]
    nodes.extend((ureg.kg / density.units, None) for density in assay._densities())
    if assay.molar_mass is not None:
        nodes.append((ureg.mol, None))
    for value in (assay.gcv, assay.ncv):
        if value is not None:
            nodes.append((ureg.MJ / value.units, None))
    return nodes


def _efficiency(value: Any) -> Decimal:
    if isinstance(value, pint.Quantity):
        value = value.to(ureg.dimensionless).magnitude
    result = _number(value, "Efficiency", positive=True)
    if result > 1:
        raise ValueError("Efficiency must be in (0, 1]; use a percent quantity for percentages.")
    return result


def efficiency_from_heat_rate(heat_rate: pint.Quantity) -> Decimal:
    """Heat rate is fuel energy / output energy, e.g. Btu_it/kWh."""
    rate = _quantity(heat_rate, "Heat rate").to(ureg.dimensionless).magnitude
    return _efficiency(1 / rate)


def heat_rate_from_efficiency(efficiency: Any, unit: str | pint.Unit = "btu_it/kWh") -> pint.Quantity:
    return ureg.Quantity(1 / _efficiency(efficiency), ureg.dimensionless).to(parse_unit(unit))


def thermal_to_output(energy: pint.Quantity, *, efficiency: Any = None,
                      heat_rate: pint.Quantity | None = None, target_unit: str | pint.Unit = "MWh") -> pint.Quantity:
    """Fuel thermal energy to output; use exactly one efficiency or heat rate."""
    return _thermal_parity(energy, efficiency, heat_rate, target_unit, inverse=False)


def output_to_thermal(energy: pint.Quantity, *, efficiency: Any = None,
                      heat_rate: pint.Quantity | None = None, target_unit: str | pint.Unit = "MWh") -> pint.Quantity:
    return _thermal_parity(energy, efficiency, heat_rate, target_unit, inverse=True)


def _thermal_parity(energy: pint.Quantity, efficiency: Any, heat_rate: pint.Quantity | None,
                    target_unit: str | pint.Unit, *, inverse: bool) -> pint.Quantity:
    if (efficiency is None) == (heat_rate is None):
        raise ValueError("Supply exactly one efficiency or heat_rate.")
    eta = _efficiency(efficiency) if heat_rate is None else efficiency_from_heat_rate(heat_rate)
    if energy.dimensionality not in (ureg.joule.dimensionality, ureg.watt.dimensionality):
        raise ValueError("Thermal parity requires energy or power units.")
    _number(energy.magnitude, "Energy")
    result = ureg.Quantity(energy.magnitude * _scale_factor(energy.magnitude, 1 / eta if inverse else eta), energy.units)
    return result.to(parse_unit(target_unit))
