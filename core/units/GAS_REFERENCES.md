# Contractual gas-reference conditions

`GasReferenceConditions` separates metering temperature/absolute pressure from
combustion reference temperature. ISO 6976 makes this distinction explicitly:
[definitions 3.11–3.12](https://www.iso.org/obp/ui?_escaped_fragment_=iso%3Astd%3Aiso%3A6976%3Aed-3%3Av1%3Aen).
The converter below implements gas-volume parity, not the ISO calorific-value
calculation or a combustion-temperature correction.

```python
from decimal import Decimal
from core.units import Q_, GasReferenceConditions, convert

# Illustrative values; use validated contract/composition data for Z and water.
metered = GasReferenceConditions(
    temperature=Q_(Decimal("310"), "K"),
    absolute_pressure=Q_(Decimal("200"), "kPa"),
    compressibility_factor=Decimal("0.98"),
    water_mole_fraction=Decimal("0.02"),
    combustion_reference_temperature=Q_(Decimal("25"), "degC"),
    standard_id="counterparty-A/metering-schedule-7",
)
delivery = GasReferenceConditions(
    temperature=Q_(Decimal("288.15"), "K"),
    absolute_pressure=Q_(Decimal("101.325"), "kPa"),
    compressibility_factor=Decimal("0.995"),
    water_mole_fraction=Decimal("0"),
    combustion_reference_temperature=Q_(Decimal("25"), "degC"),
    standard_id="counterparty-B/delivery-schedule-3",
)
delivered_volume = convert(
    Q_(Decimal("1000"), "meter**3"), "meter**3",
    source_reference=metered, target_reference=delivery,
)
```

No units are registered for these contracts. `convert_gas_reference_volume()` is
also available directly with `source`, `target` and `conserve` arguments. Carry
the resulting geometric quantity together with its destination conditions and
material identity; a bare `m³` value does not encode either.

## Conservation and required data

Default `dry_gas` conservation holds non-water gas moles fixed:

`V2 = V1 * (P1/P2) * (T2/T1) * (Z2/Z1) * ((1-y1)/(1-y2))`.

Both water mole fractions are required. Zero explicitly means dry; `None` means
unknown. Relative humidity is not a mole fraction and must first be converted
using an appropriate water-vapor/equilibrium model.

Explicit `total_gas` conservation uses the fixed-total-mole equation without the
water factor. It does not perform drying and rejects conflicting known water
fractions. The master router calls this option `gas_conservation`.

Both states require positive finite Z values at conversion time. Pass `1`
explicitly for an ideal-gas assumption. The generic converter never estimates
natural-gas relative density or chooses DAK for hydrogen, helium, CO2 or another
gas. An external EOS must be appropriate for the composition and phase at each
state, including water when present.

Kelvin and equivalent absolute temperature scales must be above absolute zero;
Celsius/Fahrenheit state temperatures are normalized to Kelvin. Temperature
intervals and gauge-pressure labels are rejected. Neutral pressure units in the
explicit `absolute_pressure` field mean absolute pressure; use
`convert_pressure()` first when the measured input is gauge pressure.

The calculation assumes a single gas phase and unchanged non-water composition.
It does not predict condensation, saturation, removal of other components or
phase equilibrium. Different combustion temperatures remain separate metadata;
they never change this metering-volume calculation.

## Existing standards

`GasReferenceConditions.from_standard()` provides `ferc`, `tx`, `iso15`, `din0`,
`gost20` and `atm60f` conveniences. Z is still required explicitly. Named dry
volume units can be bridged only when the supplied temperature, absolute
pressure and dry condition agree with their definitions, independently of the
free-text `standard_id` label.

The existing actual-to-standard and standard-to-actual helpers now delegate the
final volume calculation to this engine. They retain their existing natural-gas
EOS options and explicitly assume dry gas at both endpoints.

Groningen calorific-equivalent volume stays on the calorific assay route. A
P/T/Z-only calculation cannot establish its energy equivalence.

Conditions are frozen dataclasses with defensive copies of temperature and
pressure. Pint quantities themselves remain mutable; the converter revalidates
copies at each call and never mutates the inputs. Decimal inputs produce Decimal
results; floating inputs produce floating results.
