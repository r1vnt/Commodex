# Semantic units and commodity identity

Pint compatibility is not evidence of economic fungibility. A matching unit
does not establish a matching commodity, specification, delivery location,
origin, vintage, certification scheme, or settlement contract. Those identifiers
belong on the instrument or physical-lot model, outside `Quantity` and
`QuoteComponents`. Validate them before adding positions or comparing prices.

Carbon programs and renewable-certificate schemes have separate dimensions.
There is no automatic exchange ratio between schemes, no automatic surrender
conversion into emissions mass, and no certificate-to-electricity conversion.
Such operations require an explicit policy outside the unit registry.

## Breaking bulk-mass migration

`wmt`/`WMT`/`wet_metric_tonne` have dimension `[wet_bulk_mass]`.
`dmt`/`DMT`/`dry_metric_tonne` have dimension `[dry_bulk_mass]`.
`contained_metal_tonne` and `payable_metal_tonne` have independent contained and
payable mass dimensions. Each represents one tonne on its declared basis.
`dmtu` is 0.01 contained-metal tonnes. Existing physical `mtu` remains 10 kg;
it is not repurposed as wet or dry ore.

Plain Pint conversion, arithmetic and compatibility checks cannot equate these
bases, even when moisture happens to be zero or grade happens to be 100%.
Base-unit reduction retains the semantic dimension.

Use `convert(..., assay=ore)` or `convert_bulk_ore(...)` for explicit conversions.
Wet/dry changes require moisture; contained-metal changes require the component
grade; payable metal requires `OreAssay` smelter terms. Rates and currency quotes
use the same explicit factors. Missing information is an error.

Ore and battery assay wet/dry and contained/payable methods now return quantities
tagged with their resulting basis. `extract_contained_metal()` also defaults to
`contained_metal_tonne`, rather than an unqualified physical tonne. Preserve these
units in serialized data and migrate persisted records deliberately; old records
that discarded their basis cannot be reconstructed from the unit string alone.

For backward compatibility, explicit assay methods can accept unqualified physical
mass in the role named by the method. They reject an explicitly conflicting basis.
Internally, validated assay operations can project a semantic tonne to 1,000 kg
for physical scaling. That projection supplies no moisture, grade, or payability
factor and is not available as an ordinary Pint conversion. It never projects
carbon instruments or certificates into physical commodities.
