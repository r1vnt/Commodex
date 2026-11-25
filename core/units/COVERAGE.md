# Registry convention audit

This is a coverage audit, not certification of every contract or conversion model.
`tests/core/units/test_registry_coverage.py` checks all 32 categories and published
factor precision. Existing Pint definitions are reused rather than duplicated.

**V**: conventional units verified. **Q**: qualified by an explicit basis,
reference, assay or contract. **G**: a requested general model remains a gap;
unit aliases cannot implement it.

| Category | Status | Coverage and convention boundary |
|---|---|---|
| Mass | V | SI, metric, short/long tons, US/UK hundredweights, grain. |
| Precious metals | Q | Troy ounce, pennyweight, tola, region-named taels, Thai bullion baht, Korean don, carat. Tael and baht factors are not universal product specifications. |
| Liquid volume | V | US liquid and imperial gallons, petroleum barrel, litres and prefixes, cubic metres. |
| Dry volume | V | US dry gallon/bushel, imperial bushel, US peck, hectolitre. Crop-weight bushels are separate mass units. |
| Gas volume | G | Existing named dry reference states and explicit Z-factor conversions work. Arbitrary destination states and humidity correction require an additional reference-state model; the current dispatcher is not a universal humid-gas converter. |
| Energy | V | Joule/Wh, IT/thermochemical/59-F Btu, IT/thermochemical calorie, explicit therm/dekatherm conventions. Unqualified calorie and `cal_it` are not interchangeable. |
| Energy equivalents | Q/G | Named BOE, toe and tce conventions exist. Natural-gas energy needs measured calorific data. GGE/DGE remain specification-dependent; no universal gallon-equivalent factor was invented. |
| Power | Q | Active, reactive and apparent labels; thermal/electrical labels. Power factor and efficiency require models, not unit scaling. |
| Capacity | Q | MW/MWth/MWe, MVA/MVAr, Ah/kAh, MWh. Ah-to-energy requires voltage and operating assumptions. |
| Pressure | Q | Absolute/gauge routing, differential-pressure API and vacuum as negative gauge pressure. Conventional water/mercury pressure units exist; actual hydraulic head requires fluid density and gravity. |
| Temperature | V | K, Celsius, Fahrenheit, Rankine and differences. Offset temperatures are not temperature intervals. |
| Flow | V/Q | Composable mass, volume, molar and energy per time. Semantic bulk-mass rates retain their basis. |
| Viscosity | Q/G | cP, cSt and Pa s verified. SUS and SFS are isolated readings; no ASTM D2161 empirical converter is implemented. |
| Density | Q | SI density units and existing assay/reference-temperature methods. True, apparent, bulk and tapped describe measurement conventions; tapped-density procedures are not implemented by aliases. |
| Concentration | Q | ppm, ppb, explicit parts-per-trillion, isolated mass/mole/volume fractions and percentages, mg/L, mol/L. Changing fraction basis needs composition/density. |
| Chemical strength | Q/G | Molarity is mol/L. Normality, Brix, Plato and heavy/light Baumé are isolated readings. Reaction equivalence and empirical calibration converters remain outside the registry. |
| Quality | Q | API-gravity reading and existing gravity/assay APIs; fineness, purity, moisture and ash retain their commodity and analytical basis. |
| Electrical storage | Q | Ah, kAh, C-rate, Wh/kg and Wh/L. Capacity ratings do not establish usable energy, SOC or degradation. |
| Solar resources | V/Q | W/m² and kWh/m² are directly composable; time integration and plant-output models are separate. |
| Wind resources | V/Q | Speed, power density and fraction/percentage capacity factor. No wind-to-power curve is implied. |
| Weather | V/Q | Degree days/hours, precipitation depth and wind speed. Baselines, observations and accumulation periods belong in the model. |
| Agriculture | Q | Crop-specific mass bushels, named bags/bales, crate/sack counts and dozens. Unspecified package counts never imply a fixed mass. |
| Land | Q | International-foot acre, explicitly named US-survey acre, hectare, Chinese mu, named 4200-m² feddan convention, yield/area. No generic regional area aliases. |
| Forestry | Q | Board feet, cords, cubic metres, ADMT/BDMT. Pulp moisture transformations require an assay. |
| Mining | Q | Separate wet, dry, contained and payable dimensions; DMTU. See SEMANTICS.md. |
| Nuclear | Q | kgU/mtU, SWU, GWd/mtU. U3O8/UF6/UO2 material conversions use UraniumAssay; no fictitious universal compound-to-element unit factor. |
| Freight | Q | TEU/FEU, metric tonne–statute-mile, short-ton–statute-mile, tonne–nautical-mile and Worldscale reading. Route/year flat rates remain external. |
| Shipping | Q | Deadweight and displacement mass; separate gross/net tonnage indices. Vessel geometry and statutory formulas determine indices. |
| Environmental | Q | Separate allowance/offset/certificate/fuel-credit schemes; physical emissions mass stays separate. Eligibility and vintage are external metadata. |
| Financial | Q | Currency/minor units, basis points, percentages and tick counts. Tick value and price grid require the instrument specification. |
| Contracts | Q | Distinct lot/contract/delivery counts. Multipliers, nominal versus actual delivery and product identity remain contract metadata. |
| Water | V/Q | International-foot and explicitly named US-survey acre-feet, m³, ML and flow expressions. Abstraction allocation periods and permits remain external. |

## Ambiguous input and breaking conventions

Market parsing rejects bare `ppt`, `mu` and `ton_mile`. Use `parts_per_trillion`
(or an explicit per-thousand ratio), `mu_china`, `tonne_mile`,
`short_ton_mile` or `tonne_nautical_mile`. The legacy misleading `ton_mile`
alias for a nautical-distance unit was removed. Raw Pint may still recognize
prefixed scientific spellings; market input must use `parse_unit`/`parse_quantity`.

`dwt` retains the registry's historical deadweight meaning. Use `pennyweight`
for precious metals. `acre` and `acre_foot` retain international-foot definitions;
SP 811's explicitly US-survey-foot factors are tested against the new
`acre_us_survey` and `acre_foot_us_survey` names. Do not reinterpret historical
survey data through the unqualified names.

## Source register

- [NIST SP 811 Appendix B](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors): exact factors are distinguished from rounded published factors; tests use published precision rather than treating rounded digits as exact.
- [NIST B.8](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8) and [B.9](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9): conventional factors and quantity classification.
- [ASTM D2161](https://store.astm.org/standards/d2161): Saybolt conversion depends on scale and temperature. Readings have no fixed cSt multiplier.
- [IMO Tonnage Convention](https://www.imo.org/en/about/conventions/pages/international-convention-on-tonnage-measurement-of-ships.aspx): gross and net tonnage are independently calculated from vessel information, not fixed masses or volumes.
- [USDA China Grain and Feed Annual](https://apps.fas.usda.gov/newgainapi/api/Report/DownloadReportByFileName?fileName=Grain+and+Feed+Annual_Beijing_China+-+People%27s+Republic+of_CH2025-0064): Chinese mu convention is one fifteenth hectare.
- [FAO Sudan agricultural report](https://www.fao.org/4/x5536e/x5536e09.htm): the explicitly named feddan convention used here is 4200 m²; it is not a claim about every survey definition.
- [Samsung Gold Exchange's own unit FAQ](https://ssgold.co.kr/shopinfo/faq.html): Korean bullion-market don convention, 3.75 g. This is a merchant convention source, not an SI definition or purity guarantee.

Remaining G items are explicit model gaps, not completed features. The registry
stores qualified readings safely; future conversion models must specify sources,
validity ranges and required reference metadata before bridging their dimensions.
