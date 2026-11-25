"""commodities.core.units.registry.

Institutional-grade unit registry for physical and financial commodities.
Configures independent currency dimensions, thermodynamic gas volume standards,
petroleum volumes, flow rates, environmental instruments, and strict serialisation.
"""

from __future__ import annotations

from contextlib import contextmanager
import copy
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from functools import cache, wraps
import re
import threading
from typing import Any, Final, Iterator
import pint
from ._safety import _SEMANTIC_SCOPE, authorized_conversion

# ---------------------------------------------------------------------------
# ISO-4217 Currency Codes
# ---------------------------------------------------------------------------

CURRENCY_CODES: Final[tuple[str, ...]] = tuple(
    """AED AFN ALL AMD AOA ARS AUD AWG AZN BAM BBD BDT BHD BIF BMD BND BOB BRL BSD BTN BWP BYN BZD CAD CDF CHF CLP CNY COP CRC CUP CVE CZK DJF DKK DOP DZD
       EGP ERN ETB EUR FJD FKP GBP GEL GHS GIP GMD GNF GTQ GYD HKD HNL HTG HUF IDR ILS INR IQD IRR ISK JMD JOD JPY KES KGS KHR KMF KPW KRW KWD KYD KZT LAK
       LBP LKR LRD LSL LYD MAD MDL MGA MKD MMK MNT MOP MRU MUR MVR MWK MXN MYR MZN NAD NGN NIO NOK NPR NZD OMR PAB PEN PGK PHP PKR PLN PYG QAR RON RSD RUB
       RWF SAR SBD SCR SDG SEK SGD SHP SLE SOS SRD SSP STN SVC SYP SZL THB TJS TMT TND TOP TRY TTD TWD TZS UAH UGX USD UYU UZS VES VND VUV WST XAF XCD XCG
       XOF XPF YER ZAR ZMW ZWG""".split()
)

# Display aliases and vendor symbols mapped to ISO codes.
CURRENCY_SYMBOL_TO_CODE: Final[dict[str, str]] = {
    "kr": "SEK", "Rs": "INR", "\u20a8": "INR", "\u00a5": "JPY",
    # Pricing agency and market wire compound symbols.
    "US$": "USD",
    "AU$": "AUD",
    "CA$": "CAD",
    "NZ$": "NZD",
    "SG$": "SGD",
    "CN\u00a5": "CNY",
    "JP\u00a5": "JPY",
    # --- AED (United Arab Emirates Dirham) ---
    # #UAE #MiddleEast #Petrodollar #CrudeOil #RefinedProducts #USDPegged #DubaiGold
    "د.إ": "AED",
    # --- AFN (Afghan Afghani) ---
    # #Afghanistan #CentralAsia #MiningPotential #LithiumReserves #PreciousGems
    "؋": "AFN",
    # --- ALL (Albanian Lek) ---
    # #Albania #Balkans #ChromiumMining #Hydroelectric #AgriExports
    "L": "ALL",
    # --- AMD (Armenian Dram) ---
    # #Armenia #Caucasus #CopperMining #Molybdenum #PreciousMetals
    "֏": "AMD",
    # --- AOA (Angolan Kwanza) ---
    # #Angola #CentralAfrica #OPEC #CrudeOilExporter #Diamonds #PetroEconomy
    "Kz": "AOA",
    # --- ARS (Argentine Peso) ---
    # #Argentina #SouthAmerica #AgriCommodities #Soybeans #Corn #LithiumTriangle
    "$_ARS": "ARS",
    # --- AUD (Australian Dollar) ---
    # #Australia #CommodityCurrency #IronOre #LNGExporter #Coal #Gold #G10
    "A$": "AUD",
    # --- AWG (Aruban Florin) ---
    # #Aruba #Caribbean #TourismEconomy #USDPegged #RefineryHistory
    "ƒ": "AWG",
    # --- AZN (Azerbaijani Manat) ---
    # #Azerbaijan #CaspianSea #CrudeOil #NaturalGas #SouthernGasCorridor
    "₼": "AZN",
    # --- BAM (Bosnia and Herzegovina Convertible Mark) ---
    # #BosniaHerzegovina #Balkans #BauxiteAluminium #Hydropower #EuroPegged
    "KM": "BAM",
    # --- BBD (Barbadian Dollar) ---
    # #Barbados #Caribbean #SugarCane #Tourism #USDPegged
    "Bds$": "BBD",
    # --- BDT (Bangladeshi Taka) ---
    # #Bangladesh #SouthAsia #Textiles #JuteExport #AgriImports #LNGImporter
    "৳": "BDT",
    # --- BHD (Bahraini Dinar) ---
    # #Bahrain #GulfStates #AluminiumSmelting #OilRefining #PetroEconomy #USDPegged
    "د.ب": "BHD",
    # --- BIF (Burundian Franc) ---
    # #Burundi #EastAfrica #CoffeeTeaExport #NickelReserves #ZeroDecimal
    "FBu": "BIF",
    # --- BMD (Bermudian Dollar) ---
    # #Bermuda #OffshoreFinance #InsuranceReinsurance #USDPegged #MaritimeRegistry
    "BD$": "BMD",
    # --- BND (Brunei Dollar) ---
    # #Brunei #SoutheastAsia #PetroleumExporter #LNGExporter #SGDPegged
    "B$": "BND",
    # --- BOB (Bolivian Boliviano) ---
    # #Bolivia #SouthAmerica #LithiumReserves #NaturalGas #TinMining #Zinc
    "Bs.": "BOB",
    # --- BRL (Brazilian Real) ---
    # #Brazil #SouthAmerica #CommodityPowerhouse #Soybeans #IronOre #CrudeOil #Sugar #Coffee
    "R$": "BRL",
    # --- BSD (Bahamian Dollar) ---
    # #Bahamas #Caribbean #Tourism #OffshoreBanking #USDPegged #Aragonite
    "BS$": "BSD",
    # --- BTN (Bhutanese Ngultrum) ---
    # #Bhutan #SouthAsia #Hydroelectricity #INRPegged #Forestry #Himalayas
    "Nu.": "BTN",
    # --- BWP (Botswana Pula) ---
    # #Botswana #SouthernAfrica #DiamondMining #CopperNickel #BeefExport
    "P": "BWP",
    # --- BYN (Belarusian Ruble) ---
    # #Belarus #EasternEurope #PotashFertiliser #HeavyMachinery #OilRefining
    "Br": "BYN",
    # --- BZD (Belize Dollar) ---
    # #Belize #CentralAmerica #SugarCane #CitrusExports #CrudePetroleum #USDPegged
    "BZ$": "BZD",
    # --- CAD (Canadian Dollar) ---
    # #Canada #CommodityCurrency #CrudeOil #OilSands #NaturalGas #Wheat #Potash #G10
    "C$": "CAD",
    # --- CDF (Congolese Franc) ---
    # #DRC #CentralAfrica #CobaltCapital #CopperMining #Coltan #CriticalMinerals
    "FC": "CDF",
    # --- CHF (Swiss Franc) ---
    # #Switzerland #SafeHaven #GoldRefiningHub #CommodityTradingHub #G10
    "Fr": "CHF",
    # --- CLP (Chilean Peso) ---
    # #Chile #CopperPowerhouse #LithiumProducer #LatinAmerica #ZeroDecimal
    "$_CLP": "CLP",
    # --- CNY (Chinese Yuan Renminbi) ---
    # #China #MajorEconomy #LargestCommodityImporter #Petroyuan #MetalsConsumer
    "¥_CNY": "CNY",
    # --- COP (Colombian Peso) ---
    # #Colombia #SouthAmerica #CrudeOil #ThermalCoal #ArabicaCoffee #Nickel
    "$_COP": "COP",
    # --- CRC (Costa Rican Colón) ---
    # #CostaRica #CentralAmerica #Pineapples #Bananas #Coffee #RenewableEnergy
    "₡": "CRC",
    # --- CUP (Cuban Peso) ---
    # #Cuba #Caribbean #NickelMining #SugarCane #CobaltReserves #Tobacco
    "$_CUP": "CUP",
    # --- CVE (Cape Verdean Escudo) ---
    # #CapeVerde #WestAfrica #MaritimeLogistics #Fisheries #EuroPegged
    "Esc": "CVE",
    # --- CZK (Czech Koruna) ---
    # #Czechia #CentralEurope #IndustrialManufacturing #AutomotiveSupply #CoalMining
    "Kč": "CZK",
    # --- DJF (Djiboutian Franc) ---
    # #Djibouti #HornOfAfrica #PortLogistics #RedSeaTrade #USDPegged #ZeroDecimal
    "Fdj": "DJF",
    # --- DKK (Danish Krone) ---
    # #Denmark #Nordic #OffshoreWind #GreenEnergy #EuroPegged #MaritimeShipping
    "kr_DKK": "DKK",
    # --- DOP (Dominican Peso) ---
    # #DominicanRepublic #Caribbean #Ferronickel #GoldMining #SugarCane #Cacao
    "RD$": "DOP",
    # --- DZD (Algerian Dinar) ---
    # #Algeria #NorthAfrica #OPEC #NaturalGasExporter #CrudeOil #PipelinesToEurope
    "دج": "DZD",
    # --- EGP (Egyptian Pound) ---
    # #Egypt #NorthAfrica #SuezCanal #LNGExporter #Fertilisers #WheatImporter
    "E£": "EGP",
    # --- ERN (Eritrean Nakfa) ---
    # #Eritrea #HornOfAfrica #PotashMining #GoldCopper #RedSeaCoast
    "Nfk": "ERN",
    # --- ETB (Ethiopian Birr) ---
    # #Ethiopia #EastAfrica #ArabicaOrigin #HydroelectricPower #Oilseeds #Floriculture
    "Br_ETB": "ETB",
    # --- EUR (Euro) ---
    # #Eurozone #MajorFX #ReserveCurrency #CarbonAllowances #EUA #IndustrialHub
    "€": "EUR",
    # --- FJD (Fijian Dollar) ---
    # #Fiji #Pacific #SugarExport #GoldMining #Fisheries #Tourism
    "FJ$": "FJD",
    # --- FKP (Falkland Islands Pound) ---
    # #FalklandIslands #SouthAtlantic #Fisheries #SquidExport #OffshoreHydrocarbons #GBPPegged
    "£_FKP": "FKP",
    # --- GBP (British Pound Sterling) ---
    # #UnitedKingdom #MajorFX #G10 #BrentBenchmark #LME #NBPGas #UKA
    "£": "GBP",
    # --- GEL (Georgian Lari) ---
    # #Georgia #Caucasus #TransitCorridor #CopperOres #Ferroalloys #WineExport
    "₾": "GEL",
    # --- GHS (Ghanaian Cedi) ---
    # #Ghana #WestAfrica #GoldProducer #CocoaExporter #CrudeOil #Bauxite
    "₵": "GHS",
    # --- GIP (Gibraltar Pound) ---
    # #Gibraltar #Mediterranean #BunkeringHub #MaritimeLogistics #GBPPegged
    "£_GIP": "GIP",
    # --- GMD (Gambian Dalasi) ---
    # #Gambia #WestAfrica #GroundnutExport #Fisheries #AgriTrade
    "D": "GMD",
    # --- GNF (Guinean Franc) ---
    # #Guinea #WestAfrica #BauxiteGiant #AluminiumSupply #IronOreSimandou #ZeroDecimal
    "FG": "GNF",
    # --- GTQ (Guatemalan Quetzal) ---
    # #Guatemala #CentralAmerica #Cardamom #Coffee #SugarCane #NickelMining
    "Q": "GTQ",
    # --- GYD (Guyanese Dollar) ---
    # #Guyana #SouthAmerica #OffshoreOilBoom #StabroekBlock #GoldMining #Bauxite
    "GY$": "GYD",
    # --- HKD (Hong Kong Dollar) ---
    # #HongKong #EastAsia #CommodityTrading #GoldTrading #USDPeggedCurrencyBoard
    "HK$": "HKD",
    # --- HNL (Honduran Lempira) ---
    # #Honduras #CentralAmerica #CoffeeExport #Bananas #PalmOil #ZincLead
    "L_HNL": "HNL",
    # --- HTG (Haitian Gourde) ---
    # #Haiti #Caribbean #EssentialOils #Vetiver #AgriImports
    "G": "HTG",
    # --- HUF (Hungarian Forint) ---
    # #Hungary #CentralEurope #Manufacturing #AutomotiveSupply #NaturalGasTransit
    "Ft": "HUF",
    # --- IDR (Indonesian Rupiah) ---
    # #Indonesia #SoutheastAsia #NickelPowerhouse #ThermalCoal #PalmOil #Tin
    "Rp": "IDR",
    # --- ILS (Israeli New Shekel) ---
    # #Israel #MiddleEast #OffshoreGas #LeviathanField #PolishedDiamonds #HighTech
    "₪": "ILS",
    # --- INR (Indian Rupee) ---
    # #India #SouthAsia #MajorCommodityConsumer #CrudeOilImporter #GoldConsumer #CoalBuyer
    "₹": "INR",
    # --- IQD (Iraqi Dinar) ---
    # #Iraq #MiddleEast #OPEC #BasrahCrude #MegaOilFields #PetroState
    "ع.د": "IQD",
    # --- IRR (Iranian Rial) ---
    # #Iran #MiddleEast #CrudeOilReserves #NaturalGasGiant #Petrochemicals #Copper
    "﷼_IRR": "IRR",
    # --- ISK (Icelandic Króna) ---
    # #Iceland #Nordic #AluminiumSmelting #Geothermal #Fisheries #ZeroDecimal
    "kr_ISK": "ISK",
    # --- JMD (Jamaican Dollar) ---
    # #Jamaica #Caribbean #BauxiteAlumina #Coffee #SugarCane #Tourism
    "J$": "JMD",
    # --- JOD (Jordanian Dinar) ---
    # #Jordan #MiddleEast #PotashMining #PhosphateRock #Fertilisers #USDPegged
    "د.ا": "JOD",
    # --- JPY (Japanese Yen) ---
    # #Japan #MajorFX #G10 #ReserveCurrency #ZeroDecimal #MajorCommodityImporter #LNGImporter
    "¥_JPY": "JPY",
    # --- KES (Kenyan Shilling) ---
    # #Kenya #EastAfrica #TeaExporter #CutFlowers #GeothermalPower #RegionalHub
    "KSh": "KES",
    # --- KGS (Kyrgyzstani Som) ---
    # #Kyrgyzstan #CentralAsia #KumtorGoldMine #Hydropower #RareEarths
    "сом": "KGS",
    # --- KHR (Cambodian Riel) ---
    # #Cambodia #SoutheastAsia #RiceExport #Rubber #Garments #AgriCommodities
    "៛": "KHR",
    # --- KMF (Comorian Franc) ---
    # #Comoros #IndianOcean #VanillaExporter #YlangYlang #Cloves #ZeroDecimal
    "CF": "KMF",
    # --- KPW (North Korean Won) ---
    # #NorthKorea #EastAsia #AnthraciteCoal #IronOre #Magnesite #MiningResources
    "₩_KPW": "KPW",
    # --- KRW (South Korean Won) ---
    # #SouthKorea #EastAsia #RefiningPowerhouse #Petrochemicals #LNGImporter #ZeroDecimal
    "₩_KRW": "KRW",
    # --- KWD (Kuwaiti Dinar) ---
    # #Kuwait #MiddleEast #OPEC #CrudeOilExporter #SovereignWealth #HighestValuedCurrency
    "د.ك": "KWD",
    # --- CI$ (Cayman Islands Dollar) ---
    # #CaymanIslands #Caribbean #OffshoreFunds #FinancialServices #USDPegged
    "CI$": "KYD",
    # --- KZT (Kazakhstani Tenge) ---
    # #Kazakhstan #CentralAsia #UraniumGiant #CPCBlend #CrudeOil #Ferroalloys #Copper
    "₸": "KZT",
    # --- LAK (Lao Kip) ---
    # #Laos #SoutheastAsia #HydropowerToThailand #CopperMining #Gold #Potash
    "₭": "LAK",
    # --- LBP (Lebanese Pound) ---
    # #Lebanon #MiddleEast #CommercialTrading #GoldReserves #AgriImports
    "ل.ل": "LBP",
    # --- LKR (Sri Lankan Rupee) ---
    # #SriLanka #SouthAsia #CeylonTea #Rubber #Cinnamon #Apparel
    "Rs_LKR": "LKR",
    # --- LRD (Liberian Dollar) ---
    # #Liberia #WestAfrica #IronOre #NaturalRubber #Gold #MaritimeFlagRegistry
    "L$": "LRD",
    # --- LSL (Lesotho Loti) ---
    # #Lesotho #SouthernAfrica #DiamondMining #WaterExportToSouthAfrica #ZARPpegged
    "L_LSL": "LSL",
    # --- LYD (Libyan Dinar) ---
    # #Libya #NorthAfrica #OPEC #LightSweetCrude #MellitahGas #PetroEconomy
    "ل.د": "LYD",
    # --- MAD (Moroccan Dirham) ---
    # #Morocco #NorthAfrica #PhosphateCapital #OCPGroup #FertiliserPowerhouse #Automotive
    "د.م.": "MAD",
    # --- MDL (Moldovan Leu) ---
    # #Moldova #EasternEurope #AgriExports #WineIndustry #Sunflowers #GrainTrade
    "L_MDL": "MDL",
    # --- MGA (Malagasy Ariary) ---
    # #Madagascar #IndianOcean #VanillaMonopoly #NickelCobalt #GraphiteMining
    "Ar": "MGA",
    # --- MKD (Macedonian Denar) ---
    # #NorthMacedonia #Balkans #Ferroalloys #ZincLead #AgriExports #EuroPegged
    "ден": "MKD",
    # --- MMK (Myanmar Kyat) ---
    # #Myanmar #SoutheastAsia #JadePreciousStones #RareEarthElements #NaturalGas #Teak
    "K": "MMK",
    # --- MNT (Mongolian Tögrög) ---
    # #Mongolia #EastAsia #OyuTolgoiCopper #CokingCoal #Gold #Cashmere
    "₮": "MNT",
    # --- MOP (Macanese Pataca) ---
    # #Macau #EastAsia #GamingTourism #HKDPegged #ServiceEconomy
    "MOP$": "MOP",
    # --- MRU (Mauritanian Ouguiya) ---
    # #Mauritania #WestAfrica #IronOreSNIM #PelagicFisheries #OffshoreGasGTA #GoldMining
    "UM": "MRU",
    # --- MUR (Mauritian Rupee) ---
    # #Mauritius #IndianOcean #SugarRefining #Textiles #FreeportLogistics #FinancialHub
    "₨_MUR": "MUR",
    # --- MVR (Maldivian Rufiyaa) ---
    # #Maldives #IndianOcean #SkipjackTuna #TourismEconomy #AgriFuelImports
    "Rf": "MVR",
    # --- MWK (Malawian Kwacha) ---
    # #Malawi #SouthernAfrica #BurleyTobacco #TeaExport #SugarCane #UraniumReserves
    "MK": "MWK",
    # --- MXN (Mexican Peso) ---
    # #Mexico #NorthAmerica #SilverMiningGiant #MayaCrude #ManufacturingUSMCA #Automotive
    "Mex$": "MXN",
    # --- MYR (Malaysian Ringgit) ---
    # #Malaysia #SoutheastAsia #PalmOilPowerhouse #PetronasLNG #Electronics #Tin
    "RM": "MYR",
    # --- MZN (Mozambican Metical) ---
    # #Mozambique #SouthernAfrica #LNGProjectsCaboDelgado #CokingCoal #Aluminium #TitaniumSands
    "MT": "MZN",
    # --- NAD (Namibian Dollar) ---
    # #Namibia #SouthernAfrica #UraniumProducer #Diamonds #Zinc #ZARPpegged
    "N$": "NAD",
    # --- NGN (Nigerian Naira) ---
    # #Nigeria #WestAfrica #OPEC #BonnyLight #LNGExporter #DangoteRefinery #PetroEconomy
    "₦": "NGN",
    # --- NIO (Nicaraguan Córdoba) ---
    # #Nicaragua #CentralAmerica #GoldMining #Coffee #BeefExports #SugarCane
    "C$_NIO": "NIO",
    # --- NOK (Norwegian Krone) ---
    # #Norway #Nordic #PetroleumExporter #Equinor #NaturalGasToEurope #SalmonExport #HydroAluminium
    "kr_NOK": "NOK",
    # --- NPR (Nepalese Rupee) ---
    # #Nepal #SouthAsia #HydropowerPotential #INRPegged #AgriImports #RemittanceEconomy
    "रू": "NPR",
    # --- NZD (New Zealand Dollar) ---
    # #NewZealand #CommodityCurrency #DairyPowerhouse #Fonterra #MeatExport #Timber #G10
    "NZ$": "NZD",
    # --- OMR (Omani Rial) ---
    # #Oman #MiddleEast #OmanCrudeBenchmark #LNGExporter #SoharPort #USDPegged
    "ر.ع.": "OMR",
    # --- PAB (Panamanian Balboa) ---
    # #Panama #CentralAmerica #PanamaCanal #CopperMining #FreeTradeZone #USDPegged
    "B/.": "PAB",
    # --- PEN (Peruvian Sol) ---
    # #Peru #SouthAmerica #CopperTitan #SilverProducer #Zinc #Fishmeal #Gold
    "S/": "PEN",
    # --- PGK (Papua New Guinean Kina) ---
    # #PapuaNewGuinea #Oceania #LNGExport #GoldMining #CopperPorgera #PalmOil #Coffee
    "K_PGK": "PGK",
    # --- PHP (Philippine Peso) ---
    # #Philippines #SoutheastAsia #NickelOreExporter #CopperGold #Semiconductors #CoconutOil
    "₱": "PHP",
    # --- PKR (Pakistani Rupee) ---
    # #Pakistan #SouthAsia #TextileCotton #BasmatiRice #LNGBuyer #AgriEconomy
    "Rs_PKR": "PKR",
    # --- PLN (Polish Złoty) ---
    # #Poland #CentralEurope #KGHMCopper #SilverProducer #ThermalCoal #Manufacturing
    "zł": "PLN",
    # --- PYG (Paraguayan Guaraní) ---
    # #Paraguay #SouthAmerica #SoybeanExport #ItaipuHydropower #BeefTrade #ZeroDecimal
    "₲": "PYG",
    # --- QAR (Qatari Riyal) ---
    # #Qatar #MiddleEast #LNGSuperpower #NorthField #HeliumProducer #USDPegged
    "ر.ق": "QAR",
    # --- RON (Romanian Leu) ---
    # #Romania #EasternEurope #BlackSeaGas #CrudeRefining #GrainExporter #Automotive
    "lei": "RON",
    # --- RSD (Serbian Dinar) ---
    # #Serbia #Balkans #BorCopperMine #LithiumJadar #AgriGrains #Manufacturing
    "дин.": "RSD",
    # --- RUB (Russian Ruble) ---
    # #Russia #Eurasia #UralsCrude #NaturalGas #WheatExporter #NickelNorilsk #Palladium #Fertilisers
    "₽": "RUB",
    # --- RWF (Rwandan Franc) ---
    # #Rwanda #EastAfrica #ColtanMining #TinTungsten #CoffeeTea #ZeroDecimal
    "FRw": "RWF",
    # --- SAR (Saudi Riyal) ---
    # #SaudiArabia #MiddleEast #OPECLeader #SaudiAramco #ArabLight #Petrochemicals #USDPegged
    "ر.س": "SAR",
    # --- SBD (Solomon Islands Dollar) ---
    # #SolomonIslands #Oceania #TimberExport #RoughPalmOil #Fisheries #BauxiteNickel
    "SI$": "SBD",
    # --- SCR (Seychellois Rupee) ---
    # #Seychelles #IndianOcean #TunaCanningHub #MaritimeLogistics #TourismEconomy
    "₨_SCR": "SCR",
    # --- SDG (Sudanese Pound) ---
    # #Sudan #NorthAfrica #NileBlendCrude #GoldMining #GumArabicMonopoly #AgriPotential
    "ج.س.": "SDG",
    # --- SEK (Swedish Krona) ---
    # #Sweden #Nordic #IronOreLKAB #ForestryTimber #SteelManufacturing #HydroPower #G10
    "kr_SEK": "SEK",
    # --- SGD (Singapore Dollar) ---
    # #Singapore #SoutheastAsia #OilRefiningHub #BunkeringCapital #CommodityTradingHub #GoldStorage
    "S$": "SGD",
    # --- SHP (Saint Helena Pound) ---
    # #SaintHelena #SouthAtlantic #Fisheries #Coffee #GBPPegged #RemoteEconomy
    "£_SHP": "SHP",
    # --- SLE (Sierra Leonean Leone) ---
    # #SierraLeone #WestAfrica #DiamondMining #RutileTitanium #Bauxite #IronOre
    "Le": "SLE",
    # --- SOS (Somali Shilling) ---
    # #Somalia #HornOfAfrica #LivestockExport #Frankincense #FisheriesReserves #RedSeaCoast
    "Sh": "SOS",
    # --- SRD (Surinamese Dollar) ---
    # #Suriname #SouthAmerica #OffshoreOilBasin #GoldMining #AluminaHistory #Forestry
    "SR$": "SRD",
    # --- SSP (South Sudanese Pound) ---
    # #SouthSudan #EastAfrica #DarBlendCrude #OilPipelinesToPortSudan #AgriPotential
    "£_SSP": "SSP",
    # --- STN (São Tomé and Príncipe Dobra) ---
    # #SaoTomePrincipe #CentralAfrica #CocoaExport #OffshoreOilExploration #JDZ
    "Db": "STN",
    # --- SVC (Salvadoran Colón) ---
    # #ElSalvador #CentralAmerica #CoffeeExport #SugarCane #GeothermalPower #USDDollarised
    "₡_SVC": "SVC",
    # --- SYP (Syrian Pound) ---
    # #Syria #MiddleEast #LightCrudeFields #PhosphateRock #AgriCotton #Wheat
    "£S": "SYP",
    # --- SZL (Eswatini Lilangeni) ---
    # #Eswatini #SouthernAfrica #SugarCaneGiant #WoodPulp #ZARPpegged #AgriTrade
    "E": "SZL",
    # --- THB (Thai Baht) ---
    # #Thailand #SoutheastAsia #NaturalRubberLeader #RiceExporter #SugarCane #Petrochemicals
    "฿": "THB",
    # --- TJS (Tajikistani Somoni) ---
    # #Tajikistan #CentralAsia #AluminiumTalco #HydropowerRogun #Cotton #GoldAntimony
    "ЅМ": "TJS",
    # --- TMT (Turkmenistan Manat) ---
    # #Turkmenistan #CentralAsia #GalkynyshGasField #PipelineToChina #CottonHarvest #Petrochemicals
    "m": "TMT",
    # --- TND (Tunisian Dinar) ---
    # #Tunisia #NorthAfrica #PhosphateMining #OliveOilExporter #CrudeHydrocarbons #Dates
    "د.ت": "TND",
    # --- TOP (Tongan Paʻanga) ---
    # #Tonga #Pacific #VanillaExport #RootCrops #Fisheries #RemittanceEconomy
    "T$": "TOP",
    # --- TRY (Turkish Lira) ---
    # #Turkey #Eurasia #EnergyTransitHub #BoronMineralsMonopoly #SteelExporter #Agriculture
    "₺": "TRY",
    # --- TTD (Trinidad and Tobago Dollar) ---
    # #TrinidadTobago #Caribbean #LNGExporter #AmmoniaMethanolCapital #PetroChemicals
    "TT$": "TTD",
    # --- TWD (New Taiwan Dollar) ---
    # #Taiwan #EastAsia #SemiconductorFoundry #HighTechConsumer #EnergyImporter #LNGBuyer
    "NT$": "TWD",
    # --- TZS (Tanzanian Shilling) ---
    # #Tanzania #EastAfrica #GoldMining #GraphiteReserves #LNGProjectLindi #Cashews #Coffee
    "TSh": "TZS",
    # --- UAH (Ukrainian Hryvnia) ---
    # #Ukraine #EasternEurope #BreadbasketOfEurope #WheatCorn #SunflowerOil #IronOrePellets
    "₴": "UAH",
    # --- UGX (Ugandan Shilling) ---
    # #Uganda #EastAfrica #RobustaCoffee #TilengaKingfisherOil #EACOP #ZeroDecimal
    "USh": "UGX",
    # --- USD (United States Dollar) ---
    # #UnitedStates #GlobalReserveCurrency #WTIReference #HenryHub #CBOTGrains #Petrodollar
    "$": "USD",
    # --- UYU (Uruguayan Peso) ---
    # #Uruguay #SouthAmerica #BeefExporter #CellulosePulp #Soybeans #WindHydropower
    "$U": "UYU",
    # --- UZS (Uzbekistani Soʻm) ---
    # #Uzbekistan #CentralAsia #MuruntauGoldMine #UraniumProducer #CottonExport #NaturalGas
    "soʻm": "UZS",
    # --- VES (Venezuelan Bolívar Soberano) ---
    # #Venezuela #SouthAmerica #OPECFounding #OrinocoHeavyCrude #LargestOilReserves #Refineries
    "Bs._VES": "VES",
    # --- VND (Vietnamese Đồng) ---
    # #Vietnam #SoutheastAsia #RobustaCoffeeLeader #RiceExporter #Electronics #ZeroDecimal
    "₫": "VND",
    # --- VUV (Vanuatu Vatu) ---
    # #Vanuatu #Oceania #CopraExport #BeefProduction #MaritimeLogistics #ZeroDecimal
    "VT": "VUV",
    # --- WST (Samoan Tālā) ---
    # #Samoa #Pacific #FishExports #TaroRoot #CoconutOil #TourismEconomy
    "WS$": "WST",
    # --- XAF (Central African CFA Franc) ---
    # #CEMAC #CentralAfrica #CrudeOilExporters #TimberForestry #EuroPegged #ZeroDecimal
    "FCFA": "XAF",
    # --- XCD (East Caribbean Dollar) ---
    # #OECS #Caribbean #CurrencyUnion #BananasNutmeg #TourismLogistics #USDPegged
    "EC$": "XCD",
    # --- XCG (Caribbean Guilder) ---
    # #CuracaoSintMaarten #Caribbean #RefiningHistory #BunkeringLogistics #Tourism #USDPegged
    "Cg.": "XCG",
    # --- XOF (West African CFA Franc) ---
    # #WAEMU #WestAfrica #CocoaPowerhouse #GoldMining #BCEAO #EuroPegged #ZeroDecimal
    "CFA": "XOF",
    # --- XPF (CFP Franc) ---
    # #FrenchPacific #NewCaledonia #NickelMiningPowerhouse #FrenchPolynesia #EuroPegged #ZeroDecimal
    "CFP": "XPF",
    # --- YER (Yemeni Rial) ---
    # #Yemen #MiddleEast #MaribLNG #RasIsaTerminal #BabElMandebTransit #CoffeeMocha
    "﷼_YER": "YER",
    # --- ZAR (South African Rand) ---
    # #SouthAfrica #SouthernAfrica #PGMSuperpower #PlatinumPalladium #GoldMining #CoalExporter
    "R": "ZAR",
    # --- ZMW (Zambian Kwacha) ---
    # #Zambia #SouthernAfrica #CopperbeltPowerhouse #CobaltProducer #Emeralds #Hydropower
    "ZK": "ZMW",
    # --- ZWG (Zimbabwe Gold) ---
    # #Zimbabwe #SouthernAfrica #GoldBacked #PGMReserves #LithiumMining #TobaccoExport
    "ZiG": "ZWG",
}

# ---------------------------------------------------------------------------
# Minor Unit Sets & Exponent Mappings
# ---------------------------------------------------------------------------

_ZERO_MINOR_UNIT: Final[frozenset[str]] = frozenset({
    # Immutable set of ISO 4217 currency codes that have zero decimal places in active financial settlement and commercial circulation.
    # --- BIF (Burundian Franc) ---
    # #Burundi #EastAfrica #Franc #ZeroDecimal #GreatLakesRegion #CoffeeTeaExport
    "BIF",
    # --- CLP (Chilean Peso) ---
    # #Chile #CopperCurrency #LatinAmerica #CommodityCurrency #ZeroDecimal #CentavoPhasedOut
    "CLP",
    # --- DJF (Djiboutian Franc) ---
    # #Djibouti #HornOfAfrica #PortLogistics #USDPegged #ZeroDecimal #RedSeaTrade
    "DJF",
    # --- GNF (Guinean Franc) ---
    # #Guinea #WestAfrica #BauxiteMining #AluminiumSupply #ZeroDecimal #MiningEconomy
    "GNF",
    # --- ISK (Icelandic Króna) ---
    # #Iceland #Nordic #Fisheries #GeothermalEnergy #ZeroDecimal #AluminiumSmelting
    "ISK",
    # --- JPY (Japanese Yen) ---
    # #Japan #MajorFX #G10 #ReserveCurrency #ZeroDecimal #CommodityImporter #SafeHaven
    "JPY",
    # --- KMF (Comorian Franc) ---
    # #Comoros #IndianOcean #VanillaProducer #FrancZone #ZeroDecimal #IslandEconomy
    "KMF",
    # --- KRW (South Korean Won) ---
    # #SouthKorea #EastAsia #IndustrialDemand #RefiningHub #TechExports #ZeroDecimal
    "KRW",
    # --- PYG (Paraguayan Guaraní) ---
    # #Paraguay #SouthAmerica #SoybeanExport #HydroelectricPower #BeefTrade #ZeroDecimal
    "PYG",
    # --- RWF (Rwandan Franc) ---
    # #Rwanda #EastAfrica #ColtanMining #TinTungsten #ZeroDecimal #AgriCommodities
    "RWF",
    # --- UGX (Ugandan Shilling) ---
    # #Uganda #EastAfrica #RobustaCoffee #CrudeOilExplorer #ZeroDecimal #AgriTrade
    "UGX",
    # --- VND (Vietnamese Đồng) ---
    # #Vietnam #SoutheastAsia #CoffeeExporter #ManufacturingHub #HighDenomination #ZeroDecimal
    "VND",
    # --- VUV (Vanuatu Vatu) ---
    # #Vanuatu #Oceania #PacificIslands #CopraExport #ZeroDecimal #MaritimeEconomy
    "VUV",
    # --- XAF (Central African CFA Franc) ---
    # #CEMAC #CentralAfrica #CrudeOilExporters #BEAC #EuroPegged #ZeroDecimal
    "XAF",
    # --- XOF (West African CFA Franc) ---
    # #WAEMU #WestAfrica #CocoaCapital #GoldMining #BCEAO #ZeroDecimal
    "XOF",
    # --- XPF (CFP Franc) ---
    # #FrenchPacific #NewCaledonia #NickelMining #FrenchPolynesia #EuroPegged #ZeroDecimal
    "XPF",
})

_THREE_MINOR_UNIT: Final[frozenset[str]] = frozenset({
    # Immutable set of ISO 4217 currency codes that have 3 decimal places in active financial settlement and commercial circulation (e.g. 1 BHD = 1,000 fils).
    # --- BHD (Bahraini Dinar) ---
    # #Bahrain #GulfStates #ThreeDecimal #AluminiumSmelting #OilRefining #PetroEconomy #USDPegged #Fils
    "BHD",
    # --- IQD (Iraqi Dinar) ---
    # #Iraq #MiddleEast #OPEC #ThreeDecimal #BasrahCrude #MegaOilFields #PetroState #Fils
    "IQD",
    # --- JOD (Jordanian Dinar) ---
    # #Jordan #MiddleEast #ThreeDecimal #PotashMining #PhosphateRock #Fertilisers #USDPegged #Fils
    "JOD",
    # --- KWD (Kuwaiti Dinar) ---
    # #Kuwait #MiddleEast #OPEC #ThreeDecimal #CrudeOilExporter #SovereignWealth #HighestValuedCurrency #Fils
    "KWD",
    # --- LYD (Libyan Dinar) ---
    # #Libya #NorthAfrica #OPEC #ThreeDecimal #LightSweetCrude #MellitahGas #PetroEconomy #Dirham
    "LYD",
    # --- OMR (Omani Rial) ---
    # #Oman #MiddleEast #ThreeDecimal #OmanCrudeBenchmark #LNGExporter #SoharPort #USDPegged #Baisa
    "OMR",
    # --- TND (Tunisian Dinar) ---
    # #Tunisia #NorthAfrica #ThreeDecimal #PhosphateMining #OliveOilExporter #CrudeHydrocarbons #Millime
    "TND",
})

_NON_DECIMAL_MINOR_UNIT: Final[frozenset[str]] = frozenset({
    # Immutable set of ISO 4217 currency codes with non-decimal base-5 subdivisions
    # (divided into 5 base units rather than 100).
    # --- MGA (Malagasy Ariary) ---
    # #Madagascar #IndianOcean #NonDecimal #Base5 #Iraimbilanja #VanillaMonopoly #NickelCobalt #GraphiteMining
    "MGA",
    # --- MRU (Mauritanian Ouguiya) ---
    # #Mauritania #WestAfrica #NonDecimal #Base5 #Khoums #IronOreSNIM #PelagicFisheries #OffshoreGasGTA
    "MRU",
})

# Constructs a lookup table for the number of decimal places (minor units) for each currency code.
def _build_minor_unit_exponents() -> dict[str, int]: # Function to build a dictionary mapping currency codes to their minor unit exponents (number of decimal places).
    exponents: dict[str, int] = {} # Initialises empty dictionary to hold the exponents for each currency code.
    for code in CURRENCY_CODES: # Iterates over each currency code in the predefined set of currency codes.
        if code in _ZERO_MINOR_UNIT or code in _NON_DECIMAL_MINOR_UNIT: # If the currency code is in the set of zero minor units or non-decimal minor units, it assigns an exponent of 0 (indicating no decimal places).
            exponents[code] = 0 # Assigns an exponent of 0 for currencies with no decimal places.
        elif code in _THREE_MINOR_UNIT: # If the currency code is in the set of three minor units, assigns an exponent of 3 (indicating three decimal places).
            exponents[code] = 3 # Assigns an exponent of 3 for currencies with three decimal places.
        else: # For all other currencies, assigns an exponent of 2 (indicating two decimal places).
            exponents[code] = 2 # Assigns an exponent of 2 for currencies with two decimal places.
    return exponents # Returns the constructed dictionary mapping currency codes to their respective minor unit exponents.

CURRENCY_MINOR_UNIT_EXPONENT: Final[dict[str, int]] = _build_minor_unit_exponents() # Builds final dictionary mapping currency codes to their minor unit exponents using the helper function.

# Master downstream metadata registry for settlement and cash rounding engines populated from exponents
CURRENCY_MINOR_DIGITS: dict[str, int] = dict(CURRENCY_MINOR_UNIT_EXPONENT)

# Global lock to ensure thread-safe dynamic registration of custom units and currencies
_REGISTRY_LOCK = threading.RLock()
_REGISTRY_FROZEN: bool = False

# Constructs Pint to treat every ISO 4217 currency as a completely independent, non-interchangeable dimension of measurement.
CURRENCY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = tuple(f"{code} = [currency_{code}]" for code in CURRENCY_CODES)

# Constructs Pint definitions for the minor units of each ISO 4217 currency, allowing for precise representation of fractional currency amounts in calculations.
CURRENCY_MINOR_UNIT_DEFINITIONS: Final[tuple[str, ...]] = (
    # --- AED (United Arab Emirates Dirham) ---
    # #UAE #MiddleEast #Petrodollar #CrudeOil #RefinedProducts #USDPegged #DubaiGold
    "AED_fils = 0.01 * AED = fils_AED = AED_minor",
    # --- AFN (Afghan Afghani) ---
    # #Afghanistan #CentralAsia #MiningPotential #LithiumReserves #PreciousGems
    "AFN_pul = 0.01 * AFN = pul_AFN = AFN_minor",
    # --- ALL (Albanian Lek) ---
    # #Albania #Balkans #ChromiumMining #Hydroelectric #AgriExports
    "ALL_qindarke = 0.01 * ALL = qindarke_ALL = qintar_ALL = ALL_minor",
    # --- AMD (Armenian Dram) ---
    # #Armenia #Caucasus #CopperMining #Molybdenum #PreciousMetals
    "AMD_luma = 0.01 * AMD = luma_AMD = AMD_minor",
    # --- AOA (Angolan Kwanza) ---
    # #Angola #CentralAfrica #OPEC #CrudeOilExporter #Diamonds #PetroEconomy
    "AOA_centimo = 0.01 * AOA = centimo_AOA = AOA_minor",
    # --- ARS (Argentine Peso) ---
    # #Argentina #SouthAmerica #AgriCommodities #Soybeans #Corn #LithiumTriangle
    "ARS_centavo = 0.01 * ARS = centavo_ARS = ARS_minor",
    # --- AUD (Australian Dollar) ---
    # #Australia #CommodityCurrency #IronOre #LNGExporter #Coal #Gold #G10
    "AUD_cent = 0.01 * AUD = cent_AUD = c_AUD = AUD_minor",
    # --- AWG (Aruban Florin) ---
    # #Aruba #Caribbean #TourismEconomy #USDPegged #RefineryHistory
    "AWG_cent = 0.01 * AWG = cent_AWG = AWG_minor",
    # --- AZN (Azerbaijani Manat) ---
    # #Azerbaijan #CaspianSea #CrudeOil #NaturalGas #SouthernGasCorridor
    "AZN_qepik = 0.01 * AZN = qepik_AZN = gapik_AZN = AZN_minor",
    # --- BAM (Bosnia and Herzegovina Convertible Mark) ---
    # #BosniaHerzegovina #Balkans #BauxiteAluminium #Hydropower #EuroPegged
    "BAM_fening = 0.01 * BAM = fening_BAM = BAM_minor",
    # --- BBD (Barbadian Dollar) ---
    # #Barbados #Caribbean #SugarCane #Tourism #USDPegged
    "BBD_cent = 0.01 * BBD = cent_BBD = BBD_minor",
    # --- BDT (Bangladeshi Taka) ---
    # #Bangladesh #SouthAsia #Textiles #JuteExport #AgriImports #LNGImporter
    "BDT_poisha = 0.01 * BDT = poisha_BDT = paisa_BDT = BDT_minor",
    # --- BHD (Bahraini Dinar - 3 Decimals) ---
    # #Bahrain #GulfStates #ThreeDecimal #AluminiumSmelting #OilRefining #PetroEconomy #USDPegged #Fils
    "BHD_fils = 0.001 * BHD = fils_BHD = BHD_minor",
    # --- BIF (Burundian Franc - Nominal/Historical) ---
    # #Burundi #EastAfrica #CoffeeTeaExport #NickelReserves #ZeroDecimal
    "BIF_centime = 0.01 * BIF = centime_BIF = BIF_minor",
    # --- BMD (Bermudian Dollar) ---
    # #Bermuda #OffshoreFinance #InsuranceReinsurance #USDPegged #MaritimeRegistry
    "BMD_cent = 0.01 * BMD = cent_BMD = BMD_minor",
    # --- BND (Brunei Dollar) ---
    # #Brunei #SoutheastAsia #PetroleumExporter #LNGExporter #SGDPegged
    "BND_sen = 0.01 * BND = sen_BND = cent_BND = BND_minor",
    # --- BOB (Bolivian Boliviano) ---
    # #Bolivia #SouthAmerica #LithiumReserves #NaturalGas #TinMining #Zinc
    "BOB_centavo = 0.01 * BOB = centavo_BOB = BOB_minor",
    # --- BRL (Brazilian Real) ---
    # #Brazil #SouthAmerica #CommodityPowerhouse #Soybeans #IronOre #CrudeOil #Sugar #Coffee
    "BRL_centavo = 0.01 * BRL = centavo_BRL = BRL_minor",
    # --- BSD (Bahamian Dollar) ---
    # #Bahamas #Caribbean #Tourism #OffshoreBanking #USDPegged #Aragonite
    "BSD_cent = 0.01 * BSD = cent_BSD = BSD_minor",
    # --- BTN (Bhutanese Ngultrum) ---
    # #Bhutan #SouthAsia #Hydroelectricity #INRPegged #Forestry #Himalayas
    "BTN_chhertum = 0.01 * BTN = chhertum_BTN = chetrum_BTN = BTN_minor",
    # --- BWP (Botswana Pula) ---
    # #Botswana #SouthernAfrica #DiamondMining #CopperNickel #BeefExport
    "BWP_thebe = 0.01 * BWP = thebe_BWP = BWP_minor",
    # --- BYN (Belarusian Ruble) ---
    # #Belarus #EasternEurope #PotashFertiliser #HeavyMachinery #OilRefining
    "BYN_kopeck = 0.01 * BYN = kopeck_BYN = kapiejka_BYN = BYN_minor",
    # --- BZD (Belize Dollar) ---
    # #Belize #CentralAmerica #SugarCane #CitrusExports #CrudePetroleum #USDPegged
    "BZD_cent = 0.01 * BZD = cent_BZD = BZD_minor",
    # --- CAD (Canadian Dollar) ---
    # #Canada #CommodityCurrency #CrudeOil #OilSands #NaturalGas #Wheat #Potash #G10
    "CAD_cent = 0.01 * CAD = cent_CAD = c_CAD = CAD_minor",
    # --- CDF (Congolese Franc) ---
    # #DRC #CentralAfrica #CobaltCapital #CopperMining #Coltan #CriticalMinerals
    "CDF_centime = 0.01 * CDF = centime_CDF = CDF_minor",
    # --- CHF (Swiss Franc) ---
    # #Switzerland #SafeHaven #GoldRefiningHub #CommodityTradingHub #G10
    "CHF_rappen = 0.01 * CHF = rappen_CHF = centime_CHF = centesimo_CHF = rap_CHF = CHF_minor",
    # --- CLP (Chilean Peso - Nominal/Historical) ---
    # #Chile #CopperPowerhouse #LithiumProducer #LatinAmerica #ZeroDecimal
    "CLP_centavo = 0.01 * CLP = centavo_CLP = CLP_minor",
    # --- CNY (Chinese Yuan Renminbi - Multi-Tier) ---
    # #China #MajorEconomy #LargestCommodityImporter #Petroyuan #MetalsConsumer
    "CNY_jiao = 0.1 * CNY = jiao_CNY = mao_CNY",
    "CNY_fen = 0.01 * CNY = fen_CNY = CNY_minor",
    # --- COP (Colombian Peso) ---
    # #Colombia #SouthAmerica #CrudeOil #ThermalCoal #ArabicaCoffee #Nickel
    "COP_centavo = 0.01 * COP = centavo_COP = COP_minor",
    # --- CRC (Costa Rican Colón) ---
    # #CostaRica #CentralAmerica #Pineapples #Bananas #Coffee #RenewableEnergy
    "CRC_centimo = 0.01 * CRC = centimo_CRC = CRC_minor",
    # --- CUP (Cuban Peso) ---
    # #Cuba #Caribbean #NickelMining #SugarCane #CobaltReserves #Tobacco
    "CUP_centavo = 0.01 * CUP = centavo_CUP = CUP_minor",
    # --- CVE (Cape Verdean Escudo) ---
    # #CapeVerde #WestAfrica #MaritimeLogistics #Fisheries #EuroPegged
    "CVE_centavo = 0.01 * CVE = centavo_CVE = CVE_minor",
    # --- CZK (Czech Koruna - Nominal/Historical) ---
    # #Czechia #CentralEurope #IndustrialManufacturing #AutomotiveSupply #CoalMining
    "CZK_haler = 0.01 * CZK = haler_CZK = heller_CZK = CZK_minor",
    # --- DJF (Djiboutian Franc - Nominal/Historical) ---
    # #Djibouti #HornOfAfrica #PortLogistics #RedSeaTrade #USDPegged #ZeroDecimal
    "DJF_centime = 0.01 * DJF = centime_DJF = DJF_minor",

    # --- DKK (Danish Krone) ---
    # #Denmark #Nordic #OffshoreWind #GreenEnergy #EuroPegged #MaritimeShipping
    "DKK_ore = 0.01 * DKK = ore_DKK = DKK_minor",

    # --- DOP (Dominican Peso) ---
    # #DominicanRepublic #Caribbean #Ferronickel #GoldMining #SugarCane #Cacao
    "DOP_centavo = 0.01 * DOP = centavo_DOP = DOP_minor",

    # --- DZD (Algerian Dinar) ---
    # #Algeria #NorthAfrica #OPEC #NaturalGasExporter #CrudeOil #PipelinesToEurope
    "DZD_santeem = 0.01 * DZD = santeem_DZD = centime_DZD = DZD_minor",

    # --- EGP (Egyptian Pound - Multi-Tier) ---
    # #Egypt #NorthAfrica #SuezCanal #LNGExporter #Fertilisers #WheatImporter
    "EGP_piastre = 0.01 * EGP = piastre_EGP = qirsh_EGP = EGP_minor",
    "EGP_millieme = 0.001 * EGP = millieme_EGP",

    # --- ERN (Eritrean Nakfa) ---
    # #Eritrea #HornOfAfrica #PotashMining #GoldCopper #RedSeaCoast
    "ERN_cent = 0.01 * ERN = cent_ERN = ERN_minor",

    # --- ETB (Ethiopian Birr) ---
    # #Ethiopia #EastAfrica #ArabicaOrigin #HydroelectricPower #Oilseeds #Floriculture
    "ETB_santim = 0.01 * ETB = santim_ETB = ETB_minor",

    # --- EUR (Euro) ---
    # #Eurozone #MajorFX #ReserveCurrency #CarbonAllowances #EUA #IndustrialHub
    "EUR_cent = 0.01 * EUR = euro_cent = cent_EUR = cents_EUR = c_EUR = EUR_minor",

    # --- FJD (Fijian Dollar) ---
    # #Fiji #Pacific #SugarExport #GoldMining #Fisheries #Tourism
    "FJD_cent = 0.01 * FJD = cent_FJD = FJD_minor",

    # --- FKP (Falkland Islands Pound) ---
    # #FalklandIslands #SouthAtlantic #Fisheries #SquidExport #OffshoreHydrocarbons #GBPPegged
    "FKP_penny = 0.01 * FKP = penny_FKP = pence_FKP = FKP_minor",

    # --- GBP (British Pound Sterling) ---
    # #UnitedKingdom #MajorFX #G10 #BrentBenchmark #LME #NBPGas #UKA
    "GBP_penny = 0.01 * GBP = GBp = penny_GBP = pence_GBP = p_GBP = pence = GBP_minor = GBP_pence",

    # --- GEL (Georgian Lari) ---
    # #Georgia #Caucasus #TransitCorridor #CopperOres #Ferroalloys #WineExport
    "GEL_tetri = 0.01 * GEL = tetri_GEL = GEL_minor",

    # --- GHS (Ghanaian Cedi) ---
    # #Ghana #WestAfrica #GoldProducer #CocoaExporter #CrudeOil #Bauxite
    "GHS_pesewa = 0.01 * GHS = pesewa_GHS = GHS_minor",

    # --- GIP (Gibraltar Pound) ---
    # #Gibraltar #Mediterranean #BunkeringHub #MaritimeLogistics #GBPPegged
    "GIP_penny = 0.01 * GIP = penny_GIP = pence_GIP = GIP_minor",

    # --- GMD (Gambian Dalasi) ---
    # #Gambia #WestAfrica #GroundnutExport #Fisheries #AgriTrade
    "GMD_butut = 0.01 * GMD = butut_GMD = GMD_minor",

    # --- GNF (Guinean Franc - Nominal/Historical) ---
    # #Guinea #WestAfrica #BauxiteGiant #AluminiumSupply #IronOreSimandou #ZeroDecimal
    "GNF_centime = 0.01 * GNF = centime_GNF = GNF_minor",

    # --- GTQ (Guatemalan Quetzal) ---
    # #Guatemala #CentralAmerica #Cardamom #Coffee #SugarCane #NickelMining
    "GTQ_centavo = 0.01 * GTQ = centavo_GTQ = GTQ_minor",

    # --- GYD (Guyanese Dollar) ---
    # #Guyana #SouthAmerica #OffshoreOilBoom #StabroekBlock #GoldMining #Bauxite
    "GYD_cent = 0.01 * GYD = cent_GYD = GYD_minor",

    # --- HKD (Hong Kong Dollar - Multi-Tier) ---
    # #HongKong #EastAsia #CommodityTrading #GoldTrading #USDPeggedCurrencyBoard
    "HKD_hou = 0.1 * HKD = hou_HKD = ho_HKD",
    "HKD_cent = 0.01 * HKD = cent_HKD = c_HKD = HKD_minor",

    # --- HNL (Honduran Lempira) ---
    # #Honduras #CentralAmerica #CoffeeExport #Bananas #PalmOil #ZincLead
    "HNL_centavo = 0.01 * HNL = centavo_HNL = HNL_minor",

    # --- HTG (Haitian Gourde) ---
    # #Haiti #Caribbean #EssentialOils #Vetiver #AgriImports
    "HTG_centime = 0.01 * HTG = centime_HTG = HTG_minor",

    # --- HUF (Hungarian Forint - Nominal/Historical) ---
    # #Hungary #CentralEurope #Manufacturing #AutomotiveSupply #NaturalGasTransit
    "HUF_filler = 0.01 * HUF = filler_HUF = HUF_minor",

    # --- IDR (Indonesian Rupiah - Nominal/Historical) ---
    # #Indonesia #SoutheastAsia #NickelPowerhouse #ThermalCoal #PalmOil #Tin
    "IDR_sen = 0.01 * IDR = sen_IDR = IDR_minor",

    # --- ILS (Israeli New Shekel) ---
    # #Israel #MiddleEast #OffshoreGas #LeviathanField #PolishedDiamonds #HighTech
    "ILS_agora = 0.01 * ILS = agora_ILS = ILS_minor",

    # --- INR (Indian Rupee) ---
    # #India #SouthAsia #MajorCommodityConsumer #CrudeOilImporter #GoldConsumer #CoalBuyer
    "INR_paisa = 0.01 * INR = paisa_INR = INR_minor",

    # --- IQD (Iraqi Dinar - 3 Decimals) ---
    # #Iraq #MiddleEast #OPEC #ThreeDecimal #BasrahCrude #MegaOilFields #PetroState #Fils
    "IQD_fils = 0.001 * IQD = fils_IQD = IQD_minor",

    # --- IRR (Iranian Rial) ---
    # #Iran #MiddleEast #CrudeOilReserves #NaturalGasGiant #Petrochemicals #Copper
    "IRR_dinar = 0.01 * IRR = dinar_IRR = IRR_minor",
    "IRR_toman = 10 * IRR = toman_IRR",

    # --- ISK (Icelandic Króna - Nominal/Historical) ---
    # #Iceland #Nordic #AluminiumSmelting #Geothermal #Fisheries #ZeroDecimal
    "ISK_eyrir = 0.01 * ISK = eyrir_ISK = ISK_minor",

    # --- JMD (Jamaican Dollar) ---
    # #Jamaica #Caribbean #BauxiteAlumina #Coffee #SugarCane #Tourism
    "JMD_cent = 0.01 * JMD = cent_JMD = JMD_minor",

    # --- JOD (Jordanian Dinar - Multi-Tier / 3 Decimals) ---
    # #Jordan #MiddleEast #ThreeDecimal #PotashMining #PhosphateRock #Fertilisers #USDPegged #Fils
    "JOD_piastre = 0.01 * JOD = piastre_JOD = qirsh_JOD",
    "JOD_fils = 0.001 * JOD = fils_JOD = JOD_minor",

    # --- JPY (Japanese Yen - Nominal/Historical Subunits) ---
    # #Japan #MajorFX #G10 #ReserveCurrency #ZeroDecimal #MajorCommodityImporter #LNGImporter
    "JPY_sen = 0.01 * JPY = sen_JPY = JPY_minor",
    "JPY_rin = 0.001 * JPY = rin_JPY",

    # --- KES (Kenyan Shilling) ---
    # #Kenya #EastAfrica #TeaExporter #CutFlowers #GeothermalPower #RegionalHub
    "KES_cent = 0.01 * KES = cent_KES = KES_minor",

    # --- KGS (Kyrgyzstani Som) ---
    # #Kyrgyzstan #CentralAsia #KumtorGoldMine #Hydropower #RareEarths
    "KGS_tyiyn = 0.01 * KGS = tyiyn_KGS = KGS_minor",

    # --- KHR (Cambodian Riel - Nominal/Historical) ---
    # #Cambodia #SoutheastAsia #RiceExport #Rubber #Garments #AgriCommodities
    "KHR_kak = 0.1 * KHR = kak_KHR",
    "KHR_sen = 0.01 * KHR = sen_KHR = KHR_minor",

    # --- KMF (Comorian Franc - Nominal/Historical) ---
    # #Comoros #IndianOcean #VanillaExporter #YlangYlang #Cloves #ZeroDecimal
    "KMF_centime = 0.01 * KMF = centime_KMF = KMF_minor",

    # --- KPW (North Korean Won) ---
    # #NorthKorea #EastAsia #AnthraciteCoal #IronOre #Magnesite #MiningResources
    "KPW_chon = 0.01 * KPW = chon_KPW = KPW_minor",

    # --- KRW (South Korean Won - Accounting Subunit) ---
    # #SouthKorea #EastAsia #RefiningPowerhouse #Petrochemicals #LNGImporter #ZeroDecimal
    "KRW_jeon = 0.01 * KRW = jeon_KRW = KRW_minor",

    # --- KWD (Kuwaiti Dinar - 3 Decimals) ---
    # #Kuwait #MiddleEast #OPEC #ThreeDecimal #CrudeOilExporter #SovereignWealth #HighestValuedCurrency #Fils
    "KWD_fils = 0.001 * KWD = fils_KWD = KWD_minor",

    # --- KYD (Cayman Islands Dollar) ---
    # #CaymanIslands #Caribbean #OffshoreFunds #FinancialServices #USDPegged
    "KYD_cent = 0.01 * KYD = cent_KYD = KYD_minor",

    # --- KZT (Kazakhstani Tenge) ---
    # #Kazakhstan #CentralAsia #UraniumGiant #CPCBlend #CrudeOil #Ferroalloys #Copper
    "KZT_tiyn = 0.01 * KZT = tiyn_KZT = KZT_minor",

    # --- LAK (Lao Kip - Nominal/Historical) ---
    # #Laos #SoutheastAsia #HydropowerToThailand #CopperMining #Gold #Potash
    "LAK_att = 0.01 * LAK = att_LAK = LAK_minor",

    # --- LBP (Lebanese Pound) ---
    # #Lebanon #MiddleEast #CommercialTrading #GoldReserves #AgriImports
    "LBP_piastre = 0.01 * LBP = piastre_LBP = qirsh_LBP = LBP_minor",

    # --- LKR (Sri Lankan Rupee) ---
    # #SriLanka #SouthAsia #CeylonTea #Rubber #Cinnamon #Apparel
    "LKR_cent = 0.01 * LKR = cent_LKR = LKR_minor",

    # --- LRD (Liberian Dollar) ---
    # #Liberia #WestAfrica #IronOre #NaturalRubber #Gold #MaritimeFlagRegistry
    "LRD_cent = 0.01 * LRD = cent_LRD = LRD_minor",

    # --- LSL (Lesotho Loti) ---
    # #Lesotho #SouthernAfrica #DiamondMining #WaterExportToSouthAfrica #ZARPpegged
    "LSL_sente = 0.01 * LSL = sente_LSL = lisente_LSL = LSL_minor",

    # --- LYD (Libyan Dinar - 3 Decimals) ---
    # #Libya #NorthAfrica #OPEC #ThreeDecimal #LightSweetCrude #MellitahGas #PetroEconomy #Dirham
    "LYD_dirham = 0.001 * LYD = dirham_LYD = LYD_minor",

    # --- MAD (Moroccan Dirham) ---
    # #Morocco #NorthAfrica #PhosphateCapital #OCPGroup #FertiliserPowerhouse #Automotive
    "MAD_santim = 0.01 * MAD = santim_MAD = centime_MAD = MAD_minor",

    # --- MDL (Moldovan Leu) ---
    # #Moldova #EasternEurope #AgriExports #WineIndustry #Sunflowers #GrainTrade
    "MDL_ban = 0.01 * MDL = ban_MDL = MDL_minor",

    # --- MGA (Malagasy Ariary - Non-Decimal Base 5) ---
    # #Madagascar #IndianOcean #NonDecimal #Base5 #Iraimbilanja #VanillaMonopoly #NickelCobalt #GraphiteMining
    "MGA_iraimbilanja = 0.2 * MGA = iraimbilanja_MGA = MGA_minor",

    # --- MKD (Macedonian Denar - Nominal/Historical) ---
    # #NorthMacedonia #Balkans #Ferroalloys #ZincLead #AgriExports #EuroPegged
    "MKD_deni = 0.01 * MKD = deni_MKD = MKD_minor",

    # --- MMK (Myanmar Kyat - Nominal/Historical) ---
    # #Myanmar #SoutheastAsia #JadePreciousStones #RareEarthElements #NaturalGas #Teak
    "MMK_pya = 0.01 * MMK = pya_MMK = MMK_minor",

    # --- MNT (Mongolian Tögrög) ---
    # #Mongolia #EastAsia #OyuTolgoiCopper #CokingCoal #Gold #Cashmere
    "MNT_mongo = 0.01 * MNT = mongo_MNT = mongu_MNT = MNT_minor",

    # --- MOP (Macanese Pataca - Multi-Tier) ---
    # #Macau #EastAsia #GamingTourism #HKDPegged #ServiceEconomy
    "MOP_ho = 0.1 * MOP = ho_MOP",
    "MOP_avo = 0.01 * MOP = avo_MOP = MOP_minor",

    # --- MRU (Mauritanian Ouguiya - Non-Decimal Base 5) ---
    # #Mauritania #WestAfrica #NonDecimal #Base5 #Khoums #IronOreSNIM #PelagicFisheries #OffshoreGasGTA
    "MRU_khoums = 0.2 * MRU = khoums_MRU = MRU_minor",

    # --- MUR (Mauritian Rupee) ---
    # #Mauritius #IndianOcean #SugarRefining #Textiles #FreeportLogistics #FinancialHub
    "MUR_cent = 0.01 * MUR = cent_MUR = MUR_minor",

    # --- MVR (Maldivian Rufiyaa) ---
    # #Maldives #IndianOcean #SkipjackTuna #TourismEconomy #AgriFuelImports
    "MVR_laari = 0.01 * MVR = laari_MVR = MVR_minor",

    # --- MWK (Malawian Kwacha) ---
    # #Malawi #SouthernAfrica #BurleyTobacco #TeaExport #SugarCane #UraniumReserves
    "MWK_tambala = 0.01 * MWK = tambala_MWK = MWK_minor",

    # --- MXN (Mexican Peso) ---
    # #Mexico #NorthAmerica #SilverMiningGiant #MayaCrude #ManufacturingUSMCA #Automotive
    "MXN_centavo = 0.01 * MXN = centavo_MXN = MXN_minor",

    # --- MYR (Malaysian Ringgit) ---
    # #Malaysia #SoutheastAsia #PalmOilPowerhouse #PetronasLNG #Electronics #Tin
    "MYR_sen = 0.01 * MYR = sen_MYR = MYR_minor",

    # --- MZN (Mozambican Metical) ---
    # #Mozambique #SouthernAfrica #LNGProjectsCaboDelgado #CokingCoal #Aluminium #TitaniumSands
    "MZN_centavo = 0.01 * MZN = centavo_MZN = MZN_minor",

    # --- NAD (Namibian Dollar) ---
    # #Namibia #SouthernAfrica #UraniumProducer #Diamonds #Zinc #ZARPpegged
    "NAD_cent = 0.01 * NAD = cent_NAD = NAD_minor",

    # --- NGN (Nigerian Naira) ---
    # #Nigeria #WestAfrica #OPEC #BonnyLight #LNGExporter #DangoteRefinery #PetroEconomy
    "NGN_kobo = 0.01 * NGN = kobo_NGN = NGN_minor",

    # --- NIO (Nicaraguan Córdoba) ---
    # #Nicaragua #CentralAmerica #GoldMining #Coffee #BeefExports #SugarCane
    "NIO_centavo = 0.01 * NIO = centavo_NIO = NIO_minor",

    # --- NOK (Norwegian Krone) ---
    # #Norway #Nordic #PetroleumExporter #Equinor #NaturalGasToEurope #SalmonExport #HydroAluminium
    "NOK_ore = 0.01 * NOK = ore_NOK = NOK_minor",

    # --- NPR (Nepalese Rupee) ---
    # #Nepal #SouthAsia #HydropowerPotential #INRPegged #AgriImports #RemittanceEconomy
    "NPR_paisa = 0.01 * NPR = paisa_NPR = NPR_minor",

    # --- NZD (New Zealand Dollar) ---
    # #NewZealand #CommodityCurrency #DairyPowerhouse #Fonterra #MeatExport #Timber #G10
    "NZD_cent = 0.01 * NZD = cent_NZD = c_NZD = NZD_minor",

    # --- OMR (Omani Rial - 3 Decimals) ---
    # #Oman #MiddleEast #ThreeDecimal #OmanCrudeBenchmark #LNGExporter #SoharPort #USDPegged #Baisa
    "OMR_baisa = 0.001 * OMR = baisa_OMR = OMR_minor",

    # --- PAB (Panamanian Balboa) ---
    # #Panama #CentralAmerica #PanamaCanal #CopperMining #FreeTradeZone #USDPegged
    "PAB_centesimo = 0.01 * PAB = centesimo_PAB = PAB_minor",

    # --- PEN (Peruvian Sol) ---
    # #Peru #SouthAmerica #CopperTitan #SilverProducer #Zinc #Fishmeal #Gold
    "PEN_centimo = 0.01 * PEN = centimo_PEN = PEN_minor",

    # --- PGK (Papua New Guinean Kina) ---
    # #PapuaNewGuinea #Oceania #LNGExport #GoldMining #CopperPorgera #PalmOil #Coffee
    "PGK_toea = 0.01 * PGK = toea_PGK = PGK_minor",

    # --- PHP (Philippine Peso) ---
    # #Philippines #SoutheastAsia #NickelOreExporter #CopperGold #Semiconductors #CoconutOil
    "PHP_sentimo = 0.01 * PHP = sentimo_PHP = centavo_PHP = PHP_minor",

    # --- PKR (Pakistani Rupee) ---
    # #Pakistan #SouthAsia #TextileCotton #BasmatiRice #LNGBuyer #AgriEconomy
    "PKR_paisa = 0.01 * PKR = paisa_PKR = PKR_minor",

    # --- PLN (Polish Złoty) ---
    # #Poland #CentralEurope #KGHMCopper #SilverProducer #ThermalCoal #Manufacturing
    "PLN_grosz = 0.01 * PLN = grosz_PLN = PLN_minor",

    # --- PYG (Paraguayan Guaraní - Nominal/Historical) ---
    # #Paraguay #SouthAmerica #SoybeanExport #ItaipuHydropower #BeefTrade #ZeroDecimal
    "PYG_centimo = 0.01 * PYG = centimo_PYG = PYG_minor",

    # --- QAR (Qatari Riyal) ---
    # #Qatar #MiddleEast #LNGSuperpower #NorthField #HeliumProducer #USDPegged
    "QAR_dirham = 0.01 * QAR = dirham_QAR = QAR_minor",

    # --- RON (Romanian Leu) ---
    # #Romania #EasternEurope #BlackSeaGas #CrudeRefining #GrainExporter #Automotive
    "RON_ban = 0.01 * RON = ban_RON = RON_minor",

    # --- RSD (Serbian Dinar) ---
    # #Serbia #Balkans #BorCopperMine #LithiumJadar #AgriGrains #Manufacturing
    "RSD_para = 0.01 * RSD = para_RSD = RSD_minor",

    # --- RUB (Russian Ruble) ---
    # #Russia #Eurasia #UralsCrude #NaturalGas #WheatExporter #NickelNorilsk #Palladium #Fertilisers
    "RUB_kopeck = 0.01 * RUB = kopeck_RUB = kopeyka_RUB = RUB_minor",

    # --- RWF (Rwandan Franc - Nominal/Historical) ---
    # #Rwanda #EastAfrica #ColtanMining #TinTungsten #CoffeeTea #ZeroDecimal
    "RWF_centime = 0.01 * RWF = centime_RWF = RWF_minor",

    # --- SAR (Saudi Riyal) ---
    # #SaudiArabia #MiddleEast #OPECLeader #SaudiAramco #ArabLight #Petrochemicals #USDPegged
    "SAR_halala = 0.01 * SAR = halala_SAR = SAR_minor",

    # --- SBD (Solomon Islands Dollar) ---
    # #SolomonIslands #Oceania #TimberExport #RoughPalmOil #Fisheries #BauxiteNickel
    "SBD_cent = 0.01 * SBD = cent_SBD = SBD_minor",

    # --- SCR (Seychellois Rupee) ---
    # #Seychelles #IndianOcean #TunaCanningHub #MaritimeLogistics #TourismEconomy
    "SCR_cent = 0.01 * SCR = cent_SCR = SCR_minor",

    # --- SDG (Sudanese Pound) ---
    # #Sudan #NorthAfrica #NileBlendCrude #GoldMining #GumArabicMonopoly #AgriPotential
    "SDG_piastre = 0.01 * SDG = piastre_SDG = qirsh_SDG = SDG_minor",

    # --- SEK (Swedish Krona) ---
    # #Sweden #Nordic #IronOreLKAB #ForestryTimber #SteelManufacturing #HydroPower #G10
    "SEK_ore = 0.01 * SEK = ore_SEK = SEK_minor",

    # --- SGD (Singapore Dollar) ---
    # #Singapore #SoutheastAsia #OilRefiningHub #BunkeringCapital #CommodityTradingHub #GoldStorage
    "SGD_cent = 0.01 * SGD = cent_SGD = c_SGD = SGD_minor",

    # --- SHP (Saint Helena Pound) ---
    # #SaintHelena #SouthAtlantic #Fisheries #Coffee #GBPPegged #RemoteEconomy
    "SHP_penny = 0.01 * SHP = penny_SHP = pence_SHP = SHP_minor",

    # --- SLE (Sierra Leonean Leone) ---
    # #SierraLeone #WestAfrica #DiamondMining #RutileTitanium #Bauxite #IronOre
    "SLE_cent = 0.01 * SLE = cent_SLE = SLE_minor",

    # --- SOS (Somali Shilling) ---
    # #Somalia #HornOfAfrica #LivestockExport #Frankincense #FisheriesReserves #RedSeaCoast
    "SOS_senti = 0.01 * SOS = senti_SOS = cent_SOS = SOS_minor",

    # --- SRD (Surinamese Dollar) ---
    # #Suriname #SouthAmerica #OffshoreOilBasin #GoldMining #AluminaHistory #Forestry
    "SRD_cent = 0.01 * SRD = cent_SRD = SRD_minor",

    # --- SSP (South Sudanese Pound) ---
    # #SouthSudan #EastAfrica #DarBlendCrude #OilPipelinesToPortSudan #AgriPotential
    "SSP_piastre = 0.01 * SSP = piastre_SSP = SSP_minor",

    # --- STN (São Tomé and Príncipe Dobra) ---
    # #SaoTomePrincipe #CentralAfrica #CocoaExport #OffshoreOilExploration #JDZ
    "STN_centimo = 0.01 * STN = centimo_STN = STN_minor",

    # --- SVC (Salvadoran Colón - Nominal/Historical) ---
    # #ElSalvador #CentralAmerica #CoffeeExport #SugarCane #GeothermalPower #USDDollarised
    "SVC_centavo = 0.01 * SVC = centavo_SVC = SVC_minor",

    # --- SYP (Syrian Pound) ---
    # #Syria #MiddleEast #LightCrudeFields #PhosphateRock #AgriCotton #Wheat
    "SYP_piastre = 0.01 * SYP = piastre_SYP = qirsh_SYP = SYP_minor",

    # --- SZL (Eswatini Lilangeni) ---
    # #Eswatini #SouthernAfrica #SugarCaneGiant #WoodPulp #ZARPpegged #AgriTrade
    "SZL_cent = 0.01 * SZL = cent_SZL = SZL_minor",

    # --- THB (Thai Baht) ---
    # #Thailand #SoutheastAsia #NaturalRubberLeader #RiceExporter #SugarCane #Petrochemicals
    "THB_satang = 0.01 * THB = satang_THB = THB_minor",

    # --- TJS (Tajikistani Somoni) ---
    # #Tajikistan #CentralAsia #AluminiumTalco #HydropowerRogun #Cotton #GoldAntimony
    "TJS_diram = 0.01 * TJS = diram_TJS = TJS_minor",

    # --- TMT (Turkmenistan Manat) ---
    # #Turkmenistan #CentralAsia #GalkynyshGasField #PipelineToChina #CottonHarvest #Petrochemicals
    "TMT_tennesi = 0.01 * TMT = tennesi_TMT = tenge_TMT = TMT_minor",

    # --- TND (Tunisian Dinar - 3 Decimals) ---
    # #Tunisia #NorthAfrica #ThreeDecimal #PhosphateMining #OliveOilExporter #CrudeHydrocarbons #Millime
    "TND_millime = 0.001 * TND = millime_TND = TND_minor",

    # --- TOP (Tongan Paʻanga) ---
    # #Tonga #Pacific #VanillaExport #RootCrops #Fisheries #RemittanceEconomy
    "TOP_seniti = 0.01 * TOP = seniti_TOP = TOP_minor",

    # --- TRY (Turkish Lira) ---
    # #Turkey #Eurasia #EnergyTransitHub #BoronMineralsMonopoly #SteelExporter #Agriculture
    "TRY_kurus = 0.01 * TRY = kurus_TRY = TRY_minor",

    # --- TTD (Trinidad and Tobago Dollar) ---
    # #TrinidadTobago #Caribbean #LNGExporter #AmmoniaMethanolCapital #PetroChemicals
    "TTD_cent = 0.01 * TTD = cent_TTD = TTD_minor",

    # --- TWD (New Taiwan Dollar - Multi-Tier) ---
    # #Taiwan #EastAsia #SemiconductorFoundry #HighTechConsumer #EnergyImporter #LNGBuyer
    "TWD_jiao = 0.1 * TWD = jiao_TWD",
    "TWD_cent = 0.01 * TWD = cent_TWD = fen_TWD = TWD_minor",

    # --- TZS (Tanzanian Shilling) ---
    # #Tanzania #EastAfrica #GoldMining #GraphiteReserves #LNGProjectLindi #Cashews #Coffee
    "TZS_senti = 0.01 * TZS = senti_TZS = cent_TZS = TZS_minor",

    # --- UAH (Ukrainian Hryvnia) ---
    # #Ukraine #EasternEurope #BreadbasketOfEurope #WheatCorn #SunflowerOil #IronOrePellets
    "UAH_kopiyka = 0.01 * UAH = kopiyka_UAH = UAH_minor",

    # --- UGX (Ugandan Shilling - Nominal/Historical) ---
    # #Uganda #EastAfrica #RobustaCoffee #TilengaKingfisherOil #EACOP #ZeroDecimal
    "UGX_cent = 0.01 * UGX = cent_UGX = UGX_minor",

    # --- USD (United States Dollar - Multi-Tier) ---
    # #UnitedStates #GlobalReserveCurrency #WTIReference #HenryHub #CBOTGrains #Petrodollar
    "USD_cent = 0.01 * USD = US_cent = cent_USD = cents_USD = c_USD = cent = cents = USD_minor",
    "USD_mill = 0.001 * USD = mill_USD = mil_USD = mill = mil",

    # --- UYU (Uruguayan Peso) ---
    # #Uruguay #SouthAmerica #BeefExporter #CellulosePulp #Soybeans #WindHydropower
    "UYU_centesimo = 0.01 * UYU = centesimo_UYU = UYU_minor",

    # --- UZS (Uzbekistani Soʻm) ---
    # #Uzbekistan #CentralAsia #MuruntauGoldMine #UraniumProducer #CottonExport #NaturalGas
    "UZS_tiyin = 0.01 * UZS = tiyin_UZS = UZS_minor",

    # --- VES (Venezuelan Bolívar Soberano) ---
    # #Venezuela #SouthAmerica #OPECFounding #OrinocoHeavyCrude #LargestOilReserves #Refineries
    "VES_centimo = 0.01 * VES = centimo_VES = VES_minor",

    # --- VND (Vietnamese Đồng - Multi-Tier Nominal/Historical) ---
    # #Vietnam #SoutheastAsia #RobustaCoffeeLeader #RiceExporter #Electronics #ZeroDecimal
    "VND_hao = 0.1 * VND = hao_VND",
    "VND_xu = 0.01 * VND = xu_VND = VND_minor",

    # --- VUV (Vanuatu Vatu - Nominal/Historical) ---
    # #Vanuatu #Oceania #CopraExport #BeefProduction #MaritimeLogistics #ZeroDecimal
    "VUV_centime = 0.01 * VUV = centime_VUV = VUV_minor",

    # --- WST (Samoan Tālā) ---
    # #Samoa #Pacific #FishExports #TaroRoot #CoconutOil #TourismEconomy
    "WST_sene = 0.01 * WST = sene_WST = WST_minor",

    # --- XAF (Central African CFA Franc - Nominal/Historical) ---
    # #CEMAC #CentralAfrica #CrudeOilExporters #TimberForestry #EuroPegged #ZeroDecimal
    "XAF_centime = 0.01 * XAF = centime_XAF = XAF_minor",

    # --- XCD (East Caribbean Dollar) ---
    # #OECS #Caribbean #CurrencyUnion #BananasNutmeg #TourismLogistics #USDPegged
    "XCD_cent = 0.01 * XCD = cent_XCD = XCD_minor",

    # --- XCG (Caribbean Guilder) ---
    # #CuracaoSintMaarten #Caribbean #RefiningHistory #BunkeringLogistics #Tourism #USDPegged
    "XCG_cent = 0.01 * XCG = cent_XCG = XCG_minor",

    # --- XOF (West African CFA Franc - Nominal/Historical) ---
    # #WAEMU #WestAfrica #CocoaPowerhouse #GoldMining #BCEAO #EuroPegged #ZeroDecimal
    "XOF_centime = 0.01 * XOF = centime_XOF = XOF_minor",

    # --- XPF (CFP Franc - Nominal/Historical) ---
    # #FrenchPacific #NewCaledonia #NickelMiningPowerhouse #FrenchPolynesia #EuroPegged #ZeroDecimal
    "XPF_centime = 0.01 * XPF = centime_XPF = XPF_minor",

    # --- YER (Yemeni Rial) ---
    # #Yemen #MiddleEast #MaribLNG #RasIsaTerminal #BabElMandebTransit #CoffeeMocha
    "YER_fils = 0.01 * YER = fils_YER = YER_minor",

    # --- ZAR (South African Rand) ---
    # #SouthAfrica #SouthernAfrica #PGMSuperpower #PlatinumPalladium #GoldMining #CoalExporter
    "ZAR_cent = 0.01 * ZAR = cent_ZAR = ZAR_minor",

    # --- ZMW (Zambian Kwacha) ---
    # #Zambia #SouthernAfrica #CopperbeltPowerhouse #CobaltProducer #Emeralds #Hydropower
    "ZMW_ngwee = 0.01 * ZMW = ngwee_ZMW = ZMW_minor",

    # --- ZWG (Zimbabwe Gold) ---
    # #Zimbabwe #SouthernAfrica #GoldBacked #PGMReserves #LithiumMining #TobaccoExport
    "ZWG_cent = 0.01 * ZWG = cent_ZWG = ZWG_minor",
)


# ---------------------------------------------------------------------------
# Master Commodity Definitions
# ---------------------------------------------------------------------------

# Mutable hook for desk configurations or application settings prior to registry build
CUSTOM_COMMODITY_UNIT_DEFINITIONS: tuple[str, ...] = ()

COMMODITY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = (
    *CURRENCY_UNIT_DEFINITIONS,
    *CURRENCY_MINOR_UNIT_DEFINITIONS,

    # =========================================================================
    # Financial Margins, Quality Adjustments & Basis Points
    # =========================================================================
    "basis_point = 0.0001 = bp = bps = bips",  # 0.01 percent
    "percent = 0.01 = pct = percentage = %",

    # --- Time Constants & Operational Year Baselines ---
    # Standard 365-day operational calendar year (8,760 hours) universally adopted in
    # mining technical feasibility studies (JORC Code, NI 43-101), LNG Sale and Purchase
    # Agreements (GIIGNL SPA standards), and bulk terminal port allocation rules.
    "trading_year_365 = 365 * day = year_365 = annum_365 = op_year",                                   # Operational Year Basis (365 days / 8,760 hours nameplate capacity standard)


    # =========================================================================
    # Thermodynamic Energy & Gas Billing Units
    # =========================================================================

    # --- Energy Units: International Table Basis (Default Cash/Futures Anchor) ---
    # Defined under ISO 80000-5 and ISO 31-4 via the 1956 London International Steam
    # Table Conference (1 cal_IT = 4.1868 J exact), the International Table Btu
    # establishes an algebraically fixed Joule equivalent free from empirical
    # temperature variance, providing the reference energy basis for physical gas
    # delivery under CME Rulebook Ch. 220 (60°F, 14.73 psia) and ISO 13443 (15°C, 101.325 kPa).
    "btu_it = 1055.05585262 * joule = btu = Btu = BTU = btu_IT = Btu_IT = BTU_IT",  # British Thermal Unit (Exact IT definition, 1055.05585262 J)
    "kBtu_it = 1e3 * btu_it = kBtu = kbtu = KBTU = kBtu_it = KBtu_it",              # Kilo-British Thermal Unit (1,000 Btu)
    "MMBtu_it = 1e6 * btu_it = MMBtu = mmbtu = MMBTU = mmBtu = MMBtu_IT",            # Million British Thermal Units (NYMEX/ICE natural gas benchmark)
    "BBtu_it = 1e9 * btu_it = BBtu = bbtu = BBTU = BBtu_IT",                        # Billion British Thermal Units (1,000 MMBtu)
    "TBtu_it = 1e12 * btu_it = TBtu = tbtu = TBTU = TBtu_IT",                      # Trillion British Thermal Units (LNG cargo / regional balance scale)
    "quad_it = 1e15 * btu_it = quad = Quad = QUAD = quadBtu = QuadBtu = QUADBTU",    # Quadrillion British Thermal Units (National/global energy balances)

    # US Utility & Retail Gas Contracts (100,000 Btu)
    "therm_it = 1e5 * btu_it = therm = th = Therm = THERM = th_it",                 # Therm (100,000 Btu, utility retail gas billing standard)
    "dth_it = 10 * therm_it = dth = Dth = DTH = dekatherm = Dekatherm = DEKATHERM",  # Dekatherm (10 Therms or 1 MMBtu, US pipeline capacity tariff standard)

    # --- Energy Units: US Regulatory & Custody Transfer Basis (59°F / 15°C) ---
    # Defined under NIST SP 811, ASTM E380, and AGA standards as the heat required
    # to raise 1 lb of air-free water from 58.5°F to 59.5°F at 1 atm (14.696 psia /
    # 101.325 kPa), the 59°F Btu (1054.804 J, corresponding to the 15°C calorie)
    # serves as the legal calorimetric baseline for legacy North American utility
    # rate tariffs, statutory gas billing, and state pipeline transport contracts.
    "btu_us59 = 1054.804 * joule = btu_59F = Btu_59F = BTU_59F = btu_us = Btu_us",  # British Thermal Unit (59°F / 15°C water-heat basis, 1054.804 J)
    "kBtu_us59 = 1e3 * btu_us59 = kBtu_59F = kbtu_59F = KBTU_59F = kBtu_us",        # Kilo-British Thermal Unit (59°F basis)
    "MMBtu_us59 = 1e6 * btu_us59 = MMBtu_59F = mmbtu_59F = MMBTU_59F = MMBtu_us",    # Million British Thermal Units (59°F basis)
    "BBtu_us59 = 1e9 * btu_us59 = BBtu_59F = bbtu_59F = BBTU_59F = BBtu_us",        # Billion British Thermal Units (59°F basis)
    "TBtu_us59 = 1e12 * btu_us59 = TBtu_59F = tbtu_59F = TBTU_59F = TBtu_us",      # Trillion British Thermal Units (59°F basis)
    "quad_us59 = 1e15 * btu_us59 = quad_59F = Quad_59F = QUAD_59F = quad_us",      # Quadrillion British Thermal Units (59°F basis)

    # US 59°F Utility & Pipeline Billing
    "therm_us59 = 1e5 * btu_us59 = therm_59F = th_59F = therm_us = th_us",          # Therm (59°F basis, 100,000 Btu_59F)
    "dth_us59 = 10 * therm_us59 = dth_59F = Dth_59F = dekatherm_59F = dth_us",      # Dekatherm (59°F basis, 10 Therms_59F)

    # =========================================================================
    # Macro Energy, Upstream Reserves & Coal Thermodynamics
    # =========================================================================

    # --- Upstream Hydrocarbon Reserves & Production (BOE Framework) ---
    # Defined under US IRS § 45K / § 29 thermal equivalence (5.80 MMBtu = 5,800 scf @ 1,000 Btu/scf)
    # as the baseline physical commodity conversion, alongside the statutory US SEC Rule 4-10
    # (Regulation S-X) reporting standard (6.00 MMBtu = 6,000 scf @ 1,000 Btu/scf) used in public E&P
    # 10-K reserves disclosures and Society of Petroleum Engineers (SPE-PRMS) guidelines.
    "boe_58 = 5.8 * MMBtu_it = boe = BOE = barrel_of_oil_equivalent",                  # Barrel of Oil Equivalent (IRS standard thermal basis, 5.8 MMBtu / ~5,800 scf)
    "boe_sec = 6.0 * MMBtu_it = boe_6 = BOE_SEC = boe_S_X",                            # Barrel of Oil Equivalent (US SEC Reg S-X Rule 4-10 statutory basis, 6.0 MMBtu / 6,000 scf)
    "kboe = 1e3 * boe = KBOE = kBOE",                                                  # Thousand Barrels of Oil Equivalent (kBOE, 1,000 boe)
    "MMboe = 1e6 * boe = MMBOE = mmboe",                                               # Million Barrels of Oil Equivalent (MMBOE, field/basin reserve scale)
    "Bboe = 1e9 * boe = BBOE = bboe",                                                  # Billion Barrels of Oil Equivalent (BBOE, national/global reserve benchmark)

    # Upstream Production Flow Rates
    "boed = boe / day = BOED = boe_per_day = boe_d",                                   # Barrels of Oil Equivalent per Day (BOED, field-level production rate)
    "kboed = 1e3 * boed = KBOED = kboe_d = Mboed",                                     # Thousand Barrels of Oil Equivalent per Day (kBOED / Mboed, basin production rate)
    "MMboed = 1e6 * boed = MMBOED = mmboed = MMboe_d",                                 # Million Barrels of Oil Equivalent per Day (MMBOED, IEA/OPEC global liquids aggregate)

    # --- Macro Thermodynamic Balances (IEA / OECD / UN Standards) ---
    # Defined by the International Energy Agency (IEA), OECD, and Eurostat as exactly 10 Gcal_IT
    # (41.868 GJ net calorific value), the Tonne of Oil Equivalent (toe) normalises cross-commodity
    # macro energy statistics. The Tonne of Coal Equivalent (tce), codified under UN Statistics
    # Division standards and DIN 54115, represents exactly 7 Gcal_IT (29.3076 GJ net calorific value),
    # corresponding to the energy yield of one metric tonne of standard hard coal (7,000 kcal/kg).
    "toe = 41.868 * gigajoule = TOE = tonne_of_oil_equivalent",                        # Tonne of Oil Equivalent (TOE, IEA/OECD standard net calorific basis, 41.868 GJ / 10 Gcal)
    "ktoe = 1e3 * toe = KTOE = ktoe_it",                                               # Kilotonne of Oil Equivalent (kTOE, 1,000 toe, national energy balance accounting)
    "Mtoe = 1e6 * toe = MTOE = mtoe = Mtoe_it",                                        # Megatonne of Oil Equivalent (MTOE, 1,000,000 toe, global primary energy demand)

    "tce = 29.3076 * gigajoule = TCE = tonne_of_coal_equivalent",                      # Tonne of Coal Equivalent (TCE, UN/DIN 54115 standard, 29.3076 GJ / 7 Gcal)
    "ktce = 1e3 * tce = KTCE = ktce_it",                                               # Kilotonne of Coal Equivalent (kTCE, 1,000 tce, regional coal balance accounting)
    "Mtce = 1e6 * tce = MTCE = mtce = Mtce_it",                                        # Megatonne of Coal Equivalent (MTCE, 1,000,000 tce, global thermal coal demand)

    # --- International Table Caloric Hierarchy (Coal & Power Generation) ---
    # Defined via the 1956 International Steam Table Conference standard (1 cal_IT = 4.1868 J exact)
    # used in international seaborne coal swap settlement (API 2, API 4, Newcastle FOB) and
    # Japanese / Korean bilateral LNG import pricing (Customs Cleared LNG / JLC calorific tests).
    "cal_it = 4.1868 * joule = cal = Cal = cal_IT",                                    # Small Calorie (Exact 1956 International Steam Table definition, 4.1868 J)
    "kcal_it = 1e3 * cal_it = kcal = Kcal = KCAL = kcal_IT",                           # Kilocalorie (kcal, seaborne coal pricing benchmark, e.g. 6,000 kcal/kg NAR)
    "mcal_it = 1e6 * cal_it = mcal = Mcal = MCAL = megacalorie",                       # Megacalorie (Mcal, 1,000 kcal / 1e6 cal, district heating & industrial heat)
    "gcal_it = 1e9 * cal_it = gcal = Gcal = GCAL = gigacalorie",                       # Gigacalorie (Gcal, 1e9 cal / 4.1868 GJ, Eastern European gas & power tariffs)

    # =========================================================================
    # Gas Reference-Volume Units (Thermodynamic Standard Conditions)
    # =========================================================================

    # NGL gas-equivalent volume at 60 degrees F and one standard atmosphere.
    # Kept dimensionally distinct from FERC's 14.73-psia custody-transfer basis.
    "scf_60F_14696psi_dry = [standard_gas_volume_60F_1atm_dry] = scf_60F_1atm_dry",

    # GOST reference volume: dry gas at 20 degrees C and 101.325 kPa.
    "m3_20C_101325kPa_dry = [standard_gas_volume_20C_101325kPa_dry] = Sm3_20C",

    # --- US Standard Gas Volumes: FERC / Interstate / NYMEX (60°F, 14.73 psia) ---
    # Defined under AGA Report No. 3 (API MPMS Ch. 14.3), ASTM D3588, and GPA 2172 at 60°F (288.71 K)
    # and 14.73 psia (30.00 inHg / 101.5598 kPa absolute, dry), representing the legal custody-transfer
    # benchmark for US interstate pipeline tariffs (FERC) and CME/NYMEX Henry Hub physical delivery.
    "scf_60F_14_73psia = [standard_gas_volume_60F_14_73psia] = scf = standard_cubic_foot = scf_14_73", # Standard Cubic Foot (60°F, 14.73 psia, US pipeline & NYMEX delivery basis)
    "ccf_60F_14_73psia = 100 * scf_60F_14_73psia = ccf = Ccf = CCF",                                   # Hundred Standard Cubic Feet (100 scf, US retail utility volumetric billing)
    "Mscf_60F_14_73psia = 1e3 * scf_60F_14_73psia = Mscf = mscf = Mcf = MCF",                          # Thousand Standard Cubic Feet (1,000 scf, standard US upstream well volume)
    "MMscf_60F_14_73psia = 1e6 * scf_60F_14_73psia = MMscf = mmscf = MMcf = MMCF",                     # Million Standard Cubic Feet (1,000 Mcf, US pipeline dispatch & processing scale)
    "Bscf_60F_14_73psia = 1e9 * scf_60F_14_73psia = Bscf = bscf = Bcf = BCF",                          # Billion Standard Cubic Feet (1,000 MMcf, US storage capacity & regional flow)
    "Tscf_60F_14_73psia = 1e12 * scf_60F_14_73psia = Tscf = tscf = Tcf = TCF",                         # Trillion Standard Cubic Feet (1,000 Bcf, national basin & EIA reserve accounting)

    # --- US Statutory Wellhead Gas Volumes: Texas / Midcontinent (60°F, 14.65 psia) ---
    # Statutory base pressure defined by the Texas Railroad Commission (Rule 34), Oklahoma, Louisiana,
    # and Kansas state conservation statutes (60°F / 288.71 K and 14.65 psia / 101.008 kPa absolute, dry).
    # Accounting for this 0.543% expansion difference is essential when reconciling gathering meters to interstate lines.
    "scf_60F_14_65psia = [standard_gas_volume_60F_14_65psia] = scf_texas = scf_tx = scf_14_65",       # Standard Cubic Foot (Texas RRC / Gulf Coast state statutory wellhead basis, 14.65 psia)
    "Mscf_60F_14_65psia = 1e3 * scf_60F_14_65psia = Mscf_tx = Mcf_tx = MCF_TX",                       # Thousand Standard Cubic Feet (14.65 psia statutory wellhead basis)
    "MMscf_60F_14_65psia = 1e6 * scf_60F_14_65psia = MMscf_tx = MMcf_tx = MMCF_TX",                   # Million Standard Cubic Feet (14.65 psia statutory wellhead basis)
    "Bscf_60F_14_65psia = 1e9 * scf_60F_14_65psia = Bscf_tx = Bcf_tx = BCF_TX",                       # Billion Standard Cubic Feet (14.65 psia statutory wellhead basis)

    # --- Normal Gas Volumes: Continental European Industrial Basis (0°C, 101.325 kPa) ---
    # Defined under DIN 1343, ISO 2533, and IUPAC reference conditions at 0°C (273.15 K) and 101.325 kPa
    # (1 atm absolute, dry), universally applied across Continental European industrial gas transport,
    # German pipeline network balancing, Central European underground gas storage (UGS), and refinery gas processing.
    "Nm3_0C_101325kPa_dry = [normal_gas_volume_0C_101325kPa_dry] = Nm3 = normal_m3 = nm3",            # Normal Cubic Metre (DIN 1343 basis, 0°C, 101.325 kPa, dry)
    "kNm3_0C_101325kPa_dry = 1e3 * Nm3_0C_101325kPa_dry = kNm3 = knm3 = kNm3_dry",                    # Thousand Normal Cubic Metres (1,000 Nm3)
    "MNm3_0C_101325kPa_dry = 1e6 * Nm3_0C_101325kPa_dry = MNm3 = mnm3 = MMnm3 = MNm3_dry",            # Million Normal Cubic Metres (1,000 kNm3, Continental European industrial scale)
    "BNm3_0C_101325kPa_dry = 1e9 * Nm3_0C_101325kPa_dry = BNm3 = bnm3 = Bnm3 = Bcm_normal",           # Billion Normal Cubic Metres (1,000 MNm3, Continental European macro balance)

    # --- Standard Metric Gas Volumes: ISO / UK / Global LNG Basis (15°C, 101.325 kPa) ---
    # Defined under ISO 13443 ("Natural gas - Standard reference conditions") at 15°C (288.15 K) and
    # 101.325 kPa (1 atm absolute, dry). The global physical benchmark adopted by ENTSOG, the UK National
    # Grid (NTS), Norwegian offshore pipelines (Gassco), Australian LNG plants, and seaborne GIIGNL trade.
    "Sm3_15C_101325kPa_dry = [standard_gas_volume_15C_101325kPa_dry] = Sm3 = standard_m3 = sm3",      # Standard Cubic Metre (ISO 13443 basis, 15°C, 101.325 kPa, dry)
    "kSm3_15C_101325kPa_dry = 1e3 * Sm3_15C_101325kPa_dry = kSm3 = ksm3 = kSm3_dry",                  # Thousand Standard Cubic Metres (1,000 Sm3, Norwegian continental shelf flow metric)
    "MSm3_15C_101325kPa_dry = 1e6 * Sm3_15C_101325kPa_dry = MSm3 = msm3 = MMsm3 = MSm3_dry",          # Million Standard Cubic Metres (1,000 kSm3, European border interconnect daily rate)
    "BSm3_15C_101325kPa_dry = 1e9 * Sm3_15C_101325kPa_dry = BSm3 = bsm3 = Bcm = BCM = bcm",           # Billion Standard Cubic Metres (1,000 MSm3, global LNG & European gas benchmark)

    # --- Groningen Caloric Standard Gas Volumes: Dutch GTS / TTF Standard ---
    # Standardised low-calorific gas volume (G-gas) defined by Gasunie Transport Services (GTS) with a
    # standardised gross calorific value (GCV) of exactly 35.17 MJ/m3 (9.7694 kWh/m3) at 0°C / 101.325 kPa.
    # Serves as the volumetric reference for physical TTF title transfers and domestic Dutch/German L-gas billing.
    "m3_Groningen = [groningen_gas_volume] = m3_G = Nm3_G = G_m3",                                     # Groningen Standard Cubic Metre (GTS low-cal G-gas benchmark, 35.17 MJ/Nm3)
    "thousand_m3_Groningen = 1e3 * m3_Groningen = km3_G = kNm3_G",                                     # Thousand Groningen Cubic Metres (1,000 m3_G)
    "million_m3_Groningen = 1e6 * m3_Groningen = MMm3_G = MNm3_G",                                     # Million Groningen Cubic Metres (1,000 km3_G, GTS daily balancing standard)
    "billion_m3_Groningen = 1e9 * m3_Groningen = Bcm_G = BCM_G = BNm3_G",                              # Billion Groningen Cubic Metres (1,000 MMm3_G, annual Dutch capacity & storage)

    # =========================================================================
    # Gas Flow Rates, Pipeline Nominations & Capacity Bookings
    # =========================================================================

    # --- US Volumetric Flow Rates: Interstate / FERC Basis (60°F, 14.73 psia) ---
    # Standard volumetric flow rates governing US interstate pipeline electronic bulletin boards (EBB),
    # NAESB nomination cycles, gathering system dispatch, and EIA production reporting.
    "scfd_60F = scf_60F_14_73psia / day = scfd = SCFD = scf_d",                                       # Standard Cubic Feet per Day (scf / day)
    "scfh_60F = scf_60F_14_73psia / hour = scfh = SCFH = scf_h",                                     # Standard Cubic Feet per Hour (scf / hour)
    "Mscfd_60F = Mscf_60F_14_73psia / day = Mscfd = mscfd = Mcfd = mcfd = MCFD = Mcf_d",              # Thousand Standard Cubic Feet per Day (1,000 scfd)
    "Mscfh_60F = Mscf_60F_14_73psia / hour = Mscfh = mscfh = Mcfh = mcfh = MCFH = Mcf_h",            # Thousand Standard Cubic Feet per Hour (1,000 scfh)
    "MMscfd_60F = MMscf_60F_14_73psia / day = MMscfd = mmscfd = MMcfd = mmcfd = MMCFD = MMcf_d",     # Million Standard Cubic Feet per Day (1,000 Mcfd, standard pipeline nomination unit)
    "MMscfh_60F = MMscf_60F_14_73psia / hour = MMscfh = mmscfh = MMcfh = mmcfh = MMCFH = MMcf_h",   # Million Standard Cubic Feet per Hour (1,000 Mcfh, CCGT hourly swing burn)
    "Bcfd_60F = Bscf_60F_14_73psia / day = Bcfd = bcfd = BCFD = Bscfd = Bcf_d",                      # Billion Standard Cubic Feet per Day (1,000 MMcfd, major basin/national takeaway)

    # --- US Volumetric Flow Rates: Texas Statutory Wellhead Basis (60°F, 14.65 psia) ---
    # Statutory wellhead production and gathering flow rates compliant with Texas RRC Rule 34 and Midcontinent laws.
    "scfd_60F_14_65psia = scf_60F_14_65psia / day = scfd_tx = scf_d_tx",                              # Standard Cubic Feet per Day (Texas 14.65 psia basis)
    "Mscfd_60F_14_65psia = Mscf_60F_14_65psia / day = Mscfd_tx = Mcfd_tx = Mcf_d_tx",                # Thousand Standard Cubic Feet per Day (Texas 14.65 psia basis)
    "MMscfd_60F_14_65psia = MMscf_60F_14_65psia / day = MMscfd_tx = MMcfd_tx = MMcf_d_tx",          # Million Standard Cubic Feet per Day (Texas 14.65 psia basis)
    "Bcfd_60F_14_65psia = Bscf_60F_14_65psia / day = Bcfd_tx = Bcf_d_tx",                            # Billion Standard Cubic Feet per Day (Texas 14.65 psia basis)

    # --- Continental European Normal Volumetric Rates: DIN 1343 Basis (0°C, 101.325 kPa) ---
    # Used across Central European TSO dispatch, underground storage (UGS) injection/withdrawal, and industrial metering.
    "Nm3_per_hour_0C = Nm3_0C_101325kPa_dry / hour = Nm3_h = nm3_h = Nm3h",                          # Normal Cubic Metres per Hour (Nm3 / hour, industrial metering standard)
    "Nm3_per_day_0C = Nm3_0C_101325kPa_dry / day = Nm3_d = nm3_d = Nm3d",                            # Normal Cubic Metres per Day (Nm3 / day)
    "kNm3_per_hour_0C = kNm3_0C_101325kPa_dry / hour = kNm3_h = knm3_h",                             # Thousand Normal Cubic Metres per Hour (1,000 Nm3/h)
    "kNm3_per_day_0C = kNm3_0C_101325kPa_dry / day = kNm3_d = knm3_d",                               # Thousand Normal Cubic Metres per Day (1,000 Nm3/d)
    "MNm3_per_day_0C = MNm3_0C_101325kPa_dry / day = MNm3_d = mnm3_d = MMnm3_d",                     # Million Normal Cubic Metres per Day (1,000 kNm3/d, TSO grid dispatch)

    # --- Standard Metric Volumetric Rates: ISO 13443 Basis (15°C, 101.325 kPa) ---
    # Operating standard for Norwegian offshore pipe dispatch (Gassco), UK NTS terminal deliveries, and seaborne LNG.
    "Sm3_per_hour_15C = Sm3_15C_101325kPa_dry / hour = Sm3_h = sm3_h = Sm3h",                        # Standard Cubic Metres per Hour (Sm3 / hour)
    "Sm3_per_day_15C = Sm3_15C_101325kPa_dry / day = Sm3_d = sm3_d = Sm3d",                          # Standard Cubic Metres per Day (Sm3 / day)
    "kSm3_per_day_15C = kSm3_15C_101325kPa_dry / day = kSm3_d = ksm3_d",                             # Thousand Standard Cubic Metres per Day (1,000 Sm3/d)
    "MSm3_per_day_15C = MSm3_15C_101325kPa_dry / day = MSm3_d = msm3_d = MMsm3_d = mcm_d",           # Million Standard Cubic Metres per Day (1,000 kSm3/d, European cross-border flows)
    "Bcm_per_day_15C = BSm3_15C_101325kPa_dry / day = Bcm_d = bcm_d = BCMD",                         # Billion Standard Cubic Metres per Day (1,000 MSm3/d)

    # Annualised International Pipeline & LNG Terminal Capacity
    "bcma = BSm3_15C_101325kPa_dry / trading_year_365 = BCMA = bcm_a = Bcm_year",                     # Billion Cubic Metres per Annum (Long-term cross-border transit treaty standard)

    # --- Dutch Groningen Volumetric Rates: GTS Standard (Low-Cal G-Gas Basis) ---
    # Governs Dutch national transport, Gasunie Transport Services (GTS) entry/exit bookings, and physical TTF delivery.
    "m3_Groningen_per_hour = m3_Groningen / hour = m3_G_h = Nm3_G_h",                                  # Groningen Cubic Metres per Hour (m3_G / hour)
    "million_m3_Groningen_per_day = million_m3_Groningen / day = MMm3_G_d = MNm3_G_d",                 # Million Groningen Cubic Metres per Day (GTS system daily balancing capacity)

    # --- US Energy Nomination Rates: Thermal Basis (International Table) ---
    # Standard contractual nomination units for US interstate pipeline firm transportation (FT) and interruptible (IT) contracts.
    "MMBtu_per_day = MMBtu_it / day = mmbtu_d = MMBtud = MMBTU_d = mmBtu_d",                          # Million British Thermal Units per Day (US pipeline primary nomination unit)
    "MMBtu_per_hour = MMBtu_it / hour = mmbtu_h = MMBtuh = MMBTU_h = mmBtu_h",                        # Million British Thermal Units per Hour (Power plant hourly burn capability)
    "dth_per_day = dth_it / day = dth_d = Dth_d = DTHD = Dthd",                                       # Dekatherms per Day (10 Therms/day = 1 MMBtu/day, legal tariff standard)
    "dth_per_hour = dth_it / hour = dth_h = Dth_h = DTHH = Dthh",                                     # Dekatherms per Hour (Hourly meter nomination capacity)
    "therms_per_day = therm_it / day = th_d = therm_d = Therms_d",                                    # Therms per Day (US commercial & municipal utility nominations)

    # --- European & UK Energy Nomination Rates: CAM NC / PRISMA Standard ---
    # Defined under EU Regulation 2017/459 (CAM NC) and UK Uniform Network Code (UNC).
    # All European TSO auctions (PRISMA, RBP, GSA) and UK NTS daily entry/exit capacity operate strictly in electrical energy rates.
    "kWh_per_hour = kilowatt_hour / hour = kWh_h = kWh_per_h",                                        # Kilowatt Hours per Hour (PRISMA standard bundled capacity auction unit)
    "kWh_per_day = kilowatt_hour / day = kWh_d = kWh_per_d",                                          # Kilowatt Hours per Day (UK NTS primary daily entry/exit capacity unit)
    "MWh_per_hour = megawatt_hour / hour = MWh_h = MWh_per_h",                                        # Megawatt Hours per Hour (European interconnector hourly capacity nomination)
    "MWh_per_day = megawatt_hour / day = MWh_d = MWh_per_d",                                          # Megawatt Hours per Day (European border point aggregate daily flow)
    "GWh_per_day = gigawatt_hour / day = GWh_d = GWh_per_d",                                          # Gigawatt Hours per Day (UK/Continental high-pressure transmission line capacity)

    # =========================================================================
    # Petroleum Volumes, Refinery Throughput & Tanker Logistics
    # =========================================================================

    # --- Petroleum Volumetric Hierarchy: API MPMS Ch. 11 / ASTM D1250 Basis ---
    # Standard 42-US-gallon blue barrel defined at 60°F (15.56°C) and equilibrium vapor pressure
    # under API Manual of Petroleum Measurement Standards (MPMS) Chapter 11.1, ASTM D1250, and IP 200.
    # Enforces oilfield Roman numeral conventions (Mbbl = 1,000 bbl; MMbbl = 1,000,000 bbl).
    "bbl = oil_barrel = BBL = barrel_oil = bbls",                                                      # Standard Petroleum Barrel (42 US gallons / 158.9873 litres @ 60°F)
    "kbbl = 1e3 * bbl = Mbbl = MBBL = kbbls",                                                          # Thousand Barrels of Oil (1,000 bbl, standard oilfield M-notation)
    "MMbbl = 1e6 * bbl = mmbbl = MMBBL = MMbbls",                                                      # Million Barrels of Oil (1,000 Mbbl, commercial field and commercial tank farm scale)
    "Gbbl = 1e9 * bbl = gbbl = GBBL = Bbbl = BBBL",                                                    # Billion Barrels of Oil (1,000 MMbbl, national proved crude reserves benchmark)

    # --- Wholesale Refined Products & Rack Distribution ---
    # Standard volumetric measures for finished products (Gasoline, Diesel, Jet A-1, Gasoil).
    # NYMEX RBOB and Heating Oil contracts clear physically on 42,000-gallon lots (1,000 bbl).
    "gallon_petroleum = gallon = gal = Gal = us_gal",                                                  # US Liquid Gallon (231 cu in / 3.78541 litres, NYMEX refined products basis)
    "kgal = 1e3 * gallon_petroleum = Kgal = Mgal_us",                                                  # Thousand US Gallons (Wholesale truck rack and terminal rack dispatch)

    # --- Petroleum Production & Refinery Flow Rates ---
    # Daily production rates and refinery processing runs, distinguishing between calendar-day
    # averages (bpcd) and unconstrained nameplate stream-day operating capacities (bpsd).
    "bpd = bbl / day = bbl_d = BPD = bbl_per_day = barrels_per_day",                                   # Barrels per Day (Standard calendar-day production / throughput, bpcd)
    "bpcd = bpd = BPCD = barrels_per_calendar_day",                                                   # Barrels per Calendar Day (Throughput accounting for downtime and maintenance)
    "bpsd = bbl / day = BPSD = barrels_per_stream_day",                                                # Barrels per Stream Day (Maximum sustainable 24-hour CDU capacity)
    "kbpd = 1e3 * bpd = Mbpd = mbpd = KBPD = kbbl_d = Mbbl_d = kbd = Mbd",                            # Thousand Barrels per Day (1,000 bpd, regional pipeline flows and refinery throughput)
    "MMbpd = 1e6 * bpd = mmbpd = MMBPD = MMbbl_d = mmbbl_d = MMbd = mmbd",                            # Million Barrels per Day (1,000 kbpd, OPEC/IEA global liquids supply and demand)

    # --- Tanker Loading, Pipeline Pumping & Marine Laytime Logistics ---
    # Hourly transfer rates governing terminal nominations, NOR laytime, and demurrage calculations
    # under standard tanker charterparties (ASBATANKVOY, BPVOY5, ShellVOY6).
    "bph = bbl / hour = BPH = bbl_h = bbl_per_hour = barrels_per_hour",                                # Barrels per Hour (Marine loading-arm and pipeline pumping rate)
    "kbph = 1e3 * bph = Kbph = Mbph",                                                                  # Thousand Barrels per Hour (1,000 bph, VLCC / Suezmax loading terminal manifold capacity)

    # --- Metric Petroleum Volumes & Marine Bunkering Flows ---
    # Metric cubic metres and hourly discharge rates universally adopted across European, Asian,
    # and Middle Eastern tank storage terminals (ARA, Jurong, Fujairah) and bunker deliveries.
    "cbm = meter ** 3 = cbm_oil = CBM = m3_oil",                                                       # Cubic Metre (1,000 litres / ~6.2898 bbl, global metric liquid storage standard)
    "kcbm = 1e3 * cbm = kCBM = km3_oil",                                                               # Thousand Cubic Metres (1,000 cbm, commercial shore-tank terminal capacity)
    "Mcbm = 1e6 * cbm = MCBM = MMm3_oil",                                                              # Million Cubic Metres (1,000 kcbm, macro strategic petroleum reserve capacity)
    "cbm_h = cbm / hour = m3_per_hour = cbm_per_hour = m3_h = CBM_h",                                  # Cubic Metres per Hour (Terminal pumping rate and marine bunker manifold delivery)
    "cbm_d = cbm / day = m3_per_day = cbm_per_day = m3_d = CBM_d",                                     # Cubic Metres per Day (Refinery daily intake and cross-border liquid pipeline flows)

    # =========================================================================
    # Power Baseload Blocks, PPA Delivery & Environmental Attributes
    # =========================================================================

    # --- PPA Baseload Delivery Blocks & Contractual Strips ---
    # Standardised delivery blocks used in long-term Power Purchase Agreements (PPAs), tolling agreements,
    # and exchange-cleared power futures (EEX, Nord Pool, ICE Endex, CME PJM/ERCOT). Represents continuous
    # block generation over standard calendar intervals (1 MW baseload strip).
    "MWa = 8760 * megawatt_hour = MWA = megawatt_average = MW_year = MWyr",                            # Megawatt-Average (Standard 8,760 hours non-leap year baseload delivery, 8.76 GWh)
    "MWa_leap = 8784 * megawatt_hour = MWA_leap = MW_year_leap = MWyr_leap",                           # Megawatt-Average Leap Year (366 days / 8,784 hours baseload delivery, 8.784 GWh)
    "MW_day = 24 * megawatt_hour = MWd = MWD = megawatt_day",                                          # Megawatt Day (24-hour continuous baseload contract block, 24 MWh)
    "MW_peak_16 = 16 * megawatt_hour = MW_pk16 = peak_block_16h",                                      # Megawatt Peak Block (Standard 16-hour diurnal peak delivery block, 16 MWh)

    # --- North American Renewable Energy Certificates (EAC Framework) ---
    # Issued under US state Renewable Portfolio Standards (RPS) and tracked via regional registries
    # (PJM-GATS, NEPOOL-GIS, WREGIS, ERCOT, M-RETS). Exactly 1 MWh of certified renewable generation.
    "REC = [instrument_REC] = rec = renewable_energy_certificate = recs = RECs",                      # Renewable Energy Certificate (Standard US/Canada compliance & voluntary certificate, 1 MWh)
    "SREC = [instrument_SREC] = srec = solar_renewable_energy_certificate = srecs",                    # Solar Renewable Energy Certificate (US state RPS solar carve-out compliance unit, 1 MWh)

    # --- European & UK Guarantees of Origin (AIB / Ofgem) ---
    # Distinct legal attribute instruments post-Brexit. European GOs operate under Directive (EU) 2018/2001
    # (RED II) and the AIB EECS standard. UK REGOs are administered separately by Ofgem under statutory GB rules.
    "GO = [instrument_GO] = Guarantee_of_Origin = guarantees_of_origin = go = eu_go = EU_GO",        # Guarantee of Origin (EU RED II / AIB EECS standard tracking instrument, 1 MWh)
    "REGO = [instrument_REGO] = rego = renewable_energy_guarantee_of_origin = regos = REGOs",          # Renewable Energy Guarantee of Origin (UK Ofgem statutory compliance certificate, 1 MWh)

    # --- International & Multi-Regional Attribute Standards ---
    # Global tracking systems for corporate Scope 2 greenhouse gas reporting (GHG Protocol Scope 2 Guidance)
    # across non-AIB jurisdictions (Latin America, Middle East, Africa, Asia-Pacific).
    "IREC = [instrument_IREC] = irec = I_REC = international_rec = irecs",                             # International REC (I-REC Standard Foundation global tracking instrument, 1 MWh)
    "TIGR = 1 * megawatt_hour = tigr = tigrs = tradable_instrument_global_renewables",                 # Tradable Instrument for Global Renewables (APX international tracking standard, 1 MWh)
    "GEC = [instrument_GEC] = gec = green_electricity_certificate = lvezheng = GECs",                 # Green Electricity Certificate (China NDRC / NEA national renewable attribute, 1 MWh)

    # --- Australian Renewable Energy Target (Clean Energy Regulator) ---
    # Statutory instruments created under the Australian Renewable Energy (Electricity) Act 2000.
    # Encompasses both wholesale grid generation (LRET) and deemed small-scale installations (SRES).
    "LGC = [instrument_LGC] = lgc = large_scale_generation_certificate = lgcs = LGCs",                # Large-scale Generation Certificate (Australian LRET wholesale compliance unit, 1 MWh)
    "STC = [instrument_STC] = stc = small_scale_technology_certificate = stcs = STCs",                # Small-scale Technology Certificate (Australian SRES solar/deemed generation unit, 1 MWh)

    # --- Capacity Markets & Diurnal Block Shaping ---
    # Used in PJM Reliability Pricing Model (RPM), ISO-NE Forward Capacity Market, and grid balancing.
    "MW_month = 730 * megawatt_hour = MW_mo = MWm = megawatt_month",                                   # Megawatt Month (Average 730 hours monthly power capacity reservation block, 730 MWh)
    "kW_month = 0.001 * MW_month = kW_mo = kWm = kilowatt_month",                                      # Kilowatt Month (Standard US utility capacity tariff reservation unit)
    "MW_offpeak_8 = 8 * megawatt_hour = MW_op8 = offpeak_block_8h",                                    # Megawatt Off-Peak Block (Standard 8-hour nocturnal off-peak delivery block, 8 MWh)

    # =========================================================================
    # Mass, Bulk Metals, Precious Bullion & Mining Concentrates
    # =========================================================================

    # --- Imperial & Metric Mass Benchmarks (Exchange Clearance & Physical Trade) ---
    # Governed by the 1959 International Yard and Pound Agreement (1 lb = 0.45359237 kg exact)
    # and ISO 80000-4. Forms the primary settlement delivery units across the London Metal Exchange
    # (LME - metric tonnes), CME Group / COMEX (short tons, cwt, pounds), and SHFE (metric tonnes).
    "mt = 1e3 * kilogram = metric_tonne = tonne = Tonne = TONNE = MT",                                 # Metric Tonne (1,000 kg, global LME/SHFE base metals & bulk shipping standard)
    "kmt = 1e3 * mt = kt = KT = kMT = thousand_metric_tonnes",                                         # Thousand Metric Tonnes (1,000 mt, macro bulk cargo & mining production)
    "Mmt = 1e6 * mt = Mt = MT_macro = MMT = million_metric_tonnes",                                    # Million Metric Tonnes (1,000 kmt, national trade balances & annual mining output)

    "commodity_short_hundredweight = 100 * pound = cwt = CWT = short_cwt = hundredweight",            # Short Hundredweight (100 lb / 45.3592 kg, CME Live Cattle & agricultural benchmark)
    "commodity_long_hundredweight = 112 * pound = long_cwt = Long_cwt = cwt_uk",                       # Long Hundredweight (112 lb / 50.8023 kg, UK/Commonwealth scrap & imperial bulk)

    "commodity_short_ton = 2000 * pound = short_ton = ston = STON = s_ton = Short_Ton",               # Short Ton (2,000 lb / 907.1847 kg, US domestic steel, coal & scrap benchmark)
    "commodity_long_ton = 2240 * pound = long_ton = lton = LTON = l_ton = Long_Ton",                  # Long Ton (2,240 lb / 1,016.0469 kg, UK/Commonwealth marine bulk & legacy tin)

    # --- Bulk Mining, Iron Ore & Mineral Concentrates (Moisture & Grade Framework) ---
    # Governed by ISO 3082 (Iron ores - Sampling and sample preparation) and commercial smelter
    # contracts. Bulk ores trade on a dry basis (DMT = WMT * [1 - moisture fraction]), whilst contained
    # metal pricing benchmarks are normalised via dry metric tonne units (1 DMTU = 1% metal per tonne).
    "dmt = [dry_bulk_mass] = DMT = dry_metric_tonne = dry_tonne",                                       # Dry Metric Tonne (1,000 kg bone-dry ore mass excluding free moisture)
    "wmt = [wet_bulk_mass] = WMT = wet_metric_tonne = wet_tonne",                                       # Wet Metric Tonne (1,000 kg gross run-of-mine ore mass including free moisture)
    "contained_metal_tonne = [contained_metal_mass] = contained_mt",
    "payable_metal_tonne = [payable_metal_mass] = payable_mt",
    "dmtu = 0.01 * contained_metal_tonne = DMTU = dry_metric_tonne_unit",                                             # Dry Metric Tonne Unit (10 kg pure contained element / 1% Fe per DMT, 62% Fe CFR China benchmark)
    "mtu = 10 * kilogram = MTU = metric_tonne_unit",                                                   # Metric Tonne Unit (10 kg contained WO3/Li2O/Mn, minor metals & ferro-alloys benchmark)
    "stu = 20 * pound = STU = short_ton_unit",                                                         # Short Ton Unit (20 lb / 9.0718 kg pure element / 1% of a short ton, US tungsten & lithium standard)
    "dltu = 22.4 * pound = DLTU = dry_long_ton_unit",                                                 # Dry Long Ton Unit (22.4 lb pure metal / 1% of a long ton, legacy UK/Indian iron ore contracts)

    # --- Precious Bullion, Loco-London Standards & Regional Mandis ---
    # Governed by the London Bullion Market Association (LBMA) Good Delivery Rules, Bureau of Indian
    # Standards (BIS), and Chinese Gold & Silver Exchange Society (CGSE). Based on the troy ounce
    # (480 grains = 31.1034768 g exact), with distinct regional physical trading units.
    "toz = 31.1034768 * gram = troy_ounce = troy_oz = ozt = OZT = TOZ",                               # Troy Ounce (31.1034768 g exact, LBMA/COMEX gold, silver, platinum & palladium benchmark)
    "fine_troy_ounce = toz = ftoz = FTOZ = fine_ounce = fine_toz",                                    # Fine Troy Ounce (Pure 24-karat metal mass = gross toz * assay fineness)
    "tola = 11.6638038 * gram = tola_gold = TOLA",                                                     # Tola (180 troy grains / 11.6638 g, Indian/Pakistani bullion bazaar & MCX physical standard)
    "tael_hk = 37.79936375 * gram = tael = hk_tael = HK_tael = Tael",                                  # Hong Kong Tael (1-1/3 oz avoirdupois / 37.7994 g, CGSE Loco-Hong Kong gold bar standard)
    "tael_china = 37.429 * gram = mainland_tael = shanghai_tael",                                      # Shanghai Tael (37.429 g, traditional Chinese mainland precious metals trading weight)
    "tael_taiwan = 37.5 * gram = taiwan_tael = ryo = Ryo",                                             # Taiwan Tael / Japanese Ryo (37.5 g exact, Taiwanese & Japanese physical gold bullion bar standard)
    "baht_gold = 15.244 * gram = baht_bullion = thai_baht_gold",                                       # Thai Gold Baht (15.244 g @ 96.5% purity, Bangkok Gold Traders Association standard)
    "carat = 200 * milligram = ct = CT = metric_carat = Carat",                                        # Metric Carat (200 mg / 0.2 g, international diamond and gemstone valuation standard)

    # --- Speciality, Minor Metals & Nuclear Energy Fuels ---
    # Unique industry pricing standards for commercial industrial compounds and elemental liquid metals.
    "flask_mercury = 76 * pound = flask = FLASK = mercury_flask = Flask",                              # Flask of Mercury (76 lb / 34.473 kg, European/US benchmark elemental mercury unit)
    "pound_U3O8 = pound = lb_U3O8 = U3O8_lb = U3O8 = u3o8",                                            # Pound of Triuranium Octoxide (UxC / TradeTech spot uranium yellowcake pricing benchmark)

    # --- Battery Metals & Chemical Equivalents ---
    # Universal standard for lithium reserves, chemical supply contracts, and exchange pricing
    # (Fastmarkets, S&P Platts, CME). Normalises spodumene, lithium carbonate, and hydroxide to Li2CO3.
    "tonne_LCE = 1 * mt = LCE = t_LCE = mt_LCE = tonne_lce",                                           # Metric Tonne of Lithium Carbonate Equivalent (Global battery chemical benchmark)

    # --- Troy System Sub-Divisions (LBMA / Assay Standards) ---
    # Exact subdivisions of the troy ounce (480 grains = 1 toz exact).
    "grain = (1 / 480) * toz = gr = grain_troy",                                                       # Troy Grain (64.79891 mg exact, fundamental mass bridge between avoirdupois and troy)
    "pennyweight = 24 * grain = pwt = PWT = dwt_troy = dwt_bullion",                                   # Pennyweight (1/20 toz / 1.55517 g, jewellery, dental alloy & scrap bullion)

    # --- Ferrous Scrap & Steel Mill Imperial Mass ---
    # US domestic scrap steel benchmarks (e.g. No. 1 Heavy Melting Steel [HMS 1/2], Shredded Scrap).
    "gross_ton = commodity_long_ton = GT = gross_ton_steel = gt",                                      # Gross Ton (2,240 lb / 1,016.0469 kg, US scrap steel & metallurgical coke benchmark)

    # =========================================================================
    # Maritime Shipping, Freight Derivatives & Vessel Chartering
    # =========================================================================

    # --- Vessel Capacity, Container Slots & Freight Demand ---
    # Governed by the Baltic Exchange, BIMCO, and International Maritime Organization (IMO)
    # standards. Deadweight tonnage (dwt) reflects total vessel carrying capacity (cargo, bunkers,
    # fresh water, ballast, crew) in metric tonnes. Container freight follows ISO 668 standards
    # based on standard 20-foot and 40-foot intermodal equivalent units (TEU / FEU).
    "dwt = 1e3 * kilogram = DWT = deadweight_ton = deadweight_tonne",                                  # Deadweight Tonnage (1,000 kg vessel cargo, fuel, ballast & stores capacity)
    "TEU = [container_capacity] = teu = twenty_foot_equivalent_unit",                                  # Twenty-Foot Equivalent Unit (ISO 668 standard 20-foot container slot capacity)
    "FEU = 2 * TEU = feu = forty_foot_equivalent_unit",                                                # Forty-Foot Equivalent Unit (2 TEU, standard 40-foot container slot capacity)
    "tonne_nautical_mile = tonne * nautical_mile = tonne_nm = ton_nm",                      # Tonne-Nautical Mile (Bulk shipping transport demand: cargo mass multiplied by sea distance)

    # --- Tanker Freight Benchmarks & Gas Carrier Capacities ---
    # Managed by the Worldscale Association (London and New York). Nominal flat rates in USD/mt
    # are reissued annually; market fixtures are negotiated as percentage points of the flat rate (WS).
    "worldscale_point = [tanker_freight_index] = WS = worldscale = ws_point",                          # Worldscale Point (1% of published annual nominal tanker flat rate for a given route)
    "cbm_lng = meter ** 3 = cbm_liquid_lng = m3_lng",                                                  # Cubic Metre of Liquefied Natural Gas (-162°C cryogenic cargo capacity, e.g. 174,000 cbm)
    "cbm_lpg = meter ** 3 = cbm_liquid_lpg = m3_lpg",                                                  # Cubic Metre of Liquefied Petroleum Gas (Pressurised/refrigerated liquid parcel volume)

    # --- Charterparty Hire & Laytime Demurrage Rates ---
    # Governed by standard charterparty forms (NYPE, ASBATANKVOY, BPVOY5, ShellVOY6).
    # Essential for calculating Time Charter Equivalent (TCE) earnings and port laytime demurrage.
    "day_charter = day = charter_day = hire_day",                                                      # Charter Day (24-hour operational hire period for daily TCE earnings settlement)
    "hour_demurrage = hour = demurrage_hour",                                                          # Demurrage Hour (Hourly penalty rate incurred past agreed laytime allowances)

    # =========================================================================
    # Annualised Mining, LNG Liquefaction & Bulk Freight Flow Rates
    # =========================================================================

    # --- Hourly & Daily Industrial Processing Rates (Mines, Mills & Ports) ---
    # Hourly throughput for crushing circuits, conveyor systems, and bulk shiploader nameplates,
    # alongside daily mill run rates for concentrators, smelters, and leach pads.
    "tph = tonne / hour = TPH = tonne_per_hour = t_h",                                                 # Tonnes per Hour (Bulk shiploader, conveyor belt, and crusher circuit capacity)
    "tpd = tonne / day = TPD = tonne_per_day = t_d",                                                   # Tonnes per Day (Standard mining mill throughput and concentrator intake)
    "ktpd = 1e3 * tpd = KTPD = kt_d = ktonne_per_day",                                                 # Thousand Tonnes per Day (1,000 tpd, tier-1 open-pit copper/gold mill capacity)

    # --- Annual Metric Flow Rates (LNG Export Trains & Base Metals) ---
    # Standard annualised capacity units across international LNG liquefaction trains,
    # regasification terminals, alumina refineries, and primary metal smelters.
    "tpa = tonne / trading_year_365 = TPA = tpy = TPY = tonne_per_year = tonnes_per_annum",            # Tonnes per Annum (Baseline annualised metric flow rate)
    "ktpa = 1e3 * tpa = KTPA = ktpy = KTPY = thousand_tonnes_per_annum",                               # Thousand Tonnes per Annum (1,000 tpa, cathode copper & nickel chemical production)
    "mtpa = 1e6 * tpa = MTPA = mtpy = MTPY = million_tonnes_per_annum",                               # Million Tonnes per Annum (1,000 ktpa, global LNG liquefaction & macro bulk capacity)

    # --- Bulk Mining Moisture-Adjusted Annual Capacities (Iron Ore & Bauxite) ---
    # Specific dry-basis and wet-basis annual production volumes governing off-take agreements
    # and port throughput declarations under ISO 3082.
    "dmtpa = dmt / trading_year_365 = DMTPA = dmtpy = DMTPY = dry_metric_tonnes_per_annum",            # Dry Metric Tonnes per Annum (Bone-dry ore production excluding moisture)
    "wmtpa = wmt / trading_year_365 = WMTPA = wmtpy = WMTPY = wet_metric_tonnes_per_annum",            # Wet Metric Tonnes per Annum (Run-of-mine physical shipping volume including moisture)

    # --- Imperial Mining & Quarrying Flow Rates (US Domestic Coal & Aggregates) ---
    # Standard US short-ton flow rates for domestic coal preparation plants, trona refining,
    # and quarrying operations.
    "stpd = commodity_short_ton / day = STPD = short_ton_per_day = st_d",                              # Short Tons per Day (US coal prep plant and domestic mill throughput)
    "stpa = commodity_short_ton / trading_year_365 = STPA = stpy = STPY = short_ton_per_year",         # Short Tons per Annum (US domestic mining and aggregate production capacity)

    # =========================================================================
    # Agricultural Grains, Tropical Softs & Regional Cash Mandis
    # =========================================================================

    # --- Statutory Grain Bushels: Exchange Mass Benchmarks ---
    # Under USDA, CME Group (CBOT), and Kansas City Board of Trade rules, commercial grain
    # contracts settle strictly on statutory mass (weight), not volumetric Winchester bushels
    # (2,150.42 cu in). Binding statutory crop weights prevents severe volume-to-mass conversion errors.
    "bushel_corn = 56 * pound = bu_corn = BU_CORN = bushel_sorghum = bushel_flaxseed",                  # Standard Corn/Sorghum Bushel (56 lb / 25.4012 kg, CBOT Corn standard; 1 mt ~ 39.368 bu)
    "bushel_wheat = 60 * pound = bu_wheat = BU_WHEAT = bushel_soybeans = bu_soy = BU_SOY",             # Standard Wheat/Soybean Bushel (60 lb / 27.2155 kg, CBOT Wheat/Soybean standard; 1 mt ~ 36.744 bu)
    "bushel_barley = 48 * pound = bu_barley = BU_BARLEY",                                              # Standard Barley Bushel (48 lb / 21.7724 kg, 1 mt ~ 45.930 bu)
    "bushel_oats = 32 * pound = bu_oats = BU_OATS",                                                    # Standard Oats Bushel (32 lb / 14.5150 kg, CBOT Oats standard; 1 mt ~ 68.894 bu)
    "bushel_rye = 56 * pound = bu_rye = BU_RYE",                                                      # Standard Rye Bushel (56 lb / 25.4012 kg)

    # --- European & International Grain Trade Benchmarks ---
    # Governed by Euronext (Paris Milling Wheat / Rapeseed) and GAFTA / FOSFA contract rules.
    # The metric quintal serves as the cash quotation standard across France, Germany, and India.
    "hectolitre = 100 * liter = hl = HL = hectoliter",                                                 # Hectolitre (100 litres, European grain specific weight / test weight volume denominator)
    "kg_per_hl = kilogram / hectolitre = kg_hl",                                                       # Kilograms per Hectolitre (European grain density / milling quality benchmark, e.g. 76 kg/hl)
    "quintal = 100 * kilogram = q = metric_quintal = Quintal",                                         # Metric Quintal (100 kg, universal European and Indian cash grain pricing standard)
    "spanish_quintal = 46 * kilogram = quintal_ag = qq = QQ",                                          # Agricultural Spanish Quintal (100 Castilian lb / 46 kg, Central American coffee/sugar cash trade)

    # --- Tropical Softs Packaging & Commercial Export Bags ---
    # Standardised commercial bags and packing units codified under International Coffee Organization (ICO),
    # Federation of Cocoa Commerce (FCC), and ICE Futures Europe / US contract rules.
    "coffee_bag_60kg = 60 * kilogram = coffee_bag = bag_coffee = sac_coffee_brazil",                   # Coffee Bag 60kg (ICO / ICE Arabica & Robusta standard export bag; 1 mt ~ 16.667 bags)
    "coffee_bag_colombia = 70 * kilogram = bag_coffee_co",                                             # Colombian Federation Coffee Bag (70 kg standard export pack)
    "cocoa_bag_62_5kg = 62.5 * kilogram = cocoa_bag = bag_cocoa",                                     # Cocoa Bag 62.5kg (West African FCC/CMA export standard; exactly 16 bags = 1 mt)
    "sugar_bag_50kg = 50 * kilogram = sugar_bag = bag_sugar",                                         # Refined White Sugar Bag (50 kg, ICE White Sugar No. 407 physical parcel unit; 20 bags = 1 mt)

    # --- Cotton Bales, Asian Cash Mandis & Latin American Mass ---
    # Governed by ICE Cotton No. 2 contract specifications, B3 Brazil standards, and South Asian
    # agricultural mandi regulatory boards.
    "cotton_bale = 480 * pound = cotton_bale_net = bale_cotton = Cotton_Bale",                         # Cotton Net Export Bale (480 lb net / 217.724 kg, ICE Cotton No. 2 contract basis)
    "cotton_bale_gross = 500 * pound = cotton_bale_gr",                                                # Cotton Gross Bale (500 lb gross including bagging and wire ties / 226.796 kg)
    "candy = 355.62 * kilogram = khandi = candy_cotton",                                               # Indian Cotton Candy (355.62 kg / ~784 lb, MCX / domestic Indian physical lint pricing unit)
    "arroba = 15 * kilogram = arroba_beef = ARROBA = arroba_net",                                      # Brazilian Arroba (15 kg net carcass mass, B3 Live Cattle & domestic raw cotton benchmark)
    "maund_metric = 40 * kilogram = metric_maund = mnd_metric",                                        # Metric Maund (40 kg exact, modern Indian Mandi statutory delivery standard)
    "maund_imperial = 82.2857 * pound = maund = mnd = imperial_maund",                                  # Imperial Maund (82.2857 lb / ~37.3242 kg, legacy Pakistani/North Indian ag trade)
    "kati = 0.6 * kilogram = catty = Kati = Catty",                                                    # Southeast Asian Kati / Catty (0.6 kg / 600 g, domestic palm oil and rubber cash weight)
    "picul = 100 * kati = pikul = Picul",                                              # Southeast Asian Picul (60 kg, Malaysian / Indonesian palm oil, rubber & tin cash unit)

    # --- Legacy Imperial Volumetric Agricultural Measures ---
    # Dry capacity measures derived from the Winchester bushel (2,150.42 cu in / 35.2391 litres).
    "commodity_peck = 0.25 * bushel = peck = Peck = pk",                                               # Imperial/US Dry Peck (0.25 dry bushel / 2 dry gallons / 8.8098 litres)

    # =========================================================================
    # Forestry, Timber & Structural Building Materials
    # =========================================================================

    # --- Lumber & Wood Products (CME Group / Physical Yard Standards) ---
    # Governed by the American Softwood Lumber Standard (NIST PS 20) and CME Lumber (LBR)
    # contract rules. 1 board foot represents the nominal volume of a 1 ft x 1 ft x 1 in board
    # (1/12 cu ft / 144 cu in / ~0.00236 m3). Quoted in dollars per thousand board feet ($/MBF).
    # Pulp basis labels each denote one physical tonne. Fibre/moisture basis
    # changes require PulpLumberAssay methods, rather than plain .to() scaling.
    "admt = metric_ton = ADMT = air_dry_metric_tonne",
    "bdmt = metric_ton = BDMT = bone_dry_metric_tonne = odmt = ODMT",
    "board_foot = foot ** 2 * inch = bf = BF = fbm = FBM = board_feet",                     # Board Foot (Nominal timber volume: 144 cu in / 2.3597 litres, CME Lumber basis)
    "thousand_board_feet = 1e3 * board_foot = mbf = MBF = kbf = kFBM",                                 # Thousand Board Feet (1,000 bf, standard wholesale and CME trading lot unit)
    "million_board_feet = 1e6 * board_foot = MMbf = MMBF = MMfbm",                                     # Million Board Feet (1,000 MBF, sawmill annual capacity & commercial forestry harvest)
    "cord_wood = 128 * foot ** 3 = cord = Cord = firewood_cord",                                       # Cord of Wood (128 cu ft / 3.6246 m3 stacked roundwood / firewood standard volume)

    # =========================================================================
    # Environmental Attributes, Compliance Emissions & Carbon Offsets
    # =========================================================================

    # --- Greenhouse Gas Mass Equivalents: IPCC / UNFCCC Basis ---
    # Codified under UNFCCC reporting guidelines and IPCC Fifth/Sixth Assessment Reports (AR5/AR6).
    # Metric tonnes of direct CO2 and Global Warming Potential (GWP)-weighted greenhouse gas equivalents
    # (covering CH4, N2O, SF6, HFCs, PFCs, and NF3 normalised over a 100-year timescale).
    "tonne_CO2 = tonne = tCO2 = t_CO2 = MT_CO2 = mt_CO2",                                               # Metric Tonne of Carbon Dioxide (1,000 kg direct fossil CO2 mass benchmark)
    "tonne_CO2e = tonne = tCO2e = t_CO2e = MT_CO2e = mt_CO2e",                                         # Metric Tonne of CO2 Equivalent (1,000 kg GWP-adjusted greenhouse gas emissions)

    # Compliance instruments have independent dimensions. Surrender rules,
    # vintages and eligibility require an explicit instrument-level policy.
    "EU_allowance = [instrument_EU_allowance] = EUA = eua = eu_allowance = EUAs",                                  # European Union Allowance (EU ETS Directive 2003/87/EC Phase 4 compliance instrument, 1 tCO2e)
    "EU_aviation_allowance = [instrument_EU_aviation_allowance] = EUAA = euaa = eu_aviation_allowance",                    # EU Aviation Allowance (EU ETS aircraft operator compliance unit; merged into EUA Phase 4)
    "UK_allowance = [instrument_UK_allowance] = UKA = uka = uk_allowance = UKAs",                                  # UK Allowance (UK ETS Statutory Instrument 2020 No. 1265 post-Brexit compliance unit, 1 tCO2e)
    "California_carbon_allowance = [instrument_California_carbon_allowance] = CCA = cca = california_allowance = CCAs",           # California Carbon Allowance (CARB Cap-and-Trade / WCI linked market compliance unit, 1 tCO2e)
    "ACCU = [instrument_ACCU] = accu = ACCUs = australian_carbon_credit_unit",                             # Australian Carbon Credit Unit (Clean Energy Regulator Safeguard Mechanism compliance unit, 1 tCO2e)
    "NZU = [instrument_NZU] = nzu = NZUs = new_zealand_unit",                                             # New Zealand Unit (NZ ETS Climate Change Response Act statutory compliance unit, 1 tCO2e)
    "CEA = [instrument_CEA] = cea = CEAs = china_emission_allowance",                                     # China Emission Allowance (National ETS compliance quota, Shanghai Environment and Energy Exchange)
    "KAU = [instrument_KAU] = kau = KAUs = korea_allowance_unit",                                         # Korea Allowance Unit (K-ETS statutory allocation unit, Korea Exchange - KRX)

    # --- Sub-National Short-Ton Compliance (RGGI) ---
    # The Regional Greenhouse Gas Initiative (US Northeast power sector cap-and-invest programme)
    # is legally denominated in US short tons of CO2 rather than metric tonnes (1 short ton = 0.907185 tCO2e).
    "RGGI_allowance = [instrument_RGGI_allowance] = RGGI = rggi = rggi_allowance = RGGIs",                  # RGGI Allowance (US Northeast cap-and-invest power sector compliance unit, 1 short ton CO2)

    # --- Voluntary Carbon Markets (VCM) & Article 6 Mitigation Instruments ---
    # Project-based carbon crediting mechanisms verified against international standards
    # (Verra VCS, Gold Standard, American Carbon Registry, Climate Action Reserve) and UNFCCC Article 6.
    "VCU = [instrument_VCU] = vcu = VCUs = verified_carbon_unit",                                         # Verified Carbon Unit (Verra Verified Carbon Standard [VCS] offset credit, 1 tCO2e)
    "GS_VER = [instrument_GS_VER] = gs_ver = GS_VERs = gold_standard_ver = VER = ver",                      # Gold Standard Verified Emission Reduction (Gold Standard Foundation carbon credit, 1 tCO2e)
    "ACR_offset = [instrument_ACR_offset] = ACR = acr = american_carbon_registry = ACRs",                        # American Carbon Registry Emission Reduction Measurement (Winrock International standard, 1 tCO2e)
    "CAR_CRT = [instrument_CAR_CRT] = CRT = crt = climate_reserve_tonne = CRTs",                              # Climate Reserve Tonne (Climate Action Reserve North American offset protocol, 1 tCO2e)
    "ITMO = [instrument_ITMO] = itmo = ITMOs = article_6_itmo",                                            # Internationally Transferred Mitigation Outcome (UNFCCC Paris Agreement Article 6.2/6.4 unit)
    "CORSIA_eligible_unit = [instrument_CORSIA_eligible_unit] = ceu = CEU = CEUs = CORSIA",                                # CORSIA Eligible Unit (ICAO global aviation carbon offsetting and reduction scheme, 1 tCO2e)

    # --- US EPA Renewable Fuel Standard (RFS2) D-Codes ---
    # Mandated under 40 CFR Part 80. Renewable Identification Numbers (RINs) are generated
    # per physical gallon of renewable fuel based on its statutory Equivalence Value (EV).
    "RIN = [renewable_identification_number] = rin = rins = RFS_RIN",                                  # Renewable Identification Number (US EPA RFS2 compliance credit, 1 ethanol-equivalent gallon)
    "RIN_D3 = RIN = rin_d3 = cellulosic_biofuel_rin",                                                 # D3 Cellulosic Biofuel RIN (Cellulosic ethanol, renewable natural gas [RNG])
    "RIN_D4 = RIN = rin_d4 = biomass_based_diesel_rin",                                               # D4 Biomass-Based Diesel RIN (Biodiesel, Renewable Diesel; EV 1.5–1.7 RINs/gal)
    "RIN_D5 = RIN = rin_d5 = advanced_biofuel_rin",                                                    # D5 Advanced Biofuel RIN (Sugarcane ethanol, non-corn starch biofuels)
    "RIN_D6 = RIN = rin_d6 = renewable_fuel_rin",                                                      # D6 Conventional Renewable Fuel RIN (Corn starch ethanol)

    # --- Sub-National Clean Fuel Standards & European GHG Mandates ---
    # Performance-based lifecycle carbon-intensity compliance mechanisms.
    "LCFS_credit = [instrument_LCFS_credit] = lcfs = LCFS = carb_lcfs = clean_fuel_credit",                      # LCFS Credit (California CARB / Oregon CFP deficit offset, 1 tCO2e reduction)
    "THG_quota = [instrument_THG_quota] = thg = THG = thg_quote = ghg_quota",                                  # German THG-Quote (38. BImSchV transport emission compliance unit, 1 tCO2e reduction)

    # =========================================================================
    # Bespoke, Desk-Specific & Proprietary Commodity Extensions
    # =========================================================================

    # --- Proprietary Desk Overrides & User-Defined Extensions ---
    # Injected tuple for environment-specific, OTC counterparty-bespoke, or proprietary
    # firm-internal units (e.g. internal pipeline zone thermal factors, non-standard
    # storage terminal tank parcels, or proprietary crack spread formulations).
    # =========================================================================
    # Weather Derivatives & Thermal Degree Days
    # =========================================================================
    # Temperature differences accumulated over time, not affine temperatures.
    # HDD/CDD aliases identify index units; baseline and daily observation logic
    # belong in the weather model rather than the registry.
    "degree_day_F = delta_degF * day = degree_day = HDD = CDD = DD = degree_days",
    "degree_day_C = delta_degC * day = degree_day_c = HDD_c = CDD_c = DD_c",
    "degree_hour_F = delta_degF * hour = degree_hour = HDH = CDH",
    "degree_hour_C = delta_degC * hour = degree_hour_c = HDH_c = CDH_c",

    # =========================================================================
    # Power Generation Heat Rates & Spark/Dark Spread Quoting Scales
    # =========================================================================
    # Energy/energy is dimensionless in Pint. Named units retain the quoting
    # scale; they do not enforce fuel-versus-electrical energy semantics.
    "heat_rate_us = btu_it / kilowatt_hour = heat_rate = HR = btu_per_kwh = Btu_kWh",
    "heat_rate_macro = MMBtu_it / megawatt_hour = MMBtu_per_MWh = mmbtu_mwh = MMBtu_MWh",
    "heat_rate_metric = megajoule / kilowatt_hour = MJ_per_kWh = MJ_kWh",
    "heat_rate_gas_turbine = gigajoule / megawatt_hour = GJ_per_MWh = GJ_MWh",
    # Thermal/electrical labels carry the same dimensions. Efficiency is
    # applied explicitly by heat_content.thermal_to_output/output_to_thermal.
    "MWh_thermal = megawatt_hour = MWh_th",
    "MWh_electric = megawatt_hour = MWh_e",
    "MW_thermal = megawatt = MW_th = MWth",
    "MW_electric = megawatt = MW_e = MWe",

    # =========================================================================
    # Hydrogen, Clean Fuel Carriers & Energy Transition
    # =========================================================================
    # Commodity and calorific-basis labels retain the underlying dimensions.
    # Heating-value factors and LHV/HHV transformations belong in conversions.py.
    "kg_H2 = kilogram = kg_h2 = KG_H2",
    "tonne_H2 = 1e3 * kg_H2 = t_H2 = mt_h2",
    "Nm3_H2 = Nm3_0C_101325kPa_dry = nm3_h2",
    "MWh_lhv = megawatt_hour = MWh_LHV = mwh_lhv",
    "MWh_hhv = megawatt_hour = MWh_HHV = mwh_hhv",
    "tonne_ammonia = 1e3 * kilogram = t_NH3 = mt_NH3 = tonne_NH3",
    "gallon_SAF = gallon_petroleum = gal_saf = SAF_gal",

    # =========================================================================
    # Pipeline Operating Pressure Labels
    # =========================================================================
    # Gauge and absolute labels share pressure dimensions. Datum shifts require
    # explicit handling in conversions.py; these aliases apply no pressure offset.
    "psia = psi = PSIA = lb_per_sq_in_abs",
    "psig = psi = PSIG = lb_per_sq_in_gauge",
    "bara = bar = BARA = bar_absolute",
    "barg = bar = BARG = bar_gauge",

    # =========================================================================
    # Clean Fuel Lifecycle Carbon Intensity
    # =========================================================================
    "gCO2e = 1e-6 * tonne_CO2e = g_CO2e = g_co2e",
    "kgCO2e = 1e-3 * tonne_CO2e = kg_CO2e = kg_co2e",
    "CI_unit = gCO2e / megajoule = CI = ci_score = CI_score = gCO2e_per_MJ",
    "kgCO2e_per_MWh = kgCO2e / megawatt_hour = kg_per_MWh_CO2e",

    # =========================================================================
    # Physical Bullion Bars & Nominal Vault Specifications
    # =========================================================================
    # Good Delivery units are nominal reference masses, not actual bar assays.
    "kilobar_gold = kilogram = kilobar = kilo_gold = k_bar",
    "ten_tola_bar = 10 * tola = tt_bar = TT_bar = TT_Bar",
    "good_delivery_gold = 400 * toz = lgd_gold = bar_gold_400oz",
    "good_delivery_silver = 1000 * toz = lgd_silver = bar_silver_1000oz",

    # =========================================================================
    # Grid Dynamics, Battery Storage & Frequency Regulation
    # =========================================================================
    "MW_per_minute = megawatt / minute = MW_min = MW_per_min = ramp_rate_mw",
    "kW_per_minute = kilowatt / minute = kW_min = kW_per_min",
    "C_rate = 1 / hour = c_rate = C_Rate",

    # Crop-weight bushels, not generic dry-volume bushels.
    "bushel_canola = 50 * pound = bu_canola",
    "bushel_rice = 45 * pound = bu_rice",

    # Nuclear-fuel labels are case-sensitive: mtU must never alias mtu (10 kg).
    "SWU = [separative_work] = swu = separative_work_unit",
    "kgU = kilogram = kg_uranium",
    "mtU = 1e3 * kgU = tU = tonne_uranium",

    # Metric pressure labels retain datum metadata for convert_pressure().
    "kpaa = kilopascal = kPa_a = kpa_a = kilopascal_absolute",
    "kpag = kilopascal = kPa_g = kpa_g = kilopascal_gauge",
    "mpaa = megapascal = MPa_a = mpa_a = megapascal_absolute",
    "mpag = megapascal = MPa_g = mpa_g = megapascal_gauge",
    "paa = pascal = Pa_a = pa_a = pascal_absolute",
    "pag = pascal = Pa_g = pa_g = pascal_gauge",

    # Apparent/reactive power labels share electrical power dimensions.
    # Power factor and active/reactive dispatch require explicit grid models.
    "kilovolt_ampere = 1e3 * volt * ampere = kVA",
    "megavolt_ampere = 1e6 * volt * ampere = MVA",
    "kilovolt_ampere_reactive = 1e3 * volt * ampere = kVAR = kVAr",
    "megavolt_ampere_reactive = 1e6 * volt * ampere = MVAR = MVAr",

    # Coverage audit: reuse Pint's SI/prefixed units instead of redefining them.
    # SP 811 B.8/B.9: thermochemical Btu is distinct from IT and 59-F Btu.
    "btu_th = thermochemical_british_thermal_unit",
    # Korean bullion-market convention (not an assay or purity specification).
    "don_korea = 3.75 * gram = korean_don",

    # Explicit concentration bases. Changing fraction bases needs composition.
    "parts_per_billion = 1e-9 = ppb",
    "parts_per_trillion = 1e-12",
    "mass_fraction = [mass_fraction_basis]",
    "mole_fraction = [mole_fraction_basis]",
    "volume_fraction = [volume_fraction_basis]",
    "mass_percent = 0.01 * mass_fraction",
    "mole_percent = 0.01 * mole_fraction",
    "volume_percent = 0.01 * volume_fraction",
    "molarity = mole / liter",
    "normality = [reaction_equivalent_concentration]",

    # Conditional empirical readings: deliberately not fixed density/viscosity factors.
    "degree_Brix = [brix_reading] = degBrix",
    "degree_Plato = [plato_reading] = degPlato",
    "degree_Baume_heavy = [baume_heavy_reading]",
    "degree_Baume_light = [baume_light_reading]",
    "saybolt_universal_second = [saybolt_universal_reading] = SUS",
    "saybolt_furol_second = [saybolt_furol_reading] = SFS",

    # Region-qualified agricultural-area conventions; bare mu is ambiguous in Pint.
    "acre_us_survey = 43560 * survey_foot ** 2",
    "acre_foot_us_survey = acre_us_survey * survey_foot",
    "mu_china = hectare / 15 = chinese_mu",
    "feddan_4200 = 4200 * meter ** 2",
    "tonne_mile = tonne * mile",
    "short_ton_mile = commodity_short_ton * mile",
    "displacement_tonne = metric_ton",
    # IMO tonnage indices require vessel geometry; neither is displacement or mass.
    "gross_tonnage = [gross_tonnage_index]",
    "net_tonnage = [net_tonnage_index]",

    "gigawatt_day = gigawatt * day = GWd",
    # Counts require contract/packing metadata before physical or monetary delivery.
    "item = [item_count]",
    "dozen = 12 * item = doz",
    "api_gravity = [api_gravity_reading] = degree_API",
    "contract = [contract_count]",
    "lot = [lot_count]",
    "delivery_unit = [delivery_count]",
    "crate = [crate_count]",
    "sack = [sack_count]",
    "price_tick = [price_tick_count]",

    *CUSTOM_COMMODITY_UNIT_DEFINITIONS,

    # Examples of bespoke physical parcel and container overrides:
    # "bespoke_iso_tank_24kl = 24000 * liter = iso_tank_24k",                         # 24,000L Chemical/Bunkering ISO Tank Container
    # "flexitank_24kl = 24000 * liter = flexibag",                                    # 24,000L Bulk Liquid Foodgrade/Chemical Flexitank
    # "drum_steel_200l = 200 * liter = drum_200l = standard_drum",                    # 200L / 55 US Gal Standard Industrial Steel Drum
    # "pipeline_batch_nymex = 10000 * oil_barrel = nymex_batch",                      # 10,000 bbl Minimum NYMEX Pipeline Delivery Batch
)


# ---------------------------------------------------------------------------
# Registry Singleton & Serialisation Bindings
# ---------------------------------------------------------------------------

# Only these observed Pint-default collisions are deliberate. New collisions
# (including desk extensions) must fail instead of silently replacing a unit.
_BOOTSTRAP_COLLISIONS: Final[dict[str, frozenset[str]]] = {
    "USD_mill": frozenset({"mil"}), "basis_point": frozenset({"bp", "bps"}),
    "percent": frozenset({"%", "percent"}), "btu_it": frozenset({"BTU", "Btu"}),
    "quad_it": frozenset({"quad"}), "therm_it": frozenset({"th", "therm"}),
    "toe": frozenset({"toe", "tonne_of_oil_equivalent"}),
    "cal_it": frozenset({"cal", "cal_it"}), "bbl": frozenset({"bbl"}),
    "gallon_petroleum": frozenset({"Gal", "gal"}), "mt": frozenset({"tonne"}),
    "kmt": frozenset({"kt"}), "commodity_short_hundredweight": frozenset({"cwt", "hundredweight"}),
    "commodity_short_ton": frozenset({"short_ton"}), "commodity_long_ton": frozenset({"long_ton"}),
    "toz": frozenset({"ozt", "toz", "troy_ounce"}), "carat": frozenset({"carat", "ct"}),
    "grain": frozenset({"gr", "grain"}), "pennyweight": frozenset({"pennyweight"}),
    "dwt": frozenset({"dwt"}), "commodity_peck": frozenset({"peck", "pk"}),
    "board_foot": frozenset({"BF", "FBM", "board_feet", "board_foot"}),
}


class _SemanticCompatibility:
    def is_compatible_with(self, other, *contexts, **kwargs):
        if not super().is_compatible_with(other, *contexts, **kwargs):
            return False
        if _SEMANTIC_SCOPE.get():
            return True
        registry = self._REGISTRY
        source = self.units if isinstance(self, pint.Quantity) else self
        if isinstance(other, str):
            target = registry.parse_expression(other).units
        elif isinstance(other, pint.Quantity):
            target = other.units
        elif isinstance(other, pint.Unit):
            target = other
        else:
            target = registry.dimensionless
        return registry._semantic_tags(source._units) == registry._semantic_tags(target._units)


class _SemanticUnit(_SemanticCompatibility, pint.Unit):
    pass


class _SemanticQuantity(_SemanticCompatibility, pint.Quantity):
    pass


class _GuardedUnitRegistry(pint.UnitRegistry):
    """Guard ordinary Pint APIs as well as the package's public wrappers."""
    _sealed = False
    Unit = _SemanticUnit
    Quantity = _SemanticQuantity

    def __setattr__(self, name, value):
        if self.__dict__.get("_sealed", False) and not name.startswith("_"):
            raise RuntimeError("UnitRegistry is frozen; conversion policies cannot be reassigned.")
        super().__setattr__(name, value)

    def _mutable(self):
        if self._sealed:
            raise RuntimeError("UnitRegistry is frozen.")

    def define(self, definition):
        with _REGISTRY_LOCK:
            self._mutable()
            result = super().define(definition)
            self.__dict__.get("_semantic_tag_cache", {}).clear()
            return result

    def load_definitions(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            return super().load_definitions(*args, **kwargs)

    def add_context(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            return super().add_context(*args, **kwargs)

    def remove_context(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            return super().remove_context(*args, **kwargs)

    def get_name(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().get_name(*args, **kwargs)

    def get_dimensionality(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().get_dimensionality(*args, **kwargs)

    def parse_units(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().parse_units(*args, **kwargs)

    def parse_expression(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().parse_expression(*args, **kwargs)

    def _semantic_tags(self, container):
        cache = self.__dict__.setdefault("_semantic_tag_cache", {})
        if container in cache:
            return cache[container]
        tags = {}
        def walk(name, exponent, seen):
            if name in seen or name not in self._units:
                return
            definition = self._units[name]
            canonical = definition.name
            basis = {"wmt": "wet", "dmt": "dry", "dmtu": "contained", "contained_metal_tonne": "contained", "payable_metal_tonne": "payable", "admt": "air_dry", "bdmt": "bone_dry"}.get(canonical)
            pressure = canonical.lower()
            if basis:
                tags[basis] = tags.get(basis, 0) + exponent
                return
            if pressure in {"psig", "barg", "kpag", "mpag", "pag"}:
                tags["gauge"] = tags.get("gauge", 0) + exponent
                return
            for child, power in definition.reference.items():
                walk(child, exponent * power, seen | {name})
        for name, exponent in container.items():
            walk(name, exponent, set())
        result = {key: value for key, value in tags.items() if value}
        cache[container] = result
        return result

    def _project_mass_units(self, container):
        """Internal physical scaling for an already-authorized assay operation."""
        factor, roots = self._get_root_units(container)
        projected = dict(roots)
        for name in ("wmt", "dmt", "contained_metal_tonne", "payable_metal_tonne"):
            exponent = projected.pop(name, 0)
            if exponent:
                factor *= 1000 ** exponent
                projected["kilogram"] = projected.get("kilogram", 0) + exponent
        return factor, self.UnitsContainer({name: power for name, power in projected.items() if power})

    def _convert_mass_scale(self, value, src, dst):
        source_factor, source = self._project_mass_units(src)
        target_factor, target = self._project_mass_units(dst)
        number = Decimal(str(source_factor)) if isinstance(value, Decimal) else float(source_factor)
        result = super()._convert(value * number, source, target)
        divisor = Decimal(str(target_factor)) if isinstance(value, Decimal) else float(target_factor)
        return result / divisor

    def _convert(self, value, src, dst, *args, **kwargs):
        with _REGISTRY_LOCK:
            if not _SEMANTIC_SCOPE.get() and self._semantic_tags(src) != self._semantic_tags(dst):
                raise pint.DimensionalityError(src, dst, extra_msg="; pressure datum or mass-basis changes require a validated conversion/assay API")
            if _SEMANTIC_SCOPE.get():
                semantic = {"wet", "dry", "contained", "payable"}
                source_tags = set(self._semantic_tags(src)) & semantic
                target_tags = set(self._semantic_tags(dst)) & semantic
                if source_tags and target_tags and source_tags != target_tags:
                    raise pint.DimensionalityError(src, dst, extra_msg="; use explicit moisture/grade/payability conversion")
                if source_tags or target_tags:
                    return self._convert_mass_scale(value, src, dst)
            return super()._convert(value, src, dst, *args, **kwargs)


def _new_registry() -> pint.UnitRegistry:
    registry = _GuardedUnitRegistry(autoconvert_offset_to_baseunit=True, on_redefinition="ignore")
    for definition in COMMODITY_UNIT_DEFINITIONS:
        before = dict(registry._units)
        registry.define(definition)
        changed = {name for name, previous in before.items() if registry._units.get(name) is not previous}
        canonical = definition.split("=", 1)[0].strip()
        unexpected = changed if definition in CUSTOM_COMMODITY_UNIT_DEFINITIONS else changed - _BOOTSTRAP_COLLISIONS.get(canonical, frozenset())
        if unexpected:
            raise pint.errors.RedefinitionError(", ".join(sorted(unexpected)), "unit")
    registry._on_redefinition = "raise"
    return registry


@cache
def get_registry() -> pint.UnitRegistry:
    registry = _new_registry()
    pint.set_application_registry(registry)
    return registry


def _locked_registry(function):
    @wraps(function)
    def locked(*args, **kwargs):
        with _REGISTRY_LOCK:
            return function(*args, **kwargs)
    return locked


# Module-level singleton instance and Quantity constructor export
ureg: Final[pint.UnitRegistry] = get_registry()
Q_ = ureg.Quantity

# ---------------------------------------------------------------------------
# Parsing and Normalisation Infrastructure
# ---------------------------------------------------------------------------

def _has_special_symbol_chars(symbol: str) -> bool:
    """Returns True if symbol contains characters other than letters, numbers, or underscores."""
    for char in symbol: # For every character in the symbol
        if not char.isalnum() and char != "_": # Check if it's non-acceptable characters
            return True # Return True if it is
    return False # Return False if all characters are valid identifier characters

_SAFE_SYMBOLS: Final[dict[str, str]] = {} # Initialise empty dictionary

for symbol, code in CURRENCY_SYMBOL_TO_CODE.items(): # Unpacks pairs from master currency dictionary (e.g. "$", "USD")
    if _has_special_symbol_chars(symbol) or symbol in {"kr", "Rs", "Fr", "R"}: # True for raw currency glyphs (e.g. $, £, C$), False for alphanumeric identifiers (USD, GBP)
        _SAFE_SYMBOLS[symbol] = code # Symbol code pair into _SAFE_SYMBOLS

_SORTED_SAFE_SYMBOLS: tuple[str, ...] = tuple(sorted(_SAFE_SYMBOLS.keys(), key=len, reverse=True)) # Sorts symbols by length descending to match multi-character prefixes before single glyphs

_CURRENCY_CODES_SET: Final[set[str]] = set(CURRENCY_CODES)
_CURRENCY_CANONICAL_CODES: Final[dict[str, str]] = {code: code for code in CURRENCY_CODES}

_CURRENCY_DEFINITION_LINES: Final[frozenset[str]] = (frozenset(CURRENCY_UNIT_DEFINITIONS) | frozenset(CURRENCY_MINOR_UNIT_DEFINITIONS)) #Composite set that merges all primary currency dimensional definitions with their corresponding fractional minor unit definitions

#translation table mapping colloquial trading screen shorthands, sub-unit glyphs, and vendor ticker conventions to their canonical Pint minor currency units
minor_mapping: Final[dict[str, str]] = {
    "zac": "ZAR_cent",
    "zarc": "ZAR_cent",
    "c": "USD_cent",
    "¢": "USD_cent",
    "usc": "USD_cent",
    "us¢": "USD_cent",
    "p": "GBP_penny",
    "penny": "GBP_penny",
    "gbx": "GBP_penny",
    "€c": "EUR_cent",
    "eurc": "EUR_cent",
    "ac": "AUD_cent",
    "a¢": "AUD_cent",
    "cc": "CAD_cent",
    "c¢": "CAD_cent",
}


def _register_minor_unit_definition(definition: str) -> None:
    """Populates minor_mapping from a Pint currency minor-unit definition string."""
    parts = [part.strip() for part in definition.split("=")] #Splits input definition string e.g. "USD_cent = 0.01 * USD = US_cent = cent", parts becomes ["USD_cent", "0.01 * USD", "US_cent", "cent"]
    if not parts or not parts[0]: #evaluates whether the list is empty or if the leading element is empty string
        return #terminates execution if the definition contains no valid primary identifier
    canonical_identifier = parts[0] #Assigns first segment to the variable canonical
    minor_mapping[canonical_identifier] = canonical_identifier #Self-referencing key-value entry into minor_mapping
    minor_mapping[canonical_identifier.lower()] = canonical_identifier #Inserts lowercase representation of the canonical name
    for alias in parts[2:]: #Iterates over all segments from index 2 onward skipping 1 intentionally
        if not alias or not alias.isidentifier():  #Evaluates whethoer alias is empty or fails isidentifier
            continue
        minor_mapping[alias] = canonical_identifier #Maps exact-cased alias to canonical minor unit identifier
        if alias.upper() not in _CURRENCY_CODES_SET: #Checks whether the uppercase form of the alias colides with sovereign currency code
            minor_mapping[alias.lower()] = canonical_identifier #Maps lowercase form of the alias

# Establishing the two foundational lookup tables required to decompose and normalise <numerator>/<denominator> broker market quotes
for _minor_definition in CURRENCY_MINOR_UNIT_DEFINITIONS:
    _register_minor_unit_definition(_minor_definition)
irregular_denominators: Final[dict[str, str]] = {}

def _register_commodity_denominator_definition(definition: str) -> None:
    """This function extract valid unit identifiers and aliases from Pint definition strings
    (e.g. COMMODITY_UNIT_DEFINITIONS) and populates irregular_denominators so any commodity
    product denominator dynamically maps to its canonical unit"""
    parts = [part.strip() for part in definition.split("=")]  # Splits definition into components across '=' signs and trims whitespace
    if not parts or not parts[0]:  # Guard clause: returns early if definition string is empty or missing a leading token
        return  # Aborts processing for invalid or blank definition strings
    canonical_identifier = parts[0].split(":")[0].strip()  # Isolates base unit token, stripping away any dimensional metadata trailing a colon
    if not canonical_identifier.isidentifier():  # Ensures the extracted canonical token is a syntactically valid Python identifier
        return  # Skips registration if the token contains mathematical operators or invalid characters
    irregular_denominators[canonical_identifier] = canonical_identifier  # Maps exact-cased canonical identifier to itself for O(1) direct lookup
    irregular_denominators.setdefault(canonical_identifier.lower(), canonical_identifier)  # Registers lowercase variant without overwriting existing entries
    for alias in parts[2:]:  # Iterates through aliases starting at index 2, skipping index 1 (the mathematical conversion ratio)
        if not alias or not alias.isidentifier():  # Validates that alias token exists and is a clean alphanumeric identifier
            continue  # Discards blank entries and complex arithmetic expressions containing exponents or operators
        irregular_denominators[alias] = canonical_identifier  # Maps exact-cased alias directly back to canonical identifier
        irregular_denominators.setdefault(alias.lower(), canonical_identifier)  # Preserves case-insensitive lookup for broker screen shorthand inputs

for _commodity_definition in COMMODITY_UNIT_DEFINITIONS:  # Iterates sequentially through every unit definition string in the master commodity definitions tuple
    if _commodity_definition not in _CURRENCY_DEFINITION_LINES:  # Filters out currency and minor currency definitions using the composite set
        _register_commodity_denominator_definition(_commodity_definition)  # Parses and registers valid commodity identifiers and aliases into irregular_denominators

# Canonical market denominator overrides aligned with _QUOTE_REPLACEMENTS and _UNIT_TOKEN_REPLACEMENTS
for _denomination_aliases, _canonical_denomination in (  # Iterates through pairs of colloquial broker aliases and their canonical Pint units
    (("bbl", "BBL", "barrel_oil", "bbls", "barrel", "barrels"), "oil_barrel"),  # Standardises 42-gallon petroleum barrel variants, preventing Pint's dry-barrel collision
    (("gal", "Gal", "us_gal", "gallon_petroleum", "gallons"), "gallon"),  # Maps US liquid gallon screen abbreviations and plurals to standard volume
    (("MMBTU", "mmbtu", "MMBtu", "mmBtu", "MMBtu_IT"), "MMBtu_it"),  # Anchors natural gas benchmark tickers to the International Table energy basis
    (("MWh", "mwh"), "megawatt_hour"),  # Normalises wholesale electricity baseload contract blocks
    (("kWh", "kwh"), "kilowatt_hour"),  # Normalises retail and small-scale power generation units
    (("GWh", "gwh"), "gigawatt_hour"),  # Normalises transmission grid and utility-scale balance units
    (("th", "therm", "Therm", "THERM", "th_it"), "therm_it"),  # Maps retail/pipeline gas therm abbreviations directly to IT therms
    (("MT", "mt", "t", "metric_tonne", "tonne", "Tonne", "TONNE", "tonnes"), "tonne"),  # Normalises metric tonne bulk freight, grain, and metals mass units
    (("bu", "bushel", "bushels"), "bushel"),  # Maps standard agricultural grain contract volume/mass abbreviations
    (("lb", "lbs", "pound", "pounds"), "pound"),  # Maps imperial avoirdupois pound mass tickers and plurals
    (("toz", "oz", "ozt", "OZT", "TOZ", "troy_ounce", "troy_oz"), "toz"),  # Forces generic bullion 'oz' and troy ounce variants to LBMA/COMEX troy ounces
    (("cbm", "CBM", "cbm_oil", "m3_oil", "m3"), "meter ** 3"),  # Expands maritime freight and metric storage volume shorthand into valid algebraic syntax
    (("kg", "kgs"), "kilogram"),  # Maps SI base mass abbreviations and informal plural forms
    (("g", "gm"), "gram"),  # Normalises precious/minor metals gram abbreviations
    (("l", "ltr", "litre", "litres", "liters"), "liter"),  # Normalises metric liquid volume abbreviations across US and UK spellings
):
    for _alias in _denomination_aliases:  # Iterates across each individual alias string defined for the current commodity
        irregular_denominators[_alias] = _canonical_denomination  # Maps exact-cased alias token directly to canonical Pint unit for O(1) lookup
        if _alias.lower() != "mt" or _alias in ("mt", "MT"):  # Case collision guard: prevents mixed-case 'Mt' (megatonnes/million tonnes) from corrupting metric tonnes ('mt')
            irregular_denominators[_alias.lower()] = _canonical_denomination  # Registers lowercase version to enable case-insensitive screen quote matching



# ---------------------------------------------------------------------------
# Dynamic Quote Resolution Engine
# ---------------------------------------------------------------------------

def _resolve_quote_expression(quote: str) -> str | None:
    """Converts a quote suchas "USD/BBL" or "gbx per therm" into a parsable Pint unit expression,
      such as "USD / barrel". If it cannot identify or validate the quote it returns "None" """
    
    clean_quote = quote.strip() #Removes leading/trailing whitespaces from the string

    # Split quote into numerator and denominator across '/' or 'per'
    if "/" in clean_quote: 
        raw_numerator, raw_denomination = clean_quote.split("/", 1) #If "/" is in the quote, split it by "/" and set first section as number, second as denomination
    elif " per " in clean_quote.lower(): #If it uses "per", same thing as "/"
        parts = re.split(r"\s+per\s+", clean_quote, maxsplit=1, flags=re.I) #splits clean_quote wherever it finds the word per with whitespace on both sides, ignoring capitalisation.
        raw_numerator, raw_denomination = parts[0], parts[1] #Sets numerator and denominator
    else:
        return None

    num = raw_numerator.strip() #Strips incase of any whitespaces
    denomination = raw_denomination.strip() #Strips incase of any whitespaces

    # -----------------------------------------------------------------------
    # 1. Resolve Currency Numerator
    # -----------------------------------------------------------------------
    num_norm: str | None = None
    num_lower = num.lower()
    if num in minor_mapping:
        # Exact match in sub-unit mapping (e.g. "gbx" -> "GBP_penny")
        num_norm = minor_mapping[num]
    elif num_lower in minor_mapping:
        # Case-insensitive sub-unit match (e.g. "cents" -> "USD_cent")
        num_norm = minor_mapping[num_lower]
    elif num in CURRENCY_SYMBOL_TO_CODE:
        # Single-glyph or composite currency symbol (e.g. "$" -> "USD", "€" -> "EUR")
        num_norm = CURRENCY_SYMBOL_TO_CODE[num]
    elif num.upper() in _CURRENCY_CODES_SET:
        # Standard ISO-4217 3-letter currency code (e.g. "usd" -> "USD")
        num_norm = _CURRENCY_CANONICAL_CODES[num.upper()]
    elif num == "¥":
        # Disambiguate generic yen/yuan glyph to Japanese Yen by default
        num_norm = "JPY"
    elif num.endswith(("_cent", "_penny", "_fils", "_minor")):
        # Pre-canonicalised minor currency units pass through directly
        num_norm = num
    # Abort resolution if the numerator token cannot be identified
    if num_norm is None:
        return None
    # -----------------------------------------------------------------------
    # 2. Resolve Commodity Denominator
    # -----------------------------------------------------------------------
    # Lookup in irregular denominations table with case-insensitive fallback
    denomination_norm = irregular_denominators.get(denomination, irregular_denominators.get(denomination.lower(), denomination),)
    for pattern, replacement in _UNIT_TOKEN_REPLACEMENTS:
        denomination_norm = pattern.sub(replacement, denomination_norm)
    # Reject injection of invalid characters outside standard Pint identifiers,
    # spaces, and exponent/multiplication operators (e.g. "meter ** 3")
    sanitised_check = (denomination_norm.replace("_", "").replace(" ", "").replace("*", ""))

    if not sanitised_check.isalnum():
        return None
    expression = f"{num_norm} / {denomination_norm}"
    try:
        ureg.Unit(expression)
    except (pint.PintError, ValueError, TypeError):
        return None
    return expression

class QuoteReplacementsMap(dict[str, str]):
    """Dynamic, self-caching dictionary that can recognise quote keys it hasn’t stored yet. 
    It uses _resolve_quote_expression from the previous example to parse and normalise them

    __missing__: Resolves and caches a quote when its key is not already in the dictionary.

    __contains__: Checks whether a quote is stored or can be resolved.

    get: Returns a quote’s value, or a default if it cannot be resolved.
    """

    def __missing__(self, key: str) -> str:
        with _REGISTRY_LOCK: #Locks the registry while it resolves and possibly stores the quote, helping prevent simultaneous lookups from interfering.
            if dict.__contains__(self, key): #Checks whether the key is stored. Another lookup may have added it while this call was waiting for the lock.
                return dict.__getitem__(self, key) #If it is now stored, returns its value directly.
            resolved = _resolve_quote_expression(key) #Tries to parse and normalise the quote key.
            if resolved is not None: #Checks whether parsing succeeded.
                self[key] = resolved #Stores the resolved expression in the dictionary so it can be reused.
                return resolved #Returns the resolved expression.
            raise KeyError(key) 

    def __contains__(self, key: object) -> bool:
        if super().__contains__(key): #Checks whether the key is already stored in the dictionary.
            return True #Returns True if it is stored.
        return isinstance(key, str) and _resolve_quote_expression(key) is not None 
        #Otherwise, checks that the key is a string and tries to resolve it. 
        #Returns True if it can be resolved. This check does not store the result.

    def get(self, key: str, default: Any = None) -> Any:
        try: #Looks up the key. If it is missing, this invokes __missing__, which tries to resolve and cache it.
            return self[key]
        except KeyError:
            return default #Returns the supplied default value instead of raising KeyError.


# Master lookup map populated with baseline institutional benchmarks
_QUOTE_REPLACEMENTS: Final[QuoteReplacementsMap] = QuoteReplacementsMap({
    # Liquids & Refined Products
    "$/BBL": "USD / oil_barrel",
    "$/bbl": "USD / oil_barrel",
    "$/gal": "USD / gallon",
    "c/gal": "USD_cent / gallon",
    "cents/gal": "USD_cent / gallon",
    # Gas & Power
    "$/MMBTU": "USD / MMBtu_it",
    "$/mmbtu": "USD / MMBtu_it",
    "$/MWh": "USD / megawatt_hour",
    "p/th": "GBP_penny / therm_it",
    "p/therm": "GBP_penny / therm_it",
    "£/MWh": "GBP / megawatt_hour",
    "€/MWh": "EUR / megawatt_hour",
    # Agriculture & Freight
    "$/MT": "USD / tonne",
    "$/mt": "USD / tonne",
    "c/bu": "USD_cent / bushel",
    "cents/bu": "USD_cent / bushel",
    "c/lb": "USD_cent / pound",
    "cents/lb": "USD_cent / pound",
    "$/lb": "USD / pound",
    # Bullion & Precious Metals
    "$/toz": "USD / toz",
    "$/oz": "USD / toz",
    # Compliance Carbon & Emissions
    "€/t": "EUR / tonne_CO2",
    "€/EUA": "EUR / EU_allowance",
    "£/t": "GBP / tonne_CO2",
    "£/UKA": "GBP / UK_allowance",
    "$/CCA": "USD / California_carbon_allowance",
})


# ---------------------------------------------------------------------------
# Pre-Parsing Token Interventions & Regex Matchers
# ---------------------------------------------------------------------------
# Strict negative lookbehind/lookahead patterns to bypass Pint's prefix engine
# without corrupting compound variable names (e.g. daily_bbl_rate)
_UNIT_TOKEN_REPLACEMENTS: Final[list[tuple[re.Pattern[str], str]]] = [
    (re.compile(r"(?<![A-Za-z0-9_])psia(?![A-Za-z0-9_])", re.I), "psia"),
    (re.compile(r"(?<![A-Za-z0-9_])psig(?![A-Za-z0-9_])", re.I), "psig"),
    # Composite expressions must precede their individual token replacements.
    (re.compile(r"(?<![A-Za-z0-9_])MMBtu/MWh(?![A-Za-z0-9_])", re.I), "MMBtu_it / megawatt_hour"),
    (re.compile(r"(?<![A-Za-z0-9_])Btu/kWh(?![A-Za-z0-9_])", re.I), "btu_it / kilowatt_hour"),
    (re.compile(r"(?<![A-Za-z0-9_])bbl(?![A-Za-z0-9_])", re.I), "oil_barrel"),
    (re.compile(r"(?<![A-Za-z0-9_])MMbbl(?![A-Za-z0-9_])", re.I), "MMbbl"),
    (re.compile(r"(?<![A-Za-z0-9_])Mbbl(?![A-Za-z0-9_])", re.I), "kbbl"),
    (re.compile(r"(?<![A-Za-z0-9_])MMBtu(?![A-Za-z0-9_])", re.I), "MMBtu_it"),
    (re.compile(r"(?<![A-Za-z0-9_])MMcf(?![A-Za-z0-9_])", re.I), "MMscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Mcf(?![A-Za-z0-9_])", re.I), "Mscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Bcf(?![A-Za-z0-9_])", re.I), "Bscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Tcf(?![A-Za-z0-9_])", re.I), "Tscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])MMboe(?![A-Za-z0-9_])", re.I), "MMboe"),
    (re.compile(r"(?<![A-Za-z0-9_])Mboe(?![A-Za-z0-9_])", re.I), "kboe"),
    (re.compile(r"(?<![A-Za-z0-9_])cwt(?![A-Za-z0-9_])", re.I), "commodity_short_hundredweight"),
    (re.compile(r"(?<![A-Za-z0-9_])cbm(?![A-Za-z0-9_])", re.I), "meter ** 3"),
]

# Require digits and valid three-digit comma groups, retaining decimal/exponent support.
#This builds a regular expression for matching a numeric magnitude, such as 1234, 1,234.50 or 2.5e-3.
_MAGNITUDE_PATTERN = r"(?:(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


#Builds and stores a regular expression for matching quotes that start with a 
#currency symbol, such as -$1,234.50/bbl or ($1,234.50/bbl). 
#It prepares the matcher; it does not match a quote itself
# Matches "$1,250.50/bbl", "($45.20/bbl)", and "$.50/gal".
# Longer currency symbols must precede their single-character suffixes.
def _recompile_prefix_currency_pattern() -> None:
    """Refresh prefix matching after symbol additions; caller holds the registry lock."""
    global _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
    symbols = tuple(sorted(_SAFE_SYMBOLS, key=len, reverse=True)) #Takes the symbols in _SAFE_SYMBOLS sorting them from longest to shortest and stores the rorted results in tuple called symbols
    #Compiles a regex that splits a currency quote into its optional brackets, sign, symbol amount and denomination
    pattern = re.compile(
        rf"^\s*(?P<bracket>\()?\s*(?P<sign>[+-]?)\s*"
        rf"(?P<symbol>{'|'.join(re.escape(s) for s in symbols)})\s*"
        rf"(?P<magnitude>{_MAGNITUDE_PATTERN})\s*"
        r"(?:/\s*(?P<denomination>.+?))?\s*(?(bracket)\))\s*$"
    )
    _SORTED_SAFE_SYMBOLS = symbols
    _PREFIX_CURRENCY_PATTERN = pattern

_recompile_prefix_currency_pattern()

# Matches "1,500 MWh", "(2,500) bbl", "2.5%", and "15 bps".
_NUMERIC_MAGNITUDE_PATTERN = re.compile(
    r"^\s*(?P<bracket>\()?\s*"
    rf"(?P<magnitude>[+-]?{_MAGNITUDE_PATTERN})"
    r"(?(bracket)\s*\))"
    r"(?:\s+|(?=[^\d\s,.)+\-]))(?P<unit>.+?)\s*$"
)

# Separate batch magnitudes must not silently discard an embedded unit.
_BARE_MAGNITUDE_PATTERN = re.compile(
    r"^\s*(?P<bracket>\()?\s*"
    rf"(?P<magnitude>[+-]?{_MAGNITUDE_PATTERN})"
    r"\s*(?(bracket)\))\s*$"
)


# Symbols that can represent several currencies.
# Add every shared glyph or market alias your registry accepts.
_SHARED_SYMBOL_CURRENCIES: dict[str, frozenset[str]] = {
    "$": frozenset({
        "USD", "ARS", "AUD", "BBD", "BMD", "BND", "BSD", "BZD",
        "CAD", "CLP", "COP", "CUP", "FJD", "GYD", "HKD", "JMD",
        "KYD", "LRD", "MXN", "NAD", "NZD", "SGD", "SRD", "TTD",
        "TWD", "UYU", "XCD", "WST",
    }),
    "¥": frozenset({"JPY", "CNY"}),
    "£": frozenset({"GBP", "FKP", "GIP", "SHP", "SSP"}),
    "kr": frozenset({"SEK", "NOK", "DKK", "ISK"}),
    "Rs": frozenset({"INR", "PKR", "LKR", "NPR", "MUR", "SCR"}),
    "\u20a8": frozenset({"INR", "PKR", "LKR", "NPR", "MUR", "SCR"}),
    "Fr": frozenset({"CHF", "CDF", "DJF", "GNF", "RWF", "BIF", "KMF", "XAF", "XOF", "XPF"}),
    "R": frozenset({"ZAR"}),
    "\u00a5": frozenset({"JPY", "CNY"}),
}

def _resolve_currency_symbol(symbol: str, default_currency: str | None,) -> str:
    """Resolve any registered currency symbol, using context for shared symbols."""
    try:
        code = CURRENCY_SYMBOL_TO_CODE[symbol]
    except KeyError:
        raise ValueError(f"Unknown currency symbol: {symbol!r}") from None
    
    if default_currency is None:
        return code
    if default_currency.upper() not in _CURRENCY_CANONICAL_CODES:
        raise ValueError(f"Unknown ISO 4217 currency code: {default_currency}")
    
    default_currency = _CURRENCY_CANONICAL_CODES[default_currency.upper()]
    valid_codes = _SHARED_SYMBOL_CURRENCIES.get(symbol, frozenset({code}))
    
    if default_currency not in valid_codes:
        raise ValueError(f"{symbol!r} is not a valid currency symbol for {default_currency}")
    return default_currency
# ---------------------------------------------------------------------------
# Public Sanitisation & Ingestion APIs
# ---------------------------------------------------------------------------

@contextmanager
def sandbox_registry() -> Iterator[pint.UnitRegistry]:
    """Yield an independent registry containing the baseline master definitions.

    Runtime desk additions and master freeze state are not copied. Define test
    units directly on the yielded registry; module APIs still use the master.
    The application registry is never rebound, even when the context raises.
    """
    sandbox = _new_registry()
    yield sandbox


@contextmanager
def _registry_transaction():
    global _REGISTRY_FROZEN, _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
    with _REGISTRY_LOCK:
        _check_frozen()
        state = copy.deepcopy(ureg.__dict__, {id(ureg): ureg, id(ureg.Quantity): ureg.Quantity, id(ureg.Unit): ureg.Unit})
        names = ("CURRENCY_MINOR_DIGITS", "CURRENCY_SYMBOL_TO_CODE", "_SAFE_SYMBOLS", "_CURRENCY_CODES_SET", "_CURRENCY_CANONICAL_CODES", "minor_mapping", "irregular_denominators", "_QUOTE_REPLACEMENTS")
        snapshots = {name: globals()[name].copy() for name in names}
        tokens = _UNIT_TOKEN_REPLACEMENTS[:]
        frozen, symbols, pattern = _REGISTRY_FROZEN, _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
        try:
            yield
        except BaseException:
            ureg.__dict__.clear()
            ureg.__dict__.update(state)
            for name, snapshot in snapshots.items():
                globals()[name].clear()
                globals()[name].update(snapshot)
            _UNIT_TOKEN_REPLACEMENTS[:] = tokens
            _REGISTRY_FROZEN, _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN = frozen, symbols, pattern
            raise


def freeze_registry() -> None:
    """Seal public registration after service startup configuration.

    Register approved currencies and commodity units during process
    initialisation, then call ``freeze_registry()`` before serving requests.
    Freezing is permanent for this process and blocks the public registration
    APIs and the registry's define/load_definitions/add_context/remove_context
    mutators. Reading units and activating an existing context remain available.
    """
    global _REGISTRY_FROZEN
    with _REGISTRY_LOCK:
        _REGISTRY_FROZEN = True
        ureg._sealed = True

def is_registry_frozen() -> bool:
    """Return whether public runtime registration has been disabled."""
    with _REGISTRY_LOCK:
        return _REGISTRY_FROZEN

def _check_frozen() -> None:
    """Check the registration guard while the caller holds _REGISTRY_LOCK."""
    if _REGISTRY_FROZEN:
        raise RuntimeError("UnitRegistry is frozen. Dynamic registration of custom units or currencies is prohibited in locked production environments.")

def register_currency(code: str, *, decimals: int, symbols: tuple[str, ...] = (), aliases: tuple[str, ...] = (), sub_units: dict[str, str] | None = None, quote_aliases: dict[str, str] | None = None,) -> None:
    """
    Registers a fiat, crypto, or central bank digital currency (CBDC) into the registry.

        --- Parameters ---

        code : str
            Canonical currency symbol (e.g. "BTC", "ETH", "USDT").

        decimals : int
            Settlement precision exponent (e.g. 8 for BTC, 18 for ETH, 6 for USDC).
            Used downstream by Decimal Banker's Rounding engines.

        symbols : tuple[str, ...]
            Graphical or compound screen symbols (e.g. a cryptocurrency glyph).

        aliases : tuple[str, ...]
            Alternative screen symbols or ticker tags (e.g. ("XBT", "bitcoin")).

        sub_units : dict[str, str] | None
            Named fractional denominations mapped to their scaling definition relative to the canonical unit.
            e.g. {"satoshi": "1e-8 * BTC = sat = sats"}

        quote_aliases : dict[str, str] | None
            Broker/exchange screen quote normalisations (e.g. {"BTC/USD": "BTC / USD"}).
    """

    cleaned_currency_code_attempt = code.strip() # Remove leading/trailing whitespace from the currency code
    if not cleaned_currency_code_attempt or not cleaned_currency_code_attempt.isalnum(): # If the cleaned code is empty or contains non-alphanumeric characters, raise an error
        raise ValueError(f"Currency code '{code}' must be a non-empty alphanumeric string.")
    if isinstance(decimals, bool) or not isinstance(decimals, int) or decimals < 0:
        raise ValueError(f"Decimals for currency '{code}' must be a non-negative integer, got {decimals}.")

    # Construct dimensional definition isolating this currency
    if aliases: # If there are aliases provided, join them into a string for the definition
        joined_aliases = " = ".join(aliases) # Join aliases with " = " to create a definition string
        alias_attempt = f" = {joined_aliases}" # Create the alias attempt string
    else:
        alias_attempt = "" # If no aliases are provided, set it to an empty string

    base_definition = f"{cleaned_currency_code_attempt} = [currency_{cleaned_currency_code_attempt}]{alias_attempt}" # Construct the base definition string for the currency, including its dimension and any aliases

    with _registry_transaction(): # Roll back partial definitions and metadata on failure
        _check_frozen()
        for symbol in symbols:
            if not isinstance(symbol, str) or not symbol or symbol != symbol.strip() or any(char.isspace() for char in symbol):
                raise ValueError("Currency symbols must be non-empty strings without whitespace.")
            if not _has_special_symbol_chars(symbol) or any(char in symbol for char in "/=()*0123456789"):
                raise ValueError(f"Invalid graphical currency symbol: {symbol!r}")
            if symbol in CURRENCY_SYMBOL_TO_CODE:
                raise ValueError(f"Currency symbol {symbol!r} is already registered.")
        if cleaned_currency_code_attempt in ureg: # If the currency code is already registered in the unit registry, raise an error
            raise ValueError(f"Currency '{cleaned_currency_code_attempt}' is already registered in the unit registry.")

        # 1. Register canonical base currency dimension
        ureg.define(base_definition)

        # 2. Register minor / sub-denominations (e.g., satoshis, gwei, wei)
        if sub_units:
            for sub_name, sub_expressions in sub_units.items():
                sub_definition = f"{sub_name} = {sub_expressions}"
                ureg.define(sub_definition)
                _register_minor_unit_definition(sub_definition)

        # 3. Synchronise downstream settlement engine metadata
        _CURRENCY_CODES_SET.add(cleaned_currency_code_attempt.upper())
        _CURRENCY_CODES_SET.update(alias.upper() for alias in aliases)
        _CURRENCY_CANONICAL_CODES[cleaned_currency_code_attempt.upper()] = cleaned_currency_code_attempt
        for alias in aliases:
            _CURRENCY_CANONICAL_CODES[alias.upper()] = cleaned_currency_code_attempt
        CURRENCY_MINOR_DIGITS[cleaned_currency_code_attempt] = decimals
        for alias in aliases:
            CURRENCY_MINOR_DIGITS[alias] = decimals

        # Publish symbols before the new regex so every regex match is resolvable.
        if symbols:
            for symbol in symbols:
                CURRENCY_SYMBOL_TO_CODE[symbol] = cleaned_currency_code_attempt
                _SAFE_SYMBOLS[symbol] = cleaned_currency_code_attempt
            _recompile_prefix_currency_pattern()

        # 4. Register vendor market data quote substitutions
        if quote_aliases:
            for raw_quote, normalised_expressions in quote_aliases.items():
                _QUOTE_REPLACEMENTS[raw_quote] = normalised_expressions


def register_commodity_unit(definition: str, *, quote_aliases: dict[str, str] | None = None,) -> None:
    """
    Registers a custom or desk-specific commodity unit into the active registry.

    --- Parameters ---

    definition : str
        A valid Pint unit definition string (e.g. "drum_200l = 200 * liter = steel_drum").

    quote_aliases : dict[str, str] | None
        Optional raw broker screen quote substitutions to map to this unit
        (e.g. {"$/drum": "USD / drum_200l"}).
    """

    if not isinstance(definition, str) or not definition.strip(): # If the definition is not a string or is empty after stripping whitespace
        raise ValueError("Unit definition must be a non-empty string.") # Raise an error

    clean_definition = definition.strip() # Extract canonical unit name (everything prior to the first '=' or space)
    canonical_name = clean_definition.split("=")[0].split(":")[0].strip() # Strip whitespace from the canonical name

    # Check if the unit is already defined in the registry
    with _registry_transaction(): # Roll back partial definitions and metadata on failure
        _check_frozen()
        if canonical_name in ureg: # If the name is already in the registry, raise an error
            raise ValueError(f"Unit '{canonical_name}' is already defined in the registry.")

        ureg.define(clean_definition) # Register into active Pint registry
        _register_commodity_denominator_definition(clean_definition)

        # Resolve custom symbols and aliases before Pint can interpret them as
        # prefixed units. Keep matching case-sensitive for desk-defined names.
        for token in clean_definition.split("#", 1)[0].split("=")[2:]:
            token = token.strip()
            if token and token != "_":
                _UNIT_TOKEN_REPLACEMENTS.insert(0, (re.compile(r"(?<![A-Za-z0-9_])" + re.escape(token) + r"(?![A-Za-z0-9_])"), canonical_name,))

        if quote_aliases: # Register vendor screen quote replacements if provided
            for screen_quote, normalised_expressions in quote_aliases.items(): # For each screen quote and its corresponding normalised expression in the provided dictionary
                _QUOTE_REPLACEMENTS[screen_quote] = normalised_expressions # Quote replacements are updated with the new mappings



@_locked_registry
def sanitize_unit_string(raw_unit: str, default_currency: str | None = None) -> str:
    """Normalises raw market quotations and tickers into Pint-compliant unit expressions."""
    text = re.sub(r"\s+per\s+", " / ", raw_unit.strip(), flags=re.I)
    ambiguous = re.search(r"(?<![A-Za-z0-9_])(ppt|mu|ton_mile)(?![A-Za-z0-9_])", text)
    if ambiguous:
        raise ValueError(f"Ambiguous market unit {ambiguous[0]!r}; use an explicit ratio, area or freight convention.")

    # 1. Fast-path lookup against cached quote mappings
    if default_currency is None and text in _QUOTE_REPLACEMENTS:
        return _QUOTE_REPLACEMENTS[text]

    # 2. Extract leading currency glyphs and convert to ISO codes (e.g. "$/barrel" -> "USD/barrel")
    for symbol in _SORTED_SAFE_SYMBOLS:
        # Ensure symbol is standalone or followed immediately by a slash or space
        if text.startswith(symbol) and (len(text) == len(symbol) or text[len(symbol)] in {"/", " "}):
            resolved_code = _resolve_currency_symbol(symbol, default_currency)
            text = resolved_code + text[len(symbol):]
            break

    if text in _QUOTE_REPLACEMENTS:
        return _QUOTE_REPLACEMENTS[text]

    # 3. Normalise high-frequency exchange minor-unit prefixes
    if text.startswith("p/"):
        text = "GBP_penny/" + text[2:]
    elif text.startswith("c/"):
        text = "USD_cent/" + text[2:]
    elif text.startswith("cents/"):
        text = "USD_cent/" + text[6:]

    # 4. Apply strict-boundary regex replacements for problematic tokens (e.g. bbl, MMBtu)
    for pattern, replacement in _UNIT_TOKEN_REPLACEMENTS:
        text = pattern.sub(replacement, text)

    return text


@_locked_registry
def parse_unit(unit: str | pint.Unit, *, default_currency: str | None = None) -> pint.Unit:
    """Parse and normalise a raw unit expression into a validated Pint Unit."""
    if isinstance(unit, pint.Unit):
        return unit
    if not isinstance(unit, str) or not unit.strip():
        raise ValueError(f"Unit expression must be a non-empty string, got {unit!r}")

    return ureg.Unit(sanitize_unit_string(unit, default_currency))


def _parse_magnitude(value, as_decimal):
    if isinstance(value, str):
        match = _BARE_MAGNITUDE_PATTERN.fullmatch(value)
        if match is None:
            raise ValueError(f"Invalid numeric magnitude: {value!r}")
        value = match.group("magnitude").replace(",", "")
        if match.group("bracket"):
            value = "-" + value.lstrip("+-")
    if isinstance(value, bool):
        raise TypeError("Boolean values are not quantity magnitudes.")
    magnitude = Decimal(str(value)) if as_decimal else float(value)
    if not Decimal(str(magnitude)).is_finite():
        raise ValueError("Quantity magnitude must be finite and representable.")
    return magnitude


@_locked_registry
def parse_quantity(value: int | float | Decimal | str, unit: str | None = None, *, default_currency: str | None = None, as_decimal: bool = True,) -> pint.Quantity:
    """Parses numeric expressions or quotation strings into validated Pint Quantities.

    Handles currency-prefixed strings ("$75.50/bbl"), spaced numbers ("75.50 USD/bbl"),
    accounting negatives ("(1,500) MWh"), percentages ("2.5%"), and scientific
    notation ("1.5e6 scfd"). Preserves high-precision arithmetic by
    defaulting to `Decimal` magnitude storage.
    """
    if isinstance(value, str) and unit is None:
        # Pattern A: Leading currency symbol expressions (e.g. "$75.50/bbl", "-£45.20 / MWh")
        prefix_match = _PREFIX_CURRENCY_PATTERN.match(value)
        if prefix_match:
            sign = "-" if prefix_match.group("bracket") else (prefix_match.group("sign") or "")
            symbol = prefix_match.group("symbol")
            mag_str = sign + prefix_match.group("magnitude").replace(",", "")
            denomination = prefix_match.group("denomination")

            currency_code = _resolve_currency_symbol(symbol, default_currency)
            # Keep the screen ticker intact so explicit mappings such as EUR
            # carbon quotes retain their commodity identity on round trips.
            unit_str = f"{symbol}/{denomination}" if denomination else currency_code
            sanitized = sanitize_unit_string(unit_str, default_currency)
            magnitude = _parse_magnitude(mag_str, as_decimal)
            return ureg.Quantity(magnitude, sanitized)

        # Pattern B: Standard numeric value with trailing unit string (e.g. "75.50 USD/bbl")
        num_match = _NUMERIC_MAGNITUDE_PATTERN.match(value)
        if num_match:
            mag_str = num_match.group("magnitude").replace(",", "")
            if num_match.group("bracket"):
                mag_str = "-" + mag_str.lstrip("+-")
            unit_str = num_match.group("unit")
            sanitized = sanitize_unit_string(unit_str, default_currency)
            magnitude = _parse_magnitude(mag_str, as_decimal)
            return ureg.Quantity(magnitude, sanitized)

        if _BARE_MAGNITUDE_PATTERN.fullmatch(value):
            return ureg.Quantity(_parse_magnitude(value, as_decimal), ureg.dimensionless)

        # Pattern C: Standalone dimensional unit strings with implicit unit magnitude (e.g. "USD/bbl")
        parsed = ureg.Quantity(sanitize_unit_string(value, default_currency))
        magnitude = _parse_magnitude(parsed.magnitude, as_decimal)
        return ureg.Quantity(magnitude, parsed.units)

    # Pattern D: Numerical value supplied without an explicit unit (dimensionless quantity)
    if unit is None:
        return ureg.Quantity(_parse_magnitude(value, as_decimal), ureg.dimensionless)

    # Pattern E: Explicit separate value and unit parameters passed
    sanitized = sanitize_unit_string(unit, default_currency)
    return ureg.Quantity(_parse_magnitude(value, as_decimal), sanitized)

# ---------------------------------------------------------------------------
# Institutional Validation & Unit Enforcement Helpers
# ---------------------------------------------------------------------------

def validate_dimensionality(quantity: pint.Quantity, expected_dimension: str | pint.Unit | pint.Quantity,) -> bool:
    """Verifies that a quantity strictly conforms to an expected dimensional context.

    Prevents execution errors such as pricing volumetric barrels against mass-denominated
    hedging curves.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, received {type(quantity).__name__}")

    if isinstance(expected_dimension, (pint.Unit, pint.Quantity)):
        return quantity.check(expected_dimension.dimensionality)

    # Resolve raw string through sanitiser prior to registry dimensionality verification
    parsed_expected = ureg(sanitize_unit_string(expected_dimension))
    return quantity.check(parsed_expected.dimensionality)

def ensure_unit(quantity: pint.Quantity, target_unit: str | pint.Unit,) -> pint.Quantity:
    """Converts a quantity to a specified unit, raising Pint's DimensionalityError on mismatch."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, received {type(quantity).__name__}")
    sanitized_target = (sanitize_unit_string(target_unit) if isinstance(target_unit, str) else target_unit)
    return quantity.to(sanitized_target)


def are_convertible(source_unit: str | pint.Unit | pint.Quantity, target_unit: str | pint.Unit | pint.Quantity, *, default_currency: str | None = None,) -> bool:
    """Return whether ordinary conversion is dimensionally and semantically safe."""
    def resolve_unit(value: str | pint.Unit | pint.Quantity) -> pint.Unit:
        if isinstance(value, pint.Quantity):
            return value.units
        return parse_unit(value, default_currency=default_currency)
    try:
        source, target = resolve_unit(source_unit), resolve_unit(target_unit)
        if source.dimensionality != target.dimensionality:
            return False
        return ureg._semantic_tags(source._units) == ureg._semantic_tags(target._units)
    except (pint.PintError, TypeError, ValueError):
        return False

# ---------------------------------------------------------------------------
# Reverse Market Formatting
# ---------------------------------------------------------------------------

_CANONICAL_TO_MARKET_TICKER: Final[dict[str, str]] = {
    "USD / oil_barrel": "$/bbl",
    "USD / gallon": "$/gal",
    "USD_cent / gallon": "c/gal",
    "USD / MMBtu_it": "$/MMBtu",
    "USD / megawatt_hour": "$/MWh",
    "GBP_penny / therm_it": "p/th",
    "GBP / megawatt_hour": "\u00a3/MWh",
    "EUR / megawatt_hour": "\u20ac/MWh",
    "USD / tonne": "$/mt",
    "USD_cent / bushel": "c/bu",
    "USD_cent / pound": "c/lb",
    "USD / pound": "$/lb",
    "USD / toz": "$/toz",
    "EUR / tonne_CO2": "\u20ac/t",
    "EUR / EU_allowance": "\u20ac/EUA",
    "GBP / tonne_CO2": "\u00a3/t",
    "GBP / UK_allowance": "\u00a3/UKA",
    "USD / California_carbon_allowance": "$/CCA",
    "USD / gigawatt_hour": "$/GWh",
    "EUR / gigawatt_hour": "\u20ac/GWh",
    "GBP / gigawatt_hour": "\u00a3/GWh",
    "USD / kilowatt_hour": "$/kWh",
    "USD_cent / kilowatt_hour": "c/kWh",
    "GBP_penny / kilowatt_hour": "p/kWh",
    "USD / Mscf_60F_14_73psia": "$/Mcf",
    "USD / dth_it": "$/Dth",
    "USD / boe_58": "$/boe",
    "USD / commodity_short_hundredweight": "$/cwt",
    "USD_cent / commodity_short_hundredweight": "c/cwt",
}
# Resolve aliases through Pint once, so its preferred spelling cannot prevent a match.
_MARKET_TICKERS: Final[dict[pint.Unit, str]] = {ureg.Unit(canonical): ticker for canonical, ticker in _CANONICAL_TO_MARKET_TICKER.items()}

_CANONICAL_UNIT_DISPLAY: Final[dict[str, str]] = {
    "gCO2e / megajoule": "CI",
    "CI_unit": "CI",
    # Render complete ratios before replacing individual energy identifiers.
    "btu_it / kilowatt_hour": "Btu/kWh",
    "MMBtu_it / megawatt_hour": "MMBtu/MWh",
    "heat_rate_us": "Btu/kWh",
    "heat_rate_macro": "MMBtu/MWh",
    "heat_rate_metric": "MJ/kWh",
    "heat_rate_gas_turbine": "GJ/MWh",
    "oil_barrel": "bbl",
    "MMBtu_it": "MMBtu",
    "meter ** 3": "cbm",
    "therm_it": "th",
    "metric_tonne": "mt",
    "tonne": "mt",
    "megawatt_hour": "MWh",
    "kilowatt_hour": "kWh",
    "gigawatt_hour": "GWh",
    "dth_it": "Dth",
    "commodity_short_hundredweight": "cwt",
    "Mscf_60F_14_73psia": "Mcf",
    "MMscf_60F_14_73psia": "MMcf",
    "Bscf_60F_14_73psia": "Bcf",
    "Tscf_60F_14_73psia": "Tcf",
    "scf_60F_14_73psia": "scf",
    "boe_58": "boe",
    "basis_point": "bps",
    "percent": "%",
}

def _format_quantity_latex(quantity: pint.Quantity, magnitude: str) -> str:
    """Typeset clean unit names using Pint's structural factors and exponents."""
    escapes = {"\\": r"\backslash{}", "{": r"\{", "}": r"\}", "_": r"\_", "%": r"\%", "$": r"\$", "&": r"\&", "#": r"\#", "^": r"\hat{}", "~": r"\sim{}",}

    def display_latex(display: str) -> str:
        if "/" in display:
            numerator, denominator = display.split("/", 1)
            return rf"\frac{{{display_latex(numerator.strip())}}}{{{display_latex(denominator.strip())}}}"
        escaped = "".join(escapes.get(char, char) for char in display)
        return rf"\mathrm{{{escaped}}}"

    # A complete display override includes composites such as carbon intensity.
    whole_display = _CANONICAL_UNIT_DISPLAY.get(str(quantity.units))
    if whole_display is not None:
        unit_latex = display_latex(whole_display)
    else:
        numerator_terms: list[str] = []
        denominator_terms: list[str] = []
        for name, exponent in quantity.to_tuple()[1]:
            if exponent == 0:
                continue
            display = _CANONICAL_UNIT_DISPLAY.get(name, name)
            term = display_latex(display)
            power = abs(exponent)
            if power != 1:
                if "/" in display:
                    term = rf"\left({term}\right)"
                term = rf"{term}^{{{power}}}"
            (numerator_terms if exponent > 0 else denominator_terms).append(term)
        numerator = r"\,".join(numerator_terms) or "1"
        denominator = r"\,".join(denominator_terms)
        unit_latex = rf"\frac{{{numerator}}}{{{denominator}}}" if denominator else numerator
        if not numerator_terms and not denominator_terms:
            unit_latex = ""

    scientific = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))[eE]([+-]?\d+)", magnitude)
    if scientific:
        magnitude = rf"{scientific[1]}\times 10^{{{int(scientific[2])}}}"
    return rf"{magnitude}\,{unit_latex}" if unit_latex else magnitude


@_locked_registry
def format_quantity(quantity: pint.Quantity, *, style: str = "market", precision: int | None = None,) -> str:
    """Format a quantity as market, standard, canonical, or LaTeX notation.

    Precision is a non-negative number of decimal places; None preserves the
    magnitude's existing representation without currency-based rounding.
    Unknown market tickers fall back to standard unit display without conversion.
    LaTeX returns a math fragment without dollar delimiters.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if style not in {"market", "standard", "canonical", "latex"}:
        raise ValueError(f"Unknown quantity format style: {style!r}")
    if precision is not None and (isinstance(precision, bool) or not isinstance(precision, int) or precision < 0):
        raise ValueError("precision must be a non-negative integer or None.")

    magnitude = str(quantity.magnitude) if precision is None else f"{quantity.magnitude:.{precision}f}"
    unit_string = str(quantity.units)
    if style == "canonical":
        return f"{magnitude} {unit_string}"
    if style == "latex":
        return _format_quantity_latex(quantity, magnitude)
    if style == "market":
        ticker = _MARKET_TICKERS.get(quantity.units)
        if ticker is not None:
            if ticker.startswith(("$", "\u00a3", "\u20ac", "\u00a5")):
                return f"{ticker[0]}{magnitude}{ticker[1:]}"
            return f"{magnitude} {ticker}"

    for canonical, display in _CANONICAL_UNIT_DISPLAY.items():
        unit_string = re.sub(r"(?<![A-Za-z0-9_])" + re.escape(canonical) + r"(?![A-Za-z0-9_])", display, unit_string,)
    return f"{magnitude} {unit_string}"

# ---------------------------------------------------------------------------
# Semantic Classification & Batch Ingestion
# ---------------------------------------------------------------------------

@cache
def _classify_dimensions(dimensions) -> tuple[bool, bool, bool, bool]:
    """Cache classification by immutable Pint dimensionality, including desk units."""
    currencies = {dimension: exponent for dimension, exponent in dimensions.items() if dimension.startswith("[currency_")}
    isolated_currency = len(dimensions) == 1 and list(currencies.values()) == [1]
    physical = {dimension: exponent for dimension, exponent in dimensions.items() if dimension not in currencies}
    commodity_price = (len(currencies) == 1 and list(currencies.values()) == [1] and any(exponent < 0 for exponent in physical.values()))
    flow_rate = not currencies and dimensions.get("[time]") == -1 and any(dimension != "[time]" and exponent > 0 for dimension, exponent in physical.items())
    fx_rate = (not physical and len(currencies) == 2 and sorted(currencies.values()) == [-1, 1])
    return isolated_currency, commodity_price, flow_rate, fx_rate

def is_currency(quantity: pint.Quantity) -> bool:
    """Return True for an isolated currency or its fractional denominations."""
    return isinstance(quantity, pint.Quantity) and _classify_dimensions(quantity.dimensionality)[0]

def is_commodity_price(quantity: pint.Quantity) -> bool:
    """Return True for one currency over physical dimensions, excluding FX rates."""
    return isinstance(quantity, pint.Quantity) and _classify_dimensions(quantity.dimensionality)[1]

def is_flow_rate(quantity: pint.Quantity) -> bool:
    """Return True for non-financial quantities per time, such as bpd or tph."""
    return isinstance(quantity, pint.Quantity) and _classify_dimensions(quantity.dimensionality)[2]

def is_fx_rate(quantity_or_unit: pint.Quantity | pint.Unit) -> bool:
    """Return True for pure currency ratios such as EUR/USD or GBP/USD."""
    return (isinstance(quantity_or_unit, (pint.Quantity, pint.Unit)) and _classify_dimensions(quantity_or_unit.dimensionality)[3])

@_locked_registry
def parse_quantities(values: list[str | int | float | Decimal], units: str | list[str] | None = None, *, default_currency: str | None = None, as_decimal: bool = True,) -> list[pint.Quantity]:
    """Batch parse expressions or magnitudes with shared or paired units.

    Repeated units are sanitised once per call. Paired lists must have equal
    lengths. Separate magnitudes support grouped numbers and accounting signs.
    """
    if units is None:
        return [parse_quantity(value, default_currency=default_currency, as_decimal=as_decimal) for value in values]
    if not isinstance(units, (str, list)):
        raise TypeError("units must be a string, a list of strings, or None.")
    if isinstance(units, list) and len(values) != len(units):
        raise ValueError("values and units must have equal lengths.")

    sanitised_units: dict[str, pint.Unit] = {}
    def resolve_unit(raw_unit: str) -> pint.Unit:
        if raw_unit not in sanitised_units:
            sanitised_units[raw_unit] = ureg.Unit(sanitize_unit_string(raw_unit, default_currency))
        return sanitised_units[raw_unit]

    shared_unit = resolve_unit(units) if isinstance(units, str) else None
    results: list[pint.Quantity] = []
    for index, value in enumerate(values):
        magnitude = _parse_magnitude(value, as_decimal)
        parsed_unit = shared_unit if shared_unit is not None else resolve_unit(units[index])
        results.append(ureg.Quantity(magnitude, parsed_unit))
    return results

# ---------------------------------------------------------------------------
# Quote Structure, Scalar Serialisation & Startup Diagnostics
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class QuoteComponents:
    """Immutable financial and physical metadata for a two-sided market quote."""
    currency: str
    commodity_unit: str
    is_minor_currency: bool
    settlement_decimals: int
    canonical_expression: str

@_locked_registry
def decompose_quote(quote: str, default_currency: str | None = None) -> QuoteComponents:
    """Validate and decompose a currency/commodity quote without price conversion.

    Currency dimensions identify registered proprietary sub-units as well as
    sovereign minor units. Denominators must be physical, not currencies or FX.
    """
    if not isinstance(quote, str):
        raise TypeError("quote must be a string.")
    expression = sanitize_unit_string(quote, default_currency)
    if "/" not in expression:
        raise ValueError(f"Expression {quote!r} is not a two-sided quote.")
    numerator_string, denomination_string = (part.strip() for part in expression.split("/", 1))
    try:
        numerator = ureg.Unit(numerator_string)
        # Derive the denominator from the complete expression so chained
        # divisions preserve Pint's associativity rather than reversing a rate.
        complete_unit = ureg.Unit(expression)
        denomination = numerator / complete_unit
        dimensions = numerator.dimensionality
        if len(dimensions) != 1:
            raise ValueError("Quote numerator must be an isolated currency.")
        dimension, exponent = next(iter(dimensions.items()))
        if exponent != 1 or not dimension.startswith("[currency_"):
            raise ValueError("Quote numerator must be an isolated currency.")
        currency = dimension[len("[currency_"):-1]
        if currency.upper() not in _CURRENCY_CANONICAL_CODES:
            raise ValueError("Quote numerator currency is not registered.")
        currency = _CURRENCY_CANONICAL_CODES[currency.upper()]
        if not denomination.dimensionality or any(name.startswith("[currency_") for name in denomination.dimensionality):
            raise ValueError("Quote denominator must represent a physical commodity unit.")
        minor = ureg.Quantity(1, numerator).to(currency).magnitude < 1
    except (pint.PintError, ValueError, TypeError) as error:
        raise ValueError(f"Invalid commodity quote {quote!r}: {error}") from error
    return QuoteComponents(currency=currency, commodity_unit=str(denomination), is_minor_currency=minor, settlement_decimals=CURRENCY_MINOR_DIGITS[currency], canonical_expression=f"{numerator} / ({denomination})" if "/" in str(denomination) else f"{numerator} / {denomination}",)

@_locked_registry
def quantity_to_dict(quantity: pint.Quantity) -> dict[str, str]:
    """Serialise a finite scalar quantity to JSON-safe magnitude/unit strings."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    magnitude = quantity.magnitude
    if isinstance(magnitude, bool) or not isinstance(magnitude, (int, float, Decimal)):
        raise TypeError("Only scalar int, float, and Decimal magnitudes are supported.")
    if not Decimal(str(magnitude)).is_finite():
        raise ValueError("Quantity magnitude must be finite.")
    return {"m": str(magnitude), "u": str(quantity.units)}

@_locked_registry
def quantity_from_dict(payload: dict[str, Any], *, as_decimal: bool = True) -> pint.Quantity:
    """Restore exact stored unit identifiers, bypassing market alias normalisation."""
    if not isinstance(payload, dict) or not {"m", "u"} <= payload.keys():
        raise ValueError("Invalid quantity payload; must contain 'm' and 'u'.")
    raw_magnitude, raw_unit = payload["m"], payload["u"]
    if isinstance(raw_magnitude, bool) or not isinstance(raw_magnitude, (str, int, float, Decimal)):
        raise ValueError("Payload magnitude must be a numeric scalar or string.")
    if not isinstance(raw_unit, str) or not raw_unit.strip():
        raise ValueError("Payload unit must be a non-empty string.")
    try:
        decimal_magnitude = Decimal(str(raw_magnitude))
        if not decimal_magnitude.is_finite():
            raise ValueError("Payload magnitude must be finite.")
        magnitude = decimal_magnitude if as_decimal else float(decimal_magnitude)
        if not Decimal(str(magnitude)).is_finite():
            raise ValueError("Payload magnitude overflows float precision.")
        return ureg.Quantity(magnitude, ureg.Unit(raw_unit))
    except (InvalidOperation, pint.PintError, TypeError, OverflowError) as error:
        raise ValueError(f"Invalid quantity payload: {error}") from error

@authorized_conversion
def verify_registry_integrity() -> None:
    """Audit static definitions, dimensional isolation, shields and cached quotes.

    Raises AssertionError with the failing invariant, including under python -O.
    This diagnostic validates existing baselines; it registers no conversions.
    """
    def check_scale(name: str, quantity: pint.Quantity, target: str, expected: float) -> None:
        actual = float(quantity.to(target).magnitude)
        if not abs(actual - expected) <= max(1e-9, abs(expected) * 1e-12):
            raise AssertionError(f"{name}: expected {expected} {target}, got {actual}.")

    check_scale("Petroleum barrel", ureg.Quantity(1, "bbl"), "gallon", 42)
    check_scale("IT Btu", ureg.Quantity(1, "btu_it"), "joule", 1055.05585262)
    try:
        ureg.Quantity(1, "USD").to("EUR")
    except pint.DimensionalityError:
        pass
    else:
        raise AssertionError("USD converted to EUR without an explicit FX context.")
    for token, expected in (("Mcf", 1e3), ("MMcf", 1e6), ("Bcf", 1e9), ("Tcf", 1e12)):
        check_scale(token, ureg.Quantity(1, sanitize_unit_string(token)), "scf_60F_14_73psia", expected)
    check_scale("Basis point", ureg.Quantity(100, "basis_point"), "percent", 1)
    weather_index = ureg.Quantity(10, "delta_degF") * ureg.Quantity(1, "day")
    check_scale("Weather degree day", weather_index, "degree_day", 10)
    check_scale("Weather degree hour", ureg.Quantity(1, "degree_day_C"), "degree_hour_C", 24)
    heat_rate = ureg.Quantity(7.5, "MMBtu_it / megawatt_hour")
    expected_dimensions = (ureg.joule / (ureg.watt * ureg.second)).dimensionality
    if heat_rate.dimensionality != expected_dimensions:
        raise AssertionError("Heat rate must retain energy-per-electrical-output dimensions.")
    check_scale("Heat rate quoting scale", ureg.Quantity(1, "heat_rate_macro"), "heat_rate_us", 1000)
    check_scale("Metric heat rate quoting scale", ureg.Quantity(1, "heat_rate_gas_turbine"), "heat_rate_metric", 1)
    check_scale("Hydrogen mass", ureg.Quantity(1, "tonne_H2"), "kg_H2", 1000)
    check_scale("Hydrogen normal volume", ureg.Quantity(1, "Nm3_H2"), "Nm3_0C_101325kPa_dry", 1)
    check_scale("SAF gallon", ureg.Quantity(1, "gallon_SAF"), "gallon", 1)
    ci_dimensions = (ureg.gram / ureg.megajoule).dimensionality
    if ureg.Quantity(35.5, "CI").dimensionality != ci_dimensions:
        raise AssertionError("Carbon intensity must have mass-per-energy dimensions.")
    check_scale("Carbon intensity mass basis", ureg.Quantity(1, "gCO2e"), "gram", 1)
    check_scale("Absolute pipeline pressure", ureg.Quantity(14.73, "psia"), "psi", 14.73)
    check_scale("Gauge pipeline pressure scale", ureg.Quantity(600, "psig"), "psi", 600)
    check_scale("Absolute bar pressure", ureg.Quantity(1, "bara"), "bar", 1)
    check_scale("Gauge bar pressure scale", ureg.Quantity(1, "barg"), "bar", 1)
    check_scale("Gold kilobar", ureg.Quantity(1, "kilobar_gold"), "toz", 1000 / 31.1034768)
    check_scale("Ten tola bar", ureg.Quantity(1, "ten_tola_bar"), "tola", 10)
    check_scale("Nominal gold delivery bar", ureg.Quantity(1, "good_delivery_gold"), "toz", 400)
    check_scale("Nominal silver delivery bar", ureg.Quantity(1, "good_delivery_silver"), "toz", 1000)
    check_scale("Power ramp scale", ureg.Quantity(1, "MW_per_minute"), "kW_per_minute", 1000)
    check_scale("Battery C-rate", ureg.Quantity(1, "C_rate"), "1 / hour", 1)
    check_scale("Air-dry pulp label", ureg.Quantity(1, "admt"), "metric_ton", 1)
    check_scale("Bone-dry pulp label", ureg.Quantity(1, "bdmt"), "metric_ton", 1)
    check_scale("Corn crop bushel", ureg.Quantity(1, "bushel_corn"), "pound", 56)
    check_scale("Canola crop bushel", ureg.Quantity(1, "bushel_canola"), "pound", 50)
    check_scale("Rice crop bushel", ureg.Quantity(1, "bushel_rice"), "pound", 45)
    check_scale("Metric absolute pressure", ureg.Quantity(101.325, "kpaa"), "kilopascal", 101.325)
    check_scale("Metric gauge pressure", ureg.Quantity(500, "kpag"), "kilopascal", 500)
    check_scale("Uranium tonne", ureg.Quantity(1, "mtU"), "kgU", 1000)
    check_scale("Metric tonne unit isolation", ureg.Quantity(1, "mtu"), "kilogram", 10)
    for unit, dimension in (
        ("m3_20C_101325kPa_dry", "[standard_gas_volume_20C_101325kPa_dry]"),
        ("scf_60F_14696psi_dry", "[standard_gas_volume_60F_1atm_dry]"),
    ):
        if ureg.Unit(unit).dimensionality != {dimension: 1}:
            raise AssertionError(f"Gas reference dimensional isolation failed for {unit}.")
    with _REGISTRY_LOCK:
        quotes = tuple(_QUOTE_REPLACEMENTS.items())
    for screen_quote, expression in quotes:
        try:
            ureg.Unit(expression)
        except (pint.PintError, ValueError, TypeError) as error:
            raise AssertionError(f"Invalid cached quote {screen_quote!r}: {expression!r}") from error

# ---------------------------------------------------------------------------
# Module Type Aliases & Public Exports
# ---------------------------------------------------------------------------

Unit: Final = ureg.Unit
Quantity: Final = ureg.Quantity

__all__: Final[tuple[str, ...]] = (
    # =========================================================================
    # Core Registry Singletons, Constructors & Types
    # =========================================================================
    "ureg",
    # The master, thread-safe pint.UnitRegistry singleton instance. It holds the pre-compiled dimensional
    # definitions for every sovereign currency, physical commodity, thermodynamic reference condition,
    # and environmental compliance instrument.
    #
    # from registry import ureg
    # # Direct unit multiplication and registry inspection
    # barrel_vol = 100 * ureg.oil_barrel
    # energy_density = ureg("MMBtu_it / megawatt_hour")
    # # Accessing base dimensions directly
    # print(ureg.USD)           # 1 <Unit('USD')>
    # print(ureg.oil_barrel)    # 1 <Unit('oil_barrel')>

    "Q_",
    # Direct alias for ureg.Quantity. Used as a concise shorthand to instantiate new quantities without
    # needing verbose method calls.
    #
    # from decimal import Decimal
    # from registry import Q_
    # # Construct quantities using standard magnitude and unit parameters
    # crude_lot = Q_(1000, "bbl")
    # power_price = Q_(Decimal("84.50"), "EUR / megawatt_hour")

    "Unit",
    # Type alias pointing to ureg.Unit. Used in function signatures, data models (such as Pydantic models or
    # dataclasses), and type annotations to enforce that an argument is a unit rather than a magnitude or quantity.
    #
    # from registry import Unit, parse_unit
    #
    # def validate_trading_book_unit(expected_unit: Unit, trade_unit: Unit) -> bool:
    #     return expected_unit.dimensionality == trade_unit.dimensionality
    #
    # book_unit: Unit = parse_unit("USD / oil_barrel")

    "Quantity",
    # Type alias pointing to ureg.Quantity. Used for strict type hints in execution algorithms, pricing
    # calculators, and risk models to guarantee inputs carry physical dimensions.
    #
    # from registry import Quantity, parse_quantity
    #
    # def calculate_gross_notional(volume: Quantity, price: Quantity) -> Quantity:
    #     notional = volume * price
    #     return notional
    #
    # vol: Quantity = parse_quantity("50000 bbl")
    # px: Quantity = parse_quantity("$75.20/bbl")
    # notional: Quantity = calculate_gross_notional(vol, px)

    "get_registry",
    # The cached factory function (@cache) that initialises, configures, and binds the master UnitRegistry
    # singleton. It handles Pint's application registry binding (pint.set_application_registry) so
    # multiprocessing and pickle operations do not lose custom dimensions.
    #
    # from registry import get_registry
    #
    # # Retrieve the active singleton (returns the identical cached instance)
    # registry = get_registry()
    # assert registry is ureg

    "verify_registry_integrity",
    # A self-audit diagnostics suite designed to run during container boot sequences and test runners.
    # Asserts critical physical invariants—confirming that a petroleum barrel resolves to 42 US gallons,
    # that IT Btus match the 1956 International Steam Table joule value, that sovereign currencies remain
    # dimensionally isolated without an FX context, and that Roman numeral gas notations scale correctly.
    #
    # from registry import verify_registry_integrity
    #
    # def test_startup_invariants():
    #     verify_registry_integrity()
    #     print("Registry verified: All dimensional invariants intact.")

    "sandbox_registry",
    # A context manager providing an isolated, ephemeral fork of the registry loaded with all commodity
    # and currency definitions, leaving the master singleton untouched and un-frozen. Ideal for test
    # fixtures, backtesting runs, or evaluating untrusted counterparty definitions.
    #
    # from registry import sandbox_registry
    #
    # with sandbox_registry() as test_reg:
    #     test_reg.define("pipeline_batch = 10000 * oil_barrel")
    #     batch = test_reg.Quantity(2, "pipeline_batch")
    #     assert batch.to("oil_barrel").magnitude == 20000

    # =========================================================================
    # Immutability & Lifecycle Governance
    # =========================================================================
    "freeze_registry",
    # A production security switch that permanently locks the registry against further runtime additions.
    # Calling this prevents rogue background tasks or malformed order tickets from redefining existing
    # units or injecting unintended currency dimensions into production.
    #
    # from registry import freeze_registry, register_commodity_unit
    #
    # freeze_registry()
    # try:
    #     register_commodity_unit("desk_tank = 500 * barrel")
    # except RuntimeError as err:
    #     print(f"Action blocked: {err}")

    "is_registry_frozen",
    # A thread-safe boolean query returning True if the registry has been sealed via freeze_registry(),
    # or False if it remains open to dynamic configuration.
    #
    # from registry import is_registry_frozen
    #
    # if not is_registry_frozen():
    #     # Safe to load desk-specific bespoke units
    #     pass

    # =========================================================================
    # Lexical Ingestion, Parsing & Deconstruction
    # =========================================================================
    "sanitize_unit_string",
    # The lexical pre-parser. Cleans broker tickers and vendor shorthands into Pint-compliant algebraic
    # unit strings. Strips leading currency symbols, resolves Roman numeral conflicts (MMcf, MMbbl), converts
    # colloquial English 'per' into '/', normalises exchange minor units (p/, c/), and maps cached market tickers.
    #
    # from registry import sanitize_unit_string
    #
    # print(sanitize_unit_string("$/bbl"))               # "USD / oil_barrel"
    # print(sanitize_unit_string("p/th"))                # "GBP_penny / therm_it"
    # print(sanitize_unit_string("MMBtu per hour"))       # "MMBtu_it / hour"
    # print(sanitize_unit_string("100 Mcfd"))            # "100 Mscf_60F_14_73psia / day"
    # print(sanitize_unit_string("$/bbl", default_currency="CAD")) # "CAD / oil_barrel"

    "parse_quantity",
    # The primary ingestion engine. Converts numbers, dirty screen strings, accounting formats, and
    # array-like magnitudes into strongly typed pint.Quantity objects. Handles grouped thousands (1,500),
    # accounting negative parentheses ((45.20)), currency prefixes ($75.50/bbl), percentages (2.5%),
    # and vectorised arrays.
    #
    # from decimal import Decimal
    # import numpy as np
    # from registry import parse_quantity
    #
    # q1 = parse_quantity("$75.50/bbl")                  # 75.50 <Unit('USD / oil_barrel')>
    # q2 = parse_quantity("($1,250.00/mt)")              # -1250.00 <Unit('USD / tonne')>
    # q3 = parse_quantity("(500) MWh")                   # -500 <Unit('megawatt_hour')>
    # q4 = parse_quantity("2.5%")                        # 2.5 <Unit('percent')>
    # q5 = parse_quantity(Decimal("100.25"), "USD/bbl")  # Decimal magnitude preserved
    # prices = parse_quantity(np.array([72.5, 73.1, 71.9]), "$/bbl")

    "parse_quantities",
    # A batch ingestion utility optimised for large datasets (e.g. forward curve tenors or tick-data feeds).
    # Parses arrays of magnitudes sharing the same unit, caching the sanitised unit once rather than
    # re-parsing it for every row.
    #
    # from registry import parse_quantities
    #
    # curve_ticks = ["75.10", "75.25", "75.40", "75.80"]
    # quantities = parse_quantities(curve_ticks, units="$/bbl")
    #
    # legs_val = ["1000", "5000", "200"]
    # legs_units = ["bbl", "MMBtu", "MWh"]
    # trade_legs = parse_quantities(legs_val, units=legs_units)

    "parse_unit",
    # Parses and normalises a raw unit string directly into a validated pint.Unit instance without
    # requiring a dummy numeric magnitude.
    #
    # from registry import parse_unit
    #
    # fuel_unit = parse_unit("c/gal")
    # gas_unit = parse_unit("$/MMBTU")
    # power_unit = parse_unit("£/MWh")
    # print(type(fuel_unit))  # <class 'pint.Unit'>
    # print(fuel_unit)        # USD_cent / gallon

    "decompose_quote",
    # Breaks down any <currency>/<commodity> quote string into its discrete financial, physical,
    # and settlement parameters.
    #
    # from registry import decompose_quote
    #
    # quote = decompose_quote("c/bu")
    # print(quote.currency)              # "USD"
    # print(quote.commodity_unit)        # "bushel"
    # print(quote.is_minor_currency)     # True
    # print(quote.settlement_decimals)   # 2
    # print(quote.canonical_expression)  # "USD_cent / bushel"

    "QuoteComponents",
    # The immutable dataclass(frozen=True, slots=True) returned by decompose_quote. Provides typed field
    # access to the decomposed elements of a quotation.
    #
    # from registry import QuoteComponents, decompose_quote
    #
    # comp: QuoteComponents = decompose_quote("p/th")
    # if comp.is_minor_currency:
    #     print(f"Settling in sub-unit of {comp.currency} with {comp.settlement_decimals} decimals.")

    # =========================================================================
    # Formatting & Lossless Serialisation
    # =========================================================================
    "format_quantity",
    # The reverse canonicaliser. Converts internal Pint quantities back into clean, institutional display
    # formats across four selectable styles: 'market' (screen tickers like $/bbl), 'standard', 'canonical',
    # and 'latex'.
    #
    # from registry import parse_quantity, format_quantity
    #
    # q = parse_quantity("75.5 USD / oil_barrel")
    # print(format_quantity(q, style="market"))             # "$75.50/bbl"
    # print(format_quantity(q, style="standard"))           # "75.5 USD/bbl"
    # print(format_quantity(q, style="canonical"))          # "75.5 USD / oil_barrel"
    # print(format_quantity(q, style="latex", precision=2)) # r"75.50\,\mathrm{\frac{USD}{bbl}}"

    "quantity_to_dict",
    # Serialises a pint.Quantity into a JSON-serialisable dictionary while retaining exact magnitude precision
    # (converting Decimal to string). Prevents binary float drift across Redis caches, message queues, and document databases.
    #
    # from registry import parse_quantity, quantity_to_dict
    #
    # trade_vol = parse_quantity("1500.50 bbl")
    # packet = quantity_to_dict(trade_vol)
    # print(packet)  # {'m': '1500.50', 'u': 'oil_barrel'}

    "quantity_from_dict",
    # The exact inverse of quantity_to_dict. Reconstructs a pint.Quantity from a dictionary packet,
    # restoring full Decimal precision by default.
    #
    # from registry import quantity_from_dict
    #
    # packet = {'m': '1500.50', 'u': 'oil_barrel'}
    # restored_q = quantity_from_dict(packet, as_decimal=True)
    # print(restored_q)  # 1500.50 oil_barrel

    # =========================================================================
    # Dynamic Registrations
    # =========================================================================
    "register_currency",
    # Thread-safely adds a new sovereign currency, cryptocurrency, or central bank digital currency (CBDC)
    # into the registry before it is frozen. Automatically registers its independent dimension, typographical
    # symbols, fractional sub-units, settlement rounding exponents, and quote aliases.
    #
    # from registry import register_currency, parse_quantity
    #
    # register_currency(
    #     "BTC",
    #     decimals=8,
    #     symbols=("₿",),
    #     aliases=("XBT", "bitcoin"),
    #     sub_units={"satoshi": "1e-8 * BTC = sat = sats"},
    #     quote_aliases={"BTC/USD": "BTC / USD"}
    # )
    # q_btc = parse_quantity("₿0.5")
    # print(q_btc.to("satoshi"))  # 50000000.0 satoshi

    "register_commodity_unit",
    # Thread-safely registers a desk-specific, bespoke, or non-standard physical container or parcel into
    # the active registry.
    #
    # from registry import register_commodity_unit, parse_quantity
    #
    # register_commodity_unit(
    #     "drum_200l = 200 * liter = standard_drum = steel_drum",
    #     quote_aliases={"$/drum": "USD / drum_200l"}
    # )
    # drum_order = parse_quantity("50 drum_200l")
    # print(drum_order.to("liter"))  # 10000 liter

    # =========================================================================
    # Dimensional Governance & Semantic Predicates
    # =========================================================================
    "validate_dimensionality",
    # Verifies whether a quantity strictly matches an expected physical dimension or reference unit string.
    # Prevents operational mismatches (such as booking volumetric barrels against a mass-denominated forward curve).
    #
    # from registry import parse_quantity, validate_dimensionality
    #
    # cargo = parse_quantity("600000 bbl")
    # assert validate_dimensionality(cargo, "gallon") == True
    # assert validate_dimensionality(cargo, "meter ** 3") == True
    # assert validate_dimensionality(cargo, "tonne") == False  # Barrels are volume, not mass!

    "ensure_unit",
    # Converts a quantity to a specified target unit, raising a clean pint.DimensionalityError if the
    # conversion violates dimensional physics.
    #
    # from registry import parse_quantity, ensure_unit
    #
    # gas_flow = parse_quantity("50000 MMBtu")
    # dth_flow = ensure_unit(gas_flow, "Dth")  # Converts cleanly within same dimension
    # try:
    #     ensure_unit(gas_flow, "mt")          # Raises DimensionalityError
    # except Exception as err:
    #     print(f"Blocked: {err}")

    "are_convertible",
    # A fast boolean predicate checking whether two units or quantities share the exact same physical
    # dimensionality without raising an unhandled exception.
    #
    # from registry import are_convertible
    #
    # print(are_convertible("bbl", "gal"))       # True (both volume)
    # print(are_convertible("bbl", "mt"))        # False (volume vs mass)
    # print(are_convertible("MWh", "MMBtu"))     # True (both energy)
    # print(are_convertible("$/bbl", "EUR/bbl")) # False (USD and EUR are isolated dimensions)

    "is_currency",
    # Returns True if a quantity or unit represents a purely financial sovereign currency or fractional
    # denomination (with no physical denominator).
    #
    # from registry import parse_quantity, is_currency
    #
    # print(is_currency(parse_quantity("100 USD")))        # True
    # print(is_currency(parse_quantity("500 GBP_penny")))  # True
    # print(is_currency(parse_quantity("$75.50/bbl")))     # False (commodity price)

    "is_fx_rate",
    # Returns True if the expression represents a pure foreign exchange currency cross rate ([currency_A] / [currency_B]).
    #
    # from registry import parse_quantity, parse_unit, is_fx_rate
    #
    # eur_usd = parse_quantity("1.0850 EUR / USD")
    # gbp_usd_unit = parse_unit("GBP / USD")
    # print(is_fx_rate(eur_usd))                           # True
    # print(is_fx_rate(gbp_usd_unit))                      # True
    # print(is_fx_rate(parse_quantity("$75.50/bbl")))      # False

    "is_commodity_price",
    # Returns True if the quantity represents a valid commodity screen quote ([currency] / [physical_dimension]).
    #
    # from registry import parse_quantity, is_commodity_price
    #
    # print(is_commodity_price(parse_quantity("$75.50/bbl")))   # True
    # print(is_commodity_price(parse_quantity("p/th")))         # True
    # print(is_commodity_price(parse_quantity("100 bbl")))      # False (volume)
    # print(is_commodity_price(parse_quantity("1.08 EUR/USD"))) # False (FX rate)

    "is_flow_rate",
    # Returns True if the quantity represents a physical commodity delivered over time ([physical] / [time]).
    #
    # from registry import parse_quantity, is_flow_rate
    #
    # print(is_flow_rate(parse_quantity("10000 bpd")))  # True
    # print(is_flow_rate(parse_quantity("500 tph")))    # True (tonnes per hour)
    # print(is_flow_rate(parse_quantity("50 MMscfd")))  # True
    # print(is_flow_rate(parse_quantity("1000 bbl")))   # False (static volume)

    # =========================================================================
    # Domain Reference Metadata Collections
    # =========================================================================
    "CURRENCY_CODES",
    # An immutable tuple[str, ...] containing every official ISO-4217 sovereign currency code recognised
    # globally (USD, EUR, GBP, JPY, AED, BHD, etc.).
    #
    # from registry import CURRENCY_CODES
    #
    # if "EUR" in CURRENCY_CODES:
    #     print("EUR is a recognised sovereign currency.")

    "CURRENCY_SYMBOL_TO_CODE",
    # A dictionary mapping raw currency glyphs, multi-character symbols, and wire abbreviations
    # ('$', '£', '€', 'د.إ', 'US$', 'A$') to their canonical ISO-4217 3-letter codes.
    #
    # from registry import CURRENCY_SYMBOL_TO_CODE
    #
    # print(CURRENCY_SYMBOL_TO_CODE["$"])    # "USD"
    # print(CURRENCY_SYMBOL_TO_CODE["£"])    # "GBP"
    # print(CURRENCY_SYMBOL_TO_CODE["A$"])   # "AUD"
    # print(CURRENCY_SYMBOL_TO_CODE["د.إ"])  # "AED"

    "CURRENCY_MINOR_DIGITS",
    # A lookup dictionary mapping each ISO-4217 currency code to its statutory minor unit settlement precision
    # (number of decimal places). Used downstream by settlement accounting and cash-rounding engines to prevent decimal drift.
    #
    # from registry import CURRENCY_MINOR_DIGITS
    #
    # print(CURRENCY_MINOR_DIGITS["USD"])  # 2 (cents)
    # print(CURRENCY_MINOR_DIGITS["BHD"])  # 3 (fils - 3 decimal places)
    # print(CURRENCY_MINOR_DIGITS["JPY"])  # 0 (yen - zero decimal places in active trade)
    # print(CURRENCY_MINOR_DIGITS["CLP"])  # 0 (Chilean Peso - zero decimal places)

    "COMMODITY_UNIT_DEFINITIONS",
    # The master, immutable tuple[str, ...] holding all raw definition strings used to build the unit registry.
    # Contains the exact mathematical specifications for energy, thermodynamic gas volumes, petroleum barrels,
    # mass standards, agricultural bushels, carbon compliance units, and grid ramp rates.
    #
    # from registry import COMMODITY_UNIT_DEFINITIONS
    #
    # for definition in COMMODITY_UNIT_DEFINITIONS[:5]:
    #     print(definition)
)
