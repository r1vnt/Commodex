"""commodities.core.units.registry.

Institutional-grade unit registry for physical and financial commodities.
Configures independent currency dimensions, thermodynamic gas volume standards,
petroleum volumes, flow rates, environmental instruments, and strict serialisation.
"""

from __future__ import annotations

from decimal import Decimal
import re
import threading
from typing import Any, Final
import pint

# ---------------------------------------------------------------------------
# ISO-4217 Currency Codes
# ---------------------------------------------------------------------------

CURRENCY_CODES: Final[tuple[str, ...]] = tuple(
    """
AED AFN ALL AMD AOA ARS AUD AWG AZN
BAM BBD BDT BHD BIF BMD BND BOB BRL BSD BTN BWP BYN BZD
CAD CDF CHF CLP CNY COP CRC CUP CVE CZK
DJF DKK DOP DZD
EGP ERN ETB EUR
FJD FKP
GBP GEL GHS GIP GMD GNF GTQ GYD
HKD HNL HTG HUF
IDR ILS INR IQD IRR ISK
JMD JOD JPY
KES KGS KHR KMF KPW KRW KWD KYD KZT
LAK LBP LKR LRD LSL LYD
MAD MDL MGA MKD MMK MNT MOP MRU MUR MVR MWK MXN MYR MZN
NAD NGN NIO NOK NPR NZD
OMR
PAB PEN PGK PHP PKR PLN PYG
QAR
RON RSD RUB RWF
SAR SBD SCR SDG SEK SGD SHP SLE SOS SRD SSP STN SVC SYP SZL
THB TJS TMT TND TOP TRY TTD TWD TZS
UAH UGX USD UYU UZS
VES VND VUV WST
XAF XCD XCG XOF XPF
YER ZAR ZMW ZWG
""".split()
)

# Display aliases and vendor symbols mapped to ISO codes.
CURRENCY_SYMBOL_TO_CODE: Final[dict[str, str]] = {
    # --- AED (United Arab Emirates Dirham) ---
    "د.إ": "AED",
    # --- AFN (Afghan Afghani) ---
    "؋": "AFN",
    # --- ALL (Albanian Lek) ---
    "L": "ALL",
    # --- AMD (Armenian Dram) ---
    "֏": "AMD",
    # --- AOA (Angolan Kwanza) ---
    "Kz": "AOA",
    # --- ARS (Argentine Peso) ---
    "$_ARS": "ARS",
    # --- AUD (Australian Dollar) ---
    "A$": "AUD",
    # --- AWG (Aruban Florin) ---
    "ƒ": "AWG",
    # --- AZN (Azerbaijani Manat) ---
    "₼": "AZN",
    # --- BAM (Bosnia and Herzegovina Convertible Mark) ---
    "KM": "BAM",
    # --- BBD (Barbadian Dollar) ---
    "Bds$": "BBD",
    # --- BDT (Bangladeshi Taka) ---
    "৳": "BDT",
    # --- BHD (Bahraini Dinar) ---
    "د.ب": "BHD",
    # --- BIF (Burundian Franc) ---
    "FBu": "BIF",
    # --- BMD (Bermudian Dollar) ---
    "BD$": "BMD",
    # --- BND (Brunei Dollar) ---
    "B$": "BND",
    # --- BOB (Bolivian Boliviano) ---
    "Bs.": "BOB",
    # --- BRL (Brazilian Real) ---
    "R$": "BRL",
    # --- BSD (Bahamian Dollar) ---
    "BS$": "BSD",
    # --- BTN (Bhutanese Ngultrum) ---
    "Nu.": "BTN",
    # --- BWP (Botswana Pula) ---
    "P": "BWP",
    # --- BYN (Belarusian Ruble) ---
    "Br": "BYN",
    # --- BZD (Belize Dollar) ---
    "BZ$": "BZD",
    # --- CAD (Canadian Dollar) ---
    "C$": "CAD",
    # --- CDF (Congolese Franc) ---
    "FC": "CDF",
    # --- CHF (Swiss Franc) ---
    "Fr": "CHF",
    # --- CLP (Chilean Peso) ---
    "$_CLP": "CLP",
    # --- CNY (Chinese Yuan Renminbi) ---
    "¥_CNY": "CNY",
    # --- COP (Colombian Peso) ---
    "$_COP": "COP",
    # --- CRC (Costa Rican Colón) ---
    "₡": "CRC",
    # --- CUP (Cuban Peso) ---
    "$_CUP": "CUP",
    # --- CVE (Cape Verdean Escudo) ---
    "Esc": "CVE",
    # --- CZK (Czech Koruna) ---
    "Kč": "CZK",
    # --- DJF (Djiboutian Franc) ---
    "Fdj": "DJF",
    # --- DKK (Danish Krone) ---
    "kr_DKK": "DKK",
    # --- DOP (Dominican Peso) ---
    "RD$": "DOP",
    # --- DZD (Algerian Dinar) ---
    "دج": "DZD",
    # --- EGP (Egyptian Pound) ---
    "E£": "EGP",
    # --- ERN (Eritrean Nakfa) ---
    "Nfk": "ERN",
    # --- ETB (Ethiopian Birr) ---
    "Br_ETB": "ETB",
    # --- EUR (Euro) ---
    "€": "EUR",
    # --- FJD (Fijian Dollar) ---
    "FJ$": "FJD",
    # --- FKP (Falkland Islands Pound) ---
    "£_FKP": "FKP",
    # --- GBP (British Pound Sterling) ---
    "£": "GBP",
    # --- GEL (Georgian Lari) ---
    "₾": "GEL",
    # --- GHS (Ghanaian Cedi) ---
    "₵": "GHS",
    # --- GIP (Gibraltar Pound) ---
    "£_GIP": "GIP",
    # --- GMD (Gambian Dalasi) ---
    "D": "GMD",
    # --- GNF (Guinean Franc) ---
    "FG": "GNF",
    # --- GTQ (Guatemalan Quetzal) ---
    "Q": "GTQ",
    # --- GYD (Guyanese Dollar) ---
    "GY$": "GYD",
    # --- HKD (Hong Kong Dollar) ---
    "HK$": "HKD",
    # --- HNL (Honduran Lempira) ---
    "L_HNL": "HNL",
    # --- HTG (Haitian Gourde) ---
    "G": "HTG",
    # --- HUF (Hungarian Forint) ---
    "Ft": "HUF",
    # --- IDR (Indonesian Rupiah) ---
    "Rp": "IDR",
    # --- ILS (Israeli New Shekel) ---
    "₪": "ILS",
    # --- INR (Indian Rupee) ---
    "₹": "INR",
    # --- IQD (Iraqi Dinar) ---
    "ع.د": "IQD",
    # --- IRR (Iranian Rial) ---
    "﷼_IRR": "IRR",
    # --- ISK (Icelandic Króna) ---
    "kr_ISK": "ISK",
    # --- JMD (Jamaican Dollar) ---
    "J$": "JMD",
    # --- JOD (Jordanian Dinar) ---
    "د.ا": "JOD",
    # --- JPY (Japanese Yen) ---
    "¥_JPY": "JPY",
    # --- KES (Kenyan Shilling) ---
    "KSh": "KES",
    # --- KGS (Kyrgyzstani Som) ---
    "сом": "KGS",
    # --- KHR (Cambodian Riel) ---
    "៛": "KHR",
    # --- KMF (Comorian Franc) ---
    "CF": "KMF",
    # --- KPW (North Korean Won) ---
    "₩_KPW": "KPW",
    # --- KRW (South Korean Won) ---
    "₩_KRW": "KRW",
    # --- KWD (Kuwaiti Dinar) ---
    "د.ك": "KWD",
    # --- CI$ (Cayman Islands Dollar) ---
    "CI$": "KYD",
    # --- KZT (Kazakhstani Tenge) ---
    "₸": "KZT",
    # --- LAK (Lao Kip) ---
    "₭": "LAK",
    # --- LBP (Lebanese Pound) ---
    "ل.ل": "LBP",
    # --- LKR (Sri Lankan Rupee) ---
    "Rs_LKR": "LKR",
    # --- LRD (Liberian Dollar) ---
    "L$": "LRD",
    # --- LSL (Lesotho Loti) ---
    "L_LSL": "LSL",
    # --- LYD (Libyan Dinar) ---
    "ل.د": "LYD",
    # --- MAD (Moroccan Dirham) ---
    "د.م.": "MAD",
    # --- MDL (Moldovan Leu) ---
    "L_MDL": "MDL",
    # --- MGA (Malagasy Ariary) ---
    "Ar": "MGA",
    # --- MKD (Macedonian Denar) ---
    "ден": "MKD",
    # --- MMK (Myanmar Kyat) ---
    "K": "MMK",
    # --- MNT (Mongolian Tögrög) ---
    "₮": "MNT",
    # --- MOP (Macanese Pataca) ---
    "MOP$": "MOP",
    # --- MRU (Mauritanian Ouguiya) ---
    "UM": "MRU",
    # --- MUR (Mauritian Rupee) ---
    "₨_MUR": "MUR",
    # --- MVR (Maldivian Rufiyaa) ---
    "Rf": "MVR",
    # --- MWK (Malawian Kwacha) ---
    "MK": "MWK",
    # --- MXN (Mexican Peso) ---
    "Mex$": "MXN",
    # --- MYR (Malaysian Ringgit) ---
    "RM": "MYR",
    # --- MZN (Mozambican Metical) ---
    "MT": "MZN",
    # --- NAD (Namibian Dollar) ---
    "N$": "NAD",
    # --- NGN (Nigerian Naira) ---
    "₦": "NGN",
    # --- NIO (Nicaraguan Córdoba) ---
    "C$_NIO": "NIO",
    # --- NOK (Norwegian Krone) ---
    "kr_NOK": "NOK",
    # --- NPR (Nepalese Rupee) ---
    "रू": "NPR",
    # --- NZD (New Zealand Dollar) ---
    "NZ$": "NZD",
    # --- OMR (Omani Rial) ---
    "ر.ع.": "OMR",
    # --- PAB (Panamanian Balboa) ---
    "B/.": "PAB",
    # --- PEN (Peruvian Sol) ---
    "S/": "PEN",
    # --- PGK (Papua New Guinean Kina) ---
    "K_PGK": "PGK",
    # --- PHP (Philippine Peso) ---
    "₱": "PHP",
    # --- PKR (Pakistani Rupee) ---
    "Rs_PKR": "PKR",
    # --- PLN (Polish Złoty) ---
    "zł": "PLN",
    # --- PYG (Paraguayan Guaraní) ---
    "₲": "PYG",
    # --- QAR (Qatari Riyal) ---
    "ر.ق": "QAR",
    # --- RON (Romanian Leu) ---
    "lei": "RON",
    # --- RSD (Serbian Dinar) ---
    "дин.": "RSD",
    # --- RUB (Russian Ruble) ---
    "₽": "RUB",
    # --- RWF (Rwandan Franc) ---
    "FRw": "RWF",
    # --- SAR (Saudi Riyal) ---
    "ر.س": "SAR",
    # --- SBD (Solomon Islands Dollar) ---
    "SI$": "SBD",
    # --- SCR (Seychellois Rupee) ---
    "₨_SCR": "SCR",
    # --- SDG (Sudanese Pound) ---
    "ج.س.": "SDG",
    # --- SEK (Swedish Krona) ---
    "kr_SEK": "SEK",
    # --- SGD (Singapore Dollar) ---
    "S$": "SGD",
    # --- SHP (Saint Helena Pound) ---
    "£_SHP": "SHP",
    # --- SLE (Sierra Leonean Leone) ---
    "Le": "SLE",
    # --- SOS (Somali Shilling) ---
    "Sh": "SOS",
    # --- SRD (Surinamese Dollar) ---
    "SR$": "SRD",
    # --- SSP (South Sudanese Pound) ---
    "£_SSP": "SSP",
    # --- STN (São Tomé and Príncipe Dobra) ---
    "Db": "STN",
    # --- SVC (Salvadoran Colón) ---
    "₡_SVC": "SVC",
    # --- SYP (Syrian Pound) ---
    "£S": "SYP",
    # --- SZL (Eswatini Lilangeni) ---
    "E": "SZL",
    # --- THB (Thai Baht) ---
    "฿": "THB",
    # --- TJS (Tajikistani Somoni) ---
    "ЅМ": "TJS",
    # --- TMT (Turkmenistan Manat) ---
    "m": "TMT",
    # --- TND (Tunisian Dinar) ---
    "د.ت": "TND",
    # --- TOP (Tongan Paʻanga) ---
    "T$": "TOP",
    # --- TRY (Turkish Lira) ---
    "₺": "TRY",
    # --- TTD (Trinidad and Tobago Dollar) ---
    "TT$": "TTD",
    # --- TWD (New Taiwan Dollar) ---
    "NT$": "TWD",
    # --- TZS (Tanzanian Shilling) ---
    "TSh": "TZS",
    # --- UAH (Ukrainian Hryvnia) ---
    "₴": "UAH",
    # --- UGX (Ugandan Shilling) ---
    "USh": "UGX",
    # --- USD (United States Dollar) ---
    "$": "USD",
    # --- UYU (Uruguayan Peso) ---
    "$U": "UYU",
    # --- UZS (Uzbekistani Soʻm) ---
    "soʻm": "UZS",
    # --- VES (Venezuelan Bolívar Soberano) ---
    "Bs._VES": "VES",
    # --- VND (Vietnamese Đồng) ---
    "₫": "VND",
    # --- VUV (Vanuatu Vatu) ---
    "VT": "VUV",
    # --- WST (Samoan Tālā) ---
    "WS$": "WST",
    # --- XAF (Central African CFA Franc) ---
    "FCFA": "XAF",
    # --- XCD (East Caribbean Dollar) ---
    "EC$": "XCD",
    # --- XCG (Caribbean Guilder) ---
    "Cg.": "XCG",
    # --- XOF (West African CFA Franc) ---
    "CFA": "XOF",
    # --- XPF (CFP Franc) ---
    "CFP": "XPF",
    # --- YER (Yemeni Rial) ---
    "﷼_YER": "YER",
    # --- ZAR (South African Rand) ---
    "R": "ZAR",
    # --- ZMW (Zambian Kwacha) ---
    "ZK": "ZMW",
    # --- ZWG (Zimbabwe Gold) ---
    "ZiG": "ZWG",
}

# ---------------------------------------------------------------------------
# Pint Currency & Minor Unit Definitions
# ---------------------------------------------------------------------------

_ZERO_MINOR_UNIT: Final[frozenset[str]] = frozenset({
    "BIF", "CLP", "DJF", "GNF", "ISK", "JPY", "KMF", "KRW",
    "PYG", "RWF", "UGX", "VND", "VUV", "XAF", "XOF", "XPF",
})
_THREE_MINOR_UNIT: Final[frozenset[str]] = frozenset({
    "BHD", "IQD", "JOD", "KWD", "LYD", "OMR", "TND",
})
_NON_DECIMAL_MINOR_UNIT: Final[frozenset[str]] = frozenset({"MGA", "MRU"})

CURRENCY_MINOR_UNIT_EXPONENT: Final[dict[str, int]] = {
    code: (
        0 if code in _ZERO_MINOR_UNIT or code in _NON_DECIMAL_MINOR_UNIT
        else 3 if code in _THREE_MINOR_UNIT
        else 2
    )
    for code in CURRENCY_CODES
}

# Isolated currency dimensions prevent unhedged cross-currency operations
CURRENCY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = tuple(
    f"{code} = [currency_{code}]" for code in CURRENCY_CODES
)

CURRENCY_MINOR_UNIT_DEFINITIONS: Final[tuple[str, ...]] = tuple(
    f"{code}_minor = {10 ** -exponent:g} * {code}"
    for code, exponent in CURRENCY_MINOR_UNIT_EXPONENT.items()
    if exponent > 0
) + (
    "USD_cent = 0.01 * USD = cent_USD = cents_USD = cent = cents",
    "EUR_cent = 0.01 * EUR = cent_EUR = cents_EUR",
    "GBP_penny = 0.01 * GBP = GBp = penny_GBP = pence_GBP = pence",
    "MGA_iraimbilanja = 0.2 * MGA",
    "MRU_khoums = 0.2 * MRU",
)

# ---------------------------------------------------------------------------
# Master Commodity Definitions
# ---------------------------------------------------------------------------

COMMODITY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = (
    *CURRENCY_UNIT_DEFINITIONS,
    *CURRENCY_MINOR_UNIT_DEFINITIONS,

    # --- Energy Units (International Table Basis) ---
    "btu_it = 1055.05585262 * joule",
    "kBtu_it = 1000 * btu_it",
    "therm_it = 100000 * btu_it = therm = th",
    "dth_it = 10 * therm_it = dth = dekatherm",
    "MMBtu_it = 1000000 * btu_it = MMBtu = mmbtu = MMBTU",
    "BBtu_it = 1000000000 * btu_it = BBtu",
    "TBtu_it = 1000000000000 * btu_it = TBtu",
    "quad_it = 1000000000000000 * btu_it = quad",

    # --- US Energy Units (59°F Basis) ---
    "btu_us59 = 1054.804 * joule",
    "kBtu_us59 = 1000 * btu_us59",
    "therm_us59 = 100000 * btu_us59 = therm_us",
    "dth_us59 = 10 * therm_us59 = dth_us",
    "MMBtu_us59 = 1000000 * btu_us59 = MMBtu_us",
    "BBtu_us59 = 1000000000 * btu_us59",
    "TBtu_us59 = 1000000000000 * btu_us59",
    "quad_us59 = 1000000000000000 * btu_us59",

    # --- Gas Reference-Volume Units ---
    # US standard: 60°F (288.7056 K), 14.73 psia
    "scf_60F_14_73psia = [standard_gas_volume_60F_14_73psia] = scf = standard_cubic_foot",
    "ccf_60F_14_73psia = 100 * scf_60F_14_73psia = ccf = Ccf",
    "Mscf_60F_14_73psia = 1000 * scf_60F_14_73psia = Mscf = mscf = Mcf",
    "MMscf_60F_14_73psia = 1000000 * scf_60F_14_73psia = MMscf = mmscf = MMcf",
    "Bscf_60F_14_73psia = 1000000000 * scf_60F_14_73psia = Bscf = Bcf = bcf",
    "Tscf_60F_14_73psia = 1000000000000 * scf_60F_14_73psia = Tscf = Tcf = tcf",
    # Normal: 0°C (273.15 K), 101.325 kPa absolute, dry
    "Nm3_0C_101325kPa_dry = [normal_gas_volume_0C_101325kPa_dry] = Nm3 = normal_m3",
    "kNm3_0C_101325kPa_dry = 1000 * Nm3_0C_101325kPa_dry",
    "MNm3_0C_101325kPa_dry = 1000000 * Nm3_0C_101325kPa_dry",
    # ISO standard: 15°C (288.15 K), 101.325 kPa absolute, dry
    "Sm3_15C_101325kPa_dry = [standard_gas_volume_15C_101325kPa_dry] = Sm3 = standard_m3",
    "kSm3_15C_101325kPa_dry = 1000 * Sm3_15C_101325kPa_dry",
    "MSm3_15C_101325kPa_dry = 1000000 * Sm3_15C_101325kPa_dry",

    # --- Gas Flow Rates ---
    "scfd_60F = scf_60F_14_73psia / day = scfd",
    "Mscfd_60F = Mscf_60F_14_73psia / day = Mscfd = Mcfd",
    "MMscfd_60F = MMscf_60F_14_73psia / day = MMscfd = MMcfd",
    "Bcfd_60F = Bscf_60F_14_73psia / day = Bcfd",
    "Nm3_per_hour_0C = Nm3_0C_101325kPa_dry / hour = Nm3_h",
    "Nm3_per_day_0C = Nm3_0C_101325kPa_dry / day = Nm3_d",
    "kNm3_per_day_0C = kNm3_0C_101325kPa_dry / day",
    "Sm3_per_hour_15C = Sm3_15C_101325kPa_dry / hour = Sm3_h",
    "Sm3_per_day_15C = Sm3_15C_101325kPa_dry / day = Sm3_d",
    "kSm3_per_day_15C = kSm3_15C_101325kPa_dry / day",

    # --- Petroleum & Liquids ---
    "kbbl = 1000 * oil_barrel = Mbbl",
    "MMbbl = 1000000 * oil_barrel",
    "bpd = oil_barrel / day = bbl_d = barrels_per_day",
    "kbpd = 1000 * bpd = Mbpd = kbbl_d",
    "mmbpd = 1000000 * bpd = MMbpd",
    "m3_per_day = meter ** 3 / day = cbm_d",
    "m3_per_hour = meter ** 3 / hour = cbm_h",

    # --- Mass, Weights & Multiples ---
    "mt = 1000 * kilogram = metric_tonne = tonne",
    "commodity_short_hundredweight = 100 * pound = cwt = short_cwt",
    "commodity_long_hundredweight = 112 * pound = long_cwt",
    "commodity_short_ton = 2000 * pound = short_ton",
    "commodity_long_ton = 2240 * pound = long_ton",
    "toz = 31.1034768 * gram = troy_ounce = troy_oz",

    # --- Annualised Mining, LNG & Freight Flow Rates ---
    "trading_year_365 = 365 * day",
    "tpa = tonne / trading_year_365",
    "ktpa = 1000 * tpa = KTPA",
    "mtpa = 1000000 * tpa = MTPA",
    "tpd = tonne / day",

    # --- Agricultural Measures ---
    "commodity_peck = 0.25 * bushel = peck",
    "coffee_bag_60kg = 60 * kilogram = coffee_bag",
    "cocoa_bag_62_5kg = 62.5 * kilogram = cocoa_bag",

    # --- Environmental & Emissions Instruments ---
    "tonne_CO2 = tonne = tCO2", 
    "tonne_CO2e = tonne = tCO2e",
    "EU_allowance = [EU_allowance] = EUA",
    "EU_aviation_allowance = [EU_aviation_allowance] = EUAA",
    "UK_allowance = [UK_allowance] = UKA",
    "California_carbon_allowance = [California_carbon_allowance] = CCA",
)

# ---------------------------------------------------------------------------
# Registry Singleton & Serialisation Bindings
# ---------------------------------------------------------------------------

class CommodityUnitRegistry:
    """Thread-safe singleton managing the master Pint UnitRegistry."""

    _instance: pint.UnitRegistry | None = None
    _lock: threading.Lock = threading.Lock()

    @classmethod
    def get_registry(cls) -> pint.UnitRegistry:
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    registry = pint.UnitRegistry(
                        autoconvert_offset_to_baseunit=True,
                        on_redefinition="raise",
                    )
                    for definition in COMMODITY_UNIT_DEFINITIONS:
                        registry.define(definition)

                    # Bind registry to global process state to ensure unpickling in
                    # multiprocessing pools restores custom dimensions seamlessly.
                    pint.set_application_registry(registry)
                    cls._instance = registry
        return cls._instance


ureg: Final[pint.UnitRegistry] = CommodityUnitRegistry.get_registry()
Q_ = ureg.Quantity

# ---------------------------------------------------------------------------
# Parsing and Normalisation Infrastructure
# ---------------------------------------------------------------------------

_SAFE_SYMBOLS: Final[dict[str, str]] = {
    symbol: code
    for symbol, code in CURRENCY_SYMBOL_TO_CODE.items()
    if any(not (char.isalnum() or char == "_") for char in symbol)
}
_SORTED_SAFE_SYMBOLS: Final[tuple[str, ...]] = tuple(
    sorted(_SAFE_SYMBOLS, key=len, reverse=True)
)

_QUOTE_REPLACEMENTS: Final[dict[str, str]] = {
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
    # Metals
    "$/toz": "USD / toz",
    "$/oz": "USD / toz",
    # Carbon & Emissions
    "€/t": "EUR / tonne_CO2",
    "€/EUA": "EUR / EU_allowance",
    "£/t": "GBP / tonne_CO2",
    "£/UKA": "GBP / UK_allowance",
    "$/CCA": "USD / California_carbon_allowance",
}

_UNIT_TOKEN_REPLACEMENTS: Final[tuple[tuple[re.Pattern[str], str], ...]] = (
    (re.compile(r"(?<![A-Za-z0-9_])bbl(?![A-Za-z0-9_])", re.I), "oil_barrel"),
    (re.compile(r"(?<![A-Za-z0-9_])MMBtu(?![A-Za-z0-9_])", re.I), "MMBtu_it"),
    (re.compile(r"(?<![A-Za-z0-9_])cbm(?![A-Za-z0-9_])", re.I), "meter ** 3"),
)

_PREFIX_CURRENCY_PATTERN = re.compile(
    r"^\s*(?P<sign>[+-]?)\s*"
    r"(?P<symbol>" + "|".join(re.escape(sym) for sym in _SORTED_SAFE_SYMBOLS) + r")\s*"
    r"(?P<magnitude>\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)\s*"
    r"(?:/\s*(?P<denom>.+))?$"
)

_NUMERIC_MAGNITUDE_PATTERN = re.compile(
    r"^\s*(?P<magnitude>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(?P<unit>.+?)\s*$"
)


def _resolve_currency_symbol(symbol: str, default_currency: str | None) -> str:
    code = _SAFE_SYMBOLS[symbol]
    if default_currency is None or symbol not in {"$", "¥", "£"}:
        return code

    if default_currency not in CURRENCY_CODES:
        raise ValueError(f"Unknown ISO 4217 code: {default_currency}")

    valid_codes = {
        "$": {
            "USD", "ARS", "AUD", "BBD", "BMD", "BND", "BSD", "BZD",
            "CAD", "CLP", "COP", "CUP", "FJD", "GYD", "HKD", "JMD",
            "KYD", "LRD", "MXN", "NAD", "NZD", "SGD", "SRD", "TTD",
            "TWD", "UYU", "XCD", "WST",
        },
        "¥": {"JPY", "CNY"},
        "£": {"GBP", "FKP", "GIP", "SHP", "SSP"},
    }[symbol]

    if default_currency not in valid_codes:
        raise ValueError(f"{symbol!r} is not a valid symbol for {default_currency}")

    return default_currency


def sanitize_unit_string(raw_unit: str, default_currency: str | None = None) -> str:
    """Normalises market quotation strings into Pint registry notation."""
    text = raw_unit.strip()

    if text in _QUOTE_REPLACEMENTS:
        return _QUOTE_REPLACEMENTS[text]

    for symbol in _SORTED_SAFE_SYMBOLS:
        if text.startswith(symbol) and (
            len(text) == len(symbol) or text[len(symbol)] in {"/", " "}
        ):
            text = _resolve_currency_symbol(symbol, default_currency) + text[len(symbol):]
            break

    if text.startswith("p/"):
        text = "GBP_penny/" + text[2:]
    elif text.startswith("c/"):
        text = "USD_cent/" + text[2:]
    elif text.startswith("cents/"):
        text = "USD_cent/" + text[6:]

    for pattern, replacement in _UNIT_TOKEN_REPLACEMENTS:
        text = pattern.sub(replacement, text)

    return text


def parse_quantity(
    value: int | float | Decimal | str,
    unit: str | None = None,
    *,
    default_currency: str | None = None,
    as_decimal: bool = True,
) -> pint.Quantity:
    """Parses numeric expressions or quotation strings into validated quantities.

    Supports currency prefixes ($75.50/bbl), standard spacing (75.50 USD/bbl),
    and scientific notation (1.5e6 scfd).
    """
    if isinstance(value, str) and unit is None:
        # Match leading currency symbols: "$75.50/bbl", "-£45.20 / MWh"
        prefix_match = _PREFIX_CURRENCY_PATTERN.match(value)
        if prefix_match:
            sign = prefix_match.group("sign") or ""
            symbol = prefix_match.group("symbol")
            mag_str = sign + prefix_match.group("magnitude")
            denom = prefix_match.group("denom")

            currency_code = _resolve_currency_symbol(symbol, default_currency)
            unit_str = f"{currency_code} / {denom}" if denom else currency_code
            sanitized = sanitize_unit_string(unit_str, default_currency)
            magnitude = Decimal(mag_str) if as_decimal else float(mag_str)
            return ureg.Quantity(magnitude, sanitized)

        # Match standard numeric magnitude with trailing unit: "75.50 USD/bbl", "1.5e6 scfd"
        num_match = _NUMERIC_MAGNITUDE_PATTERN.match(value)
        if num_match:
            mag_str = num_match.group("magnitude")
            unit_str = num_match.group("unit")
            sanitized = sanitize_unit_string(unit_str, default_currency)
            magnitude = Decimal(mag_str) if as_decimal else float(mag_str)
            return ureg.Quantity(magnitude, sanitized)

        # Dimension-only string or unit with implicit unit magnitude
        return ureg.Quantity(sanitize_unit_string(value, default_currency))

    if unit is None:
        return ureg.Quantity(value, ureg.dimensionless)

    sanitized = sanitize_unit_string(unit, default_currency)
    return ureg.Quantity(value, sanitized)


# ---------------------------------------------------------------------------
# Institutional Validation Helpers
# ---------------------------------------------------------------------------

def validate_dimensionality(
    quantity: pint.Quantity,
    expected_dimension: str | pint.Unit | pint.Quantity,
) -> bool:
    """Asserts that a quantity matches an expected physical dimension or unit context."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

    if isinstance(expected_dimension, (pint.Unit, pint.Quantity)):
        return quantity.check(expected_dimension.dimensionality)

    parsed_expected = ureg(sanitize_unit_string(expected_dimension))
    return quantity.check(parsed_expected.dimensionality)


def ensure_unit(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
) -> pint.Quantity:
    """Converts a quantity to a target unit, raising a clear DimensionalityError on mismatch."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

    sanitized_target = (
        sanitize_unit_string(target_unit) if isinstance(target_unit, str) else target_unit
    )
    return quantity.to(sanitized_target)