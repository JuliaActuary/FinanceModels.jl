"""Extract cached official EIOPA values without running macros or calculations.

Usage: python3 test/data/extract_eiopa.py path/to/official.xlsm
Requires openpyxl only for fixture maintenance, never for Julia tests.
"""
import hashlib
import json
from pathlib import Path
import sys

from openpyxl import load_workbook

source = Path(sys.argv[1])
root = Path(__file__).resolve().parents[2]
book = load_workbook(source, read_only=True, data_only=True)
basic = book["Input Data & Extrapolation"]
rows = list(basic.iter_rows(min_row=25, max_row=174, min_col=7, max_col=12, values_only=True))
selected = [row for row in rows if row[1] == 1]

def vector(sheet, column, start, length=150):
    return [row[0] for row in book[sheet].iter_rows(min_row=start, max_row=start+length-1,
                                                  min_col=column, max_col=column, values_only=True)]

data = {
    "source": {
        "url": "https://www.eiopa.europa.eu/document/download/01e0930c-d78d-4806-8e83-a27609ea57c5_en",
        "version": "RFR extrapolation and VA calculation (19 May 2026)",
        "copyright": "EIOPA, 2026. Extracted numerical data; no macros executed.",
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "ranges": ["Input Data & Extrapolation!G25:L174", "CSSR Calculation!L24:L173",
                   "VA Term Structures!L23:L172", "IRR Up Shock!J22:O171",
                   "IRR Down Shock!J22:O171", "Parameters!G17:L166",
                   "CSSR Calculation!H10:L20", "VA Term Structures!H10:L20",
                   "IRR Up Shock!O19", "IRR Down Shock!N19"],
    },
    "inputs": {
        "rates": [r[3] for r in selected], "tenors": [r[0] for r in selected],
        "llfr_weights": [r[2] for r in selected], "fsp": basic["H13"].value,
        "ufr": basic["H14"].value, "alpha": basic["H15"].value,
        "cra_bp": basic["H16"].value, "frequency": basic["H12"].value,
    },
    "va": {"rcs_currency": basic["L10"].value, "rcs_country": basic["L11"].value,
           "debt_to_assets": basic["L12"].value, "usa": book["VA Term Structures"]["H12"].value,
           "usa_selected": bool(book["VA Term Structures"]["H11"].value),
           "macro_selected": bool(book["VA Term Structures"]["H10"].value),
           "assets": book["CSSR Calculation"]["L10"].value,
           "assets_notional": book["CSSR Calculation"]["L11"].value,
           "liabilities": book["CSSR Calculation"]["L13"].value,
           "liabilities_notional": book["CSSR Calculation"]["L14"].value},
    "expected": {
        "bootstrap": [r[4] for r in rows[:50]], "llfr": basic["L22"].value,
        "cssr_value": book["CSSR Calculation"]["L16"].value,
        "total_va": book["VA Term Structures"]["L13"].value,
        "up_llfr": book["IRR Up Shock"]["L19"].value,
        "down_llfr": book["IRR Down Shock"]["K19"].value,
        "cssr_llfr": book["CSSR Calculation"]["K20"].value,
        "va_llfr": book["VA Term Structures"]["K20"].value,
        "up_va_llfr": book["IRR Up Shock"]["O19"].value,
        "down_va_llfr": book["IRR Down Shock"]["N19"].value,
        "up_unfloored_prefix": vector("IRR Up Shock", 10, 22, int(basic["H13"].value)),
        "down_unfloored_prefix": vector("IRR Down Shock", 10, 22, int(basic["H13"].value)),
        "assets_pvbp": book["CSSR Calculation"]["L12"].value,
        "liabilities_pvbp": book["CSSR Calculation"]["L15"].value,
        "notional_va": book["CSSR Calculation"]["H10"].value,
        "currency_va": book["VA Term Structures"]["L10"].value,
        "country_weight": book["VA Term Structures"]["L11"].value,
        "basic": vector("Input Data & Extrapolation", 12, 25),
        "cssr": vector("CSSR Calculation", 12, 24), "va": vector("VA Term Structures", 12, 23),
        "up": vector("IRR Up Shock", 12, 22), "up_va": vector("IRR Up Shock", 15, 22),
        "down": vector("IRR Down Shock", 12, 22), "down_va": vector("IRR Down Shock", 15, 22),
    },
    "shocks": {name: vector("Parameters", col, 17) for name, col in
               [("up_relative", 8), ("up_absolute", 9), ("down_relative", 10), ("down_absolute", 11), ("floor", 12)]},
}
with (root / "test/data/eiopa_2026.toml").open("w") as f:
    for index, (section, items) in enumerate(data.items()):
        f.write(f"[{section}]\n")
        for key, value in items.items():
            f.write(f"{key} = {json.dumps(value)}\n")
        if index + 1 < len(data):
            f.write("\n")

with (root / "src/EIOPA/ShockData.jl").open("w") as f:
    f.write("# EIOPA, 2026. Articles 165–166 as amended by EU 2026/269.\n")
    f.write("# Source: https://eur-lex.europa.eu/eli/reg_del/2026/269/oj/eng\n")
    f.write("# Cross-checked against the official workbook's Parameters!H17:K66.\n")
    f.write("# Regenerate with test/data/extract_eiopa.py; fixture records source SHA-256.\n")
    for name in ("up_relative", "up_absolute", "down_relative", "down_absolute"):
        f.write(f"const _{name.upper()} = (\n")
        values = data["shocks"][name][:50]
        for offset in range(0, 50, 10):
            f.write("    " + ", ".join(repr(x) for x in values[offset:offset+10]) + ",\n")
        f.write(")\n")
