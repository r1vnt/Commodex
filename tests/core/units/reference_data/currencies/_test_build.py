"""Compile locally downloaded SIX XML and reviewed overlays into a pinned snapshot.

Test helper: independently rebuild a temporary currency snapshot from supplied XML.
Download list-one.xml and list-three.xml from SIX over verified HTTPS first.
This build never downloads data or guesses missing historical dates/precision.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def build(directory):
    directory = Path(directory)
    overlay = json.loads((directory / "overlay.json").read_text(encoding="utf-8"))
    version = overlay["version"]
    records = {}
    sources = []
    for stem, status, tag in (("list-one", "current", "CcyNtry"), ("list-three", "historical", "HstrcCcyNtry")):
        filename = stem + ".xml"
        raw = (directory / filename).read_bytes()
        root = ET.fromstring(raw)
        sources.append({"file": filename, "url": "https://www.six-group.com/dam/download/financial-information/data-center/iso-currrency/lists/" + filename,
                        "sha256": sha256(raw).hexdigest(), "published": root.attrib.get("Pblshd")})
        for entry in root.iter(tag):
            values = {element.tag: element.text for element in entry}
            code = values.get("Ccy")
            if not code:
                continue
            numeric = values.get("CcyNbr")
            numeric = numeric.zfill(3) if numeric else None
            key = (code, numeric, values.get("CcyNm"), status)
            name_element = entry.find("CcyNm")
            fund = name_element is not None and name_element.attrib.get("IsFund") == "true"
            category = "fund" if fund else "currency"
            category = {"XDR": "unit_of_account", "XAU": "precious_metal", "XAG": "precious_metal",
                        "XPT": "precious_metal", "XPD": "precious_metal", "XXX": "no_currency", "XTS": "testing"}.get(code, category)
            record = records.setdefault(key, {
                "code": code, "numeric_code": numeric, "name": values.get("CcyNm"), "status": status,
                "minor_digits": int(values["CcyMnrUnts"]) if values.get("CcyMnrUnts", "").isdigit() else None,
                "entities": [], "withdrawal_dates": [], "valid_from": None, "valid_to": None,
                "source_file": filename, "kind": "iso4217", "category": category, "is_fund": fund,
            })
            if values.get("CtryNm") not in record["entities"]:
                record["entities"].append(values.get("CtryNm"))
            withdrawal = values.get("WthdrwlDt")
            if withdrawal and withdrawal not in record["withdrawal_dates"]:
                record["withdrawal_dates"].append(withdrawal)
    # List Three includes country/entity withdrawals of currencies that remain
    # active elsewhere (for example EUR in Serbia and Montenegro). These are
    # entity history, not withdrawn monetary denominations.
    current_keys = {(record["code"], record["numeric_code"]) for record in records.values() if record["status"] == "current"}
    entity_history, denominations = [], []
    for record in records.values():
        if record["status"] == "historical" and (record["code"], record["numeric_code"]) in current_keys:
            entity_history.append(record)
        else:
            denominations.append(record)
    for record in denominations:
        record.update(overlay.get("record_overrides", {}).get(record["code"], {}))
    result = {"schema_version": 1, "version": version, "observed_on": version,
              "records": denominations + overlay.get("additional_records", []), "entity_history": entity_history}
    (directory / "reference.json").write_text(json.dumps(result, ensure_ascii=True, indent=2) + "\n", encoding="utf-8", newline="\n")
    manifest = {"version": version, "retrieved_on": version, "sources": sources,
                "derived_sha256": {name: sha256((directory / name).read_bytes()).hexdigest()
                                   for name in ("reference.json", "overlay.json", "conventions.json")}}
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    build(parser.parse_args().directory)
