"""Fixed chemistry constants and model scaling for Stage 23."""
from __future__ import annotations

from decimal import Decimal
from typing import Any


SOIL_RAW = {
    "S": {
        "soil": "Sand (Rehovot)",
        "Cmic_mg_kg": 199.44,
        "LOI_SOM_g_g": 0.016732,
        "Corg_g_kg": 9.705,
        "core_BD_g_cm3": 1.3947,
        "Corg_mg_cm3_reported": 13.535,
        "Cmic_ug_cm3_reported": 278.16,
    },
    "V": {
        "soil": "Vertisol (Bnei Re'em)",
        "Cmic_mg_kg": 406.61,
        "LOI_SOM_g_g": 0.048627,
        "Corg_g_kg": 28.204,
        "core_BD_g_cm3": 1.3482,
        "Corg_mg_cm3_reported": 38.023,
        "Cmic_ug_cm3_reported": 548.17,
    },
    "L": {
        "soil": "Loess (Mishmar HaNegev)",
        "Cmic_mg_kg": 315.57,
        "LOI_SOM_g_g": 0.027252,
        "Corg_g_kg": 15.806,
        "core_BD_g_cm3": 1.1344,
        "Corg_mg_cm3_reported": 17.931,
        "Cmic_ug_cm3_reported": 357.99,
    },
}

MODEL_REFERENCE_OM_DENSITY = 1.0
MODEL_REFERENCE_B0_DENSITY = 0.0125


def chemistry_table() -> dict[str, dict[str, Any]]:
    """Return fixed values and derive volume concentrations and model densities."""
    volumes = {}
    for code, raw in SOIL_RAW.items():
        volumes[code] = {
            "Corg_mg_cm3": raw["Corg_g_kg"] * raw["core_BD_g_cm3"],
            "Cmic_ug_cm3": raw["Cmic_mg_kg"] * raw["core_BD_g_cm3"],
        }
    mean_corg = sum(row["Corg_mg_cm3"] for row in volumes.values()) / len(volumes)
    mean_cmic = sum(row["Cmic_ug_cm3"] for row in volumes.values()) / len(volumes)
    result = {}
    for code, raw in SOIL_RAW.items():
        row = dict(raw)
        row.update(volumes[code])
        row["OM_density"] = (
            MODEL_REFERENCE_OM_DENSITY * row["Corg_mg_cm3"] / mean_corg
        )
        row["B0_density"] = (
            MODEL_REFERENCE_B0_DENSITY * row["Cmic_ug_cm3"] / mean_cmic
        )
        result[code] = row
    result["chem_S"] = result["S"]
    result["chem_V"] = result["V"]
    result["chem_L"] = result["L"]
    result["chem_eq"] = {
        "soil": "Stage 22 equal chemistry",
        "OM_density": MODEL_REFERENCE_OM_DENSITY,
        "B0_density": MODEL_REFERENCE_B0_DENSITY,
    }
    return result


CHEMISTRY = chemistry_table()

REPORTED_MODEL_DENSITIES = {
    "S": {"OM_density": 0.5843, "B0_density": 0.008807},
    "V": {"OM_density": 1.6415, "B0_density": 0.017357},
    "L": {"OM_density": 0.7741, "B0_density": 0.011335},
}


def validate_chemistry() -> None:
    """Assert prompt-derived concentrations and scaling to four significant figures."""
    for code, raw in SOIL_RAW.items():
        row = CHEMISTRY[code]
        checks = (
            (row["Corg_mg_cm3"], raw["Corg_mg_cm3_reported"]),
            (row["Cmic_ug_cm3"], raw["Cmic_ug_cm3_reported"]),
            (row["OM_density"], REPORTED_MODEL_DENSITIES[code]["OM_density"]),
            (row["B0_density"], REPORTED_MODEL_DENSITIES[code]["B0_density"]),
        )
        for actual, expected in checks:
            expected_decimal = Decimal(str(expected))
            tolerance = Decimal("0.5") * Decimal(1).scaleb(
                expected_decimal.adjusted() - 3
            )
            if abs(Decimal(str(actual)) - expected_decimal) > tolerance:
                raise AssertionError(
                    f"{code} chemistry derivation mismatch: {actual:.12g} vs {expected:.12g}"
                )


def code_for_soil(soil: str) -> str:
    normalized = soil.lower()
    if normalized == "sand":
        return "S"
    if normalized == "vertisol":
        return "V"
    if normalized == "loess":
        return "L"
    raise ValueError(f"Unknown soil: {soil}")
