import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest

import backend_bridge
import units

KEYS = [p["key"] for p in backend_bridge.STATE_PROPERTIES]


@pytest.mark.parametrize(["kind", "display", "label", "factor"], [
    ("temperature_per_pressure", {"pressure": "bar"}, "K/bar", 1e5),
    ("temperature_per_pressure", {"pressure": "MPa"}, "K/MPa", 1e6),
    ("inverse_pressure", {"pressure": "bar"}, "1/bar", 1e5),
    ("pressure_per_temperature", {"pressure": "kPa"}, "kPa/K", 1e-3),
    ("pressure_per_density", {"pressure": "bar", "density": "kg/m3"}, "bar·m³/kg", 1e-5),
    ("pressure_per_density", {"pressure": "bar", "density": "g/cm3"}, "bar·cm³/g", 1e-2),
    ("pressure_per_density_squared", {"pressure": "MPa", "density": "g/cm3"}, "MPa·cm⁶/g²", 1.0),
    ("inverse_density", {"density": "g/cm3"}, "cm³/g", 1e3),
    ("inverse_density_squared", {"density": "kg/m3"}, "m⁶/kg²", 1.0),
    ("inverse_density_squared", {"density": "g/cm3"}, "cm⁶/g²", 1e6),
    ("pressure", {"pressure": "MPa"}, "MPa", 1e-6),
    ("speed", {}, "m/s", 1.0),
    ("inverse_temperature", {}, "1/K", 1.0),
])
def test_units_follow_the_chosen_pressure_and_density(kind, display, label, factor):
    assert units.display_label(kind, display) == label
    assert units.display_value(kind, display, 2.5) == pytest.approx(2.5 * factor, rel=1e-12)


def test_every_property_has_a_unit():
    for prop in backend_bridge.STATE_PROPERTIES + backend_bridge.SATURATION_PROPERTIES:
        assert units.display_label(prop["kind"], units.DEFAULTS)


def test_state_carries_every_property():
    nitrogen = backend_bridge.state("nitrogen", 300.0, p_Pa=101325.0)
    assert all(isinstance(nitrogen[k], float) and math.isfinite(nitrogen[k]) for k in KEYS)
    air = backend_bridge.state("air", 300.0, p_Pa=101325.0)
    assert air["fugacity_coefficient"] is None
    assert backend_bridge.defined("nitrogen", "fugacity_coefficient")
    assert not backend_bridge.defined("air", "fugacity_coefficient")


def test_saturated_pair_of_a_pure_fluid_shares_gibbs_energy_and_fugacity():
    pair = backend_bridge.saturated_pair("nitrogen", 100.0)
    liquid, vapor = pair["liquid"], pair["vapor"]
    assert all(k in liquid and k in vapor for k in KEYS)
    assert liquid["gibbs_energy_J_kg"] == pytest.approx(vapor["gibbs_energy_J_kg"], rel=1e-9)
    assert liquid["fugacity_coefficient"] == pytest.approx(vapor["fugacity_coefficient"], rel=1e-9)


def test_scan_rows_carry_every_property():
    scan = backend_bridge.scan("nitrogen", "isobar", 1e5, 200.0, 400.0)
    rows = [row for row in scan["rows"] if row["status"] == "ok"]
    assert rows and all(k in row for row in rows for k in KEYS)


def test_scan_export_carries_every_property_the_fluid_defines():
    import app

    scan = backend_bridge.scan("air", "isotherm", 300.0, 1e5, 1e6)
    columns = app._defined("air", backend_bridge.STATE_PROPERTIES)
    text = app._scan_export(scan, backend_bridge.limits("air"), units.DEFAULTS, "K", "bar",
                            columns, ";", ",")
    symbols, units_row = (line.split(";") for line in text.splitlines()[4:6])
    assert symbols == ["P"] + [p["symbol"] for p in columns] + ["Nota"]
    assert "φ" not in symbols and "μ_JT" in symbols
    assert units_row[symbols.index("μ_JT")] == "K/bar"


def test_values_hang_their_sign():
    import app

    sign, digits, _, _ = app._value_cell(-1741.326)
    assert (sign.children, digits.children) == ("-", "1741")
    sign, digits, _, _ = app._value_cell(0.5)
    assert (sign.children, digits.children) == ("", "0")
