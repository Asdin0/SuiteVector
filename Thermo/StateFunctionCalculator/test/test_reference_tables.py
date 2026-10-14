import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, os.path.join(project_root, 'general_class'))
sys.path.insert(0, os.path.join(project_root, 'fluid_list'))

import math as mth
import pytest
import nitrogen
import air

NITROGEN = nitrogen.Nitrogen()
AIR = air.Air()


def printed(text):
    value = float(text)
    digits = text.lower().split("e")[0].lstrip("+-").replace(".", "").lstrip("0")
    last_digit = 10 ** (mth.floor(mth.log10(abs(value))) - len(digits) + 1)
    return value, 0.75 * last_digit


single_phase_states = [
    ("nitrogen", 300, 0.1, "0.040098", "8717.7", "191.77", "gas"),
    ("nitrogen", 1000, 0.1, "0.012023", "30135", "228.17", "gas"),
    ("nitrogen", 80, 10, "29.214", "-3069.6", "79.486", "liquid"),
    ("nitrogen", 150, 10, "14.963", "1359.2", "118.44", "supercritical"),
    ("nitrogen", 300, 50, "14.697", "7645.6", "135.10", "supercritical"),
    ("nitrogen", 200, 200, "29.542", "6809.1", "105.72", "supercritical"),
    ("nitrogen", 500, 1000, "35.595", "37163", "124.12", "supercritical"),
    ("nitrogen", 65, 0.1, "30.690", "-4117.0", "69.556", "liquid"),
    ("nitrogen", 120, 2, "2.8387", "2571.8", "134.71", "vapor"),
    ("air", 300, 0.101325, "0.040634", "8696.5", "198.85", "gas"),
    ("air", 100, 1, "26.593", "-2413.1", "99.798", "liquid"),
    ("air", 150, 5, "6.4975", "2669.9", "137.96", "supercritical"),
    ("air", 200, 200, "30.660", "6334.8", "112.70", "supercritical"),
    ("air", 1000, 10, "1.1655", "30484", "197.27", "supercritical"),
    ("air", 2000, 2000, "32.893", "1.1424e5", "176.45", "supercritical"),
]

saturation_states = [
    ("nitrogen", 65, "0.01740", "30.685", "-4118.8", "69.569", "0.01740", "0.03259", "1863.9", "161.61"),
    ("nitrogen", 77.355, "0.101325", "28.775", "-3418.2", "79.395", "0.101325", "0.16464", "2161.5", "151.53"),
    ("nitrogen", 100, "0.77827", "24.608", "-2050.8", "94.576", "0.77827", "1.1409", "2458.6", "139.67"),
    ("nitrogen", 120, "2.51058", "18.682", "-500.60", "107.89", "2.51058", "4.4653", "2077.8", "129.38"),
    ("nitrogen", 125, "3.20687", "15.210", "179.37", "113.10", "3.20687", "7.3244", "1541.7", "124.00"),
    ("air", 80, "0.11462", "30.042", "-3589.6", "87.037", "0.08232", "0.12813", "2243.6", "161.61"),
    ("air", 100, "0.66313", "26.533", "-2416.0", "99.896", "0.56742", "0.78609", "2564.0", "150.36"),
    ("air", 120, "2.15573", "21.588", "-1021.2", "112.00", "2.00674", "2.9889", "2458.3", "141.22"),
    ("air", 132, "3.72284", "14.198", "478.83", "123.13", "3.64625", "8.1273", "1553.9", "131.33"),
]

FLUIDS = {"nitrogen": NITROGEN, "air": AIR}


@pytest.mark.parametrize(["fluid_key", "t", "p_MPa", "density_text", "enthalpy_text", "entropy_text", "phase"], single_phase_states)
def test_single_phase_against_paper_tables(fluid_key, t, p_MPa, density_text, enthalpy_text, entropy_text, phase):
    fluid = FLUIDS[fluid_key]
    molar_mass = fluid.constants.molar_mass_kg_mol
    result = fluid.calculated_functions(t, p_MPa * 1e6)
    density, density_tolerance = printed(density_text)
    enthalpy, enthalpy_tolerance = printed(enthalpy_text)
    entropy, entropy_tolerance = printed(entropy_text)
    assert result.phase == phase
    assert result.density_kg_m3 / (molar_mass * 1000) == pytest.approx(density, abs=density_tolerance)
    assert result.enthalpy_J_kg * molar_mass == pytest.approx(enthalpy, abs=enthalpy_tolerance)
    assert result.entropy_J_kg_K * molar_mass == pytest.approx(entropy, abs=entropy_tolerance)


@pytest.mark.parametrize(["fluid_key", "t", "liquid_pressure_text", "liquid_density_text", "liquid_enthalpy_text", "liquid_entropy_text",
                          "vapor_pressure_text", "vapor_density_text", "vapor_enthalpy_text", "vapor_entropy_text"], saturation_states)
def test_saturation_against_paper_tables(fluid_key, t, liquid_pressure_text, liquid_density_text, liquid_enthalpy_text, liquid_entropy_text,
                                         vapor_pressure_text, vapor_density_text, vapor_enthalpy_text, vapor_entropy_text):
    fluid = FLUIDS[fluid_key]
    molar_mass = fluid.constants.molar_mass_kg_mol
    saturation = fluid.saturation_state(t)
    sides = (
        ("liquid", saturation.liquid_pressure_Pa, saturation.liquid_density_kg_m3,
         liquid_pressure_text, liquid_density_text, liquid_enthalpy_text, liquid_entropy_text),
        ("vapor", saturation.vapor_pressure_Pa, saturation.vapor_density_kg_m3,
         vapor_pressure_text, vapor_density_text, vapor_enthalpy_text, vapor_entropy_text),
    )
    for phase, pressure_Pa, density_kg_m3, pressure_text, density_text, enthalpy_text, entropy_text in sides:
        pressure, pressure_tolerance = printed(pressure_text)
        density, density_tolerance = printed(density_text)
        enthalpy, enthalpy_tolerance = printed(enthalpy_text)
        entropy, entropy_tolerance = printed(entropy_text)
        result = fluid.calculated_functions(t, density=density_kg_m3)
        assert result.phase == phase
        assert pressure_Pa / 1e6 == pytest.approx(pressure, abs=pressure_tolerance)
        assert density_kg_m3 / (molar_mass * 1000) == pytest.approx(density, abs=density_tolerance)
        assert result.enthalpy_J_kg * molar_mass == pytest.approx(enthalpy, abs=enthalpy_tolerance)
        assert result.entropy_J_kg_K * molar_mass == pytest.approx(entropy, abs=entropy_tolerance)
