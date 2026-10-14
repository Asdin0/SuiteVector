import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, os.path.join(project_root, 'general_class'))
sys.path.insert(0, os.path.join(project_root, 'fluid_list'))

import pytest
import GeneralClass
import nitrogen
import air

NITROGEN = nitrogen.Nitrogen()
AIR = air.Air()
FLUIDS = {"nitrogen": NITROGEN, "air": AIR}


def test_constants_match_papers():
    assert NITROGEN.constants.critical_pressure_Pa == 3395800
    assert NITROGEN.constants.universal_gas_constant_J_mol_K == 8.314510
    assert AIR.constants.universal_gas_constant_J_mol_K == 8.314510
    assert NITROGEN.real_helmholtz_terms[26].n == 0.0368825891208
    assert AIR.constants.reference_pressure_Pa == 3785020
    assert NITROGEN.constants.reference_density_kg_m3 == pytest.approx(11.1839 * 28.01348, rel=1e-12)
    assert AIR.constants.reference_density_kg_m3 == pytest.approx(10.4477 * 28.9586, rel=1e-12)


def test_liquid_just_above_saturation_pressure_is_liquid():
    result = NITROGEN.calculated_functions(77.355, 101355.5)
    assert result.phase == "liquid"
    assert result.density_kg_m3 == pytest.approx(806.085, abs=0.01)
    assert result.enthalpy_J_kg == pytest.approx(-122018.3, abs=2)


def test_compressed_liquid_near_critical_point_is_liquid():
    result = NITROGEN.calculated_functions(125.9, 3.35e6)
    assert result.phase == "liquid"
    assert result.density_kg_m3 == pytest.approx(384.1614, abs=1e-3)


def test_vapor_just_below_saturation_pressure_is_vapor():
    assert NITROGEN.calculated_functions(110, 1465790).phase == "vapor"


def test_nitrogen_saturation_from_maxwell_criterion():
    saturation = NITROGEN.saturation_state(77.355)
    assert saturation.vapor_pressure_Pa == pytest.approx(101325, abs=1)
    assert saturation.liquid_pressure_Pa == saturation.vapor_pressure_Pa
    assert saturation.liquid_density_kg_m3 == pytest.approx(806.085, abs=1e-3)
    assert saturation.vapor_density_kg_m3 == pytest.approx(4.6121, abs=1e-4)


@pytest.mark.parametrize(["fluid_key", "t", "density"], [("nitrogen", 100, 300.0), ("nitrogen", 77.355, 100.0), ("air", 100, 300.0)])
def test_density_inside_two_phase_region_is_refused(fluid_key, t, density):
    with pytest.raises(GeneralClass.TwoPhaseError):
        FLUIDS[fluid_key].calculated_functions(t, density=density)


@pytest.mark.parametrize(["fluid_key", "t"], [("nitrogen", 65), ("nitrogen", 77.355), ("nitrogen", 120), ("air", 65), ("air", 70), ("air", 80), ("air", 130)])
def test_saturated_density_gives_saturation_pressure(fluid_key, t):
    fluid = FLUIDS[fluid_key]
    saturation = fluid.saturation_state(t)
    assert fluid.calculated_functions(t, density=saturation.liquid_density_kg_m3).phase == "liquid"
    assert fluid.calculated_functions(t, density=saturation.vapor_density_kg_m3).phase == "vapor"
    assert fluid.pressure(t, saturation.liquid_density_kg_m3) == pytest.approx(saturation.liquid_pressure_Pa, rel=1e-6)
    assert fluid.pressure(t, saturation.vapor_density_kg_m3) == pytest.approx(saturation.vapor_pressure_Pa, rel=1e-6)


@pytest.mark.parametrize(["fluid_key", "t", "p", "density"], [("nitrogen", 150, 800e6, None), ("nitrogen", 65, None, 880.0),
                                                               ("air", 150, 800e6, None), ("air", 200, 1500e6, None)])
def test_states_above_melting_line_are_refused(fluid_key, t, p, density):
    with pytest.raises(GeneralClass.SolidPhaseError):
        FLUIDS[fluid_key].calculated_functions(t, p, density)


def test_nitrogen_melting_pressure():
    assert NITROGEN.meltingpoint_pressure(63.151) == pytest.approx(12523, abs=1e-6)
    assert NITROGEN.meltingpoint_pressure(150) == pytest.approx(12523 * (1 + 12798.61 * ((150 / 63.151) ** 1.78963 - 1)), rel=1e-12)


def test_air_melting_pressure_above_critical_temperature():
    assert AIR.meltingpoint_pressure(150) == pytest.approx(5265 * (1 + 35493.5 * ((150 / 59.75) ** 1.78963 - 1)), rel=1e-12)


def test_high_pressure_states_inside_declared_range():
    assert NITROGEN.calculated_functions(300, 2.2e9).density_kg_m3 == pytest.approx(1308.1314, abs=0.01)
    assert NITROGEN.calculated_functions(290, 2.0e9).density_kg_m3 == pytest.approx(1285.0339, abs=0.01)
    assert AIR.calculated_functions(240, 2.0e9).density_kg_m3 / (AIR.constants.molar_mass_kg_mol * 1000) == pytest.approx(47.893, abs=1e-3)


def test_air_retrograde_region_is_two_phase():
    with pytest.raises(GeneralClass.TwoPhaseError):
        AIR.calculated_functions(132.6, 3.781e6)


@pytest.mark.parametrize("t", [132.55, 132.6, 132.63])
def test_air_saturation_up_to_maxcondentherm(t):
    saturation = AIR.saturation_state(t)
    assert saturation.vapor_pressure_Pa < saturation.liquid_pressure_Pa
    assert AIR.pressure(t, saturation.liquid_density_kg_m3) == pytest.approx(saturation.liquid_pressure_Pa, rel=1e-9)
    assert AIR.pressure(t, saturation.vapor_density_kg_m3) == pytest.approx(saturation.vapor_pressure_Pa, rel=1e-9)
    assert AIR.saturation_state(AIR.constants.saturation_maximum_temperature_K) is None


def test_phase_label_does_not_depend_on_input_mode():
    by_pressure = NITROGEN.calculated_functions(200, 5e6)
    by_density = NITROGEN.calculated_functions(200, density=by_pressure.density_kg_m3)
    assert by_pressure.phase == by_density.phase == "supercritical"


def test_derived_pressure_above_validity_is_refused():
    with pytest.raises(ValueError):
        NITROGEN.calculated_functions(1000, density=1400.0)


def test_temperatures_down_to_the_papers_lower_limit():
    assert NITROGEN.calculated_functions(63.151, 1e4).phase == "vapor"
    assert NITROGEN.calculated_functions(64, 1e5).phase == "liquid"
    assert AIR.calculated_functions(59.75, 2000).phase == "vapor"


def test_saturation_reported_with_state_comes_from_equation_of_state():
    result = NITROGEN.calculated_functions(77.355, 2e5)
    saturation = NITROGEN.saturation_state(77.355)
    assert result.bubble_point_pressure_Pa == saturation.liquid_pressure_Pa
    assert result.dew_point_density_kg_m3 == saturation.vapor_density_kg_m3
