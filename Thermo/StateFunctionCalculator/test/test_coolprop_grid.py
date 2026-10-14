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

CoolProp = pytest.importorskip("CoolProp.CoolProp")

NITROGEN = nitrogen.Nitrogen()
AIR = air.Air()
AIR_MOLAR_MASS = AIR.constants.molar_mass_kg_mol
AIR_STATE = CoolProp.AbstractState("HEOS", "Air")

nitrogen_states = [(t, p) for t in (65, 70, 77, 85, 95, 105, 112, 118, 122, 124, 125, 125.9, 126.19, 126.3, 127, 130, 140, 160, 200, 300, 500, 1000)
                   for p in (1e3, 1e5, 1e6, 3e6, 3.35e6, 3.4e6, 5e6, 1e7, 1e8, 1e9, 2.2e9)]
air_states = [(t, p) for t in (60, 65, 75, 85, 100, 115, 125, 130, 132, 132.4, 132.6, 132.7, 133, 135, 140, 160, 200, 240, 300, 600, 1000, 2000)
              for p in (1e3, 1e5, 1e6, 2e6, 3e6, 3.6e6, 3.8e6, 4e6, 1e7, 1e8, 1e9, 2e9)]


def air_coolprop(t, p, phase):
    AIR_STATE.specify_phase(CoolProp.iphase_liquid if phase == "liquid" else CoolProp.iphase_gas)
    try:
        AIR_STATE.update(CoolProp.PT_INPUTS, p, t)
        return AIR_STATE.rhomolar(), AIR_STATE.hmolar(), AIR_STATE.smolar()
    finally:
        AIR_STATE.unspecify_phase()


@pytest.mark.parametrize(["t", "p"], nitrogen_states)
def test_nitrogen_against_coolprop(t, p):
    try:
        result = NITROGEN.calculated_functions(t, p)
    except (GeneralClass.TwoPhaseError, GeneralClass.SolidPhaseError):
        pytest.skip("outside the single-phase fluid region")
    near_critical = abs(t - NITROGEN.constants.critical_temperature_K) <= 2
    assert result.density_kg_m3 == pytest.approx(CoolProp.PropsSI('D', 'T', t, 'P', p, 'Nitrogen'), rel=1e-5 if near_critical else 1e-6)
    assert result.enthalpy_J_kg == pytest.approx(CoolProp.PropsSI('H', 'T', t, 'P', p, 'Nitrogen'), abs=0.5)
    assert result.entropy_J_kg_K == pytest.approx(CoolProp.PropsSI('S', 'T', t, 'P', p, 'Nitrogen'), abs=5e-3)


@pytest.mark.parametrize(["t", "p"], air_states)
def test_air_against_coolprop(t, p):
    try:
        result = AIR.calculated_functions(t, p)
    except (GeneralClass.TwoPhaseError, GeneralClass.SolidPhaseError):
        pytest.skip("outside the single-phase fluid region")
    reference = AIR.calculated_functions(300, 101325)
    reference_density, reference_enthalpy, reference_entropy = air_coolprop(300, 101325, "gas")
    density, enthalpy, entropy = air_coolprop(t, p, result.phase)
    assert result.density_kg_m3 / AIR_MOLAR_MASS == pytest.approx(density, rel=1e-9)
    assert (result.enthalpy_J_kg - reference.enthalpy_J_kg) * AIR_MOLAR_MASS == pytest.approx(enthalpy - reference_enthalpy, abs=1e-3)
    assert (result.entropy_J_kg_K - reference.entropy_J_kg_K) * AIR_MOLAR_MASS == pytest.approx(entropy - reference_entropy, abs=1e-6)


@pytest.mark.parametrize("t", [63.151, 65, 70, 77.355, 90, 100, 110, 120, 125, 126, 126.15, 126.19, 126.1919])
def test_nitrogen_saturation_against_coolprop(t):
    saturation = NITROGEN.saturation_state(t)
    assert saturation.vapor_pressure_Pa == pytest.approx(CoolProp.PropsSI('P', 'T', t, 'Q', 0, 'Nitrogen'), rel=5e-7)
    assert saturation.liquid_density_kg_m3 == pytest.approx(CoolProp.PropsSI('D', 'T', t, 'Q', 0, 'Nitrogen'), rel=5e-7)
    assert saturation.vapor_density_kg_m3 == pytest.approx(CoolProp.PropsSI('D', 'T', t, 'Q', 1, 'Nitrogen'), rel=5e-7)


@pytest.mark.parametrize("t", [60, 65, 70, 80, 90, 100, 110, 120, 125, 130, 132, 132.5])
def test_air_saturation_densities_against_coolprop(t):
    saturation = AIR.saturation_state(t)
    liquid_density = air_coolprop(t, saturation.liquid_pressure_Pa, "liquid")[0]
    vapor_density = air_coolprop(t, saturation.vapor_pressure_Pa, "vapor")[0]
    assert saturation.liquid_density_kg_m3 / AIR_MOLAR_MASS == pytest.approx(liquid_density, rel=1e-9)
    assert saturation.vapor_density_kg_m3 / AIR_MOLAR_MASS == pytest.approx(vapor_density, rel=1e-9)
