import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
sys.path.insert(0, os.path.join(project_root, 'general_class'))
sys.path.insert(0, os.path.join(project_root, 'fluid_list'))

import math as mth
import pytest
import GeneralClass
import nitrogen
import air

NITROGEN = nitrogen.Nitrogen()
AIR = air.Air()
FLUIDS = {"nitrogen": NITROGEN, "air": AIR}

NEW_PROPERTIES = ("gibbs_energy_J_kg", "helmholtz_energy_J_kg", "isochoric_heat_capacity_J_kg_K", "isobaric_heat_capacity_J_kg_K",
                  "heat_capacity_ratio", "speed_of_sound_m_s", "pressure_density_derivative_Pa_m3_kg",
                  "pressure_density_second_derivative_Pa_m6_kg2", "pressure_temperature_derivative_Pa_K", "joule_thomson_coefficient_K_Pa",
                  "isentropic_expansion_coefficient", "isothermal_expansion_coefficient", "volume_expansivity_1_K",
                  "adiabatic_compressibility_1_Pa", "adiabatic_bulk_modulus_Pa", "isothermal_compressibility_1_Pa", "isothermal_bulk_modulus_Pa",
                  "second_virial_coefficient_m3_kg", "third_virial_coefficient_m6_kg2", "fugacity_coefficient")

states = [("nitrogen", 300, 0.1e6), ("nitrogen", 80, 10e6), ("nitrogen", 150, 10e6), ("nitrogen", 500, 1000e6),
          ("air", 300, 0.101325e6), ("air", 100, 1e6), ("air", 150, 5e6), ("air", 1000, 10e6)]

nitrogen_states = [(t, p) for t in (65, 70, 77, 85, 95, 105, 112, 118, 122, 124, 125, 125.9, 126.19, 126.3, 127, 130, 140, 160, 200, 300, 500, 1000)
                   for p in (1e3, 1e5, 1e6, 3e6, 3.35e6, 3.4e6, 5e6, 1e7, 1e8, 1e9, 2.2e9)]
air_states = [(t, p) for t in (60, 65, 75, 85, 100, 115, 125, 130, 132, 132.4, 132.6, 132.7, 133, 135, 140, 160, 200, 240, 300, 600, 1000, 2000)
              for p in (1e3, 1e5, 1e6, 2e6, 3e6, 3.6e6, 3.8e6, 4e6, 1e7, 1e8, 1e9, 2e9)]


@pytest.mark.parametrize(["fluid_key", "t", "p"], states)
def test_relations_between_printed_equations(fluid_key, t, p):
    fluid = FLUIDS[fluid_key]
    result = fluid.calculated_functions(t, p)
    density = result.density_kg_m3
    pressure = fluid.pressure(t, density)
    assert result.gibbs_energy_J_kg == pytest.approx(result.enthalpy_J_kg - t * result.entropy_J_kg_K, rel=1e-10)
    assert result.helmholtz_energy_J_kg == pytest.approx(result.internal_energy_J_kg - t * result.entropy_J_kg_K, rel=1e-10)
    assert result.heat_capacity_ratio == pytest.approx(result.isobaric_heat_capacity_J_kg_K / result.isochoric_heat_capacity_J_kg_K, rel=1e-12)
    assert result.adiabatic_bulk_modulus_Pa == pytest.approx(result.isentropic_expansion_coefficient * pressure, rel=1e-10)
    assert result.adiabatic_bulk_modulus_Pa == pytest.approx(1 / result.adiabatic_compressibility_1_Pa, rel=1e-10)
    assert result.isothermal_bulk_modulus_Pa == pytest.approx(result.isothermal_expansion_coefficient * pressure, rel=1e-10)
    assert result.isothermal_bulk_modulus_Pa == pytest.approx(1 / result.isothermal_compressibility_1_Pa, rel=1e-10)
    assert result.pressure_density_derivative_Pa_m3_kg == pytest.approx(fluid.pressure_density_derivative(t, density), rel=1e-12)

    values = fluid._helmholtz_values_(t, density)
    density_term = 1 + 2 * values.reddens_pdev_reddens_red_real_helmholtz + values.reddens_pdev2_reddens_red_real_helmholtz
    temperature_term = 1 + values.reddens_pdev_reddens_red_real_helmholtz - values.reddens_redtemp_pdev2_red_real_helmholtz
    second_temperature_derivatives = values.redtemp_pdev2_redtemp_red_ideal_helmholtz + values.redtemp_pdev2_redtemp_red_real_helmholtz
    specific_gas_constant = fluid.constants.universal_gas_constant_J_mol_K / fluid.constants.molar_mass_kg_mol
    speed_of_sound_squared = specific_gas_constant * t * (density_term - temperature_term ** 2 / second_temperature_derivatives)
    assert result.speed_of_sound_m_s ** 2 == pytest.approx(speed_of_sound_squared, rel=1e-10)


@pytest.mark.parametrize(["fluid_key", "t", "p"], states)
def test_derivatives_against_finite_differences(fluid_key, t, p):
    fluid = FLUIDS[fluid_key]
    result = fluid.calculated_functions(t, p)
    density = result.density_kg_m3

    def specific_volume(t_value, p_value):
        return 1 / fluid.calculated_functions(t_value, p_value).density_kg_m3

    dt = 1e-4 * t
    pressure_temperature_derivative = (fluid.pressure(t + dt, density) - fluid.pressure(t - dt, density)) / (2 * dt)
    assert result.pressure_temperature_derivative_Pa_K == pytest.approx(pressure_temperature_derivative, rel=1e-5)

    ddensity = 1e-4 * density
    pressure_density_second_derivative = (fluid.pressure_density_derivative(t, density + ddensity)
                                          - fluid.pressure_density_derivative(t, density - ddensity)) / (2 * ddensity)
    assert result.pressure_density_second_derivative_Pa_m6_kg2 == pytest.approx(pressure_density_second_derivative, rel=1e-5)

    volume = 1 / density
    volume_expansivity = (specific_volume(t + dt, p) - specific_volume(t - dt, p)) / (2 * dt) / volume
    assert result.volume_expansivity_1_K == pytest.approx(volume_expansivity, rel=1e-5)

    dp = 1e-4 * p
    isothermal_compressibility = -(specific_volume(t, p + dp) - specific_volume(t, p - dp)) / (2 * dp) / volume
    assert result.isothermal_compressibility_1_Pa == pytest.approx(isothermal_compressibility, rel=1e-5)

    enthalpy = result.enthalpy_J_kg

    def isenthalpic_temperature(p_value):
        previous_t, current_t = t, t * (1 + 1e-6)
        previous_h = fluid.calculated_functions(previous_t, p_value).enthalpy_J_kg - enthalpy
        current_h = fluid.calculated_functions(current_t, p_value).enthalpy_J_kg - enthalpy
        for _ in range(50):
            if current_h == previous_h:
                break
            next_t = current_t - current_h * (current_t - previous_t) / (current_h - previous_h)
            previous_t, previous_h = current_t, current_h
            current_t = next_t
            current_h = fluid.calculated_functions(current_t, p_value).enthalpy_J_kg - enthalpy
            if abs(current_t - previous_t) <= 1e-13 * t:
                break
        return current_t

    dp = 1e-3 * p
    joule_thomson_coefficient = (isenthalpic_temperature(p + dp) - isenthalpic_temperature(p - dp)) / (2 * dp)
    assert result.joule_thomson_coefficient_K_Pa == pytest.approx(joule_thomson_coefficient, rel=1e-4)


@pytest.mark.parametrize("fluid_key", ["nitrogen", "air"])
@pytest.mark.parametrize("t", [100, 200, 500, 1000])
def test_virial_coefficients_are_intercept_and_slope_at_zero_density(fluid_key, t):
    fluid = FLUIDS[fluid_key]
    reference_density = fluid.constants.reference_density_kg_m3

    def z_minus_one_over_density(density):
        return (fluid.calculated_functions(t, density=density).compressibility_factor - 1) / density

    result = fluid.calculated_functions(t, density=1e-6 * reference_density)
    second_virial_coefficient = result.second_virial_coefficient_m3_kg
    assert second_virial_coefficient == pytest.approx(z_minus_one_over_density(1e-6 * reference_density), rel=1e-5)

    density = 1e-5 * reference_density
    third_virial_coefficient = (z_minus_one_over_density(density) - second_virial_coefficient) / density
    assert result.third_virial_coefficient_m6_kg2 == pytest.approx(third_virial_coefficient, rel=1e-3)


@pytest.mark.parametrize("t", [65, 77.355, 100, 120, 125])
def test_nitrogen_saturation_has_equal_gibbs_energies_and_fugacity_coefficients(t):
    saturation = NITROGEN.saturation_state(t)
    liquid = NITROGEN.calculated_functions(t, density=saturation.liquid_density_kg_m3)
    vapor = NITROGEN.calculated_functions(t, density=saturation.vapor_density_kg_m3)
    assert saturation.liquid_pressure_Pa == saturation.vapor_pressure_Pa
    assert liquid.gibbs_energy_J_kg == pytest.approx(vapor.gibbs_energy_J_kg, rel=1e-9)
    assert liquid.fugacity_coefficient == pytest.approx(vapor.fugacity_coefficient, rel=1e-9)
    assert all(type(getattr(side, name)) is float for side in (liquid, vapor) for name in NEW_PROPERTIES)


def test_fugacity_coefficient_only_where_the_paper_defines_it():
    assert NITROGEN.calculated_functions(300, 1e5).fugacity_coefficient is not None
    assert AIR.calculated_functions(300, 101325).fugacity_coefficient is None


@pytest.mark.parametrize(["fluid_key", "t", "p"], [("nitrogen", t, p) for t, p in nitrogen_states] + [("air", t, p) for t, p in air_states])
def test_new_properties_are_finite_in_the_single_phase_region(fluid_key, t, p):
    fluid = FLUIDS[fluid_key]
    try:
        result = fluid.calculated_functions(t, p)
    except (GeneralClass.TwoPhaseError, GeneralClass.SolidPhaseError):
        pytest.skip("outside the single-phase fluid region")
    for name in NEW_PROPERTIES:
        value = getattr(result, name)
        if name == "fugacity_coefficient" and not fluid.fugacity_coefficient_defined:
            assert value is None
        else:
            assert value is not None and mth.isfinite(value), name
