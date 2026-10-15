import dataclasses
import typing


@dataclasses.dataclass(frozen=True)
class FluidConstants:
    molar_mass_kg_mol:                  typing.Optional[float]=None
    critical_temperature_K:             typing.Optional[float]=None
    critical_pressure_Pa:               typing.Optional[float]=None
    critical_density_kg_m3:             typing.Optional[float]=None
    reference_temperature_K:            typing.Optional[float]=None
    reference_pressure_Pa:              typing.Optional[float]=None
    reference_density_kg_m3:            typing.Optional[float]=None
    melting_temperature_K:              typing.Optional[float]=None
    melting_pressure_Pa:                typing.Optional[float]=None
    triplepoint_temperature_K:          typing.Optional[float]=None
    triplepoint_pressure_Pa:            typing.Optional[float]=None
    saturation_maximum_temperature_K:   typing.Optional[float]=None
    acentric_factor:                    typing.Optional[float]=None
    universal_gas_constant_J_mol_K:     float=8.314510

@dataclasses.dataclass(frozen=True)
class AcceptableParameters:
    T_maximum_value_K:typing.Optional[float]=None
    T_minimum_value_K:typing.Optional[float]=None
    P_maximum_value_Pa:typing.Optional[float]=None
    P_minimum_value_Pa:typing.Optional[float]=None


@dataclasses.dataclass(frozen=True)
class SaturationTerms:
    d_p_p:typing.Optional [tuple[float]]=None
    d_p_d:typing.Optional [tuple[float]]=None
    b_p_d:typing.Optional [tuple[float]]=None
    b_p_p:typing.Optional [tuple[float]]=None

@dataclasses.dataclass(frozen=True)
class IdealHelmholtzTerms:
    e:typing.Optional[float]=None

@dataclasses.dataclass(frozen=True)
class RealHelmholtzTerms:
    n:typing.Optional[float]=None
    i:typing.Optional[float]=None
    j:typing.Optional[float]=None
    l:typing.Optional[float]=None
    o:typing.Optional[float]=None
    k:typing.Optional[float]=None
    y:typing.Optional[float]=None


@dataclasses.dataclass(frozen=True)
class Fluid:
    name:                   typing.Optional[str]=None
    fluid_constants:        typing.Optional[FluidConstants]=None
    saturation_terms:       typing.Optional[SaturationTerms]=None
    real_helmholtz_terms:   typing.Optional[RealHelmholtzTerms]=None
    ideal_helmholtz_terms:  typing.Optional[IdealHelmholtzTerms]=None
    acceptable_parameters:  typing.Optional[AcceptableParameters]=None

@dataclasses.dataclass(frozen=True)
class HelmholtzFunctions:
    reduced_ideal_exp_helmholtz:                    typing.Optional[str]=None
    reduced_ideal_num_helmholtz:                    typing.Optional[callable]=None
    reduced_ideal_helmholtz:                        typing.Optional[float]=None

    reduced_real_exp_helmholtz:                     typing.Optional[str]=None
    reduced_real_num_helmholtz:                     typing.Optional[callable]=None
    reduced_real_helmholtz:                         typing.Optional[float]=None

    redtemp_pdev_redtemp_red_ideal_exp_helmholtz:   typing.Optional[float]=None
    redtemp_pdev_redtemp_red_ideal_num_helmholtz:   typing.Optional[callable]=None
    redtemp_pdev_redtemp_red_ideal_helmholtz:       typing.Optional[float]=None

    redtemp_pdev_redtemp_red_real_exp_helmholtz:    typing.Optional[float]=None
    redtemp_pdev_redtemp_red_real_num_helmholtz:    typing.Optional[callable]=None
    redtemp_pdev_redtemp_red_real_helmholtz:        typing.Optional[float]=None

    reddens_pdev_reddens_red_real_exp_helmholtz:    typing.Optional[float]=None
    reddens_pdev_reddens_red_real_num_helmholtz:    typing.Optional[callable]=None
    reddens_pdev_reddens_red_real_helmholtz:        typing.Optional[float]=None

    reddens_pdev2_reddens_red_real_exp_helmholtz:   typing.Optional[float]=None
    reddens_pdev2_reddens_red_real_num_helmholtz:   typing.Optional[callable]=None
    reddens_pdev2_reddens_red_real_helmholtz:       typing.Optional[float]=None

    reddens_pdev3_reddens_red_real_exp_helmholtz:   typing.Optional[float]=None
    reddens_pdev3_reddens_red_real_num_helmholtz:   typing.Optional[callable]=None
    reddens_pdev3_reddens_red_real_helmholtz:       typing.Optional[float]=None

    redtemp_pdev2_redtemp_red_ideal_exp_helmholtz:  typing.Optional[float]=None
    redtemp_pdev2_redtemp_red_ideal_num_helmholtz:  typing.Optional[callable]=None
    redtemp_pdev2_redtemp_red_ideal_helmholtz:      typing.Optional[float]=None

    redtemp_pdev2_redtemp_red_real_exp_helmholtz:   typing.Optional[float]=None
    redtemp_pdev2_redtemp_red_real_num_helmholtz:   typing.Optional[callable]=None
    redtemp_pdev2_redtemp_red_real_helmholtz:       typing.Optional[float]=None

    reddens_redtemp_pdev2_red_real_exp_helmholtz:   typing.Optional[float]=None
    reddens_redtemp_pdev2_red_real_num_helmholtz:   typing.Optional[callable]=None
    reddens_redtemp_pdev2_red_real_helmholtz:       typing.Optional[float]=None

@dataclasses.dataclass(frozen=True)
class StateFunctions:
    density_kg_m3:                                  typing.Optional[float]=None
    compressibility_factor:                         typing.Optional[float]=None
    internal_energy_J_kg:                           typing.Optional[float]=None
    enthalpy_J_kg:                                  typing.Optional[float]=None
    entropy_J_kg_K:                                 typing.Optional[float]=None
    gibbs_energy_J_kg:                              typing.Optional[float]=None
    helmholtz_energy_J_kg:                          typing.Optional[float]=None
    isochoric_heat_capacity_J_kg_K:                 typing.Optional[float]=None
    isobaric_heat_capacity_J_kg_K:                  typing.Optional[float]=None
    heat_capacity_ratio:                            typing.Optional[float]=None
    speed_of_sound_m_s:                             typing.Optional[float]=None
    pressure_density_derivative_Pa_m3_kg:           typing.Optional[float]=None
    pressure_density_second_derivative_Pa_m6_kg2:   typing.Optional[float]=None
    pressure_temperature_derivative_Pa_K:           typing.Optional[float]=None
    joule_thomson_coefficient_K_Pa:                 typing.Optional[float]=None
    isentropic_expansion_coefficient:               typing.Optional[float]=None
    isothermal_expansion_coefficient:               typing.Optional[float]=None
    volume_expansivity_1_K:                         typing.Optional[float]=None
    adiabatic_compressibility_1_Pa:                 typing.Optional[float]=None
    adiabatic_bulk_modulus_Pa:                      typing.Optional[float]=None
    isothermal_compressibility_1_Pa:                typing.Optional[float]=None
    isothermal_bulk_modulus_Pa:                     typing.Optional[float]=None
    second_virial_coefficient_m3_kg:                typing.Optional[float]=None
    third_virial_coefficient_m6_kg2:                typing.Optional[float]=None
    fugacity_coefficient:                           typing.Optional[float]=None


@dataclasses.dataclass(frozen=True)
class OutputFunctions:
    temperature_input_K:                            typing.Optional[float]=None
    pressure_input_Pa:                              typing.Optional[float]=None
    phase:                                          typing.Optional[str]=None
    density_kg_m3:                                  typing.Optional[StateFunctions]=None
    compressibility_factor:                         typing.Optional[StateFunctions]=None
    internal_energy_J_kg:                           typing.Optional[StateFunctions]=None
    enthalpy_J_kg:                                  typing.Optional[StateFunctions]=None
    entropy_J_kg_K:                                 typing.Optional[StateFunctions]=None
    gibbs_energy_J_kg:                              typing.Optional[float]=None
    helmholtz_energy_J_kg:                          typing.Optional[float]=None
    isochoric_heat_capacity_J_kg_K:                 typing.Optional[float]=None
    isobaric_heat_capacity_J_kg_K:                  typing.Optional[float]=None
    heat_capacity_ratio:                            typing.Optional[float]=None
    speed_of_sound_m_s:                             typing.Optional[float]=None
    pressure_density_derivative_Pa_m3_kg:           typing.Optional[float]=None
    pressure_density_second_derivative_Pa_m6_kg2:   typing.Optional[float]=None
    pressure_temperature_derivative_Pa_K:           typing.Optional[float]=None
    joule_thomson_coefficient_K_Pa:                 typing.Optional[float]=None
    isentropic_expansion_coefficient:               typing.Optional[float]=None
    isothermal_expansion_coefficient:               typing.Optional[float]=None
    volume_expansivity_1_K:                         typing.Optional[float]=None
    adiabatic_compressibility_1_Pa:                 typing.Optional[float]=None
    adiabatic_bulk_modulus_Pa:                      typing.Optional[float]=None
    isothermal_compressibility_1_Pa:                typing.Optional[float]=None
    isothermal_bulk_modulus_Pa:                     typing.Optional[float]=None
    second_virial_coefficient_m3_kg:                typing.Optional[float]=None
    third_virial_coefficient_m6_kg2:                typing.Optional[float]=None
    fugacity_coefficient:                           typing.Optional[float]=None
    dew_point_pressure_Pa:                          typing.Optional[float]=None
    dew_point_density_kg_m3:                        typing.Optional[float]=None
    bubble_point_pressure_Pa:                       typing.Optional[float]=None
    bubble_point_density_kg_m3:                     typing.Optional[float]=None
    melting_point_pressure_Pa:                      typing.Optional[float]=None


@dataclasses.dataclass(frozen=True)
class SaturationState:
    temperature_K:              typing.Optional[float]=None
    liquid_pressure_Pa:         typing.Optional[float]=None
    vapor_pressure_Pa:          typing.Optional[float]=None
    liquid_density_kg_m3:       typing.Optional[float]=None
    vapor_density_kg_m3:        typing.Optional[float]=None
