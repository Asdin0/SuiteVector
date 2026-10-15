import models
from abc import ABC, abstractmethod
import math as mth
import sympy as sy
import functools
import typing
import dataclasses


class TwoPhaseError(ValueError):
    pass


class SolidPhaseError(ValueError):
    pass


class RealFluid(ABC):
    two_phase_pressure_tolerance_Pa=1
    saturated_density_tolerance=1e-9
    pressure_tolerance=1e-9
    fugacity_coefficient_defined=False

    def __init__(self,
                 constants,
                 saturation_terms,
                 ideal_helmholtz_terms,
                 real_helmholtz_terms,
                 acceptable_parameters):
        self.constants=constants
        self.saturation_terms=saturation_terms
        self.ideal_helmholtz_terms=ideal_helmholtz_terms
        self.real_helmholtz_terms=real_helmholtz_terms
        self.acceptable_parameters=acceptable_parameters

    def _check_input_(self,t:typing.Optional[float]=None,p:typing.Optional[float]=None):
        if t!=None:
            if t<self.acceptable_parameters.T_minimum_value_K or t>self.acceptable_parameters.T_maximum_value_K:
                raise ValueError(f"The temperature must be between {self.acceptable_parameters.T_minimum_value_K }[k] and {self.acceptable_parameters.T_maximum_value_K }[[K]")
        if p!=None:
            if p<0 or p>self.acceptable_parameters.P_maximum_value_Pa:
                raise ValueError(f"The pressure must be between 0[Pa] and {self.acceptable_parameters.P_maximum_value_Pa }[[Pa]")

    def phases(self,t:float,p:typing.Optional[float]=None,density:typing.Optional[float]=None):
        self._check_input_(t,p)
        if p is None and density is None:
            return None
        saturation=self.saturation_state(t)

        if p is None:
            if saturation is None:
                side=None
            elif density>=saturation.liquid_density_kg_m3*(1-self.saturated_density_tolerance):
                side="liquid"
            elif density<=saturation.vapor_density_kg_m3*(1+self.saturated_density_tolerance):
                side="vapor"
            else:
                return "two-phase"
            pressure=self.pressure(t,density)
        else:
            pressure=p
            if saturation is None:
                side=None
            elif pressure>=saturation.liquid_pressure_Pa+self.two_phase_pressure_tolerance_Pa:
                side="liquid"
            elif pressure<=saturation.vapor_pressure_Pa-self.two_phase_pressure_tolerance_Pa:
                side="vapor"
            else:
                return "two-phase"

        if pressure>=self.meltingpoint_pressure(t):
            return "solid"
        if side is None or t>self.constants.critical_temperature_K:
            return "supercritical" if pressure>self.constants.critical_pressure_Pa else "gas"
        return side

    def _reduced_variables_(self,t:float,real_density:float):
        return self.constants.reference_temperature_K/t,real_density/self.constants.reference_density_kg_m3

    @functools.cache
    def _numerical_functions_(self):
        ideal_helmholtz_results=self._reduced_ideal_helmholtzenergy_()
        real_helmholtz_results=self._reduced_real_helmholtzenergy_()
        partial_derivative_results=self._partial_derivatives_(helmholtz_i_r_expression=ideal_helmholtz_results.reduced_ideal_exp_helmholtz,
                                                            helmholtz_r_r_expression=real_helmholtz_results.reduced_real_exp_helmholtz)
        return ideal_helmholtz_results,real_helmholtz_results,partial_derivative_results

    def pressure(self,t:float,real_density:float):
        partial_derivative_results=self._numerical_functions_()[2]
        reduced_temperature,reduced_density=self._reduced_variables_(t,real_density)
        delta_devpar_delta_helmholtz_r_r=reduced_density*partial_derivative_results.reddens_pdev_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        return float((real_density/self.constants.molar_mass_kg_mol)*self.constants.universal_gas_constant_J_mol_K*t*(1+delta_devpar_delta_helmholtz_r_r))

    def pressure_density_derivative(self,t:float,real_density:float):
        partial_derivative_results=self._numerical_functions_()[2]
        reduced_temperature,reduced_density=self._reduced_variables_(t,real_density)
        devpar_delta_helmholtz_r_r=partial_derivative_results.reddens_pdev_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        devpar2_delta_helmholtz_r_r=partial_derivative_results.reddens_pdev2_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        return float((self.constants.universal_gas_constant_J_mol_K*t/self.constants.molar_mass_kg_mol)*(1+2*reduced_density*devpar_delta_helmholtz_r_r+reduced_density**2*devpar2_delta_helmholtz_r_r))

    def _maxwell_saturation_(self,t:float):
        self._check_input_(t)
        if t>=self.constants.saturation_maximum_temperature_K:
            return None
        real_helmholtz_results,partial_derivative_results=self._numerical_functions_()[1:]
        helmholtz_r_r=real_helmholtz_results.reduced_real_num_helmholtz
        devpar_delta=partial_derivative_results.reddens_pdev_reddens_red_real_num_helmholtz
        devpar2_delta=partial_derivative_results.reddens_pdev2_reddens_red_real_num_helmholtz
        reduced_temperature=self.constants.reference_temperature_K/t

        def reduced_pressure(delta):
            return delta+delta**2*devpar_delta(reduced_temperature,delta)

        def reduced_gibbs(delta):
            return delta*devpar_delta(reduced_temperature,delta)+helmholtz_r_r(reduced_temperature,delta)+mth.log(delta)

        def reduced_pressure_derivative(delta):
            return 1+2*delta*devpar_delta(reduced_temperature,delta)+delta**2*devpar2_delta(reduced_temperature,delta)

        def reduced_gibbs_derivative(delta):
            return 2*devpar_delta(reduced_temperature,delta)+delta*devpar2_delta(reduced_temperature,delta)+1/delta

        delta_liquid=self.bubblepoint_density(t)/self.constants.reference_density_kg_m3
        delta_vapor=self.dewpoint_density(t)/self.constants.reference_density_kg_m3
        for _ in range(100):
            residual_pressure=reduced_pressure(delta_vapor)-reduced_pressure(delta_liquid)
            residual_gibbs=reduced_gibbs(delta_vapor)-reduced_gibbs(delta_liquid)
            if abs(residual_pressure)<=1e-10*reduced_pressure(delta_liquid) and abs(residual_gibbs)<=1e-10:
                break
            a=-reduced_pressure_derivative(delta_liquid)
            b=reduced_pressure_derivative(delta_vapor)
            c=-reduced_gibbs_derivative(delta_liquid)
            d=reduced_gibbs_derivative(delta_vapor)
            determinant=a*d-b*c
            step_liquid=(b*residual_gibbs-d*residual_pressure)/determinant
            step_vapor=(c*residual_pressure-a*residual_gibbs)/determinant
            damping=1
            while not 0<delta_vapor+damping*step_vapor<delta_liquid+damping*step_liquid:
                damping=damping/2
                if damping<1e-6:
                    raise ValueError(f"The saturation state at temperature {t} K did not converge")
            delta_liquid=delta_liquid+damping*step_liquid
            delta_vapor=delta_vapor+damping*step_vapor
        else:
            raise ValueError(f"The saturation state at temperature {t} K did not converge")

        liquid_density=delta_liquid*self.constants.reference_density_kg_m3
        vapor_density=delta_vapor*self.constants.reference_density_kg_m3
        if (liquid_density<=vapor_density*(1+1e-6)
                or self.pressure_density_derivative(t,liquid_density)<=0
                or self.pressure_density_derivative(t,vapor_density)<=0):
            raise ValueError(f"The saturation state at temperature {t} K did not converge")
        saturation_pressure=self.pressure(t,vapor_density)
        return models.SaturationState(temperature_K=t,
                                      liquid_pressure_Pa=saturation_pressure,
                                      vapor_pressure_Pa=saturation_pressure,
                                      liquid_density_kg_m3=liquid_density,
                                      vapor_density_kg_m3=vapor_density)

    def _ancillary_pressure_saturation_(self,t:float):
        self._check_input_(t)
        if t>=self.constants.saturation_maximum_temperature_K:
            return None
        liquid_pressure=self.bubblepoint_pressure(t)
        vapor_pressure=self.dewpoint_pressure(t)
        return models.SaturationState(temperature_K=t,
                                      liquid_pressure_Pa=liquid_pressure,
                                      vapor_pressure_Pa=vapor_pressure,
                                      liquid_density_kg_m3=self._branch_density_(t,liquid_pressure,self.bubblepoint_density(t)),
                                      vapor_density_kg_m3=self._branch_density_(t,vapor_pressure,self.dewpoint_density(t)))

    def _branch_density_(self,t:float,p:float,guess:float):
        low=high=guess
        for _ in range(5000):
            if self.pressure_density_derivative(t,low)<=0 or self.pressure_density_derivative(t,high)<=0:
                break
            if self.pressure(t,low)<=p<=self.pressure(t,high):
                return self._monotone_root_(t,p,low,high)
            if self.pressure(t,high)<p:
                low,high=high,high*1.001
            else:
                low,high=low/1.001,low
        raise ValueError(f"At temperature {t} K and pressure {p} Pa no stable density was found near {guess} kg/m3")

    def _monotone_root_(self,t:float,p:float,low:float,high:float):
        density=0.5*(low+high)
        for _ in range(300):
            pressure_error=self.pressure(t,density)-p
            if abs(pressure_error)<=1e-12*p:
                break
            if pressure_error>0:
                high=density
            else:
                low=density
            if high-low<=4e-16*high:
                break
            derivative=self.pressure_density_derivative(t,density)
            newton_density=density-pressure_error/derivative if derivative>0 else low
            density=newton_density if low<newton_density<high else 0.5*(low+high)
        if abs(self.pressure(t,density)-p)>self.pressure_tolerance*p or self.pressure_density_derivative(t,density)<=0:
            raise ValueError(f"At temperature {t} K and pressure {p} Pa no stable density was found")
        return density

    def _density_from_pressure_(self,t:float,p:float):
        saturation=self.saturation_state(t)
        minimum_density=1e-10*self.constants.reference_density_kg_m3
        if saturation is not None and p<saturation.vapor_pressure_Pa:
            return self._monotone_root_(t,p,minimum_density,saturation.vapor_density_kg_m3)
        low=saturation.liquid_density_kg_m3 if saturation is not None else minimum_density
        high=max(1.05*low,min(p*self.constants.molar_mass_kg_mol/(self.constants.universal_gas_constant_J_mol_K*t),2*self.constants.reference_density_kg_m3))
        for _ in range(200):
            if self.pressure_density_derivative(t,high)<=0:
                break
            if self.pressure(t,high)>=p:
                return self._monotone_root_(t,p,low,high)
            low,high=high,1.25*high
        raise ValueError(f"At temperature {t} K and pressure {p} Pa no stable density was found")

    def _helmholtz_values_(self,t:float,real_density:float):
        ideal_helmholtz_results,real_helmholtz_results,partial_derivative_results=self._numerical_functions_()
        reduced_temperature,reduced_density=self._reduced_variables_(t,real_density)

        helmholtz_i_r=ideal_helmholtz_results.reduced_ideal_num_helmholtz(reduced_temperature,reduced_density)
        helmholtz_r_r=real_helmholtz_results.reduced_real_num_helmholtz(reduced_temperature,reduced_density)

        delta_devpar_delta_helmholtz_r_r=reduced_density*partial_derivative_results.reddens_pdev_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        tao_devpar_tao_helmholtz_r_r=reduced_temperature*partial_derivative_results.redtemp_pdev_redtemp_red_real_num_helmholtz(reduced_temperature,reduced_density)
        tao_devpar_helmholtz_tao_i_r=reduced_temperature*partial_derivative_results.redtemp_pdev_redtemp_red_ideal_num_helmholtz(reduced_temperature,reduced_density)

        delta2_devpar2_delta_helmholtz_r_r=reduced_density**2*partial_derivative_results.reddens_pdev2_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        delta3_devpar3_delta_helmholtz_r_r=reduced_density**3*partial_derivative_results.reddens_pdev3_reddens_red_real_num_helmholtz(reduced_temperature,reduced_density)
        tao2_devpar2_tao_helmholtz_i_r=reduced_temperature**2*partial_derivative_results.redtemp_pdev2_redtemp_red_ideal_num_helmholtz(reduced_temperature,reduced_density)
        tao2_devpar2_tao_helmholtz_r_r=reduced_temperature**2*partial_derivative_results.redtemp_pdev2_redtemp_red_real_num_helmholtz(reduced_temperature,reduced_density)
        delta_tao_devpar2_delta_tao_helmholtz_r_r=reduced_density*reduced_temperature*partial_derivative_results.reddens_redtemp_pdev2_red_real_num_helmholtz(reduced_temperature,reduced_density)

        return models.HelmholtzFunctions(reduced_ideal_helmholtz=helmholtz_i_r,reduced_real_helmholtz=helmholtz_r_r,redtemp_pdev_redtemp_red_ideal_helmholtz=tao_devpar_helmholtz_tao_i_r,redtemp_pdev_redtemp_red_real_helmholtz=tao_devpar_tao_helmholtz_r_r,reddens_pdev_reddens_red_real_helmholtz=delta_devpar_delta_helmholtz_r_r,
                                         reddens_pdev2_reddens_red_real_helmholtz=delta2_devpar2_delta_helmholtz_r_r,reddens_pdev3_reddens_red_real_helmholtz=delta3_devpar3_delta_helmholtz_r_r,redtemp_pdev2_redtemp_red_ideal_helmholtz=tao2_devpar2_tao_helmholtz_i_r,redtemp_pdev2_redtemp_red_real_helmholtz=tao2_devpar2_tao_helmholtz_r_r,reddens_redtemp_pdev2_red_real_helmholtz=delta_tao_devpar2_delta_tao_helmholtz_r_r)

    def _state_function_(self,t:float,real_density:float):
        helmholtz_values=self._helmholtz_values_(t,real_density)
        helmholtz_i_r=helmholtz_values.reduced_ideal_helmholtz
        helmholtz_r_r=helmholtz_values.reduced_real_helmholtz
        tao_devpar_helmholtz_tao_i_r=helmholtz_values.redtemp_pdev_redtemp_red_ideal_helmholtz
        tao_devpar_tao_helmholtz_r_r=helmholtz_values.redtemp_pdev_redtemp_red_real_helmholtz
        delta_devpar_delta_helmholtz_r_r=helmholtz_values.reddens_pdev_reddens_red_real_helmholtz
        delta2_devpar2_delta_helmholtz_r_r=helmholtz_values.reddens_pdev2_reddens_red_real_helmholtz
        delta3_devpar3_delta_helmholtz_r_r=helmholtz_values.reddens_pdev3_reddens_red_real_helmholtz
        tao2_devpar2_tao_helmholtz_i_r=helmholtz_values.redtemp_pdev2_redtemp_red_ideal_helmholtz
        tao2_devpar2_tao_helmholtz_r_r=helmholtz_values.redtemp_pdev2_redtemp_red_real_helmholtz
        delta_tao_devpar2_delta_tao_helmholtz_r_r=helmholtz_values.reddens_redtemp_pdev2_red_real_helmholtz

        compressibility_factor=float(1+delta_devpar_delta_helmholtz_r_r)
        specifica_internal_energy=float(((tao_devpar_helmholtz_tao_i_r+tao_devpar_tao_helmholtz_r_r)*self.constants.universal_gas_constant_J_mol_K*t)/self.constants.molar_mass_kg_mol)
        specific_enthalpy=float((((tao_devpar_helmholtz_tao_i_r+tao_devpar_tao_helmholtz_r_r+delta_devpar_delta_helmholtz_r_r+1)*self.constants.universal_gas_constant_J_mol_K*t)/self.constants.molar_mass_kg_mol))
        specific_entropy=float(((tao_devpar_helmholtz_tao_i_r+tao_devpar_tao_helmholtz_r_r-helmholtz_i_r-helmholtz_r_r)*self.constants.universal_gas_constant_J_mol_K)/self.constants.molar_mass_kg_mol)

        specific_gas_constant=self.constants.universal_gas_constant_J_mol_K/self.constants.molar_mass_kg_mol
        density=float(real_density)
        pressure=density*specific_gas_constant*t*compressibility_factor
        density_term=float(1+2*delta_devpar_delta_helmholtz_r_r+delta2_devpar2_delta_helmholtz_r_r)
        temperature_term=float(1+delta_devpar_delta_helmholtz_r_r-delta_tao_devpar2_delta_tao_helmholtz_r_r)

        gibbs_energy=float((1+helmholtz_i_r+helmholtz_r_r+delta_devpar_delta_helmholtz_r_r)*specific_gas_constant*t)
        helmholtz_energy=float((helmholtz_i_r+helmholtz_r_r)*specific_gas_constant*t)
        isochoric_heat_capacity=float(-(tao2_devpar2_tao_helmholtz_i_r+tao2_devpar2_tao_helmholtz_r_r)*specific_gas_constant)
        pressure_density_derivative=float(specific_gas_constant*t*density_term)
        pressure_density_second_derivative=float((specific_gas_constant*t/density)*(2*delta_devpar_delta_helmholtz_r_r+4*delta2_devpar2_delta_helmholtz_r_r+delta3_devpar3_delta_helmholtz_r_r))
        pressure_temperature_derivative=float(density*specific_gas_constant*temperature_term)
        second_virial_coefficient,third_virial_coefficient=self._virial_coefficients_(t)

        isobaric_heat_capacity=heat_capacity_ratio=speed_of_sound=joule_thomson_coefficient=volume_expansivity=None
        isentropic_expansion_coefficient=adiabatic_compressibility=adiabatic_bulk_modulus=None
        isothermal_expansion_coefficient=isothermal_compressibility=isothermal_bulk_modulus=None
        fugacity_coefficient=None
        if density_term>0:
            volume_expansivity=(1/density)*pressure_temperature_derivative*(1/pressure_density_derivative)
            if pressure>0:
                isothermal_expansion_coefficient=(density/pressure)*pressure_density_derivative
                isothermal_compressibility=1/(isothermal_expansion_coefficient*pressure)
                isothermal_bulk_modulus=isothermal_expansion_coefficient*pressure
            if isochoric_heat_capacity>0:
                isobaric_heat_capacity=isochoric_heat_capacity+specific_gas_constant*temperature_term**2/density_term
                heat_capacity_ratio=isobaric_heat_capacity/isochoric_heat_capacity
                speed_of_sound_squared=specific_gas_constant*t*heat_capacity_ratio*density_term
                speed_of_sound=mth.sqrt(speed_of_sound_squared)
                joule_thomson_coefficient=(t*volume_expansivity-1)/(density*isobaric_heat_capacity)
                if pressure>0:
                    isentropic_expansion_coefficient=speed_of_sound_squared*density/pressure
                    adiabatic_compressibility=1/(isentropic_expansion_coefficient*pressure)
                    adiabatic_bulk_modulus=isentropic_expansion_coefficient*pressure
        if self.fugacity_coefficient_defined and compressibility_factor>0:
            fugacity_coefficient=mth.exp(compressibility_factor-1-mth.log(compressibility_factor)+helmholtz_r_r)

        return models.StateFunctions(density_kg_m3=real_density,compressibility_factor=compressibility_factor,internal_energy_J_kg=specifica_internal_energy,enthalpy_J_kg=specific_enthalpy,entropy_J_kg_K=specific_entropy,
                                     gibbs_energy_J_kg=gibbs_energy,
                                     helmholtz_energy_J_kg=helmholtz_energy,
                                     isochoric_heat_capacity_J_kg_K=isochoric_heat_capacity,
                                     isobaric_heat_capacity_J_kg_K=isobaric_heat_capacity,
                                     heat_capacity_ratio=heat_capacity_ratio,
                                     speed_of_sound_m_s=speed_of_sound,
                                     pressure_density_derivative_Pa_m3_kg=pressure_density_derivative,
                                     pressure_density_second_derivative_Pa_m6_kg2=pressure_density_second_derivative,
                                     pressure_temperature_derivative_Pa_K=pressure_temperature_derivative,
                                     joule_thomson_coefficient_K_Pa=joule_thomson_coefficient,
                                     isentropic_expansion_coefficient=isentropic_expansion_coefficient,
                                     isothermal_expansion_coefficient=isothermal_expansion_coefficient,
                                     volume_expansivity_1_K=volume_expansivity,
                                     adiabatic_compressibility_1_Pa=adiabatic_compressibility,
                                     adiabatic_bulk_modulus_Pa=adiabatic_bulk_modulus,
                                     isothermal_compressibility_1_Pa=isothermal_compressibility,
                                     isothermal_bulk_modulus_Pa=isothermal_bulk_modulus,
                                     second_virial_coefficient_m3_kg=second_virial_coefficient,
                                     third_virial_coefficient_m6_kg2=third_virial_coefficient,
                                     fugacity_coefficient=fugacity_coefficient)

    def _virial_coefficients_(self,t:float):
        partial_derivative_results=self._numerical_functions_()[2]
        reduced_temperature=self.constants.reference_temperature_K/t
        second_virial_coefficient=float(partial_derivative_results.reddens_pdev_reddens_red_real_num_helmholtz(reduced_temperature,0.0)/self.constants.reference_density_kg_m3)
        third_virial_coefficient=float(partial_derivative_results.reddens_pdev2_reddens_red_real_num_helmholtz(reduced_temperature,0.0)/self.constants.reference_density_kg_m3**2)
        return second_virial_coefficient,third_virial_coefficient

    @functools.cache
    def _partial_derivatives_(self,helmholtz_i_r_expression,helmholtz_r_r_expression):
        t_r,d_r=sy.symbols('t_r d_r', real=True)
        devpar_helmholtz_tao_i_r_symbol=sy.diff(helmholtz_i_r_expression,t_r)
        devpar_helmholtz_tao_i_r_numerical=sy.lambdify((t_r,d_r),devpar_helmholtz_tao_i_r_symbol)

        devpar_delta_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,d_r)
        devpar_delta_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar_delta_helmholtz_r_r_symbol)

        devpar2_delta_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,d_r,2)
        devpar2_delta_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar2_delta_helmholtz_r_r_symbol)

        devpar_tao_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,t_r)
        devpar_tao_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar_tao_helmholtz_r_r_symbol)

        devpar3_delta_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,d_r,3)
        devpar3_delta_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar3_delta_helmholtz_r_r_symbol)

        devpar2_tao_helmholtz_i_r_symbol=sy.diff(helmholtz_i_r_expression,t_r,2)
        devpar2_tao_helmholtz_i_r_numerical=sy.lambdify((t_r,d_r),devpar2_tao_helmholtz_i_r_symbol)

        devpar2_tao_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,t_r,2)
        devpar2_tao_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar2_tao_helmholtz_r_r_symbol)

        devpar2_delta_tao_helmholtz_r_r_symbol=sy.diff(helmholtz_r_r_expression,d_r,t_r)
        devpar2_delta_tao_helmholtz_r_r_numerical=sy.lambdify((t_r,d_r),devpar2_delta_tao_helmholtz_r_r_symbol)

        return models.HelmholtzFunctions(redtemp_pdev_redtemp_red_ideal_exp_helmholtz  =   devpar_helmholtz_tao_i_r_symbol,
                                redtemp_pdev_redtemp_red_ideal_num_helmholtz  =   devpar_helmholtz_tao_i_r_numerical,
                                redtemp_pdev_redtemp_red_real_exp_helmholtz   =   devpar_tao_helmholtz_r_r_symbol,
                                redtemp_pdev_redtemp_red_real_num_helmholtz   =   devpar_tao_helmholtz_r_r_numerical,
                                reddens_pdev_reddens_red_real_exp_helmholtz   =   devpar_delta_helmholtz_r_r_symbol,
                                reddens_pdev_reddens_red_real_num_helmholtz   =   devpar_delta_helmholtz_r_r_numerical,
                                reddens_pdev2_reddens_red_real_exp_helmholtz  =   devpar2_delta_helmholtz_r_r_symbol,
                                reddens_pdev2_reddens_red_real_num_helmholtz  =   devpar2_delta_helmholtz_r_r_numerical,
                                reddens_pdev3_reddens_red_real_exp_helmholtz  =   devpar3_delta_helmholtz_r_r_symbol,
                                reddens_pdev3_reddens_red_real_num_helmholtz  =   devpar3_delta_helmholtz_r_r_numerical,
                                redtemp_pdev2_redtemp_red_ideal_exp_helmholtz =   devpar2_tao_helmholtz_i_r_symbol,
                                redtemp_pdev2_redtemp_red_ideal_num_helmholtz =   devpar2_tao_helmholtz_i_r_numerical,
                                redtemp_pdev2_redtemp_red_real_exp_helmholtz  =   devpar2_tao_helmholtz_r_r_symbol,
                                redtemp_pdev2_redtemp_red_real_num_helmholtz  =   devpar2_tao_helmholtz_r_r_numerical,
                                reddens_redtemp_pdev2_red_real_exp_helmholtz  =   devpar2_delta_tao_helmholtz_r_r_symbol,
                                reddens_redtemp_pdev2_red_real_num_helmholtz  =   devpar2_delta_tao_helmholtz_r_r_numerical)

    def calculated_functions(self,t:float,p:typing.Optional[float]=None,density:typing.Optional[float]=None):
        phase=self.phases(t,p,density)
        if phase=="two-phase":
            raise TwoPhaseError(f"At temperature {t} K the fluid is in a two-phase condition")
        if phase=="solid":
            raise SolidPhaseError(f"At temperature {t} K the fluid is solid above the melting pressure {self.meltingpoint_pressure(t)} Pa")

        if density is None:
            real_density=self._density_from_pressure_(t,p)
        else:
            self._check_input_(t,self.pressure(t,density))
            real_density=density
        state_function_results=self._state_function_(t,real_density)
        saturation=self.saturation_state(t)

        return models.OutputFunctions(temperature_input_K=t,
                            pressure_input_Pa=p,
                            phase=phase,
                            **dataclasses.asdict(state_function_results),
                            dew_point_density_kg_m3=saturation.vapor_density_kg_m3 if saturation else None,
                            dew_point_pressure_Pa=saturation.vapor_pressure_Pa if saturation else None,
                            bubble_point_density_kg_m3=saturation.liquid_density_kg_m3 if saturation else None,
                            bubble_point_pressure_Pa=saturation.liquid_pressure_Pa if saturation else None,
                            melting_point_pressure_Pa=self.meltingpoint_pressure(t))

    @abstractmethod
    def dewpoint_pressure(self,t:float):
        pass

    @abstractmethod
    def dewpoint_density(self,t:float):
        pass

    @abstractmethod
    def bubblepoint_pressure(self,t:float):
        pass

    @abstractmethod
    def bubblepoint_density(self,t:float):
        pass

    @abstractmethod
    def meltingpoint_pressure(self,t:float):
        pass

    @abstractmethod
    def saturation_state(self,t:float):
        pass

    @abstractmethod
    def _reduced_ideal_helmholtzenergy_(self):
        pass

    @abstractmethod
    def _reduced_real_helmholtzenergy_(self):
        pass
