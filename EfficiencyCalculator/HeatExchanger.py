import math as mth
import typing


def _logarithmic_mean_(difference_1:float,difference_2:float):
    if difference_1<=0 or difference_2<=0:
        raise ValueError("The terminal temperature differences must be positive: the hot and cold streams cross")
    if mth.isclose(difference_1,difference_2,rel_tol=1e-9):
        return difference_1
    return (difference_1-difference_2)/mth.log(difference_1/difference_2)


def _capacity_rates_(m_h:float,m_c:float,cp_h:float,cp_c:float):
    c_h=m_h*cp_h
    c_c=m_c*cp_c
    if c_h<=0 or c_c<=0:
        raise ValueError("Mass flow rates and specific heats must be positive")
    return c_h,c_c


def _epsilon_(in_t_h,out_t_h,in_t_c,out_t_c,m_h,m_c,cp_h,cp_c):
    if None in (m_h,m_c,cp_h,cp_c,in_t_h,in_t_c) or (out_t_h is None and out_t_c is None):
        return None
    c_h,c_c=_capacity_rates_(m_h,m_c,cp_h,cp_c)
    if in_t_h<=in_t_c:
        raise ValueError("The hot inlet temperature must be higher than the cold inlet temperature")
    qmax=min(c_h,c_c)*(in_t_h-in_t_c)
    if out_t_h is not None:
        q=c_h*(in_t_h-out_t_h)
    else:
        q=c_c*(out_t_c-in_t_c)
    epsilon=q/qmax
    if not 0<=epsilon<=1:
        raise ValueError(f"The effectiveness {epsilon} is outside [0, 1]: the temperatures are not consistent with the energy balance")
    return epsilon


class EqHeatExchanger():
    @staticmethod
    def logarithmic_mean_temperature_difference(in_t_h:typing.Optional[float]=None,out_t_h:typing.Optional[float]=None,in_t_c:typing.Optional[float]=None,out_t_c:typing.Optional[float]=None):
        if None in (in_t_h,out_t_h,in_t_c,out_t_c):
            return None
        return _logarithmic_mean_(in_t_h-in_t_c,out_t_h-out_t_c)

    @staticmethod
    def epsilon(in_t_h:typing.Optional[float]=None,out_t_h:typing.Optional[float]=None,in_t_c:typing.Optional[float]=None,out_t_c:typing.Optional[float]=None,m_h:typing.Optional[float]=None,m_c:typing.Optional[float]=None,cp_h:typing.Optional[float]=None,cp_c:typing.Optional[float]=None):
        return _epsilon_(in_t_h,out_t_h,in_t_c,out_t_c,m_h,m_c,cp_h,cp_c)

    @staticmethod
    def ntu(epsilon:typing.Optional[float]=None,m_h:typing.Optional[float]=None,m_c:typing.Optional[float]=None,cp_h:typing.Optional[float]=None,cp_c:typing.Optional[float]=None):
        if None in (epsilon,m_h,m_c,cp_h,cp_c):
            return None
        c_h,c_c=_capacity_rates_(m_h,m_c,cp_h,cp_c)
        c=min(c_h,c_c)/max(c_h,c_c)
        if not 0<=epsilon<1/(1+c):
            raise ValueError(f"In parallel flow the effectiveness must be lower than 1/(1+Cr) = {1/(1+c)}")
        return -(mth.log(1-epsilon*(1+c)))/(1+c)


class CfHeatExchanger():
    @staticmethod
    def logarithmic_mean_temperature_difference(in_t_h:typing.Optional[float]=None,out_t_h:typing.Optional[float]=None,in_t_c:typing.Optional[float]=None,out_t_c:typing.Optional[float]=None):
        if None in (in_t_h,out_t_h,in_t_c,out_t_c):
            return None
        return _logarithmic_mean_(in_t_h-out_t_c,out_t_h-in_t_c)

    @staticmethod
    def epsilon(in_t_h:typing.Optional[float]=None,out_t_h:typing.Optional[float]=None,in_t_c:typing.Optional[float]=None,out_t_c:typing.Optional[float]=None,m_h:typing.Optional[float]=None,m_c:typing.Optional[float]=None,cp_h:typing.Optional[float]=None,cp_c:typing.Optional[float]=None):
        return _epsilon_(in_t_h,out_t_h,in_t_c,out_t_c,m_h,m_c,cp_h,cp_c)

    @staticmethod
    def ntu(epsilon:typing.Optional[float]=None,m_h:typing.Optional[float]=None,m_c:typing.Optional[float]=None,cp_h:typing.Optional[float]=None,cp_c:typing.Optional[float]=None):
        if None in (epsilon,m_h,m_c,cp_h,cp_c):
            return None
        c_h,c_c=_capacity_rates_(m_h,m_c,cp_h,cp_c)
        c=min(c_h,c_c)/max(c_h,c_c)
        if not 0<=epsilon<1:
            raise ValueError("In counterflow the effectiveness must be lower than 1")
        if mth.isclose(c,1,rel_tol=1e-9):
            return epsilon/(1-epsilon)
        return (1/(c-1))*mth.log((epsilon-1)/(epsilon*c-1))
