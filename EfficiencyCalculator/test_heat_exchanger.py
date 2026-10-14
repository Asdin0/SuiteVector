import math as mth
import pytest
import HeatExchanger


def counterflow_effectiveness(ntu, c):
    return (1 - mth.exp(-ntu * (1 - c))) / (1 - c * mth.exp(-ntu * (1 - c)))


def parallel_flow_effectiveness(ntu, c):
    return (1 - mth.exp(-ntu * (1 + c))) / (1 + c)


def test_counterflow_ntu_with_balanced_streams():
    assert HeatExchanger.CfHeatExchanger.ntu(epsilon=0.95, m_h=1, m_c=1, cp_h=1000, cp_c=1000) == pytest.approx(19)


def test_counterflow_ntu_is_continuous_near_balanced_streams():
    assert HeatExchanger.CfHeatExchanger.ntu(epsilon=0.95, m_h=1, m_c=1.0000001, cp_h=1000, cp_c=1000) == pytest.approx(19, rel=1e-5)


@pytest.mark.parametrize(["ntu", "c"], [(0.5, 0.3), (2.5, 0.7), (6, 0.95), (3, 1.0)])
def test_counterflow_ntu_inverts_effectiveness(ntu, c):
    epsilon = ntu / (1 + ntu) if c == 1.0 else counterflow_effectiveness(ntu, c)
    assert HeatExchanger.CfHeatExchanger.ntu(epsilon=epsilon, m_h=1, m_c=1 / c, cp_h=1000, cp_c=1000) == pytest.approx(ntu, rel=1e-9)


@pytest.mark.parametrize(["ntu", "c"], [(0.5, 0.3), (2.5, 0.7), (1.2, 1.0)])
def test_parallel_flow_ntu_inverts_effectiveness(ntu, c):
    epsilon = parallel_flow_effectiveness(ntu, c)
    assert HeatExchanger.EqHeatExchanger.ntu(epsilon=epsilon, m_h=1, m_c=1 / c, cp_h=1000, cp_c=1000) == pytest.approx(ntu, rel=1e-9)


def test_lmtd_with_equal_terminal_differences():
    assert HeatExchanger.CfHeatExchanger.logarithmic_mean_temperature_difference(in_t_h=300, out_t_h=100, in_t_c=95, out_t_c=295) == 5


def test_lmtd_accepts_zero_temperature():
    assert HeatExchanger.EqHeatExchanger.logarithmic_mean_temperature_difference(in_t_h=80, out_t_h=40, in_t_c=0, out_t_c=30) == pytest.approx(70 / mth.log(8))


def test_effectiveness_accepts_zero_outlet_temperature():
    assert HeatExchanger.CfHeatExchanger.epsilon(in_t_h=80, out_t_h=0, in_t_c=-20, m_h=1, m_c=1, cp_h=1, cp_c=1) == pytest.approx(0.8)


def test_temperature_cross_is_refused():
    with pytest.raises(ValueError):
        HeatExchanger.CfHeatExchanger.logarithmic_mean_temperature_difference(in_t_h=300, out_t_h=150, in_t_c=100, out_t_c=310)


def test_parallel_flow_effectiveness_limit_is_refused():
    with pytest.raises(ValueError):
        HeatExchanger.EqHeatExchanger.ntu(epsilon=0.7, m_h=1, m_c=1, cp_h=1, cp_c=1)


def test_effectiveness_outside_energy_balance_is_refused():
    with pytest.raises(ValueError):
        HeatExchanger.CfHeatExchanger.epsilon(in_t_h=300, out_t_h=50, in_t_c=100, m_h=1, m_c=1, cp_h=1, cp_c=1)


def test_missing_inputs_return_none():
    assert HeatExchanger.CfHeatExchanger.ntu(epsilon=None, m_h=1, m_c=1, cp_h=1, cp_c=1) is None
    assert HeatExchanger.EqHeatExchanger.logarithmic_mean_temperature_difference(in_t_h=300, out_t_h=None, in_t_c=100, out_t_c=200) is None
