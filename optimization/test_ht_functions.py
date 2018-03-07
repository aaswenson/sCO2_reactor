import math
import pytest
from random import uniform
from ht_functions import Flow, oned_flow_modeling

# parameters for test case
radius = 0.005
PD = 2
c = 0.00031
L = 0.5
N = 6


def test_set_geom():
    """Test geometry initialization.
    """
    test = Flow(radius, PD, c, L)
    test.guess = N
    # expected values
    exp_De = radius * 2.0
    exp_A_flow = math.pi * radius**2
    exp_A_fuel = math.sqrt(3) * ((radius + c) * PD * 2)**2 / 2.0 -\
        (radius + c)**2 * math.pi
    # calculated values
    test.set_geom()
    # compare
    assert exp_De == test.D_e
    assert exp_A_flow == test.A_flow
    assert exp_A_fuel == test.A_fuel


def test_characterize_flow():
    """Test flow characterization.
    """
    test = Flow(radius, PD, c, L)
    test.guess_channels = N
    # get geom
    test.set_geom()
    # expected flow velocity
    exp_G = test.fps.m_dot / (test.A_flow * N)
    exp_v = exp_G / test.fps.rho
    # expected heat transfer coefficient
    exp_Re = test.fps.rho * exp_v * test.D_e / test.fps.mu
    exp_Nu = 0.023*math.pow(exp_Re, 0.8)*math.pow(test.fps.Pr, 0.4)
    exp_f = 0.184 / math.pow(exp_Re, 0.2)
    exp_h = exp_Nu * test.fps.k_cool / test.D_e
    # calculated values
    test.characterize_flow()
    # compare
    assert exp_v == test.v
    assert exp_f == test.f
    assert exp_h == test.h_bar


def test_q():
    """Test q_per_channel calculation.
    """
    test = Flow(radius, PD, c, L)
    test.guess_channels = N
    # get geom and flow conditions
    test.set_geom()
    test.characterize_flow()
    # expected value
    exp_q_per_channel = 2.43123E+04
    test.get_q_per_channel()
    # compare
    assert abs(exp_q_per_channel - test.q_per_channel) < 1.0


def test_dp():
    """Test subchannel dp calculation.
    """
    test = Flow(radius, PD, c, L)
    test.guess_channels = N
    test.set_geom()
    test.characterize_flow()
    test.calc_dp()
    # expected dp
    exp_dp = test.f * L * test.fps.rho * test.v ** 2 / (2*test.D_e)
    # compare
    assert (exp_dp - test.dp)**2 < 1e-5

def test_flow_calc():
    """
    """
    test = Flow(radius, PD, c, L)
    oned_flow_modeling(test)
    
    assert math.ceil(test.N_channels) == 5
    
def test_rand_flow_calc():
    # random check oned_flow_modeling
    rand_r = uniform(0.005, 0.01)
    rand_PD = uniform(1.1, 2)
    rand_L = uniform(0.15, 0.5)

    obs = Flow(rand_r, rand_PD, c, rand_L)
    oned_flow_modeling(obs)
    # check result manually
    exp = Flow(rand_r, rand_PD, c, rand_L)
    exp.guess_channels = obs.N_channels
    # get geom and flow conditions
    exp.set_geom()
    exp.characterize_flow()
    # expected value
    exp.get_q_per_channel()
    
    assert abs(exp.q_per_channel - obs.q_per_channel) < 1.0
    
    
