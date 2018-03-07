from ht_functions import Flow, ParametricSweep
from physical_constants import FlowProperties

def get_min_mass(temp, press, dp, mass_flow, thermal_power):
    """Callable wrapper for MATLAB returns minimum reactor mass for a set of
    input flow conditions.

    Arguments:
    ----------
        T (float): bulk coolant temperature [K]
        P (float): bulk coolant pressture [Pa]
        dp (float): pressure drop across reactor core [Pa]
        m_dot (float): coolant flow rate [kg/s]
        Q_therm (float): core thermal power [W]
    
    Returns:
    --------
        min_mass (float): minimum reactor core mass for given flow conditions.
    """
    # store primarty and calculate secondary flow properties
    properties = FlowProperties(T=temp, P=press, m_dot=mass_flow,
            Q_therm=thermal_power, dp_limit=dp)

    sweepresults = ParametricSweep()
    sweepresults.sweep_geometric_configs(props=properties)
    
    return sweepresults.get_min_mass()

if __name__ == '__main__':
    print(get_min_mass(1100, 1.7e7, 481000, 0.75, 131000))
