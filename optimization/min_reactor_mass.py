from ht_functions import Flow, ParametricSweep
from physical_constants import FlowProperties

def get_min_mass(T, P, dp, m_dot, Q_therm):
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
    properties = FlowProperties(T, P, m_dot, Q_therm)
    properties.dP_allowed = dp

    sweepresults = ParametricSweep()
    sweepresults.sweep_geometric_configs(flowprops=properties)
    print(properties.__dict__)
    print(FlowProperties.__dict__)
    
    return sweepresults.get_min_mass()


