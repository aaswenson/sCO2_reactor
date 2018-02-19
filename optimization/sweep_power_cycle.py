from thermo.chemical import Chemical
from mass_opt import sweep_geometric_configs 
import matplotlib.pyplot as plt
import numpy as np

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

    def __init__(self, T, P, mass_flow, thermal_power):
        """Inialize FlowProperties class and load required flow property data.

        Initialized Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_allowed: (float) power-cycle constrained dp [Pa]
            mass: (float) required reactor mass
        """
        self.m_dot = mass_flow
        self.Q_therm = thermal_power
        self.T = float(T)
        self.P = float(P)
        self.secondary_properties()
        self.dP_allowed = 0

    def secondary_properties(self):
        """Calculate secondary properties from primary flow properties
        
        Modified Attributes:
        --------------------
            Cp: (float) specific heat [kJ/kg-K]
            mu: (float) dynamic viscosity [kg/m-s]
            k_cool: (float) coolant conductivity [W/m-k]
            rho: (float) coolant density [kg/m^3]
            Pr: (float) cooland Prandtl number [-]
        """
        # access CO2 data from Thermodynamic Library
        tol = Chemical('carbon dioxide')
        tol.calculate(T=self.T, P=self.P)
        # calculate Cp, mu, k, rho, Pr    
        self.Cp = tol.Cp
        self.mu = tol.mu
        self.k_cool = tol.k
        self.rho = tol.rho
        self.Pr = self.Cp * self.mu / self.k_cool

class PowerCycleSweep:
    
    def __init__(self, power_cycle_input):
        self.filename = power_cycle_input
        self.names = ['Q_therm', 'm_dot', 'T', 'P', 'mass', 'flow_props']
        self.types = ['f8']*5 + ['O', 'O']

    def load_params(self):
        """Load power cycle parameters.
        """
        N = 0
        with open(self.filename, 'r') as fp:
            next(fp)
            # get number of lines
            for line in fp:
                N += 1
            # return to second line of file
            fp.seek(1,0)
            next(fp)
            self.cycle_parameters = np.ndarray((N,),
                    dtype={'names' : self.names, 'formats' : self.types})
            # loop through file, save data as you go
            for idx, line in enumerate(fp):
                data = line.split(',')
                bulklet = self.process_params(data)        
                # store float parameters
                for key in self.names[0:3]:
                    self.cycle_parameters[key][idx] = bulklet.__dict__[key]
                    self.cycle_parameters['flow_props'][idx] = bulklet
                    
        # close file
        fp.close()

    def process_params(self, data):
        """Calculate secondary properties for inlet and outlet conditions.
        Calculate bulk-averaged properties/secondary properties.
        """
        T_in = float(data[0]); T_out = float(data[1])
        P_in = float(data[2]); P_out = float(data[3])
        m_dot = float(data[4])
        Q_therm = float(data[6])

        dP = P_in - P_out
        # bulk-averaged
        T_bulk = (T_in + T_out)/ 2.0
        P_bulk = (P_in + P_out) / 2.0
        # get bulk cond
        bulklet = FlowProperties(T_bulk, P_bulk, m_dot, Q_therm)
        bulklet.dP_allowed = dP
        
        return bulklet

    def get_minimum_mass(self, Rs, PDs, z, c, N):
        """Calculate minimum-mass reactor for every power cycle configuration.
        """
        for pc_config in self.cycle_parameters:
            sweepresults = sweep_geometric_configs(Rs, PDs, z, c, N,
                                                   pc_config['flow_props'])
            # save minimum mass
            pc_config['mass'] = sweepresults.min_mass
            

    def fit_curve(self):
        """Fit a curve to the mass, q_therm data.
        """
        res = np.polyfit(self.cycle_parameters['Q_therm'], 
                         self.cycle_parameters['mass'],
                         1)

        return res
    
    def plot(self, xkey, ykey):
        plt.scatter(self.cycle_parameters[xkey], self.cycle_parameters[ykey])
        plt.title(ykey + " vs. " + xkey + " (Thermal considerations only)")
        plt.xlabel(xkey)
        plt.ylabel(ykey)
        plt.show()

if __name__=='__main__':
    
    pc_data = PowerCycleSweep('CycleParameters.csv')
    pc_data.load_params()
    pc_data.get_minimum_mass((0.005, 0.015), (1.01, 2), 0.15, 0.00031, 5)
    res = pc_data.fit_curve()
    pc_data.plot('T', 'Q_therm')
    print(res)
