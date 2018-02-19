from thermo.chemical import Chemical
from mass_opt import sweep_geometric_configs 
import matplotlib.pyplot as plt

class FlowProperties:

    def __init__(self, T, P, mass_flow, thermal_power):
        self.m_dot = mass_flow
        self.Q_therm = thermal_power
        self.T = float(T)
        self.P = float(P)
        self.secondary_properties()
        self.dP_allowed = 0
        self.mass = 0

    def secondary_properties(self):
        self.tol = Chemical('carbon dioxide')
        self.tol.calculate(T=self.T, P=self.P)
        
        self.Cp = self.tol.Cp
        self.mu = self.tol.mu
        self.k_cool = self.tol.k
        self.rho = self.tol.rho
        self.Pr = self.Cp * self.mu / self.k_cool

class PowerCycleSweep:
    
    def __init__(self, params):
        self.filename = params

    def load_params(self):
        """Load power cycle parameters.
        """
        self.params = []

        with open(self.filename, 'r') as fp:
            next(fp)
            for line in fp:
                data = line.split()
                inlet, outlet, bulklet = self.process_params(data)        
                # store parameters
                self.params.append({'in': inlet,
                                    'out' : outlet,
                                    'bulk' : bulklet})
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
        # get inlet cond
        inlet = FlowProperties(T_in, P_out, m_dot, Q_therm)
        inlet.dP_allowed = dP
        # get outlet cond
        outlet = FlowProperties(T_out, P_out, m_dot, Q_therm)
        outlet.dP_allowed = dP
        # get bulk cond
        bulklet = FlowProperties(T_bulk, P_bulk, m_dot, Q_therm)
        bulklet.dP_allowed = dP
        
        return inlet, outlet, bulklet

    def get_minimum_mass(self, Rs, PDs, z, c, N, axial_approx='bulk'):
        """Calculate minimum-mass reactor for every power cycle configuration.
        """
        x = []
        y = []
        for pc_config in self.params:
            sweepresults = sweep_geometric_configs(Rs, PDs, z, c, N,
                                                   pc_config[axial_approx])

            pc_config['bulk'].mass = sweepresults.min_mass
            
            x.append(pc_config['bulk'].Q_therm)
            y.append(pc_config['bulk'].mass)
        
        plt.plot(x, y)
        plt.show()

if __name__=='__main__':
    
    pc_data = PowerCycleSweep('CycleParameters.txt')
    pc_data.load_params()
    print(pc_data.params)
    pc_data.get_minimum_mass((0.005, 0.015), (1.01, 2), 0.15, 0.00031, 20)
