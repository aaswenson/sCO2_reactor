from thermo.chemical import Chemical
from mass_opt import sweep_geometric_configs 
from mcnp_inputs import PinCellMCNP
from scale_inputs import PinCellSCALE
import matplotlib.pyplot as plt
import numpy as np

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

    def __init__(self, T, P, mass_flow, thermal_power, efficiency):
        """Inialize FlowProperties class and load required flow property data.

        Initialized Attributes:
        -----------------------
            eta: (float) cycle efficiency [-]
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_allowed: (float) power-cycle constrained dp [Pa]
            mass: (float) required reactor mass
        """
        self.eta = efficiency
        self.m_dot = mass_flow
        self.Q_therm = thermal_power
        self.T = float(T)
        self.P = float(P) * 1e3
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
    """Load power cycle data and calculate minimum-mass reactor for every power
    cycle configuration.
    """
    def __init__(self, power_cycle_input):
        """
        """
        self.filename = power_cycle_input
        self.names = ['Q_therm', 'm_dot', 'T', 'P', 'eta', 'mass', 'TH_results',
                'flow_props']
        self.types = ['f8']*6 + ['O', 'O']

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
                for key in self.names[0:5]:
                    self.cycle_parameters[key][idx] = bulklet.__dict__[key]
                    self.cycle_parameters['flow_props'][idx] = bulklet
                    
        # close file
        fp.close()

    def process_params(self, data):
        """Calculate secondary properties for inlet and outlet conditions.
        Calculate bulk-averaged properties/secondary properties.
        
        Modified Attributes:
        --------------------
            None
        Arguments:
        ----------
            data: (list) one line of power cycle parameter csv file
        Returns:
        --------
            bulklet: (FlowProperties class) object with primary and secondary
            flow properties required for TH analysis.
        """
        # load data from input line
        T_in = float(data[0]); T_out = float(data[1])
        P_in = float(data[2]); P_out = float(data[3])
        m_dot = float(data[4])
        eta = float(data[5])
        Q_therm = float(data[6])

        # bulk-averaged flow properties
        T_bulk = (T_in + T_out)/ 2.0
        P_bulk = (P_in + P_out) / 2.0
        # get bulk secondary conditions
        bulklet = FlowProperties(T_bulk, P_bulk, m_dot, Q_therm, eta)
        # pressure drop limit 
        bulklet.dP_allowed = P_in - P_out
        
        return bulklet

    def get_minimum_mass(self, Rs, PDs, z, c, N):
        """Calculate minimum-mass reactor for every power cycle configuration.
        """
        for pc_config in self.cycle_parameters:
            sweepresults = sweep_geometric_configs(Rs, PDs, z, c, N,
                                                   pc_config['flow_props'])
            # save minimum mass reactor for each pc configuration
            TH_min_reactor = sweepresults.get_min_mass()
            pc_config['mass'] = TH_min_reactor.mass
            pc_config['flow_props'] = pc_config
            pc_config['TH_results'] = TH_min_reactor

        # get minimum-mass configuration to create neutronics model
        pc_min_idx = list(self.cycle_parameters['mass']).\
                          index(min(self.cycle_parameters['mass']))
        
        reactor = self.cycle_parameters[pc_min_idx]['TH_results']
        print(reactor.__dict__)
        print(reactor.fps.__dict__)
        #self.write_mcnp_input(reactor)

    def write_mcnp_input(self, TH_reactor):
        """ Write the mcnp input.
        """

        input = PinCellSCALE(TH_reactor)
        input.write_fuel_string(0.5, ('Nitrogen', 1))
        input.write_mat_string('Inconel-718', 'Carbon Dioxide')
        input.write_input()

    def fit_curve(self):
        """Fit a curve to the mass, q_therm data.
        """
        res = np.polyfit(self.cycle_parameters['Q_therm'], 
                         self.cycle_parameters['mass'],
                         2)

        return res
    
    def plot(self):
        """
        """
        plots = [
      [('mass', 'mass'), ('mass', 'Q_therm'), ('mass', 'T'), ('mass', 'm_dot')],
      [('Q_therm', 'mass'), ('Q_therm', 'Q_therm'), ('Q_therm', 'T'), ('Q_therm', 'm_dot')],
      [('T', 'mass'), ('T', 'Q_therm'), ('T', 'T'), ('T', 'm_dot')],
      [('m_dot', 'mass'), ('m_dot', 'Q_therm'), ('m_dot', 'T'), ('m_dot', 'm_dot')]]
    
        axis_labels = {'mass' : 'Fuel mass [kg]',
                       'Q_therm' : 'Q therm [W]',
                       'T' : 'Coolant T [K]',
                       'm_dot' : 'Cool. flow rate [kg/s]'}

        pass_idx = [4, 8, 9, 12, 13, 14] 
        
        fig = plt.figure()

        for xidx, row in enumerate(plots):
            for yidx, plot in enumerate(row): 
                plot_id = yidx*4 + xidx
                if plot_id not in pass_idx:
                    ax = fig.add_subplot(4, 4, plot_id+1)
                    xkey = plot[0] 
                    ykey = plot[1]
                    x = self.cycle_parameters[xkey]
                    y = self.cycle_parameters[ykey]
                    ax.scatter(x, y, s=5)
                    if ykey == 'Q_therm':
                        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                    if xkey == 'Q_therm':
                        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
                    if yidx == 0:
                        plt.title(axis_labels[xkey], fontsize=12)
                    if xidx == yidx:
                        plt.ylabel(axis_labels[ykey], fontsize=12)
        #plt.tight_layout(pad=0.001)
        fig.savefig('mass_vs_all.png')
        self.plot_mass_vs_Qtherm()
        
        return plt

    def plot_mass_vs_Qtherm(self):
        """
        """
        fig = plt.figure()
        plt.scatter(self.cycle_parameters['Q_therm'],
                    self.cycle_parameters['mass'], s=12)
        plt.xlabel('Q_therm [W]')
        plt.ylabel('fuel mass [kg]')
        plt.title('mass vs. thermal power (TH considerations only)')
        fig.savefig('mass_vs_therm.png')
       
        return plt

if __name__=='__main__':
    
    pc_data = PowerCycleSweep('CycleParameters_EvenSampling.csv')
    pc_data.load_params()
    pc_data.get_minimum_mass((0.005, 0.015), (1.1, 2), 0.15, 0.00031, 5)
    res = pc_data.fit_curve()
    plt = pc_data.plot()
    plt.show()
    print(res)
