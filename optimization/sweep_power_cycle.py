from thermo.chemical import Chemical


class FlowProperties:

    def __init__(self, T, P):
        self.T = float(T)
        self.P = float(P)
        self.secondary_properties()

    def secondary_properties(self):
        self.tol = Chemical('carbon dioxide')
        self.tol.calculate(T=self.T, P=self.P)
        
        self.Cp = self.tol.Cp
        self.mu = self.tol.mu
        self.k = self.tol.k
        self.rho = self.tol.rho

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
                # get inlet cond
                inlet = FlowProperties(data[0], data[2])
                # get outlet cond
                outlet = FlowProperties(data[1], data[3])

                self.params.append({'in': inlet,
                                    'out' : outlet,
                                    'm_dot' : data[4],
                                    'Q_therm' : data[6]})
    def minimize_TH_mass(self):




if __name__=='__main__':
    
    pc_data = PowerCycleParams('CycleParameters.txt')
    pc_data.load_pc_params()

    print(pc_data.params)
