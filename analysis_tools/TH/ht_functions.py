# Matplotlib imports
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.axis
from matplotlib import cm, rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter
# Other Imports
import math
import numpy as np
import operator
import sys
from scipy.optimize import minimize, minimize_scalar
from physical_constants import *
import matplotlib.pyplot as plt
    
def _error(guess, flowiteration):
    """ Calculate squared error between guess value and N channels for all
    three guess values.
    """
    flowiteration.guess = guess
    flowiteration.mass_flux_channel()
    flowiteration.calc_nondim()
    flowiteration.get_h_bar()
    flowiteration.get_q_bar()
    flowiteration.calc_aspect_ratio()
    flowiteration.error = (flowiteration.guess - flowiteration.N_channels)**2
    
    return flowiteration.error

class FlowIteration:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel. Upon instantiation,
    the user provides:
        flow_radius (float): radius of the coolant flow channels [m]
        PD (float): fuel cell pitch / coolant channel diameter [-]
        c (float): channel cladding thickness [m]
        L (float): fuel length [m]
        guess (int): initial number of flow channels [-]
    """
    # geometric attributes
    r_channel = 0; PD = 0; c = 0; L = 0; Vol_fuel = 0; AR = 0
    guess = 0; N_channels = 0; error = 0
    A_flow = 0; A_fuel = 0;
    mass = 0
    # temperature drop
    dt = T_centerline - T_bulk
    # heat transfer and non-dimensional coefficients
    h_bar = 0; Re = 0; Nu = 0; Pr = 0
    R_fuel = 0; R_clad = 0; R_conv = 0; R_tot = 0;
    # flow parameters
    G_dot = 0; D_e = 0; v = 0; dp = 0
    # heat generation
    q_bar = 0; q_per_channel = 0; q_therm_check = 0
    
    def __init__(self, diameter, PD, c, L):
        self.r_channel = diameter / 2.0
        self.c = c
        self.pitch = (self.r_channel + self.c) * PD * 2
        self.L = L
        self.iterations = 0

    def mass_flux_channel(self):
        """Calculate coolant mass flux through 1 CERMET flow channel.
        Arguments: flow_radius (float) flow channel radius
                   self.N_channels (int) number of fuel elements
        Returns: G_dot (float): coolant mass flux
        """
        self.A_flow = self.r_channel ** 2 * math.pi
        self.A_fuel = math.sqrt(3)*self.pitch**2 / 2.0 -\
                      (self.r_channel + self.c) ** 2 * math.pi
        self.G_dot = m_dot / (self.A_flow * self.guess)
        self.D_e = 2.0 * self.r_channel
        self.v = self.G_dot / rho_cool

    def calc_nondim(self):
        """ Calculate Reynolds number
        """
        self.Re = rho_cool * self.v * self.D_e / mu
        self.Pr = Cp_cool * mu / k_cool
        # nusselt correlation is Eq (9-22) from El-Wakil Nuclear Heat Transport
        self.Nu = 0.023*math.pow(self.Re,0.8)*math.pow(self.Pr, 0.4)
    
    def get_h_bar(self):
        """Calculate average heat transfer coefficient.
        """
        self.h_bar = self.Nu * k_cool / self.D_e

    def get_q_bar(self):
        """Calculate heat transfer from fuel element. Convert q_bar max from
        cylindrical geometry to hexagonal geometry. Calculate number of required
        flow channels.
        """
        r_i = self.r_channel + self.c
        r_o = self.pitch / math.sqrt(3)
       
        # El Wakil (9-62) calculates max q''' at axial centerline
        self.R_fuel = (r_o**2 / (4*k_fuel)) * ((r_i/r_o)**2 - 2*math.log(r_i/r_o) - 1)
        self.R_clad = (r_o/2)**2 *\
        (1-(r_i/r_o)**2)*(math.log(r_i/(r_i-self.c))/k_clad)
        self.R_conv = (r_o/2)**2 *\
        (1-(r_i/r_o)**2)*(1/(self.h_bar*(r_i - self.c)))
        self.R_tot = self.R_fuel + self.R_clad + self.R_conv
        
        # calculate centerline volumetric generation
        q_trip_max = self.dt / self.R_tot
        
        # consider axial flux variation
        self.q_per_channel = 2 * q_trip_max * self.A_fuel * self.L / math.pi
        
        # calculate required number of channels
        self.N_channels = Q_therm / self.q_per_channel
        # calculate total fuel volume and q_bar
        self.Vol_fuel = self.A_fuel * self.L * self.N_channels
        self.q_bar = Q_therm / self.Vol_fuel
        
    def calc_dp(self):
        """Calculate pressure drop subchannel
        """
        # El Wakil (9-4)
        f = 0.184 / math.pow(self.Re, 0.2)
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = f * self.L * rho_cool * self.v * self.v / (2*self.D_e)

    def calc_aspect_ratio(self):
        total_area = (self.A_fuel + self.A_flow) * self.N_channels
        equivalent_radius = math.sqrt(total_area / math.pi) 
        self.AR = self.L / (2*equivalent_radius)

    def Iterate(self):
        """Perform Flow Calc Iteration
        """
        res = minimize_scalar(_error, bounds=(1, 1e9), args=(self),
        method='Bounded')
        self.calc_dp()
        self.N_channels = math.ceil(self.N_channels)

    def calc_reactor_mass(self):
        """Based on results of the iteration, calculate the reactor mass.
        """
        self.mass = self.Vol_fuel * rho_fuel

class ParametricSweep():
    """Class to store results of parametric sweeps for 1D flow channel analysis.
    
    """
    titles = {'mass' : ("Total Fuel Mass", "m [kg]"),
              'Re' : ("Reynolds Number", "Re [-]"),
              'G_dot' : ("Reactor Mass Flux", "G_dot [kg/m^2-s]"),
              'N_channels' : ("Number of Fuel Channels", "N Channels [-]"),
              'Nu' : ("Nusselt Number", "Nu [-]"),
              'dp' : ("Subchannel Pressure Drop", "dP [Pa]"),
              'h_bar' : ("Heat Transfer Coefficient", "h [W / m^2 - K]"),
              'q_per_channel' : ("Total Subchannel Generation", "q/channel [W]"),
              'q_bar' : ("Average Volumetric Generation", "q_bar [W/m^3]"),
              'v' : ("Flow Velocity", "v [m/s]"),
              'R_fuel' : ("Resistance to Conduction in Fuel", "R_fuel [K/W]"),
              'R_clad' : ("Resistance to Conduction in Clad", "R_clad [K/W]"),
              'R_conv' : ("Resistance to Convection", "R_conv [K/W]"),
              'R_tot' : ("Total Resistance to Heat Transfer", "R_tot [K/W]"),
              'AR' : ("Approximate Core Aspect Ratio", "AR [-]")
             }

    D = 0; PD = 0;
    
    # dict to save data for plotting
    data = {k: [] for k in titles.keys()}
    min_idx = None; min_jdk = None; min_mass = 0; minD = 0; min_PD = 0;
    select_AR = False
    
    def __init__(self, D, PD, N, select_AR):
        self.D = D
        self.PD = PD
        self.N = N
        self.select_AR = select_AR
        for key in self.data:
            self.data[key] = np.empty([N,N])

    def save_iteration(self, iteration, i, j):
        """ Save the data from each iteration of the parametric sweep. 
        """
        for key in self.data.keys():
            self.data[key][i][j] = iteration.__dict__[key]

    def get_min_data(self):
        """ After the parametric sweep is complete, find the minimum calculated
        fuel mass corresponding with a valid aspect ratio. Save the flowdata for that
        point, display the PD, D and mass for the optimized configuration.
        """
        
        if False not in np.isnan(self.data['mass']):
            print("The range you have selected does not satisfy pressure drop\
 requirements, consider increasing coolant flow channel diameter")
            sys.exit()

        not_AR = 1
        if self.select_AR == True:
            not_AR = np.nan

        # calculate departure from AR = 1
        aspect_err = np.subtract(self.data['AR'], np.ones([self.N, self.N]))
        aspect_err = np.abs(aspect_err)
        # select AR error no greater than 0.33 (based on NuScale AR = 1.33)
        on_off = lambda x: not_AR if x > 0.33 else 1
        vfunc = np.vectorize(on_off)

        for key in self.data:
            self.data[key] = np.multiply(vfunc(aspect_err), self.data[key])
        
        if False not in np.isnan(self.data['AR']):
            print("The range you have selected does not satisfy aspect ratio\
 requirements, consider a new geometry range")
            sys.exit()
        
        # search the results for minimum-mass configuration
        mindices_mass = np.unravel_index(
                   np.nanargmin(self.data['mass']), self.data['mass'].shape)
        # get data for optimal AR
        self.min_idx, self.min_jdx = mindices_mass
        self.min_mass = self.data['mass'][self.min_idx][self.min_jdx]

        # save the optimal configuration
        self.minD = self.D[self.min_idx][self.min_jdx]
        self.minPD = self.PD[self.min_idx][self.min_jdx]
        
        # report the optimal configuration and it's corresponding fuel mass 
        outstring =  "1D Thermal Hydraulics Optimization Results:\n"
        outstring += "Reactor minimum mass (m = " + str(round(self.min_mass,3))\
        + "[kg]) " + " occurs at D = " + str(self.minD) +  "[m] & PD = "\
        + str(self.minPD) + "[-]."
        print(outstring)

    def plot(self, D, PD, key):
        """Produce surface plot of the flow results as function of PD and coolant
        channel diameter.
        """
        # get parametric sweep data
        M = self.data[key]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(D, PD, M, 
                               cmap=cm.viridis, linewidth=0,
                               vmin=0, vmax=np.nanmax(M),
                               antialiased=False)

        # set x/y axis labels, ticks
        ax.set_xlabel("Coolant Channel Diameter [m]", fontsize=7)
        plt.xticks(rotation=25, fontsize=6)
        ax.set_ylabel("Fuel Pitch to Coolant D Ratio [-]", fontsize=7)
        plt.yticks(rotation=25, fontsize=6)
        ax.set_zlabel(self.titles[key][1], fontsize=7)
         
        # Customize the z axis.
        ax.set_zlim(np.nanmin(M),np.nanmax(M))
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        
        # edit z tick labels
        for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(6)
        niceMathTextForm = ScalarFormatter(useMathText=True)
        ax.w_zaxis.set_major_formatter(niceMathTextForm)
        ax.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
        plt.title(self.titles[key][0])
        
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5, format='%.0e')
        
        return plt
