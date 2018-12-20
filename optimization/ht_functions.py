import math
from scipy.optimize import minimize_scalar
import numpy as np
# import physical properties
import physical_properties as pp


def pipeflow_turbulent(Re, Pr, LD, relrough):
    """Turbulent pipeflow correlation from EES. This function proiveds a nusselt
    number and friction factor for turbulent flow in a pipe.

    Arguments:
    ----------
        Re (float): flow Reynold's number [-]
        Pr (float): flow Prandtl number [-]
        LD (float): length over diameter [-]
        relrough (float) : relative roughness [-]
    
    Returns:
    --------
        Nusselt_L (float): nusselt number
        f (float): friction factor
    """

    #From Li, Seem, and Li, "IRJ, 
    #"A New Explicity Equation for Accurate Friction Factor 
    # Calculation for Smooth Tubes" 2011
    f_fd=(-0.001570232/math.log(Re) +
           0.394203137/math.log(Re)**2 +
           2.534153311/math.log(Re)**3) * 4 
    
    if relrough > 1e-5:
        #Offor and Alabi, 
        #Advances in Chemical Engineering and Science, 2016, 6, 237-245
        f_fd=(-2*math.log(
                 (relrough/3.71) -
                 (1.975/Re) * math.log((relrough/3.93)**1.092 +
                 (7.627/(Re + 395.9))), 10)) ** (-2)
    
    #Gnielinski, V.,, Int. Chem. Eng., 16, 359, 1976
    Nusselt_L= ((f_fd/8)*(Re-1000)*Pr)/(1+12.7*math.sqrt(f_fd/8)*(Pr **(2/3) - 1)) 
    
    if (Pr<0.5):
        # Notter and Sleicher, Chem. Eng. Sci., Vol. 27, 1972
        Nusselt_L_lp =4.8 + 0.0156 * Re**0.85 * Pr**0.93
        if (Pr<0.1):
            Nusselt_L = Nusselt_L_lp
        else:
            Nusselt_L = Nusselt_L_lp+(Pr-0.1)*(Nusselt_L-Nusselt_L_lp)/0.4

    #account for developing flow
    f=f_fd*(1+(1/LD)**0.7) 
    Nusselt_L*=(1+(1/LD)**0.7)
    
    return Nusselt_L, f

def pipeflow_laminar(Re, Pr, LD, relrough):
    """Laminar pipeflow correlation from EES. This function provides a nusselt
    number and friction factor for laminar flow in a pipe.

    Arguments:
    ----------
        Re (float): flow Reynold's numbe [-]
        Pr (float): flow Prandtl number [-]
        LD (float): length over diameter [-]
        relrough (float) : relative roughness [-]
    
    Returns:
    --------
        Nusselt_T (float): nusselt number at constant temperature
        Nusselt_H (float): nusselt number at constant heat flux
        f (float): friction factor
    """

    Gz = Re* Pr /LD     
    x = LD / Re
    fR = 3.44 / math.sqrt(x) +\
        (1.25/(4*x) + 16 - 3.44/math.sqrt(x)) /\
        (1 + 0.00021 * x**(-2))

    f = 4 * fR / Re
    Gm = Gz**(1/3)
    Nusselt_T = 3.66 + ((0.049+0.02/Pr)*Gz**1.12)/(1+0.065*Gz**0.7)
    Nusselt_H = 4.36 + ((0.1156 +0.08569 /Pr**0.4)*Gz)/(1+0.1158*Gz**0.6) 
    
    return Nusselt_T, Nusselt_H, f
    
def pipeflow_nd(Re, Pr, LD, relrough):
    """Nusselt correlation procedure from EES. This function provides a nusselt
    number and friction factor for flow in a pipe. It accounts for laminar,
    turbulent and transistional flow.
    
    Arguments:
    ----------
        Re (float): flow Reynold's numbe [-]
        Pr (float): flow Prandtl number [-]
        LD (float): length over diameter [-]
        relrough (float) : relative roughness [-]
    
    Returns:
    --------
        Nusselt_T (float): nusselt number at constant temperature
        Nusselt_H (float): nusselt number at constant heat flux
        f (float): friction factor
    """
    # smoothing factor
    m=6

    # get turbulent
    Nusselt_T, f = pipeflow_turbulent(Re, Pr, LD, relrough)
    Nusselt_H = Nusselt_T
   
    # get laminar
    Nusselt_lam_T, Nusselt_lam_H, f_lam = pipeflow_laminar(Re, Pr, LD, relrough)
    Nusselt_T = (Nusselt_lam_T**m + Nusselt_T**m)**(1/m)
    Nusselt_H = (Nusselt_lam_H**m + Nusselt_H**m)**(1/m)
    f = (f_lam**m + f**m)**(1/m)
        
    return Nusselt_T, Nusselt_H, f

def oned_flow_modeling(analyze_flow):
    """1D calculation.
    This function produces a valid, coolable reactor design given the following
    arguments:

    Arguments:
    ----------
        analyze_flow (flow) Flow object with methods and attributes to calculate
        N_channels.
    Returns:
    --------
        None
    """
    find_n_channels(analyze_flow)
    analyze_flow.adjust_dp()
    analyze_flow.calc_reactor_mass()

def _calc_n_channels_error(guess, flowiteration):
    """Calculate squared error between guess value and N channels for all
    three guess values.

    Arguments:
    ----------
        guess: a guess value for N_channels
        flowiteration: class containing TH methods and attributes

    Returns:
    --------
        error: difference between guess fuel channels and calculated required
        N_channels (float)
    """
    return flowiteration.compute_Q_from_guess(guess)


def find_n_channels(flow):
    """Perform error minimization. Using scipy's optimization package, call
    the _error function until the error is minimized to a set tolerance.

    Arguments:
    ----------
        flow: (class) Flow object. Contains attributes and
        methods required to perform an N_channels calculation for a single
        geometry (r, PD, L, c)
    Returns:
    --------
        none

    """
    res = minimize_scalar(_calc_n_channels_error, bounds=(0.2, 1), args=(flow),
                          method='Bounded', options={'xatol': 1e-10})

class Flow:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel.
    """

    ################################
    # UNIT SYSTEM: m, kg, J, W, Pa #
    ################################
    const_r = 0.1
    # set some guess values
    fuel_frac = 0.75  # number of required fuel channels for given flow conditions
    core_r = 1 # guess core radius
    T = 0
    def __init__(self, cool_r, c, AR, power, fuel, cool, clad, refl, flowprops,
                 order=3):
        """Initialize the flow iteration class.

        Initialized Attributes:
        --------------------
            Q_therm (float): core thermal power [W]
            fuel (str): reactor fuel
            coolant (str): reactor coolant
            AR (float): core aspect ratio [-]
            c (float): clad thickness
            r_channel (float): coolant channel radius
            pitch: fuel thickness (minor axis of hexagon) [m]
            L: length of core [m]
        """
        # load core parameters
        self.Q_therm = power
        self.fuel = fuel
        self.coolant = cool
        self.AR = AR
        self.c = c
        self.r_channel = cool_r

        # load material properties
        self.fuelprops = pp.fuel_props[fuel]
        self.reflprops = pp.refl_props[refl]
        self.cladprops = pp.clad_props[clad]
        self.pv_props  = pp.pv_props['SS304']
        self.fps = flowprops
        self.m_dot = flowprops.m_dot   
        # set up geometry
        self.neutronics(order)
        self.set_geom()
        self.dT = self.fuelprops['T_center'] - self.fps.T  # temp. drop fuel -> coolant

    def set_geom(self):
        """Setup the problem geometry.

        Modified Attributes:
        --------------------
            A_flow: flow area per fuel channel. [m^2]
            A_fuel: fuel area per fuel channel. [m^2]
            A_core: area of the core. [m^2]
            L: length of the core [m]
            LD: length over diameter for a coolant channel [-]
            vol_fuel: fuel volume [m^3]
            vol_cool: cool volume [m^3]
            N_channels: number of coolant channels [-]
            D_e: equivalent flow diameter [m]
        """
        self.A_core = self.core_r**2 * math.pi 
        self.clad_frac = 1 - self.r_channel**2 / (self.r_channel + self.c)**2
        self.A_flow = self.A_core * (1 - self.fuel_frac) * (1-self.clad_frac)
        self.A_clad = self.A_core * (1-self.fuel_frac) * (self.clad_frac)
        self.A_fuel = self.A_core * self.fuel_frac
        
        self.L = 2*self.AR * self.core_r
        self.LD = self.L / (self.r_channel*2)
        self.vol_fuel = self.A_fuel * self.L
        self.vol_cool = self.A_flow * self.L
        self.vol_clad = self.A_clad * self.L

        self.N_channels = self.A_flow / (self.r_channel**2 * math.pi)
        
        # hydraulic diametre
        self.D_e = 2.0 * self.r_channel

    def characterize_flow(self):
        """Calculate important non-dim and dim flow parameters. These parameters
        are required to determine generation per fuel channel.

        Correlations and equations used are from El-Wakil's Nuclear Heat
        Transport Textbook

        Modified Attributes:
        --------------------
            relrough (float): relative roughness [-]
            G_dot (float): mass flux [kg/m^2-s]
            v (float): flow velocity [m/s]
            Re (float): reynold's number [-]
            Nu (float): nusselt number [-]
            f (float): friction factor [-]
            h_bar (float): heat transfer coefficient [W/m^2-K]
        """
        # roughness of cladding
        self.rough = self.cladprops['rough']
        # relative roughness
        self.relrough = self.rough / (self.r_channel*2)
        # calculate mass flux
        self.G_dot = self.fps.m_dot / self.A_flow
        # calculate flow velocity from mass flux
        self.v = self.G_dot / self.fps.rho
        # calculate Reynolds Number
        self.Re = self.fps.rho * self.v * self.D_e / self.fps.mu
        # Nusselt Correlations from EES
        self.Nu, NuH, self.f = pipeflow_nd(self.Re, self.fps.Pr, self.LD, self.relrough)
        # heat transfer coefficient
        self.h_bar = self.Nu * self.fps.k_cool / self.D_e

    def get_q_per_channel(self):
        """Calculate achievable average volumetric generation:
        This method uses previously set geometry and calculated flow
        parameters to determine the maximum-achievable volumetric generation
        for each fuel channel. Calculates number of fuel channels required
        for desired thermal output.

        This method uses an equation from El-Wakil's Nuclear Heat Transport
        textbook.

        Modified Attributes:
        --------------------
            radius_cond (float): average distance of conduction in fuel [m]
            R_fuel: (float): resistance term (conduction in the fuel) [W/K]
            R_conv: (float): resistance term (convection to the coolant) [W/K]
            q_trip_max (float): maximum achievable thermal power [W]
            gen_Q (float): thermal power scaled for flux shape
        """
        
        self.radius_cond = math.sqrt(self.A_fuel / self.N_channels) / 2
        self.XS_A_cond = math.pi * self.r_channel * self.L * 2 * self.N_channels
        
        self.R_cond_fuel = self.radius_cond / (self.fuelprops['k'] * self.XS_A_cond)
        self.R_cond_clad = math.log(1 + self.c/self.r_channel)\
                        /(2*math.pi*self.cladprops['k']*self.L*self.N_channels)
        self.R_conv = 1 / (self.h_bar * self.r_channel *\
                           2 * math.pi * self.L * self.N_channels)
        self.R_tot = self.R_cond_fuel + self.R_cond_clad + self.R_conv 
        # calculate centerline volumetric generation
        q_trip_max = self.dT / (self.R_tot)

        # consider axial/radial flux variation (El-Wakil 4-30a)
        self.gen_Q = q_trip_max * 0.275
        #if self.Re < 3000 and self.Re > 2300:
        #    print(self.fuel_frac, self.gen_Q - self.Q_therm)

    def calc_dp(self):
        """Calculate axial pressure drop across the reactor core.

        Modified Attributes:
        --------------------
            dp (float): core pressure drop [Pa]
        """
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = self.f * self.L * self.fps.rho * \
            self.v * self.v / (2*self.D_e)

    def adjust_dp(self):
        """Check for pressure constraint. This method calls calc_dp() to get
        the pressure drop in the current condition. It checks the dp against the
        power cycle-constrained allowable dp. If the pressure is too high, it
        adjusts N_channels to the min N_channels that satisfies the dp
        constraint.

        Modified Attributes:
        --------------------
            N_channels (int): number of fuel channels [-]
        """

        self.calc_dp()
        while self.fps.dp_limit - self.dp < -1e-4:
            # set N_channels and guess_channels
            self.fuel_frac = self.get_dp_constrained_Nchannels()
            self.constrain_radius()
            self.set_geom()
            self.characterize_flow()
            self.calc_dp()

    def get_dp_constrained_Nchannels(self):
        """Set the N_channels based on the allowable dP. This method
        calculates the required number of channels to meet the pressure drop
        constraint (set by the power cycle).

        Arguments:
        ----------
            self: Flow object [-]
        Returns:
        --------
            req_channels (int): Min N_channels required to meet dp constraint [-].
        """
        v_req = math.sqrt(2*self.D_e * self.fps.dp_limit /
                          (self.f * self.L * self.fps.rho))
        
        req_G_dot = v_req * self.fps.rho
        req_A_flow = self.fps.m_dot / req_G_dot
        
        # get total channel area
        req_A_channels = req_A_flow / (1-self.clad_frac)
        
        # get required fuel fraction    
        req_fuel_frac = 1 - (req_A_channels / self.A_core)
        
        return req_fuel_frac

    def compute_Q_from_guess(self, inp_guess):
        """Perform single 1D heat flow calculation. This method calls the
        required methods to perform one iteration of the calculation.

            Arguments:
            ----------
                inp_guess: (float) guess value for fuel fraction
            Returns:
            --------
                error: (float) squared error between the possible power
                generation and the required power generation
        """

        self.fuel_frac = inp_guess
        self.constrain_radius() 
        self.set_geom()
        self.characterize_flow()
        self.get_q_per_channel()
        
        return (self.gen_Q - self.Q_therm)**2

    def PV_thickness(self):
        """Calculate required pressure vessel thickness, volume
        """
    
        R = self.r_ref
        # from faculty.washington.edu/vkumar/me356/pv_rules.pdf
        self.r_PV = R + R*self.fps.P / (self.pv_props['strength'] + 0.6*self.fps.P)
        
        
    def calc_reactor_mass(self):
        """Based on results of the iteration, calculate the reactor mass.

        Modified Attributes:
        --------------------
            Vol_fuel: total fuel volume [m^3]
            mass: total fuel mass[kg]
        """
        self.refl_m = {'UW'  : {'CO2' : 1.211337},
                       'UO2' : {'CO2' : 1.08}
                      }
        self.r_ref = self.refl_m[self.fuel][self.coolant]*self.core_r
        
        self.PV_thickness()
        # calculate volumes
        self.vol_core = self.vol_fuel + self.vol_cool + self.vol_clad
        self.vol_refl = (self.r_ref**2 * 2*self.r_ref * math.pi) - self.vol_core
        self.vol_PV = (self.r_PV**2 - self.r_ref**2) * 2*self.r_PV * math.pi
        # calculate masses
        self.fuel_mass = self.vol_fuel * self.fuelprops['rho_fuel']
        self.cool_mass = self.vol_cool * self.fps.rho
        self.clad_mass = self.vol_clad * self.cladprops['rho']
        self.refl_mass = self.vol_refl * self.reflprops['rho']
        self.PV_mass = self.vol_PV * self.pv_props['rho']
        # total mass
        self.mass = self.fuel_mass + self.cool_mass + self.refl_mass +\
                    self.clad_mass + self.PV_mass
    
    def neutronics(self, order):
        """
        """
        fp = open('../data/{0}_{1}_core_radii.txt'.format(self.fuel,
                                                          self.coolant), 'r')
        lines = fp.readlines()
        fp.close()
        
        f = []; r = []
        for line in lines[1:]:
            f.append(float(line.split(',')[0]))
            r.append(float(line.split(',')[1]) / 100)

        self.poly = np.polyfit(f, r, order)

    def constrain_radius(self):
        """Constrain the core radius based on criticality requirements.

        Modified Attributes:
        --------------------
            core_r (float): critical core radius based on fuel fraction [m]
        """

        # Critical Radius of Buried Reactor
        coeffs = {'UW'  : {'CO2' : (-85.8968, 196.3772, -161.3156, 62.3759)},
                  'UO2' : {'CO2' : ( -63.9328861, 147.44434591, -122.88976054, 48.38235672)}
                 }
        
        self.core_r = self.const_r #np.polyval(self.poly, self.fuel_frac)
