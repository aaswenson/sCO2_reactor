from ht_functions import FlowIteration, StoreIteration
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.axis
from matplotlib import cm, rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter
import numpy as np
import argparse
import sys

def sweep_configs(D, PD, z, c, N, key):
    """Perform parametric sweep through pin cell geometric space.
    """
    # calculate appropriate step sizes given range
    D_step = (D[1] - D[0]) / N
    PD_step = (PD[1] - PD[0]) / N
    # ranges for diameter and pitch/diameter ratio
    D = np.arange(D[0], D[1], D_step)
    PD = np.arange(PD[0], PD[1], PD_step)

    # create parameter mesh
    D, PD = np.meshgrid(D, PD)
    M = np.empty([N,N])
    min_mass = 1e9
    # sweep through parameter space, calculate min mass
    for i in range(N):
        for j in range(N):
            flowdata = FlowIteration(D[i,j], PD[i,j], c, z, 10)
            flowdata.Iterate()
            flowdata.calc_reactor_mass()
            M[i][j] = flowdata.__dict__[key]
            if flowdata.mass < min_mass:
                min_mass = flowdata.mass
    print(min_mass)
    return D, PD, M

def plot_mass(D, PD, M, key):
    """Produce surface plot of the reactor mass as function of PD and coolant
    channel diameter.
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(D, PD, M, cmap=cm.viridis, linewidth=0,
            antialiased=False)
    
    ax.set_xlabel("Coolant Channel Diameter [m]", fontsize=7)
    plt.xticks(rotation=25, fontsize=6)
    ax.set_ylabel("Fuel Pitch to Coolant D Ratio [-]", fontsize=7)
    plt.yticks(rotation=25, fontsize=6)
    ax.set_zlabel(key, fontsize=7)
     
    # Customize the z axis.
    ax.set_zlim(np.min(M),np.max(M))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # edit z tick labels
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(6)
    niceMathTextForm = ScalarFormatter(useMathText=True)
    ax.w_zaxis.set_major_formatter(niceMathTextForm)
    ax.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
    plt.title(key)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5, format='%.0e')
    
    savename = key + '.png'
    plt.savefig(savename, dpi=500)
#    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("d_lower", type=float, help="channel D lower lim [m]")
    parser.add_argument("d_upper", type=float, help="channel D upper lim [m]")
    parser.add_argument("pd_lower", type=float, help="PD lower lim [m]")
    parser.add_argument("pd_upper", type=float, help="PD upper lim [m]")
    parser.add_argument("z", type=float, help="axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")
    parser.add_argument("steps", type=int, help="parameter resolution")
    parser.add_argument("plotkey", type=str, help="parameter parameter to plot")

    args = parser.parse_args()
    
    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel\
diameter! Set min PD > 1!")
        sys.exit()

    Diameters, PDs, Masses = sweep_configs((args.d_lower, args.d_upper),
                                   (args.pd_lower, args.pd_upper), 
                                   args.z, args.clad_t, args.steps, args.plotkey)
    # plot results
    plot_mass(Diameters, PDs, Masses, args.plotkey)

