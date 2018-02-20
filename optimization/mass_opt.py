""" Thermal-Hydraulic Optimization Module.
This module contains functions to perform a parametric sweep of reactor
geometric parameters. It calculates fuel mass and flow characteristics for each
set of valid geometric parameters.

Functions contained in this module:
    *oneD_flow_modeling
    *sweep_configs
"""
# Other imports
import numpy as np
import argparse
import sys
# Import TH functions
from ht_functions import Flow, ParametricSweep, oned_flow_modeling
from plot import plot


def sweep_geometric_configs(radii, pds, z, c, N, properties):
    """Perform parametric sweep through pin cell geometric space. Calculate the
    minimum required mass for TH purposes at each point.
    """
    # calculate appropriate step sizes given range
    R_step = (radii[1] - radii[0]) / N
    PD_step = (pds[1] - pds[0]) / N
    # ranges for diameter and pitch/diameter ratio
    R_array = np.arange(radii[0], radii[1], R_step)
    PD_array = np.arange(pds[0], pds[1], PD_step)

    # create parameter mesh
    R_mesh, PD_mesh = np.meshgrid(R_array, PD_array)

    # initialize object to save sweep results
    sweepresults = ParametricSweep(N)
    # sweep through parameter space, calculate min mass
    for i in range(N):
        for j in range(N):
            flowdata = Flow(R_mesh[i, j], PD_mesh[i, j], c, z, properties)
            oned_flow_modeling(flowdata)
            sweepresults.save_iteration(flowdata, i, j, N)

    return sweepresults

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("r_lower", type=float, help="channel r lower lim [m]")
    parser.add_argument("r_upper", type=float, help="channel r upper lim [m]")
    parser.add_argument("pd_lower", type=float, help="PD lower lim [m]")
    parser.add_argument("pd_upper", type=float, help="PD upper lim [m]")
    parser.add_argument("z", type=float, help="axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")
    parser.add_argument("steps", type=int, help="parameter resolution")
    parser.add_argument("-plotkey", type=str,
                        help="parameter parameter to plot")
    parser.add_argument("-i", action='store_true', dest='show',
                        default=False, help="--display plot")

    args = parser.parse_args()

    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel" +\
              "diameter! Set min PD > 1!")
        sys.exit()
    
    pc_data = PowerCycleSweep('test_params.txt')
    pc_data.load_params()

    sweepresults = sweep_geometric_configs((args.r_lower, args.r_upper),
                                           (args.pd_lower, args.pd_upper),
                                            args.z, args.clad_t, args.steps,
                                            pc_data.params[0]['bulk'])
    sweepresults.disp_min_mass()

    if args.plotkey:
        plt = plot(sweepresults, args.plotkey)
        savename = args.plotkey + '.png'
        plt.savefig(savename, dpi=500)
        if args.show:
            plt.show()
