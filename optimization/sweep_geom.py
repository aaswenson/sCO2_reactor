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
import os
from pyne.material import MaterialLibrary
# Import TH functions
from ht_functions import Flow, ParametricSweep, flow_calc
from plot import plot
from mcnp_inputs import PinCellMCNP
from physical_constants import rho_fuel
from scale_inputs import PinCellSCALE



        
def sweep_configs(R, PD, z, c, N, key, AR_select, save=False):
    """Perform parametric sweep through pin cell geometric space.

    """
    # calculate appropriate step sizes given range
    D_step = (diams[1] - diams[0]) / N
    PD_step = (pds[1] - pds[0]) / N
    # ranges for diameter and pitch/diameter ratio
    D_array = np.arange(diams[0], diams[1], D_step)
    PD_array = np.arange(pds[0], pds[1], PD_step)

    # create parameter mesh
    D_mesh, PD_mesh = np.meshgrid(D_array, PD_array)

    # initialize object to save sweep results
    sweepresults = ParametricSweep(D_mesh, PD_mesh, N)
    
    # Initialize material libraries.
    path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(path_to_compendium,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")
    pool_size = 8
    pool = Pool(pool_size)
    
    t0 = time.time()

    # sweep through parameter space, calculate min mass    
    for i in range(N):
        for j in range(N):
            flowdata = Flow(D_mesh[i, j], PD_mesh[i, j], c, z)
            oneD_flow_modeling(flowdata)
            sweepresults.save_iteration(flowdata, i, j)
            # write MCNP input file
            
            infile = PinCellMCNP('hex',D[i,j]*50 , PD[i,j], c * 100, z *
                    100)
            infile.write_fuel_string(0.96, ('Nitrogen', 1), raw_matlib)
            infile.write_mat_string('Steel, Stainless 316', 'Carbon Dioxide',
                raw_matlib)
            infile.write_input([5000, 25, 40])
    
    sweepresults.get_min_data()
    sweepresults.disp_min_mass()
    
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
    

    sweepresults = sweep_configs((args.d_lower, args.d_upper),
                                 (args.pd_lower, args.pd_upper),
                                  args.z, args.clad_t, args.steps)
    
    if args.plotkey:
        plt = plot(sweepresults, args.plotkey)
        savename = args.plotkey + '.png'
        plt.savefig(savename, dpi=500)
        if args.show:
            plt.show()
