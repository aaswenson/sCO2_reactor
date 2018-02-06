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
import time
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Lock, RLock
import os
from pyne.material import MaterialLibrary
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep
from scale_inputs import PinCellSCALE
from physical_constants import rho_fuel


def THCalc(r, pd_ratio, c, z):
    # perform flow calculation
    flowdata = FlowIteration(r, pd_ratio, c, z)
    flowdata.oneD_calc()
    flowdata.check_dp()
    flowdata.calc_reactor_mass()
    flowdata.calc_aspect_ratio()

    return flowdata

def write_scale_inputs(r, pd_ratio, c, pyne_matlib):
    """
    """
    # write scale input file
    input = PinCellSCALE(r, pd_ratio, c)
    input.write_fuel_string(0.9, ('Nitrogen', 1), pyne_matlib)
    input.write_mat_string('Inconel-718', 'Carbon Dioxide', pyne_matlib)
    infilename = input.write_input()

def run_scale(R, PD, i, j):
    """
    """
    r = R[i,j] * 100
    pitch = PD[i,j]
    ifilename = "./inputs/leakage_{0}_{1}.inp".format(round(r, 5),
                                                      round(pitch, 5))
    ofilename = "./inputs/leakage_{0}_{1}".format(round(r, 5),
                                                      round(pitch, 5))
    command = "batch6.1 " + ifilename + " -o " + ofilename

    os.system(command)

def save_keff(N, saveiterations):
    for i in range(N):
        for j in range(N):

def sweep_configs(R, PD, z, c, N, key, neutronics, save=False):
    """Perform parametric sweep through pin cell geometric space.
    """
    # calculate appropriate step sizes given range
    D_step = (R[1] - R[0]) / N
    PD_step = (PD[1] - PD[0]) / N
    # ranges for diameter and pitch/diameter ratio
    R = np.arange(R[0], R[1], D_step)
    PD = np.arange(PD[0], PD[1], PD_step)

    # create parameter mesh
    R, PD = np.meshgrid(R, PD)
    
    # initialize object to save sweep results
    sweepresults = ParametricSweep(R, PD, N)

    
    # Initialize material libraries.
    path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(path_to_compendium,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")
    
    # sweep through parameter space, calculate min mass    
    for i in range(N):
        for j in range(N):
            r = R[i,j]; pd_ratio = PD[i,j]
            TH_configuration = THCalc(r, pd_ratio, c, z)
            if neutronics == True:
                write_scale_inputs(r, pd_ratio, c, raw_matlib)

    if neutronics == True:
        pool_size = 8
        pool = Pool(pool_size) 
        for i in range(N):
            for j in range(N):
                run_scale(R, PD, i, j)
#               pool.apply_async(run_scale, (R, PD, i, j, ))
#        pool.close()
#        pool.join()

    sweepresults.save_iteration(TH_configuration, i, j)

    # get minimum data        
    sweepresults.get_min_data()
    plt = sweepresults.plot(R, PD, key)
    
    if save == True:
        savename = key + '.png'
        plt.savefig(savename, dpi=500)

#    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("r_lower", type=float, help="channel r lower lim [m]")
    parser.add_argument("r_upper", type=float, help="channel r upper lim [m]")
    parser.add_argument("pd_lower", type=float, help="PD lower lim [m]")
    parser.add_argument("pd_upper", type=float, help="PD upper lim [m]")
    parser.add_argument("z", type=float, help="axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")
    parser.add_argument("steps", type=int, help="parameter resolution")
    parser.add_argument("plotkey", type=str, help="parameter parameter to plot")
    parser.add_argument("--n", action='store_true', dest='neutronics',
            default=False, help="--write scale inputs and run reactor physics\
calc")

    args = parser.parse_args()
    
    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel\
diameter! Set min PD > 1!")
        sys.exit()

    sweep_configs((args.r_lower, args.r_upper),
                  (args.pd_lower, args.pd_upper), 
                   args.z, args.clad_t, args.steps,
                   args.plotkey, args.neutronics, True)
