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
import os
from pyne.material import MaterialLibrary
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep
from scale_inputs import PinCellSCALE
from physical_constants import rho_fuel


def worker(R, PD, c, z, i, j, sweepresults, pyne_matlib):
    r = R[i,j]; pd_ratio = PD[i,j]
    # perform flow calculation
    flowdata = FlowIteration(r, pd_ratio, c, z)
    oneD_flow_modeling(flowdata)
    # write scale input file
    input = PinCellSCALE(r, pd_ratio, c)
    input.write_fuel_string(0.9, ('Nitrogen', 1), pyne_matlib)
    input.write_mat_string('Inconel-718', 'Carbon Dioxide', pyne_matlib)
    infilename = input.write_input()
    # run scale
    command = "batch6.1 " + infilename
    os.system(command)
    
    resultsfile = infilename.split('.inp')[0] + ".out"
    
    with open(resultsfile, 'r') as results:
        for line in results:
            if 'k-eff = ' in line:
                keff = line.split()
                flowdata.keff = float(keff[2])
                break

    sweepresults.save_iteration(flowdata, i, j)

def oneD_flow_modeling(analyze_flow):
    """Conduct oneD_flow_modeling.

    Arguments:
    ----------
        analyze_flow: (class) FlowIteration object. Contains attributes and
        methods required to perform an N_channels calculation for a single
        geometry (r, PD, L, c)

    Returns:
    --------
        None
    """
    analyze_flow.oneD_calc()
    analyze_flow.check_dp()
    analyze_flow.calc_reactor_mass()
    analyze_flow.calc_aspect_ratio()
        
def sweep_configs(R, PD, z, c, N, key, AR_select, save=False):
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
    sweepresults = ParametricSweep(R, PD, N, AR_select)

    
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
            pool.apply_async(worker, (R, PD, c, z, i, j, sweepresults,
                raw_matlib, ))
    pool.close()
    pool.join()

    print(time.time() - t0)

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
    parser.add_argument("--AR", action='store_true', dest='AR_select',
            default=False, help="--selects data corresponding to valid AR's")

    args = parser.parse_args()
    
    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel\
diameter! Set min PD > 1!")
        sys.exit()

    sweep_configs((args.r_lower, args.r_upper),
                  (args.pd_lower, args.pd_upper), 
                   args.z, args.clad_t, args.steps,
                   args.plotkey, args.AR_select, True)
