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
import os
import tempfile
import subprocess
import multiprocessing as mp
import glob
import itertools
from pyne.material import MaterialLibrary
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep
from scale_inputs import PinCellSCALE
from physical_constants import rho_fuel

def load_pyne_matlib():

    # Initialize material libraries.
    path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(path_to_compendium,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")
    return raw_matlib

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
    input.write_fuel_string(0.6, ('Nitrogen', 1), pyne_matlib)
    input.write_mat_string('Inconel-718', 'Carbon Dioxide', pyne_matlib)
    infilename = input.write_input()

    return infilename

def work(in_file):
    """Defines work unit on an input file"""
    subprocess.call(['batch6.1', '{}'.format(in_file)])


def run_scale():
    """
    """
    inputs = './inputs/*inp'
    tasks = glob.glob(inputs)
    # set up parallel pool
    count = mp.cpu_count()
    pool = mp.Pool(processes=count)
    # run jobs
    pool.map(work, tasks)

def save_keff(N, saveiterations, R, PD):
    
    for i in range(N):
      for j in range(N):
          # get name of output file
          r = R[i,j] * 100; pd = PD[i,j]
          ofilename = "./inputs/leakage_{0}_{1}.out".format(round(r, 5),
                                                        round(pd, 5))
          with open(ofilename, 'r') as results:
              for line in results:
                  if 'k-eff = ' in line:
                      keff = line.split()
                      saveiterations.data['keff'][i,j] =\
                                             float(keff[2])
                      break

def sweep_configs(R, PD, z, c, N, key, neutronics, save=False):
    """Perform parametric sweep through pin cell geometric space.
    """
    # calculate appropriate step sizes given range
    D_step = (R[1] - R[0]) / N
    PD_step = (PD[1] - PD[0]) / N
    # ranges for diameter and pitch/diameter ratio
    R_save = np.arange(R[0], R[1], D_step)
    PD_save = np.arange(PD[0], PD[1], PD_step)

    # create parameter mesh
    R, PD = np.meshgrid(R_save, PD_save)
    
    # initialize object to save sweep results
    sweepresults = ParametricSweep(R, PD, N)

    raw_matlib = load_pyne_matlib() 
    scale_inputs = [] 
    # sweep through parameter space, calculate min mass    
    for i in range(N):
        for j in range(N):
            r = R[i,j]; pd_ratio = PD[i,j]
            TH_configuration = THCalc(r, pd_ratio, c, z)
            sweepresults.save_iteration(TH_configuration, i, j)
            # write scale inputs
            if neutronics == True:
                scale_inputs.append(write_scale_inputs(r, pd_ratio, c,\
                    raw_matlib))
    
    if neutronics == True:
        run_scale()
        save_keff(N, sweepresults, R, PD)

    # get minimum data        
    sweepresults.get_min_data()
    plt = sweepresults.plot(R, PD, key)
    
    if save == True:
        savename = key + '.png'
        plt.savefig(savename, dpi=500)

    plt.show()


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
