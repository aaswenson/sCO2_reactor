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
import tempfile
import subprocess
import multiprocessing as mp
import glob
import itertools
from pyne.material import MaterialLibrary
# Import TH functions
from ht_functions import Flow, ParametricSweep, flow_calc
from plot import plot
from mcnp_inputs import PinCellMCNP
from scale_inputs import PinCellSCALE
from physical_constants import rho_fuel

models = ['mcnp', 'scale']
model_objects = {'mcnp': PinCellMCNP, 'scale' : PinCellSCALE}

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

def write_inputs(r, pd_ratio, c, pyne_matlib):
    """
    """
    # write scale input file
    input = PinCellInput(r, pd_ratio, c)
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

def save_keff_scale(N, saveiterations, R, PD):
    
    for i in range(N):
      for j in range(N):
          # get name of output file
          r = R[i,j]*100; pd = PD[i,j]
          ofilename = "./inputs/leakage_{0}_{1}.out".format(round(r, 5),
                                                        round(pd, 5))
          with open(ofilename, 'r') as results:
              for line in results:
                  if 'k-eff = ' in line:
                      keff = line.split()
                      saveiterations.data['keff'][i,j] =\
                                             float(keff[2])
                      break

def sweep_configs(radii, pds, z, c, N, neutronics, run_neutronics):
    """Perform parametric sweep through pin cell geometric space.
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
    sweepresults = ParametricSweep(R_mesh, PD_mesh, N)
    
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
            # get geometry
            r = R_mesh[i,j]
            pd_ratio = PD_mesh[i,j]
            
            flowdata = Flow(r, pd_ratio, c, z)
            flow_calc(flowdata)
            sweepresults.save_iteration(flowdata, i, j)
            
    
    if neutronics == True:
                run_scale(R, PD, i, j)

    if run_neutronics == True:
        run_scale()
        save_keff_scale(N, sweepresults, R, PD)

    # get minimum data        
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
    parser.add_argument("-n", action='store_true', dest='neutronics',
            default=False, help="write scale inputs and run reactor physics\
calc")
    parser.add_argument("-r", action='store_true', dest='run_neutronics',
            default=False, help="run reactor_physics calc")
    parser.add_argument("--save", action='store_true', dest='saveplot',
            default=False, help="--save, save plot as png")
    # select the version of the pin cell
    parser.add_argument("-m", "--model", type=str, default=models[0], required=False, choices=models,
                        help="Allows selection for which model to generate an input file. ")

    args = parser.parse_args()

    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel" +\
              "diameter! Set min PD > 1!")
        sys.exit()
    
    PinCellInput = model_objects[args.model]

    sweepresults = sweep_configs((args.r_lower, args.r_upper),
                                 (args.pd_lower, args.pd_upper),
                                  args.z, args.clad_t, args.steps,
                                  args.neutronics, args.run_neutronics)
    
    if args.plotkey:
        plt = plot(sweepresults, args.plotkey)
        savename = args.plotkey + '.png'
        plt.savefig(savename, dpi=500)
        if args.show:
            plt.show()
