""" Thermal-Hydraulic Optimization Module.
This module contains functions to perform a parametric sweep of reactor
geometric parameters. It calculates fuel mass and flow characteristics for each
set of valid geometric parameters.

Functions contained in this module:
    *oneD_flow_modeling
    *sweep_configs
"""
# Other imports
import argparse
import sys
# Import TH functions
from ht_functions import Flow, ParametricSweep
from plot import plot

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("d_lower", type=float, help="channel D lower lim [m]")
    parser.add_argument("d_upper", type=float, help="channel D upper lim [m]")
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
    
    sweepresults = ParametricSweep(args.steps)
    sweepresults.sweep_geometric_configs((args.d_lower, args.d_upper),
                                         (args.pd_lower, args.pd_upper),
                                          args.z, args.clad_t)
    
    if args.plotkey:
        plt = plot(sweepresults, args.plotkey, Flow.savedata)
        savename = args.plotkey + '.png'
        plt.savefig(savename, dpi=500)
        if args.show:
            plt.show()

if __name__ == '__main__':
    main()
