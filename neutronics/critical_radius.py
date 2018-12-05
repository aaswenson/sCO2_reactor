import sys
import os
import numpy as np
from scipy.optimize import minimize_scalar, curve_fit
from subprocess import call, DEVNULL
from mcnp_inputs import HomogeneousInput
import fit_data as fd
import parse_outputs as po

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

target_keff = 1.01
# critical radius range to sweep
domain = (5, 28.75)

g_to_kg = 0.001


fig = plt.figure()

def crit_radius(config, func):
    """Get critical radius = f(fuel_frac)
    """
    f = config['fuel_frac']
    m = config['ref_mult']
    interp_keff = lambda r: (func([r, f, m])[0] - target_keff)**2
    res = minimize_scalar(interp_keff, args=(), 
                          method='bounded', bounds=(10, 50))
    
    config['keff_err'] = res.fun
    config['core_r'] = res.x
    
def fuel_frac(coolant, fuel, clad, matr, func):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    rhos = {'CO2' : 252.638e-3, 'H2O' : 141.236e-3}
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : coolant,
              'clad' : clad,
              'rho_cool' : rhos[coolant],
             }
    resname = '{0}_{1}_results.txt'.format(coolant, fuel)
    
    resfile = open(resname, '+w')
    resfile.write('fuel_frac,crit_radius\n') 
    resfile.close()
    
    results = {'frac'   : [],
               'mult'   : [],
               'radius' : []
              }

    for frac in np.linspace(0.3, 0.95, 200):
        resfile = open(resname, 'a')
        config['fuel_frac'] = frac
        config['ref_mult'] = refl_mult(config, func)
        # get critical radius
        crit_radius(config, func)
        
        #save data for plotting
        results['frac'].append(config['fuel_frac'])
        results['mult'].append(config['ref_mult'])
        results['radius'].append(config['core_r'])

        resfile.write('{0:.2f},{1:.5f},{2:.5f}\n'.format(config['fuel_frac'],
                                                         config['core_r'],
                                                         config['ref_mult']))
        resfile.close()
        
    # plot mass mult curves
    plt.legend()
    plt.savefig('mass_vs_mult.png')

    # plot results
    fig = plt.figure()
    plt.scatter(results['frac'], results['mult'])
    for f in func.grid[1]:
        if f < 0.3:
            continue
        plt.axvline(x=f)
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Optimized Reflector Multiplier [-]')
    plt.title('Optimal Reflector Multipliers')
    plt.savefig('refl_mult.png')
    # core radius
    fig = plt.figure()
    plt.scatter(results['frac'], results['radius'])
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Critical Core Radius [cm]')
    plt.title('Core Radius for keff = 1.01')
    plt.savefig('core_r.png')

def refl_mult(config, func):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """

    mults = np.linspace(0.001, 0.4, 50)
    data = {'mass' : [], 'r' : [], 'mult' : [], 'keff' : []}
    refl_res = open('refl_results.txt', 'a')
     
    for mult in mults:
        config['ref_mult'] = mult
        # get critical radius
        crit_radius(config, func)
        if config['keff_err'] > 1e-5:
            continue
        input = HomogeneousInput(config=config)
        input.homog_core()
        
        data['mass'].append(input.tot_mass/ 1000)
        data['mult'].append(mult)
    
    poly = np.polyfit(data['mult'], data['mass'], 3)
    fit_mults = np.linspace(data['mult'][0], data['mult'][-1], 1000)
    fit_mass = np.polyval(poly, fit_mults)
   

    massfunc = lambda m: np.polyval(poly, m)
    opt_mult = minimize_scalar(massfunc, method='bounded', 
                               bounds=(data['mult'][0],
                                       data['mult'][-1])).x
###############################################################


    plt.scatter(data['mult'], data['mass'], s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(data['r']))))
    
    plt.plot(fit_mults, fit_mass, label='{0:.2f}'.format(config['fuel_frac']))

    plt.title('Fuel Frac: {0}'.format(config['fuel_frac']))
    plt.xlabel('reflector mult [-]')
    plt.ylabel('reactor mass [kg]')
    
    return opt_mult

if __name__ == '__main__':
    data = po.load_from_csv('./crit_results.csv')
    fn = po.interpolate_grid(data)
    fuel_frac('CO2', 'UO2', 'Inconel-718', None,fn)
