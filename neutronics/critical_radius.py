import sys
import os
import numpy as np
from scipy.optimize import minimize_scalar, curve_fit
from subprocess import call, DEVNULL
from mcnp_inputs import HomogeneousInput
import fit_data as fd
import parse_outputs as po

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

target_keff = 1.01
# critical radius range to sweep
domain = (5, 28.75)

g_to_kg = 0.001


fig = plt.figure()

def calc_keff(config):
    """Calculate keff deviation from target.
    """
    frac = config['fuel_frac']
    radius = config['core_r']

    basename = "{0}_{1}.i".format(round(frac,5), round(radius,5))
    write_inp(basename, config)
    call(["mcnp6", "n= {0} tasks 8".format(basename)], stdout=DEVNULL)
    keff = parse_output(basename)
    os.remove('{0}r'.format(basename))
    os.remove('{0}s'.format(basename))
    os.remove('{0}o'.format(basename))
    os.remove(basename)
    
    print(frac, keff)

def parse_output(basename):
    """Parse the output for keff value.
    """

    fp = open(basename + 'o', 'r')
    lines = fp.readlines()
    fp.close()

    res_idx = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_idx.append(idx)
    keff = float(lines[res_idx[-1]].split()[2])
    stdv = float(lines[res_idx[-1]].split()[3])

    return keff

def write_inp(basename, configuration):
    """Write mcnp input for keff.
    """
    input = HomogeneousInput(config=configuration)
    homog_comp = input.homog_core()
    input.write_input(basename)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    return abs(array[idx] - value)

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
               'radius' : [],
               'r_radius' : [],
               'r_mult' : [],
               'p1' : [], 'p2' : [], 'p3' : []
              }

    for frac in np.linspace(0.2, 0.95, 100):
#    for frac in func.grid[1]:
        resfile = open(resname, 'a')
        config['fuel_frac'] = frac
        config['ref_mult'], popt= refl_mult(config, func)
        # get critical radius
        crit_radius(config, func)
        
        #save data for plotting
        results['frac'].append(config['fuel_frac'])
        results['p1'].append(popt[0])
        results['p2'].append(popt[1])
        results['p3'].append(popt[2])
        results['mult'].append(config['ref_mult'])
        results['radius'].append(config['core_r'])
#        results['r_radius'].append(find_nearest(func.grid[0], config['core_r']))
#        results['r_mult'].append(find_nearest(func.grid[2], config['ref_mult']))
        resfile.write('{0:.2f},{1:.5f},{2:.5f}\n'.format(config['fuel_frac'],
                                                         config['core_r'],
                                                         config['ref_mult']))
        resfile.close()
        
    # plot mass mult curves
    plt.legend()
    plt.savefig('mass_vs_mult.png', dpi=1000, figsize=(20,10))

    # plot results
    fig = plt.figure()
    plt.scatter(results['frac'], results['mult'], #c=results['r_radius'],
                cmap=plt.cm.get_cmap('Reds', len(set(results['r_radius']))))
#    plt.colorbar(label='Distance from nearest grid point [cm]')

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
    
    # core radius
    fig = plt.figure()
    plt.scatter(results['frac'], results['p1'])
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('p1')
    plt.title('Core Radius for keff = 1.01')
    plt.savefig('p1.png')
    
    # core radius
    fig = plt.figure()
    plt.scatter(results['frac'], results['p2'])
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('p2')
    plt.title('Core Radius for keff = 1.01')
    plt.savefig('p2.png')
    
    # core radius
    fig = plt.figure()
    plt.scatter(results['frac'], results['p3'])
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('p3')
    plt.title('Core Radius for keff = 1.01')
    plt.savefig('p3.png')

def refl_mult(config, func):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """

    mults = np.linspace(0.001, 0.15, 10)
    data = {'mass' : [], 'r' : [], 'mult' : [], 'keff' : []}
    refl_res = open('refl_results.txt', 'a')
     
#    for mult in func.grid[2]:
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

#    data['mass'] = [x / sum(data['mass']) for x in data['mass']]
    
    poly = np.polyfit(data['mult'], data['mass'], 3)
    fit_mults = np.linspace(data['mult'][0], data['mult'][-1], 1000)
    fit_mass = np.polyval(poly, fit_mults)
    
#    print('{0:.3f} {1:.3f}'.format(config['fuel_frac'],
#                           np.std(data['mass'])))
    
    massfunc = lambda m: np.polyval(poly, m)
    opt_mult = minimize_scalar(massfunc, method='bounded', 
                               bounds=(data['mult'][0],
                                       data['mult'][-1])).x
###############################################################

    if config['fuel_frac']:# == func.grid[1][-2]:
#        calc_keff(config)
        plt.scatter(data['mult'], data['mass'], s=6)
        plt.plot(fit_mults, fit_mass, label='{0:.2f}'.format(config['fuel_frac']))

    plt.title('Reactor Mass Dependence on Reflector Multiplier')
    plt.xlabel('reflector mult [-]')
    plt.ylabel('reactor mass [kg]')
#    plt.ylim((0.029, 0.0321)) 
#    opt_mult = data['mult'][data['mass'].index(min(data['mass']))]

    return opt_mult, poly

def plot_hist(data):
    """
    """
    fig = plt.figure()
    plt.title('keff stdv')
    plt.hist(data['stdv'])
    plt.savefig('hist_err.png')

def plot_err(data):
    """
    """

    fig = plt.figure()
    plt.title('keff stdv vs. fuel frac')
    plt.xlabel('fuel frac [-]')
    plt.ylabel('keff stdv [-]')
    plt.ylim((0, 0.001))
    plt.scatter(data['fuel_frac'], data['stdv'], c=data['mass'], s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(data['mass']))))
    plt.colorbar(label='Core Mass [g]')

    plt.savefig('err_v_frac.png')

if __name__ == '__main__':
    data = po.load_from_csv('./crit_results.csv')
#    plot_hist(data)
#    plot_err(data)
    fn = po.interpolate_grid(data)
    fuel_frac('CO2', 'UO2', 'Inconel-718', None,fn)
