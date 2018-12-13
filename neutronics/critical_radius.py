import sys
import os
import numpy as np
from scipy.optimize import minimize_scalar, curve_fit
import matplotlib.pyplot as plt

from mcnp_inputs import HomogeneousInput
import parse_outputs as po

target_keff = 1.01

fig = plt.figure()
ax = plt.subplot(111)

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
    """returns distance from value to nearest point on the grid in the dimension
    desired. Burden is on the user to ensure the dimension of the value and the
    array match.
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    
    return abs(array[idx] - value)

def calc_rxtr_mass(config):
    """Calculate the mass of a given reactor configuration
    """
    input = HomogeneousInput(config=config)
    input.homog_core()
    # g -> kg
    mass = input.tot_mass / 1000
    
    return mass

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
    
def fuel_frac(config, func):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """
    # prepare the results file
    resname = '{0}_{1}_results.txt'.format(config['cool'], config['fuel'])
    resfile = open(resname, '+w')
    resfile.write('fuel_frac,crit_radius,ref_mult,mass_std\n') 
    resfile.close()
    
    results = {'frac'   : [],
               'mult'   : [],
               'radius' : [],
               'mass'   : [],
               'r_radius' : [],
               'r_mult' : [],
              }

    for frac in np.linspace(0.3, 0.95, 20):
        resfile = open(resname, 'a')
        config['fuel_frac'] = frac
        config['ref_mult'], mstd = refl_mult(config, func)
        # get critical radius
        crit_radius(config, func)
        
        #save data for plotting
        results['mass'].append(calc_rxtr_mass(config))
        results['frac'].append(config['fuel_frac'])
        results['mult'].append(config['ref_mult'])
        results['radius'].append(config['core_r'])
        results['r_radius'].append(find_nearest(func.grid[0], config['core_r']))
        results['r_mult'].append(find_nearest(func.grid[2], config['ref_mult']))
        resfile.write('{0:.2f},{1:.5f},{2:.5f},{3:.5f}\n'.\
        format(config['fuel_frac'],config['core_r'],config['ref_mult'], mstd))
        
        resfile.close()
        
    # plot mass mult curves
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),ncol=5,
              fancybox=True)
    for m in results['mult']:
        ax.scatter(results['mult'], results['mass'],c='r', marker='X')
    plt.title('Reactor Mass Dependence on Reflector Multiplier')
    plt.xlabel('reflector mult [-]')
    plt.ylabel('normed reactor mass [-]')
#    plt.ylim((0.0495, 0.0507))
    plt.savefig('mass_vs_mult.png', dpi=500, figsize=(20,10))
    
    # plot results
    fig = plt.figure()
    plt.scatter(results['frac'], results['mult'], #c=results['r_mult'],
                cmap=plt.cm.get_cmap('Reds', len(set(results['r_mult']))))
#    plt.colorbar(label='Distance from nearest grid point [cm]')

    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Optimized Reflector Multiplier [-]')
    plt.title('Optimal Reflector Multipliers')
    plt.savefig('refl_mult.png')
    # core radius
    fig = plt.figure()
    plt.plot(results['frac'], results['radius'])
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Critical Core Radius [cm]')
    plt.title('Core Radius for keff = 1.01, m=0.15')
    plt.savefig('core_r.png')
    
def refl_mult(config, func):
    """Determine the optimal reflector thickness for a given reactor
    configuration.
    """

    mults = np.linspace(0.001, max(func.grid[2]), 20)
    data = {'mass' : [], 'r' : [], 'mult' : [], 'keff' : []}
     
    for mult in mults:
        config['ref_mult'] = mult
        # get critical radius
        crit_radius(config, func)
        if config['keff_err'] > 1e-5:
            print('Warning, keff target not reached for r= {0}, f= {1}'.format(
                    config['core_r'], config['fuel_frac']))
            continue
        
        data['mass'].append(calc_rxtr_mass(config))
        data['mult'].append(mult)
    
    # normalize data to 1
#    data['mass'] = [x / sum(data['mass']) for x in data['mass']]
    
    poly = np.polyfit(data['mult'], data['mass'], 2)
    fit_mults = np.linspace(data['mult'][0], data['mult'][-1], 1000)
    fit_mass = np.polyval(poly, fit_mults)
    
    massfunc = lambda m: np.polyval(poly, m)
    opt_mult = minimize_scalar(massfunc, method='bounded', 
                               bounds=(data['mult'][0],
                                       data['mult'][-1])).x
###############################################################
    if config['fuel_frac']:# == func.grid[1][-2]:
        ax.scatter(data['mult'], data['mass'], s=6)
        ax.plot(fit_mults, fit_mass, label='{0:.2f}'.format(config['fuel_frac']))

    return opt_mult, np.std(data['mass'])

if __name__ == '__main__':
    
    args = sys.argv[1:]
    matr = None
    if len(args) == 5:
        matr = args[4]
    config = {'fuel' : args[1],
              'matr' : matr,
              'cool' : args[2],
              'clad' : args[3],
             }
    
    data = po.load_from_csv(args[0])
    fn = po.interpolate_grid(data)
    fuel_frac(config, fn)
