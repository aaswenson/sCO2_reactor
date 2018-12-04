import glob
import numpy as np
import argparse
from scipy.optimize import curve_fit, minimize_scalar

cm_to_m = 0.01

def min_mult(params, func):
    coeffs = tuple(params)
    res = minimize_scalar(func, args=coeffs, method='bounded', bounds=(0.001,0.6))
    
    return res.x

def load_data(file):
    """Load the data from text file
    """

    data = pandas.read_csv(file)

    return data

def load_refl(file):
    """Load the reflector data.
    """
    lines = open(file, 'r').readlines()
    mult = []
    mass = []
    for line in lines[1:]:
        data = [float(x) for x in line.split(',')]
        mult.append(data[0])
        mass.append(data[2])

    return {'mult' : mult, 'mass' : mass}
    
def cubic(x, a, b, c, d):
    """
    """
    A = a*np.power(x,3)
    B = b*np.power(x,2)
    C = c*x

    return np.add(A,np.add(B,np.add(C,d)))

def power(x, a, b):
    """Exponential fitting function
    """

    return a*np.power(x, b)

def poly(x, a, b, c):

    return np.add(np.add(a*np.power(x,2), np.multiply(x, b)), c)

def fit_data(data, func, x, y):
    """Fit the data to function
    """    
    popt, pcov = curve_fit(func, data[x], data[y])

    return popt, pcov

def plot_results(data):
    """
    """
    fig, ax = plt.subplots()


    line_formats = {'CO2 UO2 ' : ('r', 'o'),
                    'H2O UO2 ' : ('r', 'v'),
                    'CO2 UN '  : ('b', 'o'),
                    'H2O UN '  : ('b', 'v')}

    for rxtr in data:
        res = data[rxtr]
        ax.scatter(res[0], res[1], 
                   c=line_formats[rxtr][0],
                   marker=line_formats[rxtr][1], 
                   label=rxtr)

    plt.legend()

    plt.title('Critical Radius: Buried Reactor on Mars')
    plt.xlabel('Fuel Fraction [-]')
    plt.ylabel('Critical Core Radius [cm]')
    plt.savefig('crit_radius.png', dpi=700)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", '--reflector', action='store_true', default=False, help="reflector mode")
    parser.add_argument("-c", '--crit_radius', action='store_true', default=False,
            help="radius mode")
    args = parser.parse_args()

    
    res_files = glob.glob('./*results.txt')
    plot_data = {}
    for file in res_files:
        name = " ".join(file.split('results.txt')[0].split('_')).strip('./')
        if args.reflector:
            data = load_refl(file)
            popt, cov = fit_data(data, cubic, 'mult', 'mass')
            x = min_mult(popt)
            print(name + ' {0:.4f}'.format(x))

        elif args.crit_radius:
            data = load_data(file)
            plot_data[name] = (data['fuel_frac'], data['crit_radius'])
            popt, cov = fit_data(data, power, 'fuel_frac', 'crit_radius')
            print(name + ' ' + str(popt[0]) + ' ' + str(popt[1]))

    if args.crit_radius:
        plot_results(plot_data)
