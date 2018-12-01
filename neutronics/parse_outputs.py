import operator
import math
import material_data as md
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import glob
import neutronic_sweeps as ns
import pandas

names = ns.dimensions + ['keff', 'r_mass', 'c_mass', 'p_mass', 'mass']
types = ['f8']*len(names)

def load_outputs(data_dir):
    """Load the MCNP output file
    """
    file_strings = []
    files = glob.glob(data_dir)
    for output in files:
        # load file and read the lines
        fp = open(output, 'r')
        file_strings.append(fp.readlines())
    
    return file_strings

def parse_keff(lines):
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

    return keff, stdv

def calc_fuel_mass(data):
    """
    """

    return data

def parse_header_str(lines, data):
    """
    """
    data['fuel

def save_store_data(data_dir='/mnt/sdb/calculation_results/sa_results/*.i.o'):
    """
    """
    files = glob.glob(data_dir)
    N = len(files)
    data = np.zeros(N, dtype={'names' : names, 'formats' : types})

    for idx, file in enumerate(files):
        print(file)
        fp = open(file, 'r')
        string = fp.readlines()
        fp.close()
        params = parse_header_string(string)
        


    np.savetxt("depl_results.csv", data, delimiter=',', 
           fmt='%10.5f', header=','.join(names))
  
def plot_results(data, ind, dep, colorplot=None, log=None):
    """Generate Plots
    """
    label_strings = {'core_r' : 'Core Radius [cm]',
                     'enrich' : 'U-235 Enrichment [-]',
                     'keff' : 'k-eff [-]',
                     'mass' : 'reactor fuel mass [kg]',
                     'fuel_frac' : 'volume fraction fuel [-]',
                    }
    # plot
    fig = plt.figure()
    if colorplot:
        colorsave = '_'+colorplot
        plt.scatter(data[ind], data[dep], c=data[colorplot], s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(data[colorplot]))))
        plt.colorbar(label=label_strings[colorplot])
    else:
        colorsave = ''
        plt.scatter(data[ind], data[dep], s=6)
    # titles and labels
    plt.title("{0} vs. {1}".format(dep, ind))
#    plt.title(r'EOL k$_{eff}$ vs. Fuel Mass')
#    plt.title("keff vs. mass for 0.2 < enrich < 0.3")
    plt.xlabel(label_strings[ind])
    plt.ylabel(label_strings[dep])

    if log == 'log-log':
        plt.xscale('log')
        plt.yscale('log')
    if log == 'semilogy':
        plt.yscale('log')
    if log == 'semilogx':
        plt.xscale('log')

    savename = '{0}_vs_{1}{2}.png'.format(dep, ind, colorsave)
    plt.savefig(savename, dpi=1000, format='png')

    return plt


def surf_plot(data):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X = data['core_r']
    Y = data['fuel_frac']
    Z = data['ref_mult']
    k = data['keff']

    # Plot the surface.
    ax.scatter(X,Y,Z, c=Z)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    plt.show()

def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = pandas.read_csv(datafile)
    
    return data

def filter_data(filters, data):

    """Apply useful filters on the data
    """
    opers = {'<' : operator.lt,
             '=' : operator.eq,
             '>' : operator.gt}
    
    for filter in filters:
        op= filter.split()[1]
        key = filter.split()[0]
        val = float(filter.split()[2])
        data = data[opers[op](data[key], val)]
    
    return data

if __name__ == '__main__':
    save_store_data()
    data = load_from_csv()
#    data = filter_data(filter, data)
    surf_plot(data)
#    plt = plot_results(data, 'mass', 'keff')
