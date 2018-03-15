# Matplotlib imports
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.axis
from matplotlib import cm, rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter
import os
import glob
import zipfile
import numpy as np

def parse_table_210(fp, n_steps):
    """Parse print table 210 to get depletion results.
    """
    keff = []
    time = []
    burnup = []
    
    for i in range(7):
        next(fp)
    
    for entry in range(n_steps):
        line = fp.readline()
        time.append(float(line.split()[2]))
        burnup.append(float(line.split()[8]))
        keff.append(float(line.split()[4]))
    
    return time, burnup, keff

def parse_output(filename, n_steps):
    """Parse MCNP output file for burnup results.
    """
    t = 0
    # read output file
    fp = open(filename, 'r')
    for line in fp:
        if "print table 210" in line:
            t, BU, k = parse_table_210(fp, n_steps)
    fp.close()

    return t, BU, k


def unzip_collect_outputs():
    """Extract output files and collect their info
    """
    if os.path.isdir('./results') == False:
        zip_ref = zipfile.ZipFile('depl_results.zip', 'r') 
        zip_ref.extractall("results")
    files = glob.glob("./results/*")
    
    return files

def process_filename(file):
    """Parse an output filename for core radius and fuel fraction.
    """
    name = file.split('/')[-1]
    frac = float(name.split('_')[1])
    r = float(name.split('_')[2].split('.')[0])
    
    return frac, r

def EOL_keff(outputs):
    """For file in outputs
    """
    N = len(outputs)
    data = np.zeros(N, dtype={'names' : ['r', 'frac', 'EOL_k', 'grad'],
                              'formats' : ['f8', 'f8', 'f8', 'f8']})
    results = []
    for idx, file in enumerate(outputs):
        frac, r = process_filename(file)
        t, BU, k = parse_output(file, 30)
        data[idx]['r'] = r
        data[idx]['frac'] = frac
        data[idx]['EOL_k'] = k[-1]
    data['grad'] = np.gradient(data['EOL_k'])
        
    return data

def plot_grad_keff(data):

    Nx = len(set(data['r']))
    Ny = len(set(data['frac']))
    
    # setup the figure
    fig, ax2 = plt.subplots(1,1, sharex=True, sharey=True)
    # plot keff
#    cs1 = ax1.tricontourf(data['r'],data['frac'], grad, 20)
#    ax1.plot(data['r'],data['frac'], 'ko ')
    # plot grad keff
    cs2 = ax2.tricontourf(data['r'], data['frac'], data['EOL_k'], 20)
    ax2.plot(data['r'],data['frac'], 'ko ')
#    ax1.set_title("Gradient Map")
    # Add a color bar which maps values to colors.
#    fig.colorbar(cs1, ax=ax1)
    fig.colorbar(cs2, ax=ax2)
#    ax1.set_title("Gradient Map")
    ax2.set_title("EOL Keff")
    plt.show()

def plot_contour(data):
    
    fig, ax = plt.subplots(1,1)
    cs = ax.tricontourf(data['r'], data['frac'], data['EOL_k'], 10)
    ax.scatter(data['r'],data['frac'], s=12, facecolors='none', edgecolors='k')
    fig.colorbar(cs)
    plt.xlabel("Core Radius [cm]")
    plt.ylabel("Fuel Fraction [cm]")
    plt.title("EOL Reactivity - 10y @ 99 kW")
    plt.show()

def plot_slice(slice_choice, data):
    
    slice_datakey = slice_choice[0]
    slice_val = slice_choice[1]
    ind_datakey = slice_choice[2]
    
    X = []
    Y = []

    for idx, item in enumerate(data[slice_datakey]):
        if item == slice_val:
            X.append(data[ind_datakey][idx])
            Y.append(data['EOL_k'][idx])
    
    
    fig = plt.figure()
    plt.scatter(X, Y)
    xlabels = {'r' : 'Fuel Fraction [-]',
               'frac' : 'Core Radius [cm]'
              }
    titles = {'r':'EOL Reactivity vs. Fuel Fraction at Core Radius: {0} [cm]',
              'frac': 'EOL Reactivity vs. Core Radius at Fuel Frac: {0}'
             }
    title_str = titles[slice_datakey].format(slice_val)
    plt.title(title_str)
    plt.xlabel(xlabels[slice_datakey])
    plt.ylabel('EOL Reactivity [-]')
    plt.show()


testfiles = unzip_collect_outputs()
data = EOL_keff(testfiles)

plot_slice(('frac', 0.5, 'r'), data)

