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

names = ns.dimensions + ['keff', 'ave_E', 'mass', 'q_dens', 'dU', 'dK',
                         'mass_235', 'm_dens', 'fuel_frac']
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
    """Parse the keff data from the output file.
    """
    keff = []
    err = []
    days = []
    BU = []
    
    res_loc = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_loc.append(idx)
        if 'print table 210' in line:
            burndx = idx + 8

    # skip the predictor calcs 
    save_res = res_loc[0::2]

    for line_num in save_res:
        keff.append(float(lines[line_num].split()[2]))
        err.append(float(lines[line_num].split()[3]))

    for burndata in lines[burndx:]:
        if burndata == '\n':
            break
        BU.append(float(burndata.split()[8]))
        days.append(float(burndata.split()[2]))

    return (days, BU, keff, err)

def parse_etal(tally, lines):
    """Parse energy tallies from the output file.
    """
    bins = []
    vals = []
    errs = []
    tally_locations = []
    # get number of energy bins
    bupper = lines.index(' energy bins\n')
    blower = lines.index('      total bin\n')
    nbins = blower - bupper
    
    tally_num = '{0}tally'.format(tally)

    for idx, line in enumerate(lines):
        if tally_num in line and 'nps' in line:
            tally_locations.append(idx + 11)
            
    for tally in tally_locations:
        bindata = []
        valdata = []
        errdata = []
        for idx in range(tally, tally + nbins - 1):
            bindata.append(float(lines[idx].split()[0]))
            valdata.append(float(lines[idx].split()[1]))
            errdata.append(float(lines[idx].split()[2]))
        
        bins.append(bindata)
        vals.append(valdata)
        errs.append(errdata)

    average = np.average(bins, weights=vals)

    return (bins, vals, errs, average)

def parse_header_string(string):
    """Get important parameters from MCNP6 input header string.
    """
    for line in string:
        if '1-' in line:
            data = line.split()[1]
            AR = float(data.split(',')[0])
            core_r = float(data.split(',')[1])
            cool_r = float(data.split(',')[2])
            PD = float(data.split(',')[3])
            power = float(data.split(',')[4])
            enrich = float(data.split(',')[5])
            
            break

    return [AR, PD, cool_r, core_r, enrich, power]

def calc_fuel_mass(core_r, r, PD, Q):
    """
    """
    c = 0.0031
    l = core_r
    core_v = math.pi*core_r*core_r*l
    UN_to_cermet = 0.6

    pitch = 2*r*PD
    # calculate 'volumes' for fixed length
    v_cool = (r ** 2 * math.pi)
    # clad volume fraction
    v_clad = ((r + c)**2 - r**2)*math.pi
    # fuel volume fraction
    v_cermet = (math.sqrt(3)*pitch**2 / 2.0) - (r + c) ** 2 * math.pi 

    cell_vol = v_cool + v_clad + v_cermet
    # calculate vfracs from total cell volume
    vfrac_cermet = v_cermet / cell_vol
    
    fuel_vol = core_v * vfrac_cermet * UN_to_cermet
    
    fuel_frac = fuel_vol / core_v

    power_density = Q / (core_v / 1000) # W / l or MW/m^3

    mass = (fuel_vol * md.rho_UN) / 1000
    mass_dens = mass / (core_v / 1000)

    return mass, mass_dens, power_density, fuel_frac

def m235_mass(mass, enrich):
    
    mmu235 = 235.0439
    mmu238 = 238.0508
    mmN = 14.0067
    mmU = (mmu235 / enrich) + (mmu238 / (1-enrich))

    mfrac_U_in_UN = mmU / (mmU + mmN)
    
    return mass*mfrac_U_in_UN*enrich

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
        data[idx]['AR'] =  round(params[0], 5)
        data[idx]['PD'] = round(params[1], 5)
        data[idx]['cool_r'] = round(params[2], 5)
        data[idx]['core_r'] = round(params[3], 5)
        data[idx]['enrich'] = round(params[4], 5)
        data[idx]['power'] = round(params[5], 5)
        keff_vals = parse_keff(string)[2]
        data[idx]['keff'] = keff_vals[-1]
        data[idx]['dK'] = keff_vals[-1] - keff_vals[0]
        data[idx]['ave_E'] = parse_etal('1', string)[-1]
        mass, mdens, q_dens, fuel_frac = calc_fuel_mass(params[3], params[2], params[1], params[5])
        data[idx]['mass'] = mass
        data[idx]['fuel_frac'] = fuel_frac
        data[idx]['m_dens'] = mdens
        data[idx]['q_dens'] = round(q_dens, 5)
        data[idx]['mass_235'] = m235_mass(mass, params[4])

    np.savetxt("depl_results.csv", data, delimiter=',', 
           fmt='%10.5f', header=','.join(names))
  
def plot_results(data, ind, dep, colorplot=None, log=None, title=None):
    """Generate Plots
    """
    label_strings = {'AR' : 'Core Aspect Ratio[-]',
                     'PD' : 'Fuel Pitch to Coolant Channel Diameter',
                     'cool_r' : 'Coolant Channel [cm]',
                     'core_r' : 'Core Radius [cm]',
                     'enrich' : 'U-235 Enrichment [-]',
                     'power' : 'Core Thermal Power [kW]',
                     'keff' : r'k$_{eff}$ [-]',
                     'ave_E' : 'average neutron energy [MeV]',
                     'mass' : 'reactor fuel mass [kg]',
                     'fuel_frac' : 'volume fraction fuel [-]',
                     'q_dens' : 'volumetric power density [kW/l]',
                     'm_dens' : r'fuel mass density [kg/m$^3$]',
                     'dU' : 'Depleted $^{235}$U mass [kg]',
                     'dK' : r'$\delta$k',
                     'mass_235' : r'Mass of $^{235}$U [kg]',
                     'rel_depl' : r'Percent of $^{235}$U mass depleted',
                    }

    titles = {'keff' : {'core_r' : 'EOL keff Dependence on Core Radius',
                        'power'  : 'EOL keff Dependence on Thermal Power'
                       }
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
    if title:
        plt.title(title)
    else:
        plt.title("{0} vs. {1}".format(dep, ind))

    plt.axhline(y=-0.01, color='r')
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
    plt.savefig(savename, dpi=700, bbox_inches='tight')

    return


def surf_plot(data):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X = data['core_r']
    Y = data['PD']
    Z = data['keff']

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
#    save_store_data()
    data = load_from_csv()
    filter = ['mass > 100', 'mass < 500']
#    filter = ['enrich > 0.85']
#    filter = ['fuel_frac > 0.35']
    data = filter_data(filter, data)
#    surf_plot(data)
#   plot_results(data, 'mass_235', 'dK')
#   plot_results(data, 'mass', 'keff')
#   plot_results(data, 'mass_235', 'keff', None, 'log-log')
#   plot_results(data, 'mass_235', 'keff', 'm_dens', 'log-log')
#   plot_results(data, 'mass_235', 'keff', 'fuel_frac')
#   plot_results(data, 'mass_235', 'keff', 'power')
    plot_results(data, 'mass_235', 'keff', 'm_dens', None, r'EOL k$_{eff}$ Dependence on $^{235}$U Mass')
    plot_results(data, 'mass_235', 'keff', None, None, r'EOL k$_{eff}$ Dependence on $^{235}$U Mass')
#    plot_results(data, 'mass', 'dU', 'power')
#    filter = ['mass > 100']
#    dk_data = filter_data(filter, data)
    plot_results(data, 'mass', 'dK', None, None, r'10y Reactivity Swing')
