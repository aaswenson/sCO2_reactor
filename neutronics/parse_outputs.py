import math
import operator
import material_data as md
import matplotlib.pyplot as plt
import numpy as np
import glob
import neutronic_sweeps as ns

results = ['keff', 'ave_E']

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
            data = line.split()[1].split(',')[:-1]

            data = [float(x) for x in data]
            break
    return data

def calc_fuel_mass(core_r, AR, PD, r):
    """
    """
    c = 0.0031
    l = core_r*AR
    core_v = math.pi*core_r*core_r*l

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
    
    fuel_vol = core_v * vfrac_cermet

    return (fuel_vol * md.rho_UN) / 1000

def save_store_data(data_dir='data/*i.o'):
    """
    """
    files = glob.glob(data_dir)
    names = ns.dimensions + results + ['mass']
    types = ['f8']*len(names)
    N = len(files)
    data = np.zeros(N, dtype={'names' : names, 'formats' : types})

    for idx, file in enumerate(files):
        fp = open(file, 'r')
        string = fp.readlines()
        params = parse_header_string(string)
        params.append(parse_keff(string)[2][-1])
        params.append(parse_etal('1', string)[-1])
        params.append(calc_fuel_mass(params[1], params[0],params[3],params[2]))
        # save data and close file
        data[idx] = tuple(params)
        fp.close()

    np.savetxt("depl_results.csv", data, delimiter=',', fmt='%10.5f',
               header=','.join(names))

def plot_results(data, ind, dep, colorplot=None):
    """Generate Plots
    """
    label_strings = {'AR' : 'Core Aspect Ratio[-]',
                     'PD' : 'Fuel Pitch to Coolant Channel Diameter',
                     'cool_r' : 'Coolant Channel [cm]',
                     'core_r' : 'Core Radius [cm]',
                     'enrich' : 'U-235 Enrichment [-]',
                     'power' : 'Core Thermal Power [kW]',
                     'keff' : 'k-eff [-]',
                     'ave_E' : 'average neutron energy [MeV]',
                     'mass' : 'fuel mass [kg]'
                    }
    # plot
    fig = plt.figure()
    colorname = ''
    if colorplot:
        plt.scatter(data[ind], data[dep], s=6, c=data[colorplot],
                cmap=plt.cm.get_cmap('plasma', len(set(data[colorplot]))))
        plt.colorbar(label=label_strings[colorplot])
        colorname = colorplot
    else:
        plt.scatter(data[ind], data[dep], s=6)
    # titles and labels
    plt.title("{0} vs. {1}".format(dep, ind))
    plt.xlabel(label_strings[ind])
    plt.ylabel(label_strings[dep])
    
    plt.savefig("{0}_vs_{1}_{2}.png".format(dep, ind, colorname), dpi=500)
    return plt

def filter_data(filters, data):

    """Apply useful filters on the data
    """
    opers = {'less' : operator.lt,
             'equal' : operator.eq,
             'great' : operator.gt}
    
    for op in filters:
        data = data[opers[op[1]](data[op[0]], op[2])]   
    
    return data
    
def get_min_mass(data):
    """ After the parametric sweep is complete, find the minimum calculated
    fuel mass.
    """
    # search the results for minimum-mass configuration
    min_idx = list(data['mass']).index(min(data['mass']))

    # get data for min mass config
    min_mass = data[min_idx]['mass']

    # return the min_idx to get other results at minimum config
    return min_idx


def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    fp = open(datafile, 'r')
    names = fp.readlines()[0].strip('#').split(',')
    data = np.genfromtxt(datafile, delimiter=',',
            names= names)
    
    return data

if __name__ == '__main__':
#    save_store_data()
    data = load_from_csv()
#    data = filter_data([('keff', 'great', 1)], data)
    plt = plot_results(data, 'AR', 'keff', 'mass')
    plt.show()
    idx = get_min_mass(data)
    print(data[idx]['mass'])
    print(data[idx])
