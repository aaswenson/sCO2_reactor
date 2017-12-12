from base_string import base_string
from pyne.material import Material
from pyne import mcnp
import numpy as np
import math

def write_input(radius, height, mat_card, kcode_params, bc):
    bc_options = {'refl' : '+',
                  'bare' : ''}
    templ = base_string
    file_string = templ.substitute(core_radius = radius,
                           core_height = height,
                           bc = bc_options[bc],
                           half_core_height = height / 2.0,
                           n_per_cycle = kcode_params[0],
                           non_active_cycles = kcode_params[1],
                           total_cycles = kcode_params[2],
                           fuel_mat = mat_card,
                           comm = "$")
    # write the file
    ifile = open("./inputs/leakage_{0}_{1}_{2}.i".format(radius, height, bc),'w')
    ifile.write(file_string)
    ifile.close()

def homog_fuel_comp(r, PD, enrich):
    """Create mcnp data card for fast reactor fuel comp base on vol frac of fuel
    """
    # important constants relevant to W-UO2 CERMET fuel
    vol_frac_U = 0.7
    
    # get fuel/coolant volume fractions
    pitch = r*PD*2
    total = math.sqrt(3)*pitch*pitch / 2
    fuel_frac = (total - r**2*math.pi) / total
    cool_frac = 1 - fuel_frac

    # smear densities
    rho_UN = 11.9
    rho_W = 19.3
    rho_CO2 = 87.13
    cermet_rho = (rho_UN*vol_frac_U) + (1-vol_frac_U)*rho_W
    total_rho = cermet_rho*fuel_frac + rho_CO2 * cool_frac

    # get mass fractions
    mfrac_fuel = ((rho_UN*vol_frac_U) / cermet_rho) *\
                  (cermet_rho*fuel_frac) / total_rho
    mfrac_W = (rho_W*(1-vol_frac_U) / cermet_rho) *\
               (cermet_rho*fuel_frac) / total_rho
    mfrac_CO2 = rho_CO2 * cool_frac / total_rho

    HEU = Material()
    fuel = Material()
    coolant = Material({'C':0.333, 'O':0.667}, mfrac_CO2)
    cermet_matrix = Material({'W':1}, mfrac_W)

    # create HEU fuel material
    HEU = Material({'U238':1-enrich, 'U235':enrich}, mfrac_fuel)
    afrac_fuel = HEU.to_atom_frac()
    
    # add Nitrogen to fuel and adjust atom fractions
    afrac_fuel[70140000] = 1.0
    for isotope in afrac_fuel:
        afrac_fuel[isotope] /= sum(afrac_fuel.values())    
    fuel.from_atom_frac(afrac_fuel)

    # combine three fuel constituents into mass-weighted mixture
    homogeneous_fuel = fuel + coolant + cermet_matrix
    homogeneous_fuel.metadata['mat_number'] = 1

    return homogeneous_fuel.mcnp()

def main():
    # get homogeneous fuel for r=1cm, PD=1.2 90% U235
    fuel_string = homog_fuel_comp(1, 1.2, 0.9)
    for radius in range(10,110,10):
        write_input(radius, 50, fuel_string, (15000, 15, 65), 'bare')
        write_input(radius, 50, fuel_string, (15000, 15, 65), 'refl')


if __name__ == "__main__":
    main()
    homog_fuel_comp(1,1.1,0.9)
