from pyne.material import Material, MaterialLibrary
from pyne import mcnp
from string import Template
import numpy as np
import math
import sys
sys.path.append('../')
import physical_constants as pc


class PinCellMCNP:
    """Class to write MCNP input files for pin cell modeling in an infinite
    lattice (in the x-y directions).
    """
    # standard material numbers for fuel, clad, coolant
    mat_numbers = {'fuel' : 1, 'clad' : 2, 'cool' : 3}

    # base template string modified by the methods below
    base_string = Template("""\
MCNP6 pin cell study radius:${radius} cm
c Cell Card
1 3 ${cool_rho} -1 4 -5 imp:n=1   ${comm} coolant channel
2 2 ${clad_rho} 1 -2 4 -5 imp:n=1 ${comm} cladding
3 1 ${fuel_rho} -3 2 4 -5 imp:n=1 ${comm} fuel
99 0 5:-4:(3 -5 4) imp:n=0  ${comm} outside world
 
c Surface Card
1 CZ ${radius}
2 CZ ${clad_radius}
3+ rhp 0 0 -10000 0 0 100000 ${pitch} 0
4+ PZ -0.5
5+ PZ 0.5

c Data Card
${mat}\
${VOL}
kcode ${n_per_cycle} 1 ${non_active_cycles} ${total_cycles}
ksrc  0 ${ksrc_radius} 0
      0 -${ksrc_radius} 0
      ${ksrc_radius} 0 0
      -${ksrc_radius} 0 0
mode n
print
""")

    def __init__(self, radius, PD, clad_t):
        """Initialize parameters.
        """
        self.r = radius * 100
        self.c = clad_t * 100
        self.PD = PD
        self.pitch = (self.r + self.c) * self.PD * 2.0

    def vol_card(self):
        """ Generate cell volumes.
        """
        cool_vol = self.r**2 * math.pi
        fuel_vol = math.sqrt(3)*self.pitch**2 / 2.0 -\
            (self.r + self.c) ** 2 * math.pi
        clad_vol = ((self.r + self.c)**2 - self.r**2) * math.pi
    
        self.vol_str = "VOL {0} {1} {2}".format(cool_vol, clad_vol, fuel_vol)

    def write_fuel_string(self, enrich, fuel_type, matlib):
        """Get fuel material, enrich it and write fuel string.
        """
        # get the fuel-bonded compound
        const_mat = matlib[fuel_type[0]]
        const_mat.expand_elements().to_atom_frac()
        
        # get atom fractions of U and fuel-bonded compound
        Uafrac = 1 / (1 + fuel_type[1])
        constafrac = 1 - Uafrac

        afrac = Material({'U238':1-enrich, 'U235':enrich},-1.0).to_atom_frac()
        
        # account for uranium atom fraction
        afrac.update({n: Uafrac*afrac[n] for n in afrac.keys()})
        
        # add fuel-bonded compound, adjusted for U:compound
        for isotope in const_mat:
            afrac.update({isotope:const_mat[isotope]*constafrac})

        self.fuel = Material()
        self.fuel.from_atom_frac(afrac)
        
        self.fuel.metadata['mat_number'] = self.mat_numbers['fuel']

    def write_mat_string(self, clad_mat, coolant_mat, matlib):
        """Write the material data.
        """
        kg_m3_to_g_cc = 0.001
        self.mat_string = self.fuel.mcnp(frac_type='atom')
        
        self.cool = matlib[coolant_mat]
        self.clad = matlib[clad_mat]
        
        # delete O-18 because MCNP gets cranky
        del self.fuel['8018']; del self.cool['8018']; del self.clad['8018']

        # add material numbers
        self.cool.metadata['mat_number'] = self.mat_numbers['cool']
        self.clad.metadata['mat_number'] = self.mat_numbers['clad']

        self.fuel.density = pc.rho_fuel*kg_m3_to_g_cc
        self.cool.density = pc.rho_cool*kg_m3_to_g_cc
        self.clad.density = pc.rho_W*kg_m3_to_g_cc
        
        self.mat_string += self.clad.mcnp(frac_type='atom') +\
                           self.cool.mcnp(frac_type='atom')
    
    def write_input(self, kcode_params=[10000, 15, 60]):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        self.vol_card()
        templ = self.base_string
        file_string = templ.substitute(
                               cool_rho = abs(self.cool.number_density() / 1e24),
                               clad_rho = abs(self.clad.number_density() / 1e24),
                               fuel_rho = abs(self.fuel.number_density() / 1e24),
                               radius = self.r,
                               clad_radius = self.r + self.c,
                               pitch = self.pitch,
                               ksrc_radius = self.r + self.c + 0.01,
                               mat = self.mat_string,
                               VOL = self.vol_str,
                               n_per_cycle = kcode_params[0],
                               non_active_cycles = kcode_params[1],
                               total_cycles = kcode_params[2],
                               comm = "$")
        # write the file
        ifile = open("./inputs/mcnp/leakage_{0}_{1}.i".
                format(round(self.r, 5), round(self.PD, 5)),'w')
        ifile.write(file_string)
        ifile.close()




def main():
    """An example of using the PinCellMCNP class to write an input file.
    """
    # Initialize material libraries.
    path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(path_to_compendium,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")

    test = PinCellMCNP('hex', 0.5, 1.5, 0.031, 25)
    
    test.write_fuel_string(0.9999, ('Nitrogen', 1), raw_matlib)
    test.write_mat_string('Tungsten', 'Carbon Dioxide', raw_matlib)
    test.write_input([1,1,1])

if __name__ == '__main__':
    main()
