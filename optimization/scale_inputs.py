from pyne.material import Material, MaterialLibrary
from pyne import mcnp
from string import Template
import numpy as np
import math
import sys
import physical_constants as pc

class PinCellSCALE:
    """Class to write MCNP input files for pin cell modeling in an infinite
    lattice (in the x-y directions).
    """
    # standard material numbers for fuel, clad, coolant
    mat_numbers = {'fuel' : 1, 'clad' : 2, 'cool' : 3}

    # base template string modified by the methods below
    base_string = Template("""=t-newt
Scale pin cell study radius:${radius} cm
v5-44

read comp
${mat}
end comp

read celldata
latticecell squarepitch hpitch=${hpitch} 3 fuelr=${fuelr} 1 cladr=${c_rad} 2 end
end celldata

read model

read parameters
  run=yes
  epseigen=0.001
  converge=mix
  collapse=yes
  drawit=yes
  echo=yes
  prtflux=yes
  prthmmix=yes
  timed=yes
end parameters

read hmog
500 fuel 1 end
501 nonfuel 2 end
502 allsmear 1 2 3 end
end hmog

read collapse
30r1 14r2
end collapse

read materials
  mix=1 pn=1 com='UN fuel' end
  mix=2 pn=1 com='clad' end
  mix=3 pn=3 com='sCO2 coolant' end
end materials

read geometry
global unit 1
cylinder 11 ${radius}  sides=20
cylinder 12 ${clad_radius} sides=20
rhexprism 13 ${pitch}

media 3 1 11
media 2 1 12 -11
media 1 1 13 -12

boundary 13 40 40
end geometry

end model
end
""")

    m_to_cm = 100
    
    def __init__(self, TH_reactor):
        """Initialize parameters.
        """
        self.r = TH_reactor.r_channel * self.m_to_cm
        self.c = TH_reactor.c * self.m_to_cm
        self.PD = TH_reactor.pd_ratio
        self.pitch = (self.r + self.c) * self.PD * 2.0
        self.load_pyne_matlib()
        self.fps = TH_reactor.fps
    
    def load_pyne_matlib(self):
        """Load pyne material library.
        """
        # Initialize material libraries.
        path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
        self.raw_matlib = MaterialLibrary()
        # Write entire PyNE material library.
        self.raw_matlib.from_hdf5(path_to_compendium,
                             datapath="/material_library/materials",
                             nucpath="/material_library/nucid")
    
    def equivalent_fuel_ratio(self):
        """
        """
        A_fuel = math.sqrt(3)*self.pitch**2 / 2.0 -\
                      (self.r + self.c) ** 2 * math.pi
        A_flow = self.r**2 * math.pi
        A_clad = math.pi*((self.r + self.c)**2 - self.r**2)
        A_total = A_fuel + A_flow + A_clad
        
        # get equivalent fuel ratio
        self.fuelr = math.sqrt(A_fuel / math.pi)
        self.hpitch = math.sqrt(A_total) / 2.0
        self.cladr = math.sqrt((A_clad / math.pi) + self.fuelr**2)
    
    def write_fuel_string(self, enrich, fuel_type):
        """Get fuel material, enrich it and write fuel string.
        """
        kg_m3_to_g_cc = 0.001
        # get the fuel-bonded compound
        const_mat = self.raw_matlib[fuel_type[0]]
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

    def write_mat_string(self, clad_mat, coolant_mat):
        """Write the material data.
        """
        kg_m3_to_g_cc = 0.001
        
        self.cool = self.raw_matlib[coolant_mat].collapse_elements([])
        self.clad = self.raw_matlib[clad_mat].collapse_elements([])
        
        # delete O-18 because MCNP gets cranky
        del self.fuel['8018']; del self.cool['8018']; del self.clad['8018']

        # add material numbers
        self.cool.metadata['mat_number'] = self.mat_numbers['cool']
        self.clad.metadata['mat_number'] = self.mat_numbers['clad']
        
        self.cool.density = self.fps.rho*kg_m3_to_g_cc
        self.clad.density = pc.rho_W*kg_m3_to_g_cc
        self.fuel.density = pc.rho_fuel*kg_m3_to_g_cc

        self.mat_string = self.fuel.scale(300)
        self.mat_string += self.clad.scale(300) +\
                           self.cool.scale(300)
    
    def write_input(self):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        self.equivalent_fuel_ratio()
        templ = self.base_string
        file_string = templ.substitute(
                               radius = self.r,
                               clad_radius = self.r + self.c,
                               pitch = self.pitch,
                               mat = self.mat_string,
                               hpitch = self.hpitch,
                               fuelr = self.fuelr,
                               c_rad = self.cladr,
                               comm = "$")
        # write the file
        ifile_name = "./inputs/scale/leakage_{0}_{1}.inp".format(round(self.r, 5),
                                              round(self.PD, 5))
        ifile = open(ifile_name, 'w')
        ifile.write(file_string)
        ifile.close()
        
        return ifile_name


def main():
    """An example of using the PinCellSCALE class to write an input file.
    """
    # Initialize material libraries.
    path_to_compendium = "/home/alex/.local/lib/python2.7/\
site-packages/pyne/nuc_data.h5"
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(path_to_compendium,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")

    test = PinCellSCALE(0.5, 1.5, 0.031)
    
    test.write_fuel_string(0.9, ('Nitrogen', 1), raw_matlib)
    test.write_mat_string('Inconel-718', 'Carbon Dioxide', raw_matlib)
    test.write_input()

if __name__ == '__main__':
    main()
