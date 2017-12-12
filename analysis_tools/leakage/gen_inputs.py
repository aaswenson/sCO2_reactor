from base_string import base_string
import numpy as np

def write_input(radius, height, kcode_params, bc):
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
                           comm = "$")
    # write the file
    ifile = open("leakage_{0}_{1}.i".format(radius, height),'w')
    ifile.write(file_string)
    ifile.close()

def main():
    for radius in range(10,110,10):
        write_input(radius, 50, (15000, 15, 65), 'bare')
        write_input(radius, 50, (15000, 15, 65), 'refl')


if __name__ == "__main__":
    main()
