"""This module contains functions to sweep parameter-space, creating MCNP
depletion inputs to calculate keff at EOL.
* AR
* core_z
* cool_r
* PD
* power
* enrich
"""
from pyDOE import lhs
import os
import numpy as np
import tarfile

from mcnp_inputs import HomogeneousInput

# set seed for reproducibility
np.random.seed(1324291)

parameters = {'core_r'  : (10, 50),         
              'AR'      : (0.7, 1.3),
              'cool_r'  : (0.5, 1),
              'PD'      : (1.1, 1.6),        
              'power'   : (80, 200),        
              'enrich'  : (0.2, 0.9)
             }

dim = len(parameters.keys())
samples = 2

def gen_hypercube(samples, N):
    """Generate N-dimensional latin hypercube to sample dimensional reactor
    space.

    Arguments:
    ----------
        samples (int): number of test cases to try
        N (int): number of dimensions to test
    
    Returns:
    --------
        cube (ndarray): normalized, N-D latin hypercube
    """

    np.random.seed(4654562)
    hypercube = lhs(N, samples=samples)

    return hypercube

def fill_data_array(samples, parameters, cube):
    """Fill an ndarray with the sampling set generated by lhs.
    """
    # initialize array
    test_cases = np.zeros(samples, dtype={'names' : list(parameters.keys()) +\
                                           ['keff', 'mass'],
                                    'formats' : ['f8']*8 })
    # for all samples
    for sample_idx, sample in enumerate(cube):
        # get values for every dimension
        for dim_idx, dim in enumerate(sorted(parameters.keys())):
            l_limit = parameters[dim][0]
            u_limit = parameters[dim][1]
            # uniform distribution
            a = u_limit - l_limit
            b = l_limit
            # save to ndarray
            test_cases[sample_idx][dim] = b + cube[sample_idx][dim_idx] * a
    
    return test_cases

def write_inputs(sampling_data):
    """Write MCNP depletion inputs for sampled data.
    """
    datanames = sampling_data.dtype.names
    tarputs = tarfile.open('smpl_mcnp_depl_inps.tar', 'w')
    for num, sample in enumerate(sampling_data):
        input = HomogeneousInput(sample['core_r'],
                                 sample['core_r']*sample['AR'],
                                 sample['power'])
        homog_comp = input.homog_core(sample['enrich'],
                                      sample['cool_r'],
                                      sample['PD'])
        input.write_mat_string(homog_comp)
        
        # identifying header string for post-processing
        header_str = ''
        for param in sorted(datanames):
            header_str += str(round(sample[param], 5)) + ','
        # write the input and tar it
        filename = input.write_input(num, header_str)
        tarputs.add(filename)

    tarputs.close()


if __name__=='__main__':
    cube = gen_hypercube(samples, dim)
    data = fill_data_array(samples, parameters, cube)
    write_inputs(data)
    os.system('rm *.i')
