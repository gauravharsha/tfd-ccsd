import sys
sys.path.append('../src/')

import pytest
import numpy as np
import h5py

from iofuncs import *

#####################################################################
#                                                                   #
#   Functions to test the IO subsystem of the codebase              #
#                                                                   #
#####################################################################

def test_readinput():

    # 
    # Test to check the initialization of the IO class
    # -- This involves reading the Input file and updating the parameters
    # 

    iops = IOps(inp_file='TestInput')

    # According to the InputTest
    assert iops.fn == 'hub_6x1_u2_data.h5'
    assert iops.n_elec == 6
    assert iops.beta_f == 1
    assert iops.beta_pts == 101
    assert iops.ntol == 1e-4
    assert iops.deqtol == 1e-5
    assert iops.e_nuc == 0.0

def test_readintegrals():

    #
    # Test the `loadHDF' functionality of IOps class
    #

    iops = IOps(inp_file='TestInput')

    eigs, h1, eri, attrs = iops.loadHDF()

    f1 = h5py.File(iops.fn,'r')
    true_eigs = f1['eigs'][:]
    true_h1 = f1['h1'][:]
    true_eri = f1['eri'][:]

    assert eigs.all() == true_eigs.all()
    assert h1.all() == true_h1.all()
    assert eri.all() == true_eri.all()

    for (atr_name,val) in attrs:
        if atr_name == 'nsite':
            assert val == 6
        elif atr_name == 'nalfa':
            assert val == 3
        elif atr_name == 'nbeta':
            assert val == 3
        elif atr_name == 'intrn':
            assert val == 2.0
        else:
            pass

    # Check the number of spin orbitals
    assert iops.nso == 12

def test_output():

    #
    # Test the Output File Create and Edit functions
    #

    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()

    fout = h5py.File('test_output.h5','w')

    # Create only one data set called `energy'
    dsetname = ['energy']
    max_size = 1

    iops.createh5(fout, dsetname, max_size, attrs)

    # Update the energy variable to 10
    iops.updateh5([10],0)

    # Close the file for now
    fout.close()

    # check the updates
    fcheck = h5py.File('test_output.h5','r')
    en = fcheck['energy']

    assert en[0] == 10
