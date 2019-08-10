import sys
sys.path.append('../src/')

import pytest
import numpy as np
import h5py

from iofuncs import *
from inttran import *

#####################################################################
#                                                                   #
#   Functions to test the IO subsystem of the codebase              #
#                                                                   #
#####################################################################

def test_inttran2():

    #
    # IntTran2 should return the diagonalized Fock matrix
    #   with corresponding eigenvectors
    #

    # First read the file and the input
    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()

    # Next do the IntTran2
    hdiag, evecs = IntTran2(h1)

    # Check if it is diag by reconstruction
    diag_mat = np.diag( hdiag )
    h1_recon = evecs @ diag_mat @ np.transpose(evecs)

    assert np.max( np.abs( h1_recon - h1 ) ) < 1e-8

def test_inttran4():

    #
    # IntTran4 takes in as input the 2-electron integral and
    # transforms it into the basis specified by the evecs
    #

    # First read the file and the input
    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()

    # Next do the IntTran2
    hdiag, evecs = IntTran2(h1)

    # Do the IntTran4
    eri_trans = IntTran4(eri,evecs)

    # Do the transformation by hand
    eri_trans2 = np.einsum('pqrs,sd->pqrd',eri,evecs)
    eri_trans2 = np.einsum('pqrd,rc->pqcd',eri_trans2,evecs)
    eri_trans2 = np.einsum('pqcd,qb->pbcd',eri_trans2,evecs)
    eri_trans2 = np.einsum('pbcd,pa->abcd',eri_trans2,evecs)

    assert np.amax( np.abs( eri_trans2 - eri_trans ) ) < 1e-8
