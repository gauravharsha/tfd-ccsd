import sys
sys.path.append('../src/')
sys.path.append('../fort_src/')

import pytest, h5py, numpy as np
from scipy.misc import comb

from iofuncs import *
from inttran import *
from odefuncs import *

from ExpVals import *

#########################################################################
#                                                                       #
#   Testing the various Fortran to Python functions                     #
#                                                                       #
#########################################################################

def test_eval_energ_and_num():

    #
    # Test the energy eval function -- we can check the HF and CC limit
    #

    # Load the Integrals first
    iops = IOps(inp_file='TestInput')
    h1, eri_in, attrs = iops.loadHDF()
    nso = iops.nso

    # Next do the IntTran2
    hdiag, evecs = IntTran2(h1)
    eri = IntTran4(eri_in,evecs)

    # Construct t1 and t2 = use Tom's MasterCode
    nocc = 6
    t1 = np.loadtxt('T1AMP')
    t1_pre = np.zeros((nso,nso))
    t1_pre[:nocc,nocc:] = t1[:,:]
    t1 = np.einsum('ia->ai',t1_pre)

    t2 = np.loadtxt('T2AMP')
    t2_pre = eri_in*0
    m = 0
    for i in range(nocc):
        for j in range(nocc):
            for a in range(nocc):
                t2_pre[i,j,nocc+a,nocc:] = t2[m,:]
                m += 1
    t2 = np.einsum('ijab->abij',t2_pre)

    # Transform the basis
    s1 = np.transpose(evecs) @ t1 @ evecs
    s2 = IntTran4(t2,evecs)

    # Other parameters needed for residuals
    y = np.zeros(nso)
    for i in range(nocc):
        y[i] = 1.0
    x = np.sqrt(1 - y**2)

    # Eval Energy
    energy = evalenergy(hdiag, eri, s1, s2, s1*0, s2*0, x, y)
    number = evalnum(s1, s2, s1*0, s2*0, x, y)

    # Expected Values
    en_exp = -5.408955909508
    num_exp = 6.0

    # Check
    assert np.abs( energy - en_exp ) < 5e-8
    assert np.abs( number - num_exp ) < 5e-8


# def test_res_beta_cc()
# 
# def test_res_beta_ci()
# 
# def test_res_alpha_cc()
# 
# def test_res_alpha_ci()
