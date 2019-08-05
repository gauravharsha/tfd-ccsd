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

    # Eval Energy -- directly from fortran routines
    energy_f = evalenergy(hdiag, eri, s1, s2, s1*0, s2*0, x, y)
    number_f = evalnumber(s1, s2, s1*0, s2*0, x, y)

    # Eval Energy -- from python wrappers
    cc_amps = np.concatenate(([0],np.reshape(s1,nso**2),CompressT2(s2)))
    ci_amps = cc_amps * 0.0

    energy_p = eval_energy(hdiag, eri, cc_amps, ci_amps, x, y)
    number_p = eval_number(cc_amps, ci_amps, x, y)

    # Eval Twobody expectation value for Energy and number
    energy_twobody = evaltwobodyh(np.diag(hdiag), eri, s1, s2, s1*0, s2*0, x, y)
    number_twobody = evaltwobodyh(np.eye(nso), eri*0, s1, s2, s1*0, s2*0, x, y)

    # Expected Values
    en_exp = -5.408955909508
    num_exp = 6.0

    # Check
    assert np.abs( energy_f - en_exp ) < 5e-8
    assert np.abs( number_f - num_exp ) < 5e-8
    assert np.abs( energy_p - en_exp ) < 5e-8
    assert np.abs( number_p - num_exp ) < 5e-8
    assert np.abs( energy_twobody - en_exp ) < 5e-8
    assert np.abs( number_twobody - num_exp ) < 5e-8

