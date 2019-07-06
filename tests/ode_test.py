import sys
sys.path.append('../src/')
sys.path.append('../fort_src/')

import pytest, h5py, numpy as np
from scipy.misc import comb

from iofuncs import *
from inttran import *
from odefuncs import *

from ThermalCCSD import *
from ThermalCISD import *


#########################################################################
#                                                                       #
#   Testing the Beta and Alpha ODE drivers and related functiosn        #
#                                                                       #
#########################################################################

def test_compress_decompress_t2():
    
    #
    # Test the function which compresses and de-compresses T2-like 
    # amplitudes based on the symmetries
    # 

    # size
    nso = 10

    # t2
    x2 = np.random.rand(nso,nso,nso,nso)

    # anti-symmetrize
    t2 = x2 - np.einsum('qprs->pqrs',x2) \
            - np.einsum('pqsr->pqrs',x2) \
            + np.einsum('qpsr->pqrs',x2)
    t2 /= 4.0

    # Compress
    t2_compr = CompressT2(t2)

    # Check the length of the compressed t2
    assert np.shape(t2_compr) == (int(comb(nso,2)**2), )

    # DeCompress
    t2_decompr = DecompressT2(t2_compr,nso)

    # Check the shape of the decompressed t2
    assert np.shape(t2_decompr) == (nso,nso,nso,nso)

    # de compressed version should be same as original t2
    assert np.max( np.abs( t2_decompr - t2 ) ) <= 5e-8

def test_evolution_class_attributes():

    #
    # Test the attributes of the `Evolution` class
    #

    # First read the file and the input and do the formalities
    iops = IOps(inp_file='TestInput')
    h1, eri, attrs = iops.loadHDF()
    hdiag, evecs = IntTran2(h1)
    eri_trans = IntTran4(eri,evecs)

    fug = iops.n_elec/(iops.nso-iops.n_elec)
    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    assert evol.beta_step == 0.1
    assert evol.alpha_step == 0.1
    assert evol.nso == iops.nso
    assert evol.n_elec == iops.n_elec
    assert evol.fug == fug
    assert evol.h1.all() == hdiag.all()
    assert evol.eri.all() == eri_trans.all()
    assert evol.deqtol == iops.deqtol
    assert evol.ntol == iops.ntol

    # amplitudes -- before I set anything
    assert evol.cc_amps.all() == evol.ci_amps.all()
    assert evol.cc_amps.all() == np.zeros(1 + int(iops.nso**2 + comb(iops.nso,2)**2)).all()

    # amplitudes -- after changing
    cc_amps = np.random.rand(1+int(iops.nso**2 + comb(iops.nso,2)**2))
    ci_amps = np.random.rand(1+int(iops.nso**2 + comb(iops.nso,2)**2))
    evol.setAmps(cc_amps=cc_amps,ci_amps=ci_amps)

    assert evol.cc_amps.all() == cc_amps.all()
    assert evol.ci_amps.all() == ci_amps.all()



def test_evolution_class_beta_integration():

    #
    # Test the attributes of the `Evolution` class
    #

    # First read the file and the input and do the formalities
    iops = IOps(inp_file='TestInput')
    h1, eri, attrs = iops.loadHDF()
    hdiag, evecs = IntTran2(h1)
    eri_trans = IntTran4(eri,evecs)

    fug = iops.n_elec/(iops.nso-iops.n_elec)
    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    # amplitudes -- before I set anything
    assert evol.cc_amps.all() == evol.ci_amps.all()
    assert evol.cc_amps.all() == np.zeros(1 + int(iops.nso**2 + comb(iops.nso,2)**2)).all()

    # Do Beta evolution
    evol.DoBetaIntegration()
    cc_amps0 = evol.cc_amps


def test_evolution_class_eval_number():

    #
    # Test the attributes of the `Evolution` class
    #

    # First read the file and the input and do the formalities
    iops = IOps(inp_file='TestInput')
    h1, eri, attrs = iops.loadHDF()
    hdiag, evecs = IntTran2(h1)
    eri_trans = IntTran4(eri,evecs)

    fug = iops.n_elec/(iops.nso-iops.n_elec)
    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    # amplitudes -- before I set anything
    assert evol.cc_amps.all() == evol.ci_amps.all()
    assert evol.cc_amps.all() == np.zeros(1 + int(iops.nso**2 + comb(iops.nso,2)**2)).all()
    cc_amps0 = evol.cc_amps

    x = 1/np.sqrt( 1 + np.exp(-evol.beta_in*evol.h1 + evol.alpha_in)*evol.fug )
    y = np.sqrt( 1 - x**2 )

    n0 = eval_number(evol.cc_amps, evol.ci_amps, x, y)

    assert np.max(np.abs(n0 - evol.n_elec)) < 5e-8


def test_evolution_class_alpha_integration():

    #
    # Test the attributes of the `Evolution` class
    #

    # First read the file and the input and do the formalities
    iops = IOps(inp_file='TestInput')
    h1, eri, attrs = iops.loadHDF()
    hdiag, evecs = IntTran2(h1)
    eri_trans = IntTran4(eri,evecs)

    fug = iops.n_elec/(iops.nso-iops.n_elec)
    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    # amplitudes -- before I set anything
    assert evol.cc_amps.all() == evol.ci_amps.all()
    assert evol.cc_amps.all() == np.zeros(1 + int(iops.nso**2 + comb(iops.nso,2)**2)).all()
    cc_amps0 = evol.cc_amps

    x = 1/np.sqrt( 1 + np.exp(-evol.beta_in*evol.h1 + evol.alpha_in)*evol.fug )
    y = np.sqrt( 1 - x**2 )

    # Do Beta evolution - should not do anything!!
    evol.BisectionAndAlphaIntegrate()

    # amplitudes -- after evolving -- check if it evolve at all
    assert np.max(np.abs( evol.cc_amps - cc_amps0 )) < 5e-8

