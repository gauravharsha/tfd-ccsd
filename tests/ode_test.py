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

    assert evol.beta_step == 1e-2
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



def test_cc_beta_evolve():
    
    #
    # Test the integration driver functions
    #

    # Define the evolution class instance
    evol = Evolution(inp_file='TestInput')

    beta_step = evol.beta_step
    alpha_step = evol.alpha_step
    nso = evol.nso
    n_elec =  evol.n_elec
    nso = evol.nso
    fug =  evol.fug
    h1 = evol.h1
    eri = evol.eri
    deqtol = evol.deqtol
    ntol = evol.ntol
    evecs = evol.evecs

    # Construct t1 and t2 = use Tom's MasterCode
    nocc = 6
    t1 = np.loadtxt('T1AMP')
    t1_pre = np.zeros((nso,nso))
    t1_pre[:nocc,nocc:] = t1[:,:]
    t1 = np.einsum('ia->ai',t1_pre)

    t2 = np.loadtxt('T2AMP')
    t2_pre = eri*0
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

    # Compress and concatenate
    t0 = 0.0
    t1 = np.reshape( t1, int(nso**2) )
    t2 = CompressT2( t2 )
    cc_amps = np.concatenate(([t0],t1,t2))

    beta_in = 10.0
    alpha_in = beta_in*h1[nocc]

    # Check outputs from CC beta evolution
    yout = cc_beta_evolve(beta_in, cc_amps, alpha_in, fug, h1, eri)

    r0 = yout[0]
    r1 = np.reshape(yout[1:1+nso**2],(nso,nso))
    r2 = DecompressT2(yout[1+nso**2:],nso)

    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Check outputs from CC alpha evolution
    yout = cc_alpha_evolve(alpha_in, cc_amps, beta_in, fug, h1)

    r0 = yout[0]
    r1 = np.reshape(yout[1:1+nso**2],(nso,nso))
    r2 = DecompressT2(yout[1+nso**2:],nso)

    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Check outputs from CI beta evolution
    yout = ci_beta_evolve(beta_in, cc_amps, alpha_in, fug, h1, eri)

    r0 = yout[0]
    r1 = np.reshape(yout[1:1+nso**2],(nso,nso))
    r2 = DecompressT2(yout[1+nso**2:],nso)

    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Check outputs from CI alpha evolution
    yout = ci_alpha_evolve(alpha_in, cc_amps, beta_in, fug, h1)

    r0 = yout[0]
    r1 = np.reshape(yout[1:1+nso**2],(nso,nso))
    r2 = DecompressT2(yout[1+nso**2:],nso)

    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)




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

    assert np.max( np.abs(cc_amps0) ) > 0


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

