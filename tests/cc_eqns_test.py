import sys
sys.path.append('../src/')
sys.path.append('../fort_src/')

import pytest
import numpy as np
import h5py

from iofuncs import *
from inttran import *

from ThermalCCSD import *

#####################################################################
#                                                                   #
#   Functions to test the dimensionality and the accuracy of        #
#   the Fortran Libraries to compute the "Residuals" that           #
#   drive the imaginary time evolution for Thermal CCSD             #
#                                                                   #
#####################################################################

# 
# Set Up by loading the required one and two electron integrals
# We will test on both H2 in sto-3g and 6 site Hubbard with U/t = 2.0
# 

def test_cc_beta_residuals_zeroT():

    # 
    # Test to check the Coupled Cluster Beta CC Residuals in the Zero Temperature Limit
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

    # Get the residuals
    r0, r1, r2 = betacc(hdiag, eri, s1, s2, x, y)
    r1 = np.einsum('pq,p,q->pq',r1,x,y)
    r2 = np.einsum('pqrs,p,q,r,s->pqrs',r2,x,x,y,y)

    # First check the shapes of r0, r1, r2
    assert type(r0) == float
    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Compare with expected values (generated from Tom's MasterCode)
    assert np.max( np.abs( r1[nocc:,:nocc] ) ) <= 5e-8
    assert np.max( np.abs( r2[nocc:,nocc:,:nocc,:nocc] ) ) <= 5e-8

def test_cc_mu_residuals_zeroT():

    # 
    # Test to check the Coupled Cluster Mu Residuals in the Zero Temperature Limit
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

    # Get the residuals
    r0, r1, r2 = numbercc(s1, s2, x, y)
    r1 = np.einsum('pq,p,q->pq',r1,x,y)
    r2 = np.einsum('pqrs,p,q,r,s->pqrs',r2,x,x,y,y)

    print( np.max( np.abs(r1) ) )
    print( np.max( np.abs(r2) ) )

    # First check the shapes of r0, r1, r2
    assert type(r0) == float
    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Compare with expected values -- should all be zero
    assert r0 == 0
    assert np.max( np.abs( r1[nocc:,:nocc] ) ) <= 5e-8
    assert np.max( np.abs( r2[nocc:,nocc:,:nocc,:nocc] ) ) <= 5e-8
