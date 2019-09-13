import sys
sys.path.append('../src/')
sys.path.append('../fort_src/')

import pytest
import numpy as np
import h5py

from iofuncs import *
from inttran import *

from ThermalCISD import *

#####################################################################
#                                                                   #
#   Functions to test the dimensionality and the accuracy of        #
#   the Fortran Libraries to compute the "Residuals" that           #
#   drive the imaginary time evolution for Thermal CISD or the      #
#   Linear Response equations                                       #
#                                                                   #
#####################################################################

#
# NOTE: Since we do not have an existing zero temperature CISD code, we will
# instead use a trick to Convert CC amplitudes that are returned by
# Tom's MasterCode - then construct CI amplitudes and check residuals
# -- Should be Zero
# 


# 
# Set Up by loading the required one and two electron integrals
# We will test on 6 site Hubbard with U/t = 2.0
# 

def test_ci_beta_residuals_zeroT():

    # 
    # Test to check the Coupled Cluster Residuals in the Zero Temperature Limit
    # 

    # Load the Integrals first
    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()
    nso = iops.nso

    # Construct random t1 and t2 = use Tom's MasterCode
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

    # Convert CC to CI amps
    t2 += np.einsum('ai,bj->abij',t1,t1)/2

    # Other parameters needed for residuals
    y = np.zeros(nso)
    for i in range(nocc):
        y[i] = 1.0
    x = np.sqrt(1 - y**2)

    # Get the residuals
    r0, r1, r2 = betaci(eigs, h1, eri, t1, t2, x, y)
    print(np.max(np.abs(r2)))

    # First check the shapes of r0, r1, r2
    assert type(r0) == float
    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Compare with expected values (generated from Tom's MasterCode)
    # NOTE: Because T1 = 0 for Hubbard, the R1 = 0 for CI that we have constructed here
    #       but R2 will not be zero as one needs to have up to 4th order wavefunction to
    #       get CI equivalent of CC.

    assert np.max( np.abs( r1 ) ) <= 5e-8
    # assert np.max( np.abs( r2 ) ) <= 5e-8

def test_ci_beta_residuals_zeroT():

    # 
    # Test to check the Coupled Cluster Residuals in the Zero Temperature Limit
    # 

    # Load the Integrals first
    iops = IOps(inp_file='TestInput')
    eigs, h1, eri, attrs = iops.loadHDF()
    nso = iops.nso

    # Construct random t1 and t2 = use Tom's MasterCode
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

    # Convert CC to CI amps
    t2 += np.einsum('ai,bj->abij',t1,t1)/2

    # Other parameters needed for residuals
    y = np.zeros(nso)
    for i in range(nocc):
        y[i] = 1.0
    x = np.sqrt(1 - y**2)

    # Get the residuals
    r0, r1, r2 = numberci(t1, t2, x, y)
    print(np.max(np.abs(r2)))

    # First check the shapes of r0, r1, r2
    assert type(r0) == float
    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Compare with expected values -- for number operator as hamiltonian
    #   these should all be zero

    assert r0 == 0
    assert np.max( np.abs( r1 ) ) <= 5e-8
    assert np.max( np.abs( r2 ) ) <= 5e-8
