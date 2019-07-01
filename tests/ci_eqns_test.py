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
# Set Up by loading the required one and two electron integrals
# We will test on both H2 in sto-3g and 6 site Hubbard with U/t = 2.0
# 

def test_cc_residuals_zeroT():

    # 
    # Test to check the Coupled Cluster Residuals in the Zero Temperature Limit
    # 

    # Load the Integrals first
    iops = IOps(inp_file='TestInput')
    h1, eri_in, attrs = iops.loadHDF()
    nso = iops.nso

    # Next do the IntTran2
    hdiag, evecs = IntTran2(h1)
    eri = IntTran4(eri_in,evecs)

    # Construct random t1 and t2
    t1 = np.random.rand(nso,nso)
    t2 = eri_in*0.1

    # Other parameters needed for residuals
    tau = np.random.rand()
    mu = np.randon.rand()
    x = 1.0/np.sqrt(1 + np.exp(mu - beta*hdiag))
    y = np.sqrt( 1-x**2 )

    # Get the residuals
    r0, r1, r2 = ThermalCCSD(hdiag, eri, t1, t2, x, y)

    # First check the shapes of r0, r1, r2
    assert np.shape(r0) == (1,)
    assert np.shape(r1) == (nso,nso)
    assert np.shape(r2) == (nso,nso,nso,nso)

    # Compare with expected values (generated from Tom's MasterCode or my CC_SpinOrbital code)
