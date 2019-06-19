import pytest
import numpy as np

from ../LibFortCC import *

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

def test_cc_residuals():

    # 
    # Test to check the Coupled Cluster Residuals in the Zero Temperature Limit
    # 

    f1
