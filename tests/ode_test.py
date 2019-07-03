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

