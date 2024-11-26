#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py

mol = gto.M(
        verbose = 0,
        atom = '''
        o    0    0.       0.
        h    0    -0.757   0.587
        h    0    0.757    0.587''',
        basis = 'ccpvdz') 

# Read src
with h5py.File('src_h2o_ccpvdz.h5', 'r') as ar:
    src = ar['DS1'][()]
    
print(src.shape)