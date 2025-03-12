#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py
from scipy.io import savemat

geom='''
        O  0    0.       0.
        H  0    -0.757   0.587
        H  0    0.757    0.587
        '''
        
mol_h2o = gto.M(
        verbose=7, 
        atom=geom, 
        basis='ccpvtz')

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_h2o.intor('int2e')
print(eri.shape)

# Read ERI from Hai
with h5py.File('ERI_h2o_ccpvtz_1e-3.h5', 'r') as ar:
    eri0 = ar['DS1'][()]

diff = eri - eri0
print("Maximum diff between two ERIs:")
print(np.max(np.abs(diff)))
