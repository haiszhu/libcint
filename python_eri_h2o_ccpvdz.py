#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py

# mol = gto.M(atom='O 0 0 0; H  0.757 0.587 0.000; H -0.757 0.587 0.000', basis='ccpvdz', verbose=7)
mol = gto.M(
        verbose = 0,
        atom = '''
        o    0    0.       0.
        h    0    -0.757   0.587
        h    0    0.757    0.587''',
        basis = 'ccpvdz') 

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol.intor('int2e')
print(eri)
print(eri.shape)

# Read ERI from Hai
with h5py.File('ERI_h2o_ccpvdz.h5', 'r') as ar:
    eri0 = ar['DS1'][()]

diff = eri - eri0
print("Maximum diff between two ERIs:")
print(np.max(np.abs(diff)))