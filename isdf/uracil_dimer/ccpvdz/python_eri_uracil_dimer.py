#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py
from scipy.io import savemat

geom='''
    O    -1.4663316    1.0121693    0.0000000
    C    -0.6281464    1.9142678    0.0000000
    N     0.7205093    1.6882688    0.0000000
    C     1.6367290    2.7052764    0.0000000
    C     1.2769036    4.0061763    0.0000000
    C    -0.1286005    4.3621549    0.0000000
    N    -0.9777230    3.2396433    0.0000000
    O    -0.5972229    5.4864066    0.0000000
    H     2.0103504    4.7938642    0.0000000
    H     1.0232515    0.7061820    0.0000000
    H    -1.9700268    3.4323850    0.0000000
    H     2.6690620    2.3883417    0.0000000
    O     1.4663316   -1.0121693    0.0000000
    C     0.6281464   -1.9142678    0.0000000
    N    -0.7205093   -1.6882688    0.0000000
    C    -1.6367290   -2.7052764    0.0000000
    C    -1.2769036   -4.0061763    0.0000000
    C     0.1286005   -4.3621549    0.0000000
    N     0.9777230   -3.2396433    0.0000000
    O     0.5972229   -5.4864066    0.0000000
    H    -2.0103504   -4.7938642    0.0000000
    H    -1.0232515   -0.7061820    0.0000000
    H     1.9700268   -3.4323850    0.0000000
    H    -2.6690620   -2.3883417    0.0000000
    '''
        
mol_uracil_dimer = gto.M(
        verbose=7, 
        atom=geom, 
        basis='ccpvdz')

# 
x, y, z = np.meshgrid(np.linspace(0,1,5),np.linspace(0,1,5),np.linspace(0,1,5),
                         indexing='ij')
xyz = np.column_stack([x.flatten(order='F'),y.flatten(order='F'),z.flatten(order='F')])

# Evaluate GTOs at the specified points
vals = np.array(mol_uracil_dimer.eval_gto('GTOval_sph',xyz))
print(vals.shape)
savemat("uracil_dimer_basis_check.mat", {'x': x, 'y': y, 'z': z, 'xyz': xyz, 'vals': vals})

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_uracil_dimer.intor('int2e')
print(eri)
print(eri.shape)

# Read ERI from Hai
with h5py.File('ERI_uracil_dimer_ccpvdz_1e-02.h5', 'r') as ar:
    eri0 = ar['DS1'][()]

diff = eri - eri0
print("Maximum diff between two ERIs:")
print(np.max(np.abs(diff)))
