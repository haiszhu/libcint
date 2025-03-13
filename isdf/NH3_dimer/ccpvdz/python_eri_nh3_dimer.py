#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py
from scipy.io import savemat

geom='''
        N  -1.578718  -0.046611   0.000000 
        H  -2.158621   0.136396  -0.809565
        H  -2.158621   0.136396   0.809565
        H  -0.849471   0.658193   0.000000
        N   1.578718   0.046611   0.000000
        H   2.158621  -0.136396  -0.809565
        H   0.849471  -0.658193   0.000000
        H   2.158621  -0.136396   0.809565
        '''

mol_nh3_dimer = gto.M(
        verbose=7, 
        atom=geom, 
        basis='ccpvdz')

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_nh3_dimer.intor('int2e')
print(eri)
print(eri.shape)

# 
x, y, z = np.meshgrid(np.linspace(0,1,5),np.linspace(0,1,5),np.linspace(0,1,5),
                         indexing='ij')
xyz = np.column_stack([x.flatten(order='F'),y.flatten(order='F'),z.flatten(order='F')])

# Evaluate GTOs at the specified points
vals = np.array(mol_nh3_dimer.eval_gto('GTOval_sph',xyz))
print(vals.shape)
savemat("nh3_dimer_basis_check.mat", {'x': x, 'y': y, 'z': z, 'xyz': xyz, 'vals': vals})

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_nh3_dimer.intor('int2e')
print(eri)
print(eri.shape)

# Read ERI from Hai
with h5py.File('ERI_nh3_dimer_ccpvdz_1e-3.h5', 'r') as ar:
    eri0 = ar['DS1'][()]

diff = eri - eri0
print("Maximum diff between two ERIs:")
print(np.max(np.abs(diff)))
