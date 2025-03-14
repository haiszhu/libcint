#!/usr/bin/env python

import numpy as np
# from pyscf import gto, ao2mo
from pyscf import gto
import h5py
from scipy.io import savemat

geom='''
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561
        '''
        
mol_h2o_dimer = gto.M(
        verbose=7, 
        atom=geom, 
        basis='ccpvtz')

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_h2o_dimer.intor('int2e')
print(eri)
print(eri.shape)

# 
x, y, z = np.meshgrid(np.linspace(0,1,5),np.linspace(0,1,5),np.linspace(0,1,5),
                         indexing='ij')
xyz = np.column_stack([x.flatten(order='F'),y.flatten(order='F'),z.flatten(order='F')])

# Evaluate GTOs at the specified points
vals = np.array(mol_h2o_dimer.eval_gto('GTOval_sph',xyz))
print(vals.shape)
savemat("h2o_dimer_basis_check.mat", {'x': x, 'y': y, 'z': z, 'xyz': xyz, 'vals': vals})

# Calculate ERI for a given molecule and a AO basis: 
# eri has dimenions (nao, nao, nao, nao) where nao is the number of AO basis functions
eri = mol_h2o_dimer.intor('int2e')
print(eri)
print(eri.shape)

# Read ERI from Hai
with h5py.File('ERI_h2o_dimer_ccpvtz_1e-3.h5', 'r') as ar:
    eri0 = ar['DS1'][()]

diff = eri - eri0
print("Maximum diff between two ERIs:")
print(np.max(np.abs(diff)))
