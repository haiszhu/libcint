#!/usr/bin/env python
# 
# this is fine...
# 

import numpy
import scipy.special
from pyscf import gto
from pyscf.gto.mole import format_atom

ELEMENTS = [
    'X',  # Ghost
    'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca',
    'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
]

ELEMENTS_PROTON = {}
for i, symb in enumerate(ELEMENTS):
    # e.g. ELEMENTS_PROTON["C"] = 6
    ELEMENTS_PROTON[symb] = i
    # Also store uppercase version
    ELEMENTS_PROTON[symb.upper()] = i
    
BLKSIZE = 56

# for _atm, _bas, _env
CHARGE_OF  = 0
PTR_COORD  = 1
NUC_MOD_OF = 2
PTR_ZETA   = 3
PTR_FRAC_CHARGE = 4
ATM_SLOTS  = 6
ATOM_OF    = 0
ANG_OF     = 1
NPRIM_OF   = 2
NCTR_OF    = 3
RADI_POWER = 3 # for ECP
KAPPA_OF   = 4
SO_TYPE_OF = 4 # for ECP
PTR_EXP    = 5
PTR_COEFF  = 6
BAS_SLOTS  = 8
# pointer to env
PTR_EXPCUTOFF   = 0
PTR_COMMON_ORIG = 1
PTR_RINV_ORIG   = 4
PTR_RINV_ZETA   = 7
PTR_RANGE_OMEGA = 8
PTR_F12_ZETA    = 9
PTR_GTG_ZETA    = 10
NGRIDS          = 11
PTR_GRIDS       = 12
AS_RINV_ORIG_ATOM = 17
AS_ECPBAS_OFFSET = 18
AS_NECPBAS      = 19
PTR_ENV_START   = 20
# parameters from libcint
NUC_POINT = 1
NUC_GAUSS = 2
# nucleus with fractional charges. It can be used to mimic MM particles
NUC_FRAC_CHARGE = 3
NUC_ECP = 4  # atoms with pseudo potential

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
mol = gto.M(
    verbose=7, 
    atom=geom, 
    basis='ccpvdz')

atm = numpy.asarray(mol._atm, dtype=numpy.int32, order='C')
bas = numpy.asarray(mol._bas, dtype=numpy.int32, order='C')
env = numpy.asarray(mol._env, dtype=numpy.double, order='C')
natm = atm.shape[0]
nbas = bas.shape[0]

atoms_list = []
for line in geom.strip().splitlines():
    fields = line.split()
    symb = fields[0]
    x = float(fields[1])
    y = float(fields[2])
    z = float(fields[3])
    atoms_list.append((symb, (x,y,z)))
atoms_list_bohr = format_atom(atoms_list, unit='Angstrom')

unique_elems = set(a[0] for a in atoms_list_bohr)
basis_dict = {}
for elem in unique_elems:
    basis_dict[elem] = gto.basis.load('ccpvdz', elem)
    
pre_env=numpy.zeros(PTR_ENV_START, dtype=float)

# 
_atm = []
_bas = []
_env = [pre_env]
ptr_env = len(pre_env)

for ia, atom in enumerate(atoms_list_bohr):
    symb = atom[0]

    nuclear_model = NUC_POINT
    # atm0, env0 = make_atm_env(atom, ptr_env, nuclear_model, prop)
    
    # make my version
    nuc_charge = ELEMENTS_PROTON[atom[0]]
    zeta = 0
    env0 = numpy.hstack((atom[1], zeta))
    atm0 = numpy.zeros(6, dtype=numpy.int32)        
    atm0[CHARGE_OF] = nuc_charge
    atm0[PTR_COORD] = ptr_env
    atm0[NUC_MOD_OF] = nuclear_model       
    atm0[PTR_ZETA ] = ptr_env + 3
    
    ptr_env = ptr_env + len(env0)
    _atm.append(atm0)
    _env.append(env0)

_basdic = {}
for symb, basis_add in basis_dict.items():
    
    # make my version
    ptr_env0 = ptr_env
    atom_id = 0
    bas0 = []
    env0 = []
    for b in basis_add:
        angl = b[0]
        kappa = 0
        b_coeff = numpy.array(sorted(list(b[1:]), reverse=True))
        es = b_coeff[:,0]
        cs = b_coeff[:,1:]
        nprim, nctr = cs.shape
        # cs = numpy.einsum('pi,p->pi', cs, gto_norm(angl, es))
        # my version
        n1 = (angl*2+2 + 1) * .5
        gaussian_intangles = scipy.special.gamma(n1) / (2. * (2*es)**n1)
        gto_normangles = 1/numpy.sqrt(gaussian_intangles)
        cs = numpy.einsum('pi,p->pi', cs, gto_normangles)
        # if NORMALIZE_GTO: # always true
        #     cs = _nomalize_contracted_ao(angl, es, cs)
        
        ee = es.reshape(-1,1) + es.reshape(1,-1)
        # ee = gaussian_int(angl*2+2, ee)
        # r'''int_0^inf x^n exp(-alpha x^2) dx'''
        n1 = (angl*2+2 + 1) * .5
        ee = scipy.special.gamma(n1) / (2. * ee**n1)
        s1 = 1. / numpy.sqrt(numpy.einsum('pi,pq,qi->i', cs, ee, cs))
        cs = numpy.einsum('pi,i->pi', cs, s1)
        
        env0.append(es)
        env0.append(cs.T.reshape(-1))
        ptr_exp = ptr_env0
        ptr_coeff = ptr_exp + nprim
        ptr_env0 = ptr_coeff + nprim * nctr
        bas0.append([atom_id, angl, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0])
        
    env0 = numpy.hstack(env0)
    bas0 = numpy.array(bas0, numpy.int32).reshape(-1,BAS_SLOTS)
    env0 = numpy.array(env0, numpy.double)
    
    ptr_env = ptr_env + len(env0)
    _basdic[symb] = bas0
    _env.append(env0)

for ia, atom in enumerate(atoms_list_bohr):
    symb = atom[0]
    b = _basdic[symb].copy()
    b[:,ATOM_OF] = ia
    _bas.append(b)

_atm = numpy.asarray(numpy.vstack(_atm), numpy.int32).reshape(-1, ATM_SLOTS)
_bas = numpy.asarray(numpy.vstack(_bas), numpy.int32).reshape(-1, BAS_SLOTS)
_env = numpy.asarray(numpy.hstack(_env), dtype=numpy.float64)

print("==> _atm =\n", _atm - mol._atm)
print("==> _bas =\n", _bas - mol._bas)
print("==> _env =\n", _env - mol._env)