import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio  


molname = 'h2o_dimer'
molname = 'nh3_dimer'

# basmod = 'cc-pvdz.dat'
# basmod = 'aug-cc-pvdz.dat'
# basmod = 'cc-pvtz.dat'
basmod = 'aug-cc-pvtz.dat'

# 
basname = basmod.replace('.dat','')
mat_filename = f'{molname}_{basname}.mat'

# 
data = sio.loadmat(mat_filename)
epsvals = data['epsvals'].flatten()
relerrs = data['relerrs'].flatten()
relerr2s = data['relerr2s'].flatten()
relerr2ups = data['relerr2ups'].flatten()
epsvals = epsvals[0:-1]
relerrs = relerrs[0:-1]
relerr2s = relerr2s[0:-1]
relerr2ups = relerr2ups[0:-1]

# 
plt.figure(figsize=(6,4))
# plt.loglog(epsvals, relerrs, 'o-', linewidth=2, label=r'$\max_{\tilde{T}} E_{\infty}^{(\phi_i)}$ on order $k$ grid')
# plt.loglog(epsvals, relerr2s, '^-', linewidth=2, label=r'$\max_{\tilde{T}} E_{\infty}^{(\rho_{ij})}$ on order $k$ grid')
# plt.loglog(epsvals, relerr2ups, 's-', linewidth=2, label=r'$\max_{T} E_{\infty}^{(\rho_{ij})}$ on order $1.5\times k$ grid')
plt.loglog(epsvals, relerrs, 'o-', linewidth=2, label=r'$\max_{i}\ E^{(\phi_i)}$ on $\tilde{T}$')
plt.loglog(epsvals, relerr2s, '^-', linewidth=2, label=r'$\max_{ij}\ E^{(\rho_{ij})}$ on $\tilde{T}$')
plt.loglog(epsvals, relerr2ups, 's-', linewidth=2, label=r'$\max_{ij}\ E^{(\rho_{ij})}$ on $T$')


# legend
plt.legend(fontsize=10, loc='best')

# xlabel ylabel
plt.xlabel(r'$\varepsilon_{adpt}$ ', fontsize=14)
# plt.ylabel(r'$E_{\infty}$', fontsize=14)
plt.ylabel(r'$E$', fontsize=14)

# 
title_dict = {
    ('nh3_dimer', 'aug-cc-pvtz.dat'): r'NH$_3$ dimer, aug-cc-pvtz',
    ('nh3_dimer', 'cc-pvtz.dat'):    r'NH$_3$ dimer, cc-pvtz',
    ('nh3_dimer', 'aug-cc-pvdz.dat'):r'NH$_3$ dimer, aug-cc-pvdz',
    ('nh3_dimer', 'cc-pvdz.dat'):    r'NH$_3$ dimer, cc-pvdz',
    ('h2o_dimer', 'aug-cc-pvtz.dat'):r'H$_2$O dimer, aug-cc-pvtz',
    ('h2o_dimer', 'cc-pvtz.dat'):    r'H$_2$O dimer, cc-pvtz',
    ('h2o_dimer', 'aug-cc-pvdz.dat'):r'H$_2$O dimer, aug-cc-pvdz',
    ('h2o_dimer', 'cc-pvdz.dat'):    r'H$_2$O dimer, cc-pvdz'
}

title_str = title_dict.get( (molname, basmod), f'{molname}, {basmod}')
plt.title(title_str, fontsize=14)

# 
plt.grid(True)

# 
plt.tight_layout()
plt.tight_layout()
plt.savefig(f'{molname}_{basname}.png', dpi=600, bbox_inches='tight')
plt.show()
