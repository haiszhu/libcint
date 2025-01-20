ISDF + Adapted Grid
-------------------
**Last Updated:** Dec 19, 2024

To demonstrate the applicability of the ISDF+adapted grid approach, 
we examine its convergence with respect to: 
1. the sizes of the basis sets
2. the chemical species: Consider heavier elements
3. the sizes of molecules

# Benchmark systems
To examine the basis set effects, we consider the following molecules 
from the S22 benchmark dataset:
1. NH3 dimer
2. H2O dimer
3. Uracil dimer 
4. ...
5. ...

For each system, the basis set is systematically enlarged using  
1. ccvpdz
2. aug-ccpvdz
3. ccpvtz
4. aug-ccpvtz

To examine the performance of ISDF+adapted grid as the molecule size increasese, 
we consider linear alkane C_{n}H_{2n+2} with n=1~10
