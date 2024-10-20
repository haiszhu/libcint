#ifndef CGTO_H
#define CGTO_H

#include <stdint.h>

void GTOval_sph(int ngrids, int *shls_slice, int *ao_loc,
                double *ao, double *coord, uint8_t *non0table,
                int *atm, int natm, int *bas, int nbas, double *env);

void GTOval_sph_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

void GTOval_h2o_ccpvtz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0table,
                      int *atm, int *natm, int *bas, int *nbas, double *env);

#endif // CGTO_H
