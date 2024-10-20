#include "cgto.h"
#include <stdlib.h>

void GTOval_sph(int ngrids, int *shls_slice, int *ao_loc,
                double *ao, double *coord, uint8_t *non0table,
                int *atm, int natm, int *bas, int nbas, double *env);

void GTOval_sph_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    GTOval_sph(*ngrids, shls_slice, ao_loc, ao, coord, non0tab_uint8, atm, *natm, bas, *nbas, env);

}

void GTOval_h2o_ccpvtz_mwrap(int *ngrids, int *shls_slice, int *ao_loc,
                      double *ao, double *coord, int *non0tab,
                      int *atm, int *natm, int *bas, int *nbas, double *env) {
    uint8_t non0tab_uint8[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // TODO: change to uint8_t
    GTOval_sph(*ngrids, shls_slice, ao_loc, ao, coord, non0tab_uint8, atm, *natm, bas, *nbas, env);

}
