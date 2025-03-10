/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Parameters and function signature for libcint.
 */

#define CINT_VERSION            "6.1.2"
#define CINT_SOVERSION          6

/* #undef I8 */
#ifdef I8
#include <stdint.h>
#define FINT int64_t
#else
#define FINT int
#endif

/* #undef CACHE_SIZE_I8 */
#ifdef CACHE_SIZE_I8
#include <stdint.h>
#define CACHE_SIZE_T int64_t
#else
#define CACHE_SIZE_T FINT
#endif

// global parameters in env
// Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
#define PTR_EXPCUTOFF           0
// R_C of (r-R_C) in dipole, GIAO operators
#define PTR_COMMON_ORIG         1
// R_O in 1/|r-R_O|
#define PTR_RINV_ORIG           4
// ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
#define PTR_RINV_ZETA           7
// omega parameter in range-separated coulomb operator
// LR interaction: erf(omega*r12)/r12 if omega > 0
// SR interaction: erfc(omega*r12)/r12 if omega < 0
#define PTR_RANGE_OMEGA         8
// Yukawa potential and Slater-type geminal e^{-zeta r}
#define PTR_F12_ZETA            9
// Gaussian type geminal e^{-zeta r^2}
#define PTR_GTG_ZETA            10
#define NGRIDS                  11
#define PTR_GRIDS               12
#define PTR_ENV_START           20


// slots of atm
#define CHARGE_OF       0
#define PTR_COORD       1
#define NUC_MOD_OF      2
#define PTR_ZETA        3
#define PTR_FRAC_CHARGE 4
#define RESERVE_ATMSLOT 5
#define ATM_SLOTS       6


// slots of bas
#define ATOM_OF         0
#define ANG_OF          1
#define NPRIM_OF        2
#define NCTR_OF         3
#define KAPPA_OF        4
#define PTR_EXP         5
#define PTR_COEFF       6
#define RESERVE_BASLOT  7
#define BAS_SLOTS       8

// slots of gout
#define POSX            0
#define POSY            1
#define POSZ            2
#define POS1            3
// For 2-electron integral with two spin operators
// SIGMA1X * SIGMA2X     0
// SIGMA1Y * SIGMA2X     1
// SIGMA1Z * SIGMA2X     2
// I1_2x2  * SIGMA2X     3
// SIGMA1X * SIGMA2Y     4
// SIGMA1Y * SIGMA2Y     5
// SIGMA1Z * SIGMA2Y     6
// I1_2x2  * SIGMA2Y     7
// SIGMA1X * SIGMA2Z     8
// SIGMA1Y * SIGMA2Z     9
// SIGMA1Z * SIGMA2Z     10
// I1_2x2  * SIGMA2Z     11
// SIGMA1X * I2_2x2      12
// SIGMA1Y * I2_2x2      13
// SIGMA1Z * I2_2x2      14
// I1_2x2  * I2_2x2      15
#define POSXX           0
#define POSYX           1
#define POSZX           2
#define POS1X           3
#define POSXY           4
#define POSYY           5
#define POSZY           6
#define POS1Y           7
#define POSXZ           8
#define POSYZ           9
#define POSZZ           10
#define POS1Z           11
#define POSX1           12
#define POSY1           13
#define POSZ1           14
#define POS11           15

// tensor
#define TSRX        0
#define TSRY        1
#define TSRZ        2
#define TSRXX       0
#define TSRXY       1
#define TSRXZ       2
#define TSRYX       3
#define TSRYY       4
#define TSRYZ       5
#define TSRZX       6
#define TSRZY       7
#define TSRZZ       8

// other boundaries
#define MXRYSROOTS      32 // > ANG_MAX*2+1 for 4c2e
#define ANG_MAX         15 // l = 0..15
#define LMAX1           16 // > ANG_MAX
#define CART_MAX        136 // > (ANG_MAX*(ANG_MAX+1)/2)
#define SHLS_MAX        1048576
#define NPRIM_MAX       64
#define NCTR_MAX        64

#define POINT_NUC       1
#define GAUSSIAN_NUC    2
#define FRAC_CHARGE_NUC 3

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    double rij[3];
    double eij;
    double cceij;
} PairData;
typedef struct {
    FINT **index_xyz_array; // LMAX1**4 pointers to index_xyz
    FINT **non0ctr;
    FINT **sortedidx;
    FINT nbas;
    double **log_max_coeff;
    PairData **pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
} CINTOpt;

// Add this macro def to make pyscf compatible with both v4 and v5
#define HAVE_DEFINED_CINTENVVARS_H
typedef struct {
        FINT *atm;
        FINT *bas;
        double *env;
        FINT *shls;
        FINT natm;
        FINT nbas;

        FINT i_l;
        FINT j_l;
        FINT k_l;
        FINT l_l;
        FINT nfi;  // number of cartesian components
        FINT nfj;
        // in int1e_grids, the grids_offset and the number of grids
        union {FINT nfk; FINT grids_offset;};
        union {FINT nfl; FINT ngrids;};
        FINT nf;  // = nfi*nfj*nfk*nfl;
        FINT rys_order; // = nrys_roots for regular ERIs. can be nrys_roots/2 for SR ERIs
        FINT x_ctr[4];

        FINT gbits;
        FINT ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        FINT ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint.h
        FINT ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        FINT li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        FINT lj_ceil;
        FINT lk_ceil;
        FINT ll_ceil;
        FINT g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        FINT g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        FINT g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        FINT g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        FINT nrys_roots;
        FINT g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        FINT g2d_ijmax;
        FINT g2d_klmax;
        double common_factor;
        double expcutoff;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        double *rk;
        // in int2e or int3c2e, the coordinates of the fourth shell
        // in int1e_grids, the pointer for the grids coordinates
        union {double *rl; double *grids;};

        FINT (*f_g0_2e)();
        void (*f_g0_2d4d)();
        void (*f_gout)();
        CINTOpt *opt;

        /* values are assigned during calculation */
        int *idx;
        double ai[1];
        double aj[1];
        double ak[1];
        double al[1];
        double fac[1];
        double rij[3];
        double rkl[3];
} CINTEnvVars;
#endif

FINT CINTlen_cart(const FINT l);
FINT CINTlen_spinor(const FINT bas_id, const FINT *bas);

FINT CINTcgtos_cart(const FINT bas_id, const FINT *bas);
FINT CINTcgtos_spheric(const FINT bas_id, const FINT *bas);
FINT CINTcgtos_spinor(const FINT bas_id, const FINT *bas);
FINT CINTcgto_cart(const FINT bas_id, const FINT *bas);
FINT CINTcgto_spheric(const FINT bas_id, const FINT *bas);
FINT CINTcgto_spinor(const FINT bas_id, const FINT *bas);

FINT CINTtot_pgto_spheric(const FINT *bas, const FINT nbas);
FINT CINTtot_pgto_spinor(const FINT *bas, const FINT nbas);

FINT CINTtot_cgto_cart(const FINT *bas, const FINT nbas);
FINT CINTtot_cgto_spheric(const FINT *bas, const FINT nbas);
FINT CINTtot_cgto_spinor(const FINT *bas, const FINT nbas);

void CINTshells_cart_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);
void CINTshells_spheric_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);
void CINTshells_spinor_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);

double *CINTc2s_bra_sph(double *sph, FINT nket, double *cart, FINT l);
double *CINTc2s_ket_sph(double *sph, FINT nket, double *cart, FINT l);
double *CINTc2s_ket_sph1(double *sph, double *cart, FINT lds, FINT ldc, FINT l);


double CINTgto_norm(FINT n, double a);


void CINTinit_2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);


FINT cint2e_cart(double *opijkl, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                CINTOpt *opt);
void cint2e_cart_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
FINT cint2e_sph(double *opijkl, FINT *shls,
               FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
               CINTOpt *opt);
void cint2e_sph_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env);
FINT cint2e(double *opijkl, FINT *shls,
           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
           CINTOpt *opt);
void cint2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env);

#ifndef __cplusplus
#include <complex.h>

void CINTc2s_ket_spinor_sf1(double complex *gspa, double complex *gspb, double *gcart,
                            FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_iket_spinor_sf1(double complex *gspa, double complex *gspb, double *gcart,
                             FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_ket_spinor_si1(double complex *gspa, double complex *gspb, double *gcart,
                            FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_iket_spinor_si1(double complex *gspa, double complex *gspb, double *gcart,
                             FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
#endif
