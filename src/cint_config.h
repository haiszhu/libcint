#include "cint.h"

#define HAVE_EXPL
#define HAVE_SQRTL
/* #undef HAVE_FABSL */

#define HAVE_QUADMATH_H

#define WITH_RANGE_COULOMB

#ifndef M_PI
#define M_PI            3.1415926535897932384626433832795028
#endif
#define SQRTPI          1.7724538509055160272981674833411451

// ng[*]
#define IINC            0
#define JINC            1
#define KINC            2
#define LINC            3
#define GSHIFT          4
#define POS_E1          5
#define POS_E2          6
#define SLOT_RYS_ROOTS  6
#define TENSOR          7

#define EXPCUTOFF       60
#ifndef MIN_EXPCUTOFF
// ~ 1e-15
#define MIN_EXPCUTOFF   40
#endif

#define OF_CMPLX        2

#define GRID_BLKSIZE    104
