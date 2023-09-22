#ifndef MATRIX_CHM_COMMON_H
#define MATRIX_CHM_COMMON_H

#include <Rinternals.h>
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"

#define AS_CHM_DN(x) \
as_cholmod_dense  ((CHM_DN) alloca(sizeof(cholmod_dense  )), x)
/* always checking versions : */
#define AS_CHM_FR(x) \
as_cholmod_factor3((CHM_FR) alloca(sizeof(cholmod_factor )), x,  TRUE)
#define AS_CHM_SP(x) \
as_cholmod_sparse ((CHM_SP) alloca(sizeof(cholmod_sparse )), x,  TRUE, FALSE)
#define AS_CHM_TR(x) \
as_cholmod_triplet((CHM_TR) alloca(sizeof(cholmod_triplet)), x,  TRUE)
/* never checking versions : */
#define AS_CHM_FR__(x) \
as_cholmod_factor3((CHM_FR) alloca(sizeof(cholmod_factor )), x, FALSE)
#define AS_CHM_SP__(x) \
as_cholmod_sparse ((CHM_SP) alloca(sizeof(cholmod_sparse )), x, FALSE, FALSE)
#define AS_CHM_TR__(x) \
as_cholmod_triplet((CHM_TR) alloca(sizeof(cholmod_triplet)), x, FALSE)

typedef       cholmod_common *       CHM_CM;
typedef       cholmod_dense  *       CHM_DN;
typedef const cholmod_dense  * const_CHM_DN;
typedef       cholmod_factor *       CHM_FR;
typedef const cholmod_factor * const_CHM_FR;
typedef       cholmod_sparse *       CHM_SP;
typedef const cholmod_sparse * const_CHM_SP;
typedef       cholmod_triplet*       CHM_TR;
typedef const cholmod_triplet* const_CHM_TR;

extern cholmod_common c ; /* structure for              int routines */
extern cholmod_common cl; /* structure for SuiteSparse_long routines */

/* NB: Versions of these are *EXPORTED* via ../inst/include/Matrix.h : */

SEXP get_SuiteSparse_version(void);

int R_cholmod_start(CHM_CM);

SEXP CHM_set_common_env(SEXP);
void CHM_store_common(void);
void CHM_restore_common(void);

CHM_SP as_cholmod_sparse   (CHM_SP, SEXP, Rboolean, Rboolean);
CHM_TR as_cholmod_triplet  (CHM_TR, SEXP, Rboolean);
CHM_DN as_cholmod_dense    (CHM_DN, SEXP);
CHM_FR as_cholmod_factor   (CHM_FR, SEXP);
CHM_FR as_cholmod_factor3  (CHM_FR, SEXP, Rboolean);
CHM_DN numeric_as_chm_dense(CHM_DN, double *, int, int);

SEXP chm_factor_to_SEXP (CHM_FR, int);
SEXP chm_sparse_to_SEXP (CHM_SP, int, int, int, const char *, SEXP);
SEXP chm_triplet_to_SEXP(CHM_TR, int, int, int, const char *, SEXP);

#endif /* MATRIX_CHM_COMMON_H */
