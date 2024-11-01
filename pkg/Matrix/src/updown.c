/* C implementation of methods for updown, update */

#include "cholmod-etc.h"
#include "Mdefines.h"
#include "M5.h"

SEXP sparseCholesky_updown(SEXP s_trf, SEXP s_obj, SEXP s_update)
{
	cholmod_factor *L = M2CHF(s_trf, 1);
	cholmod_sparse *A = M2CHS(s_obj, 1);

	L = cholmod_copy_factor(L, &c);
	cholmod_updown(Rf_asLogical(s_update) != 0, A, L, &c);

#define UPDOWN_FINISH \
	do { \
		SEXP dimnames = PROTECT(DIMNAMES(s_trf, 0)); \
		PROTECT(s_trf = CHF2M(L, 1)); \
		cholmod_free_factor(&L, &c); \
		if (TYPEOF(s_trf) == CHARSXP) \
			Rf_error("%s", CHAR(s_trf)); \
		SET_DIMNAMES(s_trf, 0, dimnames); \
		UNPROTECT(2); \
	} while (0)

	UPDOWN_FINISH;
	return s_trf;
}

SEXP sparseCholesky_update(SEXP s_trf, SEXP s_obj, SEXP s_beta)
{
	Rcomplex beta = Rf_asComplex(s_beta);
	if (!R_FINITE(beta.r) || !R_FINITE(beta.i))
		Rf_error(_("'%s' is not a number or not finite"), "beta");

	cholmod_factor *L = M2CHF(s_trf, 1);
	cholmod_sparse *A = M2CHS(s_obj, 1);
	double b[2]; b[0] = beta.r; b[1] = beta.i;

	/* defined in ./objects.c : */
	char Matrix_shape(SEXP);
	if (Matrix_shape(s_obj) == 's')
		A->stype = (UPLO(s_obj) == 'U') ? 1 : -1;

	L = cholmod_copy_factor(L, &c);

	c.final_asis = 0;
	c.final_ll = L->is_ll;
	c.final_super = L->is_super;
	c.final_pack = 1;
	c.final_monotonic = 1;

	cholmod_factorize_p(A, b, NULL, 0, L, &c);
	cholmod_defaults(&c);

	UPDOWN_FINISH;
	return s_trf;
}
