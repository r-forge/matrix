#include "lmer.h"
#include "pedigree.h"
#include "Matrix.h"
#include <R_ext/Rdynload.h>

static R_CallMethodDef CallEntries[] = {
    {"glmer_MCMCsamp", (DL_FUNC) &glmer_MCMCsamp, 5},
    {"glmer_PQL", (DL_FUNC) &glmer_PQL, 1},
    {"glmer_devLaplace", (DL_FUNC) &glmer_devLaplace, 2},
    {"glmer_finalize", (DL_FUNC) &glmer_finalize, 1},
    {"glmer_init", (DL_FUNC) &glmer_init, 1},
    {"lme4_rWishart", (DL_FUNC) &lme4_rWishart, 3},
    {"mer_ECMEsteps", (DL_FUNC) &mer_ECMEsteps, 3},
    {"mer_MCMCsamp", (DL_FUNC) &mer_MCMCsamp, 5},
    {"mer_coef", (DL_FUNC) &mer_coef, 2},
    {"mer_coefGets", (DL_FUNC) &mer_coefGets, 3},
    {"mer_create", (DL_FUNC) &mer_create, 10},
    {"mer_dtCMatrix", (DL_FUNC) &mer_dtCMatrix, 1},
    {"mer_dtCMatrix_inv", (DL_FUNC) &mer_dtCMatrix_inv, 1},
    {"mer_factor", (DL_FUNC) &mer_factor, 1},
    {"mer_fitted", (DL_FUNC) &mer_fitted, 1},
    {"mer_fixef", (DL_FUNC) &mer_fixef, 1},
    {"mer_gradComp", (DL_FUNC) &mer_gradComp, 1},
    {"mer_gradient", (DL_FUNC) &mer_gradient, 2},
    {"mer_hat_trace", (DL_FUNC) &mer_hat_trace, 1},
    {"mer_hat_trace2", (DL_FUNC) &mer_hat_trace2, 1},
    {"mer_initial", (DL_FUNC) &mer_initial, 1},
    {"mer_isNested", (DL_FUNC) &mer_isNested, 1},
    {"mer_postVar", (DL_FUNC) &mer_postVar, 1},
    {"mer_ranef", (DL_FUNC) &mer_ranef, 1},
    {"mer_secondary", (DL_FUNC) &mer_secondary, 1},
    {"mer_sigma", (DL_FUNC) &mer_sigma, 2},
    {"mer_simulate", (DL_FUNC) &mer_simulate, 2},
    {"mer_update_ZXy", (DL_FUNC) &mer_update_ZXy, 1},
    {"mer_update_y", (DL_FUNC) &mer_update_y, 2},
    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},
    {"Zt_create", (DL_FUNC) &Zt_create, 2},
    {NULL, NULL, 0}
};

cholmod_common c;

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    /* Get pointers to functions defined in the Matrix package */
    M_alloc_dgeMatrix = (SEXP(*)(int,int,SEXP,SEXP))
	R_GetCCallable("Matrix", "alloc_dgeMatrix");
    M_alloc_dpoMatrix = (SEXP(*)(int,char*,SEXP,SEXP))
	R_GetCCallable("Matrix", "alloc_dpoMatrix");
    M_alloc_dtrMatrix = (SEXP(*)(int,char*,char*,SEXP,SEXP))
	R_GetCCallable("Matrix", "alloc_dtrMatrix");
    M_alloc_dsCMatrix = (SEXP(*)(int,int,char*,SEXP,SEXP))
	R_GetCCallable("Matrix", "alloc_dsCMatrix");
    M_as_cholmod_dense = (cholmod_dense*(*)(SEXP))
	R_GetCCallable("Matrix", "as_cholmod_dense");
    M_as_cholmod_factor = (cholmod_factor*(*)(SEXP))
	R_GetCCallable("Matrix", "as_cholmod_factor");
    M_as_cholmod_sparse = (cholmod_sparse*(*)(SEXP))
	R_GetCCallable("Matrix", "as_cholmod_sparse");
    M_chm_factor_to_SEXP = (SEXP(*)(cholmod_factor*,int))
	R_GetCCallable("Matrix", "chm_factor_to_SEXP");
    M_cholmod_aat = (cholmod_sparse*(*)(cholmod_sparse*,int*,size_t,
					int,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_aat");
    M_cholmod_add = (cholmod_sparse*(*)(cholmod_sparse*,cholmod_sparse*,
					double*,double*,int,int,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_add");
    M_cholmod_allocate_dense = (cholmod_dense*(*)(size_t nrow, size_t ncol, size_t d,
						  int xtype, cholmod_common *Common))
	R_GetCCallable("Matrix", "cholmod_allocate_dense");
    M_cholmod_allocate_sparse = (cholmod_sparse*(*)
				 (size_t,size_t,size_t,int,int,int,int,
				  cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_allocate_sparse");
    M_cholmod_analyze = (cholmod_factor*(*)(cholmod_sparse*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_analyze");
    M_cholmod_copy = (cholmod_sparse*(*)(cholmod_sparse*,int,int,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_copy");
    M_cholmod_copy_dense = (cholmod_dense*(*)(cholmod_dense*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_copy_dense");
    M_cholmod_copy_factor = (cholmod_factor*(*)(cholmod_factor*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_copy_factor");
    M_cholmod_copy_sparse = (cholmod_sparse*(*)(cholmod_sparse*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_copy_sparse");
    M_cholmod_factor_to_sparse = (cholmod_sparse*(*)(cholmod_factor*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_factor_to_sparse");
    M_cholmod_factorize = (int(*)(cholmod_sparse*,cholmod_factor*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_factorize");
    M_cholmod_finish = (int(*)(cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_finish");
    M_cholmod_free_dense = (int(*)(cholmod_dense**,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_free_dense");
    M_cholmod_free_factor = (int(*)(cholmod_factor**,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_free_factor");
    M_cholmod_free_sparse = (int(*)(cholmod_sparse**,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_free_sparse");
    M_cholmod_nnz = (long(*)(cholmod_sparse*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_nnz");
    M_cholmod_sdmult = (int(*)(cholmod_sparse*,int,double*,double*,
			       cholmod_dense*,cholmod_dense*,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_sdmult");
    M_cholmod_solve = (cholmod_dense*(*)(int,cholmod_factor*,cholmod_dense*,
					 cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_solve");
    M_cholmod_speye = (cholmod_sparse*(*)(size_t,size_t,int,cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_speye");
    M_cholmod_spsolve = (cholmod_sparse*(*)(int,cholmod_factor*,
					    cholmod_sparse*, cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_spsolve");
    M_cholmod_start = (int(*)(cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_start");
    M_cholmod_transpose = (cholmod_sparse*(*)(cholmod_sparse*,int,
					      cholmod_common*))
	R_GetCCallable("Matrix", "cholmod_transpose");
    M_dpoMatrix_chol = (SEXP(*)(SEXP))
	R_GetCCallable("Matrix", "dpoMatrix_chol");
    M_numeric_as_chm_dense = (cholmod_dense*(*)(double*,int))
	R_GetCCallable("Matrix", "numeric_as_chm_dense");

    M_cholmod_start(&c);
    lme4_DSym = install("D");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_OmegaSym = install("Omega");
    lme4_RXXSym = install("RXX");
    lme4_RZXSym = install("RZX");
    lme4_RZXinvSym = install("RZXinv");
    lme4_XSym = install("X");
    lme4_XtXSym = install("XtX");
    lme4_XtySym = install("Xty");
    lme4_ZZpOSym = install("ZZpO");
    lme4_ZtSym = install("Zt");
    lme4_ZtXSym = install("ZtX");
    lme4_ZtZSym = install("ZtZ");
    lme4_ZtySym = install("Zty");
    lme4_bVarSym = install("bVar");
    lme4_callSym = install("call");
    lme4_cnamesSym = install("cnames");
    lme4_devCompSym = install("devComp");
    lme4_devianceSym = install("deviance");
    lme4_diagSym = install("diag");
    lme4_factorSym = install("factor");
    lme4_familySym = install("family");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_gradCompSym = install("gradComp");
    lme4_iSym = install("i");
    lme4_methodSym = install("method");
    lme4_ncSym = install("nc");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_rXySym = install("rXy");
    lme4_rZySym = install("rZy");
    lme4_ranefSym = install("ranef");
    lme4_statusSym = install("status");
    lme4_uploSym = install("uplo");
    lme4_useScaleSym = install("useScale");
    lme4_wrkresSym = install("wrkres");
    lme4_wtsSym = install("wts");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
}

void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}
