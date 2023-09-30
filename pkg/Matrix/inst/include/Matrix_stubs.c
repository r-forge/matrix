#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "Matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==== CHOLMOD stubs =============================================== */

CHM_SP attribute_hidden
M_cholmod_aat(CHM_SP A, int *fset, size_t fsize,int mode, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int *, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int *, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_aat");
	return fun(A, fset, fsize, mode, Common);
}

CHM_SP attribute_hidden
M_cholmod_add(CHM_SP A, CHM_SP B, double alpha[2], double beta[2],
              int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, double[2], double[2],
	                    int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, double[2], double[2],
		                 int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_add");
	return fun(A, B, alpha, beta, values, sorted, Common);
}

CHM_DN attribute_hidden
M_cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d,
                         int xtype, CHM_CM Common)
{
	static CHM_DN(*fun)(size_t, size_t, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(size_t, size_t, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_dense");
	return fun(nrow, ncol, d, xtype, Common);
}

CHM_SP attribute_hidden
M_cholmod_allocate_sparse(size_t nrow, size_t ncol, size_t nzmax,
                          int sorted, int packed, int stype, int xtype,
                          CHM_CM Common)
{
	static CHM_SP(*fun)(size_t, size_t, size_t, int, int, int, int,
	                    CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(size_t, size_t, size_t, int, int, int, int,
		                 CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_sparse");
	return fun(nrow, ncol, nzmax, sorted, packed, stype, xtype, Common);
}

CHM_TR attribute_hidden
M_cholmod_allocate_triplet(size_t nrow, size_t ncol, size_t nzmax,
                           int stype, int xtype, CHM_CM Common)
{
	static CHM_TR(*fun)(size_t, size_t, size_t, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_TR(*)(size_t, size_t, size_t, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_allocate_triplet");
	return fun(nrow, ncol, nzmax, stype, xtype, Common);
}

CHM_FR attribute_hidden
M_cholmod_analyze(CHM_SP A, CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_SP,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_analyze");
	return fun(A, Common);
}

CHM_FR attribute_hidden
M_cholmod_analyze_p(CHM_SP A, int *Perm, int *fset, size_t fsize,
                    CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_SP, int *, int *, size_t, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_SP, int *, int *, size_t, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_analyze_p");
	return fun(A, Perm, fset, fsize, Common);
}

int attribute_hidden
M_cholmod_band_inplace(int k1, int k2, int mode, CHM_SP A, CHM_CM Common)
{
	static int(*fun)(int, int, int, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, int, int, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_band_inplace");
	return fun(k1, k2, mode, A, Common);
}

int attribute_hidden
M_cholmod_change_factor(int to_xtype, int to_ll, int to_super, int to_packed,
                        int to_monotonic, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(int, int, int, int, int, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
	fun = (int(*)(int, int, int, int, int, CHM_FR, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_change_factor");
	return fun(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
}

CHM_SP attribute_hidden
M_cholmod_copy(CHM_SP A, int stype, int mode, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int, int, CHM_CM) = NULL;
	if (fun == NULL)
	fun = (CHM_SP(*)(CHM_SP, int, int, CHM_CM))
		R_GetCCallable("Matrix", "cholmod_copy");
	return fun(A, stype, mode, Common);
}

CHM_DN attribute_hidden
M_cholmod_copy_dense(CHM_DN  A, CHM_CM Common)
{
	static CHM_DN(*fun)(CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_dense");
	return fun(A, Common);
}

CHM_FR attribute_hidden
M_cholmod_copy_factor(CHM_FR L, CHM_CM Common)
{
	static CHM_FR(*fun)(CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_factor");
	return fun(L, Common);
}

CHM_SP attribute_hidden
M_cholmod_copy_sparse(CHM_SP A, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_copy_sparse");
	return fun(A, Common);
}

int attribute_hidden
M_cholmod_defaults(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_defaults");
	return fun(Common);
}

CHM_SP attribute_hidden
M_cholmod_dense_to_sparse(CHM_DN X, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_DN, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_DN, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_dense_to_sparse");
	return fun(X, values, Common);
}

CHM_SP attribute_hidden
M_cholmod_factor_to_sparse(CHM_FR L, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factor_to_sparse");
	return fun(L, Common);
}

int attribute_hidden
M_cholmod_factorize(CHM_SP A, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factorize");
	return fun(A, L, Common);
}

int attribute_hidden
M_cholmod_factorize_p(CHM_SP A, double beta[2], int *fset,
                      size_t fsize, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(CHM_SP, double[2], int *, size_t, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, double[2], int *, size_t, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_factorize_p");
	return fun(A, beta, fset, fsize, L, Common);
}

int attribute_hidden
M_cholmod_finish(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_finish");
	return fun(Common);
}

int attribute_hidden
M_cholmod_free_dense(CHM_DN *A, CHM_CM Common)
{
	static int(*fun)(CHM_DN *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_DN *,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_dense");
	return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_factor(CHM_FR *L, CHM_CM Common)
{
	static int(*fun)(CHM_FR *,CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_FR *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_factor");
	return fun(L, Common);
}

int attribute_hidden
M_cholmod_free_sparse(CHM_SP *A, CHM_CM Common)
{
	static int(*fun)(CHM_SP *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_sparse");
	return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_triplet(CHM_TR *T, CHM_CM Common)
{
	static int(*fun)(CHM_TR *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_TR *,CHM_CM))
			R_GetCCallable("Matrix", "cholmod_free_triplet");
	return fun(T, Common);
}

int attribute_hidden
M_cholmod_nnz(CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_nnz");
	return fun(A, Common);
}

int attribute_hidden
M_cholmod_scale(CHM_DN S, int scale, CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_DN, int, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_DN, int, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_scale");
	return fun(S, scale, A, Common);
}

#if 0 /* give PRIMME, robustlmm more time to adjust their code */
int attribute_hidden
M_cholmod_sdmult(CHM_SP A, int transpose,
                 double alpha[2], double beta[2],
                 CHM_DN X, CHM_DN Y, CHM_CM Common)
{
	static int(*fun)(CHM_SP, int, double[2], double[2],
	                 CHM_DN, CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, int, double[2], double[2],
		              CHM_DN, CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sdmult");
	return fun(A, transpose, alpha, beta, X, Y, Common);
}
#else
int attribute_hidden
M_cholmod_sdmult(const cholmod_sparse *A, int transpose,
                 const double *alpha, const double *beta,
                 const cholmod_dense *X, cholmod_dense *Y, cholmod_common *Common)
{
	static int(*fun)(const cholmod_sparse *, int, const double *, const double *,
	                 const cholmod_dense *, cholmod_dense *, cholmod_common *) = NULL;
	if (fun == NULL)
		fun = (int(*)(const cholmod_sparse *, int, const double *, const double *,
		              const cholmod_dense *, cholmod_dense *, cholmod_common *))
			R_GetCCallable("Matrix", "cholmod_sdmult");
	return fun(A, transpose, alpha, beta, X, Y, Common);
}
#endif

CHM_DN attribute_hidden
M_cholmod_solve(int sys, CHM_FR L, CHM_DN B, CHM_CM Common)
{
	static CHM_DN(*fun)(int, CHM_FR, CHM_DN, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(int, CHM_FR, CHM_DN, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_solve");
	return fun(sys, L, B, Common);
}

int attribute_hidden
M_cholmod_solve2(int sys, CHM_FR L, CHM_DN B,
                 CHM_DN *X_Handle, CHM_DN *Y_Handle, CHM_DN *E_Handle,
                 CHM_CM Common)
{
	static int(*fun)(int, CHM_FR, CHM_DN, CHM_SP,
	                 CHM_DN *, CHM_SP *, CHM_DN *, CHM_DN *, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, CHM_FR, CHM_DN, CHM_SP,
		              CHM_DN *, CHM_SP *, CHM_DN *, CHM_DN *, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_solve2");
	return fun(sys, L, B, NULL, X_Handle, NULL, Y_Handle, E_Handle, Common);
}

int attribute_hidden
M_cholmod_sort(CHM_SP A, CHM_CM Common)
{
	static int(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sort");
	return fun(A, Common);
}

CHM_DN attribute_hidden
M_cholmod_sparse_to_dense(CHM_SP A, CHM_CM Common)
{
	static CHM_DN(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sparse_to_dense");
	return fun(A, Common);
}

CHM_TR attribute_hidden
M_cholmod_sparse_to_triplet(CHM_SP A, CHM_CM Common)
{
	static CHM_TR(*fun)(CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_TR(*)(CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_sparse_to_triplet");
	return fun(A, Common);
}

CHM_SP attribute_hidden
M_cholmod_speye(size_t nrow, size_t ncol, int xtype, CHM_CM Common)
{
	static CHM_SP(*fun)(size_t, size_t, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(size_t, size_t, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_speye");
	return fun(nrow, ncol, xtype, Common);
}

CHM_SP attribute_hidden
M_cholmod_spsolve(int sys, CHM_FR L, CHM_SP B, CHM_CM Common)
{
	static CHM_SP(*fun)(int, CHM_FR, CHM_SP, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(int,CHM_FR, CHM_SP, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_spsolve");
	return fun(sys, L, B, Common);
}

CHM_SP attribute_hidden
M_cholmod_ssmult(CHM_SP A, CHM_SP B,
                 int stype, int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, int, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, int, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_ssmult");
	return fun(A, B, stype, values, sorted, Common);
}

CHM_SP attribute_hidden
M_cholmod_submatrix(CHM_SP A, int *rset, int rsize, int *cset,
                    int csize, int values, int sorted, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int *, int, int *,
	                    int, int, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int *, int, int *,
		                 int, int, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_submatrix");
	return fun(A, rset, rsize, cset, csize, values, sorted, Common);
}

CHM_SP attribute_hidden
M_cholmod_transpose(CHM_SP A, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_transpose");
	return fun(A, values, Common);
}

CHM_SP attribute_hidden
M_cholmod_triplet_to_sparse(CHM_TR T, int nzmax, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_TR, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_TR, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_triplet_to_sparse");
	return fun(T, nzmax, Common);
}

int attribute_hidden
M_cholmod_updown(int update, CHM_SP C, CHM_FR L, CHM_CM Common)
{
	static int(*fun)(int, CHM_SP, CHM_FR, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(int, CHM_SP, CHM_FR, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_updown");
	return fun(update, C, L, Common);
}

CHM_SP attribute_hidden
M_cholmod_vertcat(CHM_SP A, CHM_SP B, int values, CHM_CM Common)
{
	static CHM_SP(*fun)(CHM_SP, CHM_SP, int, CHM_CM) = NULL;
	if (fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, CHM_SP, int, CHM_CM))
			R_GetCCallable("Matrix", "cholmod_vertcat");
	return fun(A, B, values, Common);
}


/* ---- cholmod_start ----------------------------------------------- */
/* NB: keep synchronized with analogues in ../../src/chm_common.c     */

#if 0
static int attribute_hidden
M_R_cholmod_printf(const char *fmt, ...)
{
	va_list(ap);
	va_start(ap, fmt);
	Rprintf((char *) fmt, ap);
	va_end(ap);
	return 0;
}
#endif

void attribute_hidden
M_R_cholmod_error(int status, const char *file, int line,
                  const char *message)
{
	/* NB: Matrix itself uses CHM_set_common_env, CHM_store_common, and
	       CHM_restore_common to preserve settings through error calls.
	       Consider defining *your* own error handler and restoring the
	       instance of cholmod_common that *you* use.
	*/

	if (status < 0)
		error(    "Cholmod error '%s' at file '%s', line %d", message, file, line);
	else
		warning("Cholmod warning '%s' at file '%s', line %d", message, file, line);
}

int attribute_hidden
M_R_cholmod_start(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_start");
	int ans = fun(Common);
#if 0
	/* No longer, with SuiteSparse 5.7.1 : */
	Common->print_function =
# if 0
		M_R_cholmod_printf;
# else
		NULL;
# endif
#endif
	Common->error_handler = M_R_cholmod_error;
	return ans;
}

int attribute_hidden
M_R_cholmod_finish(CHM_CM Common)
{
	static int(*fun)(CHM_CM) = NULL;
	if (fun == NULL)
		fun = (int(*)(CHM_CM))
			R_GetCCallable("Matrix", "cholmod_finish");
	return fun(Common);
}


/* ==== Matrix stubs ================================================ */

bool Matrix_isclass_ge_dense(SEXP x) {
	static const char *valid[] = { MATRIX_VALID_ge_dense, ""};
	return R_check_class_etc(x, valid) >= 0;
}

bool Matrix_isclass_Csparse(SEXP x) {
	static const char *valid[] = { MATRIX_VALID_Csparse , ""};
	return R_check_class_etc(x, valid) >= 0;
}

CHM_FR attribute_hidden
M_sexp_as_cholmod_factor(CHM_FR L, SEXP from)
{
	static CHM_FR(*fun)(CHM_FR, SEXP) = NULL;
	if(fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, SEXP))
			R_GetCCallable("Matrix", "sexp_as_cholmod_factor");
	return fun(L, from);
}

CHM_SP attribute_hidden
M_sexp_as_cholmod_sparse(CHM_SP A, SEXP from,
                         Rboolean checkUnit, Rboolean sortInPlace)
{
	static CHM_SP(*fun)(CHM_SP, SEXP, Rboolean, Rboolean)= NULL;
	if(fun == NULL)
		fun = (CHM_SP(*)(CHM_SP, SEXP, Rboolean, Rboolean))
			R_GetCCallable("Matrix", "sexp_as_cholmod_sparse");
	return fun(A, from, checkUnit, sortInPlace);
}

CHM_DN attribute_hidden
M_sexp_as_cholmod_dense(CHM_DN A, SEXP from)
{
	static CHM_DN(*fun)(CHM_DN, SEXP) = NULL;
	if(fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, SEXP))
			R_GetCCallable("Matrix", "sexp_as_cholmod_dense");
	return fun(A, from);
}

CHM_DN attribute_hidden
M_numeric_as_cholmod_dense(CHM_DN A, double *data, int nrow, int ncol)
{
	static CHM_DN(*fun)(CHM_DN, double *, int, int) = NULL;
	if (fun == NULL)
		fun = (CHM_DN(*)(CHM_DN, double *, int, int))
			R_GetCCallable("Matrix", "numeric_as_cholmod_dense");
	return fun(A, data, nrow, ncol);
}

SEXP attribute_hidden
M_cholmod_factor_as_sexp(CHM_FR L, int doFree)
{
	static SEXP(*fun)(CHM_FR, int) = NULL;
	if(fun == NULL)
		fun = (SEXP(*)(CHM_FR, int))
			R_GetCCallable("Matrix", "cholmod_factor_as_sexp");
	return fun(L, doFree);
}

SEXP attribute_hidden
M_cholmod_sparse_as_sexp(CHM_SP A, int doFree,
                         int ttype, int doLogic, const char *diagString,
                         SEXP dimnames)
{
	static SEXP(*fun)(CHM_SP, int, int, int, const char *, SEXP) = NULL;
	if (fun == NULL)
		fun = (SEXP(*)(CHM_SP, int, int, int, const char *, SEXP))
			R_GetCCallable("Matrix", "cholmod_sparse_as_sexp");
	return fun(A, doFree, ttype, doLogic, diagString, dimnames);
}

double attribute_hidden
M_cholmod_factor_ldetA(CHM_FR L)
{
	static double(*fun)(CHM_FR) = NULL;
	if (fun == NULL)
		fun = (double(*)(CHM_FR))
			R_GetCCallable("Matrix", "cholmod_factor_ldetA");
	return fun(L);
}

CHM_FR attribute_hidden
M_cholmod_factor_update(CHM_FR L, CHM_SP A, double beta)
{
	static CHM_FR(*fun)(CHM_FR, CHM_SP, double) = NULL;
	if (fun == NULL)
		fun = (CHM_FR(*)(CHM_FR, CHM_SP, double))
			R_GetCCallable("Matrix", "cholmod_factor_update");
	return fun(L, A, beta);
}

#ifdef	__cplusplus
}
#endif
