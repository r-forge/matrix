#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "cholmod.h"

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif

CHM_DN attribute_hidden
M_as_cholmod_dense(CHM_DN ans, SEXP x)
{
    static CHM_DN(*fun)(CHM_DN,SEXP) = NULL;
    if(fun == NULL)
	fun = (CHM_DN(*)(CHM_DN,SEXP))
	    R_GetCCallable("Matrix", "as_cholmod_dense");
    return fun(ans, x);
}

CHM_FR attribute_hidden
M_as_cholmod_factor(CHM_FR ans, SEXP x)
{
    static CHM_FR(*fun)(CHM_FR,SEXP) = NULL;
    if(fun == NULL)
	fun = (CHM_FR(*)(CHM_FR,SEXP))
	    R_GetCCallable("Matrix", "as_cholmod_factor");
    return fun(ans, x);
}

CHM_SP attribute_hidden
M_as_cholmod_sparse(CHM_SP ans, SEXP x, Rboolean check_Udiag, Rboolean sort_in_place)
{
    static CHM_SP(*fun)(CHM_SP,SEXP,Rboolean,Rboolean)= NULL;
    if(fun == NULL)
	fun = (CHM_SP(*)(CHM_SP,SEXP,Rboolean,Rboolean))
	    R_GetCCallable("Matrix", "as_cholmod_sparse");
    return fun(ans, x, check_Udiag, sort_in_place);
}

CHM_SP attribute_hidden
M_as_cholmod_triplet(CHM_SP ans, SEXP x, Rboolean check_Udiag)
{
    static CHM_SP(*fun)(CHM_SP,SEXP,Rboolean)= NULL;
    if(fun == NULL)
	fun = (CHM_SP(*)(CHM_SP,SEXP,Rboolean))
	    R_GetCCallable("Matrix", "as_cholmod_triplet");
    return fun(ans, x, check_Udiag);
}

SEXP attribute_hidden
M_Csparse_diagU2N(SEXP x)
{
    static SEXP(*fun)(SEXP) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(SEXP))
	    R_GetCCallable("Matrix", "Csparse_diagU2N");
    return fun(x);
}

SEXP attribute_hidden
M_chm_factor_to_SEXP(const cholmod_factor *f, int dofree)
{
    static SEXP(*fun)(const cholmod_factor*,int) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(const cholmod_factor*,int))
	    R_GetCCallable("Matrix", "chm_factor_to_SEXP");
    return fun(f, dofree);
}

double attribute_hidden
M_chm_factor_ldetL2(const cholmod_factor *f)
{
    static double(*fun)(const cholmod_factor*) = NULL;
    if(fun == NULL)
	fun = (double(*)(const cholmod_factor*))
	    R_GetCCallable("Matrix", "chm_factor_ldetL2");
    return fun(f);
}

CHM_FR attribute_hidden
M_chm_factor_update(CHM_FR f, const cholmod_sparse *A, double mult)
{
    static CHM_FR(*fun)(CHM_FR,const cholmod_sparse*,double) = NULL;
    if(fun == NULL)
	fun = (CHM_FR(*)(CHM_FR,const cholmod_sparse*,double))
	    R_GetCCallable("Matrix", "chm_factor_update");
    return fun(f, A, mult);
}

SEXP attribute_hidden
M_chm_sparse_to_SEXP(const cholmod_sparse *a, int dofree,
		     int uploT, int Rkind, char *diag, SEXP dn)
{
    static SEXP(*fun)(const cholmod_sparse*,int,int,int,char*,SEXP) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(const cholmod_sparse*,int,int,int,char*,SEXP))
	    R_GetCCallable("Matrix", "chm_sparse_to_SEXP");
    return fun(a, dofree, uploT, Rkind, diag, dn);
}

SEXP attribute_hidden
M_chm_triplet_to_SEXP(const CHM_TR a, int dofree,
		      int uploT, int Rkind, const char *diag, SEXP dn)
{
    static SEXP(*fun)(const CHM_TR,int,int,int,const char*,SEXP) = NULL;
    if(fun == NULL)
	fun = (SEXP(*)(const CHM_TR,int,int,int,const char*,SEXP))
	    R_GetCCallable("Matrix", "chm_triplet_to_SEXP");
    return fun(a, dofree, uploT, Rkind, diag, dn);
}

CHM_SP attribute_hidden
M_cholmod_aat(const cholmod_sparse *A, int *fset, size_t fsize,
	      int mode, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,int*,size_t,
			int,CHM_CM) = NULL;
    if(fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,int*,size_t,
			 int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_aat");
    return fun(A, fset, fsize, mode, Common);
}

int attribute_hidden
M_cholmod_band_inplace(CHM_SP A, int k1, int k2, int mode,
		       CHM_CM Common)
{
    static int(*fun)(CHM_SP,int,int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_SP,int,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_band_inplace");
    return fun(A, k1, k2, mode, Common);
}

CHM_SP attribute_hidden
M_cholmod_add(const cholmod_sparse *A, const cholmod_sparse *B,
	      double alpha[2], double beta[2], int values,
	      int sorted, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,const cholmod_sparse*,
			double*,double*,int,int,
			CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,const cholmod_sparse*,
			 double*,double*,int,int,
			 CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_add");
    return fun(A, B, alpha, beta, values, sorted, Common);
}

CHM_DN attribute_hidden
M_cholmod_allocate_dense(size_t nrow, size_t ncol, size_t d,
			 int xtype, CHM_CM Common)
{
    static CHM_DN(*fun)(size_t,size_t,size_t,
			int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_DN(*)(size_t,size_t,size_t,
			 int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_allocate_dense");
    return fun(nrow, ncol, d, xtype, Common);
}

CHM_SP attribute_hidden
M_cholmod_allocate_sparse(size_t nrow, size_t ncol, size_t nzmax,
			  int sorted, int packed, int stype,
			  int xtype, CHM_CM Common)
{
    static CHM_SP(*fun)(size_t,size_t,size_t,int,int,
			int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)
	       (size_t,size_t,size_t,int,int,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_allocate_sparse");
    return fun(nrow,ncol,nzmax,sorted,packed,stype,xtype,Common);
}

CHM_TR attribute_hidden
M_cholmod_allocate_triplet(size_t nrow, size_t ncol, size_t nzmax,
			   int stype, int xtype, CHM_CM Common)
{
    static CHM_TR(*fun)(size_t,size_t,size_t, int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_TR(*)(size_t,size_t,size_t,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_allocate_triplet");
    return fun(nrow,ncol,nzmax,stype,xtype,Common);
}

CHM_SP attribute_hidden
M_cholmod_triplet_to_sparse(const cholmod_triplet* T, int nzmax,
			    CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_triplet*,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_triplet*,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_triplet_to_sparse");
    return fun(T, nzmax, Common);
}

CHM_TR attribute_hidden
M_cholmod_sparse_to_triplet(const cholmod_sparse *A, CHM_CM Common)
{
    static CHM_TR(*fun)(const cholmod_sparse*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_TR(*)(const cholmod_sparse*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_sparse_to_triplet");
    return fun(A, Common);
}

CHM_DN attribute_hidden
M_cholmod_sparse_to_dense(const cholmod_sparse *A, CHM_CM Common)
{
    static CHM_DN(*fun)(const cholmod_sparse*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_DN(*)(const cholmod_sparse*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_sparse_to_dense");
    return fun(A, Common);
}

CHM_FR attribute_hidden
M_cholmod_analyze(const cholmod_sparse *A, CHM_CM Common)
{
    static CHM_FR(*fun)(const cholmod_sparse*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_FR(*)(const cholmod_sparse*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_analyze");
    return fun(A, Common);
}

CHM_FR attribute_hidden
M_cholmod_analyze_p(const cholmod_sparse *A, int *Perm, int *fset,
		    size_t fsize, CHM_CM Common)
{
    static CHM_FR(*fun)(const cholmod_sparse*,int*,int*,size_t,
			CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_FR(*)(const cholmod_sparse*,int*,int*,
			 size_t,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_analyze_p");
    return fun(A, Perm, fset, fsize, Common);
}

CHM_SP attribute_hidden
M_cholmod_copy(const cholmod_sparse *A, int stype,
	       int mode, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_copy");
    return fun(A, stype, mode, Common);
}

CHM_DN attribute_hidden
M_cholmod_copy_dense(const cholmod_dense * A, CHM_CM Common)
{
    static CHM_DN(*fun)(const cholmod_dense*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_DN(*)(const cholmod_dense*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_copy_dense");
    return fun(A, Common);
}

CHM_FR attribute_hidden
M_cholmod_copy_factor(const cholmod_factor *L, CHM_CM Common)
{
    static CHM_FR(*fun)(const cholmod_factor*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_FR(*)(const cholmod_factor*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_copy_factor");
    return fun(L, Common);
}

int attribute_hidden
M_cholmod_change_factor(int to_xtype, int to_ll, int to_super, int to_packed,
			int to_monotonic, CHM_FR L, CHM_CM Common)
{
    static int(*fun)(int,int,int,int,int,CHM_FR,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(int,int,int,int,int,CHM_FR,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_change_factor");
    return fun(to_xtype, to_ll, to_super, to_packed, to_monotonic, L, Common);
}

CHM_SP attribute_hidden
M_cholmod_copy_sparse(const cholmod_sparse *A, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_copy_sparse");
    return fun(A, Common);
}

CHM_SP attribute_hidden
M_cholmod_factor_to_sparse(const cholmod_factor *L, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_factor*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_factor*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_factor_to_sparse");
    return fun(L, Common);
}

CHM_SP attribute_hidden
M_cholmod_submatrix(const cholmod_sparse *A, int *rset, int rsize, int *cset,
		    int csize, int values, int sorted, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,int*,int,int*,int,
			int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,int*,int,int*,
			 int,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_submatrix");
    return fun(A, rset, rsize, cset, csize, values, sorted, Common);
}

CHM_SP attribute_hidden
M_cholmod_dense_to_sparse(const cholmod_dense * X, int values, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_dense*,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_dense*,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_dense_to_sparse");
    return fun(X, values, Common);
}

int attribute_hidden
M_cholmod_factorize(const cholmod_sparse *A, CHM_FR L, CHM_CM Common)
{
    static int(*fun)(const cholmod_sparse*,CHM_FR,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(const cholmod_sparse*,CHM_FR,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_factorize");
    return fun(A, L, Common);
}

int attribute_hidden
M_cholmod_factorize_p(const cholmod_sparse *A, double *beta, int *fset,
		      size_t fsize, CHM_FR L,
		      CHM_CM Common)
{
    static int(*fun)(const cholmod_sparse*,double*,int*,size_t,
		     CHM_FR,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(const cholmod_sparse*,double*,int*,size_t,
		      CHM_FR,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_factorize_p");
    return fun(A, beta, fset, fsize, L, Common);
}

int attribute_hidden
M_cholmod_finish(CHM_CM Common)
{

    static int(*fun)(CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_finish");
    return fun(Common);
}

int attribute_hidden
M_cholmod_sort(CHM_SP A, CHM_CM Common)
{
    static int(*fun)(CHM_SP,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_SP,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_sort");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_dense(CHM_DN  *A, CHM_CM Common)
{
    static int(*fun)(CHM_DN*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_DN*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_free_dense");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_factor(CHM_FR *L, CHM_CM Common)
{
    static int(*fun)(CHM_FR*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_FR*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_free_factor");
    return fun(L, Common);
}

int attribute_hidden
M_cholmod_free_sparse(CHM_SP *A, CHM_CM Common)
{
    static int(*fun)(CHM_SP*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_SP*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_free_sparse");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_free_triplet(cholmod_triplet **T, CHM_CM Common)
{
    static int(*fun)(cholmod_triplet**,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(cholmod_triplet**,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_free_triplet");
    return fun(T, Common);
}

long attribute_hidden
M_cholmod_nnz(const cholmod_sparse *A, CHM_CM Common)
{
    static long(*fun)(const cholmod_sparse*,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (long(*)(const cholmod_sparse*,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_nnz");
    return fun(A, Common);
}

int attribute_hidden
M_cholmod_sdmult(const cholmod_sparse *A, int transpose,
		 const double *alpha, const double *beta,
		 const cholmod_dense *X, CHM_DN  Y,
		 CHM_CM Common)
{
    static int(*fun)(const cholmod_sparse*,int,const double*,
		     const double*,const cholmod_dense*,
		     CHM_DN,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(const cholmod_sparse*,int,const double*,
		      const double*, const cholmod_dense*,
		      CHM_DN,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_sdmult");
    return fun(A, transpose, alpha, beta, X, Y, Common);
}

CHM_SP attribute_hidden
M_cholmod_ssmult(const cholmod_sparse *A, const cholmod_sparse *B,
		 int stype, int values, int sorted,
		 CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,const cholmod_sparse*,
			int,int,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,const cholmod_sparse*,
			 int,int,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_ssmult");
    return fun(A, B, stype, values, sorted, Common);
}

CHM_DN attribute_hidden
M_cholmod_solve(int sys, const cholmod_factor *L,
		const cholmod_dense * B, CHM_CM Common)
{
    static CHM_DN(*fun)(int,const cholmod_factor*,const cholmod_dense*,
			CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_DN(*)(int,const cholmod_factor*,const cholmod_dense*,
			 CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_solve");
    return fun(sys, L, B, Common);
}

CHM_SP attribute_hidden
M_cholmod_speye(size_t nrow, size_t ncol,
		int xtype, CHM_CM Common)
{
    static CHM_SP(*fun)(size_t,size_t,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(size_t,size_t,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_speye");
    return fun(nrow, ncol, xtype, Common);
}

CHM_SP attribute_hidden
M_cholmod_spsolve(int sys, const cholmod_factor *L,
		  const cholmod_sparse *B, CHM_CM Common)
{
    static CHM_SP(*fun)(int,const cholmod_factor*,
			const cholmod_sparse*, CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(int,const cholmod_factor*,
			 const cholmod_sparse*, CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_spsolve");
    return fun(sys, L, B, Common);
}

int attribute_hidden
M_cholmod_defaults (CHM_CM Common)
{
    static int(*fun)(CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_defaults");
    return fun(Common);
}

/* extern cholmod_common c; */

void attribute_hidden
M_R_cholmod_error(int status, const char *file, int line, const char *message)
{
/* NB: keep in sync with R_cholmod_error(), ../../src/chm_common.c */

    if(status < 0) {
/* Note: Matrix itself uses CHM_set_common_env, CHM_store_common 
 *   and CHM_restore_common to preserve settings through error calls.
 *  Consider defining your own error handler, *and* possibly restoring
 *  *your* version of the cholmod_common that *you* use.
 */
	error("Cholmod error '%s' at file:%s, line %d", message, file, line);
    }
    else
	warning("Cholmod warning '%s' at file:%s, line %d",
		message, file, line);
}

#if 0  /* no longer used */
/* just to get 'int' instead of 'void' as required by CHOLMOD's print_function */
static int
R_cholmod_printf(const char* fmt, ...)
{
    va_list(ap);

    va_start(ap, fmt);
    Rprintf((char *)fmt, ap);
    va_end(ap);
    return 0;
}
#endif

int attribute_hidden
M_R_cholmod_start(CHM_CM Common)
{
    int val;
    static int(*fun)(CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_start");
    val = fun(Common);
/*-- NB: keep in sync with  R_cholmod_l_start() --> ../../src/chm_common.c */
    /* do not allow CHOLMOD printing - currently */
    Common->print_function = NULL;/* was  R_cholmod_printf; /.* Rprintf gives warning */
/* Consider using your own error handler: */
    Common->error_handler = M_R_cholmod_error;
    return val;
}

CHM_SP attribute_hidden
M_cholmod_transpose(const cholmod_sparse *A, int values, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,int,CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_transpose");
    return fun(A, values, Common);
}

CHM_SP attribute_hidden
M_cholmod_vertcat(const cholmod_sparse *A, const cholmod_sparse *B, int values, CHM_CM Common)
{
    static CHM_SP(*fun)(const cholmod_sparse*,const cholmod_sparse*,int,CHM_CM) = NULL;
    if (fun == NULL)
	fun = (CHM_SP(*)(const cholmod_sparse*,const cholmod_sparse*, int, CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_vertcat");
    return fun(A, B, values, Common);
}

SEXP attribute_hidden
M_dpoMatrix_chol(SEXP x)
{
    static SEXP(*fun)(SEXP) = NULL;
    if (fun == NULL)
	fun = (SEXP(*)(SEXP))
	    R_GetCCallable("Matrix", "dpoMatrix_chol");
    return fun(x);
}

CHM_DN attribute_hidden
M_numeric_as_chm_dense(CHM_DN ans, double *v, int nr, int nc)
{
    static CHM_DN(*fun)(CHM_DN,double*,int,int) = NULL;
    if (fun == NULL)
	fun = (CHM_DN(*)(CHM_DN,double*,int,int))
	    R_GetCCallable("Matrix", "numeric_as_chm_dense");
    return fun(ans, v, nr, nc);
}

int attribute_hidden
M_cholmod_scale(const cholmod_dense *S, int scale, CHM_SP A,
		CHM_CM Common)
{
    static int(*fun)(const cholmod_dense*,int,CHM_SP, CHM_CM) = NULL;
    if (fun == NULL)
	fun = (int(*)(const cholmod_dense*,int,CHM_SP, CHM_CM))
	    R_GetCCallable("Matrix", "cholmod_l_scale");
    return fun(S, scale, A, Common);
}

#ifdef	__cplusplus
}
#endif
