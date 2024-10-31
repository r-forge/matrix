/* C implementation of methods for t, ct */

#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

SEXP dense_transpose(SEXP from, const char *class, char op_ct)
{
	int packed = class[2] == 'p';
	char ul = '\0';

	if (class[0] != 'z')
		op_ct = '\0';

	SEXP to = PROTECT(newObject(class));

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, n, m);
	SET_DIMNAMES(to, class[1] != 's' && class[1] != 'p' && class[1] != 'o', DIMNAMES(from, 0));
	if (class[1] != 'g' && (ul = UPLO(from)) == 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);
	if (class[1] == 'o')
		COPY_SLOT(to, from, Matrix_sdSym);

	SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)),
		x1 = PROTECT(Rf_allocVector(TYPEOF(x0), XLENGTH(x0)));
	size_t m_ = (size_t) m, n_ = (size_t) n;

#define TEMPLATE(c) \
	do { \
		c##TYPE *px0 = c##PTR(x0), *px1 = c##PTR(x1); \
		if (!packed) \
			c##NAME(trans2)(px1, px0, m_, n_, op_ct); \
		else \
			c##NAME(trans1)(px1, px0, n_, ul, op_ct); \
	} while (0)

	SWITCH4((class[0] == 'c') ? 'd' : class[0], TEMPLATE);

#undef TEMPLATE

	SET_SLOT(to, Matrix_xSym, x1);

	UNPROTECT(3); /* x1, x0, to */
	return to;
}

SEXP sparse_transpose(SEXP from, const char *class, char op_ct,
                      int lazy)
{
	if (class[0] != 'z')
		op_ct = '\0';

	SEXP to;
	if (class[2] == 'T' || !lazy)
		PROTECT(to = newObject(class));
	else {
		char cl[] = "...Matrix";
		cl[0] = class[0];
		cl[1] = class[1];
		cl[2] = (class[2] == 'C') ? 'R' : 'C';
		PROTECT(to = newObject(cl));
	}

	int *pdim = DIM(from), m = pdim[0], n = pdim[1];
	SET_DIM(to, n, m);
	SET_DIMNAMES(to, class[1] != 's' && class[1] != 'p', DIMNAMES(from, 0));
	if (class[1] != 'g' && UPLO(from) == 'U')
		SET_UPLO(to);
	if (class[1] == 's' && class[0] == 'z' && TRANS(from) != 'C')
		SET_TRANS(to);
	if (class[1] == 't' && DIAG(from) != 'N')
		SET_DIAG(to);

	if (class[2] != 'T') {

		if (class[2] == 'R')
			SWAP(m, n, int, );

		SEXP iSym = (class[2] == 'C') ? Matrix_iSym : Matrix_jSym,
			jSym = (class[2] == 'C') ? Matrix_jSym : Matrix_iSym,
			p0 = PROTECT(GET_SLOT(from, Matrix_pSym)),
			i0 = PROTECT(GET_SLOT(from,        iSym));
		int *pp0 = INTEGER(p0), *pi0 = INTEGER(i0), nnz = INTEGER(p0)[n];

		if (lazy) {

		SET_SLOT(to, Matrix_pSym, p0);
		SET_SLOT(to,        jSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (op_ct != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, nnz));
				zvconj(COMPLEX(x1), COMPLEX(x0), (size_t) nnz);
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		} else {

		SEXP p1 = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) m + 1)),
			i1 = PROTECT(Rf_allocVector(INTSXP, nnz));
		int *pp1 = INTEGER(p1), *pi1 = INTEGER(i1), *iwork = NULL;
		SET_SLOT(to, Matrix_pSym, p1);
		SET_SLOT(to,        iSym, i1);
		Matrix_Calloc(iwork, m, int);

#define TEMPLATE(c) \
		do { \
			c##TYPE *px0 = NULL, *px1 = NULL; \
			c##IF_NPATTERN( \
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym)), \
				x1 = PROTECT(Rf_allocVector(c##TYPESXP, nnz)); \
			px0 = c##PTR(x0); \
			px1 = c##PTR(x1); \
			SET_SLOT(to, Matrix_xSym, x1); \
			UNPROTECT(2); /* x1, x0 */ \
			); \
			c##csptrans(pp1, pi1, px1, pp0, pi0, px0, m, n, op_ct, iwork); \
		} while (0)

		SWITCH5(class[0], TEMPLATE);

#undef TEMPLATE

		Matrix_Free(iwork, m);
		UNPROTECT(2); /* j1, p1 */

		}

		UNPROTECT(2); /* i0, p0 */

	} else {

		SEXP i0 = PROTECT(GET_SLOT(from, Matrix_iSym)),
			j0 = PROTECT(GET_SLOT(from, Matrix_jSym));
		SET_SLOT(to, Matrix_iSym, j0);
		SET_SLOT(to, Matrix_jSym, i0);

		if (class[0] != 'n') {
			SEXP x0 = PROTECT(GET_SLOT(from, Matrix_xSym));
			if (op_ct != 'C')
				SET_SLOT(to, Matrix_xSym, x0);
			else {
				SEXP x1 = PROTECT(Rf_allocVector(CPLXSXP, XLENGTH(x0)));
				zvconj(COMPLEX(x1), COMPLEX(x0), (size_t) XLENGTH(x0));
				SET_SLOT(to, Matrix_xSym, x1);
				UNPROTECT(1); /* x1 */
			}
			UNPROTECT(1); /* x0 */
		}

		UNPROTECT(2); /* j0, i0 */

	}

	UNPROTECT(1); /* to */
	return to;
}

SEXP R_dense_transpose(SEXP s_from, SEXP s_trans)
{
	const char *class = Matrix_class(s_from, valid_dense, 0, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	return dense_transpose(s_from, class, ct);
}

SEXP R_sparse_transpose(SEXP s_from, SEXP s_trans, SEXP s_lazy)
{
	const char *class = Matrix_class(s_from, valid_sparse, 0, __func__);

	char ct;
	VALID_TRANS(s_trans, ct);

	int lazy;
	VALID_LOGIC2(s_lazy, lazy);

	return sparse_transpose(s_from, class, ct, lazy);
}
