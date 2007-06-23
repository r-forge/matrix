/*------ Definition of a template for [diln]gCMatrix_colsums(...) : *
 *                       --------     ~~~~~~~~~~~~~~~~~~~~~~
 * i.e., included several times from ./dgCMatrix.c
 *                                   ~~~~~~~~~~~~~
 */


/* for all cases with an 'x' slot -- i.e. almost all cases ;
 * just redefine this in the other cases:
 */

#ifdef _dgC_

# define gCMatrix_colSums dgCMatrix_colSums
# define _DOUBLE_ans
# define _has_x_slot_
# define Type_x_  double
# define NA_x_    NA_REAL
# undef _dgC_

#elif defined (_igC_)

# define gCMatrix_colSums igCMatrix_colSums
# define _DOUBLE_ans
# define _has_x_slot_
# define Type_x_  int
# define NA_x_    NA_INTEGER
# undef _igC_

#elif defined (_lgC_)

# define gCMatrix_colSums lgCMatrix_colSums_i
# define _INT_ans
# define _has_x_slot_
# define Type_x_  int
# define NA_x_    NA_LOGICAL
# undef _lgC_

#elif defined (_lgC_mn)

# define gCMatrix_colSums lgCMatrix_colSums_d
# define _DOUBLE_ans
# define _has_x_slot_
# define Type_x_  int
# define NA_x_    NA_LOGICAL
# undef _lgC_mn

#elif defined (_ngC_)

# define gCMatrix_colSums ngCMatrix_colSums_i
# define _INT_ans
# undef  has_x_slot_ /* withOUT 'x' slot */
# undef _ngC_

#elif defined (_ngC_mn)

# define gCMatrix_colSums ngCMatrix_colSums_d
# define _DOUBLE_ans
# undef  has_x_slot_ /* withOUT 'x' slot */
# undef _ngC_mn

#elif defined (_zgC_)

#  error "zgC* not yet implemented"

#else

#  error "no valid  _[dilnz]gC_ option"

#endif

/* - - - - - - - - - - - - - - - - - - - - */

#ifdef _DOUBLE_ans

# define SparseResult_class "dsparseVector"
# define Type_ans double
# define STYP_ans REAL
# define NA_ans NA_REAL
# define SXP_ans  REALSXP
#undef _DOUBLE_ans

#elif defined (_INT_ans)

# define SparseResult_class "isparseVector"
# define Type_ans int
# define STYP_ans INTEGER
# define NA_ans NA_INTEGER
# define SXP_ans  INTSXP
#undef _INT_ans

#else
#  error "invalid macro logic"
#endif

/* - - - - - - - - - - - - - - - - - - - - */

#ifdef _has_x_slot_
# define ColSUM_column(_i1_,_i2_,_SUM_)					\
		if(mn) dnm = cx->nrow;	/* denominator for means */	\
		for(i = _i1_, _SUM_ = 0; i < _i2_; i++)			\
		    if (mn && na_rm && xx[i] == NA_x_)			\
			dnm--; /* skip NAs but decrement denominator*/	\
		    else _SUM_ += xx[i];				\
		if(mn) _SUM_ = (dnm > 0) ? _SUM_/dnm : NA_ans

#else /* no 'x' slot */
# define ColSUM_column(_i1_,_i2_,_SUM_)		\
		_SUM_ = _i2_ - _i1_;		\
		if(mn) _SUM_ /= cx->nrow
#endif

/* Now comes the template -- which depends on the above macros : */

SEXP gCMatrix_colSums(SEXP x, SEXP NArm, SEXP spRes, SEXP trans, SEXP means)
{
    int mn = asLogical(means), sp = asLogical(spRes), tr = asLogical(trans);
    cholmod_sparse *cx = as_cholmod_sparse(x);
#ifdef _has_x_slot_
    int na_rm = asLogical(NArm), i, dnm = 0/*Wall*/;
    Type_x_ *xx = (Type_x_ *)(cx -> x);
#endif
    int j, n = (tr ? cx->nrow : cx->ncol);
    int *xp = (int *)(cx->p);
    SEXP ans = PROTECT(sp ? NEW_OBJECT(MAKE_CLASS(SparseResult_class))
			  : allocVector(SXP_ans, n));

    if (tr) {
	cholmod_sparse *cxt = cholmod_transpose(cx, 1 /*values*/, &c);
	Free(cx);
	cx = cxt;
    }

    if (sp) { /* sparseResult - never allocating length-n ... */
	int nza, i1, i2, p, *ai;
	Type_ans *ax;

	for (j = 0, nza = 0; j < cx->ncol; j++)
	    if(xp[j] < xp[j + 1])
		nza++;

	ai =  INTEGER(ALLOC_SLOT(ans, Matrix_iSym, INTSXP,  nza));
	ax = STYP_ans(ALLOC_SLOT(ans, Matrix_xSym, SXP_ans, nza));

	SET_SLOT(ans, Matrix_lengthSym, ScalarInteger(n));

	i2 = xp[0];
	for (j = 1, p = 0; j <= cx->ncol; j++) {
	    /* j' =j+1, since 'i' slot will be 1-based */
	    i1 = i2; i2 = xp[j];
	    if(i1 < i2) {
		Type_ans sum;
		ColSUM_column(i1,i2, sum);

		ai[p]	= j;
		ax[p++] = sum;
	    }
	}
    }
    else { /* "numeric" (non sparse) result */
	Type_ans *a = STYP_ans(ans);
	for (j = 0; j < cx->ncol; j++) {
	    ColSUM_column(xp[j], xp[j + 1], a[j]);
	}
    }

    if (tr) cholmod_free_sparse(&cx, &c); else Free(cx);
    UNPROTECT(1);
    return ans;
}

#undef ColSUM_column

#undef NA_ans
#undef NA_x_
#undef STYP_ans
#undef SXP_ans
#undef SparseResult_class
#undef Type_ans
#undef Type_x_
#undef _has_x_slot_

#undef gCMatrix_colSums
