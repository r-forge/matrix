/* MJ: work in progress ... see also ../R/bind2.R */
#include "bind.h"

static void scanArgs(SEXP args, SEXP exprs, int margin, int level,
                     int *rdim, int *rdimnames, char *kind, char *repr)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	const char *cl;
	SEXP a, e, x, tmp;
	int nS4 = 0, nDense = 0,
		anyCsparse = 0, anyRsparse = 0, anyTsparse = 0, anyDiagonal = 0,
		anyN = 0, anyL = 0, anyI = 0, anyD = 0, anyZ = 0,
		i, ivalid, *xdim;
	R_xlen_t xlen;

	rdim[!margin] = -1;
	rdim[ margin] =  0;
	rdimnames[0] = rdimnames[1] = 0;

	for (a = args; a != R_NilValue; a = CDR(a)) {
		x = CAR(a);
		if (x == R_NilValue)
			continue;
		if (IS_S4_OBJECT(x)) {
			++nS4;
			ivalid = R_check_class_etc(x, valid);
			if (ivalid < 0) {
				if (margin == 1)
					ERROR_INVALID_CLASS(x, "cbind.Matrix");
				else
					ERROR_INVALID_CLASS(x, "rbind.Matrix");
			}
			cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

			tmp = GET_SLOT(x, Matrix_DimSym);
			xdim = INTEGER(tmp);
			if (rdim[!margin] < 0)
				rdim[!margin] = xdim[!margin];
			else if (xdim[!margin] != rdim[!margin]) {
				if (margin == 1)
					error(_("number of rows of matrices must match"));
				else
					error(_("number of columns of matrices must match"));
			}
			if (xdim[margin] > INT_MAX - rdim[margin])
				error(_("dimensions cannot exceed 2^31-1"));
			rdim[margin] += xdim[margin];

			if (!rdimnames[0] || !rdimnames[1]) {
				tmp = GET_SLOT(x, Matrix_DimNamesSym);
				if (cl[1] == 's') {
					if (VECTOR_ELT(tmp, 0) != R_NilValue ||
					    VECTOR_ELT(tmp, 1) != R_NilValue)
						rdimnames[0] = rdimnames[1] = 1;
				} else
					for (i = 0; i < 2; ++i)
						if (!rdimnames[i] &&
						    VECTOR_ELT(tmp, i) != R_NilValue)
							rdimnames[i] = 1;
			}

			switch (cl[0]) {
			case 'n':
				anyN = 1;
				break;
			case 'l':
				anyL = 1;
				break;
			case 'i':
				if (cl[2] != 'd')
				anyI = 1;
				break;
			case 'd':
				anyD = 1;
				break;
			case 'z':
				anyZ = 1;
				break;
			default:
				break;
			}

			switch (cl[2]) {
			case 'e':
			case 'y':
			case 'r':
			case 'p':
				++nDense;
				break;
			case 'C':
				anyCsparse = 1;
				break;
			case 'R':
				anyRsparse = 1;
				break;
			case 'T':
				SETCAR(a, Tsparse_aggregate(x));
				anyTsparse = 1;
				break;
			case 'i':
				anyDiagonal = 1;
				break;
			case 'd':
				PROTECT(tmp = GET_SLOT(x, Matrix_marginSym));
				if (INTEGER(tmp)[0] - 1 != margin) {
					anyN = 1;
					if (margin == 1)
						anyCsparse = 1;
					else
						anyRsparse = 1;
				}
				UNPROTECT(1);
				break;
			default:
				break;
			}
		} else {
			switch (TYPEOF(x)) {
			case LGLSXP:
				anyL = 1;
				break;
			case INTSXP:
				anyI = 1;
				break;
			case REALSXP:
				anyD = 1;
				break;
			case CPLXSXP:
				anyZ = 1;
				break;
			default:
				if (margin == 1)
					ERROR_INVALID_TYPE("object", TYPEOF(x), "cbind.Matrix");
				else
					ERROR_INVALID_TYPE("object", TYPEOF(x), "rbind.Matrix");
				break;
			}

			tmp = getAttrib(x, R_DimSymbol);
			if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2) {
				xdim = INTEGER(tmp);
				if (rdim[!margin] < 0)
					rdim[!margin] = xdim[!margin];
				else if (rdim[!margin] != xdim[!margin]) {
					if (margin == 1)
						error(_("number of rows of matrices must match"));
					else
						error(_("number of columns of matrices must match"));
				}
				if (xdim[margin] > INT_MAX - rdim[margin])
					error(_("dimensions cannot exceed 2^31-1"));
				rdim[margin] += xdim[margin];

				if (!rdimnames[0] || !rdimnames[1]) {
					tmp = getAttrib(x, R_DimNamesSymbol);
					if (tmp != R_NilValue)
						for (i = 0; i < 2; ++i)
							if (!rdimnames[i] &&
							    VECTOR_ELT(tmp, i) != R_NilValue)
								rdimnames[i] = 1;
				}
			}
		}
	}

	for (a = args, e = exprs; a != R_NilValue; a = CDR(a), e = CDR(e)) {
		x = CAR(a);
		if (x == R_NilValue || IS_S4_OBJECT(x))
			continue;
		tmp = getAttrib(x, R_DimSymbol);
		if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2)
			continue;
		xlen = XLENGTH(x);
		if (rdim[!margin] > 0 && xlen == 0)
			continue;
		if (rdim[margin] == INT_MAX)
			error(_("dimensions cannot exceed 2^31-1"));
		rdim[margin] += 1;
		if (xlen > rdim[!margin] || rdim[!margin] % (int) xlen) {
			if (margin == 1)
				warning(_("number of rows of result is not a multiple of vector length"));
			else
				warning(_("number of columns of result is not a multiple of vector length"));
		}
		if (!rdimnames[ margin]) {
			if (TAG(a) != R_NilValue ||
			    level == 2 || (level == 1 && TYPEOF(CAR(e)) == SYMSXP))
				rdimnames[ margin] = 1;
		}
		if (!rdimnames[!margin] && xlen == rdim[!margin]) {
			tmp = getAttrib(x, R_NamesSymbol);
			if (tmp != R_NilValue)
				rdimnames[!margin] = 1;
		}
	}

	if (anyZ)
		*kind = 'z';
	else if (anyD)
		*kind = 'd';
	else if (anyI)
		*kind = 'i';
	else if (anyL)
		*kind = 'l';
	else if (anyN)
		*kind = 'n';
	else
		*kind = '\0';

	if (nDense == nS4)
		*repr = 'e';
	else if (nDense == 0) {
		if (anyCsparse && anyRsparse)
			*repr = (margin == 1) ? 'C' : 'R';
		else if (anyCsparse)
			*repr = 'C';
		else if (anyRsparse)
			*repr = 'R';
		else if (anyTsparse)
			*repr = 'T';
		else if (anyDiagonal)
			*repr = (margin == 1) ? 'C' : 'R';
		else
			*repr = '\0';
	} else {
		/* The length of the result is at most INT_MAX * INT_MAX,
		   which cannot overflow Matrix_int_fast64_t as long as R
		   builds require sizeof(int) equal to 4
		 */
		Matrix_int_fast64_t nnz = 0, len = 0, xnnz, xlen;
		for (a = args; a != R_NilValue && nnz < INT_MAX; a = CDR(a)) {
			x = CAR(a);
			if (!IS_S4_OBJECT(x))
				continue;
			ivalid = R_check_class_etc(x, valid);
			cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];

			PROTECT(tmp = GET_SLOT(x, Matrix_DimSym));
			xdim = INTEGER(tmp);
			xlen = (Matrix_int_fast64_t) xdim[0] * xdim[1];

			switch (cl[2]) {
			case 'e':
			case 'y':
			case 'r':
			case 'p':
				xnnz = (cl[1] != 't') ? xlen : ((xlen + xdim[0]) / 2);
				break;
			case 'C':
			case 'R':
			{
				SEXP p = PROTECT(GET_SLOT(x, Matrix_pSym));
				int *pp = INTEGER(p), n = xdim[(cl[2] == 'C') ? 1 : 0];
				xnnz = pp[n];
				if (cl[1] == 's') {
					SEXP iSym = (cl[2] == 'C') ? Matrix_iSym : Matrix_jSym,
						i = PROTECT(GET_SLOT(x, iSym));
					int *pi = INTEGER(i), j;
					xnnz *= 2;
					if (*CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0)) == 'U') {
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j + 1] - 1] == j)
								--xnnz;
					} else {
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j]] == j)
								--xnnz;
					}
					UNPROTECT(1);
				} else if (cl[1] == 't' && *CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0)) != 'N')
					xnnz += xdim[0];
				UNPROTECT(1);
				break;
			}
			case 'T':
			{
				SEXP i = PROTECT(GET_SLOT(x, Matrix_iSym));
				xnnz = XLENGTH(i);
				if (cl[1] == 's') {
					SEXP j = PROTECT(GET_SLOT(x, Matrix_jSym));
					int *pi = INTEGER(i), *pj = INTEGER(j);
					R_xlen_t k = XLENGTH(i);
					xnnz *= 2;
					while (k--)
						if (*(pi++) == *(pj++))
							--xnnz;
					UNPROTECT(1);
				} else if (cl[1] == 't' && *CHAR(STRING_ELT(GET_SLOT(x, Matrix_diagSym), 0)) != 'N')
					xnnz += xdim[0];
				UNPROTECT(1);
				break;
			}
			case 'i':
				xnnz = xdim[0];
				break;
			case 'd':
				xnnz = XLENGTH(GET_SLOT(x, Matrix_permSym));
				break;
			default:
				break;
			}

			nnz += xnnz;
			len += xlen;
		}

		if (nnz > INT_MAX || nnz > len / 2)
			*repr = 'e';
		else if (anyCsparse && anyRsparse)
			*repr = (margin == 1) ? 'C' : 'R';
		else if (anyCsparse)
			*repr = 'C';
		else if (anyRsparse)
			*repr = 'R';
		else if (anyTsparse)
			*repr = 'T';
		else
			*repr = (margin == 1) ? 'C' : 'R';
	}

	return;
}

static SEXP bind(SEXP args, SEXP exprs, int margin, int level)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	const char *cl;
	SEXP a, e, x, tmp;
	int rdim[2], rdimnames[2], i, ivalid;
	char kind, repr;
	char rcl[] = "...Matrix";
	scanArgs(args, exprs, margin, level,
	         rdim, rdimnames, &kind, &repr);
	if (kind == '\0' || repr == '\0') {
		if (kind != repr)
			error(_("should never happen ..."));
		rcl[0] = 'i';
		rcl[1] = 'n';
		rcl[2] = 'd';
	} else {
		rcl[0] = kind;
		rcl[1] = 'g';
		rcl[2] = repr;
	}
	SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(rcl));

	if (repr == 'e') {
		if ((Matrix_int_fast64_t) rdim[0] * rdim[1] > R_XLEN_T_MAX)
			error(_("attempt to allocate vector of length exceeding R_XLEN_T_MAX"));
		/* TODO */
	} else if (repr == 'C') {
		/* TODO */
	} else if (repr == 'R') {
		/* TODO */
	} else if (repr == 'T') {
		/* TODO */
	} else {
		SEXP perm = PROTECT(allocVector(INTSXP, rdim[margin]));
		int *pperm = INTEGER(perm);
		R_xlen_t len;
		for (a = args; a != R_NilValue; a = CDR(a)) {
			x = CAR(a);
			if (x == R_NilValue)
				continue;
			tmp = GET_SLOT(x, Matrix_permSym);
			len = XLENGTH(tmp);
			Matrix_memcpy(pperm, INTEGER(tmp), len, sizeof(int));
			pperm += len;
		}
		SET_SLOT(res, Matrix_permSym, perm);
		UNPROTECT(1);
		if (margin == 1)
			INTEGER(GET_SLOT(res, Matrix_marginSym))[0] = 2;
	}

	SEXP dim = PROTECT(GET_SLOT(res, Matrix_DimSym));
	INTEGER(dim)[0] = rdim[0];
	INTEGER(dim)[1] = rdim[1];
	UNPROTECT(1);

	if (rdimnames[0] || rdimnames[1]) {
		Rprintf("rd[0] = %d, rd[1] = %d ... \n", rdimnames[0], rdimnames[1]);
		SEXP dimnames = PROTECT(GET_SLOT(res, Matrix_DimNamesSym)),
			marnames, s[2], s_;
		int len, pos = 0;
		if (rdimnames[margin])
			PROTECT(marnames = allocVector(STRSXP, rdim[margin]));

		for (a = args, e = exprs; a != R_NilValue; a = CDR(a), e = CDR(e)) {
			x = CAR(a);
			if (x == R_NilValue)
				continue;
			s[0] = s[1] = R_NilValue;
			if (IS_S4_OBJECT(x)) {
				ivalid = R_check_class_etc(x, valid);
				cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];
				tmp = GET_SLOT(x, Matrix_DimSym);
				len = INTEGER(tmp)[margin];
				tmp = GET_SLOT(x, Matrix_DimNamesSym);
				if (cl[1] == 's') {
					if ((s_ = VECTOR_ELT(tmp, 1)) != R_NilValue ||
					    (s_ = VECTOR_ELT(tmp, 0)) != R_NilValue)
						s[0] = s[1] = s_;
				} else
					for (i = 0; i < 2; ++i)
						s[i] = VECTOR_ELT(tmp, i);
			} else {
				tmp = getAttrib(x, R_DimSymbol);
				if (TYPEOF(tmp) == INTSXP && LENGTH(tmp) == 2) {
					len = INTEGER(tmp)[margin];
					tmp = getAttrib(x, R_DimNamesSymbol);
					if (tmp != R_NilValue)
						for (i = 0; i < 2; ++i)
							s[i] = VECTOR_ELT(tmp, i);
				} else if (rdim[!margin] == 0 || XLENGTH(x) > 0) {
					len = 1;
					if (rdim[!margin] > 0 && XLENGTH(x) == rdim[!margin])
						s[!margin] = getAttrib(x, R_NamesSymbol);
					if (TAG(a) != R_NilValue)
						s[margin] = coerceVector(TAG(a), STRSXP);
					else if (level == 2) {
						PROTECT(s_ = allocVector(EXPRSXP, 1));
						SET_VECTOR_ELT(s_, 0, CAR(e));
						s[margin] = coerceVector(s_, STRSXP);
						UNPROTECT(1);
					} else if (level == 1 && TYPEOF(CAR(e)) == SYMSXP)
						s[margin] = coerceVector(CAR(e), STRSXP);
				} else
					len = 0;
			}
			if (rdimnames[!margin] && s[!margin] != R_NilValue) {
				SET_VECTOR_ELT(dimnames, !margin, s[!margin]);
				if (!rdimnames[margin])
					break;
				rdimnames[!margin] = 1;
			}
			if (rdimnames[ margin] && s[ margin] != R_NilValue)
				for (i = 0; i < len; ++i)
					SET_STRING_ELT(marnames, pos + i, STRING_ELT(s[margin], i));
			pos += len;
		}

		if (rdimnames[margin]) {
			SET_VECTOR_ELT(dimnames, margin, marnames);
			UNPROTECT(1);
		}
		UNPROTECT(1);
	}

	UNPROTECT(1);
	return res;
}

SEXP R_bind(SEXP args)
{
	SEXP margin, level, exprs;
	args = CDR(args); margin = CAR(args);
	args = CDR(args);  level = CAR(args);
	args = CDR(args);  exprs = CAR(args);
	return bind(CDR(args), CDR(exprs), asInteger(margin), asInteger(level));
}
