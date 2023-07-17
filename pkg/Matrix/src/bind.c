/* MJ: work in progress ... see also ../R/bind2.R */

static void scanArgs(SEXP args, int margin, int *rdim, char *kind, char *repr)
{
	static const char *valid[] = { VALID_NONVIRTUAL_MATRIX, "" };
	const char *cl;
	SEXP a, x, tmp;
	int nS4 = 0, nDense, nEmpty = 0, ivalid, *xdim,
		anyCsparse = 0, anyRsparse = 0, anyTsparse = 0, anyDiagonal = 0,
		anyN = 0, anyL = 0, anyI = 0, anyD = 0, anyZ = 0;

	*rdim[!margin] = -1;
	*rdim[ margin] =  0;

	for (a = CDR(args); a != R_NilValue; a = CDR(args)) {
		x = CAR(a);
		if (IS_S4_OBJECT(x)) {
			++nS4;
			ivalid = R_check_class_etc(x, valid);
			if (ivalid < 0) {
				if (margin == 0)
					ERROR_INVALID_CLASS(x, "rbind.Matrix");
				else
					ERROR_INVALID_CLASS(x, "cbind.Matrix");
			}

			tmp = GET_SLOT(x, Matrix_DimSym);
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

			cl = valid[ivalid + VALID_NONVIRTUAL_SHIFT(ivalid, 1)];
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
			} else {
				if (XLENGTH(x) > 0)
					rdim[margin] += 1;
				else
					++nEmpty;
			}
		}
	}
	if (rdim[!margin] == 0 && nEmpty > 0) {
		if (nEmpty > INT_MAX - rdim[margin])
			error(_("dimensions cannot exceed 2^31-1"));
		rdim[margin] += nEmpty;
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
		Matrix_int_fast64_t nnz = 0, len = 0, xnnz, xlen;
		for (a = CDR(args); a != R_NilValue && nnz < INT_MAX; a = CDR(args)) {
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
					SEXP i = PROTECT(GET_SLOT(x, Matrix_iSym));
					int *pi = INTEGER(i), j;
					xnnz *= 2;
					if (*CHAR(STRING_ELT(GET_SLOT(x, Matrix_uploSym), 0)) == 'U')
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j + 1] - 1] == j)
								--xnnz;
					else
						for (j = 0; j < n; ++j)
							if (pp[j] < pp[j + 1] && pi[pp[j]] == j)
								--xnnz;
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
