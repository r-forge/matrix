#include "Mdefines.h"

const char *valid_dense            [] = { VALID_DENSE            , "" };
const char *valid_sparse           [] = { VALID_SPARSE           , "" };
const char *valid_sparse_compressed[] = { VALID_SPARSE_COMPRESSED, "" };
const char *valid_sparse_triplet   [] = { VALID_SPARSE_TRIPLET   , "" };
const char *valid_diagonal         [] = { VALID_DIAGONAL         , "" };
const char *valid_index            [] = { VALID_INDEX            , "" };
const char *valid_matrix           [] = { VALID_MATRIX           , "" };
const char *valid_vector           [] = { VALID_VECTOR           , "" };
const char *valid_matrix_or_vector [] = { VALID_MATRIX_OR_VECTOR , "" };

SEXP newObject(const char *what)
{
	SEXP class = PROTECT(R_do_MAKE_CLASS(what)), obj = R_do_new_object(class);
	UNPROTECT(1);
	return obj;
}

char typeToKind(SEXPTYPE type)
{
	switch (type) {
	case LGLSXP:
		return 'l';
	case INTSXP:
		return 'i';
	case REALSXP:
		return 'd';
	case CPLXSXP:
		return 'z';
	default:
		Rf_error(_("unexpected type \"%s\" in '%s'"),
		         Rf_type2char(type), __func__);
		return '\0';
	}
}

SEXPTYPE kindToType(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
		return LGLSXP;
	case 'i':
		return INTSXP;
	case 'd':
		return REALSXP;
	case 'z':
		return CPLXSXP;
	default:
		Rf_error(_("unexpected kind \"%c\" in '%s'"), kind, __func__);
		return NILSXP;
	}
}

size_t kindToSize(char kind)
{
	switch (kind) {
	case 'n':
	case 'l':
	case 'i':
		return sizeof(int);
	case 'd':
		return sizeof(double);
	case 'z':
		return sizeof(Rcomplex);
	default:
		Rf_error(_("unexpected kind \"%c\" in '%s'"), kind, __func__);
		return 0;
	}
}

const char *Matrix_superclass(const char *class, int mode)
{
	if (class[0] == 'p') {
		if (mode & 1)
			return "indMatrix";
	} else if (class[1] == 'o') {
		if (mode & 4)
			switch (class[2]) {
			case 'r': return "dsyMatrix";
			case 'p': return "dspMatrix";
			}
		if (mode & 2)
			switch (class[2]) {
			case 'r': return "dpoMatrix";
			case 'p': return "dppMatrix";
			}
	} else if (class[1] == 'p') {
		if (mode & 4) {
			if (class[0] == 'z')
			switch (class[2]) {
			case 'C': return "zsCMatrix";
			case 'R': return "zsRMatrix";
			case 'T': return "zsTMatrix";
			case 'o': return "zsyMatrix";
			case 'p': return "zspMatrix";
			}
			else
			switch (class[2]) {
			case 'C': return "dsCMatrix";
			case 'R': return "dsRMatrix";
			case 'T': return "dsTMatrix";
			case 'o': return "dsyMatrix";
			case 'p': return "dspMatrix";
			}
		}
	}
	return class;
}

const char *Matrix_class(SEXP x, const char **valid, int mode,
                         const char *caller)
{
	int i = R_check_class_etc(x, valid);
	if (i >= 0)
		return (mode <= 0) ? valid[i] : Matrix_superclass(valid[i], mode);
	else {
		if (caller)
			ERROR_INVALID_CLASS(x, caller);
		return NULL;
	}
}

char Matrix_kind(SEXP obj)
{
	if (TYPEOF(obj) != OBJSXP)
		switch (TYPEOF(obj)) {
		case LGLSXP:
			return 'l';
		case INTSXP:
			return 'i';
		case REALSXP:
			return 'd';
		case CPLXSXP:
			return 'z';
		default:
			return '\0';
		}
	const char *class = Matrix_class(obj, valid_matrix_or_vector, 7, NULL);
	if (!class)
		return '\0';
	return (class[2] == 'd') ? 'n' : class[0];
}

char Matrix_shape(SEXP obj, int mode)
{
	if (TYPEOF(obj) != OBJSXP)
		return '\0';
	const char *class = Matrix_class(obj, valid_matrix_or_vector, mode, NULL);
	if (!class)
		return '\0';
	return (class[2] == 'd') ? 'i' : (class[3] == 'M') ? class[1] : 'g';
}

char Matrix_repr(SEXP obj)
{
	if (TYPEOF(obj) != OBJSXP)
		return '\0';
	const char *class = Matrix_class(obj, valid_matrix_or_vector, 7, NULL);
	if (!class)
		return '\0';
	switch (class[2]) {
	case 'e':
	case 'y':
	case 'r':
		return 'n'; /* unpackedMatrix */
	case 'p':
		return 'p'; /*   packedMatrix */
	case 'C':
		return 'C'; /*  CsparseMatrix */
	case 'R':
		return 'R'; /*  RsparseMatrix */
	case 'T':
		return 'T'; /*  TsparseMatrix */
	case 'i':
		return 'd'; /* diagonalMatrix */
	case 'd':
		return 'i'; /*      indMatrix */
	default:
		return '\0';
	}
}

SEXP R_Matrix_class(SEXP s_obj, SEXP s_mode)
{
	const char *class = Matrix_class(s_obj, valid_matrix_or_vector, Rf_asInteger(s_mode), NULL);
	return Rf_mkString((!class) ? "" : class);
}

SEXP R_Matrix_kind(SEXP s_obj)
{
	char s[] = { Matrix_kind (s_obj), '\0' };
	return Rf_mkString(s);
}

SEXP R_Matrix_shape(SEXP s_obj, SEXP s_mode)
{
	char s[] = { Matrix_shape(s_obj, Rf_asInteger(s_mode)), '\0' };
	return Rf_mkString(s);
}

SEXP R_Matrix_repr(SEXP s_obj)
{
	char s[] = { Matrix_repr (s_obj), '\0' };
	return Rf_mkString(s);
}
