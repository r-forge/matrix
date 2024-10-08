#ifndef MATRIX_LAPACK_ETC_H
#define MATRIX_LAPACK_ETC_H

#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#ifdef PR18534fixed
# define usePR18534fix
#endif
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

#include <R_ext/Error.h>

#define ERROR_LAPACK_1(_ROUTINE_, _INFO_) \
do { \
	if ((_INFO_) < 0) \
		Rf_error(_("LAPACK routine '%s': argument %d had illegal value"), \
		         #_ROUTINE_, -(_INFO_)); \
} while (0)

#define ERROR_LAPACK_2(_ROUTINE_, _INFO_, _WARN_, _LETTER_) \
do { \
	ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
	if ((_INFO_) > 0 && (_WARN_) > 0) { \
		if (_WARN_ > 1) \
			Rf_error  (_("LAPACK routine '%s': matrix is exactly singular, %s[i,i]=0, i=%d"), \
			           #_ROUTINE_, #_LETTER_, (_INFO_)); \
		else \
			Rf_warning(_("LAPACK routine '%s': matrix is exactly singular, %s[i,i]=0, i=%d"), \
			           #_ROUTINE_, #_LETTER_, (_INFO_)); \
	} \
} while (0)

#define ERROR_LAPACK_3(_ROUTINE_, _INFO_, _WARN_) \
do { \
	ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
	if ((_INFO_) > 0 && (_WARN_) > 0) { \
		if (_WARN_ > 1) \
			Rf_error  (_("LAPACK routine '%s': leading principal minor of order %d is not positive"), \
			           #_ROUTINE_, (_INFO_)); \
		else \
			Rf_warning(_("LAPACK routine '%s': leading principal minor of order %d is not positive"), \
			           #_ROUTINE_, (_INFO_)); \
	} \
} while (0)

#define ERROR_LAPACK_4(_ROUTINE_, _INFO_, _WARN_, _RANK_) \
	do { \
		ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
		if ((_INFO_) > 0 && (_WARN_) > 0) { \
			if (_WARN_ > 1) \
				Rf_error  (_("LAPACK routine '%s': matrix is rank deficient or not positive definite, the _computed_ rank is %d"), \
				           #_ROUTINE_, (_RANK_)); \
			else \
				Rf_warning(_("LAPACK routine '%s': matrix is rank deficient or not positive definite, the _computed_ rank is %d"), \
				           #_ROUTINE_, (_RANK_)); \
		} \
	} while (0)

#define ERROR_LAPACK_5(_ROUTINE_, _INFO_, _WARN_) \
	do { \
		ERROR_LAPACK_1(_ROUTINE_, _INFO_); \
		if ((_INFO_) > 0 && (_WARN_) > 0) { \
			if (_WARN_ > 1) \
				Rf_error  (_("LAPACK routine '%s': failed to compute all eigenvalues"), \
				           #_ROUTINE_); \
			else \
				Rf_warning(_("LAPACK routine '%s': failed to compute all eigenvalues"), \
				           #_ROUTINE_); \
		} \
	} while (0)

#endif /* MATRIX_LAPACK_ETC_H */
