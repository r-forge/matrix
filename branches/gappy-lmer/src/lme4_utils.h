#ifndef LME4_UTILS_H
#define LME4_UTILS_H
#include <R_ext/Constants.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rversion.h>
#include <Rconfig.h>
#include "Matrix.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
#else
# define attr_hidden
#endif

extern
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

extern cholmod_common c;

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

#endif
