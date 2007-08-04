#ifndef LME4_LMER_H
#define LME4_LMER_H

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rversion.h>
#include <R_ext/Lapack.h>
#include <R_ext/stats_package.h>
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


/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SET_SLOT(obj, nm, allocVector(type, length));
    return GET_SLOT(obj, nm);
}


extern cholmod_common c;

SEXP ST_getPars(SEXP x);
SEXP ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP ST_setPars(SEXP x, SEXP pars);

SEXP glmer_condMode(SEXP x);
SEXP glmer_dev_resids(SEXP x);
SEXP glmer_eta(SEXP x);
SEXP glmer_linkinv(SEXP x);
SEXP glmer_validate(SEXP x);

SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal);

SEXP lmer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		   SEXP verbose, SEXP deviance);
SEXP lmer_update_dev(SEXP x);
SEXP lmer_update_effects(SEXP x);
SEXP lmer_validate(SEXP x);

SEXP mer_create_L(SEXP Vt);
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP Gp);
SEXP mer_optimize(SEXP x, SEXP verb, SEXP mtype);
SEXP mer_postVar(SEXP x, SEXP useScale);
SEXP mer_sigma(SEXP x, SEXP which);
SEXP mer_update_b(SEXP x);
SEXP mer_update_Vt(SEXP x);
SEXP mer_validate(SEXP x);

SEXP nlmer_condMode(SEXP x);
SEXP nlmer_create_Mt(SEXP Vt, SEXP sP);
SEXP nlmer_eval_model(SEXP x);
SEXP nlmer_update_Mt(SEXP x);
/* SEXP nlmer_update_wrkres(SEXP x); */
SEXP nlmer_validate(SEXP x);

SEXP pedigree_chol(SEXP x, SEXP ans);

#endif /* LME4_LMER_H */
