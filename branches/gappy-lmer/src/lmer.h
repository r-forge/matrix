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

#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#endif

extern
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif


SEXP ST_getPars(SEXP x);
SEXP ST_initialize(SEXP ST, SEXP Gp, SEXP Zt);
SEXP ST_setPars(SEXP x, SEXP pars);

SEXP glmer_dev_resids(SEXP x);
SEXP glmer_dmu_deta(SEXP x);
SEXP glmer_var(SEXP x);

SEXP lme4_rWishart(SEXP ns, SEXP dfp, SEXP scal);

SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		   SEXP verbose, SEXP deviance);
SEXP mer_condMode(SEXP x);
SEXP mer_create_L(SEXP Vt);
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP Gp);
SEXP mer_eta(SEXP x);
SEXP mer_optimize(SEXP x, SEXP verb);
SEXP mer_postVar(SEXP x, SEXP useScale);
SEXP mer_sigma(SEXP x, SEXP which);
SEXP mer_update_L(SEXP x);
SEXP mer_update_Vt(SEXP x);
SEXP mer_update_effects(SEXP x);
SEXP mer_update_eta(SEXP x);
SEXP mer_update_mu_res(SEXP x);
SEXP mer_update_lpdisc(SEXP x);
SEXP mer_validate(SEXP x);

SEXP nlmer_create_Mt(SEXP Vt, SEXP sP);

SEXP pedigree_chol(SEXP x, SEXP ans);

#endif /* LME4_LMER_H */
