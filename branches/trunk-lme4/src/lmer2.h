#ifndef LME4_LMER2_H
#define LME4_LMER2_H

#include "lme4_utils.h"
#include <R_ext/stats_package.h>

SEXP mer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		   SEXP verbose, SEXP deviance);
SEXP mer2_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp,
		SEXP nc, SEXP cnames, SEXP offset, SEXP wts);
SEXP lmer2_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp,
		  SEXP nc, SEXP cnames, SEXP offset, SEXP wts, SEXP fr,
		  SEXP terms, SEXP call, SEXP fam);
SEXP mer2_deviance(SEXP x, SEXP which);
SEXP mer2_getPars(SEXP x);
SEXP mer2_optimize(SEXP x, SEXP verb);
SEXP mer2_postVar(SEXP x);
SEXP mer2_ranef(SEXP x);
SEXP mer2_setPars(SEXP x, SEXP pars);
SEXP mer2_sigma(SEXP x, SEXP which);
SEXP mer2_update_effects(SEXP x);
SEXP mer2_validate(SEXP x);
SEXP mer2_vcov(SEXP x);

#endif /* LME4_LMER2_H */
