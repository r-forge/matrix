#ifndef LME4_LMER2_H
#define LME4_LMER2_H

#include "lme4_utils.h"
#include <R_ext/stats_package.h>

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
SEXP lmer_update_L(SEXP x);
SEXP lmer_update_effects(SEXP x);
SEXP lmer_update_dev(SEXP x);
SEXP lmer_validate(SEXP x);

SEXP mer_create_L(SEXP Vt);
SEXP mer_create_Vt(SEXP Zt, SEXP ST, SEXP Gp);
SEXP mer_optimize(SEXP x, SEXP verb, SEXP mtype);
SEXP mer_postVar(SEXP x);
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

#endif /* LME4_LMER2_H */
