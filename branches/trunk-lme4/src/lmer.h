#ifndef LME4_LMER_H
#define LME4_LMER_H

#include "lme4_utils.h"
#include "Wishart.h"
				/* positions in the deviance vector */
enum devP {ML_POS=0, REML_POS, ldZ_POS, ldX_POS, lr2_POS};
			 /* {"ML", "REML", "ldZ", "ldX", "lr2", ""} */
				/* positions in the dims vector */
enum dimP {nf_POS=0, n_POS, p_POS, q_POS, isREML_POS, isGLMM_POS};
	      /* {"nf", "n", "p", "q", "isREML", "isGLMM", ""} */

#define isREML(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isREML_POS]
#define isGLMM(x) INTEGER(GET_SLOT(x, lme4_dimsSym))[isGLMM_POS]

SEXP mer_ECMEsteps(SEXP x, SEXP nsteps, SEXP Verbp);
SEXP mer_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		  SEXP verbose, SEXP deviance);
SEXP mer2_MCMCsamp(SEXP x, SEXP savebp, SEXP nsampp, SEXP transp,
		   SEXP verbose, SEXP deviance);
SEXP mer_coef(SEXP x, SEXP pType);
SEXP mer_coefGets(SEXP x, SEXP coef, SEXP pType);
SEXP mer_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp, SEXP nc, SEXP cnames);
SEXP mer2_create(SEXP fl, SEXP Zt, SEXP Xp, SEXP yp, SEXP REMLp,
		SEXP nc, SEXP cnames, SEXP offset, SEXP wts);
SEXP mer_dtCMatrix(SEXP x);
SEXP mer_dtCMatrix_inv(SEXP x);
SEXP mer2_deviance(SEXP x, SEXP which);
SEXP mer_fitted(SEXP x);
SEXP mer_fixef(SEXP x);
SEXP mer2_getPars(SEXP x);
SEXP mer_gradient(SEXP x, SEXP pType);
SEXP mer_hat_trace(SEXP x);
SEXP mer_hat_trace2(SEXP x);
SEXP mer_initial(SEXP x);
SEXP mer_isNested(SEXP x);
SEXP mer_postVar(SEXP x);
SEXP mer2_postVar(SEXP x);
SEXP mer_ranef(SEXP x);
SEXP mer2_ranef(SEXP x);
SEXP mer2_setPars(SEXP x, SEXP pars);
SEXP mer_sigma(SEXP x, SEXP REML);
SEXP mer2_sigma(SEXP x, SEXP which);
SEXP mer_simulate(SEXP x, SEXP nsimP);
SEXP mer2_update_effects(SEXP x);
SEXP mer_update_ZXy(SEXP x);
SEXP mer_update_y(SEXP x, SEXP ynew);
SEXP mer_validate(SEXP x);
SEXP mer2_vcov(SEXP x);

SEXP Ztl_sparse(SEXP fl, SEXP Ztl);
SEXP Zt_carryOver(SEXP f, SEXP Zt, SEXP tvar, SEXP discount);

#endif /* LME4_LMER_H */
