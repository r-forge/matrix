#include "lme4_utils.h"
#include "lmer.h"
#include "pedigree.h" 
#include <R_ext/Rdynload.h>

/* Syms.h needs to be included a second time (it is already included
 * through lme4_utils.h) so the symbols are defined without extern.
 */
#include "Syms.h" 

static R_CallMethodDef CallEntries[] = {
    {"ST_getPars", (DL_FUNC) &ST_getPars, 1},
    {"ST_initialize", (DL_FUNC) &ST_initialize, 3},
    {"ST_setPars", (DL_FUNC) &ST_setPars, 2},

    {"glmer_condMode", (DL_FUNC) &glmer_condMode, 1},
    {"glmer_dev_resids", (DL_FUNC) &glmer_dev_resids, 1},
    {"glmer_eta", (DL_FUNC) &glmer_eta, 1},
    {"glmer_linkinv", (DL_FUNC) &glmer_linkinv, 1},
    {"glmer_validate", (DL_FUNC) &glmer_validate, 1},

    {"lme4_rWishart", (DL_FUNC) &lme4_rWishart, 3},

    {"lmer_MCMCsamp", (DL_FUNC) &lmer_MCMCsamp, 6},
    {"lmer_update_L", (DL_FUNC) &lmer_update_L, 1},
    {"lmer_update_effects", (DL_FUNC) &lmer_update_effects, 1},
    {"lmer_update_dev", (DL_FUNC) &lmer_update_dev, 1},
    {"lmer_validate", (DL_FUNC) &lmer_validate, 1},

    {"mer_create_L", (DL_FUNC) &mer_create_L, 1},
    {"mer_create_Vt", (DL_FUNC) &mer_create_Vt, 4},
    {"mer_optimize", (DL_FUNC) &mer_optimize, 3},
    {"mer_postVar", (DL_FUNC) &mer_postVar, 1},
    {"mer_sigma", (DL_FUNC) &mer_sigma, 2},
    {"mer_update_Vt", (DL_FUNC) &mer_update_Vt, 1},
    {"mer_validate", (DL_FUNC) &mer_validate, 1},

    {"nlmer_bhat", (DL_FUNC) &nlmer_bhat, 1},
    {"nlmer_eval_model", (DL_FUNC) &nlmer_eval_model, 2},
    {"nlmer_update_Mt", (DL_FUNC) &nlmer_update_Mt, 1},
    {"nlmer_update_Vt", (DL_FUNC) &nlmer_update_Vt, 1},
    {"nlmer_update_ranef", (DL_FUNC) &nlmer_update_ranef, 1},
    {"nlmer_update_wrkres", (DL_FUNC) &nlmer_update_wrkres, 1},
    {"nlmer_validate", (DL_FUNC) &nlmer_validate, 1},

    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},

/*     {"Zt_carryOver", (DL_FUNC) &Zt_carryOver, 4}, */

    {NULL, NULL, 0}
};

cholmod_common c;

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void R_init_lme4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);


    M_R_cholmod_start(&c);
    c.final_ll = 1;	    /* LL' form of simplicial factorization */

    lme4_DimNamesSym = install("Dimnames");
    lme4_DimSym = install("Dim");
    lme4_GpSym = install("Gp");
    lme4_LSym = install("L");
    lme4_MtSym = install("Mt");
    lme4_RVXySym = install("RVXy");
    lme4_RXySym = install("RXy");
    lme4_STSym = install("ST");
    lme4_VtSym = install("Vt");
    lme4_XSym = install("X");
    lme4_XtSym = install("Xt");
    lme4_XytXySym = install("XytXy");
    lme4_ZtSym = install("Zt");
    lme4_ZtXySym = install("ZtXy");
    lme4_cnamesSym = install("cnames");
    lme4_devResidSym = install("devResid");
    lme4_devianceSym = install("deviance");
    lme4_dimsSym = install("dims");
    lme4_envSym = install("env");
    lme4_etaSym = install("eta");
    lme4_famNameSym = install("famName");
    lme4_fixefSym = install("fixef");
    lme4_flistSym = install("flist");
    lme4_gradientSym = install("gradient");
    lme4_iSym = install("i");
    lme4_modelSym = install("model");
    lme4_muSym = install("mu");
    lme4_offsetSym = install("offset");
    lme4_pSym = install("p");
    lme4_permSym = install("perm");
    lme4_pnamesSym = install("pnames");
    lme4_pwtsSym = install("pwts");
    lme4_ranefSym = install("ranef");
    lme4_uvecSym = install("uvec");
    lme4_weightsSym = install("weights");
    lme4_wrkresSym = install("wrkres");
    lme4_wtsSym = install("wts");
    lme4_xSym = install("x");
    lme4_ySym = install("y");
    lme4_zSym = install("z");
}

void R_unload_lme4(DllInfo *dll){
    M_cholmod_finish(&c);
}
