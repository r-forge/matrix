#include "lmer.h"
#include <R_ext/Rdynload.h>

/* Syms.h needs to be included a second time (it is already included
 * through lmer.h) so the symbols are defined without the extern modifier.
 */
#include "Syms.h" 

static R_CallMethodDef CallEntries[] = {
    {"ST_getPars", (DL_FUNC) &ST_getPars, 1},
    {"ST_initialize", (DL_FUNC) &ST_initialize, 3},
    {"ST_setPars", (DL_FUNC) &ST_setPars, 2},

    {"glmer_dev_resids", (DL_FUNC) &glmer_dev_resids, 1},
    {"glmer_dmu_deta", (DL_FUNC) &glmer_dmu_deta, 1},
    {"glmer_var", (DL_FUNC) &glmer_var, 1},

    {"lme4_rWishart", (DL_FUNC) &lme4_rWishart, 3},

    {"mer_MCMCsamp", (DL_FUNC) &mer_MCMCsamp, 6},
    {"mer_condMode", (DL_FUNC) &mer_condMode, 1},
    {"mer_create_L", (DL_FUNC) &mer_create_L, 1},
    {"mer_create_Vt", (DL_FUNC) &mer_create_Vt, 3},
    {"mer_optimize", (DL_FUNC) &mer_optimize, 2},
    {"mer_postVar", (DL_FUNC) &mer_postVar, 2},
    {"mer_sigma", (DL_FUNC) &mer_sigma, 2},
    {"mer_update_L", (DL_FUNC) &mer_update_L, 1},
    {"mer_update_Vt", (DL_FUNC) &mer_update_Vt, 1},
    {"mer_update_effects", (DL_FUNC) &mer_update_effects, 1},
    {"mer_update_eta", (DL_FUNC) &mer_update_eta, 1},
    {"mer_update_mu_res", (DL_FUNC) &mer_update_mu_res, 1},
    {"mer_update_lpdisc", (DL_FUNC) &mer_update_lpdisc, 1},
    {"mer_validate", (DL_FUNC) &mer_validate, 1},

    {"nlmer_create_Mt", (DL_FUNC) &nlmer_create_Mt, 2},

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
    lme4_dmu_detaSym = install("dmu_deta");
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
    lme4_residSym = install("resid");
    lme4_uvecSym = install("uvec");
    lme4_varSym = install("var");
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
