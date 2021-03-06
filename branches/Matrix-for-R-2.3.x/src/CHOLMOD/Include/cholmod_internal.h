/* ========================================================================== */
/* === Include/cholmod_internal.h =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_internal.h.  Version 0.6.
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_internal.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* CHOLMOD internal include file.
 *
 * This file contains internal definitions for CHOLMOD, not meant to be included
 * in user code.  They define macros that are not prefixed with CHOLMOD_.  This
 * file can safely #include'd in user code if you want to make use of the
 * macros defined here, and don't mind the possible name conflicts with your
 * code, however.
 *
 * Required by all CHOLMOD routines.  Not required by any user routine that
 * uses CHOLMOMD.  Unless debugging is enabled, this file does not require any
 * CHOLMOD module (not even the Core module).
 *
 * If debugging is enabled, all CHOLMOD modules require the Check module.
 * Enabling debugging requires that this file be editted.  Debugging cannot be
 * enabled with a compiler flag.  This is because CHOLMOD is exceedingly slow
 * when debugging is enabled.  Debugging is meant for development of CHOLMOD
 * itself, not by users of CHOLMOD.
 */

#ifndef CHOLMOD_INTERNAL_H
#define CHOLMOD_INTERNAL_H

/* turn off debugging */
#ifndef NDEBUG
#define NDEBUG
#endif

/* Uncomment this line to enable debugging.  CHOLMOD will be very slow.
#undef NDEBUG
 */

#if !defined(NPRINT) || !defined(NDEBUG)
#include <stdio.h>
#endif

#include <stddef.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

/* ========================================================================== */
/* === basic definitions ==================================================== */
/* ========================================================================== */

/* Some non-conforming compilers insist on defining TRUE and FALSE. */
#undef TRUE
#undef FALSE
#define TRUE 1
#define FALSE 0
#define BOOLEAN(x) ((x) ? TRUE : FALSE)

/* NULL should already be defined, but ensure it is here. */
#ifndef NULL
#define NULL ((void *) 0)
#endif

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP (i) : (i))

/* MAX and MIN are not safe to use for NaN's */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MAX3(a,b,c) (((a) > (b)) ? (MAX (a,c)) : (MAX (b,c)))
#define MAX4(a,b,c,d) (((a) > (b)) ? (MAX3 (a,c,d)) : (MAX3 (b,c,d)))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define IMPLIES(p,q) (!(p) || (q))

/* find the sign: -1 if x < 0, 1 if x > 0, zero otherwise.
 * Not safe for NaN's */
#define SIGN(x) (((x) < 0) ? (-1) : (((x) > 0) ? 1 : 0))

/* round up an integer x to a multiple of s */
#define ROUNDUP(x,s) ((s) * (((x) + ((s) - 1)) / (s)))

#define ERROR(status,msg) \
    CHOLMOD(error) (status, __FILE__, __LINE__, msg, Common)

/* Check a pointer and return if null.  Set status to invalid, unless the
 * status is already "out of memory" */
#define RETURN_IF_NULL(A,result) \
{ \
    if ((A) == NULL) \
    { \
	if (Common->status != CHOLMOD_OUT_OF_MEMORY) \
	{ \
	    ERROR (CHOLMOD_INVALID, "argument missing") ; \
	} \
	return (result) ; \
    } \
}

/* Return if Common is NULL or invalid */
#define RETURN_IF_NULL_COMMON(result) \
{ \
    if (Common == NULL) \
    { \
	return (result) ; \
    } \
    if (Common->itype != ITYPE || Common->dtype != DTYPE) \
    { \
	Common->status = CHOLMOD_INVALID ; \
	return (result) ; \
    } \
}

#define IS_NAN(x)	CHOLMOD_IS_NAN(x)
#define IS_ZERO(x)	CHOLMOD_IS_ZERO(x)
#define IS_NONZERO(x)	CHOLMOD_IS_NONZERO(x)
#define IS_LT_ZERO(x)	CHOLMOD_IS_LT_ZERO(x)
#define IS_GT_ZERO(x)	CHOLMOD_IS_GT_ZERO(x)
#define IS_LE_ZERO(x)	CHOLMOD_IS_LE_ZERO(x)

/* ========================================================================== */
/* === int/long and double/float definitions ================================ */
/* ========================================================================== */

/* CHOLMOD is designed for 3 types of integer variables:
 *
 *	(1) all integers are int
 *	(2) most integers are int, some are long
 *	(3) all integers are long
 *
 * and two kinds of floating-point values:
 *
 *	(1) double
 *	(2) float
 *
 * the complex types (ANSI-compatible complex, and MATLAB-compatable zomplex)
 * are based on the double or float type, and are not selected here.  They
 * are typically selected via template routines.
 *
 * This gives 6 different modes in which CHOLMOD can be compiled (only the
 * first two are currently supported):
 *
 *	DINT	double, int		prefix: cholmod_
 *	DLONG	double, long		prefix: cholmod_l_
 *	DMIX	double, mixed int/long	prefix: cholmod_m_
 *	SINT	float, int		prefix: cholmod_si_
 *	SLONG	float, long		prefix: cholmod_sl_
 *	SMIX	float, mixed int/log	prefix: cholmod_sm_
 *
 * These are selected with compile time flags (-DDLONG, for example).  If no
 * flag is selected, the default is DINT.
 *
 * All six versions use the same include files.  The user-visible include files
 * are completely independent of which int/long/double/float version is being
 * used.  The integer / real types in all data structures (sparse, triplet,
 * dense, common, and triplet) are defined at run-time, not compile-time, so
 * there is only one "cholmod_sparse" data type.  Void pointers are used inside
 * that data structure to point to arrays of the proper type.  Each data
 * structure has an itype and dtype field which determines the kind of basic
 * types used.  These are defined in Include/cholmod_core.h.
 *
 * FUTURE WORK: support all six types (float, and mixed int/long)
 */

/* -------------------------------------------------------------------------- */
/* double, long */
/* -------------------------------------------------------------------------- */

#ifdef DLONG
#define Real double
#define Int long
#define Int_max LONG_MAX
#define Size_max LONG_MAX
#define CHOLMOD(name) cholmod_l_ ## name
#define LONG
#define DOUBLE
#define ITYPE CHOLMOD_LONG
#define DTYPE CHOLMOD_DOUBLE
#define ID "%ld"

/* -------------------------------------------------------------------------- */
/* double, int/long */
/* -------------------------------------------------------------------------- */

#elif defined (DMIX)
#error "mixed int/long not yet supported"

/* -------------------------------------------------------------------------- */
/* single, int */
/* -------------------------------------------------------------------------- */

#elif defined (SINT)
#error "single-precision not yet supported"

/* -------------------------------------------------------------------------- */
/* single, long */
/* -------------------------------------------------------------------------- */

#elif defined (SLONG)
#error "single-precision not yet supported"

/* -------------------------------------------------------------------------- */
/* single, int/long */
/* -------------------------------------------------------------------------- */

#elif defined (SMIX)
#error "single-precision not yet supported"

/* -------------------------------------------------------------------------- */
/* double, int: this is the default */
/* -------------------------------------------------------------------------- */

#else

#ifndef DINT
#define DINT
#endif
#define INT
#define DOUBLE

#define Real double
#define Int int
#define Int_max INT_MAX
#define Size_max INT_MAX
#define CHOLMOD(name) cholmod_ ## name
#define ITYPE CHOLMOD_INT
#define DTYPE CHOLMOD_DOUBLE
#define ID "%d"

#endif


/* ========================================================================== */
/* === real/complex arithmetic ============================================== */
/* ========================================================================== */

#include "cholmod_complexity.h"

/* ========================================================================== */
/* === Architecture and BLAS ================================================ */
/* ========================================================================== */

#include "cholmod_blas.h"

/* ========================================================================== */
/* === debugging definitions ================================================ */
/* ========================================================================== */

#ifndef NDEBUG

#include <assert.h>
#include "cholmod.h"

/* The cholmod_dump routines are in the Check module.  No CHOLMOD routine
 * calls the cholmod_check_* or cholmod_print_* routines in the Check module,
 * since they use Common workspace that may already be in use.  Instead, they
 * use the cholmod_dump_* routines defined there, which allocate their own
 * workspace if they need it. */

/* double, int */
extern int cholmod_dump ;
extern int cholmod_dump_malloc ;
long cholmod_dump_sparse (cholmod_sparse  *, char *, cholmod_common *) ;
int  cholmod_dump_factor (cholmod_factor  *, char *, cholmod_common *) ;
int  cholmod_dump_triplet (cholmod_triplet *, char *, cholmod_common *) ;
int  cholmod_dump_dense (cholmod_dense   *, char *, cholmod_common *) ;
int  cholmod_dump_subset (int *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_dump_perm (int *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_dump_parent (int *, size_t, char *, cholmod_common *) ;
void cholmod_dump_init (char *, cholmod_common *) ;
int  cholmod_dump_mem (char *, long, cholmod_common *) ;
void cholmod_dump_real (char *, Real *, long, long, int, int,
	cholmod_common *) ;
void cholmod_dump_super (long, int *, int *, int *, int *, double *, int,
	cholmod_common *) ;
int  cholmod_dump_partition (long, int *, int *, int *, int *, long,
	cholmod_common *) ;
int  cholmod_dump_work(int, int, long, cholmod_common *) ;

/* double, long */
extern int cholmod_l_dump ;
extern int cholmod_l_dump_malloc ;
long cholmod_l_dump_sparse (cholmod_sparse  *, char *, cholmod_common *) ;
int  cholmod_l_dump_factor (cholmod_factor  *, char *, cholmod_common *) ;
int  cholmod_l_dump_triplet (cholmod_triplet *, char *, cholmod_common *) ;
int  cholmod_l_dump_dense (cholmod_dense   *, char *, cholmod_common *) ;
int  cholmod_l_dump_subset (long *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_l_dump_perm (long *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_l_dump_parent (long *, size_t, char *, cholmod_common *) ;
void cholmod_l_dump_init (char *, cholmod_common *) ;
int  cholmod_l_dump_mem (char *, long, cholmod_common *) ;
void cholmod_l_dump_real (char *, Real *, long, long, int, int,
	cholmod_common *) ;
void cholmod_l_dump_super (long, long *, long *, long *, long *, double *, int,
	cholmod_common *) ;
int  cholmod_l_dump_partition (long, long *, long *, long *, long *, long,
	cholmod_common *) ;
int  cholmod_l_dump_work(int, int, long, cholmod_common *) ;

#define DEBUG_INIT(s)  { CHOLMOD(dump_init)(s, Common) ; }
#define ASSERT(expression) (assert (expression))

#define PRK(k,params) \
{ \
    if (CHOLMOD(dump) >= (k) && Common->print_function != NULL) \
    { \
	/* (Common->print_function) */ printf params ; \
    } \
}

#define PRINT0(params) PRK (0, params)
#define PRINT1(params) PRK (1, params)
#define PRINT2(params) PRK (2, params)
#define PRINT3(params) PRK (3, params)
#define PRINTM(params) PRK (CHOLMOD(dump_malloc), params)

#define DEBUG(statement) statement

#else

/* Debugging disabled (the normal case) */
#define PRK(k,params)
#define DEBUG_INIT(s)
#define PRINT0(params)
#define PRINT1(params)
#define PRINT2(params)
#define PRINT3(params)
#define PRINTM(params)
#define ASSERT(expression)
#define DEBUG(statement)
#endif

#endif
