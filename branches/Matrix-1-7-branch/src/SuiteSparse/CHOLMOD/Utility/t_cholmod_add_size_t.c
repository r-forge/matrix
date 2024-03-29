//------------------------------------------------------------------------------
// CHOLMOD/Utility/t_cholmod_add_size_t: add two size_t values
//------------------------------------------------------------------------------

// CHOLMOD/Utility Module. Copyright (C) 2023, Timothy A. Davis, All Rights
// Reserved.
// SPDX-License-Identifier: LGPL-2.1+

//------------------------------------------------------------------------------

#include "cholmod_internal.h"

size_t CHOLMOD(add_size_t) (size_t a, size_t b, int *ok)
{

    //--------------------------------------------------------------------------
    // add a and b
    //--------------------------------------------------------------------------

    size_t s = a + b ;

    //--------------------------------------------------------------------------
    // check for size_t overflow
    //--------------------------------------------------------------------------

    if (s < a || s < b)
    {
        (*ok) = FALSE ;
        s = 0 ;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    return (s) ;
}

