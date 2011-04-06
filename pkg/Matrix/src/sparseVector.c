#include "Mutils.h"

/* Define all of
 *  dsparseVector_sub(....)
 *  isparseVector_sub(....)
 *  lsparseVector_sub(....)
 *  nsparseVector_sub(....)
 *  zsparseVector_sub(....)
 */

#define _dspV_
#include "t_sparseVector.c"

#define _ispV_
#include "t_sparseVector.c"

#define _lspV_
#include "t_sparseVector.c"

#define _nspV_
#include "t_sparseVector.c"

#define _zspV_
#include "t_sparseVector.c"
