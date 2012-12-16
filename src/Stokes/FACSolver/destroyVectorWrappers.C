/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "Stokes/FACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
  namespace solv {

    /*
***********************************************************************
* Delete the vector wrappers.  Do not freeVectorComponents because    *
* we do not control their data allocation.  The user does that.       *
***********************************************************************
*/
    void Stokes::FACSolver::destroyVectorWrappers() {
      d_uv.reset();
      d_fv.reset();
    }

  }
}
