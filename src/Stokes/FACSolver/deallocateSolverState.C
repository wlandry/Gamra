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

void Stokes::FACSolver::deallocateSolverState()
{
  if (d_hierarchy) {

    d_fac_precond.deallocateSolverState();

    /*
     * Delete internally managed data.
     */
    int ln;
    for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->deallocatePatchData(s_weight_id[d_dim.getValue()
                                                                      - 1]);
    }

    d_hierarchy.reset();
    d_ln_min = -1;
    d_ln_max = -1;
    d_solver_is_initialized = false;

    destroyVectorWrappers();

  }
}
