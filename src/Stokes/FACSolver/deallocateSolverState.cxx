/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
#include <SAMRAI/pdat/CellVariable.h>
#include "Stokes/FACSolver.hxx"
#include <SAMRAI/tbox/PIO.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/tbox/StartupShutdownManager.h>

#include IOMANIP_HEADER_FILE

void Stokes::FACSolver::deallocateSolverState()
{
  if (d_hierarchy)
    {
      d_fac_precond.deallocateSolverState();
      d_hierarchy.reset();
      d_level_min = -1;
      d_level_max = -1;
      d_solver_is_initialized = false;
      destroyVectorWrappers();
    }
}
