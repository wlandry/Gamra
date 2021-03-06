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

void Stokes::FACSolver::setBcObject(const SAMRAI::solv::RobinBcCoefStrategy*
                                    bc_object)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  if (!bc_object) {
    TBOX_ERROR(d_object_name << ": NULL pointer for boundary condition\n"
               << "object.\n");
  }
#endif
  d_bc_object = bc_object;
  // d_fac_ops.setPhysicalBcCoefObject(d_bc_object);
}
