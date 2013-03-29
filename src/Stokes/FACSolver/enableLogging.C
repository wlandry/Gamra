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

/*
*************************************************************************
* Enable logging and propagate logging flag to major components.        *
*************************************************************************
*/

void Stokes::FACSolver::enableLogging(bool logging)
{
  d_enable_logging = logging;
  d_fac_precond.enableLogging(d_enable_logging);
  d_fac_ops.enableLogging(d_enable_logging);
}
