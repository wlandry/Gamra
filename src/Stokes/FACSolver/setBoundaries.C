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

    void Stokes::FACSolver::setBoundaries(const std::string& boundary_type,
                                          const int fluxes,
                                          const int flags,
                                          int* bdry_types)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_bc_object && d_bc_object != &d_simple_bc) {
        TBOX_ERROR(
                   d_object_name << ": Bad attempt to set boundary condition\n"
                   <<
                   "by using default bc object after it has been overriden.\n");
      }
#endif
      d_simple_bc.setBoundaries(boundary_type,
                                fluxes,
                                flags,
                                bdry_types);
      d_bc_object = &d_simple_bc;
      // d_fac_ops.setPhysicalBcCoefObject(d_bc_object);
    }


  }
}
