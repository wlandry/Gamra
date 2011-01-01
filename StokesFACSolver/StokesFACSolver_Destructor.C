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
#include "StokesFACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
  namespace solv {


    /*
*************************************************************************
*                                                                       *
* Destructor for StokesFACSolver.                            *
* Deallocate internal data.                                             *
*                                                                       *
*************************************************************************
*/

    StokesFACSolver::~StokesFACSolver()
    {
      s_instance_counter[d_dim.getValue() - 1]--;

      deallocateSolverState();

      if (s_instance_counter[d_dim.getValue() - 1] == 0) {
        hier::VariableDatabase::getDatabase()->
          removeInternalSAMRAIVariablePatchDataIndex(s_weight_id[d_dim.getValue() - 1]);
        s_weight_id[d_dim.getValue() - 1] = -1;
      }
    }

  }
}
