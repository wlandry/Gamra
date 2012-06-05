/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar Elastic equation. 
 *
 ************************************************************************/
#include "SAMRAI/pdat/CellVariable.h"
#include "ElasticFACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
  namespace solv {

    void ElasticFACSolver::initializeStatics() {

      for (int d = 0; d < SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
        s_weight_id[d] = -1;
        s_instance_counter[d] = -1;
      }

      s_initialized = 1;
    }

  }
}
