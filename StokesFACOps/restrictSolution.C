/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "StokesFACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "StokesHypreSolver.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

namespace SAMRAI {
  namespace solv {

    /*
********************************************************************
* FACOperatorStrategy virtual restrictSolution function.     *
* After restricting solution, update ghost cells of the affected   *
* level.                                                           *
********************************************************************
*/

    void StokesFACOps::restrictSolution(
                                        const SAMRAIVectorReal<double>& s,
                                        SAMRAIVectorReal<double>& d,
                                        int dest_ln) {

      t_restrict_solution->start();

      xeqScheduleURestriction(d.getComponentDescriptorIndex(0),
                              s.getComponentDescriptorIndex(0),
                              dest_ln);

      d_bc_helper.setHomogeneousBc(false);
      d_bc_helper.setTargetDataId(d.getComponentDescriptorIndex(0));

      if (dest_ln == d_ln_min) {
        // xeqScheduleGhostFillNoCoarse(d.getComponentDescriptorIndex(0),
        //                              dest_ln);
        abort();
      } else {
        // xeqScheduleGhostFill(d.getComponentDescriptorIndex(0),
        //                      dest_ln);
        abort();
      }

      t_restrict_solution->stop();
    }

  }
}
