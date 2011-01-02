/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#ifndef included_solv_StokesFACOps_C
#define included_solv_StokesFACOps_C

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

    void
    StokesFACOps::xeqScheduleGhostFill(
                                       int dst_id,
                                       int dest_ln)
    {
      if (!d_ghostfill_refine_schedules[dest_ln]) {
        TBOX_ERROR("Expected schedule not found.");
      }
      xfer::RefineAlgorithm refiner(d_dim);
      refiner.
        registerRefine(dst_id,
                       dst_id,
                       dst_id,
                       d_ghostfill_refine_operator);
      refiner.
        resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
      d_ghostfill_refine_schedules[dest_ln]->fillData(0.0);
      d_ghostfill_refine_algorithm->
        resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
    }

  }
}
#endif
