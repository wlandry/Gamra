/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.h"

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
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

/*
********************************************************************
* FACOperatorStrategy virtual solveCoarsestLevel             *
* function                                                         *
********************************************************************
*/

int Elastic::FACOps::solveCoarsestLevel
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int coarsest_ln)
{
  t_solve_coarsest->start();

  int return_value = 0;

  if (d_coarse_solver_choice == "Tackley"
      || d_coarse_solver_choice == "Gerya") {
    d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
    smoothError(data,
                residual,
                coarsest_ln,
                d_coarse_solver_max_iterations);
    d_residual_tolerance_during_smoothing = -1.0;
  } else {
    TBOX_ERROR(
               d_object_name << ": Bad coarse level solver choice '"
               << d_coarse_solver_choice
               <<
               "' in Elastic::FACOps::solveCoarsestLevel.");
  }

  xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
                               coarsest_ln);
  t_solve_coarsest->stop();

  return return_value;
}
