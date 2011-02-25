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

/*
********************************************************************
* FACOperatorStrategy virtual solveCoarsestLevel             *
* function                                                         *
********************************************************************
*/

int SAMRAI::solv::StokesFACOps::solveCoarsestLevel
(SAMRAIVectorReal<double>& data,
 const SAMRAIVectorReal<double>& residual,
 int coarsest_ln)
{
  t_solve_coarsest->start();

  checkInputPatchDataIndices();

  int return_value = 0;

  if (d_coarse_solver_choice == "Tackley"
      || d_coarse_solver_choice == "Gerya") {
    d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
    smoothError(data,
                residual,
                coarsest_ln,
                d_coarse_solver_max_iterations);
    d_residual_tolerance_during_smoothing = -1.0;
  } else if (d_coarse_solver_choice == "hypre") {
#ifndef HAVE_HYPRE
    TBOX_ERROR(d_object_name << ": Coarse level solver choice '"
               << d_coarse_solver_choice
               << "' unavailable in "
               << "scapStokesOps::solveCoarsestLevel.");
#else
    return_value = solveCoarsestLevel_HYPRE(data, residual, coarsest_ln);
#endif
  } else {
    TBOX_ERROR(
               d_object_name << ": Bad coarse level solver choice '"
               << d_coarse_solver_choice
               <<
               "' in scapStokesOps::solveCoarsestLevel.");
  }

  xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
                               data.getComponentDescriptorIndex(1),
                               coarsest_ln);
  t_solve_coarsest->stop();

  return return_value;
}
