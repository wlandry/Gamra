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
* FACOperatorStrategy virtual computeResidualNorm             *
* function                                                         *
********************************************************************
*/

double SAMRAI::solv::StokesFACOps::computeResidualNorm
(const SAMRAIVectorReal<double>& residual,
 int fine_ln,
 int coarse_ln)
{
  if (coarse_ln != residual.getCoarsestLevelNumber() ||
      fine_ln != residual.getFinestLevelNumber()) {
    TBOX_ERROR("StokesFACOps::computeResidualNorm() is not\n"
               << "set up to compute residual except on the range of\n"
               << "levels defining the vector.\n");
  }
  t_compute_residual_norm->start();
  /*
   * The residual vector was cloned from vectors that has
   * the proper weights associated with them, so we do not
   * have to explicitly weight the residuals.
   *
   * maxNorm: not good to use because Hypre's norm does not
   *   correspond to it.  Also maybe too sensitive to spikes.
   * L2Norm: maybe good.  But does not correspond to the
   *   scale of the quantity.
   * L1Norm: maybe good.  Correspond to scale of quantity,
   *   but may be too insensitive to spikes.
   * RMSNorm: maybe good.
   */
  double norm = residual.RMSNorm();
  t_compute_residual_norm->stop();
  return norm;
}
