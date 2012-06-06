/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "Stokes/FACOps.h"

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
#include "Stokes/HypreSolver.h"
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
* FACOperatorStrategy virtual postprocessOneCycle function.  *
********************************************************************
*/

void SAMRAI::solv::Stokes::FACOps::postprocessOneCycle
(int fac_cycle_num,
 const SAMRAIVectorReal<double>& current_soln,
 const SAMRAIVectorReal<double>& residual)
{
  NULL_USE(current_soln);
  NULL_USE(residual);

  if (d_enable_logging) {
    if (d_preconditioner) {
      /*
       * Output convergence progress.  This is probably only appropriate
       * if the solver is NOT being used as a preconditioner.
       */
      double avg_factor, final_factor;
      d_preconditioner->getConvergenceFactors(avg_factor, final_factor);
      tbox::plog
        << "iter=" << std::setw(4) << fac_cycle_num
        << " resid=" << d_preconditioner->getResidualNorm()
        << " net conv=" << d_preconditioner->getNetConvergenceFactor()
        << " final conv=" << d_preconditioner->getNetConvergenceFactor()
        << " avg conv=" << d_preconditioner->getAvgConvergenceFactor()
        << std::endl;
    }
  }
}
