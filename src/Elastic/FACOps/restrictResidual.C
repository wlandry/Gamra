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
#include "Elastic/HypreSolver.h"
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
* FACOperatorStrategy virtual restrictresidual function.     *
********************************************************************
*/

void SAMRAI::solv::Elastic::FACOps::restrictResidual
(const SAMRAIVectorReal<double>& s,
 SAMRAIVectorReal<double>& d,
 int dest_ln)
{
  t_restrict_residual->start();

  int p_src(s.getComponentDescriptorIndex(0)),
    p_dst(d.getComponentDescriptorIndex(0)),
    v_src(s.getComponentDescriptorIndex(1)),
    v_dst(d.getComponentDescriptorIndex(1));

  /* Need to do a sync because the coarsening for v uses ghost zones */
  v_coarsen_patch_strategy.setSourceDataId(v_src);
  xeqScheduleGhostFillNoCoarse(invalid_id,v_src,dest_ln+1);

  xeqScheduleRRestriction(p_dst,p_src,v_dst,v_src,dest_ln);

  t_restrict_residual->stop();
}
