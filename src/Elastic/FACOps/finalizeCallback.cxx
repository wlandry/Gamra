/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.hxx"

#include IOMANIP_HEADER_FILE

#include <SAMRAI/hier/BoundaryBoxUtils.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Index.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/pdat/CellDoubleConstantRefine.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/OutersideData.h>
#include <SAMRAI/pdat/OutersideVariable.h>
#include <SAMRAI/hier/PatchData.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/solv/FACPreconditioner.h>
#include <SAMRAI/tbox/Array.h>
#include <SAMRAI/tbox/MathUtilities.h>
#include <SAMRAI/tbox/StartupShutdownManager.h>
#include <SAMRAI/tbox/Timer.h>
#include <SAMRAI/tbox/TimerManager.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/tbox/MathUtilities.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/PatchLevelFullFillPattern.h>

void
Elastic::FACOps::finalizeCallback()
{
  for (int d = 0; d < SAMRAI::MAX_DIM_VAL; ++d) {
    s_cell_scratch_var[d].reset();
    s_side_scratch_var[d].reset();
  }
}
