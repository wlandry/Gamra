/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
#include "Stokes/FACOps.hxx"

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
#include "Stokes/HypreSolver.hxx"
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

/*
********************************************************************
********************************************************************
*/

void Stokes::FACOps::smoothError
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int level,
 int num_sweeps)
{
  t_smooth_error->start();

  if (d_smoothing_choice == "Gerya") {
    smooth_Gerya(data,residual,level,num_sweeps,
                 d_residual_tolerance_during_smoothing);
  } else if (d_smoothing_choice == "Tackley") {
    if(d_dim.getValue()==2)
      smooth_Tackley_2D(data,residual,level,num_sweeps,
                        d_residual_tolerance_during_smoothing);
    else if(d_dim.getValue()==3)
      smooth_Tackley_3D(data,residual,level,num_sweeps,
                        d_residual_tolerance_during_smoothing);
  } else {
    TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
               << d_smoothing_choice
               << "' in Stokes::FACOps.");
  }

  t_smooth_error->stop();
}
