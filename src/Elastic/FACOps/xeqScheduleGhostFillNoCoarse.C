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

void Elastic::FACOps::xeqScheduleGhostFillNoCoarse(int v_id, int dest_ln)
{
  if (!v_nocoarse_refine_schedules[dest_ln]) {
    TBOX_ERROR("Expected side schedule not found.");
  }
  SAMRAI::xfer::RefineAlgorithm refiner;
  refiner.registerRefine(v_id,v_id,v_id,
                         boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  refiner.resetSchedule(v_nocoarse_refine_schedules[dest_ln]);
  v_nocoarse_refine_schedules[dest_ln]->fillData(0.0,false);
}
