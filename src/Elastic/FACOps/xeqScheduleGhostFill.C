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

void SAMRAI::solv::Elastic::FACOps::xeqScheduleGhostFill(int v_id,
                                                         int dest_ln)
{
  /* v */
  {
    if (!v_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }
    set_boundaries(v_id,dest_ln-1);
    xfer::RefineAlgorithm refiner(d_dim);
    refiner.registerRefine(v_id,v_id,v_id,v_ghostfill_refine_operator);
    refiner.resetSchedule(v_ghostfill_refine_schedules[dest_ln]);
    v_ghostfill_refine_schedules[dest_ln]->fillData(0.0,false);
  }
}

