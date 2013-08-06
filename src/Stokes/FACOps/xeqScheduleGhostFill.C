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

void Stokes::FACOps::xeqScheduleGhostFill(int p_id, int v_id,
                                                      int dest_ln)
{
  /* p */
  {
    if (!p_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }
    SAMRAI::xfer::RefineAlgorithm refiner;
    refiner.registerRefine(p_id,p_id,p_id,p_ghostfill_refine_operator);
    refiner.resetSchedule(p_ghostfill_refine_schedules[dest_ln]);
    p_ghostfill_refine_schedules[dest_ln]->fillData(0.0,false);
  }

  /* v */
  {
    if (!v_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }
    set_boundaries(invalid_id,v_id,dest_ln-1);
    SAMRAI::xfer::RefineAlgorithm refiner;
    refiner.registerRefine(v_id,v_id,v_id,v_ghostfill_refine_operator);
    refiner.resetSchedule(v_ghostfill_refine_schedules[dest_ln]);
    v_ghostfill_refine_schedules[dest_ln]->fillData(0.0,false);
  }
}

