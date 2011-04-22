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

void SAMRAI::solv::StokesFACOps::xeqScheduleGhostFillNoCoarse(int p_id,
                                                              int v_id,
                                                              int dest_ln)
{
  /* p */
  if(p_id!=invalid_id)
  {
    if (!p_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected cell schedule not found.");
    }
    xfer::RefineAlgorithm refiner(d_dim);
    refiner.registerRefine(p_id,p_id,p_id,tbox::Pointer<xfer::RefineOperator>(0));
    refiner.resetSchedule(p_nocoarse_refine_schedules[dest_ln]);
    p_nocoarse_refine_schedules[dest_ln]->fillData(0.0,false);
  }

  /* v */
  if(v_id!=invalid_id)
  {
    if (!v_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected side schedule not found.");
    }
    xfer::RefineAlgorithm refiner(d_dim);
    refiner.registerRefine(v_id,v_id,v_id,tbox::Pointer<xfer::RefineOperator>(0));
    refiner.resetSchedule(v_nocoarse_refine_schedules[dest_ln]);
    v_nocoarse_refine_schedules[dest_ln]->fillData(0.0,false);
  }
}
