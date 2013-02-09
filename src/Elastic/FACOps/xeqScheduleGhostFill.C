/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for Elastic using FAC 
 *
 ************************************************************************/
#include "Elastic/FACOps.h"

void Elastic::FACOps::xeqScheduleGhostFill(int v_id, int dest_ln)
{
  /* v */
  if (!v_ghostfill_refine_schedules[dest_ln]) {
    TBOX_ERROR("Expected schedule not found.");
  }
  SAMRAI::xfer::RefineAlgorithm refiner(d_dim);

  refiner.registerRefine(v_id,v_id,v_id,v_ghostfill_refine_operator);
  refiner.registerRefine(dv_diagonal_id,dv_diagonal_id,dv_diagonal_id,
                         boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  refiner.registerRefine(dv_mixed_id,dv_mixed_id,dv_mixed_id,
                         boost::shared_ptr<SAMRAI::hier::RefineOperator>());
  refiner.resetSchedule(v_ghostfill_refine_schedules[dest_ln]);
  v_ghostfill_refine_schedules[dest_ln]->fillData(0.0,false);
}

