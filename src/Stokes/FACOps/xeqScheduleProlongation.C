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

void Stokes::FACOps::xeqScheduleProlongation
(int p_dst, int p_src, int p_scr, int v_dst, int v_src, int v_scr,
 int dest_ln)
{
  /* p */
  {
    if (!p_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }
    SAMRAI::xfer::RefineAlgorithm refiner(d_dim);
    refiner.registerRefine(p_dst, p_src, p_scr, p_prolongation_refine_operator);
    refiner.resetSchedule(p_prolongation_refine_schedules[dest_ln]);
    p_prolongation_refine_schedules[dest_ln]->fillData(0.0,false);
  }

  /* v */
  {
    if (!v_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }
    SAMRAI::xfer::RefineAlgorithm refiner(d_dim);
    refiner.registerRefine(v_dst, v_src, v_scr, v_prolongation_refine_operator);
    refiner.resetSchedule(v_prolongation_refine_schedules[dest_ln]);
    v_prolongation_refine_schedules[dest_ln]->fillData(0.0,false);
  }
}
