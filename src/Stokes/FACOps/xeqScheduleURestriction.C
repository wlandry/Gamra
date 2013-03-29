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

void Stokes::FACOps::xeqScheduleURestriction(int p_dst, int p_src,
                                                           int v_dst, int v_src,
                                                           int dest_ln)
{
  /* p */
  {
    if (!p_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }

    SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
    coarsener.registerCoarsen(p_dst, p_src, p_urestriction_coarsen_operator);
    coarsener.resetSchedule(p_urestriction_coarsen_schedules[dest_ln]);
    p_urestriction_coarsen_schedules[dest_ln]->coarsenData();
  }

  /* v */
  {
    if (!v_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }

    SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
    coarsener.registerCoarsen(v_dst, v_src, v_urestriction_coarsen_operator);
    coarsener.resetSchedule(v_urestriction_coarsen_schedules[dest_ln]);
    v_urestriction_coarsen_schedules[dest_ln]->coarsenData();
  }
}
