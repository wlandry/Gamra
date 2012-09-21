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

void Elastic::FACOps::xeqScheduleRRestriction(int v_dst, int v_src,
                                              int dest_ln)
{
  /* v */
  {
    if (!v_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
    }

    SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
    coarsener.registerCoarsen(v_dst,v_src,v_rrestriction_coarsen_operator);
    coarsener.resetSchedule(v_rrestriction_coarsen_schedules[dest_ln]);
    v_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
  }
}
