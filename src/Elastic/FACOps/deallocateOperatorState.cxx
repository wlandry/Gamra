/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::deallocateOperatorState()
{
  if (initialized)
    {
      d_cf_boundary.clear();
      d_ln_min = -1;
      d_ln_max = -1;

      v_prolongation_refine_schedules.setNull();
      v_urestriction_coarsen_schedules.setNull();
      v_rrestriction_coarsen_schedules.setNull();
      v_ghostfill_refine_schedules.setNull();
      v_nocoarse_refine_schedules.setNull();
    }
  initialized=false;
}
