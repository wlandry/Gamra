/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::deallocateOperatorState()
{
  d_cf_boundary.clear();
  d_ln_min = -1;
  d_ln_max = -1;

  refine_schedules.clear();
  coarsen_solution_schedules.clear();
  coarsen_resid_schedules.clear();
  ghostfill_schedules.clear();
  ghostfill_nocoarse_schedules.clear();
}
