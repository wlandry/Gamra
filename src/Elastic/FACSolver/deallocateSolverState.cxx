/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

#include "Elastic/FACSolver.hxx"

void Elastic::FACSolver::deallocateSolverState()
{
  if (d_hierarchy)
    {
      d_fac_precond.deallocateSolverState();
      for (int ln = d_ln_min; ln <= d_ln_max; ++ln)
        d_hierarchy->getPatchLevel(ln)->deallocatePatchData
          (s_weight_id[d_dim.getValue() - 1]);

      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;
      d_solver_is_initialized = false;
      destroyVectorWrappers();
    }
}
