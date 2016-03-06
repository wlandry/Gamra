/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

void Elastic::FACSolver::deallocateSolverState()
{
  if (d_hierarchy)
    {
      preconditioner.deallocateSolverState();
      d_hierarchy.reset();
      level_min = -1;
      level_max = -1;
      d_solver_is_initialized = false;
      destroyVectorWrappers();
    }
}
