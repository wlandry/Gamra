/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

int Elastic::Operators::solveCoarsestLevel
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int coarsest_ln)
{
  t_solve_coarsest->start();

  smooth(data, residual, coarsest_ln, coarse_solver_max_iterations,
         coarse_solver_tolerance);
  ghostfill_nocoarse(data.getComponentDescriptorIndex(0), coarsest_ln);
  t_solve_coarsest->stop();
  return 0;
}
