/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

int Elastic::FACOps::solveCoarsestLevel
(SAMRAI::solv::SAMRAIVectorReal<double>& data,
 const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int coarsest_ln)
{
  t_solve_coarsest->start();

  d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
  smoothError(data, residual, coarsest_ln, d_coarse_solver_max_iterations);
  d_residual_tolerance_during_smoothing = -1.0;
  xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
                               coarsest_ln);
  t_solve_coarsest->stop();
  return 0;
}
