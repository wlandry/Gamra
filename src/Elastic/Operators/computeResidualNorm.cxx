/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

double Elastic::Operators::computeResidualNorm
(const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int fine_level,
 int coarse_level)
{
  if (coarse_level != residual.getCoarsestLevelNumber() ||
      fine_level != residual.getFinestLevelNumber())
    { TBOX_ERROR("Elastic::Operators::computeResidualNorm() is not\n"
                 << "set up to compute residual except on the range of\n"
                 << "levels defining the vector.\n"); }

  t_compute_residual_norm->start();
  /// We use maxNorm because it is a definite upper bound on the error.
  double norm = residual.maxNorm();
  t_compute_residual_norm->stop();
  return norm;
}
