#include "Elastic/FACOps.hxx"

double Elastic::FACOps::computeResidualNorm
(const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int fine_ln,
 int coarse_ln)
{
  if (coarse_ln != residual.getCoarsestLevelNumber() ||
      fine_ln != residual.getFinestLevelNumber())
    TBOX_ERROR("Elastic::FACOps::computeResidualNorm() is not\n"
               << "set up to compute residual except on the range of\n"
               << "levels defining the vector.\n");

  t_compute_residual_norm->start();
  /// We use maxNorm because it is a definite upper bound on the error.
  double norm = residual.maxNorm();
  t_compute_residual_norm->stop();
  return norm;
}
