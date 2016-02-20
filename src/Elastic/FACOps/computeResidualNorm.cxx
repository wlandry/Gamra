#include "Elastic/FACOps.hxx"

double Elastic::FACOps::computeResidualNorm
(const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
 int fine_ln,
 int coarse_ln)
{
  if (coarse_ln != residual.getCoarsestLevelNumber() ||
      fine_ln != residual.getFinestLevelNumber()) {
    TBOX_ERROR("Elastic::FACOps::computeResidualNorm() is not\n"
               << "set up to compute residual except on the range of\n"
               << "levels defining the vector.\n");
  }
  t_compute_residual_norm->start();
  /*
   * The residual vector was cloned from vectors that has
   * the proper weights associated with them, so we do not
   * have to explicitly weight the residuals.
   *
   * maxNorm: not good to use because Hypre's norm does not
   *   correspond to it.  Also maybe too sensitive to spikes.
   * L2Norm: maybe good.  But does not correspond to the
   *   scale of the quantity.
   * L1Norm: maybe good.  Correspond to scale of quantity,
   *   but may be too insensitive to spikes.
   * RMSNorm: maybe good.
   *
   * We use maxNorm because it is a definite upper bound on the error.
   */
  double norm = residual.maxNorm();
  t_compute_residual_norm->stop();
  return norm;
}
