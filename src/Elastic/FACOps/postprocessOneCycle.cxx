/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

void Elastic::FACOps::postprocessOneCycle
(int fac_cycle_num,
 const SAMRAI::solv::SAMRAIVectorReal<double> &,
 const SAMRAI::solv::SAMRAIVectorReal<double> &)
{
  if (logging && d_preconditioner)
    {
      double avg_factor, final_factor;
      d_preconditioner->getConvergenceFactors(avg_factor, final_factor);
      SAMRAI::tbox::plog
        << "iter=" << std::setw(4) << fac_cycle_num
        << " resid=" << d_preconditioner->getResidualNorm()
        << " net conv=" << d_preconditioner->getNetConvergenceFactor()
        << " final conv=" << d_preconditioner->getNetConvergenceFactor()
        << " avg conv=" << d_preconditioner->getAvgConvergenceFactor()
        << std::endl;
    }
}
