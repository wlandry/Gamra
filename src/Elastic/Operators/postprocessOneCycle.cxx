/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::postprocessOneCycle
(int fac_cycle_num,
 const SAMRAI::solv::SAMRAIVectorReal<double> &,
 const SAMRAI::solv::SAMRAIVectorReal<double> &)
{
  if (logging && preconditioner)
    {
      double avg_factor, final_factor;
      preconditioner->getConvergenceFactors(avg_factor, final_factor);
      SAMRAI::tbox::plog
        << "iter=" << std::setw(4) << fac_cycle_num
        << " resid=" << preconditioner->getResidualNorm()
        << " net conv=" << preconditioner->getNetConvergenceFactor()
        << " final conv=" << preconditioner->getNetConvergenceFactor()
        << " avg conv=" << preconditioner->getAvgConvergenceFactor()
        << std::endl;
    }
}
