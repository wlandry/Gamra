/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACOps.hxx"

/// Set physical boundaries

void Elastic::FACOps::set_physical_boundaries
(const int &v_id,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel> &level, const bool &rhs)
{
  for (SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
       pi!=level->end(); ++pi)
    {
      boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;
      d_boundary_conditions.set_physical_boundary(*patch,v_id,rhs);
    }
}
