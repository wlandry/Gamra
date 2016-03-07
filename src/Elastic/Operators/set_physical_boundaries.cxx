/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::set_physical_boundaries
(const int &v_id,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel> &patch_level, const bool &rhs)
{
  for (SAMRAI::hier::PatchLevel::Iterator pi(patch_level->begin());
       pi!=patch_level->end(); ++pi)
    {
      boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;
      boundary_conditions.set_physical_boundary(*patch,v_id,rhs);
    }
}
