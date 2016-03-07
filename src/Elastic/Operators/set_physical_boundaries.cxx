/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Operators.hxx"

void Elastic::Operators::set_physical_boundaries
(const int &v_id,
 const SAMRAI::hier::PatchLevel &patch_level, const bool &rhs)
{
  for (SAMRAI::hier::PatchLevel::Iterator p(patch_level.begin());
       p!=patch_level.end(); ++p)
    { boundary_conditions.set_physical_boundary(**p,v_id,rhs); }
}
