/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"

void Elastic::FAC::resetHierarchyConfiguration
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &new_hierarchy,
 int , int )
{
  d_hierarchy = new_hierarchy;
}
