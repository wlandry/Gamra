/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"
#include "Constants.hxx"

void Elastic::Coarse_Fine_Boundary_Refine::refine
(SAMRAI::hier::Patch& fine_patch,
 const SAMRAI::hier::Patch& coarse,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::BoxOverlap& fine_overlap,
 const SAMRAI::hier::IntVector& ratio) const
{
   const SAMRAI::pdat::SideOverlap* t_overlap =
      dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);
   TBOX_ASSERT(t_overlap != NULL);
   for(Gamra::Dir ix=0; ix<fine_patch.getDim().getValue(); ++ix)
     {
       const SAMRAI::hier::BoxContainer&
         boxes = t_overlap->getDestinationBoxContainer(ix);
       for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
            b!=boxes.end(); ++b)
         {
           refine_box(fine_patch,coarse,dst_component,src_component,*b,ratio,ix);
         }
     }
}
