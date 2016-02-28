/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Stokes/V_Boundary_Refine.hxx"
#include "Stokes/set_boundary.hxx"
#include "Constants.hxx"

void Stokes::V_Boundary_Refine::refine
(SAMRAI::hier::Patch& fine,
 const SAMRAI::hier::Patch& coarse,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::BoxOverlap& fine_overlap,
 const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::pdat::SideOverlap* t_overlap =
    dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap);

   Stokes_set_boundary(coarse,invalid_id,src_component,true);

   for(Gamra::Dir axis=0; axis<fine.getDim().getValue(); ++axis)
     {
       const SAMRAI::hier::BoxContainer& boxes
         = t_overlap->getDestinationBoxContainer(axis);
       for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
            b!=boxes.end(); ++b)
         {
           refine_box(fine,coarse,dst_component,src_component,*b,ratio,axis);
         }
     }
}
