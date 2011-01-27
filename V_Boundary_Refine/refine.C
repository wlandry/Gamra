/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#include "V_Boundary_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"

void SAMRAI::geom::V_Boundary_Refine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::SideOverlap* t_overlap =
      dynamic_cast<const pdat::SideOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   for(int axis=0; axis<getDim().getValue(); ++axis)
     {
       const hier::BoxList& boxes = t_overlap->getDestinationBoxList(axis);
       for (hier::BoxList::Iterator b(boxes); b; b++)
         {
           refine(fine,coarse,dst_component,src_component,b(),ratio,axis);
         }
     }
}

void SAMRAI::geom::V_Boundary_Refine::refine(hier::Patch& fine,
                                             const hier::Patch& coarse,
                                             const int dst_component,
                                             const int src_component,
                                             const hier::Box& overlap_box,
                                             const hier::IntVector& ratio,
                                             const int &axis) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine, coarse, overlap_box, ratio);

   tbox::Pointer<pdat::SideData<double> >
   v = coarse.getPatchData(src_component);
   tbox::Pointer<pdat::SideData<double> >
   v_fine = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!v.isNull());
   TBOX_ASSERT(!v_fine.isNull());
   TBOX_ASSERT(v->getDepth() == v_fine->getDepth());
   TBOX_ASSERT(v->getDepth() == 1);
#endif

   hier::Box coarse_box=coarse.getBox();
   hier::Box fine_box=fine.getBox();

   std::cout << "VBR "
             << fine.getPatchLevelNumber() << " "
             << axis << " "
             << coarse_box.lower(0) << " "
             << coarse_box.upper(0) << " "
             << coarse_box.lower(1) << " "
             << coarse_box.upper(1) << " "
             << fine_box.lower(0) << " "
             << fine_box.upper(0) << " "
             << fine_box.lower(1) << " "
             << fine_box.upper(1) << " "

             << overlap_box.lower(0) << " "
             << overlap_box.upper(0) << " "
             << overlap_box.lower(1) << " "
             << overlap_box.upper(1) << " "
             << "\n";

   /* We have to infer where the boundary is from the boxes */
   int boundary_direction;
   bool boundary_positive(false);
   if(std::abs(overlap_box.lower(0)-overlap_box.upper(0))==(axis==0 ? 1 : 0))
     {
       boundary_direction=0;
       if(fine_box.upper(0)<overlap_box.lower(0))
         boundary_positive=true;
       else if(fine_box.lower(0)>overlap_box.upper(0))
         boundary_positive=false;
       else
         abort();
     }
   else if(std::abs(overlap_box.lower(1)-overlap_box.upper(1))==
           (axis==1 ? 1 : 0))
     {
       boundary_direction=1;
       if(fine_box.upper(1)<overlap_box.lower(1))
         boundary_positive=true;
       else if(fine_box.lower(1)>overlap_box.upper(1))
         boundary_positive=false;
       else
         abort();
     }
   else
     {
       abort();
     }

   int i_min(overlap_box.lower(0)), i_max(overlap_box.upper(0)),
     j_min(overlap_box.lower(1)), j_max(overlap_box.upper(1));
   if(axis==0)
     {
       if(boundary_direction==0)
         {
           if(boundary_positive)
             {
               i_max=i_min;
             }
           else
             {
               i_min=i_max;
             }
           j_min=std::max(j_min,fine_box.lower(1));
           j_max=std::min(j_max,fine_box.upper(1));
         }
       /* We need to shrink the box because we do not want the edges.
          Those are points that are either covered by other patches or
          are ghost points that we do not care about.  */
       else
         {
           --i_max;
           ++i_min;
         }
     }
   else if(axis==1)
     {
       if(boundary_direction==1)
         {
           if(boundary_positive)
             {
               j_max=j_min;
             }
           else
             {
               j_min=j_max;
             }
           i_min=std::max(i_min,fine_box.lower(0));
           i_max=std::min(i_max,fine_box.upper(0));
         }
       else
         {
           --j_max;
           ++j_min;
         }
     }

   for(int j=j_min; j<=j_max; ++j)
     for(int i=i_min; i<=i_max; ++i)
       {
         pdat::SideIndex fine(hier::Index(i,j),axis,pdat::SideIndex::Lower);
         hier::Index ip(1,0), jp(0,1);

         if(axis==0)
           {
             Update_V(axis,boundary_direction,boundary_positive,fine,
                      ip,jp,j,v,v_fine);
           }
         else if(axis==1)
           {
             Update_V(axis,boundary_direction,boundary_positive,fine,
                      jp,ip,i,v,v_fine);
           }
       }
}
