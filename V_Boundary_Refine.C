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

#ifndef included_geom_V_Boundary_Refine_C
#define included_geom_V_Boundary_Refine_C

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
   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = coarse.getPatchGeometry();

   double dx = geom->getDx()[0];
   double dy = geom->getDx()[1];

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
   // hier::Index lower_offset(0,0), upper_offset(0,0);
   // if(axis==boundary_direction)
   //   {
   //     if(boundary_sign==1)
   //       upper_offset[axis]=-1;
   //     else
   //       lower_offset[axis]=1;
   //   }

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
         pdat::SideIndex center(fine);
         center[0]/=2;
         center[1]/=2;

         /* Note that at boundaries that are not in the same direction
            as the axis, we store the derivative in the variable */
         if(axis==0)
           {
             if(boundary_direction==0)
               {
                 /* Interpolate in the y direction */
                 double dv_plus, dv_minus;

                 if(center[1]==coarse_box.lower(1)
                    && geom->getTouchesRegularBoundary(1,0))
                   {
                     dv_plus=(*v)(center+jp)-(*v)(center);
                     dv_minus=(*v)(center-jp)*dy;
                   }
                 else if(center[1]==coarse_box.upper(1)
                         && geom->getTouchesRegularBoundary(1,1))
                   {
                     dv_plus=(*v)(center+jp)*dy;
                     dv_minus=(*v)(center)-(*v)(center-jp);
                   }
                 else
                   {
                     dv_plus=(*v)(center+jp)-(*v)(center);
                     dv_minus=(*v)(center)-(*v)(center-jp);
                   }

                 double v_plus=(*v)(center)
                   + (5.0/32)*dv_plus - (3.0/32)*dv_minus;
                 double v_minus=(*v)(center)
                   + (5.0/32)*dv_minus - (3.0/32)*dv_plus;

                 (*v_fine)(fine)=v_minus*(2*(*v)(center))/(v_plus + v_minus);
                 (*v_fine)(fine+jp)=v_minus*(2*(*v)(center))/(v_plus + v_minus);
                 ++j;
               }
             else
               {
                 /* We are computing derivatives here */
                 if(i%2==0)
                   {
                     (*v_fine)(fine)=((*v)(center+jp) - (*v)(center))/dy;
                   }
                 else
                   {
                     (*v_fine)(fine)=
                       0.5*((*v)(center+jp) - (*v)(center)
                            + (*v)(center+jp+ip) - (*v)(center+ip))/dy;
                   }
               }
           }
         else if(axis==1)
           {
             if(boundary_direction==1)
               {
                 /* Interpolate in the x direction */
                 double dv_plus, dv_minus;

                 if(center[0]==coarse_box.lower(0)
                    && geom->getTouchesRegularBoundary(0,0))
                   {
                     dv_plus=(*v)(center+ip)-(*v)(center);
                     dv_minus=(*v)(center-ip)*dx;
                   }
                 else if(center[0]==coarse_box.upper(0)
                         && geom->getTouchesRegularBoundary(0,1))
                   {
                     dv_plus=(*v)(center+ip)*dx;
                     dv_minus=(*v)(center)-(*v)(center-ip);
                   }
                 else
                   {
                     dv_plus=(*v)(center+ip)-(*v)(center);
                     dv_minus=(*v)(center)-(*v)(center-ip);
                   }

                 double v_plus=(*v)(center)
                   + (5.0/32)*dv_plus - (3.0/32)*dv_minus;
                 double v_minus=(*v)(center)
                   + (5.0/32)*dv_minus - (3.0/32)*dv_plus;

                 (*v_fine)(fine)=v_minus*(2*(*v)(center))/(v_plus + v_minus);
                 (*v_fine)(fine+ip)=v_minus*(2*(*v)(center))/(v_plus + v_minus);
                 ++i;
               }
             else
               {
                 /* We are computing derivatives here */
                 if(j%2==0)
                   {
                     (*v_fine)(fine)=((*v)(center+ip) - (*v)(center))/dx;
                   }
                 else
                   {
                     (*v_fine)(fine)=
                       0.5*((*v)(center+ip) - (*v)(center)
                            + (*v)(center+ip+jp) - (*v)(center+jp))/dx;
                   }
               }
           }
         else
           {
             abort();
           }
       }
}
#endif
