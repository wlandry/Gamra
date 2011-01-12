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

#ifndef included_geom_V_Refine_C
#define included_geom_V_Refine_C

#include "V_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"

void SAMRAI::geom::V_Refine::refine(
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

void SAMRAI::geom::V_Refine::refine(hier::Patch& fine,
                                    const hier::Patch& coarse,
                                    const int dst_component,
                                    const int src_component,
                                    const hier::Box& fine_box,
                                    const hier::IntVector& ratio,
                                    const int &axis) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine, coarse, fine_box, ratio);

   tbox::Pointer<pdat::SideData<double> >
   v = coarse.getPatchData(src_component);
   tbox::Pointer<pdat::SideData<double> >
   v_fine = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!v.isNull());
   TBOX_ASSERT(!v_fine.isNull());
   TBOX_ASSERT(v->getDepth() == v_fine->getDepth());
#endif

   hier::Box coarse_box=coarse.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = coarse.getPatchGeometry();

   for(int j=fine_box.lower(1); j<=fine_box.upper(1); ++j)
     for(int i=fine_box.lower(0); i<=fine_box.upper(0); ++i)
       {
         pdat::SideIndex fine(hier::Index(i,j),axis,pdat::SideIndex::Lower);

         hier::Index ip(1,0), jp(0,1);
         pdat::SideIndex center(fine);
         center[0]/=2;
         center[1]/=2;

         if(axis==0)
           {
             double dvx_dy;

             if(i%2==0)
               {
                 if(center[1]==coarse_box.lower(1)
                    && geom->getTouchesRegularBoundary(1,0))
                   {
                     dvx_dy=((*v)(center+jp)-(*v)(center))/4;
                   }
                 else if(center[1]==coarse_box.upper(1)
                         && geom->getTouchesRegularBoundary(1,1))
                   {
                     dvx_dy=((*v)(center)-(*v)(center-jp))/4;
                   }
                 else
                   {
                     dvx_dy=((*v)(center+jp)-(*v)(center-jp))/8;
                   }

                 (*v_fine)(fine)=(*v)(center)
                   + ((j%2==0) ? (-dvx_dy) : dvx_dy);
               }
             else
               {
                 if(center[1]==coarse_box.lower(1)
                    && geom->getTouchesRegularBoundary(1,0))
                   {
                     dvx_dy=((*v)(center+jp)-(*v)(center)
                             + (*v)(center+ip+jp)-(*v)(center+ip))/4;
                   }
                 else if(center[1]==coarse_box.upper(1)
                         && geom->getTouchesRegularBoundary(1,1))
                   {
                     dvx_dy=((*v)(center)-(*v)(center-jp)
                             + (*v)(center+ip)-(*v)(center+ip-jp))/4;
                   }
                 else
                   {
                     dvx_dy=((*v)(center+jp)-(*v)(center-jp)
                             + (*v)(center+ip+jp)-(*v)(center+ip-jp))/8;
                   }

                 (*v_fine)(fine)=((*v)(center) + (*v)(center+ip)
                                  + ((j%2==0) ? (-dvx_dy) : dvx_dy))/2;
               }
           }
         else
           {
             double dvy_dx;

             if(j%2==0)
               {
                 if(center[0]==coarse_box.lower(0)
                    && geom->getTouchesRegularBoundary(0,0))
                   {
                     dvy_dx=((*v)(center+ip)-(*v)(center))/4;
                   }
                 else if(center[0]==coarse_box.upper(0)
                         && geom->getTouchesRegularBoundary(0,1))
                   {
                     dvy_dx=((*v)(center)-(*v)(center-ip))/4;
                   }
                 else
                   {
                     dvy_dx=((*v)(center+ip)-(*v)(center-ip))/8;
                   }

                 (*v_fine)(fine)=(*v)(center)
                   + ((i%2==0) ? (-dvy_dx) : dvy_dx);
               }
             else
               {
                 if(center[0]==coarse_box.lower(0)
                    && geom->getTouchesRegularBoundary(0,0))
                   {
                     dvy_dx=((*v)(center+ip)-(*v)(center)
                             + (*v)(center+jp+ip)-(*v)(center+jp))/4;
                   }
                 else if(center[0]==coarse_box.upper(0)
                         && geom->getTouchesRegularBoundary(0,1))
                   {
                     dvy_dx=((*v)(center)-(*v)(center-ip)
                             + (*v)(center+jp)-(*v)(center+jp-ip))/4;
                   }
                 else
                   {
                     dvy_dx=((*v)(center+ip)-(*v)(center-ip)
                             + (*v)(center+jp+ip)-(*v)(center+jp-ip))/8;
                   }

                 (*v_fine)(fine)=((*v)(center) + (*v)(center+jp)
                                  + ((i%2==0) ? (-dvy_dx) : dvy_dx))/2;
               }
           }
       }
}
#endif
