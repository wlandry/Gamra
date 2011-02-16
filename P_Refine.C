/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for cell-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef included_geom_P_Refine_C
#define included_geom_P_Refine_C

#include "P_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

void SAMRAI::geom::P_Refine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      dynamic_cast<const pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxList& boxes = t_overlap->getDestinationBoxList();
   for (hier::BoxList::Iterator b(boxes); b; b++) {
      refine(fine,
         coarse,
         dst_component,
         src_component,
         b(),
         ratio);
   }
}

void SAMRAI::geom::P_Refine::refine(
   hier::Patch& fine_patch,
   const hier::Patch& coarse_patch,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine_patch, coarse_patch,
                                   fine_box, ratio);

   tbox::Pointer<pdat::CellData<double> >
   p = coarse_patch.getPatchData(src_component);
   tbox::Pointer<pdat::CellData<double> >
   p_fine = fine_patch.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!p.isNull());
   TBOX_ASSERT(!p_fine.isNull());
   TBOX_ASSERT(p->getDepth() == p_fine->getDepth());
#endif

   hier::Box coarse_box=coarse_patch.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = coarse_patch.getPatchGeometry();

   for(int j=fine_box.lower(1); j<=fine_box.upper(1); ++j)
     for(int i=fine_box.lower(0); i<=fine_box.upper(0); ++i)
       {
         pdat::CellIndex fine(tbox::Dimension(2));
         fine[0]=i;
         fine[1]=j;

         pdat::CellIndex ip(fine), jp(fine);
         ip[0]=1;
         ip[1]=0;
         jp[0]=0;
         jp[1]=1;
         pdat::CellIndex center(hier::Index::coarsen(fine,hier::Index(2,2)));
         double dp_dx,dp_dy;

         /* Pressure is cell-centered, so prolongation is a
            linear interpolation from nearby cells. */

         /* This assumes that the levels are always properly nested,
            so that we always have an extra grid space for
            interpolation.  So we only have to have a special case for
            physical boundaries, where we do not have an extra grid
            space. */

         /* We could, in theory, use a refine_patch_strategy to fill
            in the ghost zones with extrapolations, and then use
            simple differences here. */

         if(center[0]==coarse_box.lower(0)
            && geom->getTouchesRegularBoundary(0,0))
           {
             dp_dx=((*p)(center+ip)-(*p)(center))/4;
           }
         else if(center[0]==coarse_box.upper(0)
                 && geom->getTouchesRegularBoundary(0,1))
           {
             dp_dx=((*p)(center)-(*p)(center-ip))/4;
           }
         else
           {
             dp_dx=((*p)(center+ip)-(*p)(center-ip))/8;
           }

         if(center[1]==coarse_box.lower(1)
            && geom->getTouchesRegularBoundary(1,0))
           {
             dp_dy=((*p)(center+jp)-(*p)(center))/4;
           }
         else if(center[1]==coarse_box.upper(1)
                 && geom->getTouchesRegularBoundary(1,1))
           {
             dp_dy=((*p)(center)-(*p)(center-jp))/4;
           }
         else
           {
             dp_dy=((*p)(center+jp)-(*p)(center-jp))/8;
           }

         (*p_fine)(fine)=(*p)(center)
           + ((i%2==0) ? (-dp_dx) : dp_dx)
           + ((j%2==0) ? (-dp_dy) : dp_dy);

         // tbox::plog << "P_Refine "
         //            << fine_patch.getPatchLevelNumber() << " "
         //            << i << " "
         //            << j << " "
         //            << (*p_fine)(fine) << " "
         //            << (*p)(center) << " "
         //            << dp_dx << " "
         //            << dp_dy << " "
         //            << "\n";

       }
}
#endif
