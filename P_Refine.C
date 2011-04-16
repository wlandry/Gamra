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

   for(pdat::CellIterator ci(fine_box); ci; ci++)
     {
       pdat::CellIndex fine(*ci);

       pdat::CellIndex
         center(hier::Index::coarsen(fine,hier::Index::getOneIndex(dim)*2));

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

       (*p_fine)(fine)=(*p)(center);
       for(int d=0; d<dim.getValue(); ++d)
         {
           hier::Index ip(hier::Index::getZeroIndex(dim));
           ip[d]=1;
           double dp;
           if(center[d]==coarse_box.lower(d)
              && geom->getTouchesRegularBoundary(d,0))
             {
               dp=((*p)(center+ip)-(*p)(center))/4;
             }
           else if(center[d]==coarse_box.upper(d)
                   && geom->getTouchesRegularBoundary(d,1))
             {
               dp=((*p)(center)-(*p)(center-ip))/4;
             }
           else
             {
               dp=((*p)(center+ip)-(*p)(center-ip))/8;
             }
           (*p_fine)(fine)+=((fine[d]%2==0) ? (-dp) : dp);
         }           
       }
}
#endif
