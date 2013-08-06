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

#include "P_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

void Stokes::P_Refine::refine(
   SAMRAI::hier::Patch& fine,
   const SAMRAI::hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const SAMRAI::hier::BoxOverlap& fine_overlap,
   const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::pdat::CellOverlap* t_overlap =
    dynamic_cast<const SAMRAI::pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap);

   const SAMRAI::hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
        b!=boxes.end(); ++b) {
      refine(fine,coarse,dst_component,src_component,*b,ratio);
   }
}

void Stokes::P_Refine::refine(
   SAMRAI::hier::Patch& fine_patch,
   const SAMRAI::hier::Patch& coarse_patch,
   const int dst_component,
   const int src_component,
   const SAMRAI::hier::Box& fine_box,
   const SAMRAI::hier::IntVector&) const
{
  const SAMRAI::tbox::Dimension& dim(fine_patch.getDim());

  boost::shared_ptr<SAMRAI::pdat::CellData<double> > p =
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
     (coarse_patch.getPatchData(src_component));
  boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_fine = 
    boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
     (fine_patch.getPatchData(dst_component));
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(p);
   TBOX_ASSERT(p_fine);
   TBOX_ASSERT(p->getDepth() == p_fine->getDepth());
#endif

   SAMRAI::hier::Box coarse_box=coarse_patch.getBox();
   boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
     boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
     (coarse_patch.getPatchGeometry());

   SAMRAI::pdat::CellIterator cend(fine_box,false);
   for(SAMRAI::pdat::CellIterator ci(fine_box,true); ci!=cend; ++ci)
     {
       const SAMRAI::pdat::CellIndex &fine(*ci);

       SAMRAI::pdat::CellIndex
         center(SAMRAI::hier::Index::coarsen(fine,SAMRAI::hier::Index::getOneIndex(dim)*2));

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
           SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dim));
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
