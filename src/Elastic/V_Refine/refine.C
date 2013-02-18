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

#include "Elastic/V_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"

#include "FTensor.hpp"

void Elastic::V_Refine::refine(SAMRAI::hier::Patch& fine,
                               const SAMRAI::hier::Patch& coarse,
                               const int dst_component,
                               const int src_component,
                               const SAMRAI::hier::BoxOverlap& fine_overlap,
                               const SAMRAI::hier::IntVector& ratio) const
{
   const SAMRAI::pdat::SideOverlap* t_overlap =
      dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   for(int axis=0; axis<getDim().getValue(); ++axis)
     {
       const SAMRAI::hier::BoxContainer&
         boxes = t_overlap->getDestinationBoxContainer(axis);
       for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin());
            b!=boxes.end(); b++)
         {
           refine(fine,coarse,dst_component,src_component,*b,ratio,axis);
         }
     }
}

void Elastic::V_Refine::refine
(SAMRAI::hier::Patch &fine_patch,
 const SAMRAI::hier::Patch &coarse_patch,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box &fine_box,
 const SAMRAI::hier::IntVector &,
 const int &axis) const
{
   const SAMRAI::tbox::Dimension &dimension(getDim());
   const int dim(dimension.getValue());

   boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
     boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
     (coarse_patch.getPatchData(src_component));
   SAMRAI::pdat::SideData<double> &v(*v_ptr);
   boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine_ptr =
     boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
     (fine_patch.getPatchData(dst_component));
   SAMRAI::pdat::SideData<double> &v_fine(*v_fine_ptr);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(v_ptr);
   TBOX_ASSERT(v_fine_ptr);
   TBOX_ASSERT(v.getDepth() == v_fine.getDepth());
   TBOX_ASSERT(v.getDepth() == 1);
#endif

   SAMRAI::hier::Box coarse_box=coarse_patch.getBox();
   boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> coarse_geom =
     boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
     (coarse_patch.getPatchGeometry());

   SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)),
     jp(ip), kp(ip);
   ip[0]=1;
   jp[1]=1;
   if(dim>2)
     kp[2]=1;
   SAMRAI::hier::Index pp[]={ip,jp,kp};

   SAMRAI::pdat::CellIterator cend(fine_box,false);
   for(SAMRAI::pdat::CellIterator ci(fine_box,true); ci!=cend; ci++)
     {
       SAMRAI::pdat::SideIndex fine(*ci,axis,SAMRAI::pdat::SideIndex::Lower);

       SAMRAI::pdat::SideIndex coarse(fine);
       coarse.coarsen(SAMRAI::hier::Index::getOneIndex(dimension)*2);

       if(fine[axis]%2==0)
         {
           v_fine(fine)=refine_along_line(v,axis,dim,pp,fine,coarse,coarse_box,
                                          *coarse_geom);
         }
       else
         {
           v_fine(fine)=(refine_along_line(v,axis,dim,pp,fine,coarse,coarse_box,
                                           *coarse_geom)
                         + refine_along_line(v,axis,dim,pp,fine,coarse+pp[axis],
                                             coarse_box,*coarse_geom))/2;
         }
     }
}
