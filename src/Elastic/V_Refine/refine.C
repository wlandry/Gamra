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
       const SAMRAI::hier::BoxList&
         boxes = t_overlap->getDestinationBoxList(axis);
       for (SAMRAI::hier::BoxList::Iterator b(boxes); b; b++)
         {
           refine(fine,coarse,dst_component,src_component,b(),ratio,axis);
         }
     }
}

void Elastic::V_Refine::refine
(SAMRAI::hier::Patch& fine_patch,
 const SAMRAI::hier::Patch& coarse_patch,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box& fine_box,
 const SAMRAI::hier::IntVector& ratio,
 const int &axis) const
{
   const SAMRAI::tbox::Dimension& dimension(getDim());
   const int dim(dimension.getValue());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, fine_patch, coarse_patch,
                                   fine_box, ratio);

   SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
   v_ptr = coarse_patch.getPatchData(src_component);
   SAMRAI::pdat::SideData<double> &v(*v_ptr);
   SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
   v_fine_ptr = fine_patch.getPatchData(dst_component);
   SAMRAI::pdat::SideData<double> &v_fine(*v_fine_ptr);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!v_ptr.isNull());
   TBOX_ASSERT(!v_fine_ptr.isNull());
   TBOX_ASSERT(v.getDepth() == v_fine.getDepth());
   TBOX_ASSERT(v.getDepth() == 1);
#endif

   SAMRAI::hier::Box coarse_box=coarse_patch.getBox();
   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
     coarse_geom = coarse_patch.getPatchGeometry();

   SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
     fine_geom = fine_patch.getPatchGeometry();
   const double *dx=fine_geom->getDx();

   const SAMRAI::hier::Box &fine_patch_box(fine_patch.getBox());

   SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)),
     jp(ip), kp(ip);
   ip[0]=1;
   jp[1]=1;
   if(dim>2)
     kp[2]=1;
   SAMRAI::hier::Index pp[]={ip,jp,kp};

   for(SAMRAI::pdat::CellIterator ci(fine_box); ci; ci++)
     {
       SAMRAI::pdat::SideIndex fine(*ci,axis,SAMRAI::pdat::SideIndex::Lower);

       SAMRAI::pdat::SideIndex coarse(fine);
       coarse.coarsen(SAMRAI::hier::Index::getOneIndex(dimension)*2);

       FTensor::Tensor1<double,3> offset(0,0,0);
       offset(axis)=dx[axis]/2;
       FTensor::Tensor1<double,3> xyz(0,0,0);
       for(int d=0;d<dim;++d)
         xyz(d)=fine_geom->getXLower()[d]
           + dx[d]*(fine[d]-fine_patch_box.lower()[d] + 0.5) - offset(d);

       if(fine[axis]%2==0)
         {
           v_fine(fine)=
             refine_along_line(v,axis,dim,pp,fine,coarse,coarse_box,
                               *coarse_geom,xyz,dx);
         }
       else
         {
           FTensor::Tensor1<double,3> xyz_low, xyz_high;
           FTensor::Index<'a',3> a;

           xyz_low(a)=xyz(a) - 2*offset(a);
           xyz_high(a)=xyz(a) + 2*offset(a);

           v_fine(fine)=
             (refine_along_line(v,axis,dim,pp,fine,coarse,coarse_box,
                                *coarse_geom,xyz_low,dx)
              + refine_along_line(v,axis,dim,pp,fine,coarse+pp[axis],
                                  coarse_box,*coarse_geom,xyz_high,dx))/2;
         }
     }
}
