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

#include "V_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"

#include "FTensor.hpp"

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

void SAMRAI::geom::V_Refine::refine(hier::Patch& fine_patch,
                                    const hier::Patch& coarse_patch,
                                    const int dst_component,
                                    const int src_component,
                                    const hier::Box& fine_box,
                                    const hier::IntVector& ratio,
                                    const int &axis) const
{
   const tbox::Dimension& dimension(getDim());
   const int dim(dimension.getValue());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, fine_patch, coarse_patch, fine_box, ratio);

   tbox::Pointer<pdat::SideData<double> >
   v_ptr = coarse_patch.getPatchData(src_component);
   pdat::SideData<double> &v(*v_ptr);
   tbox::Pointer<pdat::SideData<double> >
   v_fine_ptr = fine_patch.getPatchData(dst_component);
   pdat::SideData<double> &v_fine(*v_fine_ptr);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!v_ptr.isNull());
   TBOX_ASSERT(!v_fine_ptr.isNull());
   TBOX_ASSERT(v.getDepth() == v_fine.getDepth());
   TBOX_ASSERT(v.getDepth() == 1);
#endif

   hier::Box coarse_box=coarse_patch.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry>
     coarse_geom = coarse_patch.getPatchGeometry();

   tbox::Pointer<geom::CartesianPatchGeometry>
     fine_geom = fine_patch.getPatchGeometry();
   const double *dx=fine_geom->getDx();

   const hier::Box &fine_patch_box(fine_patch.getBox());

   hier::Index ip(hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
   ip[0]=1;
   jp[1]=1;
   if(dim>2)
     kp[2]=1;
   hier::Index pp[]={ip,jp,kp};

   for(pdat::CellIterator ci(fine_box); ci; ci++)
     {
       pdat::SideIndex fine(*ci,axis,pdat::SideIndex::Lower);

       pdat::SideIndex coarse(fine);
       coarse.coarsen(hier::Index::getOneIndex(dimension)*2);

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
