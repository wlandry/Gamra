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

#ifndef included_geom_P_MDPI_Refine_C
#define included_geom_P_MDPI_Refine_C

#include "P_MDPI_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "dRc_dp.h"

void SAMRAI::geom::P_MDPI_Refine::refine(
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

void SAMRAI::geom::P_MDPI_Refine::refine(
   hier::Patch& fine_patch,
   const hier::Patch& coarse_patch,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dimension(getDim());
   const int dim(dimension.getValue());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, fine_patch, coarse_patch,
                                   fine_box, ratio);

   tbox::Pointer<pdat::CellData<double> >
     p_ptr = coarse_patch.getPatchData(src_component);
   pdat::CellData<double> &p(*p_ptr);
   tbox::Pointer<pdat::CellData<double> >
     p_fine_ptr = fine_patch.getPatchData(dst_component);
   pdat::CellData<double> &p_fine(*p_fine_ptr);
   tbox::Pointer<pdat::SideData<double> >
     v_ptr = fine_patch.getPatchData(v_id);
   pdat::SideData<double> &v(*v_ptr);
   tbox::Pointer<pdat::CellData<double> >
     cell_viscosity_ptr = fine_patch.getPatchData(cell_viscosity_id);
   pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);
   tbox::Pointer<pdat::NodeData<double> > edge_viscosity2D_ptr;
   tbox::Pointer<pdat::EdgeData<double> > edge_viscosity3D_ptr;
   if(dim==2)
     edge_viscosity2D_ptr = fine_patch.getPatchData(edge_viscosity_id);
   else
     edge_viscosity3D_ptr = fine_patch.getPatchData(edge_viscosity_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!p_ptr.isNull());
   TBOX_ASSERT(!p_fine_ptr.isNull());
   TBOX_ASSERT(p.getDepth() == p_fine.getDepth());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = fine_patch.getPatchGeometry();

   const double *Dx=geom->getDx();

   hier::Box interior(p_fine.getBox());

   const hier::Box coarse_box
     (hier::Box::coarsen(fine_box,hier::IntVector::getOne(dimension)*2));

   hier::Box cell_box(hier::Index::getZeroIndex(dimension),
                      hier::Index::getOneIndex(dimension));

   for(pdat::CellIterator ci(coarse_box); ci; ci++)
     {
       pdat::CellIndex coarse(*ci);
       pdat::CellIndex fine(coarse*2);

       if(interior.contains(fine))
         {
           pdat::CellData<double>
             dRc_dp(cell_box,1,hier::Index::getZeroIndex(dimension));
           double dRc_dp_total(0);

           for(pdat::CellIterator ii(cell_box); ii; ii++)
             {
               pdat::CellIndex c_fine(fine+*ii);
               
               if(dim==2)
                 {
                   pdat::SideIndex x(c_fine,0,pdat::SideIndex::Lower),
                     y(c_fine,1,pdat::SideIndex::Lower);
                   dRc_dp(*ii)=dRc_dp_2D(fine_box,c_fine,x,y,cell_viscosity,
                                         *edge_viscosity2D_ptr,v,Dx[0],Dx[1]);
                 }
               else
                 {
                   const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
                   const hier::Index pp[]={ip,jp,kp};
                   dRc_dp(*ii)=dRc_dp_3D(fine_box,c_fine,cell_viscosity,
                                         *edge_viscosity3D_ptr,v,Dx,pp);
                 }
               dRc_dp_total+=dRc_dp(*ii,0);
             }

           for(pdat::CellIterator ii(cell_box); ii; ii++)
             {
               pdat::CellIndex c_fine(fine+*ii);
               p_fine(c_fine)=p(coarse)*dRc_dp_total
                 /(4*(dim-1)*dRc_dp(*ii));
             }
         }
       else
         {
           /* This should never be used as a real value, so we put in
            * a bad value so that we will notice when it happens */
           for(pdat::CellIterator ii(cell_box); ii; ii++)
             {
               pdat::CellIndex c_fine(fine+*ii);
               p_fine(c_fine)=boundary_value;
             }
         }
     }
}
#endif
