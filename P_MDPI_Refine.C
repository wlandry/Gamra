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
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine_patch, coarse_patch,
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
   if(dim.getValue()==2)
     edge_viscosity2D_ptr = fine_patch.getPatchData(edge_viscosity_id);
   else
     edge_viscosity3D_ptr = fine_patch.getPatchData(edge_viscosity_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!p_ptr.isNull());
   TBOX_ASSERT(!p_fine_ptr.isNull());
   TBOX_ASSERT(p.getDepth() == p_fine.getDepth());
#endif

   hier::Box coarse_box=coarse_patch.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = coarse_patch.getPatchGeometry();

   hier::Box gbox=p_fine.getGhostBox();

   double Dx[dim.getValue()];
   for(int d=0;d<dim.getValue();++d)
     Dx[d]=geom->getDx()[d];

   hier::Box gbox_interior(gbox);
   gbox_interior.grow(hier::Index::getOneIndex(dim)*(-1));

   hier::Box cell_box(hier::Index::getZeroIndex(dim),
                      hier::Index::getOneIndex(dim)*2);

   for(pdat::CellIterator ci(fine_box); ci; ci++)
     {
       pdat::CellIndex fine(*ci);
       pdat::CellIndex
         center(hier::Index::coarsen(fine,hier::Index::getOneIndex(dim)*2));

       if(gbox_interior.contains(fine))
         {
           double dRc_dp_total(0), dRc_dp_fine(0);
           /* This is horribly inefficient */

           for(pdat::CellIterator ii(cell_box); ii; ii++)
               {
                 pdat::CellIndex c_fine(center*2);
                 c_fine+=ii;
               
                 double dRc_dp_weight;
                 if(dim.getValue()==2)
                   {
                     pdat::SideIndex x(c_fine,0,pdat::SideIndex::Lower),
                       y(c_fine,1,pdat::SideIndex::Lower);
                     dRc_dp_weight=dRc_dp_2D(fine_box,c_fine,x,y,
                                             cell_viscosity,
                                             *edge_viscosity2D_ptr,v,Dx[0],Dx[1]);
                   }
                 else
                   {
                     const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
                     const hier::Index pp[]={ip,jp,kp};
                     dRc_dp_weight=dRc_dp_3D(fine_box,c_fine,cell_viscosity,
                                             *edge_viscosity3D_ptr,v,Dx,pp);
                   }

                 if(c_fine==fine)
                   dRc_dp_fine=dRc_dp_weight;
                 dRc_dp_total+=dRc_dp_weight;
               }

           p_fine(fine)=p(center)*dRc_dp_total/(4*dRc_dp_fine);
         }
       else
         {
           p_fine(fine)=p(center);
         }

       // tbox::plog << "P_MDPI_Refine "
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
