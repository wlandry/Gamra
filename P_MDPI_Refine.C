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
   tbox::Pointer<pdat::NodeData<double> >
     edge_viscosity_ptr = fine_patch.getPatchData(edge_viscosity_id);
   pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!p_ptr.isNull());
   TBOX_ASSERT(!p_fine_ptr.isNull());
   TBOX_ASSERT(p.getDepth() == p_fine.getDepth());
#endif

   hier::Box coarse_box=coarse_patch.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry>
     geom = coarse_patch.getPatchGeometry();

   hier::Box gbox=p_fine.getGhostBox();

   const double dx = geom->getDx()[0];
   const double dy = geom->getDx()[1];

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

         if(fine[0]>gbox.lower(0) && fine[0]<gbox.upper(0)
            && fine[1]>gbox.lower(1) && fine[1]<gbox.upper(1))
           {
             double dRc_dp_total(0), dRc_dp_fine(0);
             /* This is horribly inefficient */
             for(int xx=0;xx<2;++xx)
               for(int yy=0;yy<2;++yy)
                 {
                   pdat::CellIndex c_fine(center*2);
                   c_fine[0]+=xx;
                   c_fine[1]+=yy;
               
                   pdat::CellIndex left(c_fine-ip),right(c_fine+ip),
                     down(c_fine-jp), up(c_fine+jp);
                   pdat::SideIndex left_x(left,0,pdat::SideIndex::Lower),
                     right_x(right,0,pdat::SideIndex::Lower),
                     down_y(down,1,pdat::SideIndex::Lower),
                     up_y(up,1,pdat::SideIndex::Lower);

                   double dRc_dp_weight=
                     dRc_dp(fine_box,c_fine,left,right,down,up,
                            left_x,right_x,down_y,up_y,
                            cell_viscosity,edge_viscosity,v,dx,dy);

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
