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

#include "P_Boundary_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"

void SAMRAI::geom::P_Boundary_Refine::refine(hier::Patch& fine,
                                             const hier::Patch& coarse,
                                             const int dst_component,
                                             const int src_component,
                                             const hier::BoxOverlap& overlap,
                                             const hier::IntVector& ratio)
  const
{
  const pdat::CellOverlap* t_overlap =
    dynamic_cast<const pdat::CellOverlap *>(&overlap);
  const hier::BoxList& boxes = t_overlap->getDestinationBoxList();
  for (hier::BoxList::Iterator b(boxes); b; b++)
    {
      hier::Box &overlap_box=b();

      const tbox::Dimension& dimension(getDim());
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, fine, coarse, overlap_box, ratio);
      const int dim(dimension.getValue());

      tbox::Pointer<pdat::CellData<double> >
        p = coarse.getPatchData(src_component);
      tbox::Pointer<pdat::CellData<double> >
        p_fine = fine.getPatchData(dst_component);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(!p.isNull());
      TBOX_ASSERT(!p_fine.isNull());
      TBOX_ASSERT(p->getDepth() == p_fine->getDepth());
      TBOX_ASSERT(p->getDepth() == 1);
#endif

      hier::Box coarse_box=coarse.getBox();
      hier::Box fine_box=fine.getBox();
      hier::Box gbox=p_fine->getGhostBox();
      hier::Box coarse_gbox=p->getGhostBox();

      /* We have to infer where the boundary is from the boxes */
      int boundary_direction;
      bool boundary_positive;

      for(int d=0;d<dim;++d)
        {
          if(overlap_box.lower(d)-overlap_box.upper(d)==0)
            {
              boundary_direction=d;
              if(fine_box.upper(d)<overlap_box.lower(d))
                boundary_positive=true;
              else if(fine_box.lower(d)>overlap_box.upper(d))
                boundary_positive=false;
              else
                abort();
              break;
            }
        }

      hier::Index p_min(dimension), p_max(dimension);
      for(int d=0;d<dim;++d)
        {
          p_min[d]=std::max(overlap_box.lower(d),gbox.lower(d));
          p_max[d]=std::min(overlap_box.upper(d),gbox.upper(d));
        }
      hier::Box box(p_min,p_max);

      hier::Index ip(hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
      ip[0]=1;
      jp[1]=1;
      if(dim>2)
        kp[2]=1;

      if(dim==2)
        {
          for(int j=p_min[1]; j<=p_max[1]; ++j)
            for(int i=p_min[0]; i<=p_max[0]; ++i)
              {
                pdat::CellIndex fine(hier::Index(i,j));
        
                switch(boundary_direction)
                  {
                  case 0:
                    if(j<p_max[1])
                      Update_P_2D(fine,boundary_positive ? ip : -ip,jp,j,
                                  p,p_fine);
                    break;
                  case 1:
                    if(i<p_max[0])
                      Update_P_2D(fine,boundary_positive ? jp : -jp,ip,i,
                                  p,p_fine);
                    break;
                  }
              }
        }
      else
        {
          for(int k=p_min[2]; k<=p_max[2]; ++k)
            for(int j=p_min[1]; j<=p_max[1]; ++j)
              for(int i=p_min[0]; i<=p_max[0]; ++i)
                {
                  pdat::CellIndex fine(hier::Index(i,j));
        
                  switch(boundary_direction)
                    {
                      // case 0:
                      //   if(j<p_max[1] && k<p_max[2])
                      //     Update_P_3D(fine,boundary_positive ? ip : -ip,jp,kp,j,
                      //                 p,p_fine);
                      //   break;
                      // case 1:
                      //   if(i<p_max[0] && k<p_max[2])
                      //     Update_P_3D(fine,boundary_positive ? jp : -jp,kp,ip,k,
                      //                 p,p_fine);
                      //   break;
                      // case 1:
                      //   if(i<p_max[0] && j<p_max[1])
                      //     Update_P_3D(fine,boundary_positive ? kp : -kp,ip,jp,i,
                      //                 p,p_fine);
                      //   break;
                    }
                }
        }
    }
}
