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
#include "quad_offset_interpolate.h"

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

      const tbox::Dimension& dim(getDim());
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine, coarse, overlap_box, ratio);

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
      if(overlap_box.lower(0)-overlap_box.upper(0)==0)
        {
          boundary_direction=0;
          if(fine_box.upper(0)<overlap_box.lower(0))
            boundary_positive=true;
          else if(fine_box.lower(0)>overlap_box.upper(0))
            boundary_positive=false;
          else
            abort();
        }
      else if(overlap_box.lower(1)-overlap_box.upper(1)==0)
        {
          boundary_direction=1;
          if(fine_box.upper(1)<overlap_box.lower(1))
            boundary_positive=true;
          else if(fine_box.lower(1)>overlap_box.upper(1))
            boundary_positive=false;
          else
            abort();
        }
      else
        {
          abort();
        }

      /* We have to infer where the boundary is from the boxes */
      tbox::plog << "PBR "
                 << fine.getPatchLevelNumber() << " "
                 << boundary_direction << " "
                 << std::boolalpha
                 << boundary_positive << " "
                 << coarse_box.lower(0) << " "
                 << coarse_box.upper(0) << " "
                 << coarse_box.lower(1) << " "
                 << coarse_box.upper(1) << " "
                 << fine_box.lower(0) << " "
                 << fine_box.upper(0) << " "
                 << fine_box.lower(1) << " "
                 << fine_box.upper(1) << " "

                 << overlap_box.lower(0) << " "
                 << overlap_box.upper(0) << " "
                 << overlap_box.lower(1) << " "
                 << overlap_box.upper(1) << " "
                 << "\n";

      int i_min(std::max(overlap_box.lower(0),gbox.lower(0))),
        i_max(std::min(overlap_box.upper(0),gbox.upper(0))),
        j_min(std::max(overlap_box.lower(1),gbox.lower(1))),
        j_max(std::min(overlap_box.upper(1),gbox.upper(1)));

      for(int j=j_min; j<=j_max; ++j)
        for(int i=i_min; i<=i_max; ++i)
          {
            pdat::CellIndex fine(hier::Index(i,j));
            hier::Index ip(1,0), jp(0,1);

            switch(boundary_direction)
              {
              case 0:
                if(j<j_max)
                  Update_P(fine,boundary_positive ? ip : -ip,jp,j,p,p_fine);
                break;
              case 1:
                if(i<i_max)
                  Update_P(fine,boundary_positive ? jp : -jp,ip,i,p,p_fine);
              }
          }
    }
}

/* Interpolate the pressure from the coarse (C) to the fine (f+/-)
   coordinates using the intermediate fine points (F+/-, F_).

   i-1      i       i+1

   ------- -------
   |       |       |
   j-1 |  C    |   C   |   C
   |       |       |
   ------- -------
   |       |f- F-  |
   j   |   C   |F_ C   |   C
   |       |f+ F+  |
   ------- -------
   |       |       |
   j+1 |   C   |   C   |   C
   |       |       |
   ------- -------

   F+/- = interpolate(C(i,j-1),C(i,j),C(i,j+1))
   F_   = interpolate(C(i-1,j),C(i,j),C(i+1,j))

   then

   f+/- = F+/- + F_ - C(i,j) 
   + (C(i-1,j-1) + C(i,j) + C(i-1,j) + C(i,j-1))/16

   This example show a boundary in the positive x direction.  To reverse
   the direction, pass in ip -> -ip.  To do the y direction, switch ip
   and jp, and replace j with i.

*/
     

void SAMRAI::geom::P_Boundary_Refine::Update_P
(const pdat::CellIndex &fine,
 const hier::Index &ip, const hier::Index &jp,
 int &j,
 tbox::Pointer<SAMRAI::pdat::CellData<double> > &p,
 tbox::Pointer<SAMRAI::pdat::CellData<double> > &p_fine) const
{
  double p_plus, p_minus, p_offset;
  pdat::CellIndex center(fine);
  center.coarsen(hier::Index(2,2));

  quad_offset_interpolate((*p)(center+jp),(*p)(center),(*p)(center-jp),
                          p_plus,p_minus);
  p_offset=
    quad_offset_interpolate((*p)(center-ip),(*p)(center),(*p)(center+ip));

  const double p_low=p_minus + p_offset - (*p)(center)
    + ((*p)(center-ip-jp) + (*p)(center)
       - (*p)(center-ip) - (*p)(center-jp))/16;

  const double p_high=p_plus + p_offset - (*p)(center)
    + ((*p)(center-ip+jp) + (*p)(center)
       - (*p)(center-ip) - (*p)(center+jp))/16;


  /* If we are at an even index, update both of the elements in the cell */
  if(j%2==0)
    {
      (*p_fine)(fine)=p_low;

      (*p_fine)(fine+jp)=p_high;

      /* Since we update two points on j at once, we increment j again.
         This is ok, since the box in the 'i' direction is defined to be
         only one cell wide */
      ++j;
    }
  else
    {
      (*p_fine)(fine)=p_high;
    }

  tbox::plog << "p bc "
             << fine[0] << " "
             << fine[1] << " "
             << center[0] << " "
             << center[1] << " "
             << jp[0] << " "
             << jp[1] << " "
             << (*p_fine)(fine) << " "
             << (*p_fine)(fine+jp) << " "
             << p_minus << " "
             << p_plus << " "
             << p_offset << " "
             << (*p)(center+jp) << " "
             << (*p)(center) << " "
             << (*p)(center-jp) << " "
             << (*p)(center+ip) << " "
             << (*p)(center-ip) << " "
             << (*p)(center-ip+jp) << " "
             << (*p)(center-ip-jp) << " "
             << "\n";
}
