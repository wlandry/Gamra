/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#include "Elastic/V_Coarsen_Patch_Strategy.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "Constants.h"

#include "FTensor.hpp"

inline void
coarsen_point_2D(const SAMRAI::pdat::SideIndex& coarse,
                 const SAMRAI::pdat::SideIndex& fine,
                 const SAMRAI::hier::Index& ip, const SAMRAI::hier::Index& jp,
                 SAMRAI::pdat::SideData<double>& v,
                 const SAMRAI::pdat::SideData<double>& v_fine )
{
  v(coarse)=(v_fine(fine) + v_fine(fine+jp))/4
    + (v_fine(fine-ip) + v_fine(fine-ip+jp)
       + v_fine(fine+ip) + v_fine(fine+jp+ip))/8;
}

inline double
coarsen_correction_2D(const SAMRAI::pdat::SideIndex& fine,
                      const int& axis,
                      const SAMRAI::hier::Index& ip,
                      const SAMRAI::hier::Index& jp,
                      const SAMRAI::pdat::CellData<double>& dv_diagonal,
                      const SAMRAI::pdat::SideData<double>& dv_mixed)
{
  SAMRAI::pdat::CellIndex cell(fine);
  return (dv_diagonal(cell-ip,axis) - dv_diagonal(cell,axis)
          + dv_diagonal(cell-ip+jp,axis) - dv_diagonal(cell+jp,axis))/8
    + (dv_mixed(fine,0) + dv_mixed(fine+jp,1))/2;
}

void Elastic::V_Coarsen_Patch_Strategy::coarsen_2D
(SAMRAI::pdat::SideData<double>& v,
 const SAMRAI::pdat::SideData<double>& v_fine,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
 const boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal,
 const SAMRAI::geom::CartesianPatchGeometry& coarse_geom,
 const SAMRAI::hier::Box& coarse_box) const
{
  /* Numbering of v nodes is

     x--x--x--x--x  Fine
     0  1  2  3  4

     x-----x-----x Coarse
     0     1     2

     So the i'th coarse point is affected by the i*2-1,
     i*2, and i*2+1 fine points.  So, for example, i_fine=3
     affects both i_coarse=1 and i_coarse=2.

     |---------------------------------------------------------------|
     |               |               |               |               |
     f       f       f       f       f       f       f       f       f
     |               |               |               |               |
     c               c               c               c               c
     |               |               |               |               |
     f       f       f       f       f       f       f       f       f
     |               |               |               |               |
     |---------------------------------------------------------------|
     |               |               |               |               |
     f       f       f       f       f       f       f       f       f
     |               |               |               |               |
     c               c               c               c               c
     |               |               |               |               |
     f       f       f       f       f       f       f       f       f
     |               |               |               |               |
     |---------------------------------------------------------------|

     In 2D, a coarse point depends on six points.  In this
     case, (i*2,j*2), (i*2,j*2+1), (i*2-1,j*2),
     (i*2-1,j*2+1), (i*2+1,j*2), (i*2+1,j*2+1).

     The coarse/fine boundaries get fixed up later in
     V_Coarsen_Patch_Strategy::postprocessCoarsen.
  */
  SAMRAI::hier::Index ip(1,0), jp(0,1);
  const SAMRAI::hier::Index unit[]={ip,jp};
  int ijk[2];


  /* From reading CoarsenSchedule::coarsenSourceData in
     SAMRAI/source/SAMRAI/xfer/CoarsenSchedule.C, it seems that the
     coarse box is created from the fine box.  So the coarse box is
     always covered by the fine box, meaning we do not have to do an
     intersection. */
  for(ijk[1]=coarse_box.lower(1); ijk[1]<=coarse_box.upper(1)+1; ++ijk[1])
    for(ijk[0]=coarse_box.lower(0); ijk[0]<=coarse_box.upper(0)+1; ++ijk[0])
      {
        for(int axis=0;axis<2;++axis)
          {
            const int off_axis((axis+1)%2);
            if(ijk[off_axis]!=coarse_box.upper(off_axis)+1)
              {
                SAMRAI::pdat::SideIndex
                  coarse(SAMRAI::hier::Index(ijk[0],ijk[1]),axis,
                         SAMRAI::pdat::SideIndex::Lower);
                SAMRAI::pdat::SideIndex fine(coarse*2);
                if((ijk[axis]==coarse_box.lower(axis)
                    && coarse_geom.getTouchesRegularBoundary(axis,0))
                   || (ijk[axis]==coarse_box.upper(axis)+1
                       && coarse_geom.getTouchesRegularBoundary(axis,1)))
                  {
                    v(coarse)=
                      (v_fine(fine) + v_fine(fine+unit[off_axis]))/2;
                    if(have_faults() && !is_residual)
                      v(coarse)+=((*dv_mixed)(fine,0)
                                  + (*dv_mixed)(fine+unit[off_axis],1))/2;
                  }
                else
                  {
                    coarsen_point_2D(coarse,fine,unit[axis],unit[off_axis],
                                     v,v_fine);
                    if(have_faults() && !is_residual)
                      {
                        v(coarse)+=
                          coarsen_correction_2D(fine,axis,unit[axis],
                                                unit[off_axis],*dv_diagonal,
                                                *dv_mixed);
                      }
                  }
              }
          }
      }
}

