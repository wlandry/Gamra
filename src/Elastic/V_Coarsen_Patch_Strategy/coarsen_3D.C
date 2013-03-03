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

#ifndef included_geom_V_Coarsen_C
#define included_geom_V_Coarsen_C

#include "Elastic/V_Coarsen_Patch_Strategy.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "Constants.h"

void Elastic::V_Coarsen_Patch_Strategy::coarsen_3D
(SAMRAI::pdat::SideData<double>& v,
 const SAMRAI::pdat::SideData<double>& v_fine,
 const SAMRAI::pdat::SideData<double>& dv_mixed,
 const SAMRAI::pdat::CellData<double>& dv_diagonal,
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


          --------------------
         /                   /|
        /                   / |
       /                   /  |
      /                   /   |
      -------------------     |
     |                   |    |
     |    f        f     |    |
     |                   |    |
     |        C          |   /
     |                   |  /
     |    f        f     | /
     |                   |/
     -------------------

     In 3D, a coarse point depend on 12 points
     (i*2,j*2,k*2), (i*2,j*2+1,k*2), (i*2,j*2,k*2+1), (i*2,j*2+1,k*2+1),
     (i*2+1,j*2,k*2), (i*2+1,j*2+1,k*2), (i*2+1,j*2,k*2+1), (i*2+1,j*2+1,k*2+1),
     (i*2-1,j*2,k*2), (i*2-1,j*2+1,k*2), (i*2-1,j*2,k*2+1), (i*2-1,j*2+1,k*2+1)

     The coarse/fine boundaries get fixed up later in
     V_Coarsen_Patch_Strategy::postprocessCoarsen.
  */
  SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const SAMRAI::hier::Index unit[]={ip,jp,kp};
  const int dim(3);
  int ijk[3];
  for(ijk[2]=coarse_box.lower(2); ijk[2]<=coarse_box.upper(2)+1; ++ijk[2])
    for(ijk[1]=coarse_box.lower(1); ijk[1]<=coarse_box.upper(1)+1; ++ijk[1])
      for(ijk[0]=coarse_box.lower(0); ijk[0]<=coarse_box.upper(0)+1; ++ijk[0])
        {
          for(int ix=0;ix<dim;++ix)
            {
              const int iy((ix+1)%dim), iz((ix+2)%dim);
              if(ijk[iy]!=coarse_box.upper(iy)+1
                 && ijk[iz]!=coarse_box.upper(iz)+1)
                {
                  SAMRAI::pdat::SideIndex
                    coarse(SAMRAI::hier::Index(ijk[0],ijk[1],ijk[2]),ix,
                           SAMRAI::pdat::SideIndex::Lower);
                  SAMRAI::pdat::SideIndex fine(coarse*2);
                  if((ijk[ix]==coarse_box.lower(ix)
                      && coarse_geom.getTouchesRegularBoundary(ix,0))
                     || (ijk[ix]==coarse_box.upper(ix)+1
                         && coarse_geom.getTouchesRegularBoundary(ix,1)))
                    {
                      v(coarse)=coarsen_plane(v_fine,fine,unit[iy],unit[iz])
                        + coarsen_plane_correction(dv_mixed,fine,unit[iy],
                                                   unit[iz]);
                    }
                  else
                    {
                      v(coarse)=coarsen_point_3D(v_fine,dv_diagonal,dv_mixed,
                                                 fine,ix,unit[ix],unit[iy],
                                                 unit[iz]);
                    }
                }
            }
        }
}

#endif
