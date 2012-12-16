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

#include "Elastic/V_Coarsen.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "Constants.h"

inline void
coarsen_point_3D(const SAMRAI::pdat::SideIndex &coarse,
                 const SAMRAI::hier::Index &ip,
                 const SAMRAI::hier::Index &jp,
                 const SAMRAI::hier::Index &kp,
                 boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v,
                 boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v_fine)
{
  SAMRAI::pdat::SideIndex center(coarse*2);
  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp)
                + (*v_fine)(center+kp) + (*v_fine)(center+jp+kp))/8
    + ((*v_fine)(center-ip) + (*v_fine)(center-ip+jp)
       + (*v_fine)(center-ip+kp) + (*v_fine)(center-ip+jp+kp)
       + (*v_fine)(center+ip) + (*v_fine)(center+jp+ip)
       + (*v_fine)(center+ip+kp) + (*v_fine)(center+jp+ip+kp))/16;
}


void Elastic::V_Coarsen::coarsen_3D(SAMRAI::hier::Patch& coarse,
                                    const SAMRAI::hier::Patch& fine,
                                    const int dst_component,
                                    const int src_component,
                                    const SAMRAI::hier::Box& coarse_box,
                                    const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::tbox::Dimension& dim(getDim());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(src_component));
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse.getPatchData(dst_component));

  TBOX_ASSERT(v);
  TBOX_ASSERT(v_fine);
  TBOX_ASSERT(v_fine->getDepth() == v->getDepth());
  TBOX_ASSERT(v->getDepth() == 1);

  const SAMRAI::hier::IntVector& directions(v->getDirectionVector());

  TBOX_ASSERT(directions ==
              SAMRAI::hier::IntVector::min(directions, v_fine->getDirectionVector()));

  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> cgeom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (coarse.getPatchGeometry());

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
  for(int k=coarse_box.lower(2); k<=coarse_box.upper(2)+1; ++k)
    for(int j=coarse_box.lower(1); j<=coarse_box.upper(1)+1; ++j)
      for(int i=coarse_box.lower(0); i<=coarse_box.upper(0)+1; ++i)
        {
          if(directions(0) && j!=coarse_box.upper(1)+1
             && k!=coarse_box.upper(2)+1)
            {
              SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j,k),0,
                                             SAMRAI::pdat::SideIndex::Lower);
              SAMRAI::pdat::SideIndex center(coarse*2);
              if((i==coarse_box.lower(0)
                  && cgeom->getTouchesRegularBoundary(0,0))
                 || (i==coarse_box.upper(0)+1
                     && cgeom->getTouchesRegularBoundary(0,1)))
                {
                  (*v)(coarse)=
                    ((*v_fine)(center) + (*v_fine)(center+jp)
                     + (*v_fine)(center+kp) + (*v_fine)(center+jp+kp))/4;
                }
              else
                {
                  coarsen_point_3D(coarse,ip,jp,kp,v,v_fine);
                }
            }
          if(directions(1) && i!=coarse_box.upper(0)+1
             && k!=coarse_box.upper(2)+1)
            {
              SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j,k),1,
                                             SAMRAI::pdat::SideIndex::Lower);
              SAMRAI::pdat::SideIndex center(coarse*2);
              if((j==coarse_box.lower(1)
                  && cgeom->getTouchesRegularBoundary(1,0))
                 || (j==coarse_box.upper(1)+1
                     && cgeom->getTouchesRegularBoundary(1,1)))
                {
                  (*v)(coarse)=
                    ((*v_fine)(center) + (*v_fine)(center+ip)
                     + (*v_fine)(center+kp) + (*v_fine)(center+ip+kp))/4;
                }
              else
                {
                  coarsen_point_3D(coarse,jp,kp,ip,v,v_fine);
                }
            }
          if(directions(2) && i!=coarse_box.upper(0)+1
             && j!=coarse_box.upper(1)+1)
            {
              SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j,k),2,
                                             SAMRAI::pdat::SideIndex::Lower);
              SAMRAI::pdat::SideIndex center(coarse*2);
              if((k==coarse_box.lower(2)
                  && cgeom->getTouchesRegularBoundary(2,0))
                 || (k==coarse_box.upper(2)+1
                     && cgeom->getTouchesRegularBoundary(2,1)))
                {
                  (*v)(coarse)=
                    ((*v_fine)(center) + (*v_fine)(center+ip)
                     + (*v_fine)(center+jp) + (*v_fine)(center+ip+jp))/4;
                }
              else
                {
                  coarsen_point_3D(coarse,kp,ip,jp,v,v_fine);
                }
            }
        }
}

#endif
