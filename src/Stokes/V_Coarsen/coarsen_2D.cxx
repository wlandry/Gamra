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

#include "Stokes/V_Coarsen.hxx"

#include <float.h>
#include <math.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Index.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/tbox/Utilities.h>
#include "Constants.hxx"

inline void coarsen_point_2D(const SAMRAI::pdat::SideIndex &coarse,
                             const SAMRAI::hier::Index &ip, const SAMRAI::hier::Index &jp,
                             boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v,
                             boost::shared_ptr<SAMRAI::pdat::SideData<double> > &v_fine )
{
  SAMRAI::pdat::SideIndex center(coarse*2);
  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/4
    + ((*v_fine)(center-ip) + (*v_fine)(center-ip+jp)
       + (*v_fine)(center+ip) + (*v_fine)(center+jp+ip))/8;
}


void Stokes::V_Coarsen::coarsen_2D
(SAMRAI::hier::Patch& coarse,
 const SAMRAI::hier::Patch& fine,
 const int dst_component,
 const int src_component,
 const SAMRAI::hier::Box& coarse_box,
 const SAMRAI::hier::IntVector&) const
{
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

     The coarse/fine boundaries get fixed up later in
     V_Coarsen_Patch_Strategy::postprocessCoarsen.
  */
  SAMRAI::hier::Index ip(1,0), jp(0,1);
  for(int j=coarse_box.lower(1); j<=coarse_box.upper(1)+1; ++j)
    for(int i=coarse_box.lower(0); i<=coarse_box.upper(0)+1; ++i)
      {
        if(directions(0) && j!=coarse_box.upper(1)+1)
          {
            SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),0,
                                           SAMRAI::pdat::SideIndex::Lower);
            SAMRAI::pdat::SideIndex center(coarse*2);
            if((i==coarse_box.lower(0)
                && cgeom->getTouchesRegularBoundary(0,0))
               || (i==coarse_box.upper(0)+1
                   && cgeom->getTouchesRegularBoundary(0,1)))
              {
                (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/2;
              }
            else
              {
                coarsen_point_2D(coarse,ip,jp,v,v_fine);
              }
          }
        if(directions(1) && i!=coarse_box.upper(0)+1)
          {
            SAMRAI::pdat::SideIndex coarse(SAMRAI::hier::Index(i,j),1,
                                           SAMRAI::pdat::SideIndex::Lower);
            SAMRAI::pdat::SideIndex center(coarse*2);
            if((j==coarse_box.lower(1)
                && cgeom->getTouchesRegularBoundary(1,0))
               || (j==coarse_box.upper(1)+1
                   && cgeom->getTouchesRegularBoundary(1,1)))
              {
                (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+ip))/2;
              }
            else
              {
                coarsen_point_2D(coarse,jp,ip,v,v_fine);
              }
          }
      }
}
