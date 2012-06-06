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

#include "FTensor.hpp"

using namespace SAMRAI;

inline void coarsen_point_2D(const pdat::SideIndex &coarse,
                             const hier::Index &ip, const hier::Index &jp,
                             tbox::Pointer<pdat::SideData<double> > &v,
                             tbox::Pointer<pdat::SideData<double> > &v_fine )
{
  pdat::SideIndex center(coarse*2);
  (*v)(coarse)=((*v_fine)(center) + (*v_fine)(center+jp))/4
    + ((*v_fine)(center-ip) + (*v_fine)(center-ip+jp)
       + (*v_fine)(center+ip) + (*v_fine)(center+jp+ip))/8;
}


void SAMRAI::geom::V_Coarsen::coarsen_2D(hier::Patch& coarse,
                                         const hier::Patch& fine,
                                         const int dst_component,
                                         const int src_component,
                                         const hier::Box& coarse_box,
                                         const hier::IntVector& ratio) const
{
  const tbox::Dimension& Dim(getDim());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(Dim, coarse, fine, coarse_box, ratio);

  tbox::Pointer<pdat::SideData<double> >
    v_fine = fine.getPatchData(src_component);
  tbox::Pointer<pdat::SideData<double> >
    v = coarse.getPatchData(dst_component);

  TBOX_ASSERT(!v.isNull());
  TBOX_ASSERT(!v_fine.isNull());
  TBOX_ASSERT(v_fine->getDepth() == v->getDepth());
  TBOX_ASSERT(v->getDepth() == 1);

  const hier::IntVector& directions(v->getDirectionVector());

  TBOX_ASSERT(directions ==
              hier::IntVector::min(directions, v_fine->getDirectionVector()));

  const tbox::Pointer<CartesianPatchGeometry> cgeom =
    coarse.getPatchGeometry();
  const double *dx=cgeom->getDx();

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
  hier::Index ip(1,0), jp(0,1);
  for(int j=coarse_box.lower(1); j<=coarse_box.upper(1)+1; ++j)
    for(int i=coarse_box.lower(0); i<=coarse_box.upper(0)+1; ++i)
      {
        if(directions(0) && j!=coarse_box.upper(1)+1)
          {
            pdat::SideIndex coarse(hier::Index(i,j),0,
                                   pdat::SideIndex::Lower);
            pdat::SideIndex fine(coarse*2);
            if((i==coarse_box.lower(0)
                && cgeom->getTouchesRegularBoundary(0,0))
               || (i==coarse_box.upper(0)+1
                   && cgeom->getTouchesRegularBoundary(0,1)))
              {
                (*v)(coarse)=((*v_fine)(fine) + (*v_fine)(fine+jp))/2;
              }
            else
              {
                const int axis=0;
                FTensor::Tensor1<double,3> offset(0,0,0);
                offset(axis)=dx[axis]/2;
                FTensor::Tensor1<double,3> xyz(0,0,0);
                for(int d=0;d<Dim.getValue();++d)
                  xyz(d)=cgeom->getXLower()[d]
                    + dx[d]*(coarse[d]-coarse_box.lower()[d] + 0.5) - offset(d);

                coarsen_point_2D(coarse,ip,jp,v,v_fine);


                // tbox::pout << "coarsen "
                //            << coarse << " "
                //            << xyz(0) << " "
                //            << xyz(1) << " "
                //            << dx[0] << " "
                //            << fine << " "
                //            << (*v)(coarse) << " "
                //            << (*v_fine)(fine) << " "
                //            << (*v_fine)(fine+ip) << " "
                //            << (*v_fine)(fine-ip) << " "
                //            << (*v_fine)(fine+jp) << " "
                //            << (*v_fine)(fine+jp+ip) << " "
                //            << (*v_fine)(fine+jp-ip) << " "
                //            << "\n";
              }
          }
        if(directions(1) && i!=coarse_box.upper(0)+1)
          {
            pdat::SideIndex coarse(hier::Index(i,j),1,
                                   pdat::SideIndex::Lower);
            pdat::SideIndex fine(coarse*2);
            if((j==coarse_box.lower(1)
                && cgeom->getTouchesRegularBoundary(1,0))
               || (j==coarse_box.upper(1)+1
                   && cgeom->getTouchesRegularBoundary(1,1)))
              {
                (*v)(coarse)=((*v_fine)(fine) + (*v_fine)(fine+ip))/2;
              }
            else
              {
                const int axis=1;
                FTensor::Tensor1<double,3> offset(0,0,0);
                offset(axis)=dx[axis]/2;
                FTensor::Tensor1<double,3> xyz(0,0,0);
                for(int d=0;d<Dim.getValue();++d)
                  xyz(d)=cgeom->getXLower()[d]
                    + dx[d]*(coarse[d]-coarse_box.lower()[d] + 0.5) - offset(d);

                if(xyz(0)-dx[0]<0.5 && xyz(0)+dx[0]>0.5)
                  {
                    /* Interface between coarse and fine+1 */
                    if((xyz(1)+dx[1]/2>0.6 && xyz(1)<0.6)
                       || (xyz(1)+dx[1]/2>0.4 && xyz(1)<0.4))
                      {
                        // tbox::pout << "coarsen m ";
                        (*v)(coarse)=(((*v_fine)(fine) + (*v_fine)(fine+ip))*2
                                      + (*v_fine)(fine-jp) + (*v_fine)(fine+ip-jp))/3;
                      }
                    /* Interface between coarse and fine-1 */
                    else if((xyz(1)-dx[1]/2<0.6 && xyz(1)>0.6)
                            || (xyz(1)-dx[1]/2<0.4 && xyz(1)>0.4))
                      {
                        // tbox::pout << "coarsen p ";
                        (*v)(coarse)=(((*v_fine)(fine) + (*v_fine)(fine+ip))*2
                                      + (*v_fine)(fine+jp) + (*v_fine)(fine+ip+jp))/3;
                      }
                    else
                      {
                        // tbox::pout << "coarsen z ";
                        coarsen_point_2D(coarse,jp,ip,v,v_fine);
                      }
                  }
                else
                  {
                    // tbox::pout << "coarsen   ";
                    coarsen_point_2D(coarse,jp,ip,v,v_fine);
                  }

                // /* Interface between coarse and fine-1 */
                // if(xyz(0)+dx[0]/4>0.5 && xyz(0)<0.5)
                //   {
                //     tbox::pout << "coarsen m ";
                //     (*v)(coarse)=(*v_fine)(fine)/2
                //       + ((*v_fine)(fine-jp) + (*v_fine)(fine+jp))/4;
                //   }
                // /* Interface between coarse and fine+1 */
                // else if(xyz(0)-dx[0]/4<0.5 && xyz(0)>0.5)
                //   {
                //     tbox::pout << "coarsen p ";
                //     (*v)(coarse)=(*v_fine)(fine+ip)/2
                //       + ((*v_fine)(fine-jp+ip) + (*v_fine)(fine+jp+ip))/4;
                //   }
                // else
                //   {
                //     tbox::pout << "coarsen   ";
                //     coarsen_point_2D(coarse,jp,ip,v,v_fine);
                //   }

                // tbox::pout << coarse << " "
                //            << xyz(0) << " "
                //            << xyz(1) << " "
                //            << dx[0] << " "
                //            << fine << " "
                //            << (*v)(coarse) << " "
                //            << (*v_fine)(fine) << " "
                //            << (*v_fine)(fine+jp) << " "
                //            << (*v_fine)(fine-jp) << " "
                //            << (*v_fine)(fine+ip) << " "
                //            << (*v_fine)(fine+ip+jp) << " "
                //            << (*v_fine)(fine+ip-jp) << " "
                //            << "\n";
              }
          }
      }
}

#endif
