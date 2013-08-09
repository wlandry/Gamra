/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#include "V_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"

void Stokes::V_Refine::refine(SAMRAI::hier::Patch& fine,
                                            const SAMRAI::hier::Patch& coarse,
                                            const int dst_component,
                                            const int src_component,
                                            const SAMRAI::hier::BoxOverlap& fine_overlap,
                                            const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::pdat::SideOverlap* t_overlap =
    dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

  TBOX_ASSERT(t_overlap);

  for(int axis=0; axis<fine.getDim().getValue(); ++axis)
    {
      const SAMRAI::hier::BoxContainer&
        boxes = t_overlap->getDestinationBoxContainer(axis);
      SAMRAI::hier::BoxContainer::const_iterator bend(boxes.end());
      for (SAMRAI::hier::BoxContainer::const_iterator b(boxes.begin()); b!=bend; ++b)
        {
          refine(fine,coarse,dst_component,src_component,*b,ratio,axis);
        }
    }
}

void Stokes::V_Refine::refine(SAMRAI::hier::Patch& fine,
                                            const SAMRAI::hier::Patch& coarse,
                                            const int dst_component,
                                            const int src_component,
                                            const SAMRAI::hier::Box& fine_box,
                                            const SAMRAI::hier::IntVector&,
                                            const int &axis) const
{
  const SAMRAI::tbox::Dimension& dimension(fine.getDim());
  const int dim(dimension.getValue());

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (coarse.getPatchData(src_component));
  SAMRAI::pdat::SideData<double> &v(*v_ptr);
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_fine_ptr = 
    boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (fine.getPatchData(dst_component));
  SAMRAI::pdat::SideData<double> &v_fine(*v_fine_ptr);

#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(v_ptr);
  TBOX_ASSERT(v_fine_ptr);
  TBOX_ASSERT(v.getDepth() == v_fine.getDepth());
  TBOX_ASSERT(v.getDepth() == 1);
#endif

  SAMRAI::hier::Box coarse_box=coarse.getBox();
  boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (coarse.getPatchGeometry());

  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(dim>2)
    kp[2]=1;
  SAMRAI::hier::Index pp[]={ip,jp,kp};

  SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(fine_box));
  for(SAMRAI::pdat::CellIterator
        ci(SAMRAI::pdat::CellGeometry::begin(fine_box)); ci!=cend; ++ci)
    {
      SAMRAI::pdat::SideIndex fine(*ci,axis,SAMRAI::pdat::SideIndex::Lower);

      SAMRAI::pdat::SideIndex center(fine);
      center.coarsen(SAMRAI::hier::Index::getOneIndex(dimension)*2);

      /* This assumes that the levels are always properly nested, so
         that we always have an extra grid space for interpolation.
         So we only have to have a special case for physical
         boundaries, where we do not have an extra grid space. */

      double dvx_dy;

      if(fine[axis]%2==0)
        {
          /* Maybe this has to be fixed when dvx/dy != 0 on the
             outer boundary because the approximation to the
             derivative is not accurate enough? */

          /* Maybe in 3D we should include cross derivatives? */

          v_fine(fine)=v(center);

          for(int d=(axis+1)%dim;d!=axis;d=(d+1)%dim)
            {
              if(center[d]==coarse_box.lower(d)
                 && geom->getTouchesRegularBoundary(d,0))
                {
                  dvx_dy=(v(center+pp[d])-v(center))/4;
                }
              else if(center[d]==coarse_box.upper(d)
                      && geom->getTouchesRegularBoundary(d,1))
                {
                  dvx_dy=(v(center)-v(center-pp[d]))/4;
                }
              else
                {
                  dvx_dy=(v(center+pp[d])-v(center-pp[d]))/8;
                }
              v_fine(fine)+=((fine[d]%2==0) ? (-dvx_dy) : dvx_dy);
            }
        }
      else
        {
          v_fine(fine)=(v(center) + v(center+pp[axis]))/2;

          for(int d=(axis+1)%dim;d!=axis;d=(d+1)%dim)
            {
              if(center[d]==coarse_box.lower(d)
                 && geom->getTouchesRegularBoundary(d,0))
                {
                  dvx_dy=(v(center+pp[d])-v(center)
                          + v(center+pp[d]+pp[axis])-v(center+pp[axis]))/8;
                }
              else if(center[d]==coarse_box.upper(d)
                      && geom->getTouchesRegularBoundary(d,1))
                {
                  dvx_dy=(v(center)-v(center-pp[d])
                          + v(center+pp[axis])-v(center-pp[d]+pp[axis]))/8;
                }
              else
                {
                  dvx_dy=
                    (v(center+pp[d])-v(center-pp[d])
                     + v(center+pp[d]+pp[axis])-v(center-pp[d]+pp[axis]))/16;
                }
              v_fine(fine)+=((fine[d]%2==0) ? (-dvx_dy) : dvx_dy);
            }
        }
    }
}
