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

#include "V_Boundary_Refine.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"


/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void SAMRAI::geom::V_Boundary_Refine::Update_V
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const pdat::SideIndex &fine,
 const hier::Index &ip, const hier::Index &jp,
 int &j,
 tbox::Pointer<SAMRAI::pdat::SideData<double> > &v,
 tbox::Pointer<SAMRAI::pdat::SideData<double> > &v_fine) const
{
  pdat::SideIndex center(fine);
  center/=2;
  const int off_axis= (axis==0) ? 1 : 0;
  if(boundary_direction==axis)
    {
      /* Interpolate in the y direction */
      const double dv_plus=(*v)(center+jp)-(*v)(center);
      const double dv_minus=(*v)(center)-(*v)(center-jp);

      /* Quadratic interpolation */
      double v_plus=(*v)(center)
        + (5.0/32)*dv_plus - (3.0/32)*dv_minus;
      double v_minus=(*v)(center)
        - (5.0/32)*dv_minus + (3.0/32)*dv_plus;

      (*v_fine)(fine)=v_minus*(2*(*v)(center))/(v_plus + v_minus);
      (*v_fine)(fine+jp)=v_plus*(2*(*v)(center))/(v_plus + v_minus);

      /* Set the point outside of the boundary to be max_double.  This
         give us a marker for whether the current value is a boundary
         condition. */
      hier::Index offset(ip);
      if(!boundary_positive)
        {
          offset[axis]=-1;
        }
      (*v_fine)(fine+offset)=(*v_fine)(fine+offset+jp)=
        std::numeric_limits<double>::max();
      ++j;
    }
  else
    {
      double dv;
      /* Compute derivatives and use that to set the new vx */
      if(fine[axis]%2==0)
        {
          dv=((*v)(center+jp) - (*v)(center))/2;
        }
      else
        {
          dv=((*v)(center+jp) - (*v)(center)
              + (*v)(center+jp+ip) - (*v)(center+ip))/4;
        }
      hier::Index offset(jp);
      if(!boundary_positive)
        {
          offset[off_axis]=-1;
        }
      (*v_fine)(fine)=(*v_fine)(fine-offset)+dv;
    }
}
