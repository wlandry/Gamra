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
#include "quad_offset_interpolate.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"

#include "Boundary.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void SAMRAI::geom::V_Boundary_Refine::Update_V
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const pdat::SideIndex &fine,
 const hier::Index &ip, const hier::Index &jp,
 int &i,
 int &j,
 const int &i_max,
 const int &j_max,
 tbox::Pointer<SAMRAI::pdat::SideData<double> > &v,
 tbox::Pointer<SAMRAI::pdat::SideData<double> > &v_fine) const
{
  pdat::SideIndex center(fine);
  center.coarsen(hier::Index(2,2));
  const int off_axis= (axis==0) ? 1 : 0;

  // tbox::plog << "Boundary Update V "
  //            << axis << " "
  //            << fine[0] << " "
  //            << fine[1] << " "
  //            << center[0] << " "
  //            << center[1] << " ";
  /* Set the derivative for the normal direction */
  if(boundary_direction==axis)
    {
      /* Return early if we are at j==j_max, because that is a corner
         that we do not care about. */
      if(j==j_max)
        return;
      /* Compute the derivative at the nearest three coarse points and
         then interpolate */

      const double dv_plus=(*v)(center+jp+ip)-(*v)(center+jp-ip);
      const double dv_minus=(*v)(center-jp+ip)-(*v)(center-jp-ip);
      const double dv_center=(*v)(center+ip)-(*v)(center-ip);

      double dv_fine_minus, dv_fine_plus;

      quad_offset_interpolate(dv_plus,dv_minus,dv_center,
                              dv_fine_plus,dv_fine_minus);

      hier::Index offset(ip);

      if(boundary_positive)
        {
          offset[axis]=2;
        }
      else
        {
          offset[axis]=-2;
          dv_fine_minus=-dv_fine_minus;
          dv_fine_plus=-dv_fine_plus;
        }
        
      if(j%2==0)
        {
          (*v_fine)(fine)=(*v_fine)(fine-offset) + dv_fine_minus/2;
          (*v_fine)(fine+jp)=(*v_fine)(fine-offset+jp) + dv_fine_plus/2;
          /* Since we update two points on j at once, we increment j
             again.  This is ok, since the box in the 'i' direction is
             defined to be only one cell wide */
          ++j;
        }
      else
        {
          (*v_fine)(fine)=(*v_fine)(fine-offset) + dv_fine_plus/2;
        }          
      // tbox::plog << offset[0] << " "
      //            << offset[1] << " "
      //            << (*v_fine)(fine) << " "
      //            << (*v_fine)(fine+jp) << " "
      //            << (*v_fine)(fine-offset) << " "
      //            << (*v_fine)(fine-offset+jp) << " "
      //            << (*v)(center) << " "
      //            << (*v)(center+jp) << " "
      //            << (*v)(center-jp) << " "
      //            // << (&(*v)(center-jp)) << " "
      //            // << v_plus << " "
      //            // << v_minus << " "
      //            << dv_fine_plus << " "
      //            << dv_fine_minus << " "
      //            << "\n";

    }
  /* Set the value for the tangential direction */
  else
    {
      double v_center, v_plus;

      v_center=
        quad_offset_interpolate((*v)(center-jp),(*v)(center),(*v)(center+jp));

      if(i%2==0)
        {
          (*v_fine)(fine)=v_center;

          // tbox::plog << (*v_fine)(fine) << " ";
          // // << (*v_fine)(fine+jp) << " "
          // // << (*v_fine)(fine-offset) << " "
          // // << (*v_fine)(fine-offset+jp) << " "
          // // << (*v)(center) << " "
          // // << (*v)(center+jp) << " "
          // // << (*v)(center-jp) << " "
          // // << dv_fine_plus << " "
          // // << dv_fine_minus << " "

          if(i<i_max)
            {
              v_plus=quad_offset_interpolate((*v)(center+ip-jp),(*v)(center+ip),
                                             (*v)(center+ip+jp));
              (*v_fine)(fine+ip)=(v_center+v_plus)/2;

              // tbox::plog << (*v_fine)(fine+ip) << " ";

              /* Since we update two points on 'i' at once, we increment 'i' again.
                 This is ok, since the box in the 'j' direction is defined to be
                 only one cell wide */
              ++i;
            }
        }
      else
        {
          v_plus=quad_offset_interpolate((*v)(center+ip-jp),(*v)(center+ip),
                                         (*v)(center+ip+jp));
          (*v_fine)(fine)=(v_center+v_plus)/2;
          // tbox::plog << (*v_fine)(fine) << " ";
        }
      // tbox::plog << "\n";
    }
}
