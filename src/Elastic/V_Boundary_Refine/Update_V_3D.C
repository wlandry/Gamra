#include "Elastic/V_Boundary_Refine.h"
#include "Constants.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void SAMRAI::geom::V_Boundary_Refine::Update_V_3D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const pdat::SideIndex &fine,
 const hier::Index pp[],
 const hier::Index &ijk,
 const hier::Index &p_min, const hier::Index &p_max,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_fine) const
{
  /* Set the derivative for the normal direction.

     We set the derivative on the i=constant plane.  If we look at a
     slice in that plane.

         k-1      k       k+1

        ------- -------
       |       |       |
   j-1 |   D   |   D   |   D
       |       |       |
        ------- -------
       |       |d-- d+-|
   j   |   D   |   D   |   D
       |       |d-+ d++|
        ------- -------
       |       |       |
   j+1 |   D   |   D   |   D
       |       |       |
        ------- -------
               |
               |
               |
        Coarse-Fine Boundary

  where D are the coarse derivatives, and d are the fine derivatives.
  This picture is the same as what is seen in P_Boundary_Refine in 2D.
  So we can use the formula there to compute d--.

   d(-,-) = a - b/4 + c/16 - d/4 + e/16 + f/16
          = D(-,-)/16 + (15/16)*D(0,0)
            + (3/32)*(-D(+,0) - D(0,+) + D(-,0) + D(0,-))

  */

  if(boundary_direction==axis)
    {
      /* Return early if we are at j==j_max, because that is a corner
         that we do not care about.  We also skip if j==j_min as long
         as we do not have to do j_min+1. We have to skip these even
         though they are not used because otherwise we could end up
         reading past the end of the array.  */
      const int axis2((axis+1)%3), axis3((axis+2)%3);
      if(ijk[axis2]==p_max[axis2] || (ijk[axis2]==p_min[axis2] && ijk[axis2]%2!=0)
         || ijk[axis3]==p_max[axis3] 
         || (ijk[axis3]==p_min[axis3] && ijk[axis3]%2!=0))
        return;
      /* Compute the derivative at all of the interpolation points.  */

      const hier::Index ip(boundary_positive ? pp[axis] : -pp[axis]),
        jp(pp[(axis+1)%3]), kp(pp[(axis+2)%3]);

      pdat::SideIndex center(fine-ip);
      center.coarsen(hier::Index(2,2,2));

      const double dv_mm=v(center-jp-kp+ip) - v(center-jp-kp-ip);
      const double dv_m0=v(center-jp+ip) - v(center-jp-ip);
      const double dv_mp=v(center-jp+kp+ip) - v(center-jp+kp-ip);

      const double dv_0m=v(center-kp+ip) - v(center-kp-ip);
      const double dv_00=v(center+ip) - v(center-ip);
      const double dv_0p=v(center+kp+ip) - v(center+kp-ip);

      const double dv_pm=v(center+jp-kp+ip) - v(center+jp-kp-ip);
      const double dv_p0=v(center+jp+ip) - v(center+jp-ip);
      const double dv_pp=v(center+jp+kp+ip) - v(center+jp+kp-ip);

      const double dv_fine_mm=dv_mm/16 + (15.0/16)*dv_00
        + (3/32)*(-dv_p0 - dv_0p + dv_m0 + dv_0m);

      const double dv_fine_mp=dv_mp/16 + (15.0/16)*dv_00
        + (3/32)*(-dv_p0 - dv_0m + dv_m0 + dv_0p);

      const double dv_fine_pm=dv_pm/16 + (15.0/16)*dv_00
        + (3/32)*(-dv_m0 - dv_0p + dv_p0 + dv_0m);

      const double dv_fine_pp=dv_pp/16 + (15.0/16)*dv_00
        + (3/32)*(-dv_m0 - dv_0m + dv_p0 + dv_0p);

      hier::Index offset(ip*2);

      /* Be careful about using the right interpolation if the fine
       * points are not aligned with the coarse points. */
      if(ijk[axis2]%2==0)
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=v_fine(fine-offset) + dv_fine_mm/2;
              if(ijk[axis2]<p_max[axis2])
                v_fine(fine+jp)=v_fine(fine-offset+jp) + dv_fine_pm/2;
              if(ijk[axis3]<p_max[axis3])
                v_fine(fine+kp)=v_fine(fine-offset+kp) + dv_fine_mp/2;
              if(ijk[axis2]<p_max[axis2] && ijk[axis3]<p_max[axis3])
                v_fine(fine+jp+kp)=v_fine(fine-offset+jp+kp) + dv_fine_pp/2;
            }
          else
            {
              v_fine(fine)=v_fine(fine-offset) + dv_fine_mp/2;
              if(ijk[axis2]<p_max[axis2])
                v_fine(fine+jp)=v_fine(fine-offset+jp) + dv_fine_pp/2;
            }
        }
      else
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=v_fine(fine-offset) + dv_fine_pm/2;
              if(ijk[axis3]<p_max[axis3])
                v_fine(fine+kp)=v_fine(fine-offset+kp) + dv_fine_pp/2;
            }
          else
            {
              v_fine(fine)=v_fine(fine-offset) + dv_fine_pp/2;
            }
        }          
    }
  /* Set the value for the tangential direction.

     Again, if we look at a slice in the i=constant plane.

          j-1      j      j+1

        ------- -------
       |       |       |
   k-1 |   V   |   V   |   V
       |       |       |
        ------- -------
       |       |v--    |
   k   |   V   |   V   |   V
       |       |v-+    |
        ------- -------
       |       |       |
   k+1 |   V   |   V   |   V
       |       |       |
        ------- -------
               |
               |
               |
        Coarse-Fine Boundary

  where V are the coarse velocities, and v are the fine velocities.
  This picture is the same as what is seen in P_Boundary_Refine in 2D.
  So we can use the formula there to compute v--.

   v(-,-) = V(-,-)/16 + (15/16)*V(0,0)
            + (3/32)*(-V(+,0) - V(0,+) + V(-,0) + V(0,-))

 */
  else
    {
      const int axis3((axis+1)%3 != boundary_direction ? (axis+1)%3 : (axis+2)%3);
      const hier::Index ip(pp[axis]),
        jp(boundary_positive ? pp[boundary_direction] : -pp[boundary_direction]),
        kp(pp[axis3]);

      pdat::SideIndex center(fine);
      center.coarsen(hier::Index(2,2,2));

      double v_minus=v(center-jp-kp)/16 + (15.0/16)*v(center)
        + (3.0/32)*(-v(center+jp) - v(center+kp) + v(center-jp) + v(center-kp));
      
      double v_plus=v(center-jp+kp)/16 + (15.0/16)*v(center)
        + (3.0/32)*(-v(center+jp) - v(center-kp) + v(center-jp) + v(center+kp));


      /* Be careful about using the right interpolation if the fine
       * points are not aligned with the coarse points. */
      if(ijk[axis]%2==0)
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=v_minus;
              if(ijk[axis3]<p_max[axis3])
                v_fine(fine+kp)=v_plus;
              if(ijk[axis]<p_max[axis])
                {
                  double v_minus_off=v(center-jp-kp+ip)/16
                    + (15.0/16)*v(center+ip)
                    + (3.0/32)*(-v(center+jp+ip) - v(center+kp+ip)
                                + v(center-jp+ip) + v(center-kp+ip));
      
                  double v_plus_off=v(center-jp+kp+ip)/16
                    + (15.0/16)*v(center+ip)
                    + (3.0/32)*(-v(center+jp+ip) - v(center-kp+ip)
                                + v(center-jp+ip) + v(center+kp+ip));

                  v_fine(fine+ip)=(v_minus+v_minus_off)/2;
                  if(ijk[axis3]<p_max[axis3])
                    v_fine(fine+ip+kp)=(v_plus+v_plus_off)/2;
                }
            }
          else
            {
              v_fine(fine)=v_plus;
              if(ijk[axis]<p_max[axis])
                {
                  double v_plus_off=v(center-jp+kp+ip)/16
                    + (15.0/16)*v(center+ip)
                    + (3.0/32)*(-v(center+jp+ip) - v(center-kp+ip)
                                + v(center-jp+ip) + v(center+kp+ip));

                  v_fine(fine+ip)=(v_plus+v_plus_off)/2;
                }
            }
        }
      else
        {
          if(ijk[axis3]%2==0)
            {
              double v_minus_off=v(center-jp-kp+ip)/16
                + (15.0/16)*v(center+ip)
                + (3.0/32)*(-v(center+jp+ip) - v(center+kp+ip)
                            + v(center-jp+ip) + v(center-kp+ip));
      
              double v_plus_off=v(center-jp+kp+ip)/16
                + (15.0/16)*v(center+ip)
                + (3.0/32)*(-v(center+jp+ip) - v(center-kp+ip)
                            + v(center-jp+ip) + v(center+kp+ip));

              v_fine(fine)=(v_minus+v_minus_off)/2;
              if(ijk[axis3]<p_max[axis3])
                v_fine(fine+kp)=(v_plus+v_plus_off)/2;
            }
          else
            {
              double v_plus_off=v(center-jp+kp+ip)/16
                + (15.0/16)*v(center+ip)
                + (3.0/32)*(-v(center+jp+ip) - v(center-kp+ip)
                            + v(center-jp+ip) + v(center+kp+ip));

              v_fine(fine)=(v_plus+v_plus_off)/2;
            }
        }
    }
}
