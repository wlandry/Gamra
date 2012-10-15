#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::V_Boundary_Refine::Update_V_3D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::Index pp[],
 const SAMRAI::hier::Index &ijk,
 const SAMRAI::hier::Index &p_max,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_fine) const
{
  /* Quadratic interpolation involving both coarse and fine grids for
     the normal direction

      i-1      i      i+1

        ------- -------
       |   f   f   F   |
   j-1 C       C       C
       |   f   f   F   |
        ------- -------
       |   f   f   F   |
   j   C       C       C
       |   f   f   F   |
        ------- -------
       |   f   f   F   |
   j+1 C       C       C
       |   f   f   F   |
        ------- -------
               |
               |
               |
        Coarse-Fine Boundary

      Interpolate to F.

      Note that F is offset out of the plane

           --------------------
          /                   /|
         /                   / |
        /                   /  |
       /            C      /   |
      /      F     /      /    |
      -------------------      |
     |     /     /       |     |
     |    f     /  f     |     |
     |         /         |     |
     |        C          |    /
     |                   |   /
     |                   |  /
     |    f        f     | /
     |                   |/
     -------------------


     So need to do two interpolations of coarse values on the face and
     then an interpolation using coarse and fine values to get inside
     the cube.

  */

  if(boundary_direction==axis)
    {
      const int axis2((axis+1)%3), axis3((axis+2)%3);
      const SAMRAI::hier::Index ip_s(boundary_positive ? pp[axis] : -pp[axis]),
        jp(pp[(axis+1)%3]), kp(pp[(axis+2)%3]);

      SAMRAI::pdat::SideIndex center(fine-ip_s);
      center.coarsen(SAMRAI::hier::Index(2,2,2));

      double v_pp, v_pm, v_mp, v_mm;
      double v_p[3], v_m[3];
      quad_offset_interpolate(v(center+ip_s+jp+kp),v(center+ip_s+kp),
                              v(center+ip_s-jp+kp),v_p[0],v_m[0]);
      quad_offset_interpolate(v(center+ip_s+jp),v(center+ip_s),
                              v(center+ip_s-jp),v_p[1],v_m[1]);
      quad_offset_interpolate(v(center+ip_s+jp-kp),v(center+ip_s-kp),
                              v(center+ip_s-jp-kp),v_p[2],v_m[2]);

      quad_offset_interpolate(v_p[0],v_p[1],v_p[2],v_pp,v_pm);
      quad_offset_interpolate(v_m[0],v_m[1],v_m[2],v_mp,v_mm);

      if(ijk[axis2]%2==0)
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=v_fine(fine-ip_s)
                + (v_mm - v_fine(fine-ip_s-ip_s))/3;
              // if(ijk[axis2]<p_max[axis2])
              //   v_fine(fine+jp)=v_fine(fine-ip_s+jp)
              //     + (v_pm - v_fine(fine-ip_s-ip_s+jp))/3;
              // if(ijk[axis3]<p_max[axis3])
              //   v_fine(fine+kp)=v_fine(fine-ip_s+kp)
              //     + (v_mp - v_fine(fine-ip_s-ip_s+kp))/3;
              // if(ijk[axis2]<p_max[axis2] && ijk[axis3]<p_max[axis3])
              //   v_fine(fine+jp+kp)=v_fine(fine-ip_s+jp+kp)
              //     + (v_pp - v_fine(fine-ip_s-ip_s+jp+kp))/3;
            }
          else
            {
              v_fine(fine)=v_fine(fine-ip_s)
                + (v_mp - v_fine(fine-ip_s-ip_s))/3;
              // if(ijk[axis2]<p_max[axis2])
              //   v_fine(fine+jp)=v_fine(fine-ip_s+jp)
              //     + (v_pp - v_fine(fine-ip_s-ip_s+jp))/3;
            }
        }
      else
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=v_fine(fine-ip_s)
                + (v_pm - v_fine(fine-ip_s-ip_s))/3;
              // if(ijk[axis3]<p_max[axis3])
              //   v_fine(fine+kp)=v_fine(fine-ip_s+kp)
              //     + (v_pp - v_fine(fine-ip_s-ip_s+kp))/3;
            }
          else
            {
              v_fine(fine)=v_fine(fine-ip_s)
                + (v_pp - v_fine(fine-ip_s-ip_s))/3;
            }
        }          
    }
  /* Set the value for the tangential direction.

     If we look at a slice in the i=constant plane.

          j-1      j      j+1

        ------- -------
       | f   f | F     |
   k-1 |   C   |   C   |   C
       | f   f | F     |
        ------- -------
       | f   f | F     |
   k   |   C   |   C   |   C
       | f   f | F     |
        ------- -------
       | f   f | F     |
   k+1 |   C   |   C   |   C
       | f   f | F     |
        ------- -------
               |
               |
               |
        Coarse-Fine Boundary

  where C are the coarse velocities, f are the fine velocities, and we
  interpolate to the F velocities.

 */
  else
    {
      const int axis3((axis+1)%3 != boundary_direction ? (axis+1)%3 : (axis+2)%3);
      const SAMRAI::hier::Index ip(pp[axis]),
        jp(boundary_positive ? pp[boundary_direction] : -pp[boundary_direction]),
        kp(pp[axis3]);

      SAMRAI::pdat::SideIndex center(fine);
      center.coarsen(SAMRAI::hier::Index(2,2,2));

      double v_m, v_p;
      quad_offset_interpolate(v(center+kp),v(center),v(center-kp),v_p,v_m);

      /* Be careful about using the right interpolation if the fine
       * points are not aligned with the coarse points.  There is some
       * double calls to quad_offset_interpolate going on, but fixing
       * that would require mucking with the iteration order in an
       * annoying way. */
      if(ijk[axis]%2==0)
        {
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=(8*v_m + 10*v_fine(fine-jp)
                            - 3*v_fine(fine-jp-jp))/15;
            }
          else
            {
              v_fine(fine)=(8*v_p + 10*v_fine(fine-jp)
                            - 3*v_fine(fine-jp-jp))/15;
            }
        }
      else
        {
          double vv_m, vv_p;
          quad_offset_interpolate(v(center+kp+ip),v(center+ip),
                                  v(center-kp+ip),vv_p,vv_m);
          if(ijk[axis3]%2==0)
            {
              v_fine(fine)=(4*(v_m+vv_m) + 10*v_fine(fine-jp)
                            - 3*v_fine(fine-jp-jp))/15;
            }
          else
            {
              v_fine(fine)=(4*(v_p+vv_p) + 10*v_fine(fine-jp)
                            - 3*v_fine(fine-jp-jp))/15;
            }
        }
    }
}
