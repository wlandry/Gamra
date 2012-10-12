#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

namespace {
  void update_offset_tangent_V_2D(const SAMRAI::pdat::SideIndex &fine,
                                  SAMRAI::pdat::SideIndex &center,
                                  const SAMRAI::hier::Index &ip,
                                  SAMRAI::hier::Index &jp_s,
                                  SAMRAI::pdat::SideData<double> &v,
                                  SAMRAI::pdat::SideData<double> &v_fine)
  {
    double v_coarse;
    if(v(center+ip+ip)==boundary_value)
      {
        v_coarse=(-v(center-ip) + 6*v(center)
                  + 3*v(center+ip))/8;
      }
    else if(v(center-ip)==boundary_value)
      {
        v_coarse=(-v(center+ip+ip) + 6*v(center+ip)
                  + 3*v(center))/8;
      }
    else
      {
        v_coarse=(-v(center-ip) + 9*v(center)
                  + 9*v(center+ip) - v(center+ip+ip))/16;
      }
    v_fine(fine)=(8*v_coarse + 10*v_fine(fine-jp_s)
                  - 3*v_fine(fine-jp_s-jp_s))/15;
  }
}

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::V_Boundary_Refine::Update_V_2D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::Index &ip, const SAMRAI::hier::Index &jp,
 int &i, int &j,
 const int &i_max,
 const int &j_min,
 const int &j_max,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_fine) const
{
  /* Neumann'ish conditions for the normal direction

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
  */

  if(boundary_direction==axis)
    {
      SAMRAI::hier::Index ip_s(boundary_positive ? ip : -ip);
      SAMRAI::pdat::SideIndex center(fine-ip_s);
      center.coarsen(SAMRAI::hier::Index(2,2));

      double v_p, v_m;
      quad_offset_interpolate(v(center+ip_s+jp),v(center+ip_s),
                              v(center+ip_s-jp),v_p,v_m);

      if(j%2==0)
        {
          v_fine(fine)=v_fine(fine-ip_s) + (v_m - v_fine(fine-ip_s-ip_s))/3;
          if(j<j_max)
            {
              v_fine(fine+jp)=v_fine(fine-ip_s+jp)
                + (v_p - v_fine(fine-ip_s-ip_s+jp))/3;
            }
          ++j;
        }
      else
        {
          v_fine(fine)=v_fine(fine-ip_s) + (v_p - v_fine(fine-ip_s-ip_s))/3;
        }
    }
  /* Neumann'ish conditions for the tangential direction

      i-1      i      i+1

        ------- -------
       |       |       |
   j-1 C       C       C
       F   F   F   F   F
        ------- -------    --- Coarse-Fine Boundary
       f   f   f   f   f
   j   C       C       C
       f   f   f   f   f
        ------- -------
       f   f   f   f   f
   j+1 C       C       C
       f   f   f   f   f
        ------- -------

    C are the coarse velocities, f are the interior fine velocities,
    and F are the boundary fine velocities that we need to set.
 */


  else
    {
      SAMRAI::pdat::SideIndex center(fine);
      center.coarsen(SAMRAI::hier::Index(2,2));

      SAMRAI::hier::Index jp_s(boundary_positive ? jp : -jp);

      if(i%2==0)
        {
          v_fine(fine)=(8*v(center) + 10*v_fine(fine-jp_s)
                        - 3*v_fine(fine-jp_s-jp_s))/15;

          if(i<i_max)
            {
              update_offset_tangent_V_2D(fine+ip,center,ip,jp_s,v,v_fine);
              /* Since we update two points on 'i' at once, we
                 increment 'i' again.  This is ok, since the box in
                 the 'j' direction is defined to be only one cell
                 wide */
              ++i;
            }
        }
      else
        {
          update_offset_tangent_V_2D(fine,center,ip,jp_s,v,v_fine);
        }
    }
}
