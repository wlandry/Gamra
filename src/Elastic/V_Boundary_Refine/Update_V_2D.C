#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

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
  /* Dirichlet conditions for the normal direction

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

      const double v_center=
        (-v(center-ip_s) + 6*v(center) + 3*v(center+ip_s))/8;
      const double v_plus=
        (-v(center-ip_s+jp) + 6*v(center+jp) + 3*v(center+ip_s+jp))/8;
      const double v_minus=
        (-v(center-ip_s-jp) + 6*v(center-jp) + 3*v(center+ip_s-jp))/8;
                            
      double v_p, v_m;
      quad_offset_interpolate(v_plus,v_center,v_minus,v_p,v_m);

      if(j%2==0)
        {
          v_fine(fine)=v_m;
          if(j<j_max)
            v_fine(fine+jp)=v_p;
          ++j;
        }
      else
        {
          v_fine(fine)=v_p;
        }

      SAMRAI::tbox::plog << "Update V "
                         << fine << " "
                         << center << " "
                         << ip_s << " "
                         << jp << " "
                         << v_fine(fine) << " "
                         << v(center-ip_s) << " "
                         << v(center) << " "
                         << v(center+ip_s) << " "
                         << v_center << " "
                         << v_plus << " "
                         << v_minus << " "
                         << "\n";
    }
  /* Dirichlet conditions for the tangential direction

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
          v_fine(fine)=(8*v(center+jp_s) + 10*v_fine(fine-jp_s)
                        - 3*v_fine(fine-jp_s-jp_s))/15;

          if(i<i_max)
            {
              double v_coarse;
              if(v(center+ip+ip+jp_s)==boundary_value)
                {
                  v_coarse=(-v(center+jp_s-ip) + 6*v(center+jp_s)
                            + 3*v(center+jp_s+ip))/8;
                }
              else if(v(center-ip+jp_s)==boundary_value)
                {
                  v_coarse=(-v(center+jp_s+ip+ip) + 6*v(center+jp_s+ip)
                            + 3*v(center+jp_s))/8;
                }
              else
                {
                  v_coarse=(-v(center+jp_s-ip) + 9*v(center+jp_s)
                            + 9*v(center+jp_s+ip) - v(center+jp_s+ip+ip))/16;
                }
              v_fine(fine+ip)=(8*v_coarse + 10*v_fine(fine-jp_s+ip)
                               - 3*v_fine(fine-jp_s-jp_s+ip))/15;


              /* Since we update two points on 'i' at once, we
                 increment 'i' again.  This is ok, since the box in
                 the 'j' direction is defined to be only one cell
                 wide */
              ++i;
            }
        }
      else
        {
          double v_coarse;
          if(v(center+ip+ip+jp_s)==boundary_value)
            {
              v_coarse=(-v(center+jp_s-ip) + 6*v(center+jp_s)
                        + 3*v(center+jp_s+ip))/8;
            }
          else if(v(center-ip+jp_s)==boundary_value)
            {
              v_coarse=(-v(center+jp_s+ip+ip) + 6*v(center+jp_s+ip)
                        + 3*v(center+jp_s))/8;
            }
          else
            {
              v_coarse=(-v(center+jp_s-ip) + 9*v(center+jp_s)
                        + 9*v(center+jp_s+ip) - v(center+jp_s+ip+ip))/16;
            }
          v_fine(fine)=(8*v_coarse + 10*v_fine(fine-jp_s)
                        - 3*v_fine(fine-jp_s-jp_s))/15;
        }
    }
}
