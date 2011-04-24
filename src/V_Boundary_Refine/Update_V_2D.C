#include "V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Boundary.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void SAMRAI::geom::V_Boundary_Refine::Update_V_2D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const pdat::SideIndex &fine,
 const hier::Index &ip, const hier::Index &jp,
 int &i, int &j,
 const int &i_max,
 const int &j_min,
 const int &j_max,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_fine) const
{
  pdat::SideIndex center(fine);
  center.coarsen(hier::Index(2,2));

  /* Set the derivative for the normal direction

         i-1      i       i+1

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

    So we need to set F such that the derivative at the coarse-fine
    boundary is coorrect.  We can compute the derivative at the coarse
    points on the boundary, and then use quadratic interpolation to
    get the derivative at the fine points on the boundary.

 */
  if(boundary_direction==axis)
    {
      /* Return early if we are at j==j_max, because that is a corner
         that we do not care about.  We could also skip if j==j_min as
         long as we do not have to do j_min+1. We do not really have
         to skip these since we are guaranteed to have valid data for
         those "past the end" points since they are needed for
         pressure refinement.  */
      if(j==j_max || (j==j_min && j%2!=0))
        return;
      /* Compute the derivative at the nearest three coarse points and
         then interpolate */

      hier::Index ip_s(boundary_positive ? ip : -ip);

      const double dv_plus=v(center+jp+ip_s)-v(center+jp-ip_s);
      const double dv_minus=v(center-jp+ip_s)-v(center-jp-ip_s);
      const double dv_center=v(center+ip_s)-v(center-ip_s);

      double dv_fine_minus, dv_fine_plus;

      quad_offset_interpolate(dv_plus,dv_minus,dv_center,
                              dv_fine_plus,dv_fine_minus);

      hier::Index offset(ip_s*2);

      if(j%2==0)
        {
          v_fine(fine)=v_fine(fine-offset) + dv_fine_minus/2;
          v_fine(fine+jp)=v_fine(fine-offset+jp) + dv_fine_plus/2;
          /* Since we update two points on j at once, we increment j
             again.  This is ok, since the box in the 'i' direction is
             defined to be only one cell wide */
          ++j;
        }
      else
        {
          v_fine(fine)=v_fine(fine-offset) + dv_fine_plus/2;
        }          
    }
  /* Set the value for the tangential direction

         i-1      i       i+1

        -f-C-f- -F-C---    C
       |       |       |
   j-1 |       |       |    
       |       |       |
        -f-C-f- -F-C---    C
       |       |       |
   j   |       |       |    
       |       |       |
        -f-C-f- -F-C---    C
       |       |       |
   j+1 |       |       |    
       |       |       |
        -f-C-f- -F-C---    C
               |
               |
               |
        Coarse-Fine Boundary

    C are the coarse velocities, f are the interior fine velocities,
    and F are the boundary fine velocities that we need to set.  So we
    use quadratic interpolation from C to F.
 */
  else
    {
      double v_center, v_plus;
      hier::Index jp_s(boundary_positive ? jp : -jp);

      v_center=
        quad_offset_interpolate(v(center-jp_s),v(center),v(center+jp_s));

      if(i%2==0)
        {
          v_fine(fine)=v_center;

          if(i<i_max)
            {
              /* This is a bit inefficient, because we compute v_plus
               * twice.  Once for the in-between point, and again
               * later for the actual point. */

              v_plus=quad_offset_interpolate(v(center+ip-jp_s),v(center+ip),
                                             v(center+ip+jp_s));
              v_fine(fine+ip)=(v_center+v_plus)/2;

              /* Since we update two points on 'i' at once, we increment 'i' again.
                 This is ok, since the box in the 'j' direction is defined to be
                 only one cell wide */
              ++i;
            }
        }
      else
        {
          v_plus=quad_offset_interpolate(v(center+ip-jp_s),v(center+ip),
                                         v(center+ip+jp_s));
          v_fine(fine)=(v_center+v_plus)/2;
        }
    }
}
