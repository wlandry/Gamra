#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void SAMRAI::geom::Elastic::V_Boundary_Refine::Update_V_2D
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
  /* Dirichlet conditions for the normal direction

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

      Set F to the sentinel value, and the value just inside is
      interpolated.
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

      v_fine(fine)=boundary_value;
      if(j%2==0)
        {
          v_fine(fine+jp)=boundary_value;
          ++j;
        }
    }
  /* Neumann conditions for the tangential direction

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
    and F are the boundary fine velocities that we need to set.
 */


  else
    {
      pdat::SideIndex center(fine);
      center.coarsen(hier::Index(2,2));

      hier::Index jp_s(boundary_positive ? jp : -jp);

      const double dv_minus=v(center-jp_s)-v(center);

      if(i%2==0)
        {
          v_fine(fine)=v_fine(fine-jp_s) - dv_minus/2;

          if(i<i_max)
            {
              /* This is a bit inefficient, because we compute v_plus
               * twice.  Once for the in-between point, and again
               * later for the actual point. */

              const double dv_plus=v(center+ip-jp_s)-v(center+ip);
              v_fine(fine+ip)=v_fine(fine+ip-jp_s) - (dv_minus + dv_plus)/4;

              /* Since we update two points on 'i' at once, we
                 increment 'i' again.  This is ok, since the box in
                 the 'j' direction is defined to be only one cell
                 wide */
              ++i;
            }
        }
      else
        {
          const double dv_plus=v(center+ip-jp_s)-v(center+ip);
          v_fine(fine)=v_fine(fine-jp_s) - (dv_minus + dv_plus)/4;
        }
    }
}
