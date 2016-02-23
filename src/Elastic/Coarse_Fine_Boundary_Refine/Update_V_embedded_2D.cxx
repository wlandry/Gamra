#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"
#include "quad_offset_interpolate.hxx"
#include "Constants.hxx"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::Coarse_Fine_Boundary_Refine::Update_V_embedded_2D
(const int &axis,
 const int &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::IntVector &ip,
 const SAMRAI::hier::IntVector &jp,
 const int &i,
 const int &j,
 const SAMRAI::pdat::SideData<double> &level_set,
 const SAMRAI::pdat::SideData<double> &,
 const SAMRAI::pdat::SideData<double> &v,
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
  */

  if(boundary_direction==axis)
    {
      SAMRAI::hier::IntVector ip_s(boundary_positive ? ip : -ip);
      SAMRAI::pdat::SideIndex coarse(fine-ip_s);
      coarse.coarsen(SAMRAI::hier::Index(2,2));

      double v_coarse_center;
      if(level_set(coarse+ip_s)>=0)
        {
          v_coarse_center=v(coarse+ip_s);
        }
      else
        {
          /* FIXME: We should do a proper interpolation. Assuming
           * v=0.  */
          v_coarse_center=0;
        }

      double v_coarse_p;
      if(level_set(coarse+ip_s+jp)>=0)
        {
          v_coarse_p=v(coarse+ip_s+jp);
        }
      else
        {
          /* FIXME: We should do a proper interpolation. Assuming
           * d/dx=0.  */
          v_coarse_p=v_coarse_center;
        }

      double v_coarse_m;
      if(level_set(coarse+ip_s-jp)>=0)
        {
          v_coarse_m=v(coarse+ip_s-jp);
        }
      else
        {
          v_coarse_m=v_coarse_center;
        }

      double v_p, v_m;
      quad_offset_interpolate(v_coarse_p,v_coarse_center,v_coarse_m,v_p,v_m);

      double v_coarse(j%2==0 ? v_m : v_p);
      /* FIXME: Need to do something about the fine points */
      v_fine(fine)=v_fine(fine-ip_s) + (v_coarse - v_fine(fine-ip_s-ip_s))/3;
    }
  /* Quadratic interpolation involving both coarse and fine grids for
     the tangential direction

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
      SAMRAI::pdat::SideIndex coarse(fine);
      coarse.coarsen(SAMRAI::hier::Index(2,2));
      SAMRAI::hier::Index jp_s(boundary_positive ? jp : -jp);

      double v_m;
      if(level_set(coarse)>=0)
        {
          v_m=v(coarse);
        }
      else
        {
          /* FIXME: We should do a proper interpolation. Assuming
           * v=0.  */
          v_m=0;
        }
      double v_coarse;
      if(i%2==0)
        {
          v_coarse=v_m;
        }
      else
        {
          double v_p;
          if(level_set(coarse+ip)>=0)
            {
              v_p=v(coarse+ip);
            }
          else
            {
              v_p=0;
            }
          v_coarse=(v_m+v_p)/2;
        }

      /* FIXME: Need to do something about the fine points */
      v_fine(fine)=(8*v_coarse + 10*v_fine(fine-jp_s)
                    - 3*v_fine(fine-jp_s-jp_s))/15;
    }
}
