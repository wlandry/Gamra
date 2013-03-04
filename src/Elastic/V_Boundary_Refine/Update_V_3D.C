#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::V_Boundary_Refine::Update_V_3D
(const int &ix,
 const int &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::Index unit[],
 const SAMRAI::hier::Index &ijk,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::geom::CartesianPatchGeometry &geom,
 const boost::shared_ptr<SAMRAI::pdat::CellData<double> > &dv_diagonal,
 const boost::shared_ptr<SAMRAI::pdat::CellData<double> > &dv_diagonal_fine,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_fine,
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


     So need to do a diagonal interpolation of coarse values on the
     face and then an interpolation using coarse and fine values to
     get inside the cube.

  */

  if(boundary_direction==ix)
    {
      const int iy((ix+1)%3), iz((ix+2)%3);
      const SAMRAI::hier::Index ip_s(boundary_positive ? unit[ix] : -unit[ix]);

      SAMRAI::pdat::SideIndex coarse(fine-ip_s);
      coarse.coarsen(SAMRAI::hier::Index(2,2,2));
      coarse+=ip_s;

      SAMRAI::hier::Index jp(unit[iy]),kp(unit[iz]);
      const int ijk_mod_y(ijk[iy]%2), ijk_mod_z(ijk[iz]%2);
      int lower_y, upper_y, lower_z, upper_z;
      if(ijk_mod_y==0)
        {
          lower_y=pbox.lower(iy);
          upper_y=pbox.upper(iy);
        }
      else
        {
          lower_y=pbox.upper(iy);
          upper_y=pbox.lower(iy);
          jp=-jp;
        }

      if(ijk_mod_z==0)
        {
          lower_z=pbox.lower(iz);
          upper_z=pbox.upper(iz);
        }
      else
        {
          lower_z=pbox.upper(iz);
          upper_z=pbox.lower(iz);
          kp=-kp;
        }

      int plus,minus;
      if(ijk_mod_y==0)
        {
          if(ijk_mod_z==0)
            {
              plus=4;
              minus=6;
            }
          else
            {
              plus=5;
              minus=7;
            }
        }
      else
        {
          if(ijk_mod_z==0)
            {
              plus=7;
              minus=5;
            }
          else
            {
              plus=6;
              minus=4;
            }
        }
          
      double v_coarse;
      if(ijk[iy]==lower_y && geom.getTouchesRegularBoundary(iy,ijk_mod_y)
         && ijk[iz]==lower_z && geom.getTouchesRegularBoundary(iz,ijk_mod_z))
        {
          v_coarse=(5*v(coarse) - v(coarse+jp+kp))/4;
          if(!is_residual)
            v_coarse-=((*dv_mixed)(coarse+jp+kp,minus) - (*dv_mixed)(coarse,plus))/4;
        }
      else if(ijk[iy]==upper_y
              && geom.getTouchesRegularBoundary(iy,(ijk_mod_y+1)%2)
              && ijk[iz]==upper_z
              && geom.getTouchesRegularBoundary(iz,(ijk_mod_z+1)%2))
        {
          v_coarse=(3*v(coarse) + v(coarse-jp-kp))/4;
          if(!is_residual)
            v_coarse+=((*dv_mixed)(coarse-jp-kp,plus) - (*dv_mixed)(coarse,minus))/4;
        }
      else
        {
          v_coarse=
            quad_offset_interpolate(v(coarse-jp-kp),v(coarse),v(coarse+jp+kp));
          if(!is_residual)
            v_coarse+=quad_offset_correction((*dv_mixed)(coarse-jp-kp,plus),
                                             (*dv_mixed)(coarse,minus),
                                             (*dv_mixed)(coarse,plus),
                                             (*dv_mixed)(coarse+jp+kp,minus));
        }
      if(!is_residual)
        {
          SAMRAI::pdat::CellIndex cell(coarse);
          v_coarse-=(*dv_mixed_fine)(fine-ip_s,plus)
            - (boundary_positive ? -(*dv_diagonal)(cell-ip_s,ix) :
               (*dv_diagonal)(cell,ix));
        }

      v_fine(fine)=v_fine(fine-ip_s) + (v_coarse - v_fine(fine-ip_s-ip_s))/3;
      if(!is_residual)
        {
          SAMRAI::pdat::CellIndex cell(fine);
          if(boundary_positive)
            {
              v_fine(fine)+=(*dv_diagonal_fine)(cell-ip_s,ix)
                - (*dv_diagonal_fine)(cell-ip_s-ip_s,ix)/3;
            }
          else
            {
              v_fine(fine)-=(*dv_diagonal_fine)(cell,ix)
                - (*dv_diagonal_fine)(cell-ip_s,ix)/3;
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
      const int iz((ix+1)%3 != boundary_direction ? (ix+1)%3 : (ix+2)%3);
      const SAMRAI::hier::Index ip(unit[ix]),
        jp(boundary_positive ? unit[boundary_direction] : -unit[boundary_direction]),
        kp(ijk[iz]%2==0 ? -unit[iz] : unit[iz]);

      SAMRAI::pdat::SideIndex coarse(fine);
      coarse.coarsen(SAMRAI::hier::Index(2,2,2));

      double v_pm(quad_offset_interpolate(v(coarse+kp),v(coarse),v(coarse-kp)));

      /* Be careful about using the right interpolation if the fine
       * points are not aligned with the coarse points.  There is some
       * double calls to quad_offset_interpolate going on, but fixing
       * that would require mucking with the iteration order in an
       * annoying way. */
      if(ijk[ix]%2==0)
        {
          v_fine(fine)=(8*v_pm + 10*v_fine(fine-jp)
                        - 3*v_fine(fine-jp-jp))/15;
        }
      else
        {
          double vv_pm(quad_offset_interpolate(v(coarse+kp+ip),v(coarse+ip),
                                               v(coarse-kp+ip)));
          v_fine(fine)=(4*(v_pm + vv_pm) + 10*v_fine(fine-jp)
                        - 3*v_fine(fine-jp-jp))/15;
        }
    }
}
