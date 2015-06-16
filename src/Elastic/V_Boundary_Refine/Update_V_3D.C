#include "Elastic/V_Boundary_Refine.h"
#include "quad_offset_interpolate.h"
#include "Constants.h"

/* This is written from the perspective of axis==x.  For axis==y, we
   switch i and j and everything works out. */
void Elastic::V_Boundary_Refine::Update_V_3D
(const Gamra::Dir &ix,
 const Gamra::Dir &boundary_direction,
 const bool &boundary_positive,
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::hier::IntVector unit[],
 const SAMRAI::hier::Index &ijk,
 const SAMRAI::hier::Box &coarse_box,
 const SAMRAI::hier::Index &fine_min,
 const SAMRAI::hier::Index &fine_max,
 const SAMRAI::geom::CartesianPatchGeometry &coarse_geom,
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
      const Gamra::Dir iy(ix.next(3));
      const Gamra::Dir iz(iy.next(3));
      const SAMRAI::hier::IntVector ip_s(boundary_positive ? unit[ix]
                                         : -unit[ix]);

      SAMRAI::pdat::SideIndex coarse(fine-ip_s);
      coarse.coarsen(SAMRAI::hier::Index(2,2,2));
      coarse+=ip_s;

      SAMRAI::hier::IntVector jp(unit[iy]),kp(unit[iz]);
      const int ijk_mod_y(ijk[iy]%2), ijk_mod_z(ijk[iz]%2);
      int lower_y, upper_y, lower_z, upper_z;
      if(ijk_mod_y==0)
        {
          lower_y=coarse_box.lower(iy);
          upper_y=coarse_box.upper(iy);
        }
      else
        {
          lower_y=coarse_box.upper(iy);
          upper_y=coarse_box.lower(iy);
          jp=-jp;
        }

      if(ijk_mod_z==0)
        {
          lower_z=coarse_box.lower(iz);
          upper_z=coarse_box.upper(iz);
        }
      else
        {
          lower_z=coarse_box.upper(iz);
          upper_z=coarse_box.lower(iz);
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
          
      /* We need to check when interpolating whether the stencil goes
         off a boundary.  Fault corrections are not defined outside
         the boundary.  Also, boundary values are not defined at the
         outside corner.  So we use a simpler interpolation there. */
      double v_coarse;
      if((coarse[iy]==lower_y
          && coarse_geom.getTouchesRegularBoundary(iy,ijk_mod_y))
         || (coarse[iz]==lower_z
             && coarse_geom.getTouchesRegularBoundary(iz,ijk_mod_z)))
        {
          v_coarse=(5*v(coarse) - v(coarse+jp+kp))/4;
          /// Correct v(coarse+jp+kp) to v(coarse)
          if(have_faults() && !is_residual)
            v_coarse-=
              ((*dv_mixed)(coarse+jp+kp,minus) - (*dv_mixed)(coarse,plus))/4;
        }
      else if((coarse[iy]==upper_y
               && coarse_geom.getTouchesRegularBoundary(iy,(ijk_mod_y+1)%2))
              || (coarse[iz]==upper_z
                  && coarse_geom.getTouchesRegularBoundary(iz,(ijk_mod_z+1)%2)))
        {
          v_coarse=(3*v(coarse) + v(coarse-jp-kp))/4;
          /// Correct v(coarse-jp-kp) to v(coarse)
          if(have_faults() && !is_residual)
            v_coarse+=
              ((*dv_mixed)(coarse-jp-kp,plus) - (*dv_mixed)(coarse,minus))/4;
        }
      else
        {
          v_coarse=
            quad_offset_interpolate(v(coarse-jp-kp),v(coarse),v(coarse+jp+kp));
          /// Correct to v(coarse)
          if(have_faults() && !is_residual)
            v_coarse+=quad_offset_correction((*dv_mixed)(coarse-jp-kp,plus),
                                             (*dv_mixed)(coarse,minus),
                                             (*dv_mixed)(coarse,plus),
                                             (*dv_mixed)(coarse+jp+kp,minus));
        }

      /// Previous corrections translated the result to v(coarse)
      /// This last bit translates to v(coarse-ip_s) and then to v(fine-ip_s);
      if(have_faults() && !is_residual)
        {
          SAMRAI::pdat::CellIndex cell(coarse);
          v_coarse+=(boundary_positive ? -(*dv_diagonal)(cell-ip_s,ix) :
                     (*dv_diagonal)(cell,ix))
            - (*dv_mixed_fine)(fine-ip_s,plus);
        }
      v_fine(fine)=v_fine(fine-ip_s) + (v_coarse - v_fine(fine-ip_s*2))/3;
      if(have_faults() && !is_residual)
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
      const Gamra::Dir iz(ix.next(3) != boundary_direction
                          ? ix.next(3) : ix.next(3).next(3));
      const SAMRAI::hier::Index ip(unit[ix]),
        jp(boundary_positive ? unit[boundary_direction]
           : -unit[boundary_direction]),
        kp(ijk[iz]%2==0 ? -unit[iz] : unit[iz]);

      SAMRAI::pdat::SideIndex coarse(fine);
      coarse.coarsen(SAMRAI::hier::Index(2,2,2));

      double v_coarse(quad_offset_interpolate(v(coarse+kp),v(coarse),
                                              v(coarse-kp)));
      const Gamra::Dir dim(3);
      const Gamra::Dir ix_iz(index_map(ix,iz,dim));
      if(have_faults() && !is_residual)
        {
          if(ijk[iz]%2==0)
            {
              v_coarse+=quad_offset_correction((*dv_mixed)(coarse+kp,ix_iz),
                                               (*dv_mixed)(coarse,ix_iz+1),
                                               (*dv_mixed)(coarse,ix_iz),
                                               (*dv_mixed)(coarse-kp,ix_iz+1));
            }
          else
            {
              v_coarse+=quad_offset_correction((*dv_mixed)(coarse+kp,ix_iz+1),
                                               (*dv_mixed)(coarse,ix_iz),
                                               (*dv_mixed)(coarse,ix_iz+1),
                                               (*dv_mixed)(coarse-kp,ix_iz));
            }
        }      
      double v_m(v_fine(fine-jp)), v_mm(v_fine(fine-jp-jp));

      if(have_faults() && !is_residual)
        {
          int ix_iy_in(index_map(ix,boundary_direction,dim)
                       + (boundary_positive ? 1 : 0));
          int ix_iy_out(index_map(ix,boundary_direction,dim)
                        + (boundary_positive ? 0 : 1));

          double v_m_correction(-(*dv_mixed_fine)(fine,ix_iy_in)
                                + (*dv_mixed_fine)(fine-jp,ix_iy_out));
          v_m+= v_m_correction;
          v_mm+= -(*dv_mixed_fine)(fine-jp,ix_iy_in)
            + (*dv_mixed_fine)(fine-jp-jp,ix_iy_out)
            + v_m_correction;
        }

      /* Numbering determined by Elastic::FAC::compute_intersections_3D */

      int directions[2][2][2]={{{4,5},{7,6}},{{4,7},{5,6}}};
      int direction=directions[boundary_direction==(ix+1)%dim ? 0 : 1]
        [boundary_positive ? 0 : 1][ijk[iz]%2];

      /* Be careful about using the right interpolation if the fine
       * points are not aligned with the coarse points.  There is some
       * double calls to quad_offset_interpolate going on, but fixing
       * that would require mucking with the iteration order in an
       * annoying way. */

      if(ijk[ix]%2==0)
        {
          if(have_faults() && !is_residual)
            v_coarse+= -(*dv_mixed_fine)(fine,direction);

          v_fine(fine)=(8*v_coarse + 10*v_m - 3*v_mm)/15;
        }
      else
        {
          double v_coarse_p(quad_offset_interpolate
                            (v(coarse+kp+ip),v(coarse+ip),v(coarse-kp+ip)));
                             
          if(have_faults() && !is_residual)
            {
              SAMRAI::pdat::CellIndex cell(fine), coarse_cell(coarse);

              double coarse_to_fine[]={0,0};
              if(fine_min(ix)!=fine(ix))
                {
                  coarse_to_fine[0]=-(*dv_mixed_fine)(fine-ip,direction)
                    + (*dv_diagonal_fine)(cell-ip,ix);
                }

              if(fine_max(ix)!=fine(ix))
                {
                  coarse_to_fine[1]=-(*dv_mixed_fine)(fine+ip,direction)
                    - (*dv_diagonal_fine)(cell,ix);
                }

              /// If the fine mesh does not cover the sides where the
              /// coarse points are, then we use the correction from
              /// the other coarse point and then correct with
              /// dv_diagonal.
              if(fine_min(ix)==fine(ix))
                coarse_to_fine[0]=coarse_to_fine[1]
                  + (*dv_diagonal)(coarse_cell,ix);

              if(fine_max(ix)==fine(ix))
                coarse_to_fine[1]=coarse_to_fine[0]
                  - (*dv_diagonal)(coarse_cell,ix);

              v_coarse+=coarse_to_fine[0];
                
              v_coarse_p+=coarse_to_fine[1];
              if(ijk[iz]%2==0)
                {
                  v_coarse_p+=
                    quad_offset_correction((*dv_mixed)(coarse+kp+ip,ix_iz),
                                           (*dv_mixed)(coarse+ip,ix_iz+1),
                                           (*dv_mixed)(coarse+ip,ix_iz),
                                           (*dv_mixed)(coarse-kp+ip,ix_iz+1));
                }
              else
                {
                  v_coarse_p+=
                    quad_offset_correction((*dv_mixed)(coarse+kp+ip,ix_iz+1),
                                           (*dv_mixed)(coarse+ip,ix_iz),
                                           (*dv_mixed)(coarse+ip,ix_iz+1),
                                           (*dv_mixed)(coarse-kp+ip,ix_iz));
                }
            }
          v_fine(fine)=(4*(v_coarse + v_coarse_p) + 10*v_m - 3*v_mm)/15;
        }
    }
}
