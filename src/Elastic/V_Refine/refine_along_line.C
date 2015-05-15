#include "Elastic/V_Refine.h"

/* This assumes that the levels are always properly nested, so that we
   always have an extra grid space for interpolation.  So we only have
   to have a special case for physical boundaries, where we do not
   have an extra grid space. */

/* Maybe this has to be fixed when dvx/dy != 0 on the outer boundary
   because the approximation to the derivative is not accurate
   enough? */

/* Maybe in 3D we should include cross derivatives? */

double Elastic::V_Refine::refine_along_line
(SAMRAI::pdat::SideData<double> &v,
 const Gamra::Dir &ix,
 const Gamra::Dir &dim,
 const SAMRAI::hier::Index unit[],
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::pdat::SideIndex &coarse,
 const SAMRAI::hier::Box &coarse_box,
 const SAMRAI::geom::CartesianPatchGeometry &coarse_geom) const
{
  double result=v(coarse);

  for(Gamra::Dir d=ix.next(dim);d!=ix;d=d.next(dim))
    {
      const int sgn(fine[d]%2==0 ? -1 : 1);

      double dvx_dy;
      if(coarse[d]==coarse_box.lower(d)
         && coarse_geom.getTouchesRegularBoundary(d,0))
        {
          dvx_dy=sgn*(v(coarse+unit[d])-v(coarse))/4;
        }
      else if(coarse[d]==coarse_box.upper(d)
              && coarse_geom.getTouchesRegularBoundary(d,1))
        {
          dvx_dy=sgn*(v(coarse)-v(coarse-unit[d]))/4;
        }
      else
        {
          dvx_dy=sgn*(v(coarse+unit[d])-v(coarse-unit[d]))/8;
        }
      result+=dvx_dy;
    }
  return result;
}

/* Version for embedded boundaries */

double Elastic::V_Refine::refine_along_line
(SAMRAI::pdat::SideData<double> &v,
 const Gamra::Dir &ix,
 const Gamra::Dir &dim,
 const SAMRAI::hier::Index unit[],
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::pdat::SideIndex &coarse,
 const SAMRAI::pdat::SideData<double> &level_set_coarse) const
{
  if(level_set_coarse(coarse)<0)
    {
      /* FIXME: This should calculate the Dirichlet or Neumann
       * conditions properly as well as the proper weights.  For now,
       * hard code v=0 bc's. */
      double result(0);
      // double result(std::numeric_limits<double>::max());
      for(int d=(ix+1)%dim;d!=ix;d=(d+1)%dim)
        {
          if(level_set_coarse(coarse+unit[d])>=0)
            {
              result=v(coarse+unit[d]);
            }
          else if(level_set_coarse(coarse-unit[d])>=0)
            {
              result=v(coarse-unit[d]);
            }
        }
      return result;
    }

  double result=v(coarse);

  for(int d=(ix+1)%dim;d!=ix;d=(d+1)%dim)
    {
      const int sgn(fine[d]%2==0 ? -1 : 1);

      /* FIXME: This should calculate the Dirichlet or Neumann
       * conditions properly as well as the proper weights.  For now,
       * hard code d/dx=0 bc's. */
      // double v_plus(v(coarse)), v_minus(v(coarse));
      double dvx_dy;
      if(level_set_coarse(coarse+unit[d])<0)
        {
          dvx_dy=(v(coarse)-v(coarse-unit[d]))/4;
        }
      else if(level_set_coarse(coarse-unit[d])<0)
        {
          dvx_dy=(v(coarse+unit[d])-v(coarse))/4;
        }
      else
        {
          dvx_dy=(v(coarse+unit[d])-v(coarse-unit[d]))/8;
        }
      result+=sgn*dvx_dy;
    }
  return result;
}

