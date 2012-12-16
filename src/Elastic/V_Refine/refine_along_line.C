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
 const int &axis,
 const int &dim,
 const SAMRAI::hier::Index pp[],
 const SAMRAI::pdat::SideIndex &fine,
 const SAMRAI::pdat::SideIndex &coarse,
 const SAMRAI::hier::Box &coarse_box,
 const SAMRAI::geom::CartesianPatchGeometry &coarse_geom) const
{
  double result=v(coarse);

  for(int d=(axis+1)%dim;d!=axis;d=(d+1)%dim)
    {
      const int sgn(fine[d]%2==0 ? -1 : 1);

      double dvx_dy;
      if(coarse[d]==coarse_box.lower(d)
         && coarse_geom.getTouchesRegularBoundary(d,0))
        {
          dvx_dy=sgn*(v(coarse+pp[d])-v(coarse))/4;
        }
      else if(coarse[d]==coarse_box.upper(d)
              && coarse_geom.getTouchesRegularBoundary(d,1))
        {
          dvx_dy=sgn*(v(coarse)-v(coarse-pp[d]))/4;
        }
      else
        {
          dvx_dy=sgn*(v(coarse+pp[d])-v(coarse-pp[d]))/8;
        }
      result+=dvx_dy;
    }
  return result;
}

