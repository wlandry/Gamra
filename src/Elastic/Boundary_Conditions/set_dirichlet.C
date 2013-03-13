#include "Constants.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_dirichlet
(SAMRAI::pdat::SideData<double> &v,
 const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_ptr,
 const SAMRAI::hier::Index unit[],
 const int &dim,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::hier::Box &gbox,
 const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom,
 const double *dx,
 const bool &homogeneous)
{
  for(int ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;
      SAMRAI::pdat::SideIterator s_end(gbox,ix,false);
      for(SAMRAI::pdat::SideIterator si(gbox,ix,true); si!=s_end; si++)
        {
          SAMRAI::pdat::SideIndex x(*si);

          for(int d=0;d<dim;++d)
            coord[d]=geom->getXLower()[d]
              + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

          /* For normal BC's, for the point just outside the
             boundary, set a sentinel value for normal dirichlet
             BC or the derivative for normal traction BC. */

          /* A point can be affected by multiple boundaries.  So we
             have to fall through for every point in case a point is
             on a corner but the normal direction is not dirichlet. */

          if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0)
             && is_dirichlet[ix][ix][0])
            {
              v(x)=boundary_value;
            }
          else if(x[ix]>pbox.upper(ix)+1
                  && geom->getTouchesRegularBoundary(ix,1)
                  && is_dirichlet[ix][ix][1])
            {
              v(x)=boundary_value;
            }
          /* If at the boundary line, set values for normal
           * components. */
          else if(x[ix]==pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0)
                  && is_dirichlet[ix][ix][0])
            {
              if(!homogeneous)
                v(x)=expression[ix][ix][0].eval(coord);
              else
                v(x)=0;
            }
          else if(x[ix]==pbox.upper(ix)+1
                  && geom->getTouchesRegularBoundary(ix,1)
                  && is_dirichlet[ix][ix][1])
            {
              if(!homogeneous)
                v(x)=expression[ix][ix][1].eval(coord);
              else
                v(x)=0;
            }
          /* Set tangential components. */
          else
            {
              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  const int ix_iy(index_map(ix,iy,dim));
                  if(x[iy]<pbox.lower(iy)
                     && geom->getTouchesRegularBoundary(iy,0))
                    {
                      if(is_dirichlet[ix][iy][0])
                        {
                          v(x)=-v(x+unit[iy]);
                          if(!homogeneous)
                            {
                              double coord_save(geom->getXLower()[iy]);
                              std::swap(coord[iy],coord_save);
                              if(have_faults)
                                v(x)+= -(*dv_mixed_ptr)(x+unit[iy],ix_iy+1);
                              v(x)+=2*expression[ix][iy][0].eval(coord);
                              std::swap(coord[iy],coord_save);
                            }
                        }
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      if(is_dirichlet[ix][iy][1])
                        {
                          v(x)=-v(x-unit[iy]);
                          if(!homogeneous)
                            {
                              double coord_save(geom->getXUpper()[iy]);
                              std::swap(coord[iy],coord_save);
                              if(have_faults)
                                v(x)+= -(*dv_mixed_ptr)(x-unit[iy],ix_iy);
                              v(x)+=2*expression[ix][iy][1].eval(coord);
                              std::swap(coord[iy],coord_save);
                            }
                        }
                    }
                }
            }
        }
    }
}
