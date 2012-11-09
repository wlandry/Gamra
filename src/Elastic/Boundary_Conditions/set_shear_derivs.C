#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_shear_derivs
(SAMRAI::pdat::SideData<double> &v,
 SAMRAI::hier::Index pp[], const int &dim,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::hier::Box &gbox,
 const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry> geom,
 const double *dx,
 const bool &homogeneous)
{
  for(int ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;
      for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
        {
          SAMRAI::pdat::SideIndex x(*si);
              
          /* Set the derivative. */
          if(x[ix]>=pbox.lower(ix) && x[ix]<=pbox.upper(ix)+1)
            {
              for(int d=0;d<dim;++d)
                coord[d]=geom->getXLower()[d]
                  + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  if(x[iy]<pbox.lower(iy)
                     && geom->getTouchesRegularBoundary(iy,0))
                    {
                      if(!is_dirichlet[ix][iy][0])
                        {
                          double coord_save(geom->getXLower()[iy]);
                          std::swap(coord[iy],coord_save);
                          SAMRAI::pdat::SideIndex
                            y(x+pp[iy],iy,SAMRAI::pdat::SideIndex::Lower);
                          const double duyx((v(y)-v(y-pp[ix]))/dx[ix]);
                          v(x)=v(x+pp[iy]) + duyx*dx[iy];
                          if(!homogeneous)
                            {
                              int iz=(iy+1)%dim;
                              if(iz==ix)
                                iz=(iz+1)%dim;
                              v(x)-=shear_derivs[ix][iy][0].Eval()*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      if(!is_dirichlet[ix][iy][1])
                        {
                          double coord_save(geom->getXUpper()[iy]);
                          std::swap(coord[iy],coord_save);
                          SAMRAI::pdat::SideIndex
                            y(x,iy,SAMRAI::pdat::SideIndex::Lower);
                          const double duyx((v(y)-v(y-pp[ix]))/dx[ix]);
                          v(x)=v(x-pp[iy]) - duyx*dx[iy];
                          if(!homogeneous)
                            {
                              int iz=(iy+1)%dim;
                              if(iz==ix)
                                iz=(iz+1)%dim;
                              v(x)+=shear_derivs[ix][iy][1].Eval()*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                    }
                }
            }
        }
    }
}

