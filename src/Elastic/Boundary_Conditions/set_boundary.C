#include "Constants.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_boundary
(const SAMRAI::hier::Patch& patch, const int &v_id, const bool &rhs)
{
  SAMRAI::hier::Box pbox=patch.getBox();

  SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
    geom = patch.getPatchGeometry();
  const SAMRAI::tbox::Dimension Dim(patch.getDim());
  const int dim(patch.getDim().getValue());

  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(Dim));
  SAMRAI::hier::Index pp[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    pp[i][i]=1;

  /* This should really get read from the input file. */
  const double *dx=geom->getDx();

  SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
    v_ptr = patch.getPatchData(v_id);
  SAMRAI::pdat::SideData<double> &v(*v_ptr);

  SAMRAI::hier::Box gbox=v.getGhostBox();
  for(int ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;
      try
        {
          for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
            {
              SAMRAI::pdat::SideIndex x(*si);

              for(int d=0;d<dim;++d)
                coord[d]=geom->getXLower()[d]
                  + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

              /* For normal BC's, for the point just outside the
                 boundary, set a sentinel value for normal dirichlet
                 BC or the derivative for normal neumann BC. */
              if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
                {
                  if(is_dirichlet[ix][ix][0])
                    v(x)=boundary_value;
                  else
                    v(x)=v(x+pp[ix]*2) - neumann[ix][ix][0].Eval()*2*dx[ix];
                }
              else if(x[ix]>pbox.upper(ix)+1
                      && geom->getTouchesRegularBoundary(ix,1))
                {
                  if(is_dirichlet[ix][ix][1])
                    v(x)=boundary_value;
                  else
                    v(x)=v(x-pp[ix]*2) + neumann[ix][ix][1].Eval()*2*dx[ix];
                }
              /* If at the boundary line, set values for normal
               * components. */
              else if(x[ix]==pbox.lower(ix)
                      && geom->getTouchesRegularBoundary(ix,0)
                      && !rhs && is_dirichlet[ix][ix][0])
                {
                  v(x)=dirichlet[ix][ix][0].Eval();
                }
              else if(x[ix]==pbox.upper(ix)+1
                      && geom->getTouchesRegularBoundary(ix,1)
                      && !rhs && is_dirichlet[ix][ix][1])
                {
                  v(x)=dirichlet[ix][ix][1].Eval();
                }
              /* Set tangential components. */
              else
                {
                  for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                    {
                      if(x[iy]<pbox.lower(iy)
                         && geom->getTouchesRegularBoundary(iy,0))
                        {
                          double coord_save(geom->getXLower()[iy]);
                          std::swap(coord[iy],coord_save);
                          if(is_dirichlet[ix][iy][0])
                            {
                              if(!rhs)
                                v(x)=2*dirichlet[ix][iy][0].Eval()
                                  - v(x+pp[iy]);
                              else
                                abort();
                              /* We need to coarsen the dirichlet
                               * value and store it somewhere. */
                            }
                          else
                            {
                              v(x)=v(x+pp[iy])
                                - neumann[ix][iy][0].Eval()*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                      else if(x[iy]>pbox.upper(iy)
                              && geom->getTouchesRegularBoundary(iy,1))
                        {
                          double coord_save(geom->getXUpper()[iy]);
                          std::swap(coord[iy],coord_save);
                          if(is_dirichlet[ix][iy][1])
                            {
                              if(!rhs)
                                v(x)=2*dirichlet[ix][iy][1].Eval()
                                  - v(x-pp[iy]);
                              else
                                abort();
                              /* We need to coarsen the dirichlet
                               * value and store it somewhere. */
                            }
                          else
                            {
                              v(x)=v(x-pp[iy])
                                + neumann[ix][iy][1].Eval()*dx[iy];
                            }
                          std::swap(coord[iy],coord_save);
                        }
                    }
                }
            }
          /* Fix up the edges.  This has to be done in a different
             loop, because the values on the faces will be used to
             compute values in the edges and corners.  In 2D, I am not
             sure that this is needed.  It is certainly not needed if
             we have dirichlet conditions on the boundary.  It is
             definitely needed in 3D.  We use the d/dx conditions
             here, though we could use the d/dy conditions and get the
             same number, as long as we have pure Neumann boundary
             conditions, and not Robin. */

          for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
            {
              SAMRAI::pdat::SideIndex x(*si);
              if((x[ix]<pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0)
                  && !is_dirichlet[ix][ix][0])
                 || (x[ix]>pbox.upper(ix)+1
                     && geom->getTouchesRegularBoundary(ix,1)
                     && !is_dirichlet[ix][ix][1]))
                {
                  for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                    {
                      for(int d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);
                      if(x[iy]<pbox.lower(iy) 
                         && geom->getTouchesRegularBoundary(iy,0))
                        {
                          double coord_save(geom->getXLower()[iy]);
                          std::swap(coord[iy],coord_save);
                          v(x)=v(x+pp[iy]) - neumann[ix][iy][0].Eval()*dx[iy];
                          std::swap(coord[iy],coord_save);
                        }
                      else if(x[iy]>pbox.upper(iy) 
                              && geom->getTouchesRegularBoundary(iy,1))
                        {
                          double coord_save(geom->getXUpper()[iy]);
                          std::swap(coord[iy],coord_save);
                          v(x)=v(x-pp[iy]) + neumann[ix][iy][1].Eval()*dx[iy];
                          std::swap(coord[iy],coord_save);
                        }
                    }
                }
            }
          /* In 3D, fix the corners.  I am not sure that this is
             needed. */
          if(dim==3)
            {
              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  const int iz((iy+1)%dim);
                  for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
                    {
                      SAMRAI::pdat::SideIndex x(*si);
                      for(int d=0;d<dim;++d)
                        coord[d]=geom->getXLower()[d]
                          + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

                      if(((x[ix]<pbox.lower(ix)
                           && geom->getTouchesRegularBoundary(ix,0)
                           && !is_dirichlet[ix][ix][0])
                          || (x[ix]>pbox.upper(ix)+1
                              && geom->getTouchesRegularBoundary(ix,1)
                              && !is_dirichlet[ix][ix][1]))
                         && ((x[iy]<pbox.lower(iy)
                              && geom->getTouchesRegularBoundary(iy,0))
                             || (x[iy]>pbox.upper(iy) 
                                 && geom->getTouchesRegularBoundary(iy,1))))
                        {
                          if(x[iz]<pbox.lower(iz) 
                             && geom->getTouchesRegularBoundary(iz,0))
                            {
                              double coord_save(geom->getXLower()[iz]);
                              std::swap(coord[iz],coord_save);
                              v(x)=v(x+pp[iz])
                                - neumann[ix][iz][0].Eval()*dx[iz];
                              std::swap(coord[iz],coord_save);
                            }
                          else if(x[iz]>pbox.upper(iz) 
                                  && geom->getTouchesRegularBoundary(iz,1))
                            {
                              double coord_save(geom->getXUpper()[iz]);
                              std::swap(coord[iz],coord_save);
                              v(x)=v(x-pp[iz])
                                + neumann[ix][iz][1].Eval()*dx[iz];
                              std::swap(coord[iz],coord_save);
                            }
                        }
                    }
                }
            }
        }
      catch(mu::Parser::exception_type &e)
        {
          TBOX_ERROR("Error in input formula\n"
                     << "Message:  " << e.GetMsg() << "\n"
                     << "Formula:  " << e.GetExpr() << "\n"
                     << "Token:    " << e.GetToken() << "\n"
                     << "Position: " << e.GetPos() << "\n"
                     << "Errc:     " << e.GetCode() << "\n");
        }
    }
}
