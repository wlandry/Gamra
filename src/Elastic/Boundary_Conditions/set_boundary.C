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
  bool lower_dirichlet[]={true,true,true};
  bool upper_dirichlet[]={true,true,true};

  const double *dx=geom->getDx();

  SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
    v_ptr = patch.getPatchData(v_id);
  SAMRAI::pdat::SideData<double> &v(*v_ptr);

  SAMRAI::hier::Box gbox=v.getGhostBox();
  for(int ix=0; ix<dim; ++ix)
    {
      double offset[]={0.5,0.5,0.5};
      offset[ix]=0;
      for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
        {
          SAMRAI::pdat::SideIndex x(*si);

          for(int d=0;d<dim;++d)
            coord[d]=geom->getXLower()[d]
              + dx[d]*(x[d]-pbox.lower()[d]+offset[d]);

          /* Set a sentinel value for normal components */
          if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
            {
              if(is_dirichlet[ix][0][ix])
                v(x)=boundary_value;
              else
                v(x)=v(x+pp[ix]*2);
            }
          else if(x[ix]>pbox.upper(ix)+1
                  && geom->getTouchesRegularBoundary(ix,1))
            {
              if(is_dirichlet[ix][1][ix])
                v(x)=boundary_value;
              else
                v(x)=v(x-pp[ix]*2);
            }
          /* Set values for normal components */
          else if(x[ix]==pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0)
                  && !rhs && is_dirichlet[ix][0][ix])
            {
              v(x)=dirichlet[ix][0][ix].Eval();
            }
          else if(x[ix]==pbox.upper(ix)+1
                  && geom->getTouchesRegularBoundary(ix,1)
                  && !rhs && is_dirichlet[ix][1][ix])
            {
              v(x)=dirichlet[ix][1][ix].Eval();
            }
          /* Set derivatives for tangential component.  The edges
             and corners are incorrect for now.  */
          else
            {
              for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                {
                  if(x[iy]<pbox.lower(iy)
                     && geom->getTouchesRegularBoundary(iy,0))
                    {
                      if(rhs)
                        {
                          v(x)=v(x+pp[iy]);
                        }
                      else
                        {
                          coord[iy]=geom->getXLower()[iy];
                          v(x)=v(x+pp[iy])
                            - neumann[ix][0][iy].Eval()*dx[iy];
                        }
                    }
                  else if(x[iy]>pbox.upper(iy)
                          && geom->getTouchesRegularBoundary(iy,1))
                    {
                      if(rhs)
                        {
                          v(x)=v(x-pp[iy]);
                        }
                      else
                        {
                          coord[iy]=geom->getXUpper()[iy];
                          v(x)=v(x-pp[iy])
                            - neumann[ix][1][iy].Eval()*dx[iy];
                        }
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

      // for(pdat::SideIterator si(gbox,ix); si; si++)
      //   {
      //     pdat::SideIndex x(*si);
      //     if((x[ix]<pbox.lower(ix)
      //         && geom->getTouchesRegularBoundary(ix,0)
      //         && !lower_dirichlet[ix])
      //        || (x[ix]>pbox.upper(ix)+1
      //            && geom->getTouchesRegularBoundary(ix,1)
      //            && !upper_dirichlet[ix]))
      //       {
      //         for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
      //           {
      //             if(x[iy]<pbox.lower(iy) 
      //                && geom->getTouchesRegularBoundary(iy,0))
      //               {
      //                 v(x)=v(x+pp[iy]);
      //               }
      //             else if(x[iy]>pbox.upper(iy) 
      //                     && geom->getTouchesRegularBoundary(iy,1))
      //               {
      //                 v(x)=v(x-pp[iy]);
      //               }
      //           }
      //       }
      //   }
      /* In 3D, fix the corners.  I am not sure that this is
         needed. */
      if(dim==3)
        {
          const int iy((ix+1)%dim), iz((iy+1)%dim);
          for(SAMRAI::pdat::SideIterator si(gbox,ix); si; si++)
            {
              SAMRAI::pdat::SideIndex x(*si);
              if(((x[ix]<pbox.lower(ix)
                   && geom->getTouchesRegularBoundary(ix,0)
                   && !lower_dirichlet[ix])
                  || (x[ix]>pbox.upper(ix)+1
                      && geom->getTouchesRegularBoundary(ix,1)
                      && !upper_dirichlet[ix]))
                 && ((x[iy]<pbox.lower(iy)
                      && geom->getTouchesRegularBoundary(iy,0))
                     || (x[iy]>pbox.upper(iy) 
                         && geom->getTouchesRegularBoundary(iy,1))))
                {
                  if(x[iz]<pbox.lower(iz) 
                     && geom->getTouchesRegularBoundary(iz,0))
                    {
                      v(x)=v(x+pp[iz]);
                    }
                  else if(x[iz]>pbox.upper(iz) 
                          && geom->getTouchesRegularBoundary(iz,1))
                    {
                      v(x)=v(x-pp[iz]);
                    }
                }
            }
        }
    }
}
