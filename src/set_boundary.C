#include "set_boundary.h"
#include "Constants.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

void set_boundary(const SAMRAI::hier::Patch& patch, const int &p_id,
                    const int &v_id, const bool &rhs)
{
  hier::Box pbox=patch.getBox();

  tbox::Pointer<geom::CartesianPatchGeometry> geom = patch.getPatchGeometry();
  const tbox::Dimension Dim(patch.getDim());
  const int dim(patch.getDim().getValue());

  const hier::Index zero(hier::Index::getZeroIndex(Dim));
  hier::Index pp[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    pp[i][i]=1;
  /* This should really get read from the input file. */
  double lower_boundary[]={0,0,0};
  bool lower_dirichlet[]={true,true,true};

  // double upper_boundary[]={-6.94444444444e4,0,0};
  // bool upper_dirichlet[]={true,false,true};

  // bool upper_dirichlet[]={true,true,true};
  // double upper_boundary[]={-1,1,0};

  bool upper_dirichlet[]={true,true,true};
  double upper_boundary[]={0,0,0};

  if(p_id!=invalid_id)
    {
      tbox::Pointer<pdat::CellData<double> > p_ptr = patch.getPatchData(p_id);
      pdat::CellData<double> &p(*p_ptr);

      hier::Box gbox=p.getGhostBox();
      
      for(pdat::CellIterator ci(gbox); ci; ci++)
        {
          pdat::CellIndex center(*ci);
          const int ix(0), iy(1), iz(2);
          hier::Index ip(zero), jp(zero), kp(zero);

          /* Check if we are at a boundary.  If it is a Neumann
           * boundary and not on a corner, then use Neumann
           * conditions.  Otherwise use extrapolation */

          bool neumann_boundary(false), at_boundary(false);
          if(center[ix]<pbox.lower(ix)
             && geom->getTouchesRegularBoundary(ix,0))
            {
              ip[ix]=1;
              if(!lower_dirichlet[ix])
                neumann_boundary=true;
              else
                at_boundary=true;
            }
          else if(center[ix]>pbox.upper(ix)
                  && geom->getTouchesRegularBoundary(ix,1))
            {
              ip[ix]=-1;
              if(!upper_dirichlet[ix])
                neumann_boundary=true;
              else
                at_boundary=true;
            }

          if(center[iy]<pbox.lower(iy)
             && geom->getTouchesRegularBoundary(iy,0))
            {
              jp[iy]=1;
              if(!lower_dirichlet[iy])
                neumann_boundary=true;
              else
                at_boundary=true;
            }
          else if(center[iy]>pbox.upper(iy)
                  && geom->getTouchesRegularBoundary(iy,1))
            {
              jp[iy]=-1;
              if(!upper_dirichlet[iy])
                neumann_boundary=true;
              else
                at_boundary=true;
            }

          if(dim>2)
            {
              if(center[iz]<pbox.lower(iz)
                 && geom->getTouchesRegularBoundary(iz,0))
                {
                  kp[iz]=1;
                  if(!lower_dirichlet[iz])
                    neumann_boundary=true;
                  else
                    at_boundary=true;
                }
              else if(center[iz]>pbox.upper(iz)
                      && geom->getTouchesRegularBoundary(iz,1))
                {
                  kp[iz]=-1;
                  if(!upper_dirichlet[iz])
                    neumann_boundary=true;
                  else
                    at_boundary=true;
                }
            }

          /* Actually apply the BC */
          if(neumann_boundary && !at_boundary)
            {
              p(center)=p(center+ip+jp+kp);
            }
          else if(at_boundary)
            {
              p(center)=2*p(center+ip+jp+kp)-p(center+(ip+jp+kp)*2);
            }
        }
    }


  if(v_id!=invalid_id)
    {
      tbox::Pointer<pdat::SideData<double> > v_ptr = patch.getPatchData(v_id);
      pdat::SideData<double> &v(*v_ptr);

      hier::Box gbox=v.getGhostBox();
      for(int ix=0; ix<dim; ++ix)
        {
          for(pdat::SideIterator si(gbox,ix); si; si++)
            {
              pdat::SideIndex x(*si);

              /* Set a sentinel value for normal components */
              if(x[ix]<pbox.lower(ix) && geom->getTouchesRegularBoundary(ix,0))
                {
                  if(lower_dirichlet[ix])
                    v(x)=boundary_value;
                  else
                    v(x)=v(x+pp[ix]*2);
                }
              else if(x[ix]>pbox.upper(ix)+1
                      && geom->getTouchesRegularBoundary(ix,1))
                {
                  if(upper_dirichlet[ix])
                    v(x)=boundary_value;
                  else
                    v(x)=v(x-pp[ix]*2);
                }
              /* Set values for normal components */
              else if(x[ix]==pbox.lower(ix)
                      && geom->getTouchesRegularBoundary(ix,0)
                      && !rhs && lower_dirichlet[ix])
                {
                  v(x)=lower_boundary[ix];
                }
              else if(x[ix]==pbox.upper(ix)+1
                      && geom->getTouchesRegularBoundary(ix,1)
                      && !rhs && upper_dirichlet[ix])
                {
                  v(x)=upper_boundary[ix];
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
                          v(x)=v(x+pp[iy]);
                        }
                      else if(x[iy]>pbox.upper(iy)
                              && geom->getTouchesRegularBoundary(iy,1))
                        {
                          v(x)=v(x-pp[iy]);
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

          for(pdat::SideIterator si(gbox,ix); si; si++)
            {
              pdat::SideIndex x(*si);
              if((x[ix]<pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0)
                  && !lower_dirichlet[ix])
                 || (x[ix]>pbox.upper(ix)+1
                     && geom->getTouchesRegularBoundary(ix,1)
                     && !upper_dirichlet[ix]))
                {
                  for(int iy=(ix+1)%dim; iy!=ix; iy=(iy+1)%dim)
                    {
                      if(x[iy]<pbox.lower(iy) 
                         && geom->getTouchesRegularBoundary(iy,0))
                        {
                          v(x)=v(x+pp[iy]);
                        }
                      else if(x[iy]>pbox.upper(iy) 
                              && geom->getTouchesRegularBoundary(iy,1))
                        {
                          v(x)=v(x-pp[iy]);
                        }
                    }
                }
            }
          /* In 3D, fix the corners.  I am not sure that this is
             needed. */
          if(dim==3)
            {
              const int iy((ix+1)%dim), iz((iy+1)%dim);
              for(pdat::SideIterator si(gbox,ix); si; si++)
                {
                  pdat::SideIndex x(*si);
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
}
