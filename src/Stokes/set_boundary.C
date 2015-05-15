#include "Stokes/set_boundary.h"
#include "Constants.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

void Stokes_set_boundary(const SAMRAI::hier::Patch& patch, const int &p_id,
                         const int &v_id, const bool &rhs)
{
  SAMRAI::hier::Box pbox=patch.getBox();

  boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
    (patch.getPatchGeometry());
  const SAMRAI::tbox::Dimension Dim(patch.getDim());
  const Gamra::Dir dim(patch.getDim().getValue());

  const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(Dim));
  SAMRAI::hier::Index pp[]={zero,zero,zero};
  for(int i=0;i<dim;++i)
    pp[i][i]=1;
  /* This should really get read from the input file. */
  double lower_boundary[]={0,0,0};
  bool lower_dirichlet[]={true,true,true};

  // double upper_boundary[]={-6.94444444444e4,0,0};
  // bool upper_dirichlet[]={true,false,true};

  // bool upper_dirichlet[]={true,true,true};
  // double upper_boundary[]={-1,1,0};

  double p_lower[]={0,0,0};
  double p_upper[]={0,0,0};

  bool upper_dirichlet[]={true,true,true};
  double upper_boundary[]={0,0,0};

  if(p_id!=invalid_id)
    {
      boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
        (patch.getPatchData(p_id));
      SAMRAI::pdat::CellData<double> &p(*p_ptr);

      SAMRAI::hier::Box gbox=p.getGhostBox();

      SAMRAI::pdat::CellIterator
        cend(SAMRAI::pdat::CellGeometry::end(gbox));      
      for(SAMRAI::pdat::CellIterator
            ci(SAMRAI::pdat::CellGeometry::begin(gbox)); ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex &center(*ci);
          SAMRAI::hier::Index ip(zero), jp(zero), kp(zero);
          SAMRAI::hier::Index pp[]={ip,jp,kp};

          /* Check if we are at a boundary.  If it is a Neumann
           * boundary and not on a corner, then use Neumann
           * conditions.  Otherwise use extrapolation */

          bool dirichlet_boundary(false);

          for(Gamra::Dir ix=0;ix<dim;++ix)
            {
              if(center[ix]<pbox.lower(ix)
                 && geom->getTouchesRegularBoundary(ix,0))
                {
                  pp[ix][ix]=1;
                  if(!lower_dirichlet[ix])
                    {
                      p(center)=-p(center+pp[ix]);
                      if(rhs)
                        p(center)+=2*p_lower[ix];
                    }
                  else
                    dirichlet_boundary=true;
                }
              else if(center[ix]>pbox.upper(ix)
                      && geom->getTouchesRegularBoundary(ix,1))
                {
                  pp[ix][ix]=-1;
                  if(!upper_dirichlet[ix])
                    {
                      p(center)=-p(center+pp[ix]);
                      if(rhs)
                        p(center)+=2*p_upper[ix];
                    }
                  else
                    dirichlet_boundary=true;
                }
            }

          /* Actually apply the BC for Dirichlet velocities. */
          if(dirichlet_boundary)
            {
              p(center)=
                2*p(center+pp[0]+pp[1]+pp[2])-p(center+(pp[0]+pp[1]+pp[2])*2);
            }
        }
    }


  if(v_id!=invalid_id)
    {
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(v_id));
      SAMRAI::pdat::SideData<double> &v(*v_ptr);

      SAMRAI::hier::Box gbox=v.getGhostBox();
      for(Gamra::Dir ix=0; ix<dim; ++ix)
        {
          SAMRAI::pdat::SideIterator
            send(SAMRAI::pdat::SideGeometry::end(gbox,ix));
          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(gbox,ix)); si!=send; ++si)
            {
              const SAMRAI::pdat::SideIndex &x(*si);

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
                  for(Gamra::Dir iy=ix.next(dim); iy!=ix; iy=iy.next(dim))
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

          for(SAMRAI::pdat::SideIterator
                si(SAMRAI::pdat::SideGeometry::begin(gbox,ix)); si!=send; ++si)
            {
              const SAMRAI::pdat::SideIndex &x(*si);
              if((x[ix]<pbox.lower(ix)
                  && geom->getTouchesRegularBoundary(ix,0)
                  && !lower_dirichlet[ix])
                 || (x[ix]>pbox.upper(ix)+1
                     && geom->getTouchesRegularBoundary(ix,1)
                     && !upper_dirichlet[ix]))
                {
                  for(Gamra::Dir iy=ix.next(dim); iy!=ix; iy=iy.next(dim))
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
              const Gamra::Dir iy(ix.next(dim));
              const Gamra::Dir iz(iy.next(dim));
              SAMRAI::pdat::SideIterator
                send(SAMRAI::pdat::SideGeometry::end(gbox,ix));
              for(SAMRAI::pdat::SideIterator
                    si(SAMRAI::pdat::SideGeometry::begin(gbox,ix));
                  si!=send; ++si)
                {
                  const SAMRAI::pdat::SideIndex &x(*si);
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
