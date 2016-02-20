#include "Stokes/FACOps.hxx"
#include "Constants.hxx"

void Stokes::FACOps::residual_3D
(SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::CellData<double> &cell_viscosity,
 SAMRAI::pdat::CellData<double> &p_rhs,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::CellData<double> &p_resid,
 SAMRAI::pdat::SideData<double> &v_resid,
 SAMRAI::hier::Patch &patch,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::geom::CartesianPatchGeometry &geom)
{
  boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_viscosity_ptr = 
    boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
    (patch.getPatchData(edge_viscosity_id));
  SAMRAI::pdat::EdgeData<double> &edge_viscosity(*edge_viscosity_ptr);

  const double *Dx = geom.getDx();
  const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const SAMRAI::hier::Index pp[]={ip,jp,kp};

  SAMRAI::pdat::CellIterator cend(SAMRAI::pdat::CellGeometry::end(pbox));
  for(SAMRAI::pdat::CellIterator ci(SAMRAI::pdat::CellGeometry::begin(pbox));
      ci!=cend; ++ci)
    {
      const SAMRAI::pdat::CellIndex &center(*ci);
      SAMRAI::pdat::CellIndex up(center), down(center), right(center),
        left(center), front(center), back(center);

      ++right[0];
      --left[0];
      ++up[1];
      --down[1];
      ++front[2];
      --back[2];

      /* p */
      if(center[0]!=pbox.upper(0) && center[1]!=pbox.upper(1)
         && center[2]!=pbox.upper(2))
        {
          const SAMRAI::pdat::SideIndex
            x(center,0,SAMRAI::pdat::SideIndex::Lower),
            y(center,1,SAMRAI::pdat::SideIndex::Lower),
            z(center,2,SAMRAI::pdat::SideIndex::Lower);

          double dvx_dx=(v(x+ip) - v(x))/Dx[0];
          double dvy_dy=(v(y+jp) - v(y))/Dx[1];
          double dvz_dz=(v(z+kp) - v(z))/Dx[2];
          p_resid(center)=p_rhs(center) - dvx_dx - dvy_dy - dvz_dz;
        }

      for(Gamra::Dir ix=0;ix<3;++ix)
        {
          const Gamra::Dir iy(ix.next(3));
          const Gamra::Dir iz(iy.next(3));
          const SAMRAI::pdat::SideIndex
            x(center,ix,SAMRAI::pdat::SideIndex::Lower),
            y(center,iy,SAMRAI::pdat::SideIndex::Lower),
            z(center,iz,SAMRAI::pdat::SideIndex::Lower);
          const SAMRAI::pdat::EdgeIndex
            edge_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
            edge_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);

          if(center[iy]!=pbox.upper(iy) && center[iz]!=pbox.upper(iz))
            {
              if((center[ix]==pbox.lower(ix) && v(x-pp[ix])==boundary_value)
                 || (center[ix]==pbox.upper(ix) && v(x+pp[ix])==boundary_value))
                {
                  v_resid(x)=0;
                }
              else
                {
                  v_resid(x)=v_rhs(x)
                    - v_operator_3D(v,p,cell_viscosity,edge_viscosity,
                                    center,edge_y,edge_z,x,y,z,
                                    pp[ix],pp[iy],pp[iz],Dx[ix],Dx[iy],Dx[iz]);
                }
            }
        }          
    }
}

