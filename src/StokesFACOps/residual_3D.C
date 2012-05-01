#include "StokesFACOps.h"
#include "Constants.h"

void SAMRAI::solv::StokesFACOps::residual_3D
(pdat::CellData<double> &p,
 pdat::SideData<double> &v,
 pdat::CellData<double> &cell_moduli,
 pdat::CellData<double> &p_rhs,
 pdat::SideData<double> &v_rhs,
 pdat::CellData<double> &p_resid,
 pdat::SideData<double> &v_resid,
 hier::Patch &patch,
 const hier::Box &pbox,
 const geom::CartesianPatchGeometry &geom)
{
  tbox::Pointer<pdat::EdgeData<double> >
    edge_moduli_ptr = patch.getPatchData(edge_moduli_id);
  pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);

  const double *Dx = geom.getDx();
  const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const hier::Index pp[]={ip,jp,kp};

  for(pdat::CellIterator ci(pbox); ci; ci++)
    {
      pdat::CellIndex center(*ci);
      pdat::CellIndex up(center), down(center), right(center),
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
          p_resid(center)=0;
        }

      for(int ix=0;ix<3;++ix)
        {
          const int iy((ix+1)%3), iz((ix+2)%3);
          const pdat::SideIndex
            x(center,ix,pdat::SideIndex::Lower),
            y(center,iy,pdat::SideIndex::Lower),
            z(center,iz,pdat::SideIndex::Lower);
          const pdat::EdgeIndex
            edge_y(center,iy,pdat::EdgeIndex::LowerLeft),
            edge_z(center,iz,pdat::EdgeIndex::LowerLeft);

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
                    - v_operator_3D(v,cell_moduli,edge_moduli,
                                    center,edge_y,edge_z,x,y,z,
                                    pp[ix],pp[iy],pp[iz],Dx[ix],Dx[iy],Dx[iz]);
                }
            }
        }          
    }
}

