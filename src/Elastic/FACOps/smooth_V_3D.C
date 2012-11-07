#include "Elastic/FACOps.h"
#include "Constants.h"
#include "Elastic/dRm_dv.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void Elastic::FACOps::smooth_V_3D
(const int &ix,
 const SAMRAI::hier::Box &pbox,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::CellData<double> &cell_moduli,
 SAMRAI::pdat::EdgeData<double> &edge_moduli,
 const SAMRAI::pdat::CellIndex &center,
 const double Dx[3],
 const double &theta_momentum,
 const SAMRAI::hier::Index pp[3],
 double &maxres)
{
  const int iy((ix+1)%3), iz((ix+2)%3);
  const SAMRAI::pdat::SideIndex x(center,ix,SAMRAI::pdat::SideIndex::Lower),
    y(center,iy,SAMRAI::pdat::SideIndex::Lower),
    z(center,iz,SAMRAI::pdat::SideIndex::Lower);
  const SAMRAI::pdat::EdgeIndex
    edge_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
    edge_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);
    
  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((center[ix]==pbox.lower(ix) && v(x-pp[ix])==boundary_value)
       || (center[ix]==pbox.upper(ix)+1 && v(x+pp[ix])==boundary_value)))
    {
      double C_vx=dRm_dv_3D(cell_moduli,edge_moduli,center,center-pp[ix],
                            edge_y+pp[iz],edge_y,edge_z+pp[iy],edge_z,
                            Dx[ix],Dx[iy],Dx[iz]);

      double delta_Rx=v_rhs(x)
        - v_operator_3D(v,cell_moduli,edge_moduli,center,edge_y,edge_z,
                        x,y,z,pp[ix],pp[iy],pp[iz],Dx[ix],Dx[iy],Dx[iz]);

      /* No scaling here, though there should be. */
      maxres=std::max(maxres,std::fabs(delta_Rx));

      v(x)+=delta_Rx*theta_momentum/C_vx;
    }
}

