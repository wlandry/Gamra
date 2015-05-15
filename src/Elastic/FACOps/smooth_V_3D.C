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
(const Gamra::Dir &ix,
 const SAMRAI::hier::Box &pbox,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::CellData<double> &cell_moduli,
 SAMRAI::pdat::EdgeData<double> &edge_moduli,
 const SAMRAI::pdat::CellIndex &cell,
 const double Dx[3],
 const double &theta_momentum,
 const SAMRAI::hier::Index unit[3],
 double &maxres)
{
  const Gamra::Dir iy(ix.next(3));
  const Gamra::Dir iz(iy.next(3));
  const SAMRAI::pdat::SideIndex x(cell,ix,SAMRAI::pdat::SideIndex::Lower),
    y(cell,iy,SAMRAI::pdat::SideIndex::Lower),
    z(cell,iz,SAMRAI::pdat::SideIndex::Lower);
  const SAMRAI::pdat::EdgeIndex
    edge_y(cell,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
    edge_z(cell,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);
    
  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((cell[ix]==pbox.lower(ix) && v(x-unit[ix])==boundary_value)
       || (cell[ix]==pbox.upper(ix)+1 && v(x+unit[ix])==boundary_value)))
    {
      double C_vx=dRm_dv_3D(cell_moduli,edge_moduli,cell,cell-unit[ix],
                            edge_y+unit[iz],edge_y,edge_z+unit[iy],edge_z,
                            Dx[ix],Dx[iy],Dx[iz]);

      double delta_Rx=v_rhs(x)
        - v_operator_3D(v,cell_moduli,edge_moduli,cell,edge_y,edge_z,
                        x,y,z,unit[ix],unit[iy],unit[iz],Dx[ix],Dx[iy],Dx[iz]);

      /* No scaling here, though there should be. */
      maxres=std::max(maxres,std::fabs(delta_Rx));

      v(x)+=delta_Rx*theta_momentum/C_vx;
    }
}

