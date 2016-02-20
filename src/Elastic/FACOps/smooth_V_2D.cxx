#include "Elastic/FACOps.hxx"
#include "Constants.hxx"
#include "Elastic/dRm_dv.hxx"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void Elastic::FACOps::smooth_V_2D
(const Gamra::Dir &axis,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::pdat::CellIndex &cell,
 const SAMRAI::hier::Index &ip,
 const SAMRAI::hier::Index &jp,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_rhs,
 double &maxres,
 const double &dx,
 const double &dy,
 SAMRAI::pdat::CellData<double> &cell_moduli,
 SAMRAI::pdat::NodeData<double> &edge_moduli,
 const double &theta_momentum)
{
  const Gamra::Dir off_axis=(axis==0) ? 1 : 0;

  const SAMRAI::pdat::SideIndex x(cell,axis,SAMRAI::pdat::SideIndex::Lower),
    y(cell,off_axis,SAMRAI::pdat::SideIndex::Lower);
  const SAMRAI::pdat::NodeIndex edge(cell,SAMRAI::pdat::NodeIndex::LowerLeft);

  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((cell[axis]==pbox.lower(axis) && v(x-ip)==boundary_value)
       || (cell[axis]==pbox.upper(axis)+1 && v(x+ip)==boundary_value)))
    {
      double C_vx=dRm_dv_2D(cell_moduli,edge_moduli,cell,cell-ip,
                            edge+jp,edge,dx,dy);

      double delta_Rx=v_rhs(x)
        - v_operator_2D(v,cell_moduli,edge_moduli,cell,
                        edge,x,y,ip,jp,dx,dy);

      maxres=std::max(maxres,std::fabs(delta_Rx));
      v(x)+=delta_Rx*theta_momentum/C_vx;
    }
}

