#include "Elastic/FACOps.h"
#include "Constants.h"
#include "Elastic/dRm_dv.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void Elastic::FACOps::smooth_V_2D
(const int &axis,
 const SAMRAI::hier::Box &pbox,
 SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry> &geom,
 const SAMRAI::pdat::CellIndex &center,
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
  const int off_axis=(axis==0) ? 1 : 0;

  const SAMRAI::pdat::SideIndex x(center,axis,SAMRAI::pdat::SideIndex::Lower),
    y(center,off_axis,SAMRAI::pdat::SideIndex::Lower);
  const SAMRAI::pdat::NodeIndex edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((center[axis]==pbox.lower(axis) && v(x-ip)==boundary_value)
       || (center[axis]==pbox.upper(axis)+1 && v(x+ip)==boundary_value)))
    {
      // /* If at the boundary, set things up so that the derivative does
      //    not change. */
      // bool set_lower_boundary(false), set_upper_boundary(false);
      // double dv_lower(0), dv_upper(0);

      // SAMRAI::hier::Index offset(0,0);
      // offset[axis]=2;

      // /* Need to fix this for Neumann conditions at a physical boundary */

      // if(center[off_axis]==pbox.lower(off_axis)+1)
      //   {
      //     set_lower_boundary=true;
      //     dv_lower=v(x-offset) - v(x);
      //   }
      // if(center[off_axis]==pbox.upper(off_axis)-1)
      //   {
      //     set_upper_boundary=true;
      //     dv_upper=v(x+offset) - v(x);
      //   }

      double C_vx=dRm_dv_2D(cell_moduli,edge_moduli,center,center-ip,
                            edge+jp,edge,dx,dy);

      double delta_Rx=v_rhs(x)
        - v_operator_2D(v,cell_moduli,edge_moduli,center,
                        edge,x,y,ip,jp,dx,dy);

      SAMRAI::tbox::plog << "smooth "
                         << axis << " "
                         << x << " "
                         << v(x) << " "
                         << v_rhs(x) << " "
                         << delta_Rx << " "
                         << v(x+ip) << " "
                         << v(x-ip) << " "
                         << v(x+jp) << " "
                         << v(x-jp) << " "
                         << v(y) << " "
                         << v(y+jp) << " "
                         << v(y-ip) << " "
                         << v(y-ip+jp) << " ";

        maxres=std::max(maxres,std::fabs(delta_Rx));
        v(x)+=delta_Rx*theta_momentum/C_vx;

        SAMRAI::tbox::plog << v(x) << " "
                         << "\n";

      // /* Set the boundary elements so that the derivative is
      //    unchanged. */
      // if(set_lower_boundary)
      //   {
      //     v(x-offset)=v(x) + dv_lower;
      //   }
      // if(set_upper_boundary)
      //   {
      //     v(x+offset)=v(x) + dv_upper;
      //   }
    }
}

