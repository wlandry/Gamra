#include "StokesFACOps.h"
#include "Constants.h"
#include "dRm_dv.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void SAMRAI::solv::StokesFACOps::smooth_V_2D
(const int &axis,
 const hier::Box &pbox,
 tbox::Pointer<geom::CartesianPatchGeometry> &geom,
 const pdat::CellIndex &center,
 const hier::Index &ip,
 const hier::Index &jp,
 pdat::CellData<double> &p,
 pdat::SideData<double> &v,
 pdat::SideData<double> &v_rhs,
 double &maxres,
 const double &dx,
 const double &dy,
 pdat::CellData<double> &cell_viscosity,
 pdat::NodeData<double> &edge_viscosity,
 const double &theta_momentum)
{
  const int off_axis=(axis==0) ? 1 : 0;

  const pdat::SideIndex x(center,axis,pdat::SideIndex::Lower),
    y(center,off_axis,pdat::SideIndex::Lower);
  const pdat::NodeIndex edge(center,pdat::NodeIndex::LowerLeft);
    
  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((center[axis]==pbox.lower(axis) && v(x-ip)==boundary_value)
       || (center[axis]==pbox.upper(axis)+1 && v(x+ip)==boundary_value)))
    {
      /* If at the boundary, set things up so that the derivative does
         not change. */
      hier::Index offset(0,0);
      offset[axis]=2;
      bool set_lower_boundary(false), set_upper_boundary(false);
      double dv_lower(0), dv_upper(0);
      if(center[axis]==pbox.lower(axis)+1
         && !geom->getTouchesRegularBoundary(axis,0))
        {
          set_lower_boundary=true;
          dv_lower=v(x-offset) - v(x);
        }
      if(center[axis]==pbox.upper(axis)
         && !geom->getTouchesRegularBoundary(axis,1))
        {
          set_upper_boundary=true;
          dv_upper=v(x+offset) - v(x);
        }

      double C_vx=dRm_dv_2D(cell_viscosity,edge_viscosity,center,center-ip,
                            edge+jp,edge,dx,dy);

      double delta_Rx=v_rhs(x)
        - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                        edge,x,y,ip,jp,dx,dy);

      /* No scaling here, though there should be. */
      maxres=std::max(maxres,std::fabs(delta_Rx));

      v(x)+=delta_Rx*theta_momentum/C_vx;

      /* Set the boundary elements so that the derivative is
         unchanged. */
      if(set_lower_boundary)
        {
          v(x-offset)=v(x) + dv_lower;
        }
      if(set_upper_boundary)
        {
          v(x+offset)=v(x) + dv_upper;
        }
    }
}

