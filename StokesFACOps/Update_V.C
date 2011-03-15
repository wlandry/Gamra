#include "StokesFACOps.h"
#include "Boundary.h"
#include "dRc_dp.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void SAMRAI::solv::StokesFACOps::Update_V
(const int &axis,
 const int j,
 const hier::Box &pbox,
 tbox::Pointer<geom::CartesianPatchGeometry> &geom,
 const pdat::CellIndex &center,
 const pdat::CellIndex &left,
 const pdat::CellIndex &right, 
 const pdat::CellIndex &down,
 const pdat::CellIndex &up,
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
  hier::Index ip(0,0), jp(0,0);
  ip[axis]=1;
  jp[off_axis]=1;

  const pdat::SideIndex
    center_x(center,axis,pdat::SideIndex::Lower),
    left_x(left,axis,pdat::SideIndex::Lower),
    right_x(right,axis,pdat::SideIndex::Lower),
    down_x(down,axis,pdat::SideIndex::Lower),
    up_x(up,axis,pdat::SideIndex::Lower),
    center_y(center,off_axis,pdat::SideIndex::Lower),
    up_y(up,off_axis,pdat::SideIndex::Lower);
  const pdat::NodeIndex
    center_e(center,pdat::NodeIndex::LowerLeft),
    up_e(up,pdat::NodeIndex::LowerLeft);
    
  /* Update vx */
  if(j<pbox.upper(off_axis)+1)
    {
      /* If at the 'x' boundaries, leave vx as is */
      if(!((center[axis]==pbox.lower(axis) && v(left_x)==boundary_value)
           || (center[axis]==pbox.upper(axis)+1 && v(right_x)==boundary_value)))
        {
          double C_vx;
          /* If y==0 */
          hier::Index offset(0,0);
          bool set_boundary(false);
          if(center[axis]==pbox.lower(axis)+1
             && !geom->getTouchesRegularBoundary(axis,0))
            {
              offset[axis]=-2;
              set_boundary=true;
            }
          else if(center[axis]==pbox.upper(axis)
                  && !geom->getTouchesRegularBoundary(axis,1))
            {
              offset[axis]=2;
              set_boundary=true;
            }

          double dv(0);
          if(set_boundary)
            {
              dv=v(center_x+offset) - v(center_x);
            }

          C_vx=dRm_dv(cell_viscosity,edge_viscosity,center,left,up_e,center_e,
                      dx,dy);

          double delta_Rx=v_rhs(center_x)
            - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                            center_e,center_x,center_y,ip,jp,dx,dy);

          /* No scaling here, though there should be. */
          maxres=std::max(maxres,std::fabs(delta_Rx));

          v(center_x)+=delta_Rx*theta_momentum/C_vx;

          /* Set the boundary elements so that the
             derivative is zero. */
          if(set_boundary)
            {
              v(center_x+offset)=v(center_x) + dv;
            }
        }
    }
}
