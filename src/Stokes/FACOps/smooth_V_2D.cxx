#include "Stokes/FACOps.hxx"
#include "Constants.hxx"
#include "Stokes/dRm_dv.hxx"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void Stokes::FACOps::smooth_V_2D
(const Gamra::Dir &axis,
 const SAMRAI::hier::Box &pbox,
 boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> &geom,
 const SAMRAI::pdat::CellIndex &center,
 const SAMRAI::hier::Index &ip,
 const SAMRAI::hier::Index &jp,
 SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_rhs,
 double &maxres,
 const double &dx,
 const double &dy,
 SAMRAI::pdat::CellData<double> &cell_viscosity,
 SAMRAI::pdat::NodeData<double> &edge_viscosity,
 const double &theta_momentum)
{
  const Gamra::Dir off_axis=axis.next(2);

  const SAMRAI::pdat::SideIndex x(center,axis,SAMRAI::pdat::SideIndex::Lower),
    y(center,off_axis,SAMRAI::pdat::SideIndex::Lower);
  const SAMRAI::pdat::NodeIndex edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);
    
  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((center[axis]==pbox.lower(axis) && v(x-ip)==boundary_value)
       || (center[axis]==pbox.upper(axis)+1 && v(x+ip)==boundary_value)))
    {
      /* If at the boundary, set things up so that the derivative does
         not change. */
      SAMRAI::hier::Index offset(0,0);
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

      double C_vx=Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,center,center-ip,
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

