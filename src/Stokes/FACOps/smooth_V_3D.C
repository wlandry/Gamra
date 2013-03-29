#include "Stokes/FACOps.h"
#include "Constants.h"
#include "Stokes/dRm_dv.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void Stokes::FACOps::smooth_V_3D
(const int &ix,
 const SAMRAI::hier::Box &pbox,
 boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> &geom,
 SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::CellData<double> &cell_viscosity,
 SAMRAI::pdat::EdgeData<double> &edge_viscosity,
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
      /* If at the boundary, set things up so that the derivative does
         not change. */
      SAMRAI::hier::Index offset(0,0,0);
      offset[ix]=2;
      bool set_lower_boundary(false), set_upper_boundary(false);
      double dv_lower(0), dv_upper(0);
      if(center[ix]==pbox.lower(ix)+1
         && !geom->getTouchesRegularBoundary(ix,0))
        {
          set_lower_boundary=true;
          dv_lower=v(x-offset) - v(x);
        }
      if(center[ix]==pbox.upper(ix) && !geom->getTouchesRegularBoundary(ix,1))
        {
          set_upper_boundary=true;
          dv_upper=v(x+offset) - v(x);
        }

      double C_vx=Stokes_dRm_dv_3D(cell_viscosity,edge_viscosity,center,center-pp[ix],
                                   edge_y+pp[iz],edge_y,edge_z+pp[iy],edge_z,
                                   Dx[ix],Dx[iy],Dx[iz]);

      double delta_Rx=v_rhs(x)
        - v_operator_3D(v,p,cell_viscosity,edge_viscosity,center,edge_y,edge_z,
                        x,y,z,pp[ix],pp[iy],pp[iz],Dx[ix],Dx[iy],Dx[iz]);

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

