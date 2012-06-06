#include "Elastic/FACOps.h"
#include "Constants.h"
#include "Elastic/dRc_dp.h"
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void SAMRAI::solv::Elastic::FACOps::smooth_V_3D
(const int &ix,
 const hier::Box &pbox,
 tbox::Pointer<geom::CartesianPatchGeometry> &geom,
 pdat::SideData<double> &v,
 pdat::SideData<double> &v_rhs,
 pdat::CellData<double> &cell_moduli,
 pdat::EdgeData<double> &edge_moduli,
 const pdat::CellIndex &center,
 const double Dx[3],
 const double &theta_momentum,
 const hier::Index pp[3],
 double &maxres)
{
  const int iy((ix+1)%3), iz((ix+2)%3);
  const pdat::SideIndex x(center,ix,pdat::SideIndex::Lower),
    y(center,iy,pdat::SideIndex::Lower),
    z(center,iz,pdat::SideIndex::Lower);
  const pdat::EdgeIndex edge_y(center,iy,pdat::EdgeIndex::LowerLeft),
    edge_z(center,iz,pdat::EdgeIndex::LowerLeft);
    
  /* If at a Dirichlet 'x' boundary, leave vx as is */
  if(!((center[ix]==pbox.lower(ix) && v(x-pp[ix])==boundary_value)
       || (center[ix]==pbox.upper(ix)+1 && v(x+pp[ix])==boundary_value)))
    {
      /* If at the boundary, set things up so that the derivative does
         not change. */
      hier::Index offset(0,0,0);
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

      double C_vx=dRm_dv_3D(cell_moduli,edge_moduli,center,center-pp[ix],
                            edge_y+pp[iz],edge_y,edge_z+pp[iy],edge_z,
                            Dx[ix],Dx[iy],Dx[iz]);

      double delta_Rx=v_rhs(x)
        - v_operator_3D(v,cell_moduli,edge_moduli,center,edge_y,edge_z,
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

