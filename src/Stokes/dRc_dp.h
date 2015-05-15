#ifndef GAMR_STOKES_DRC_DP_H
#define GAMR_STOKES_DRC_DP_H

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/SideData.h"
#include "Stokes/dRm_dv.h"
#include "Constants.h"

/* These two functions should really have a more similar API */

/* The derivative of the continuity equation with respect to
   pressure.  Note that pressure does not appear in the continuity
   equation, so we use Tackley's method to chain together
   derivatives */

inline double Stokes_dRc_dp_2D(const SAMRAI::hier::Box &pbox,
                               const SAMRAI::pdat::CellIndex &center,
                               const SAMRAI::pdat::SideIndex &x,
                               const SAMRAI::pdat::SideIndex &y,
                               SAMRAI::pdat::CellData<double> &cell_viscosity,
                               SAMRAI::pdat::NodeData<double> &edge_viscosity,
                               SAMRAI::pdat::SideData<double> &v,
                               const double &dx,
                               const double &dy)
{
  const SAMRAI::hier::Index ip(1,0), jp(0,1);
  const SAMRAI::pdat::NodeIndex
    center_e(center,SAMRAI::pdat::NodeIndex::LowerLeft),
    up_e(center_e+jp),right_e(center_e+ip);
    
  const double dRm_dp_xp(1/dx), dRm_dp_xm(-1/dx),
    dRm_dp_yp(1/dy), dRm_dp_ym(-1/dy),
    dRc_dvx_p(-1/dx), dRc_dvx_m(1/dx),
    dRc_dvy_p(-1/dy), dRc_dvy_m(1/dy);

  double result(0);

  if(!(center[0]==pbox.lower(0) && v(x-ip)==boundary_value))
    result+=dRc_dvx_m * dRm_dp_xm/Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,
                                                   center,center-ip,up_e,center_e,
                                                   dx,dy);

  if(!(center[0]==pbox.upper(0) && v(x+ip*2)==boundary_value))
    result+=dRc_dvx_p * dRm_dp_xp/Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,
                                                   center+ip,center,up_e+ip,
                                                   center_e+ip,dx,dy);
  if(!(center[1]==pbox.lower(1) && v(y-jp)==boundary_value))
    result+=dRc_dvy_m * dRm_dp_ym/Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,
                                                   center,center-jp,right_e,center_e,
                                                   dy,dx);

  if(!(center[1]==pbox.upper(1) && v(y+jp*2)==boundary_value))
    result+=dRc_dvy_p * dRm_dp_yp/Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,
                                                   center+jp,center,right_e+jp,
                                                   center_e+jp,dy,dx);
                                            
  return result;
}


inline double Stokes_dRc_dp_3D(const SAMRAI::hier::Box &pbox,
                               const SAMRAI::pdat::CellIndex &center,
                               SAMRAI::pdat::CellData<double> &cell_viscosity,
                               SAMRAI::pdat::EdgeData<double> &edge_viscosity,
                               SAMRAI::pdat::SideData<double> &v,
                               const double Dx[3],
                               const SAMRAI::hier::Index pp[3])
{
  double result(0);
  for(Gamra::Dir ix=0;ix<3;++ix)
    {
      const Gamra::Dir iy(ix.next(3));
      const Gamra::Dir iz(iy.next(3));
      const SAMRAI::pdat::SideIndex x(center,ix,SAMRAI::pdat::SideIndex::Lower);
      const SAMRAI::pdat::EdgeIndex
        center_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
        front_y(center_y+pp[iz]),
        center_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft),
        up_z(center_z+pp[iy]);

      if(!(center[ix]==pbox.lower(ix) && v(x-pp[ix])==boundary_value))
        result+=-1.0/(Stokes_dRm_dv_3D(cell_viscosity,edge_viscosity,
                                       center,center-pp[ix],front_y,center_y,
                                       up_z,center_z,Dx[ix],Dx[iy],Dx[iz])
                      * Dx[ix]*Dx[ix]);

      if(!(center[ix]==pbox.upper(ix) && v(x+pp[ix]*2)==boundary_value))
        result+=-1.0/(Stokes_dRm_dv_3D(cell_viscosity,edge_viscosity,
                                       center+pp[ix],center,front_y+pp[ix],center_y+pp[ix],
                                       up_z+pp[ix],center_z+pp[ix],
                                       Dx[ix],Dx[iy],Dx[iz])
                      * Dx[ix]*Dx[ix]);
    }
  return result;
}
#endif
