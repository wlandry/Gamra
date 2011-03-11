#ifndef GAMR_DRC_DP
#define GAMR_DRC_DP

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/SideData.h"
#include "dRm_dv.h"

/* The derivative of the continuity equation with respect to
   pressure.  Note that pressure does not appear in the continuity
   equation, so we use Tackley's method to chain together
   derivatives */

inline double dRc_dp(const SAMRAI::hier::Box &pbox,
                     const SAMRAI::pdat::CellIndex &center,
                     const SAMRAI::pdat::CellIndex &left,
                     const SAMRAI::pdat::CellIndex &right, 
                     const SAMRAI::pdat::CellIndex &down,
                     const SAMRAI::pdat::CellIndex &up,
                     const SAMRAI::pdat::SideIndex &left_x,
                     const SAMRAI::pdat::SideIndex &right_x,
                     const SAMRAI::pdat::SideIndex &down_y,
                     const SAMRAI::pdat::SideIndex &up_y,
                     SAMRAI::pdat::CellData<double> &cell_viscosity,
                     SAMRAI::pdat::NodeData<double> &edge_viscosity,
                     SAMRAI::pdat::SideData<double> &v,
                     const double &dx,
                     const double &dy)
{
  const SAMRAI::pdat::NodeIndex
    center_e(center,SAMRAI::pdat::NodeIndex::LowerLeft),
    up_e(up,SAMRAI::pdat::NodeIndex::LowerLeft),
    right_e(right,SAMRAI::pdat::NodeIndex::LowerLeft);
  const SAMRAI::hier::Index ip(1,0), jp(0,1);
  const double dRm_dp_xp(1/dx), dRm_dp_xm(-1/dx),
    dRm_dp_yp(1/dy), dRm_dp_ym(-1/dy),
    dRc_dvx_p(-1/dx), dRc_dvx_m(1/dx),
    dRc_dvy_p(-1/dy), dRc_dvy_m(1/dy);

  double result(0);

  if(!(center[0]==pbox.lower(0) && v(left_x)==boundary_value))
    result+=dRc_dvx_p * dRm_dp_xp/dRm_dv(cell_viscosity,edge_viscosity,right,
                                         center,up_e+ip,center_e+ip,dx,dy);

  if(!(center[0]==pbox.upper(0)+1 && v(right_x)==boundary_value))
    result+=dRc_dvx_m * dRm_dp_xm/dRm_dv(cell_viscosity,edge_viscosity,
                                         center,left,up_e,center_e,dx,dy);

  if(!(center[1]==pbox.lower(1) && v(down_y)==boundary_value))
    result+=dRc_dvy_p * dRm_dp_yp/dRm_dv(cell_viscosity,edge_viscosity,up,
                                         center,right_e+jp,center_e+jp,dy,dx);

  if(!(center[1]==pbox.upper(1)+1 && v(up_y)==boundary_value))
    result+=dRc_dvy_m * dRm_dp_ym/dRm_dv(cell_viscosity,edge_viscosity,
                                         center,down,right_e,center_e,dy,dx);

  return result;
}


#endif
