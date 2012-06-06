#ifndef GAMR_DRC_DP
#define GAMR_DRC_DP

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/SideData.h"
#include "dRm_dv.h"
#include "Constants.h"

/* These two functions should really have a more similar API */

/* The derivative of the continuity equation with respect to
   pressure.  Note that pressure does not appear in the continuity
   equation, so we use Tackley's method to chain together
   derivatives */

inline double dRc_dp_2D(const SAMRAI::hier::Box &pbox,
                        const SAMRAI::pdat::CellIndex &center,
                        const SAMRAI::pdat::SideIndex &x,
                        const SAMRAI::pdat::SideIndex &y,
                        SAMRAI::pdat::CellData<double> &cell_moduli,
                        SAMRAI::pdat::NodeData<double> &edge_moduli,
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
    result+=dRc_dvx_m * dRm_dp_xm/dRm_dv_2D(cell_moduli,edge_moduli,
                                            center,center-ip,up_e,center_e,
                                            dx,dy);

  if(!(center[0]==pbox.upper(0) && v(x+ip*2)==boundary_value))
    result+=dRc_dvx_p * dRm_dp_xp/dRm_dv_2D(cell_moduli,edge_moduli,
                                            center+ip,center,up_e+ip,
                                            center_e+ip,dx,dy);
  if(!(center[1]==pbox.lower(1) && v(y-jp)==boundary_value))
    result+=dRc_dvy_m * dRm_dp_ym/dRm_dv_2D(cell_moduli,edge_moduli,
                                            center,center-jp,right_e,center_e,
                                            dy,dx);

  if(!(center[1]==pbox.upper(1) && v(y+jp*2)==boundary_value))
    result+=dRc_dvy_p * dRm_dp_yp/dRm_dv_2D(cell_moduli,edge_moduli,
                                            center+jp,center,right_e+jp,
                                            center_e+jp,dy,dx);
                                            
  return result;
}

#endif
