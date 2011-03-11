#ifndef GAMR_DRM_DV
#define GAMR_DRM_DV

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"

/* The derivative of the momentum equation w/respect to velocity. It
   is written from the perspective of vx(center_x), but pass in
   different values for center etc. to get vy or vx(!center_x). */

inline double dRm_dv(SAMRAI::pdat::CellData<double> &cell_viscosity,
                     SAMRAI::pdat::NodeData<double> &edge_viscosity,
                     const SAMRAI::pdat::CellIndex &center,
                     const SAMRAI::pdat::CellIndex &left,
                     const SAMRAI::pdat::NodeIndex &up_e,
                     const SAMRAI::pdat::NodeIndex &center_e,
                     const double &dx,
                     const double &dy)
{
  return -2*(cell_viscosity(center) + cell_viscosity(left))/(dx*dx)
    - (edge_viscosity(up_e) + edge_viscosity(center_e))/(dy*dy);
}




#endif
