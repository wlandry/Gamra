#ifndef GAMR_DRM_DV
#define GAMR_DRM_DV

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/EdgeData.h"

/* The derivative of the momentum equation w/respect to velocity. It
   is written from the perspective of vx(center_x), but pass in
   different values for center etc. to get vy or vx(!center_x). */

template<class E_data, class E_index>
double dRm_dv_2D(SAMRAI::pdat::CellData<double> &cell_moduli,
                 E_data &edge_moduli,
                 const SAMRAI::pdat::CellIndex &center,
                 const SAMRAI::pdat::CellIndex &left,
                 const E_index &up_e,
                 const E_index &center_e,
                 const double &dx,
                 const double &dy)
{
  return -2*(cell_moduli(center) + cell_moduli(left))/(dx*dx)
    - (edge_moduli(up_e) + edge_moduli(center_e))/(dy*dy);
}

#endif
