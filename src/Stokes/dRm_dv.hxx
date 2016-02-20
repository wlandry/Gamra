#ifndef GAMR_STOKES_DRM_DV_H
#define GAMR_STOKES_DRM_DV_H

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/EdgeData.h>

/* The derivative of the momentum equation w/respect to velocity. It
   is written from the perspective of vx(center_x), but pass in
   different values for center etc. to get vy or vx(!center_x). */

template<class E_data, class E_index>
double Stokes_dRm_dv_2D(SAMRAI::pdat::CellData<double> &cell_viscosity,
                        E_data &edge_viscosity,
                        const SAMRAI::pdat::CellIndex &center,
                        const SAMRAI::pdat::CellIndex &left,
                        const E_index &up_e,
                        const E_index &center_e,
                        const double &dx,
                        const double &dy)
{
  return -2*(cell_viscosity(center) + cell_viscosity(left))/(dx*dx)
    - (edge_viscosity(up_e) + edge_viscosity(center_e))/(dy*dy);
}

inline double Stokes_dRm_dv_3D(SAMRAI::pdat::CellData<double> &cell_viscosity,
                               SAMRAI::pdat::EdgeData<double> &edge_viscosity,
                               const SAMRAI::pdat::CellIndex &center,
                               const SAMRAI::pdat::CellIndex &left,
                               const SAMRAI::pdat::EdgeIndex &front_y,
                               const SAMRAI::pdat::EdgeIndex &center_y,
                               const SAMRAI::pdat::EdgeIndex &up_z,
                               const SAMRAI::pdat::EdgeIndex &center_z,
                               const double &dx,
                               const double &dy,
                               const double &dz)
{
  return Stokes_dRm_dv_2D(cell_viscosity,edge_viscosity,center,left,front_y,center_y,dx,dy)
    - (edge_viscosity(up_z) + edge_viscosity(center_z))/(dz*dz);
}




#endif
