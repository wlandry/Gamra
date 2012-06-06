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
  return -(   (cell_moduli(center,0) + cell_moduli(left,0))
           +2*(cell_moduli(center,1) + cell_moduli(left,1)))/(dx*dx)
         -(edge_moduli(up_e,1) + edge_moduli(center_e,1))/(dy*dy);
}

inline double dRm_dv_3D(SAMRAI::pdat::CellData<double> &cell_moduli,
		        SAMRAI::pdat::EdgeData<double> &edge_moduli,
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
	  return dRm_dv_2D(cell_moduli,edge_moduli,center,left,front_y,center_y,dx,dy)
		      - (edge_moduli(up_z,1) + edge_moduli(center_z,1))/(dz*dz);
}

#endif
