#include "StokesFACOps.h"

#include IOMANIP_HEADER_FILE

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "StokesHypreSolver.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

#include "Boundary.h"
#include "dRc_dp.h"

void SAMRAI::solv::StokesFACOps::residual_2D
(pdat::CellData<double> &p,
 pdat::SideData<double> &v,
 pdat::CellData<double> &cell_viscosity,
 pdat::CellData<double> &p_rhs,
 pdat::SideData<double> &v_rhs,
 pdat::CellData<double> &p_resid,
 pdat::SideData<double> &v_resid,
 hier::Patch &patch,
 const hier::Box &pbox,
 const geom::CartesianPatchGeometry &geom)
{
  const hier::Index ip(1,0), jp(0,1);

  tbox::Pointer<pdat::NodeData<double> >
    edge_viscosity_ptr = patch.getPatchData(edge_viscosity_id);
  pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];

  for(pdat::CellIterator ci(pbox); ci; ci++)
    {
      pdat::CellIndex center(*ci);
      pdat::CellIndex up(center), down(center), right(center),
        left(center);

      ++up[1];
      --down[1];
      ++right[0];
      --left[0];

      const pdat::SideIndex
        center_x(center,0,pdat::SideIndex::Lower),
        left_x(left,0,pdat::SideIndex::Lower),
        right_x(right,0,pdat::SideIndex::Lower),
        down_x(down,0,pdat::SideIndex::Lower),
        up_x(up,0,pdat::SideIndex::Lower),
        center_y(center,1,pdat::SideIndex::Lower),
        left_y(left,1,pdat::SideIndex::Lower),
        right_y(right,1,pdat::SideIndex::Lower),
        down_y(down,1,pdat::SideIndex::Lower),
        up_y(up,1,pdat::SideIndex::Lower);
      const pdat::NodeIndex
        center_e(center,pdat::NodeIndex::LowerLeft),
        up_e(up,pdat::NodeIndex::LowerLeft),
        right_e(right,pdat::NodeIndex::LowerLeft);

      /* p */
      if(center[0]!=pbox.upper(0) && center[1]!=pbox.upper(1))
        {
          double dvx_dx=(v(right_x) - v(center_x))/dx;
          double dvy_dy=(v(up_y) - v(center_y))/dy;
          p_resid(center)=p_rhs(center) - dvx_dx - dvy_dy;
        }

      /* vx */
      if(center[1]!=pbox.upper(1))
        {
          /* If x==0 */
          if((center[0]==pbox.lower(0) && v(left_x)==boundary_value)
             || (center[0]==pbox.upper(0)
                 && v(right_x)==boundary_value))
                       
            {
              v_resid(center_x)=0;
            }
          else
            {
              v_resid(center_x)=v_rhs(center_x)
                - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                                left,center_x,right_x,left_x,up_x,down_x,
                                center_y,up_y,center_e,up_e,ip,dx,dy);
            }
        }

      /* vy */
      if(center[0]!=pbox.upper(0))
        {
          /* If y==0 */
          if((center[1]==pbox.lower(1) && v(down_y)==boundary_value)
             || (center[1]==pbox.upper(1) && v(up_y)==boundary_value))
            {
              v_resid(center_y)=0;
            }
          else
            {
              v_resid(center_y)=v_rhs(center_y)
                - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                                down,center_y,up_y,down_y,right_y,left_y,
                                center_x,right_x,center_e,right_e,jp,
                                dy,dx);
            }
        }
    }
}

