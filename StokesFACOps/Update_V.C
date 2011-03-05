/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Stokes using FAC 
 *
 ************************************************************************/
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
/*
********************************************************************
* Updates one component of the velocity during a red-black *
* Gauss-Seidel iteration.  *
********************************************************************
*/
void SAMRAI::solv::StokesFACOps::Update_V
(const int &axis,
 const int j,
 const hier::Box &pbox,
 tbox::Pointer<geom::CartesianPatchGeometry> &geom,
 const pdat::CellIndex &center,
 const pdat::CellIndex &left,
 const pdat::CellIndex &right, 
 const pdat::CellIndex &down,
 const pdat::CellIndex &up,
 pdat::CellData<double> &p,
 pdat::SideData<double> &v,
 pdat::SideData<double> &v_rhs,
 double &maxres,
 const double &dx,
 const double &dy,
 pdat::CellData<double> &cell_viscosity,
 pdat::NodeData<double> &edge_viscosity,
 const double &theta_momentum)
{
  const int off_axis=(axis==0) ? 1 : 0;
  hier::Index ip(0,0);
  ip[axis]=1;

  const pdat::SideIndex
    center_x(center,axis,pdat::SideIndex::Lower),
    left_x(left,axis,pdat::SideIndex::Lower),
    right_x(right,axis,pdat::SideIndex::Lower),
    down_x(down,axis,pdat::SideIndex::Lower),
    up_x(up,axis,pdat::SideIndex::Lower),
    center_y(center,off_axis,pdat::SideIndex::Lower),
    up_y(up,off_axis,pdat::SideIndex::Lower);
  const pdat::NodeIndex
    center_e(center,pdat::NodeIndex::LowerLeft),
    up_e(up,pdat::NodeIndex::LowerLeft);
    
  /* Update vx */
  if(j<pbox.upper(off_axis)+1)
    {
      /* If at the 'x' boundaries, leave vx as is */
      if(!((center[axis]==pbox.lower(axis) && v(left_x)==boundary_value)
           || (center[axis]==pbox.upper(axis)+1 && v(right_x)==boundary_value)))
        {
          double d2vx_dxx, d2vx_dyy, C_vx;
          /* If y==0 */
          hier::Index offset(0,0);
          bool set_boundary(false);
          if(center[axis]==pbox.lower(axis)+1
             && !geom->getTouchesRegularBoundary(axis,0))
            {
              offset[axis]=-2;
              set_boundary=true;
            }
          else if(center[axis]==pbox.upper(axis)
                  && !geom->getTouchesRegularBoundary(axis,1))
            {
              offset[axis]=2;
              set_boundary=true;
            }

          double dv(0);
          if(set_boundary)
            {
              dv=v(center_x+offset) - v(center_x);
            }

          C_vx=dRm_dv(cell_viscosity,edge_viscosity,center,left,up_e,center_e,
                      dx,dy);

          double delta_Rx=v_rhs(center_x)
            - v_operator(v,p,cell_viscosity,edge_viscosity,center,left,center_x,
                         right_x,left_x,up_x,down_x,center_y,up_y,center_e,up_e,
                         ip,dx,dy);

          /* No scaling here, though there should be. */
          maxres=std::max(maxres,std::fabs(delta_Rx));


          {
    const double dtau_xx_dx =
      2*((v(right_x)-v(center_x))*cell_viscosity(center)
         - (v(center_x)-v(left_x))*cell_viscosity(left))/(dx*dx);

    const double dtau_xy_dy = 
      edge_viscosity(up_e)*((v(up_x)-v(center_x))/(dy*dy)
                            + (v(up_y)-v(up_y-ip))/(dx*dy))
      - edge_viscosity(center_e)*((v(center_x)-v(down_x))/(dy*dy)
                                  + (v(center_y)-v(center_y-ip))/(dx*dy));
    const double dp_dx=(p(center)-p(left))/dx;

          // tbox::plog << "v " << axis << " "
          //            << center[0] << " "
          //            << center[1] << " "
          //            << v(center_x) << " "
          //            << v(left_x) << " "
          //            << v(right_x) << " "
          //            << v(up_x) << " "
          //            << v(down_x) << " "
          //            << v(up_y) << " "
          //            << v(up_y-ip) << " "
          //            << v(center_y) << " "
          //            << v(center_y-ip) << " "
          //            << v_rhs(center_x) << " "
          //            << delta_Rx << " "
          //            << p(center) << " "
          //            << p(left) << " "
          //            << edge_viscosity(up_e) << " "
          //            << edge_viscosity(center_e) << " "
          //            << cell_viscosity(center) << " "
          //            << cell_viscosity(left) << " "
          //            << dx << " "
          //            << dy << " "
          //            << dtau_xx_dx << " "
          //            << dtau_xy_dy << " "
          //            << dp_dx << " "
          //            << "\n";
          }

          v(center_x)+=delta_Rx*theta_momentum/C_vx;

          /* Set the boundary elements so that the
             derivative is zero. */
          if(set_boundary)
            {
              v(center_x+offset)=v(center_x) + dv;
                
              // tbox::plog << "setbc "
              //            << (center+offset)(0) << " "
              //            << (center+offset)(1) << " "
              //            // << center(0) << " "
              //            // << center(1) << " "
              //            // << offset(0) << " "
              //            // << offset(1) << " "
              //            // << pbox.lower(0) << " "
              //            // << pbox.upper(0) << " "
              //            // << pbox.lower(1) << " "
              //            // << pbox.upper(1) << " "
              //            << v(pdat::SideIndex(center+offset,axis,
              //                                    pdat::SideIndex::Lower)) << " "
              //            << v(pdat::SideIndex(center,axis,pdat::SideIndex::Lower)) << " "
              //            << dv << " ";
            }
        }
    }
}
