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
 tbox::Pointer<pdat::CellData<double> > &p,
 tbox::Pointer<pdat::SideData<double> > &v,
 tbox::Pointer<pdat::SideData<double> > &v_rhs,
 double &maxres,
 const double &dx,
 const double &dy,
 tbox::Pointer<pdat::CellData<double> > &cell_viscosity,
 tbox::Pointer<pdat::NodeData<double> > &node_viscosity,
 const double &theta_momentum)
{
  const int off_axis=(axis==0) ? 1 : 0;
  /* Update vx */
  if(j<pbox.upper(off_axis)+1)
    {
      /* If at the 'x' boundaries, leave vx as is */
      if(!((center[axis]==pbox.lower(axis)
            && (*v)(pdat::SideIndex(left,
                                    axis,
                                    pdat::SideIndex::Lower))
            ==boundary_value)
           || (center[axis]==pbox.upper(axis)+1
               && (*v)(pdat::SideIndex
                       (right,
                        axis,
                        pdat::SideIndex::Lower))
               ==boundary_value)))
        {
          double dp_dx, d2vx_dxx, d2vx_dyy, C_vx;
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
              dv=(*v)(pdat::SideIndex
                      (center+offset,
                       axis,
                       pdat::SideIndex::Lower))
                - (*v)(pdat::SideIndex
                       (center,axis,
                        pdat::SideIndex::Lower));
            }

          d2vx_dyy=
            ((*v)(pdat::SideIndex(up,axis,
                                  pdat::SideIndex::Lower))
             - 2*(*v)(pdat::SideIndex
                      (center,axis,
                       pdat::SideIndex::Lower))
             + (*v)(pdat::SideIndex
                    (down,axis,
                     pdat::SideIndex::Lower)))
            /(dy*dy);

          C_vx=-2*(*cell_viscosity)(center)*(1/(dx*dx) + 1/(dy*dy));
          // C_vx=-2*((*cell_viscosity)(center) + (*cell_viscosity)(left))/(dx*dx)
            // - ((*node_viscosity)(up) + (*node_viscosity)(center))/(dx*dx)/(dy*dy));

          d2vx_dxx=((*v)(pdat::SideIndex
                         (left,axis,
                          pdat::SideIndex::Lower))
                    - 2*(*v)(pdat::SideIndex
                             (center,axis,
                              pdat::SideIndex::Lower))
                    + (*v)(pdat::SideIndex
                           (right,axis,
                            pdat::SideIndex::Lower)))
            /(dx*dx);

          dp_dx=((*p)(center)-(*p)(left))/dx;
                              
          double delta_Rx=
            (*v_rhs)(pdat::SideIndex(center,
                                     axis,
                                     pdat::SideIndex::Lower))
            - (*cell_viscosity)(center)*(d2vx_dxx + d2vx_dyy) + dp_dx;


          /* No scaling here, though there should be. */
          maxres=std::max(maxres,std::fabs(delta_Rx));

          // tbox::plog << "v " << axis << " "
          //            // << (*v)(pdat::SideIndex(center,
          //            //                         axis,
          //            //                         pdat::SideIndex::Lower))
          //            // << " "
          //            // << maxres << " "
          //            // << (*v_rhs)(pdat::SideIndex(center,
          //            //                             axis,
          //            //                             pdat::SideIndex::Lower))
          //            // << " "
          //            // << &(*v_rhs)(pdat::SideIndex(center,
          //            //                             axis,
          //            //                             pdat::SideIndex::Lower))
          //            // << " ";
          //            << delta_Rx << " ";
          //            // << d2vx_dxx << " "
          //            // << d2vx_dyy << " "
          //            // << dp_dx << " "
          //            // << (*v)(pdat::SideIndex(left,axis,
          //            //                         pdat::SideIndex::Lower)) << " "
          //            // << (*v)(pdat::SideIndex
          //            //         (center,axis,
          //            //          pdat::SideIndex::Lower)) << " "
          //            // << (*v)(pdat::SideIndex
          //            //         (right,axis,
          //            //          pdat::SideIndex::Lower)) << " "
          //            // << (*v)(pdat::SideIndex(up,axis,
          //            //                         pdat::SideIndex::Lower)) << " "
          //            // << (*v)(pdat::SideIndex
          //            //         (down,axis,
          //            //          pdat::SideIndex::Lower)) << " ";

          //            // << (*p)(center) << " "
          //            // << (*p)(left) << " ";

          //            // << &(*v_rhs)(pdat::SideIndex(center,
          //            //                             axis,
          //            //                             pdat::SideIndex::Lower))
          //            // << " "
          //            // << std::boolalpha
          //            // << set_boundary << " ";

          (*v)(pdat::SideIndex(center,axis,
                               pdat::SideIndex::Lower))+=
            delta_Rx*theta_momentum/C_vx;

          /* Set the boundary elements so that the
             derivative is zero. */
          if(set_boundary)
            {
              (*v)(pdat::SideIndex(center+offset,axis,
                                   pdat::SideIndex::Lower))=
                (*v)(pdat::SideIndex(center,axis,pdat::SideIndex::Lower)) + dv;
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
              //            << (*v)(pdat::SideIndex(center+offset,axis,
              //                                    pdat::SideIndex::Lower)) << " "
              //            << (*v)(pdat::SideIndex(center,axis,pdat::SideIndex::Lower)) << " "
              //            << dv << " ";
            }
        }
    }
}
