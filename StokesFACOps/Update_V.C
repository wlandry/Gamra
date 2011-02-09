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
 const double &viscosity,
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
          if(center[axis]==pbox.lower(axis)+1)
            {
              offset[axis]=-1;
              set_boundary=true;
            }
          else if(center[axis]==pbox.upper(axis)-1)
            {
              offset[axis]=1;
              set_boundary=true;
            }

          double dv;
          if(set_boundary)
            {
              dv=(*v)(pdat::SideIndex
                      (center-offset,
                       axis,
                       pdat::SideIndex::Lower))
                - (*v)(pdat::SideIndex
                       (center+offset,axis,
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

          C_vx=-2*viscosity*(1/(dx*dx) + 1/(dy*dy));

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
            - viscosity*(d2vx_dxx + d2vx_dyy) + dp_dx;


          tbox::plog << "v " << axis << " "
                     << (*v)(pdat::SideIndex(center,
                                             axis,
                                             pdat::SideIndex::Lower))
                     << " "
                     << (*v_rhs)(pdat::SideIndex(center,
                                                 axis,
                                                 pdat::SideIndex::Lower))
                     << " ";

          // tbox::plog << "Update "
          //            << axis << " "
          //            << off_axis << " "
          //            << j << " "
          //            << center[axis] << " "
          //            << pbox.lower(axis) << " "
          //            << pbox.upper(axis) << " "
          //            << pbox.lower(off_axis) << " "
          //            << pbox.upper(off_axis) << " "
          //            << (*v_rhs)(pdat::SideIndex(center,
          //                                        axis,
          //                                        pdat::SideIndex::Lower)) << " "
          //            << right[0] << " "
          //            << right[1] << " "
          //            << (*v)(pdat::SideIndex
          //                    (right,
          //                     axis,
          //                     pdat::SideIndex::Lower)) << " "
          //            << delta_Rx << " "
          //            << (theta_momentum/C_vx) << " "
          //            << "\n";
            

          /* No scaling here, though there should be. */
          maxres=std::max(maxres,delta_Rx);

          (*v)(pdat::SideIndex(center,axis,
                               pdat::SideIndex::Lower))+=
            delta_Rx*theta_momentum/C_vx;

          /* Set the boundary elements so that the
             derivative is zero. */
          if(set_boundary)
            {
              (*v)(pdat::SideIndex
                   (center+offset,axis,
                    pdat::SideIndex::Lower))=
                (*v)(pdat::SideIndex(center-offset,axis,
                                     pdat::SideIndex::Lower))
                -dv;
            }
        }
    }
}
