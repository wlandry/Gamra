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

#include <unistd.h>

/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void SAMRAI::solv::StokesFACOps::relax(SAMRAIVectorReal<double>& solution,
                                       const SAMRAIVectorReal<double>& residual,
                                       int ln,
                                       int num_sweeps,
                                       double residual_tolerance)
{
  const int p_id(solution.getComponentDescriptorIndex(0)),
    p_rhs_id(residual.getComponentDescriptorIndex(0)),
    v_id(solution.getComponentDescriptorIndex(1)),
    v_rhs_id(residual.getComponentDescriptorIndex(1));

  checkInputPatchDataIndices();

#ifdef DEBUG_CHECK_ASSERTIONS
  if (solution.getPatchHierarchy() != d_hierarchy
      || residual.getPatchHierarchy() != d_hierarchy) {
    TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
               "internal hierarchy.");
  }
#endif
  tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  // if (ln > d_ln_min) {
  //   /*
  //    * Perform a one-time transfer of data from coarser level,
  //    * to fill ghost boundaries that will not change through
  //    * the smoothing loop.
  //    */
  //   xeqScheduleGhostFill(data_id, ln);
  // }

  double viscosity=1;
  double theta_momentum=1.2;
  double theta_continuity=0.3;

  /*
   * Smooth the number of sweeps specified or until
   * the convergence is satisfactory.
   */
  int isweep(0);
  double red_maxres, blk_maxres, maxres = 0;
  red_maxres = blk_maxres = residual_tolerance + 1;
  /*
   * Instead of checking residual convergence globally,
   * we check the not_converged flag.  This avoids possible
   * round-off errors affecting different processes differently,
   * leading to disagreement on whether to continue smoothing.
   */
  bool converged = false;
  for (; isweep < num_sweeps && !converged; ++isweep) {
    red_maxres = blk_maxres = 0;

    for(int rb=0;rb<2;++rb)
      {
        // Need to sync
        tbox::plog << "syncing\n";
        xeqScheduleGhostFillNoCoarse(p_id,v_id,ln);
        for (hier::PatchLevel::Iterator pi(*level); pi; pi++) {
          tbox::Pointer<hier::Patch> patch = *pi;

          tbox::Pointer<pdat::CellData<double> >
            p = patch->getPatchData(p_id);
          tbox::Pointer<pdat::CellData<double> >
            p_rhs = patch->getPatchData(p_rhs_id);
          tbox::Pointer<pdat::SideData<double> >
            v = patch->getPatchData(v_id);
          tbox::Pointer<pdat::SideData<double> >
            v_rhs = patch->getPatchData(v_rhs_id);

          hier::Box pbox=patch->getBox();
          tbox::Pointer<geom::CartesianPatchGeometry>
            geom = patch->getPatchGeometry();
          double dx = *(geom->getDx());
          double dy = *(geom->getDx());

          for(int j=pbox.lower(1); j<=pbox.upper(1)+1; ++j)
            {
              /* Do the red-black skip */
              int i_min=pbox.lower(0) + (abs(pbox.lower(0) + j + rb))%2;
              for(int i=i_min; i<=pbox.upper(0)+1; i+=2)
                {
                  pdat::CellIndex center(tbox::Dimension(2));
                  center[0]=i;
                  center[1]=j;

                  pdat::CellIndex up(center), down(center), right(center),
                    left(center);

                  ++up[1];
                  --down[1];
                  ++right[0];
                  --left[0];

                  /* Update p */
                  if(i!=pbox.upper(0)+1 && j!=pbox.upper(1)+1)
                    {
                      double dvx_dx=
                        ((*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                              pdat::SideIndex::Upper))
                         - (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                pdat::SideIndex::Lower)))/dx;
                      double dvy_dy=
                        ((*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                              pdat::SideIndex::Upper))
                         - (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                pdat::SideIndex::Lower)))/dy;

                      double delta_R_continuity=
                        (*p_rhs)(center) - dvx_dx - dvy_dy;

                      (*p)(center)+=
                        viscosity*delta_R_continuity*theta_continuity;
                    }

                  /* Update vx */
                  if(j!=pbox.upper(1)+1)
                    {
                      /* If x==0 */
                      if((center[0]==pbox.lower(0)
                          && geom->getTouchesRegularBoundary(0,0))
                         || (center[0]==pbox.upper(0)+1
                             && geom->getTouchesRegularBoundary(0,1)))
                        {
                          (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                               pdat::SideIndex::Lower))=0;

                          // tbox::plog << "vx x0 boundary ";
                        }
                      else
                        {
                          double dp_dx, d2vx_dxx, d2vx_dyy, C_vx;
                          /* If y==0 */
                          if(center[1]==pbox.lower(1)
                             && geom->getTouchesRegularBoundary(1,0))
                            {
                              d2vx_dyy=
                                ((*v)(pdat::SideIndex(up,pdat::SideIndex::X,
                                                      pdat::SideIndex::Lower))
                                 - (*v)(pdat::SideIndex(center,
                                                        pdat::SideIndex::X,
                                                        pdat::SideIndex::Lower)))
                                /(dy*dy);
                              C_vx=-viscosity*(2/(dx*dx) + 1/(dy*dy));
                            }
                          /* If y==max_y */
                          else if(center[1]==pbox.upper(1)
                                  && geom->getTouchesRegularBoundary(1,1))
                            {
                              d2vx_dyy=
                                (-(*v)(pdat::SideIndex(center,
                                                       pdat::SideIndex::X,
                                                       pdat::SideIndex::Lower))
                                 + (*v)(pdat::SideIndex(down,
                                                        pdat::SideIndex::X,
                                                        pdat::SideIndex::Lower)))
                                /(dy*dy);
                              C_vx=-viscosity*(2/(dx*dx) + 1/(dy*dy));
                              // tbox::plog << "vx y1 boundary ";
                            }
                          else
                            {
                              d2vx_dyy=
                                ((*v)(pdat::SideIndex(up,pdat::SideIndex::X,
                                                      pdat::SideIndex::Lower))
                                 - 2*(*v)(pdat::SideIndex(center,
                                                          pdat::SideIndex::X,
                                                          pdat::SideIndex::Lower))
                                 + (*v)(pdat::SideIndex(down,
                                                        pdat::SideIndex::X,
                                                        pdat::SideIndex::Lower)))
                                /(dy*dy);

                              C_vx=-2*viscosity*(1/(dx*dx) + 1/(dy*dy));
                            }
                          d2vx_dxx=((*v)(pdat::SideIndex(left,
                                                         pdat::SideIndex::X,
                                                         pdat::SideIndex::Lower))
                                    - 2*(*v)(pdat::SideIndex(center,
                                                             pdat::SideIndex::X,
                                                             pdat::SideIndex::Lower))
                                    + (*v)(pdat::SideIndex(right,
                                                           pdat::SideIndex::X,
                                                           pdat::SideIndex::Lower)))
                            /(dx*dx);

                          dp_dx=((*p)(center)-(*p)(left))/dx;
                              
                          double delta_Rx=
                            (*v_rhs)(pdat::SideIndex(center,
                                                     pdat::SideIndex::X,
                                                     pdat::SideIndex::Lower))
                            - viscosity*(d2vx_dxx + d2vx_dyy) + dp_dx;
                          (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                               pdat::SideIndex::Lower))+=
                            delta_Rx*theta_momentum/C_vx;
                        }
                    }

                  /* Update vy */
                  if(i!=pbox.upper(0)+1)
                    {
                      /* If y==0 */
                      if((center[1]==pbox.lower(1)
                          && geom->getTouchesRegularBoundary(1,0))
                         || (center[1]==pbox.upper(1)+1
                             && geom->getTouchesRegularBoundary(1,1)))
                        {
                          (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                               pdat::SideIndex::Lower))=0;
                        }
                      else
                        {
                          double dp_dy, d2vy_dxx, d2vy_dyy, C_vy;
                          /* If x==0 */
                          if(center[0]==pbox.lower(0)
                             && geom->getTouchesRegularBoundary(0,0))
                            {
                              d2vy_dxx=
                                ((*v)(pdat::SideIndex(right,pdat::SideIndex::Y,
                                                      pdat::SideIndex::Lower))
                                 - (*v)(pdat::SideIndex(center,
                                                        pdat::SideIndex::Y,
                                                        pdat::SideIndex::Lower)))
                                /(dx*dx);
                              C_vy=-viscosity*(1/(dx*dx) + 2/(dy*dy));
                            }
                          /* If x==max_x */
                          else if(center[0]==pbox.upper(0)
                                  && geom->getTouchesRegularBoundary(0,1))
                            {
                              d2vy_dxx=
                                ((*v)(pdat::SideIndex(left,pdat::SideIndex::Y,
                                                      pdat::SideIndex::Lower))
                                 - (*v)(pdat::SideIndex(center,
                                                        pdat::SideIndex::Y,
                                                        pdat::SideIndex::Lower)))
                                /(dx*dx);
                              C_vy=-viscosity*(1/(dx*dx) + 2/(dy*dy));
                            }
                          else
                            {
                              d2vy_dxx=
                                ((*v)(pdat::SideIndex(left,pdat::SideIndex::Y,
                                                      pdat::SideIndex::Lower))
                                 - 2*(*v)(pdat::SideIndex(center,
                                                          pdat::SideIndex::Y,
                                                          pdat::SideIndex::Lower))
                                 + (*v)(pdat::SideIndex(right,
                                                        pdat::SideIndex::Y,
                                                        pdat::SideIndex::Lower)))
                                /(dx*dx);

                              C_vy=-2*viscosity*(1/(dx*dx) + 1/(dy*dy));
                            }
                          d2vy_dyy=((*v)(pdat::SideIndex(up,pdat::SideIndex::Y,
                                                         pdat::SideIndex::Lower))
                                    - 2*(*v)(pdat::SideIndex(center,
                                                             pdat::SideIndex::Y,
                                                             pdat::SideIndex::Lower))
                                    + (*v)(pdat::SideIndex(down,
                                                           pdat::SideIndex::Y,
                                                           pdat::SideIndex::Lower)))
                            /(dy*dy);

                          dp_dy=((*p)(center)-(*p)(down))/dy;
                              
                          double delta_Ry=
                            (*v_rhs)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                     pdat::SideIndex::Lower))
                            - viscosity*(d2vy_dxx + d2vy_dyy) + dp_dy;
                          (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                               pdat::SideIndex::Lower))+=
                            delta_Ry*theta_momentum/C_vy;
                        }
                    }
                }
            }
        }
      }
  }

  // if (residual_tolerance >= 0.0) {
  //   /*
  //    * Check for early end of sweeps due to convergence
  //    * only if it is numerically possible (user gave a
  //    * non negative value for residual tolerance).
  //    */
  //   maxres = tbox::MathUtilities<double>::Max(red_maxres, blk_maxres);
  //   not_converged = maxres > residual_tolerance;
  //   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getDomainMappedBoxLevel().getMPI());
  //   if (mpi.getSize() > 1) {
  //     mpi.AllReduce(&not_converged, 1, MPI_MAX);
  //   }
  // }

  // if (d_enable_logging) tbox::plog
  //                         << d_object_name << " RBGS smoothing maxres = " << maxres << "\n"
  //                         << "  after " << isweep << " sweeps.\n";

}

