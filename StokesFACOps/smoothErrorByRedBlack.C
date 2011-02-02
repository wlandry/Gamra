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
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/

void SAMRAI::solv::StokesFACOps::smoothErrorByRedBlack
(SAMRAIVectorReal<double>& solution,
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
      || residual.getPatchHierarchy() != d_hierarchy)
    {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                 "internal hierarchy.");
    }
#endif
  tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

  {
    hier::PatchLevel::Iterator pi(*level);
    for (pi.initialize(*level); pi; pi++) {
      tbox::Pointer<hier::Patch> patch = *pi;
      tbox::Pointer<pdat::CellData<double> > p_rhs_data =
        patch->getPatchData(p_rhs_id);
      for(pdat::CellIterator ci(p_rhs_data->getBox()); ci; ci++)
        {
          pdat::CellIndex cc=ci();
          tbox::plog << "p_rhs "
                     << cc[0] << " "
                     << cc[1] << " "
                     << (*p_rhs_data)(cc) << " "
                     // << (&(*p_rhs_data)(cc)) << " "
                     << "\n";
        }
    }
  }
  /* Only need to sync the rhs once. This sync is needed because
     calculating a new pressure update requires computing in the ghost
     region so that the update for the velocity inside the box will be
     correct. */
  v_refine_patch_strategy.setTargetDataId(v_id);
  xeqScheduleGhostFillNoCoarse(p_rhs_id,v_rhs_id,ln);
  set_boundaries(v_id,level);

  if (ln > d_ln_min) {
    /*
     * Perform a one-time transfer of data from coarser level,
     * to fill ghost boundaries that will not change through
     * the smoothing loop.
     */
    xeqScheduleGhostFill(p_id, v_id, ln);
  }

  double viscosity=1;
  double theta_momentum=1.2;
  double theta_continuity=0.3;

  /*
   * Smooth the number of sweeps specified or until
   * the convergence is satisfactory.
   */
  double maxres;
  /*
   * Instead of checking residual convergence globally, we check the
   * converged flag.  This avoids possible round-off errors affecting
   * different processes differently, leading to disagreement on
   * whether to continue smoothing.
   */
  bool converged = false;
  for (int sweep=0; sweep < num_sweeps*(1<<(d_ln_max-ln)) && !converged; ++sweep)
    {
      maxres=0;
      for(int rb=0;rb<2;++rb)
        {
          // Need to sync
          xeqScheduleGhostFillNoCoarse(p_id,v_id,ln);
          for (hier::PatchLevel::Iterator pi(*level); pi; pi++)
            {
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

                      tbox::plog << "smooth "
                                 << ln << " "
                                 << i << " "
                                 << j << " "
                                 << pbox.lower(0) << " "
                                 << pbox.upper(0) << " "
                                 << pbox.lower(1) << " "
                                 << pbox.upper(1) << " ";

                      /* Update p */
                      if(!((*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                pdat::SideIndex::Lower))
                           ==boundary_value
                           || (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
                                                   pdat::SideIndex::Upper))
                           ==boundary_value
                           || (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                   pdat::SideIndex::Lower))
                           ==boundary_value
                           || (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
                                                   pdat::SideIndex::Upper))
                           ==boundary_value))
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

                          /* No scaling here, though there should be. */
                          maxres=std::max(maxres,delta_R_continuity);

                          (*p)(center)+=
                            viscosity*delta_R_continuity*theta_continuity;

                          tbox::plog << "p "
                                     << (*p)(center) << " "
                                     << (*p_rhs)(center) << " "
                                     << dvx_dx << " "
                                     << dvy_dy << " "
                                     << delta_R_continuity << " ";
                        }
                      /* Update v */
                      Update_V(0,j,pbox,center,left,right,down,up,p,v,v_rhs,
                               maxres,dx,dy,viscosity,theta_momentum);
                      Update_V(1,i,pbox,center,down,up,left,right,p,v,v_rhs,
                               maxres,dy,dx,viscosity,theta_momentum);

                                // << (*v)(pdat::SideIndex
                                //         (center,
                                //          pdat::SideIndex::X,
                                //          pdat::SideIndex::Lower)) << " "
                                // << (*v)(pdat::SideIndex
                                //         (center,
                                //          pdat::SideIndex::Y,
                                //          pdat::SideIndex::Lower)) << " "
                      tbox::plog << "\n";
                    }
                }
            }
        }
      if (residual_tolerance >= 0.0) {
        /*
         * Check for early end of sweeps due to convergence
         * only if it is numerically possible (user gave a
         * non negative value for residual tolerance).
         */
        converged = maxres < residual_tolerance;
        const tbox::SAMRAI_MPI&
          mpi(d_hierarchy->getDomainMappedBoxLevel().getMPI());
        int tmp= converged ? 1 : 0;
        if (mpi.getSize() > 1)
          {
            mpi.AllReduce(&tmp, 1, MPI_MIN);
          }
        converged=(tmp==1);
        if (d_enable_logging)
          tbox::plog
            << d_object_name << "\n"
            << " RBGS sweep #" << sweep << " : " << maxres << "\n";
      }
    }
}

