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
#include "set_V_boundary.h"

/* Set the physical boundaries for the velocity. */

void SAMRAI::solv::StokesFACOps::set_boundaries
(const int &v_id, tbox::Pointer<hier::PatchLevel> &level)
{
  for (hier::PatchLevel::Iterator pi(*level); pi; pi++)
    {
      tbox::Pointer<hier::Patch> patch = *pi;
      set_V_boundary(*patch,v_id);

      // tbox::Pointer<pdat::SideData<double> >
      //   v = patch->getPatchData(v_id);

      // hier::Box pbox=patch->getBox();
      // tbox::Pointer<geom::CartesianPatchGeometry>
      //   geom = patch->getPatchGeometry();

      // hier::Box gbox=v->getGhostBox();

      // tbox::plog << "boundary "
      //            << gbox.lower(0) << " "
      //            << gbox.upper(0) << " "
      //            << gbox.lower(1) << " "
      //            << gbox.upper(1) << " "
      //            << pbox.lower(0) << " "
      //            << pbox.upper(0) << " "
      //            << pbox.lower(1) << " "
      //            << pbox.upper(1) << " "
      //            << "\n";
      // for(int j=gbox.lower(1); j<=gbox.upper(1)+1; ++j)
      //   for(int i=gbox.lower(0); i<=gbox.upper(0)+1; ++i)
      //     {
      //       pdat::CellIndex center(tbox::Dimension(2));
      //       center[0]=i;
      //       center[1]=j;
      //       hier::Index ip(1,0), jp(0,1);

      //       /* vx */
      //       if(j<=gbox.upper(1))
      //         {
      //           /* Set a sentinel value */
      //           if((i<pbox.lower(0) && geom->getTouchesRegularBoundary(0,0))
      //              || (i>pbox.upper(0)+1
      //                  && geom->getTouchesRegularBoundary(0,1)))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
      //                                    pdat::SideIndex::Lower))=
      //                 boundary_value;
      //             }
      //           /* Set the value so the derivative=0 */
      //           else if(j<pbox.lower(0)
      //                   && geom->getTouchesRegularBoundary(1,0))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
      //                                    pdat::SideIndex::Lower))=
      //                 (*v)(pdat::SideIndex(center+jp,pdat::SideIndex::X,
      //                                      pdat::SideIndex::Lower));
      //             }
      //           else if(j>pbox.upper(0)
      //                   && geom->getTouchesRegularBoundary(1,1))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
      //                                    pdat::SideIndex::Lower))=
      //                 (*v)(pdat::SideIndex(center-jp,pdat::SideIndex::X,
      //                                      pdat::SideIndex::Lower));
      //             }
      //           // tbox::plog << "set bc vx "
      //           //            << level->getLevelNumber() << " "
      //           //            << i << " "
      //           //            << j << " "
      //           //            << (*v)(pdat::SideIndex(center,pdat::SideIndex::X,
      //           //                                    pdat::SideIndex::Lower))
      //           //            << " "
      //           //            << "\n";
      //         }
      //       /* vy */
      //       if(i<=gbox.upper(0))
      //         {
      //           if((j<pbox.lower(1) && geom->getTouchesRegularBoundary(1,0))
      //              || (j>pbox.upper(1)+1
      //                  && geom->getTouchesRegularBoundary(1,1)))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
      //                                    pdat::SideIndex::Lower))=
      //                 boundary_value;
      //             }
      //           else if(i<pbox.lower(0)
      //                   && geom->getTouchesRegularBoundary(0,0))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
      //                                    pdat::SideIndex::Lower))=
      //                 (*v)(pdat::SideIndex(center+ip,pdat::SideIndex::Y,
      //                                      pdat::SideIndex::Lower));
      //             }
      //           else if(i>pbox.upper(0)
      //                   && geom->getTouchesRegularBoundary(0,1))
      //             {
      //               (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
      //                                    pdat::SideIndex::Lower))=
      //                 (*v)(pdat::SideIndex(center-ip,pdat::SideIndex::Y,
      //                                      pdat::SideIndex::Lower));
      //             }
      //           // tbox::plog << "set bc vy "
      //           //            << level->getLevelNumber() << " "
      //           //            << i << " "
      //           //            << j << " "
      //           //            << (*v)(pdat::SideIndex(center,pdat::SideIndex::Y,
      //           //                                    pdat::SideIndex::Lower))
      //           //            << " "
      //           //            << "\n";
      //         }
      //     }
    }

}
