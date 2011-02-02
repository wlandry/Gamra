/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "FACStokes.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "StokesSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

extern "C" {
  void F77_FUNC(setexactandrhs2d, SETEXACTANDRHS2D) (const int & ifirst0,
                                                     const int & ilast0,
                                                     const int & ifirst1,
                                                     const int & ilast1,
                                                     double* exact,
                                                     double* rhs,
                                                     const double* dx,
                                                     const double* xlower);
}

namespace SAMRAI {

  /*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
  void FACStokes::initializeLevelData
  (const tbox::Pointer<hier::BasePatchHierarchy> patch_hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer<hier::BasePatchLevel> old_level,
   const bool allocate_data)
  {

    (void)init_data_time;
    (void)can_be_refined;
    (void)initial_time;
    (void)old_level;

    tbox::Pointer<hier::PatchHierarchy> hierarchy = patch_hierarchy;
    tbox::Pointer<geom::CartesianGridGeometry> grid_geom =
      hierarchy->getGridGeometry();

    tbox::Pointer<hier::PatchLevel> level =
      hierarchy->getPatchLevel(level_number);

    if (allocate_data) {
      level->allocatePatchData(p_id);
      level->allocatePatchData(p_rhs_id);
      level->allocatePatchData(p_exact_id);
      level->allocatePatchData(v_id);
      level->allocatePatchData(v_rhs_id);
    }

    /*
     * Initialize data in all patches in the level.
     */
    hier::PatchLevel::Iterator pi(*level);
    for (pi.initialize(*level); pi; pi++) {

      tbox::Pointer<hier::Patch> patch = *pi;
      if (patch.isNull()) {
        TBOX_ERROR(d_object_name
                   << ": Cannot find patch.  Null patch pointer.");
      }
      tbox::Pointer<pdat::CellData<double> > p_rhs_data =
        patch->getPatchData(p_rhs_id);

      p_rhs_data->fill(0.0);

      // for(pdat::CellIterator ci(p_rhs_data->getBox()); ci; ci++)
      //   {
      //     pdat::CellIndex cc=ci();
      //     tbox::plog << "p_rhs "
      //                << cc[0] << " "
      //                << cc[1] << " "
      //                << (*p_rhs_data)(cc) << " "
      //                << (&(*p_rhs_data)(cc)) << " "
      //                << "\n";
      //   }

      tbox::Pointer<pdat::SideData<double> > v_rhs_data =
        patch->getPatchData(v_rhs_id);

      hier::Box pbox = v_rhs_data->getBox();

      for(pdat::SideIterator si(pbox,0); si; si++)
        {
          pdat::SideIndex s=si();
          (*v_rhs_data)(s)=0;
        }
      for(pdat::SideIterator si(pbox,1); si; si++)
        {
          pdat::SideIndex s=si();

          (*v_rhs_data)(s)=10;
        }
    }    // End patch loop.
  }

}
