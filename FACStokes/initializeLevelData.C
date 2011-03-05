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

/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void SAMRAI::FACStokes::initializeLevelData
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
    level->allocatePatchData(cell_viscosity_id);
    level->allocatePatchData(edge_viscosity_id);
    level->allocatePatchData(dp_id);
    level->allocatePatchData(p_rhs_id);
    level->allocatePatchData(p_exact_id);
    level->allocatePatchData(v_id);
    level->allocatePatchData(v_rhs_id);
  }

  const double inclusion_radius=0.5;
  const double inclusion_viscosity=10;
  const double background_viscosity=1;

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
    tbox::Pointer<geom::CartesianPatchGeometry>
      geom = patch->getPatchGeometry();

    tbox::Pointer<pdat::CellData<double> > cell_viscosity_ptr =
      patch->getPatchData(cell_viscosity_id);
    pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);

    hier::Box visc_box = cell_viscosity.getBox();
    for(pdat::CellIterator ci(cell_viscosity.getGhostBox()); ci; ci++)
      {
        pdat::CellIndex c=ci();
        double x=geom->getXLower()[0]
          + geom->getDx()[0]*(c[0]-visc_box.lower()[0] + 0.5);
        double y=geom->getXLower()[1]
          + geom->getDx()[1]*(c[1]-visc_box.lower()[1] + 0.5);

        if(x*x + y*y < inclusion_radius*inclusion_radius)
          cell_viscosity(c)=inclusion_viscosity;
        else
          cell_viscosity(c)=background_viscosity;

        // tbox::plog << "cell "
        //            << c[0] << " "
        //            << c[1] << " "
        //            << x << " "
        //            << y << " "
        //            << geom->getXLower()[0] << " "
        //            << geom->getXLower()[1] << " "
        //            << geom->getDx()[0] << " "
        //            << geom->getDx()[1] << " "
        //            << std::boolalpha
        //            << (x*x + y*y < inclusion_radius*inclusion_radius) << " "
        //            << cell_viscosity(c) << " "
        //            << "\n";

      }

    tbox::Pointer<pdat::NodeData<double> > edge_viscosity_ptr =
      patch->getPatchData(edge_viscosity_id);
    pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

    for(pdat::NodeIterator ei(edge_viscosity.getGhostBox()); ei; ei++)
      {
        pdat::NodeIndex e=ei();
        double x=geom->getXLower()[0]
          + geom->getDx()[0]*(e[0]-visc_box.lower()[0]);
        double y=geom->getXLower()[1]
          + geom->getDx()[1]*(e[1]-visc_box.lower()[1]);
        if(x*x + y*y < inclusion_radius*inclusion_radius)
          edge_viscosity(e)=inclusion_viscosity;
        else
          edge_viscosity(e)=background_viscosity;

        // tbox::plog << "edge "
        //            << e[0] << " "
        //            << e[1] << " "
        //            << x << " "
        //            << y << " "
        //            << geom->getXLower()[0] << " "
        //            << geom->getXLower()[1] << " "
        //            << geom->getDx()[0] << " "
        //            << geom->getDx()[1] << " "
        //            << std::boolalpha
        //            << (x*x + y*y < inclusion_radius*inclusion_radius) << " "
        //            << edge_viscosity(e) << " "
        //            << "\n";
      }


      // tbox::Pointer<pdat::CellData<double> > cell_viscosity_data =
      //   patch->getPatchData(cell_viscosity_id);

      // /* At some point this needs to do the proper interpolation for
      //    lower levels */
      // cell_viscosity_data->fill(1.0);

      // tbox::Pointer<pdat::NodeData<double> > edge_viscosity_data =
      //   patch->getPatchData(edge_viscosity_id);

      // /* At some point this needs to do the proper interpolation for
      //    lower levels */
      // edge_viscosity_data->fill(1.0);


    tbox::Pointer<pdat::CellData<double> > dp_data =
      patch->getPatchData(dp_id);

    /* This is mostly so that the outer boundaries are set properly. */
    dp_data->fill(0.0);

    tbox::Pointer<pdat::CellData<double> > p_rhs_data =
      patch->getPatchData(p_rhs_id);

    p_rhs_data->fill(0.0);

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
        // (*v_rhs_data)(s)=10;
        (*v_rhs_data)(s)=0;
      }
  }    // End patch loop.
}
