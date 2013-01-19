/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "FTensor.hpp"

/* Initialize data on a level, including allocating memory.  It does
   not set the terms due to the faults, since the moduli have to be
   fixed first. */

void Elastic::FAC::initializeLevelData
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& patch_hierarchy,
 const int level_number,
 const double ,
 const bool ,
 const bool ,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel>& ,
 const bool allocate_data)
{
  boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
    hierarchy = patch_hierarchy;
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());

  boost::shared_ptr<SAMRAI::hier::PatchLevel> level =
    hierarchy->getPatchLevel(level_number);
  const int dim=d_dim.getValue();

  if (allocate_data) {
    level->allocatePatchData(cell_moduli_id);
    level->allocatePatchData(edge_moduli_id);
    level->allocatePatchData(v_id);
    level->allocatePatchData(v_rhs_id);
    level->allocatePatchData(dv_diagonal_id);
    level->allocatePatchData(dv_mixed_id);
  }

  /*
   * Initialize data in all patches in the level.
   */
  SAMRAI::hier::PatchLevel::Iterator p_i(level->begin());
  for (; p_i!=level->end(); p_i++) {

    boost::shared_ptr<SAMRAI::hier::Patch> patch = *p_i;
    if (!patch) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find patch.  Null patch pointer.");
    }
    boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
      boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>
      (patch->getPatchGeometry());
    const double *dx=geom->getDx();

    /* Initialize cell moduli */
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
      (patch->getPatchData(cell_moduli_id));

    SAMRAI::hier::Box cell_moduli_box = cell_moduli->getBox();

    SAMRAI::pdat::CellIterator cend(cell_moduli->getGhostBox(),false);
    for(SAMRAI::pdat::CellIterator ci(cell_moduli->getGhostBox(),true);
        ci!=cend; ci++)
      {
        SAMRAI::pdat::CellIndex c=*ci;
        double xyz[3];
        for(int d=0;d<dim;++d)
          xyz[d]=geom->getXLower()[d]
            + dx[d]*(c[d]-cell_moduli_box.lower()[d] + 0.5);

        (*cell_moduli)(c,0)=lambda.eval(xyz);
        (*cell_moduli)(c,1)=mu.eval(xyz);
      }

    /* v_rhs */
    boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_data =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
      (patch->getPatchData(v_rhs_id));

    v_rhs_data->fillAll(0);
    /* FIXME: need to add in the v_rhs from the input file */
  }
}
