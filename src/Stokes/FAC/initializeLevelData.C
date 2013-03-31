/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#include "Stokes/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

/*
*************************************************************************
* Initialize data on a level.                                           *
*                                                                       *
* Allocate the solution, exact solution and rhs memory.                 *
* Fill the rhs and exact solution.                                      *
*************************************************************************
*/
void Stokes::FAC::initializeLevelData
(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& patch_hierarchy,
 const int level_number,
 const double ,
 const bool ,
 const bool ,
 const boost::shared_ptr<SAMRAI::hier::PatchLevel>& ,
 const bool allocate_data)
{
  boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy = patch_hierarchy;
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geom =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>(hierarchy->getGridGeometry());

  boost::shared_ptr<SAMRAI::hier::PatchLevel> level =
    hierarchy->getPatchLevel(level_number);
  const int dim=d_dim.getValue();

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

  /*
   * Initialize data in all patches in the level.
   */
  SAMRAI::hier::PatchLevel::Iterator pi(level->begin());
  for (; pi!=level->end(); ++pi) {

    boost::shared_ptr<SAMRAI::hier::Patch> patch = *pi;
    if (!patch) {
      TBOX_ERROR(d_object_name
                 << ": Cannot find patch.  Null patch pointer.");
    }
    boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> geom =
      boost::dynamic_pointer_cast<SAMRAI::geom::CartesianPatchGeometry>(patch->getPatchGeometry());
    const double *dx=geom->getDx();

    /* Initialize cell viscosity */
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_viscosity =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >(patch->getPatchData(cell_viscosity_id));

    SAMRAI::hier::Box cell_visc_box = cell_viscosity->getBox();
    SAMRAI::pdat::CellIterator cend(cell_viscosity->getGhostBox(),false);
    for(SAMRAI::pdat::CellIterator ci(cell_viscosity->getGhostBox(),true); ci!=cend; ++ci)
      {
        SAMRAI::pdat::CellIndex c=*ci;
        double xyz[dim];
        for(int d=0;d<dim;++d)
          xyz[d]=geom->getXLower()[d]
            + dx[d]*(c[d]-cell_visc_box.lower()[d] + 0.5);

        int ijk(0), factor(1);
        for(int d=0;d<dim;++d)
          {
            int i=static_cast<int>(xyz[d]*(viscosity_ijk[d]-1)
                                   /(viscosity_xyz_max[d]-viscosity_xyz_min[d]));
            i=std::max(0,std::min(viscosity_ijk[d]-1,i));
            ijk+=i*factor;
            factor*=viscosity_ijk[d];
          }
        (*cell_viscosity)(c)=viscosity[ijk];
      }

    /* I do not think this is actually necessary. */
    boost::shared_ptr<SAMRAI::pdat::CellData<double> > dp_data =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >(patch->getPatchData(dp_id));
    dp_data->fill(0.0);

    boost::shared_ptr<SAMRAI::pdat::CellData<double> > p_rhs_data =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >(patch->getPatchData(p_rhs_id));
    p_rhs_data->fill(0.0);

    /* v_rhs */
    boost::shared_ptr<SAMRAI::pdat::SideData<double> > v_rhs_data =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >(patch->getPatchData(v_rhs_id));

    if(v_rhs.empty())
      {
        v_rhs_data->fill(0,0);
      }
    else
      {
        SAMRAI::hier::Box pbox = v_rhs_data->getBox();
        int ix_offset(0);
        for(int ix=0;ix<dim;++ix)
          {
            double offset[]={0.5,0.5,0.5};
            offset[ix]=0;

            SAMRAI::pdat::SideIterator send(pbox,ix,false);
            for(SAMRAI::pdat::SideIterator si(pbox,ix,true); si!=send; ++si)
              {
                SAMRAI::pdat::SideIndex s=*si;
                double xyz[dim];
                for(int d=0;d<dim;++d)
                  xyz[d]=geom->getXLower()[d]
                    + dx[d]*(s[d]-pbox.lower()[d]+offset[d]);
            
                int ijk(0), factor(1);
                for(int d=0;d<dim;++d)
                  {
                    int i=static_cast<int>(xyz[d]*(v_rhs_ijk[d]-1)
                                           /(v_rhs_xyz_max[d]-v_rhs_xyz_min[d]));
                    i=std::max(0,std::min(v_rhs_ijk[d]-1,i));
                    ijk+=i*factor;
                    factor*=v_rhs_ijk[d];
                  }
                (*v_rhs_data)(s)=v_rhs[ijk+ix_offset];
              }
            int i=1;
            for(int d=0;d<dim;++d)
              i*=v_rhs_ijk[d];
            ix_offset+=i;
          }
      }
  }    // End patch loop.
}
