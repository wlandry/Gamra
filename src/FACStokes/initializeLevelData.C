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
#include "SAMRAI/geom/CartesianGridGeometry.h"

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
 const double ,
 const bool ,
 const bool ,
 const tbox::Pointer<hier::BasePatchLevel> ,
 const bool allocate_data)
{
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

  // const double inclusion_radius=0.5;
  const double inclusion_viscosity=1e2;
  const double background_viscosity=1;

  // const double background_density(1), block_density(1);
  const double background_density(1), block_density(1.03);

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

    hier::Box cell_visc_box = cell_viscosity.getBox();
    for(pdat::CellIterator ci(cell_viscosity.getGhostBox()); ci; ci++)
      {
        pdat::CellIndex c=ci();
        double x=geom->getXLower()[0]
          + geom->getDx()[0]*(c[0]-cell_visc_box.lower()[0] + 0.5);
        double y=geom->getXLower()[1]
          + geom->getDx()[1]*(c[1]-cell_visc_box.lower()[1] + 0.5);

        if(!viscosity.empty())
          {
            int i(static_cast<int>(x*(viscosity_ijk[0]-1)
                                   /(viscosity_xyz_max[0]-viscosity_xyz_min[0]))),
              j(static_cast<int>(y*(viscosity_ijk[1]-1)
                                 /(viscosity_xyz_max[1]-viscosity_xyz_min[1])));
            i=std::max(0,std::min(viscosity_ijk[0]-1,i));
            j=std::max(0,std::min(viscosity_ijk[1]-1,j));

            if(d_dim.getValue()==2)
              cell_viscosity(c)=viscosity[i+viscosity_ijk[0]*j];
            else
              {
                double z=geom->getXLower()[2]
                  + geom->getDx()[2]*(c[2]-cell_visc_box.lower()[2] + 0.5);
                int k(static_cast<int>(z*(viscosity_ijk[2]-1)
                                       /(viscosity_xyz_max[2]
                                         - viscosity_xyz_min[2])));
                k=std::max(0,std::min(viscosity_ijk[2]-1,k));
                cell_viscosity(c)=
                  viscosity[i+viscosity_ijk[0]*(j+viscosity_ijk[1]*k)];
              }
          }
        else
          {
            // if(x*x + y*y < inclusion_radius*inclusion_radius)
            //   cell_viscosity(c)=inclusion_viscosity;
            // else
            //   cell_viscosity(c)=background_viscosity;

            if(d_dim.getValue()==2)
              {
                if(x<1.0/3 || x>2.0/3 || y<1.0/3 || y>2.0/3)
                  cell_viscosity(c)=background_viscosity;
                else
                  cell_viscosity(c)=inclusion_viscosity;
              }
            else
              {
                double z=geom->getXLower()[2]
                  + geom->getDx()[2]*(c[2]-cell_visc_box.lower()[2] + 0.5);
                if(x<1.0/3 || x>2.0/3 || y<1.0/3 || y>2.0/3 || z<1.0/3 || z>2.0/3)
                  cell_viscosity(c)=background_viscosity;
                else
                  cell_viscosity(c)=inclusion_viscosity;
              }
          }
      }

    tbox::Pointer<pdat::CellData<double> > dp_data =
      patch->getPatchData(dp_id);

    /* I do not think this is actually necessary. */
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

        double x=geom->getXLower()[0]
          + geom->getDx()[0]*(s[0]-pbox.lower()[0]+0.5);
        double y=geom->getXLower()[1]
          + geom->getDx()[1]*(s[1]-pbox.lower()[1]);

        if(d_dim.getValue()==2)
          {
            if(x<1.0/3 || x>2.0/3 || y<1.0/3 || y>2.0/3)
              (*v_rhs_data)(s)=background_density;
            else
              (*v_rhs_data)(s)=block_density;
          }
        else
          {
            double z=geom->getXLower()[2]
              + geom->getDx()[2]*(s[2]-pbox.lower()[2]);
            if(x<1.0/3 || x>2.0/3 || y<1.0/3 || y>2.0/3 || z<1.0/3 || z>2.0/3)
              (*v_rhs_data)(s)=background_density;
            else
              (*v_rhs_data)(s)=block_density;
          }
            
        // (*v_rhs_data)(s)=10;
        // (*v_rhs_data)(s)=0;
      }
    if(d_dim.getValue()==3)
      for(pdat::SideIterator si(pbox,2); si; si++)
        {
          pdat::SideIndex s=si();
          (*v_rhs_data)(s)=0;
        }
  }    // End patch loop.
}
