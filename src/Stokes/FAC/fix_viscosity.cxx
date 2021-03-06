#include "Stokes/FAC.hxx"
#include <SAMRAI/geom/CartesianGridGeometry.h>

/* Fix the viscosity on the coarse grids by coarsening from the finer
   grids. */

void Stokes::FAC::fix_viscosity()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  {
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
      cell_viscosity_coarsen_operator;
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
      cell_viscosity_coarsen_algorithm(new SAMRAI::xfer::CoarsenAlgorithm(d_dim));
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
      cell_viscosity_coarsen_schedules;

    SAMRAI::hier::VariableDatabase*
      vdb = SAMRAI::hier::VariableDatabase::getDatabase();
    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
      boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
      (d_hierarchy->getGridGeometry());
    boost::shared_ptr<SAMRAI::hier::Variable> variable;
    vdb->mapIndexToVariable(cell_viscosity_id, variable);
    cell_viscosity_coarsen_operator =
      geometry->lookupCoarsenOperator(variable,
                                      "CONSERVATIVE_COARSEN");

    if (!cell_viscosity_coarsen_operator)
      { TBOX_ERROR("Stokes::FAC: Cannot find cell viscosity coarsening "
                   "operator"); }

    cell_viscosity_coarsen_schedules.resize(ln_max + 1);
    cell_viscosity_coarsen_algorithm->
      registerCoarsen(cell_viscosity_id,cell_viscosity_id,
                      cell_viscosity_coarsen_operator);

    for (int dest_ln = 0; dest_ln < ln_max; ++dest_ln) {
      cell_viscosity_coarsen_schedules[dest_ln] =
        cell_viscosity_coarsen_algorithm->
        createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                       d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!cell_viscosity_coarsen_schedules[dest_ln])
        { TBOX_ERROR("Stokes::FAC: Cannot create a coarsen schedule for "
                     "cell viscosity restriction!\n"); }
    }

    for(int dest_ln=ln_max-1; dest_ln>=0; --dest_ln)
      {
        SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
        coarsener.registerCoarsen(cell_viscosity_id, cell_viscosity_id,
                                  cell_viscosity_coarsen_operator);
        coarsener.resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
        cell_viscosity_coarsen_schedules[dest_ln]->coarsenData();
        cell_viscosity_coarsen_algorithm->
          resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
      }
  }

  /* Compute edge_viscosity by averaging the cell viscosities. */

  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(d_dim)),
    jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(d_dim.getValue()>2)
    kp[2]=1;
  SAMRAI::hier::Index pp[]={ip,jp,kp};

  for (int level = 0; level <= d_hierarchy->getFinestLevelNumber(); ++level)
    {
      SAMRAI::hier::PatchLevel &patch_level = *d_hierarchy->getPatchLevel(level);
      for (SAMRAI::hier::PatchLevel::Iterator patch_iter(patch_level.begin());
           patch_iter!=patch_level.end(); ++patch_iter)
        {
          SAMRAI::hier::Patch &patch = **patch_iter;
          boost::shared_ptr<SAMRAI::pdat::CellData<double> >cell_viscosity_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(cell_viscosity_id));
          SAMRAI::pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);
          if(d_dim.getValue()==2)
            {
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> >
                edge_viscosity_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch.getPatchData(edge_viscosity_id));
              SAMRAI::pdat::NodeData<double>
                &edge_viscosity(*edge_viscosity_ptr);

              SAMRAI::pdat::NodeIterator
                nend(SAMRAI::pdat::NodeGeometry::end(edge_viscosity.getBox()));
              for(SAMRAI::pdat::NodeIterator
                    ni(SAMRAI::pdat::NodeGeometry::begin(edge_viscosity.getBox()));
                  ni!=nend; ++ni)
                {
                  const SAMRAI::pdat::NodeIndex &e(*ni);
                  SAMRAI::pdat::CellIndex c(e);
                  edge_viscosity(e)=
                    pow(cell_viscosity(c)*cell_viscosity(c-ip)
                        *cell_viscosity(c-jp)*cell_viscosity(c-ip-jp),0.25);
                }
            }
          else
            {
              boost::shared_ptr<SAMRAI::pdat::EdgeData<double> >
                edge_viscosity_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
                (patch.getPatchData(edge_viscosity_id));
              SAMRAI::pdat::EdgeData<double>
                &edge_viscosity(*edge_viscosity_ptr);
              for(Gamra::Dir axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  SAMRAI::hier::Box pbox=patch.getBox();
                  pbox.grow(axis,edge_viscosity.getGhostCellWidth()[axis]);
                  
                  SAMRAI::pdat::EdgeIterator
                    nend(SAMRAI::pdat::EdgeGeometry::end(pbox,axis));
                  for(SAMRAI::pdat::EdgeIterator
                        ni(SAMRAI::pdat::EdgeGeometry::begin(pbox,axis));
                      ni!=nend; ++ni)
                    {
                      const SAMRAI::pdat::EdgeIndex &e(*ni);
                      SAMRAI::pdat::CellIndex c(e);
                      edge_viscosity(e)=
                        pow(cell_viscosity(c)*cell_viscosity(c-pp[axis2])
                            *cell_viscosity(c-pp[axis3])
                            *cell_viscosity(c-pp[axis2]-pp[axis3]),0.25);
                    }
                }
            }
        }

      /* Ghost fill */
      SAMRAI::xfer::RefineAlgorithm refiner;
      refiner.registerRefine(edge_viscosity_id,edge_viscosity_id,
                             edge_viscosity_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());

      boost::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule=
        refiner.createSchedule(d_hierarchy->getPatchLevel(level));
        
      schedule->fillData(0.0,false);
    }
}
