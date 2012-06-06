#include "Stokes/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

/* Fix the viscosity on the coarse grids by coarsening from the finer
   grids. */

void SAMRAI::Stokes::FAC::fix_viscosity()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  tbox::Pointer<xfer::CoarsenOperator> cell_viscosity_coarsen_operator;
  tbox::Pointer<xfer::CoarsenAlgorithm> cell_viscosity_coarsen_algorithm;
  tbox::Array<tbox::Pointer<xfer::CoarsenSchedule> >
    cell_viscosity_coarsen_schedules;

  hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
  tbox::Pointer<geom::CartesianGridGeometry> geometry =
    d_hierarchy->getGridGeometry();
  tbox::Pointer<hier::Variable> variable;
  vdb->mapIndexToVariable(cell_viscosity_id, variable);
  cell_viscosity_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CONSERVATIVE_COARSEN");
                                    // "CELL_VISCOSITY_COARSEN");

  if (!cell_viscosity_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find cell viscosity coarsening operator");
  }

  cell_viscosity_coarsen_schedules.resizeArray(ln_max + 1);
  cell_viscosity_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  cell_viscosity_coarsen_algorithm->
    registerCoarsen(cell_viscosity_id,cell_viscosity_id,
                    cell_viscosity_coarsen_operator);

  for (int dest_ln = 0; dest_ln < ln_max; ++dest_ln) {
    cell_viscosity_coarsen_schedules[dest_ln] =
      cell_viscosity_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!cell_viscosity_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for cell viscosity restriction!\n");
    }
  }

  for(int dest_ln=ln_max-1; dest_ln>=0; --dest_ln)
    {
      xfer::CoarsenAlgorithm coarsener(d_dim);
      coarsener.registerCoarsen(cell_viscosity_id, cell_viscosity_id,
                                cell_viscosity_coarsen_operator);
      coarsener.resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
      cell_viscosity_coarsen_schedules[dest_ln]->coarsenData();
      cell_viscosity_coarsen_algorithm->
        resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
    }

  cell_viscosity_coarsen_algorithm.setNull();
  cell_viscosity_coarsen_schedules.setNull();

  /* Compute edge_viscosity by averaging the cell viscosities. */

  hier::Index ip(hier::Index::getZeroIndex(d_dim)), jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(d_dim.getValue()>2)
    kp[2]=1;
  hier::Index pp[]={ip,jp,kp};

  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
      tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator i_p(*level);
      for ( ; i_p; i_p++)
        {
          tbox::Pointer<hier::Patch> patch = *i_p;
          tbox::Pointer<pdat::CellData<double> >
            cell_viscosity_ptr = patch->getPatchData(cell_viscosity_id);
          pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);
          if(d_dim.getValue()==2)
            {
              tbox::Pointer<pdat::NodeData<double> >
                edge_viscosity_ptr = patch->getPatchData(edge_viscosity_id);
              pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

              for(pdat::NodeIterator ni(edge_viscosity.getBox()); ni; ni++)
                {
                  pdat::NodeIndex e=ni();
                  pdat::CellIndex c(e);
                  cell_viscosity(c);
                  cell_viscosity(c-ip);
                  cell_viscosity(c-jp);
                  cell_viscosity(c-ip-jp);
                  edge_viscosity(e)=
                    pow(cell_viscosity(c)*cell_viscosity(c-ip)
                        *cell_viscosity(c-jp)*cell_viscosity(c-ip-jp),0.25);
                }
            }
          else
            {
              tbox::Pointer<pdat::EdgeData<double> >
                edge_viscosity_ptr = patch->getPatchData(edge_viscosity_id);
              pdat::EdgeData<double> &edge_viscosity(*edge_viscosity_ptr);
              for(int axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  hier::Box pbox=patch->getBox();
                  pbox.grow(axis,edge_viscosity.getGhostCellWidth()[axis]);

                  for(pdat::EdgeIterator ni(pbox,axis); ni; ni++)
                    {
                      pdat::EdgeIndex e=ni();
                      pdat::CellIndex c(e);
                      edge_viscosity(e)=
                        pow(cell_viscosity(c)*cell_viscosity(c-pp[axis2])
                            *cell_viscosity(c-pp[axis3])
                            *cell_viscosity(c-pp[axis2]-pp[axis3]),0.25);
                    }
                }
            }
        }

      /* Ghost fill */
      xfer::RefineAlgorithm refiner(d_dim);
      refiner.registerRefine(edge_viscosity_id,edge_viscosity_id,
                             edge_viscosity_id,
                             tbox::Pointer<xfer::RefineOperator>(0));

      tbox::Pointer<xfer::RefineSchedule> schedule=
        refiner.createSchedule(d_hierarchy->getPatchLevel(ln));
        
      schedule->fillData(0.0,false);
    }
}
