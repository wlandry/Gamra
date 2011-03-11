#include "FACStokes.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

/* Fix the viscosity on the coarse grids by coarsening from the finer
   grids. */

void SAMRAI::FACStokes::fix_viscosity()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  tbox::Pointer<xfer::CoarsenOperator> cell_viscosity_coarsen_operator;
  tbox::Pointer<xfer::CoarsenAlgorithm> cell_viscosity_coarsen_algorithm;
  tbox::Array<tbox::Pointer<xfer::CoarsenSchedule> >
    cell_viscosity_coarsen_schedules;

  tbox::Pointer<xfer::CoarsenOperator> edge_viscosity_coarsen_operator;
  tbox::Pointer<xfer::CoarsenAlgorithm> edge_viscosity_coarsen_algorithm;
  tbox::Array<tbox::Pointer<xfer::CoarsenSchedule> >
    edge_viscosity_coarsen_schedules;

  hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
  tbox::Pointer<geom::CartesianGridGeometry> geometry =
    d_hierarchy->getGridGeometry();
  tbox::Pointer<hier::Variable> variable;
  vdb->mapIndexToVariable(cell_viscosity_id, variable);
  cell_viscosity_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CELL_VISCOSITY_COARSEN");

  vdb->mapIndexToVariable(edge_viscosity_id, variable);
  edge_viscosity_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "EDGE_VISCOSITY_COARSEN");

  if (!cell_viscosity_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find cell viscosity coarsening operator");
  }

  if (!edge_viscosity_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find edge viscosity coarsening operator");
  }

  cell_viscosity_coarsen_schedules.resizeArray(ln_max + 1);
  cell_viscosity_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  cell_viscosity_coarsen_algorithm->
    registerCoarsen(cell_viscosity_id,cell_viscosity_id,
                    cell_viscosity_coarsen_operator);

  edge_viscosity_coarsen_schedules.resizeArray(ln_max + 1);
  edge_viscosity_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  edge_viscosity_coarsen_algorithm->
    registerCoarsen(edge_viscosity_id,edge_viscosity_id,
                    edge_viscosity_coarsen_operator);

  for (int dest_ln = 0; dest_ln < ln_max; ++dest_ln) {
    cell_viscosity_coarsen_schedules[dest_ln] =
      cell_viscosity_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!cell_viscosity_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for cell viscosity restriction!\n");
    }
    edge_viscosity_coarsen_schedules[dest_ln] =
      edge_viscosity_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!edge_viscosity_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for edge viscosity restriction!\n");
    }
  }

  for(int dest_ln=ln_max-1; dest_ln>=0; --dest_ln)
    {
      {
        xfer::CoarsenAlgorithm coarsener(d_dim);
        coarsener.registerCoarsen(cell_viscosity_id, cell_viscosity_id,
                                  cell_viscosity_coarsen_operator);
        coarsener.resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
        cell_viscosity_coarsen_schedules[dest_ln]->coarsenData();
        cell_viscosity_coarsen_algorithm->
          resetSchedule(cell_viscosity_coarsen_schedules[dest_ln]);
      }
      {
        xfer::CoarsenAlgorithm coarsener(d_dim);
        coarsener.registerCoarsen(edge_viscosity_id, edge_viscosity_id,
                                  edge_viscosity_coarsen_operator);
        coarsener.resetSchedule(edge_viscosity_coarsen_schedules[dest_ln]);
        edge_viscosity_coarsen_schedules[dest_ln]->coarsenData();
        edge_viscosity_coarsen_algorithm->
          resetSchedule(edge_viscosity_coarsen_schedules[dest_ln]);
      }
    }

  cell_viscosity_coarsen_algorithm.setNull();
  cell_viscosity_coarsen_schedules.setNull();

  edge_viscosity_coarsen_algorithm.setNull();
  edge_viscosity_coarsen_schedules.setNull();
}
