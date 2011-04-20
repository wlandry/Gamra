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
                                    "CONSERVATIVE_COARSEN");
                                    // "CELL_VISCOSITY_COARSEN");

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

  hier::Index ip(hier::Index::getZeroIndex(d_dim)), jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(d_dim.getValue()>2)
    kp[2]=1;

  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
    tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
    hier::PatchLevel::Iterator i_p(*level);
    for ( ; i_p; i_p++) {
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
                1/(1/cell_viscosity(c)
                   + 1/cell_viscosity(c-ip)
                   + 1/cell_viscosity(c-jp)
                   + 1/cell_viscosity(c-ip-jp));
            }
        }
      else
        {
          tbox::Pointer<pdat::EdgeData<double> >
            edge_viscosity_ptr = patch->getPatchData(edge_viscosity_id);
          pdat::EdgeData<double> &edge_viscosity(*edge_viscosity_ptr);
          for(int axis=0;axis<3;++axis)
            {
              for(pdat::EdgeIterator ni(edge_viscosity.getBox(),axis); ni; ni++)
                {
                  pdat::EdgeIndex e=ni();
                  pdat::CellIndex c(e);
                  edge_viscosity(e)=
                    1/(1/cell_viscosity(c)
                       + 1/cell_viscosity(c-ip)
                       + 1/cell_viscosity(c-jp)
                       + 1/cell_viscosity(c-ip-jp)
                       + 1/cell_viscosity(c-kp)
                       + 1/cell_viscosity(c-ip-kp)
                       + 1/cell_viscosity(c-jp-kp)
                       + 1/cell_viscosity(c-ip-jp-kp));
                }
            }
        }
    }
  }


}
