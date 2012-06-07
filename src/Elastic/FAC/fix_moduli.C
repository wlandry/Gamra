#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

/* Fix the moduli on the coarse grids by coarsening from the finer
   grids. */

void SAMRAI::Elastic::FAC::fix_moduli()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  tbox::Pointer<xfer::CoarsenOperator> cell_moduli_coarsen_operator;
  tbox::Pointer<xfer::CoarsenAlgorithm> cell_moduli_coarsen_algorithm;
  tbox::Array<tbox::Pointer<xfer::CoarsenSchedule> >
    cell_moduli_coarsen_schedules;

  hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
  tbox::Pointer<geom::CartesianGridGeometry> geometry =
    d_hierarchy->getGridGeometry();
  tbox::Pointer<hier::Variable> variable;
  vdb->mapIndexToVariable(cell_moduli_id, variable);
  cell_moduli_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,
                                    "CONSERVATIVE_COARSEN");
                                    // "CELL_VISCOSITY_COARSEN");

  if (!cell_moduli_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find cell moduli coarsening operator");
  }

  cell_moduli_coarsen_schedules.resizeArray(ln_max + 1);
  cell_moduli_coarsen_algorithm = new xfer::CoarsenAlgorithm(d_dim);
  cell_moduli_coarsen_algorithm->
    registerCoarsen(cell_moduli_id,cell_moduli_id,
                    cell_moduli_coarsen_operator);

  for (int dest_ln = 0; dest_ln < ln_max; ++dest_ln) {
    cell_moduli_coarsen_schedules[dest_ln] =
      cell_moduli_coarsen_algorithm->
      createSchedule(d_hierarchy->getPatchLevel(dest_ln),
                     d_hierarchy->getPatchLevel(dest_ln + 1));
    if (!cell_moduli_coarsen_schedules[dest_ln]) {
      TBOX_ERROR(d_object_name
                 << ": Cannot create a coarsen schedule for cell moduli restriction!\n");
    }
  }

  for(int dest_ln=ln_max-1; dest_ln>=0; --dest_ln)
    {
      xfer::CoarsenAlgorithm coarsener(d_dim);
      coarsener.registerCoarsen(cell_moduli_id, cell_moduli_id,
                                cell_moduli_coarsen_operator);
      coarsener.resetSchedule(cell_moduli_coarsen_schedules[dest_ln]);
      cell_moduli_coarsen_schedules[dest_ln]->coarsenData();
      cell_moduli_coarsen_algorithm->
        resetSchedule(cell_moduli_coarsen_schedules[dest_ln]);
    }

  cell_moduli_coarsen_algorithm.setNull();
  cell_moduli_coarsen_schedules.setNull();

  /* Compute edge_moduli by averaging the cell moduli. */

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
            cell_moduli_ptr = patch->getPatchData(cell_moduli_id);
          pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
          if(2==d_dim.getValue())
            {
              tbox::Pointer<pdat::NodeData<double> >
                edge_moduli_ptr = patch->getPatchData(edge_moduli_id);
              pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

              for(pdat::NodeIterator ni(edge_moduli.getBox()); ni; ni++)
                {
                  for (int m=0;m<2;++m)
                    {
                      pdat::NodeIndex e=ni();
                      pdat::CellIndex c(e);
                      cell_moduli(c,m);
                      cell_moduli(c-ip,m);
                      cell_moduli(c-jp,m);
                      cell_moduli(c-ip-jp,m);
                      edge_moduli(e,m)=
                        pow(cell_moduli(c,m)*cell_moduli(c-ip,m)
                            *cell_moduli(c-jp,m)*cell_moduli(c-ip-jp,m),0.25);
                    }
		}
            }
          else
            {
              tbox::Pointer<pdat::EdgeData<double> >
                edge_moduli_ptr = patch->getPatchData(edge_moduli_id);
              pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);
              for(int axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  hier::Box pbox=patch->getBox();
                  pbox.grow(axis,edge_moduli.getGhostCellWidth()[axis]);

                  for(pdat::EdgeIterator ni(pbox,axis); ni; ni++)
                    {
                      pdat::EdgeIndex e=ni();
                      pdat::CellIndex c(e);
		      for (int m=0;m<2;++m)
                        {
                          edge_moduli(e,m)=
                            pow(cell_moduli(c,m)*cell_moduli(c-pp[axis2],m)
                                *cell_moduli(c-pp[axis3],m)
                                *cell_moduli(c-pp[axis2]-pp[axis3],m),0.25);
                        }
                    }
                }
            }
        }

      /* Ghost fill */
      xfer::RefineAlgorithm refiner(d_dim);
      refiner.registerRefine(edge_moduli_id,edge_moduli_id,
                             edge_moduli_id,
                             tbox::Pointer<xfer::RefineOperator>(0));

      tbox::Pointer<xfer::RefineSchedule> schedule=
        refiner.createSchedule(d_hierarchy->getPatchLevel(ln));
        
      schedule->fillData(0.0,false);
    }
}
