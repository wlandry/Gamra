#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

/* Fix the moduli on the coarse grids by coarsening from the finer
   grids. */

void Elastic::FAC::fix_moduli()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
    cell_moduli_coarsen_operator;
  SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm>
    cell_moduli_coarsen_algorithm;
  SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule> >
    cell_moduli_coarsen_schedules;

  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();
  SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry> geometry =
    d_hierarchy->getGridGeometry();
  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable> variable;
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
  cell_moduli_coarsen_algorithm = new SAMRAI::xfer::CoarsenAlgorithm(d_dim);
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
      SAMRAI::xfer::CoarsenAlgorithm coarsener(d_dim);
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

  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(d_dim)),
    jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(d_dim.getValue()>2)
    kp[2]=1;
  SAMRAI::hier::Index pp[]={ip,jp,kp};

  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel>
        level = d_hierarchy->getPatchLevel(ln);
      SAMRAI::hier::PatchLevel::Iterator i_p(*level);
      for ( ; i_p; i_p++)
        {
          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch> patch = *i_p;
          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<double> >
            cell_moduli_ptr = patch->getPatchData(cell_moduli_id);
          SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
          if(2==d_dim.getValue())
            {
              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<double> >
                edge_moduli_ptr = patch->getPatchData(edge_moduli_id);
              SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

              for(SAMRAI::pdat::NodeIterator ni(edge_moduli.getBox()); ni; ni++)
                {
                  for (int m=0;m<2;++m)
                    {
                      SAMRAI::pdat::NodeIndex e=ni();
                      SAMRAI::pdat::CellIndex c(e);
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<double> >
                edge_moduli_ptr = patch->getPatchData(edge_moduli_id);
              SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);
              for(int axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  SAMRAI::hier::Box pbox=patch->getBox();
                  pbox.grow(axis,edge_moduli.getGhostCellWidth()[axis]);

                  for(SAMRAI::pdat::EdgeIterator ni(pbox,axis); ni; ni++)
                    {
                      SAMRAI::pdat::EdgeIndex e=ni();
                      SAMRAI::pdat::CellIndex c(e);
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
      SAMRAI::xfer::RefineAlgorithm refiner(d_dim);
      refiner.registerRefine(edge_moduli_id,edge_moduli_id,
                             edge_moduli_id,
                             SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>(0));

      SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule> schedule=
        refiner.createSchedule(d_hierarchy->getPatchLevel(ln));
        
      schedule->fillData(0.0,false);
    }
}
