#include "Elastic/FAC.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/NodeGeometry.h"

/* Fix the moduli on the coarse grids by coarsening from the finer
   grids geometrically averaging the cell moduli to get the edge
   moduli. */

void Elastic::FAC::fix_moduli()
{
  const int ln_max(d_hierarchy->getFinestLevelNumber());

  boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    cell_moduli_coarsen_operator;
  boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
    cell_moduli_coarsen_algorithm;
  SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    cell_moduli_coarsen_schedules;

  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (d_hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;
  vdb->mapIndexToVariable(cell_moduli_id, variable);
  cell_moduli_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,"CONSERVATIVE_COARSEN");

  if (!cell_moduli_coarsen_operator) {
    TBOX_ERROR(d_object_name
               << ": Cannot find cell moduli coarsening operator");
  }

  cell_moduli_coarsen_schedules.resizeArray(ln_max + 1);
  cell_moduli_coarsen_algorithm =
    boost::make_shared<SAMRAI::xfer::CoarsenAlgorithm >(d_dim,true);
  cell_moduli_coarsen_algorithm->
    registerCoarsen(cell_moduli_id,cell_moduli_id,
                    cell_moduli_coarsen_operator);
  if(have_embedded_boundary())
    cell_moduli_coarsen_algorithm->
      registerCoarsen(level_set_id,level_set_id,
                      boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());

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
      cell_moduli_coarsen_algorithm->
        resetSchedule(cell_moduli_coarsen_schedules[dest_ln]);
      cell_moduli_coarsen_schedules[dest_ln]->coarsenData();
    }

  cell_moduli_coarsen_algorithm.reset();
  cell_moduli_coarsen_schedules.setNull();

  /* Compute edge_moduli by averaging the cell moduli. */

  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(d_dim)),
    jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(d_dim.getValue()>2)
    kp[2]=1;
  SAMRAI::hier::Index unit[]={ip,jp,kp};

  for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
      boost::shared_ptr<SAMRAI::hier::PatchLevel>
        level = d_hierarchy->getPatchLevel(ln);
      
      for (SAMRAI::hier::PatchLevel::Iterator i_p(level->begin());
           i_p!=level->end(); ++i_p)
        {
          boost::shared_ptr<SAMRAI::hier::Patch> patch = *i_p;
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch->getPatchData(cell_moduli_id));
          SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
          if(2==d_dim.getValue())
            {
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> >
                edge_moduli_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch->getPatchData(edge_moduli_id));
              SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

              SAMRAI::pdat::NodeIterator
                nend(SAMRAI::pdat::NodeGeometry::end(edge_moduli.getBox()));
              for(SAMRAI::pdat::NodeIterator
                    ni(SAMRAI::pdat::NodeGeometry::begin(edge_moduli.getBox()));
                  ni!=nend; ni++)
                {
                  for (int m=0;m<2;++m)
                    {
                      const SAMRAI::pdat::NodeIndex &e(*ni);
                      SAMRAI::pdat::CellIndex c(e);
                      edge_moduli(e,m)=
                        pow(cell_moduli(c,m)*cell_moduli(c-ip,m)
                            *cell_moduli(c-jp,m)*cell_moduli(c-ip-jp,m),0.25);
                    }
		}
            }
          else
            {
              boost::shared_ptr<SAMRAI::pdat::EdgeData<double> >
                edge_moduli_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
                (patch->getPatchData(edge_moduli_id));
              SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);
              for(int axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  SAMRAI::hier::Box pbox=patch->getBox();
                  /* Grow in axis direction only, because the
                     cell_moduli neighbors are not available on the
                     corners */
                  pbox.grow(axis,edge_moduli.getGhostCellWidth()[axis]);

                  SAMRAI::pdat::EdgeIterator
                    nend(SAMRAI::pdat::EdgeGeometry::end(pbox,axis));
                  for(SAMRAI::pdat::EdgeIterator
                        ni(SAMRAI::pdat::EdgeGeometry::begin(pbox,axis));
                      ni!=nend; ++ni)
                    {
                      const SAMRAI::pdat::EdgeIndex &e(*ni);
                      SAMRAI::pdat::CellIndex c(e);
		      for (int m=0;m<2;++m)
                        {
                          edge_moduli(e,m)=
                            pow(cell_moduli(c,m)*cell_moduli(c-unit[axis2],m)
                                *cell_moduli(c-unit[axis3],m)
                                *cell_moduli(c-unit[axis2]-unit[axis3],m),0.25);
                        }
                    }
                }
            }
        }

      /* Ghost fill */
      SAMRAI::xfer::RefineAlgorithm refiner;
      refiner.registerRefine(edge_moduli_id,edge_moduli_id,
                             edge_moduli_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());

      boost::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule=
        refiner.createSchedule(d_hierarchy->getPatchLevel(ln));
        
      schedule->fillData(0.0,false);
    }
}
