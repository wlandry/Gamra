/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FAC.hxx"
#include "Constants.hxx"

#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/pdat/NodeData.h>
#include <SAMRAI/pdat/EdgeData.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/pdat/NodeGeometry.h>

/// Fix the moduli on the coarse grids by coarsening from the finer
/// grids geometrically averaging the cell moduli to get the edge
/// moduli.

void Elastic::FAC::fix_moduli()
{
  const int level_max(hierarchy->getFinestLevelNumber());

  boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    cell_moduli_coarsen_operator;
  boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm>
    cell_moduli_coarsen_algorithm;
  std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    cell_moduli_coarsen_schedules (level_max + 1);

  SAMRAI::hier::VariableDatabase*
    vdb = SAMRAI::hier::VariableDatabase::getDatabase();
  boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> geometry =
    boost::dynamic_pointer_cast<SAMRAI::geom::CartesianGridGeometry>
    (hierarchy->getGridGeometry());
  boost::shared_ptr<SAMRAI::hier::Variable> variable;
  vdb->mapIndexToVariable(cell_moduli_id, variable);
  cell_moduli_coarsen_operator =
    geometry->lookupCoarsenOperator(variable,"CONSERVATIVE_COARSEN");

  if (!cell_moduli_coarsen_operator)
    { TBOX_ERROR("Elastic::FAC: Cannot find cell moduli coarsening operator"); }

  cell_moduli_coarsen_algorithm =
    boost::make_shared<SAMRAI::xfer::CoarsenAlgorithm >(dimension);
  cell_moduli_coarsen_algorithm->
    registerCoarsen(cell_moduli_id,cell_moduli_id,
                    cell_moduli_coarsen_operator);
  if(have_embedded_boundary())
    cell_moduli_coarsen_algorithm->
      registerCoarsen(level_set_id,level_set_id,
                      boost::shared_ptr<SAMRAI::hier::CoarsenOperator>());

  for (int dest_level = 0; dest_level < level_max; ++dest_level)
    {
      cell_moduli_coarsen_schedules[dest_level] =
        cell_moduli_coarsen_algorithm->
        createSchedule(hierarchy->getPatchLevel(dest_level),
                       hierarchy->getPatchLevel(dest_level + 1));
      if (!cell_moduli_coarsen_schedules[dest_level])
        { TBOX_ERROR("Elastic::FAC: Cannot create a coarsen schedule for cell "
                     "moduli restriction!\n"); }
    }

  for(int dest_level=level_max-1; dest_level>=0; --dest_level)
    {
      cell_moduli_coarsen_algorithm->
        resetSchedule(cell_moduli_coarsen_schedules[dest_level]);
      cell_moduli_coarsen_schedules[dest_level]->coarsenData();
    }

  cell_moduli_coarsen_algorithm.reset();

  /// Compute edge_moduli by averaging the cell moduli.

  SAMRAI::hier::Index ip(SAMRAI::hier::Index::getZeroIndex(dimension)),
    jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(dimension.getValue()>2)
    kp[2]=1;
  SAMRAI::hier::Index unit[]={ip,jp,kp};

  for (int level = 0; level <= hierarchy->getFinestLevelNumber(); ++level)
    {
      SAMRAI::hier::PatchLevel &patch_level = *hierarchy->getPatchLevel(level);
      
      for (SAMRAI::hier::PatchLevel::Iterator i_p(patch_level.begin());
           i_p!=patch_level.end(); ++i_p)
        {
          SAMRAI::hier::Patch &patch = **i_p;
          boost::shared_ptr<SAMRAI::pdat::CellData<double> > cell_moduli_ptr =
            boost::dynamic_pointer_cast<SAMRAI::pdat::CellData<double> >
            (patch.getPatchData(cell_moduli_id));
          SAMRAI::pdat::CellData<double> &cell_moduli(*cell_moduli_ptr);
          if(2==dimension.getValue())
            {
              boost::shared_ptr<SAMRAI::pdat::NodeData<double> >
                edge_moduli_ptr =
                boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
                (patch.getPatchData(edge_moduli_id));
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
                (patch.getPatchData(edge_moduli_id));
              SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);
              for(Gamra::Dir axis=0;axis<3;++axis)
                {
                  const int axis2((axis+1)%3), axis3((axis+2)%3);
                  SAMRAI::hier::Box pbox=patch.getBox();
                  /// Grow in axis direction only, because the
                  /// cell_moduli neighbors are not available on the
                  /// corners
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

      /// Ghost fill
      SAMRAI::xfer::RefineAlgorithm refiner;
      refiner.registerRefine(edge_moduli_id,edge_moduli_id,
                             edge_moduli_id,
                             boost::shared_ptr<SAMRAI::hier::RefineOperator>());

      boost::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule=
        refiner.createSchedule(hierarchy->getPatchLevel(level));
        
      schedule->fillData(0.0,false);
    }
}
