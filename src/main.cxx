/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/tbox/InputManager.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>

#include "Stokes/FAC.hxx"
#include "Stokes/P_Refine.hxx"
#include "Stokes/V_Refine.hxx"
#include "Stokes/P_Boundary_Refine.hxx"
#include "Stokes/V_Boundary_Refine.hxx"
#include "Stokes/V_Coarsen.hxx"
#include "Stokes/Resid_Coarsen.hxx"
#include "Elastic/FAC.hxx"

#include "solve_system.hxx"

int main(int argc, char* argv[])
{
  SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  SAMRAI::tbox::SAMRAIManager::initialize();
  SAMRAI::tbox::SAMRAIManager::startup();

  bool converged;
  /// Create block to force pointer deallocation.  If this is not done
  /// then there will be memory leaks reported.
  // FIXME: Really???
  {
    std::string input_filename;
    if (argc != 2)
      {
        TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n");
      }
    else
      {
        input_filename = argv[1];
      }

    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);

    boost::shared_ptr<SAMRAI::tbox::InputDatabase>
      input_db(SAMRAI::tbox::InputManager::getManager()->parseInputFile
               (input_filename));
    if(input_db->isDatabase("TimerManager"))
      SAMRAI::tbox::TimerManager::createManager(input_db->getDatabase
                                                ("TimerManager"));
    boost::shared_ptr<SAMRAI::tbox::Database>
      main_db(input_db->getDatabase("Main"));
    const SAMRAI::tbox::Dimension
      dim(static_cast<unsigned short>(main_db->getInteger("dim")));

    std::string base_name(main_db->getString("base_name"));

    const std::string log_file_name = base_name + ".log";
    bool log_all_nodes(main_db->getBoolWithDefault("log_all_nodes",false));
                                                   
    if (log_all_nodes)
      { SAMRAI::tbox::PIO::logAllNodes(log_file_name); }
    else
      { SAMRAI::tbox::PIO::logOnlyNodeZero(log_file_name); }

    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry>
      grid_geometry(new SAMRAI::geom::CartesianGridGeometry
                    (dim, "CartesianGridGeometry",
                     input_db->getDatabase("CartesianGridGeometry")));

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
      patch_hierarchy(new SAMRAI::hier::PatchHierarchy
                      ("PatchHierarchy", grid_geometry,
                       input_db->getDatabase("PatchHierarchy")));

    if(input_db->isDatabase("Stokes"))
      {
        Stokes::FAC stokes(dim, *input_db->getDatabase("Stokes"));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::CellVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Stokes::P_Refine()));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Stokes::V_Refine()));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::CellVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Stokes::P_Boundary_Refine()));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Stokes::V_Boundary_Refine()));
        grid_geometry->addCoarsenOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
           (new Stokes::V_Coarsen()));
        grid_geometry->addCoarsenOperator
          (typeid(SAMRAI::pdat::CellVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
           (new Stokes::Resid_Coarsen
            (stokes.cell_viscosity_id)));

        converged=solve_system(stokes,main_db,input_db,patch_hierarchy,
                               base_name,dim);
      }
    else
      {
        Elastic::FAC elastic(dim, *input_db->getDatabase("Elastic"));

        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Elastic::V_Refine()));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Elastic::Coarse_Fine_Boundary_Refine()));
        grid_geometry->addCoarsenOperator
          (typeid(SAMRAI::pdat::CellVariable<double>).name(),
           boost::make_shared<SAMRAI::geom::CartesianCellDoubleWeightedAverage>
           ());

        converged=solve_system(elastic,main_db,input_db,patch_hierarchy,
                               base_name,dim);
      }
  }
  if(converged)
    SAMRAI::tbox::pout << "PASSED\n";
  else
    SAMRAI::tbox::pout << "FAILED\n";


  SAMRAI::tbox::SAMRAIManager::shutdown();
  SAMRAI::tbox::SAMRAIManager::finalize();
  SAMRAI::tbox::SAMRAI_MPI::finalize();

  return (converged ? 0 : 1);
}
