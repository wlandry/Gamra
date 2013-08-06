#include "SAMRAI/SAMRAI_config.h"

#include <string>

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "Stokes/P_Refine.h"
#include "Stokes/V_Refine.h"
#include "Stokes/P_Boundary_Refine.h"
#include "Stokes/V_Boundary_Refine.h"
#include "Stokes/V_Coarsen.h"
#include "Stokes/Resid_Coarsen.h"
#include "Elastic/V_Refine.h"
#include "Elastic/V_Boundary_Refine.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"

#include "Stokes/FAC.h"
#include "Elastic/FAC.h"

#include "solve_system.h"

int main(int argc, char* argv[])
{
  SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  SAMRAI::tbox::SAMRAIManager::initialize();
  SAMRAI::tbox::SAMRAIManager::startup();

  bool converged;
  /*
   * Create block to force pointer deallocation.  If this is not done
   * then there will be memory leaks reported.
   */
  {
    /*
     * Process command line arguments.  For each run, the input
     * filename must be specified.  Usage is:
     *
     *    executable <input file name>
     *
     */
    std::string input_filename;

    if (argc != 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                 << "  options:\n"
                 << "  none at this time\n");
    } else {
      input_filename = argv[1];
    }

    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
    SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);

    /*
     * Create input database and parse all data in input file.
     */

    boost::shared_ptr<SAMRAI::tbox::InputDatabase>
      input_db(new SAMRAI::tbox::InputDatabase("input_db"));
    SAMRAI::tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                             input_db);

    /*
     * Set up the timer manager.
     */
    if(input_db->isDatabase("TimerManager"))
      SAMRAI::tbox::TimerManager::createManager(input_db->getDatabase
                                                ("TimerManager"));

    /*
     * Retrieve "Main" section from input database.
     * The main database is used only in main().
     * The base_name variable is a base name for
     * all name strings in this program.
     */

    boost::shared_ptr<SAMRAI::tbox::Database>
      main_db(input_db->getDatabase("Main"));

    const SAMRAI::tbox::Dimension
      dim(static_cast<unsigned short>(main_db->getInteger("dim")));

    std::string base_name = "unnamed";
    base_name = main_db->getStringWithDefault("base_name", base_name);

    const std::string log_file_name = base_name + ".log";
    bool log_all_nodes(main_db->getBoolWithDefault("log_all_nodes",false));
                                                   
    if (log_all_nodes)
      {
        SAMRAI::tbox::PIO::logAllNodes(log_file_name);
      }
    else
      {
        SAMRAI::tbox::PIO::logOnlyNodeZero(log_file_name);
      }

    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry>
      grid_geometry(new SAMRAI::geom::CartesianGridGeometry
                    (dim, base_name + "CartesianGridGeometry",
                     input_db->getDatabase("CartesianGridGeometry")));

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy>
      patch_hierarchy(new SAMRAI::hier::PatchHierarchy
                      (base_name + "::PatchHierarchy",
                       grid_geometry,
                       input_db->getDatabase("PatchHierarchy")));

    if(input_db->isDatabase("Stokes"))
      {
        Stokes::FAC fac_stokes(base_name + "::Stokes::FAC", dim,
                               input_db->getDatabase("Stokes"));
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
            (fac_stokes.cell_viscosity_id)));

        converged=solve_system(fac_stokes,main_db,input_db,patch_hierarchy,
                               base_name,dim);
      }
    else
      {
        Elastic::FAC fac_elastic(base_name + "::Elastic::FAC", dim,
                                 input_db->getDatabase("Elastic"));

        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Elastic::V_Refine()));
        grid_geometry->addRefineOperator
          (typeid(SAMRAI::pdat::SideVariable<double>).name(),
           boost::shared_ptr<SAMRAI::hier::RefineOperator>
           (new Elastic::V_Boundary_Refine()));
        grid_geometry->addCoarsenOperator
          (typeid(SAMRAI::pdat::CellVariable<double>).name(),
           boost::make_shared<SAMRAI::geom::CartesianCellDoubleWeightedAverage>
           ());

        converged=solve_system(fac_elastic,main_db,input_db,patch_hierarchy,
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
