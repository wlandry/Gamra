/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Main program
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <string>
using namespace std;

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "Stokes/P_Refine.h"
#include "Stokes/V_Refine.h"
#include "Stokes/P_Boundary_Refine.h"
#include "Stokes/V_Boundary_Refine.h"
#include "Stokes/V_Coarsen.h"
#include "Stokes/Resid_Coarsen.h"
#include "Elastic/V_Refine.h"
#include "Elastic/V_Boundary_Refine.h"
#include "Elastic/V_Coarsen.h"
#include "Elastic/Resid_Coarsen.h"

#include "Stokes/FAC.h"
#include "Elastic/FAC.h"

#include "solve_system.h"

using namespace SAMRAI;

/*
************************************************************************
*                                                                      *
* This is the driver program to demonstrate                            *
* how to use the Stokes and Elastic FAC solver.                                   *
*                                                                      *
* We set up the simple problem                                         *
*          u + div(grad(u)) = sin(x)*sin(y)                            *
* in the domain [0:1]x[0:1], with u=0 on the                           *
* boundary.                                                            *
*                                                                      *
* Stokes::FAC and Elastic::FAC are the primary objects used to         *
* set up and solve the system.  It maintains                           *
* the data for the computed solution u, the                            *
* exact solution, and the right hand side.                             *
*                                                                      *
* The hierarchy created to solve this problem                          *
* has only one level.  (The Stokes::FAC and Elastic::FAC solver        *
* is a single-level solver.)                                           *
*                                                                      *
*************************************************************************
*/

int main(
         int argc,
         char* argv[])
{
  /*
   * Initialize MPI, SAMRAI.
   */

  tbox::SAMRAI_MPI::init(&argc, &argv);
  tbox::SAMRAIManager::initialize();
  tbox::SAMRAIManager::startup();

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
    string input_filename;

    if (argc != 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                 << "  options:\n"
                 << "  none at this time" << endl);
    } else {
      input_filename = argv[1];
    }

    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit(true);
    tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort(true);

    /*
     * Create input database and parse all data in input file.
     */

    tbox::Pointer<tbox::Database> input_db(new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

    /*
     * Set up the timer manager.
     */
    if (input_db->isDatabase("TimerManager")) {
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    }

    /*
     * Retrieve "Main" section from input database.
     * The main database is used only in main().
     * The base_name variable is a base name for
     * all name strings in this program.
     */

    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

    string base_name = "unnamed";
    base_name = main_db->getStringWithDefault("base_name", base_name);

    /*
     * Start logging.
     */
    const string log_file_name = base_name + ".log";
    bool log_all_nodes = false;
    log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
                                                log_all_nodes);
    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
    } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
    }

    /*
     * Create major algorithm and data objects which comprise application.
     * Each object will be initialized either from input data or restart
     * files, or a combination of both.  Refer to each class constructor
     * for details.  For more information on the composition of objects
     * for this application, see comments at top of file.
     */

    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry>
      grid_geometry(new SAMRAI::geom::CartesianGridGeometry
                    (dim, base_name + "CartesianGridGeometry",
                     input_db->getDatabase("CartesianGridGeometry")));

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy>
      patch_hierarchy(new SAMRAI::hier::PatchHierarchy
                      (base_name + "::PatchHierarchy",
                       grid_geometry,
                       input_db->getDatabase("PatchHierarchy")));

    /*
     * The Stokes::FAC and Elastic::FAC objects is the main user
     * object specific to the problem being solved.  It provides the
     * implementations for setting up the grid and plotting data.  It
     * also wraps up the solve process that includes making the
     * initial guess, specifying the boundary conditions and call the
     * solver.
     */

    if(input_db->isDatabase("Stokes"))
      {
        Stokes::FAC fac_stokes(base_name + "::Stokes::FAC", dim,
                               input_db->getDatabase("Stokes"));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new SAMRAI::geom::Stokes::P_Refine(dim)));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new SAMRAI::geom::Stokes::V_Refine(dim)));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new SAMRAI::geom::Stokes::P_Boundary_Refine(dim)));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new SAMRAI::geom::Stokes::V_Boundary_Refine(dim)));
        grid_geometry->addSpatialCoarsenOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
           (new SAMRAI::geom::Stokes::V_Coarsen(dim)));
        grid_geometry->addSpatialCoarsenOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
           (new SAMRAI::geom::Stokes::Resid_Coarsen(dim,fac_stokes.cell_viscosity_id)));

        solve_system(fac_stokes,main_db,input_db,patch_hierarchy,
                     base_name,dim);
      }
    else
      {
        Elastic::FAC fac_elastic(base_name + "::Elastic::FAC", dim,
                                 input_db->getDatabase("Elastic"));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new Elastic::V_Refine(dim)));
        grid_geometry->addSpatialRefineOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
           (new Elastic::V_Boundary_Refine(dim)));
        grid_geometry->addSpatialCoarsenOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
           (new Elastic::V_Coarsen(dim)));
        grid_geometry->addSpatialCoarsenOperator
          (SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
           (new Elastic::Resid_Coarsen(dim,fac_elastic.cell_moduli_id)));

        solve_system(fac_elastic,main_db,input_db,patch_hierarchy,
                     base_name,dim);
      }

    /*
     * Deallocate objects when done.
     */

  }
  /*
   * This print is for the SAMRAI testing framework.  Passing here
   * means application ran.  A better test would actually test the
   * results.
   */
  tbox::pout << "\nPASSED:  FAC" << endl;

  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();

  return 0;
}
