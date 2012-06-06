/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/

#include "SAMRAI/pdat/CellVariable.h"
#include "Stokes/FACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
  namespace solv {


    /*
********************************************************************
* Set state from database                                          *
*                                                                  *
* Do not allow FAC preconditioner and Stokes FAC operators to be  *
* set from database, as that may cause them to be inconsistent     *
* with this object if user does not coordinate the inputs          *
* correctly.  This is also why we don't allow direct access to     *
* those objects.  The responsibility for maintaining consistency   *
* lies in the public functions to set parameters, so use them      *
* instead of setting the parameters directly in this function.     *
********************************************************************
*/

    void Stokes::FACSolver::getFromInput(tbox::Pointer<tbox::Database> database)
    {
      if (database) {
        if (database->isBool("enable_logging")) {
          bool logging = database->getBool("enable_logging");
          enableLogging(logging);
        }
        if (database->isInteger("max_cycles")) {
          int max_cycles = database->getInteger("max_cycles");
          setMaxCycles(max_cycles);
        }
        if (database->isDouble("residual_tol")) {
          double residual_tol = database->getDouble("residual_tol");
          setResidualTolerance(residual_tol);
        }
        if (database->isInteger("num_pre_sweeps")) {
          int num_pre_sweeps = database->getInteger("num_pre_sweeps");
          setPresmoothingSweeps(num_pre_sweeps);
        }
        if (database->isInteger("num_post_sweeps")) {
          int num_post_sweeps = database->getInteger("num_post_sweeps");
          setPostsmoothingSweeps(num_post_sweeps);
        }
        if (database->isString("coarse_fine_discretization")) {
          std::string s = database->getString("coarse_fine_discretization");
          setCoarseFineDiscretization(s);
        }
        if (database->isString("p_prolongation_method")) {
          std::string s = database->getString("p_prolongation_method");
          set_P_ProlongationMethod(s);
        }
        if (database->isString("v_prolongation_method")) {
          std::string s = database->getString("v_prolongation_method");
          set_V_ProlongationMethod(s);
        }
        if (database->isString("coarse_solver_choice")) {
          std::string s = database->getString("coarse_solver_choice");
          setCoarsestLevelSolverChoice(s);
        }
        if (database->isDouble("coarse_solver_tolerance")) {
          double tol = database->getDouble("coarse_solver_tolerance");
          setCoarsestLevelSolverTolerance(tol);
        }
        if (database->isInteger("coarse_solver_max_iterations")) {
          int itr = database->getInteger("coarse_solver_max_iterations");
          setCoarsestLevelSolverMaxIterations(itr);
        }
#ifdef HAVE_HYPRE
        if (database->isBool("use_smg")) {
          bool smg = database->getBool("use_smg");
          setUseSMG(smg);
        }
#endif
      }
    }


  }
}
