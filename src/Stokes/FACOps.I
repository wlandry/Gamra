/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for solving scalar Stokes using FAC 
 *
 ************************************************************************/
namespace Stokes {

  inline
  void FACOps::setPreconditioner(
                                 const SAMRAI::solv::FACPreconditioner* preconditioner) {
    d_preconditioner = preconditioner;
  }

#ifdef HAVE_HYPRE
  inline
  void FACOps::setUseSMG(
                         bool use_smg)
  {
    if (d_hierarchy) {
      TBOX_ERROR(
                 d_object_name << ": setUseSMG(bool) may NOT be called\n"
                 <<
                 "while the solver state is initialized, as that\n"
                 << "would lead to a corrupted solver state.\n");
    }
    d_hypre_solver.setUseSMG(use_smg);
  }
#endif

  /*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

  // inline
  // void FACOps::setPhysicalBcCoefObject(
  //    const RobinBcCoefStrategy* physical_bc_coef)
  // {
  //    d_physical_bc_coef = physical_bc_coef;
  //    // d_bc_helper.setCoefImplementation(physical_bc_coef);
  // #ifdef HAVE_HYPRE
  //    d_hypre_solver.setPhysicalBcCoefObject(d_physical_bc_coef);
  // #endif
  // }

  /*
********************************************************************
********************************************************************
*/

  inline
  void FACOps::enableLogging(
                             bool enable_logging)
  {
    d_enable_logging = enable_logging;
  }

  /*
********************************************************************
* Set the choice for smoothing algorithm.                          *
********************************************************************
*/

  inline
  void FACOps::setSmoothingChoice(
                                  const std::string& smoothing_choice)
  {
#ifdef DEBUG_CHECK_ASSERTIONS
    if (smoothing_choice != "Tackley" && smoothing_choice != "Gerya") {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                 << smoothing_choice
                 << "' in Stokes::FACOps::setSmoothingChoice.");
    }
#endif
    d_smoothing_choice = smoothing_choice;
  }

  /*
********************************************************************
* Set the choice for the coarse level solver.                      *
********************************************************************
*/

  inline
  void FACOps::setCoarsestLevelSolverChoice(
                                            const std::string& choice) {
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef HAVE_HYPRE
    if (choice == "hypre") {
      TBOX_ERROR(d_object_name << ": HYPRe library is not available.\n");
    }
#endif
#endif
    if (choice == "Tackley"
        || choice == "Gerya"
        || choice == "hypre") {
      d_coarse_solver_choice = choice;
    } else {
      TBOX_ERROR(
                 d_object_name << ": Bad coarse level solver choice '"
                 << choice
                 <<
                 "' in Stokes::FACOps::setCoarseLevelSolverChoice.");
    }
  }

  /*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

  inline
  void FACOps::setCoarsestLevelSolverTolerance(
                                               double tol) {
    d_coarse_solver_tolerance = tol;
  }

  /*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

  inline
  void FACOps::setCoarsestLevelSolverMaxIterations(
                                                   int max_iterations) {
#ifdef DEBUG_CHECK_ASSERTIONS
    if (max_iterations < 0) {
      TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
    }
#endif
    d_coarse_solver_max_iterations = max_iterations;
  }

  /*
********************************************************************
* Set the prolongation method                                      *
********************************************************************
*/

  inline
  void FACOps::set_P_ProlongationMethod(
                                        const std::string& prolongation_method) {
#ifdef DEBUG_CHECK_ASSERTIONS
    if (d_hierarchy) {
      TBOX_ERROR(
                 d_object_name << ": Cannot change p prolongation method\n"
                 <<
                 "while operator state is initialized because that\n"
                 << "causes a corruption in the state.\n");
    }
#endif
    p_prolongation_method = prolongation_method;
  }
}
