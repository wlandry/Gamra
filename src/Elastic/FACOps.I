/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for solving scalar Elastic using FAC 
 *
 ************************************************************************/

inline
void Elastic::FACOps::setPreconditioner(
   const SAMRAI::solv::FACPreconditioner* preconditioner) {
   d_preconditioner = preconditioner;
}

/*
 ********************************************************************
 * Set the physical boundary condition object.                      *
 ********************************************************************
 */

// inline
// void Elastic::FACOps::setPhysicalBcCoefObject(
//    const RobinBcCoefStrategy* physical_bc_coef)
// {
//    d_physical_bc_coef = physical_bc_coef;
//    // d_bc_helper.setCoefImplementation(physical_bc_coef);
// }

/*
 ********************************************************************
 ********************************************************************
 */

inline
void Elastic::FACOps::enableLogging(
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
void Elastic::FACOps::setSmoothingChoice(
   const std::string& smoothing_choice)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (smoothing_choice != "Tackley" && smoothing_choice != "Gerya") {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                               << smoothing_choice
                               << "' in Elastic::FACOps::setSmoothingChoice.");
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
void Elastic::FACOps::setCoarsestLevelSolverChoice(
   const std::string& choice) {
   if (choice == "Tackley") {
      d_coarse_solver_choice = choice;
   } else {
      TBOX_ERROR(
         d_object_name << ": Bad coarse level solver choice '"
         << choice
         <<
         "' in Elastic::FACOps::setCoarseLevelSolver.");
   }
}

/*
 ********************************************************************
 * Set the tolerance for the coarse level solver.                   *
 ********************************************************************
 */

inline
void Elastic::FACOps::setCoarsestLevelSolverTolerance(
   double tol) {
   d_coarse_solver_tolerance = tol;
}

/*
 ********************************************************************
 * Set the tolerance for the coarse level solver.                   *
 ********************************************************************
 */

inline
void Elastic::FACOps::setCoarsestLevelSolverMaxIterations(
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
 * Set the coarse-fine discretization method.                       *
 ********************************************************************
 */

inline
void Elastic::FACOps::setCoarseFineDiscretization(
   const std::string& coarsefine_method) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_hierarchy) {
      TBOX_ERROR(
         d_object_name << ": Cannot change coarse-fine\n"
         <<
         "discretization method while operator state\n"
         << "is initialized because that causes a\n"
         << "corruption in the state.\n");
   }
#endif
   d_cf_discretization = coarsefine_method;
}

/*
 ********************************************************************
 * Set the prolongation method                                      *
 ********************************************************************
 */

inline
void Elastic::FACOps::set_V_ProlongationMethod(
   const std::string& prolongation_method) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_hierarchy) {
      TBOX_ERROR(
         d_object_name << ": Cannot change v prolongation method\n"
         <<
         "while operator state is initialized because that\n"
         << "causes a corruption in the state.\n");
   }
#endif
   v_prolongation_method = prolongation_method;
}

