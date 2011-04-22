/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for solving scalar Stokes using FAC 
 *
 ************************************************************************/
namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setPreconditioner(
   const FACPreconditioner* preconditioner) {
   d_preconditioner = preconditioner;
}

#ifdef HAVE_HYPRE
SAMRAI_INLINE_KEYWORD
void StokesFACOps::setUseSMG(
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

// SAMRAI_INLINE_KEYWORD
// void StokesFACOps::setPhysicalBcCoefObject(
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

SAMRAI_INLINE_KEYWORD
void StokesFACOps::enableLogging(
   bool enable_logging)
{
   d_enable_logging = enable_logging;
}

/*
 ********************************************************************
 * Set the choice for smoothing algorithm.                          *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setSmoothingChoice(
   const std::string& smoothing_choice)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (smoothing_choice != "Tackley" && smoothing_choice != "Gerya") {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                               << smoothing_choice
                               << "' in StokesFACOps::setSmoothingChoice.");
   }
#endif
   d_smoothing_choice = smoothing_choice;
}

/*
 ********************************************************************
 * Set the choice for the coarse level solver.                      *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setCoarsestLevelSolverChoice(
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
         "' in scapStokesOpsX::setCoarseLevelSolver.");
   }
}

/*
 ********************************************************************
 * Set the tolerance for the coarse level solver.                   *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setCoarsestLevelSolverTolerance(
   double tol) {
   d_coarse_solver_tolerance = tol;
}

/*
 ********************************************************************
 * Set the tolerance for the coarse level solver.                   *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setCoarsestLevelSolverMaxIterations(
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

SAMRAI_INLINE_KEYWORD
void StokesFACOps::setCoarseFineDiscretization(
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

SAMRAI_INLINE_KEYWORD
void StokesFACOps::set_P_ProlongationMethod(
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

SAMRAI_INLINE_KEYWORD
void StokesFACOps::set_V_ProlongationMethod(
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

}
}
