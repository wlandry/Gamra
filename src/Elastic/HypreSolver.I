/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Level solver for diffusion-like elliptic problems. 
 *
 ************************************************************************/
namespace SAMRAI {
namespace solv {

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setSolnIdDepth(
   const int depth) {
   d_soln_depth = depth;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setRhsIdDepth(
   const int depth) {
   d_rhs_depth = depth;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setUseSMG(
   bool use_smg) {
   d_use_smg = use_smg;
}

/*
 ********************************************************************
 * Specify bc using the default internal bc coefficient object.     *
 * Clear up data supporting external bc coefficient setter.         *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setBoundaries(
   const std::string& boundary_type,
   const int fluxes,
   const int flags,
   int* bdry_types)
{
   d_physical_bc_simple_case.setBoundaries(boundary_type,
      fluxes,
      flags,
      bdry_types);
   d_physical_bc_coef_strategy = &d_physical_bc_simple_case;
   d_physical_bc_variable.setNull();
}

/*
 ********************************************************************
 * Set the physical boundary condition object.                      *
 ********************************************************************
 */

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setPhysicalBcCoefObject(
   const RobinBcCoefStrategy* physical_bc_coef_strategy,
   const tbox::Pointer<hier::Variable> variable)
{
   d_physical_bc_coef_strategy = physical_bc_coef_strategy;
   d_physical_bc_variable = variable;
}

SAMRAI_INLINE_KEYWORD
int Elastic::HypreSolver::getNumberOfIterations() const
{
   return d_number_iterations;
}

SAMRAI_INLINE_KEYWORD
double Elastic::HypreSolver::getRelativeResidualNorm() const
{
   return d_relative_residual_norm;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setNumPreRelaxSteps(
   const int steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
#endif
   d_num_pre_relax_steps = steps;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setNumPostRelaxSteps(
   const int steps)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
#endif
   d_num_post_relax_steps = steps;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setPrintSolverInfo(
   const bool print)
{
   d_print_solver_info = print;
}

SAMRAI_INLINE_KEYWORD
void Elastic::HypreSolver::setStoppingCriteria(
   const int max_iterations,
   const double residual_tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_iterations >= 0);
   TBOX_ASSERT(residual_tol >= 0.0);
#endif
   d_max_iterations = max_iterations;
   d_relative_residual_tol = residual_tol;
}

}
}
