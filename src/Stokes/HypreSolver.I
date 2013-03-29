/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Level solver for diffusion-like elliptic problems. 
 *
 ************************************************************************/
namespace Stokes {

  inline
  void HypreSolver::setSolnIdDepth(
                                   const int depth) {
    d_soln_depth = depth;
  }

  inline
  void HypreSolver::setRhsIdDepth(
                                  const int depth) {
    d_rhs_depth = depth;
  }

  inline
  void HypreSolver::setUseSMG(
                              bool use_smg) {
    d_use_smg = use_smg;
  }

  /*
********************************************************************
* Specify bc using the default internal bc coefficient object.     *
* Clear up data supporting external bc coefficient setter.         *
********************************************************************
*/

  inline
  void HypreSolver::setBoundaries(
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
    d_physical_bc_variable.reset();
  }

  /*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

  inline
  void HypreSolver::setPhysicalBcCoefObject(
                                            const RobinBcCoefStrategy* physical_bc_coef_strategy,
                                            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable> variable)
  {
    d_physical_bc_coef_strategy = physical_bc_coef_strategy;
    d_physical_bc_variable = variable;
  }

  inline
  int HypreSolver::getNumberOfIterations() const
  {
    return d_number_iterations;
  }

  inline
  double HypreSolver::getRelativeResidualNorm() const
  {
    return d_relative_residual_norm;
  }

  inline
  void HypreSolver::setNumPreRelaxSteps(
                                        const int steps)
  {
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy);
#endif
    d_num_pre_relax_steps = steps;
  }

  inline
  void HypreSolver::setNumPostRelaxSteps(
                                         const int steps)
  {
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy);
#endif
    d_num_post_relax_steps = steps;
  }

  inline
  void HypreSolver::setPrintSolverInfo(
                                       const bool print)
  {
    d_print_solver_info = print;
  }

  inline
  void HypreSolver::setStoppingCriteria(
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
