/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar elastic equation. 
 *
 ************************************************************************/

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::set_V_ProlongationMethod(
   const std::string& v_prolongation_method)
{
   d_fac_ops.set_V_ProlongationMethod(v_prolongation_method);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setCoarsestLevelSolverChoice(
   const std::string& choice)
{
   d_fac_ops.setCoarsestLevelSolverChoice(choice);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setCoarsestLevelSolverTolerance(
   double tol)
{
   d_fac_ops.setCoarsestLevelSolverTolerance(tol);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setCoarsestLevelSolverMaxIterations(
   int max_iterations)
{
   d_fac_ops.setCoarsestLevelSolverMaxIterations(max_iterations);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setCoarseFineDiscretization(
   const std::string& coarsefine_method)
{
   d_fac_ops.setCoarseFineDiscretization(coarsefine_method);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setPresmoothingSweeps(
   int num_pre_sweeps) {
   d_fac_precond.setPresmoothingSweeps(num_pre_sweeps);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setPostsmoothingSweeps(
   int num_post_sweeps) {
   d_fac_precond.setPostsmoothingSweeps(num_post_sweeps);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setMaxCycles(
   int max_cycles) {
   d_fac_precond.setMaxCycles(max_cycles);
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::setResidualTolerance(
   double residual_tol) {
   d_fac_precond.setResidualTolerance(residual_tol);
}

SAMRAI_INLINE_KEYWORD
int Elastic::FACSolver::getNumberOfIterations() const
{
   return d_fac_precond.getNumberOfIterations();
}

SAMRAI_INLINE_KEYWORD
double Elastic::FACSolver::getResidualNorm() const
{
   return d_fac_precond.getResidualNorm();
}

SAMRAI_INLINE_KEYWORD
void Elastic::FACSolver::getConvergenceFactors(
   double& avg_factor,
   double& final_factor)
const
{
   d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
}

