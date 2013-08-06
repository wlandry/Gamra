/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar elastic equation. 
 *
 ************************************************************************/

inline
void Elastic::FACSolver::set_V_ProlongationMethod(
   const std::string& v_prolongation_method)
{
   d_fac_ops->set_V_ProlongationMethod(v_prolongation_method);
}

inline
void Elastic::FACSolver::setCoarsestLevelSolverTolerance(
   double tol)
{
   d_fac_ops->setCoarsestLevelSolverTolerance(tol);
}

inline
void Elastic::FACSolver::setCoarsestLevelSolverMaxIterations(
   int max_iterations)
{
   d_fac_ops->setCoarsestLevelSolverMaxIterations(max_iterations);
}

inline
void Elastic::FACSolver::setCoarseFineDiscretization(
   const std::string& coarsefine_method)
{
   d_fac_ops->setCoarseFineDiscretization(coarsefine_method);
}

inline
int Elastic::FACSolver::getNumberOfIterations() const
{
   return d_fac_precond.getNumberOfIterations();
}

inline
double Elastic::FACSolver::getResidualNorm() const
{
   return d_fac_precond.getResidualNorm();
}

inline
void Elastic::FACSolver::getConvergenceFactors(
   double& avg_factor,
   double& final_factor)
const
{
   d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
}

