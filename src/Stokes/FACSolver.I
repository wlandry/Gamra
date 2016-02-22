/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar stokes equation. 
 *
 ************************************************************************/
namespace Stokes {
  inline
  void FACSolver::setCoarsestLevelSolverChoice(
                                               const std::string& choice)
  {
    d_fac_ops->setCoarsestLevelSolverChoice(choice);
  }

  inline
  void FACSolver::setCoarsestLevelSolverTolerance(
                                                  double tol)
  {
    d_fac_ops->setCoarsestLevelSolverTolerance(tol);
  }

  inline
  void FACSolver::setCoarsestLevelSolverMaxIterations(
                                                      int max_iterations)
  {
    d_fac_ops->setCoarsestLevelSolverMaxIterations(max_iterations);
  }

#ifdef HAVE_HYPRE
  inline
  void FACSolver::setUseSMG(
                            bool use_smg)
  {
    if (d_solver_is_initialized) {
      TBOX_ERROR(
                 d_object_name << ": setUseSMG(bool) may NOT be called\n"
                 <<
                 "while the solver state is initialized, as that\n"
                 << "would lead to a corrupted solver state.\n");
    }
    d_fac_ops->setUseSMG(use_smg);
  }
#endif

  inline
  int FACSolver::getNumberOfIterations() const
  {
    return d_fac_precond.getNumberOfIterations();
  }

  inline
  double FACSolver::getResidualNorm() const
  {
    return d_fac_precond.getResidualNorm();
  }

  inline
  void FACSolver::getConvergenceFactors(
                                        double& avg_factor,
                                        double& final_factor)
    const
  {
    d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
  }

}
