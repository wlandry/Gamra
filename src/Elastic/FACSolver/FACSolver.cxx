/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

Elastic::FACSolver::FACSolver
(const SAMRAI::tbox::Dimension& dim,
 const std::string& object_name,
 boost::shared_ptr<SAMRAI::tbox::Database> database,
 Boundary_Conditions &bc):
  d_dim(dim),
  d_object_name(object_name),
  d_boundary_conditions(bc),
  d_fac_ops(boost::make_shared<FACOps>(d_dim,database,bc)),
  d_fac_precond(object_name + "::fac_precond",d_fac_ops,database),
  d_ln_min(-1),
  d_ln_max(-1),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::CONTEXT")),
  d_solver_is_initialized(false)
{
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);

  d_fac_ops->setPreconditioner(&d_fac_precond);
  if (database)
    { getFromInput(*database); }
}
