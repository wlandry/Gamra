/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

bool Elastic::FACSolver::s_initialized = false;

Elastic::FACSolver::FACSolver
(const SAMRAI::tbox::Dimension& dim,
 const std::string& object_name,
 boost::shared_ptr<SAMRAI::tbox::Database> database,
 Boundary_Conditions &bc):
  d_dim(dim),
  d_object_name(object_name),
  d_boundary_conditions(bc),
  d_fac_ops(boost::make_shared<FACOps>(d_dim, object_name + "::fac_ops",
                                       database,bc)),
  d_fac_precond(object_name + "::fac_precond",d_fac_ops,database),
  d_hierarchy(),
  d_ln_min(-1),
  d_ln_max(-1),
  d_context(SAMRAI::hier::VariableDatabase::getDatabase()
            ->getContext(object_name + "::CONTEXT")),
  d_uv(),
  d_fv(),
  d_solver_is_initialized(false),
  d_enable_logging(false)
{
  if (!s_initialized)
    initializeStatics();

  /// FIXME: Does this do anything?
  setCoarseFineDiscretization("Ewing");
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);

  d_fac_ops->setPreconditioner(&d_fac_precond);
  if (database)
    getFromInput(*database);

}
