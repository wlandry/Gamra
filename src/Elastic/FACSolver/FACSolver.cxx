/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include "Elastic/FACSolver.hxx"

Elastic::FACSolver::FACSolver
(const SAMRAI::tbox::Dimension& dim,
 const std::string& object_name,
 boost::shared_ptr<SAMRAI::tbox::Database> database,
 Boundary_Conditions &bc):
  dimension(dim),
  boundary_conditions(bc),
  operators(boost::make_shared<FACOps>(dimension,database,bc)),
  preconditioner("Elastic::FACSolver::FACPreconditioner",operators,database),
  level_min(-1),
  level_max(-1),
  context(SAMRAI::hier::VariableDatabase::getDatabase()
          ->getContext(object_name + "::CONTEXT"))
{
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);

  operators->setPreconditioner(&preconditioner);
  if (database)
    { getFromInput(*database); }
}
