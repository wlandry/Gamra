/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

#include "Elastic/FACSolver.hxx"

bool Elastic::FACSolver::s_initialized = 0;
int Elastic::FACSolver::s_weight_id[SAMRAI::MAX_DIM_VAL];
int Elastic::FACSolver::s_instance_counter[SAMRAI::MAX_DIM_VAL];

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

  setCoarseFineDiscretization("Ewing");
  setCoarsestLevelSolverTolerance(1e-8);
  setCoarsestLevelSolverMaxIterations(10);

  SAMRAI::hier::VariableDatabase*
    var_db = SAMRAI::hier::VariableDatabase::getDatabase();

  {
    static std::string cell_weight_name("Elastic::FACSolver_cell_weight");

    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > weight =
      boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
      (var_db->getVariable(cell_weight_name));
    if (!weight)
      weight = boost::make_shared<SAMRAI::pdat::CellVariable<double> >
        (d_dim, cell_weight_name, 1);

    if (s_weight_id[d_dim.getValue() - 1] < 0)
      s_weight_id[d_dim.getValue() - 1] = var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
  }

  {
    static std::string side_weight_name("Elastic::FACSolver_side_weight");

    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > weight =
      boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
      (var_db->getVariable(side_weight_name));
    if (!weight)
      weight = boost::make_shared<SAMRAI::pdat::SideVariable<double> >
        (d_dim,side_weight_name,SAMRAI::hier::IntVector::getOne(d_dim),1);

    if (s_weight_id[d_dim.getValue() - 2] < 0)
      s_weight_id[d_dim.getValue() - 2] = var_db->registerInternalSAMRAIVariable
        (weight,SAMRAI::hier::IntVector::getZero(d_dim));
  }

  d_fac_ops->setPreconditioner(&d_fac_precond);
  if (database)
    getFromInput(*database);

  s_instance_counter[d_dim.getValue() - 1]++;
}
