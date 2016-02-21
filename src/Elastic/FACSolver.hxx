#pragma once

/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

#include "FACOps.hxx"
#include "Boundary_Conditions.hxx"

namespace Elastic
{
  class FACSolver
  {
  public:
    FACSolver(const SAMRAI::tbox::Dimension& dim,
              const std::string& object_name,
              boost::shared_ptr<SAMRAI::tbox::Database> database,
              Boundary_Conditions &bc);

    ~FACSolver(void);

    void enableLogging(bool logging)
    {
      d_enable_logging = logging;
      d_fac_ops->enableLogging(d_enable_logging);
    }
    bool solveSystem(const int v, const int v_rhs)
    {
      createVectorWrappers(v, v_rhs);
      return d_fac_precond.solveSystem(*d_uv, *d_fv);
    }
    void setCoarsestLevelSolverTolerance(double tol)
    {
      d_fac_ops->setCoarsestLevelSolverTolerance(tol);
    }
    void setCoarsestLevelSolverMaxIterations(int max_iterations)
    {
      d_fac_ops->setCoarsestLevelSolverMaxIterations(max_iterations);
    }

    /// These functions pass std::string's to
    /// xfer::Geometry::lookupRefineOperator() to get refinement and
    /// coarsening operators
    void setCoarseFineDiscretization(const std::string& coarsefine_method)
    {
      d_fac_ops->setCoarseFineDiscretization(coarsefine_method);
    }
    void set_V_ProlongationMethod(const std::string& v_prolongation_method)
    {
      d_fac_ops->set_V_ProlongationMethod(v_prolongation_method);
    }
    void initializeSolverState
    (const int cell_moduli_id,
     const int edge_moduli_id,
     const int dv_diagonal_id,
     const int dv_mixed_id, 
     const int level_set_id,
     const int v_id,
     const int v_rhs_id,
     boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
     const int coarse_level = -1,
     const int fine_level = -1);

    void deallocateSolverState();
    int getNumberOfIterations() const
    {
      return d_fac_precond.getNumberOfIterations();
    }
    void getConvergenceFactors(double& avg_factor, double& final_factor) const
    {
      d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
    }
    double getResidualNorm() const
    {
      return d_fac_precond.getResidualNorm();
    }

    void set_boundaries(const int &v_id,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
                        const bool &homogeneous)
    {
      d_fac_ops->set_boundaries(v_id,level,homogeneous);
    }
  private:
    const SAMRAI::tbox::Dimension d_dim;
    std::string d_object_name;
    Boundary_Conditions &d_boundary_conditions;
    /// The FACOps has to be a shared_ptr so that we can pass it to
    /// FACPreconditioner's constructor.
    boost::shared_ptr<FACOps> d_fac_ops;
    SAMRAI::solv::FACPreconditioner d_fac_precond;
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;
    int d_ln_min;
    int d_ln_max;
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_uv;
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_fv;
    bool d_solver_is_initialized;
    bool d_enable_logging;

    static bool s_initialized;
    static int s_weight_id[SAMRAI::MAX_DIM_VAL];
    static int s_instance_counter[SAMRAI::MAX_DIM_VAL];

    void getFromInput(SAMRAI::tbox::Database &database);
    void createVectorWrappers(int v, int v_rhs);
    void destroyVectorWrappers()
    {
      d_uv.reset();
      d_fv.reset();
    }
    static void initializeStatics();
  };
}


