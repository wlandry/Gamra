#pragma once

/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

#include <SAMRAI/SAMRAI_config.h>
#include <SAMRAI/solv/FACPreconditioner.h>
#include <SAMRAI/solv/SimpleCellRobinBcCoefs.h>
#include <SAMRAI/tbox/Database.h>

#include "Elastic/FACOps.hxx"
#include "Elastic/Boundary_Conditions.hxx"

namespace Elastic {
  class FACSolver
  {

  public:
    FACSolver(const SAMRAI::tbox::Dimension& dim,
              const std::string& object_name,
              boost::shared_ptr<SAMRAI::tbox::Database> database,
              Boundary_Conditions &bc);

    ~FACSolver(void);

    void enableLogging(bool logging);
    bool solveSystem(const int v, const int v_rhs);
    void setCoarsestLevelSolverTolerance(double tol);
    void setCoarsestLevelSolverMaxIterations(int max_iterations);

    /// These functions pass std::string's to
    /// xfer::Geometry::lookupRefineOperator() to get refinement and
    /// coarsening operators
    void setCoarseFineDiscretization(const std::string& coarsefine_method);
    void set_V_ProlongationMethod(const std::string& prolongation_method);
    void initializeSolverState
    (const int cell_moduli_id,
     const int edge_moduli_id,
     const int dv_diagonal_id,
     const int dv_mixed_id, 
     const int level_set_id,
     const int v_id,
     const int v_rhs_id,
     boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
     const int coarse_level = -1,
     const int fine_level = -1);

    void deallocateSolverState();
    int getNumberOfIterations() const;
    void getConvergenceFactors(double& avg_factor, double& final_factor) const;
    double getResidualNorm() const;

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

    void getFromInput(boost::shared_ptr<SAMRAI::tbox::Database> database);
    void createVectorWrappers(int v, int v_rhs);
    void destroyVectorWrappers();
    static void initializeStatics();
  };

}

#include "Elastic/FACSolver.I"


