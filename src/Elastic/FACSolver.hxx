#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

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

    bool solveSystem(const int v, const int v_rhs)
    {
      createVectorWrappers(v, v_rhs);
      return preconditioner.solveSystem(*uv, *fv);
    }
    void setCoarsestLevelSolverTolerance(double tol)
    {
      operators->setCoarsestLevelSolverTolerance(tol);
    }
    void setCoarsestLevelSolverMaxIterations(int max_iterations)
    {
      operators->setCoarsestLevelSolverMaxIterations(max_iterations);
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
      return preconditioner.getNumberOfIterations();
    }
    void getConvergenceFactors(double& avg_factor, double& final_factor) const
    {
      preconditioner.getConvergenceFactors(avg_factor, final_factor);
    }
    double getResidualNorm() const
    {
      return preconditioner.getResidualNorm();
    }

    void set_physical_boundaries
    (const int &v_id,
     boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
     const bool &homogeneous)
    {
      operators->set_physical_boundaries(v_id,level,homogeneous);
    }
  private:
    const SAMRAI::tbox::Dimension dimension;
    Boundary_Conditions &boundary_conditions;
    /// The FACOps has to be a shared_ptr so that we can pass it to
    /// FACPreconditioner's constructor.
    boost::shared_ptr<FACOps> operators;
    SAMRAI::solv::FACPreconditioner preconditioner;
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy;
    int level_min;
    int level_max;
    boost::shared_ptr<SAMRAI::hier::VariableContext> context;
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > uv;
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > fv;

    void getFromInput(SAMRAI::tbox::Database &database);
    void createVectorWrappers(int v, int v_rhs);
    void destroyVectorWrappers()
    {
      uv.reset();
      fv.reset();
    }
    static void initializeStatics();
  };
}


