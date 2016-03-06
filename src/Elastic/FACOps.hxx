#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/solv/FACOperatorStrategy.h>
#include <SAMRAI/solv/FACPreconditioner.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>

#include "Elastic/V_Refine_Patch_Strategy.hxx"
#include "Elastic/V_Coarsen_Patch_Strategy.hxx"
#include "Elastic/Coarse_Fine_Boundary_Refine.hxx"
#include "Elastic/V_Refine.hxx"

namespace Elastic
{
  class FACOps: public SAMRAI::solv::FACOperatorStrategy
  {
  public:
    FACOps(const SAMRAI::tbox::Dimension& dim,
           const boost::shared_ptr<SAMRAI::tbox::Database> &database,
           const Boundary_Conditions &bc);
    void setCoarsestLevelSolverTolerance(double tol)
    {
      d_coarse_solver_tolerance = tol;
    }
    void setCoarsestLevelSolverMaxIterations(int max_iterations)
    {
      if (max_iterations < 0)
        { TBOX_ERROR(__FILE__ << "Invalid number of max iterations"); }
      d_coarse_solver_max_iterations = max_iterations;
    }

    void set_extra_ids(const int &Cell_moduli_id, const int &Edge_moduli_id,
                       const int &Dv_diagonal_id, const int &Dv_mixed_id,
                       const int &Level_set_id)
    {
      cell_moduli_id=Cell_moduli_id;
      edge_moduli_id=Edge_moduli_id;
      dv_diagonal_id=Dv_diagonal_id;
      dv_mixed_id=Dv_mixed_id;
      level_set_id=Level_set_id;

      Elastic::Coarse_Fine_Boundary_Refine::dv_diagonal_id=dv_diagonal_id;
      Elastic::Coarse_Fine_Boundary_Refine::dv_mixed_id=dv_mixed_id;
      Elastic::Coarse_Fine_Boundary_Refine::level_set_id=level_set_id;
      Elastic::V_Refine::level_set_id=level_set_id;
      v_coarsen_patch_strategy.set_extra_ids(dv_diagonal_id,dv_mixed_id,
                                             level_set_id);
    }

    bool have_faults() const
    {
      return dv_diagonal_id!=invalid_id;
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }

    void setPreconditioner(const SAMRAI::solv::FACPreconditioner* preconditioner)
    {
      d_preconditioner = preconditioner;
    }
    virtual void restrictSolution
    (const SAMRAI::solv::SAMRAIVectorReal<double>& source,
     SAMRAI::solv::SAMRAIVectorReal<double>& dest,
     int dest_ln);
    virtual void restrictResidual
    (const SAMRAI::solv::SAMRAIVectorReal<double>& source,
     SAMRAI::solv::SAMRAIVectorReal<double>& dest,
     int dest_ln);

    virtual void prolongErrorAndCorrect
    (const SAMRAI::solv::SAMRAIVectorReal<double>& source,
     SAMRAI::solv::SAMRAIVectorReal<double>& dest,
     int dest_ln);

    virtual void smoothError
    (SAMRAI::solv::SAMRAIVectorReal<double>& error,
     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     int ln,
     int num_sweeps)
    {
      smooth(error,residual,ln,num_sweeps,0);
    }
    
    void smooth(SAMRAI::solv::SAMRAIVectorReal<double>& data,
                const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                const int &ln,
                const int &num_sweeps,
                const double &tolerance);

    virtual int solveCoarsestLevel
    (SAMRAI::solv::SAMRAIVectorReal<double>& error,
     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     int coarsest_ln);

    virtual void
    computeCompositeResidualOnLevel
    (SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
     const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
     int ln,
     bool error_equation_indicator);

    virtual double
    computeResidualNorm(const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                        int fine_ln,
                        int coarse_ln);

    virtual void initializeOperatorState
    (const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
     const SAMRAI::solv::SAMRAIVectorReal<double>& rhs);

    virtual void deallocateOperatorState();

    virtual void
    postprocessOneCycle(int fac_cycle_num,
                        const SAMRAI::solv::SAMRAIVectorReal<double> &,
                        const SAMRAI::solv::SAMRAIVectorReal<double> &);

    void set_physical_boundaries(const int &v_id,
                                 const SAMRAI::hier::PatchHierarchy &hierarchy,
                                 const int &l)
    {
      set_physical_boundaries(v_id,hierarchy,l,true);
    }
    void set_physical_boundaries(const int &v_id,
                                 const SAMRAI::hier::PatchHierarchy &hierarchy,
                                 const int &l, const bool &rhs)
    {
      set_physical_boundaries(v_id,hierarchy.getPatchLevel(l),rhs);
    }
    void set_physical_boundaries
    (const int &v_id,
     const boost::shared_ptr<SAMRAI::hier::PatchLevel> &level)
    {
      set_physical_boundaries(v_id,level,true);
    }
    void set_physical_boundaries
    (const int &v_id, 
     const boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
     const bool &rhs);
  private:
    void Gauss_Seidel_red_black_2D
    (SAMRAI::solv::SAMRAIVectorReal<double>& error,
     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     int ln, int num_sweeps, double residual_tolerance);

    void Gauss_Seidel_red_black_3D
    (SAMRAI::solv::SAMRAIVectorReal<double>& solution,
     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     int ln, int num_sweeps,
     double residual_tolerance);

    void update_V_2D(const Gamra::Dir &axis,
                     const SAMRAI::hier::Box &pbox,
                     const SAMRAI::pdat::CellIndex &cell,
                     const SAMRAI::hier::Index &ip,
                     const SAMRAI::hier::Index &jp,
                     SAMRAI::pdat::SideData<double> &v,
                     SAMRAI::pdat::SideData<double> &v_rhs,
                     double &maxres,
                     const double &dx,
                     const double &dy,
                     SAMRAI::pdat::CellData<double> &cell_moduli,
                     SAMRAI::pdat::NodeData<double> &edge_moduli,
                     const double &theta_momentum);

    void update_V_3D(const Gamra::Dir &ix,
                     const SAMRAI::hier::Box &pbox,
                     SAMRAI::pdat::SideData<double> &v,
                     SAMRAI::pdat::SideData<double> &v_rhs,
                     SAMRAI::pdat::CellData<double> &cell_moduli,
                     SAMRAI::pdat::EdgeData<double> &edge_moduli,
                     const SAMRAI::pdat::CellIndex &cell,
                     const double Dx[3],
                     const double &theta_momentum,
                     const SAMRAI::hier::Index pp[3],
                     double &maxres);

    void refine(int v_dst, int v_src, int v_scr, int dest_ln);
    void coarsen_solution(int v_dst, int v_src, int dest_ln);
    void coarsen_resid(int v_dst, int v_src, int dest_ln);
    void ghostfill(int v_id, int dest_ln);
    void ghostfill_nocoarse(int v_id, int dest_ln);

    static void finalizeCallback();

    const SAMRAI::tbox::Dimension d_dim;
    int d_ln_min;
    int d_ln_max;

    std::vector<boost::shared_ptr<SAMRAI::hier::CoarseFineBoundary> >
    d_cf_boundary;
    double d_coarse_solver_tolerance;
    int d_coarse_solver_max_iterations;
    int cell_moduli_id, edge_moduli_id, dv_diagonal_id, dv_mixed_id,
      level_set_id;

    static boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;

    int d_side_scratch_id;
    boost::shared_ptr<SAMRAI::hier::RefineOperator> refine_operator;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    refine_schedules;

    boost::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsen_solution_operator;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    coarsen_solution_schedules;

    boost::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsen_resid_operator;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    coarsen_resid_schedules;

    boost::shared_ptr<SAMRAI::hier::RefineOperator> ghostfill_operator;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    ghostfill_schedules;

    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    ghostfill_nocoarse_schedules;

    V_Refine_Patch_Strategy v_refine_patch_strategy;
    V_Coarsen_Patch_Strategy v_coarsen_patch_strategy;

  public:
    bool logging;
  private:
    const SAMRAI::solv::FACPreconditioner* d_preconditioner;
    const Boundary_Conditions &d_boundary_conditions;
        
    boost::shared_ptr<SAMRAI::tbox::Timer> t_restrict_solution;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_restrict_residual;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_prolong;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_smooth_error;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_solve_coarsest;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_compute_composite_residual;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_compute_residual_norm;

    static SAMRAI::tbox::StartupShutdownManager::Handler s_finalize_handler;
  };
}
