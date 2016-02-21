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
#include "Elastic/V_Boundary_Refine.hxx"
#include "Elastic/V_Refine.hxx"

namespace Elastic
{
  class FACOps: public SAMRAI::solv::FACOperatorStrategy
  {
  public:
    FACOps(const SAMRAI::tbox::Dimension& dim,
           const std::string& object_name,
           const boost::shared_ptr<SAMRAI::tbox::Database> &database,
           const Boundary_Conditions &bc);
    void enableLogging(bool enable_logging)
    {
      d_enable_logging = enable_logging;
    }
    void setCoarsestLevelSolverTolerance(double tol)
    {
      d_coarse_solver_tolerance = tol;
    }
    void setCoarsestLevelSolverMaxIterations(int max_iterations)
    {
      if (max_iterations < 0)
        TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
      d_coarse_solver_max_iterations = max_iterations;
    }
    void setCoarseFineDiscretization(const std::string& coarsefine_method)
    {
      if (initialized)
        TBOX_ERROR(d_object_name << ": Cannot change coarse-fine\n"
                   << "discretization method while operator state\n"
                   << "is initialized because that causes a\n"
                   << "corruption in the state.\n");
      d_cf_discretization = coarsefine_method;
    }
    void set_V_ProlongationMethod(const std::string& prolongation_method)
    {
      if (initialized)
        TBOX_ERROR(d_object_name << ": Cannot change v prolongation method\n"
                   << "while operator state is initialized because that\n"
                   << "causes a corruption in the state.\n");
      v_prolongation_method = prolongation_method;
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

      Elastic::V_Boundary_Refine::dv_diagonal_id=dv_diagonal_id;
      Elastic::V_Boundary_Refine::dv_mixed_id=dv_mixed_id;
      Elastic::V_Boundary_Refine::level_set_id=level_set_id;
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

    void computeVectorWeights(const SAMRAI::hier::PatchHierarchy &hierarchy,
                              int weight_id,
                              int coarsest_ln = -1,
                              int finest_ln = -1) const;

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
     int num_sweeps);
    
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

    void residual_2D
    (SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::CellData<double> &cell_moduli,
     SAMRAI::pdat::SideData<double> &v_rhs,
     SAMRAI::pdat::SideData<double> &v_resid,
     SAMRAI::hier::Patch &patch,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::geom::CartesianPatchGeometry &geom);

    void residual_3D
    (SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::CellData<double> &cell_moduli,
     SAMRAI::pdat::SideData<double> &v_rhs,
     SAMRAI::pdat::SideData<double> &v_resid,
     SAMRAI::hier::Patch &patch,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::geom::CartesianPatchGeometry &geom);

    virtual double
    computeResidualNorm(const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                        int fine_ln,
                        int coarse_ln);

    virtual void initializeOperatorState
    (const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
     const SAMRAI::solv::SAMRAIVectorReal<double>& rhs);

    virtual void
    deallocateOperatorState();

    virtual void
    postprocessOneCycle(int fac_cycle_num,
                        const SAMRAI::solv::SAMRAIVectorReal<double>& current_soln,
                        const SAMRAI::solv::SAMRAIVectorReal<double>& residual);

    void set_boundaries(const int &v_id,
                        const SAMRAI::hier::PatchHierarchy &hierarchy,
                        const int &l)
    {
      set_boundaries(v_id,hierarchy,l,true);
    }
    void set_boundaries(const int &v_id,
                        const SAMRAI::hier::PatchHierarchy &hierarchy,
                        const int &l, const bool &rhs)
    {
      set_boundaries(v_id,hierarchy.getPatchLevel(l),rhs);
    }
    void set_boundaries(const int &v_id,
                        const boost::shared_ptr<SAMRAI::hier::PatchLevel> &level)
    {
      set_boundaries(v_id,level,true);
    }
    void set_boundaries(const int &v_id, 
                        const boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
                        const bool &rhs);
  private:
    void smooth_2D(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                   const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                   int ln, int num_sweeps, double residual_tolerance = -1.0);

    void smooth_3D(SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                   const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                   int ln, int num_sweeps,
                   double residual_tolerance = -1.0);

    void smooth_V_2D(const Gamra::Dir &axis,
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

    void smooth_V_3D(const Gamra::Dir &ix,
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

    /* The mixed derivative of the stress.  We have to use a template
       because 2D uses Node's for the edge moduli, while 3D uses
       Edge's.  Written as if it is dtau_xy_dy. */

    template<class E_data, class E_index>
    double shear_noncell(const SAMRAI::pdat::SideData<double> &v,
                         const E_data &edge_moduli,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::pdat::SideIndex &y,
                         const E_index &edge,
                         const SAMRAI::hier::Index &ip,
                         const SAMRAI::hier::Index &jp,
                         const double &dx,
                         const double &dy)
    {
      return 
        edge_moduli(edge+jp,1)*(v(x+jp)-v(x   ))/(dy*dy)
        -edge_moduli(edge   ,1)*(v(x   )-v(x-jp))/(dy*dy)
        +edge_moduli(edge+jp,1)*(v(y+jp)-v(y+jp-ip))/(dx*dy) 
        -edge_moduli(edge   ,1)*(v(y   )-v(y-ip   ))/(dx*dy);
    }

    /* The action of the velocity operator. It is written from the
       perspective of vx, but pass in different values for cell_x
       etc. to get vy. */

    double aligned_terms(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::CellIndex &cell,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::hier::Index &ip,
                         const double &dx)
    {
      return (( v(x+ip)-v(x   ))
              *(cell_moduli(cell   ,0)+2*cell_moduli(cell   ,1))
              -(v(x   )-v(x-ip))
              *(cell_moduli(cell-ip,0)+2*cell_moduli(cell-ip,1)))/(dx*dx);
    }

    double lame_mixed(const SAMRAI::pdat::SideData<double> &v,
                      const SAMRAI::pdat::CellData<double> &cell_moduli,
                      const SAMRAI::pdat::CellIndex &cell,
                      const SAMRAI::pdat::SideIndex &y,
                      const SAMRAI::hier::Index &ip,
                      const SAMRAI::hier::Index &jp,
                      const double &dx,
                      const double &dy)
    {
      return (+ cell_moduli(cell   ,0)*(v(y+jp   )-v(y   ))/dy
              - cell_moduli(cell-ip,0)*(v(y+jp-ip)-v(y-ip))/dy)/dx;
    }

    double v_operator_2D(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::NodeData<double> &edge_moduli,
                         const SAMRAI::pdat::CellIndex &cell,
                         const SAMRAI::pdat::NodeIndex &edge,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::pdat::SideIndex &y,
                         const SAMRAI::hier::Index &ip,
                         const SAMRAI::hier::Index &jp,
                         const double &dx,
                         const double &dy)
    {
      return aligned_terms(v,cell_moduli,cell,x,ip,dx)
        +lame_mixed(v,cell_moduli,cell,y,ip,jp,dx,dy)
        +shear_noncell(v,edge_moduli,x,y,edge,ip,jp,dx,dy);
    }

    double v_level_set_operator_2D
    (const SAMRAI::pdat::SideData<double> &level_set,
     const SAMRAI::pdat::SideData<double> &v,
     const SAMRAI::pdat::CellData<double> &cell_moduli,
     const SAMRAI::pdat::NodeData<double> &edge_moduli,
     const SAMRAI::pdat::CellIndex &cell,
     const SAMRAI::pdat::NodeIndex &edge,
     const SAMRAI::pdat::SideIndex &x,
     const SAMRAI::pdat::SideIndex &y,
     const SAMRAI::hier::Index &ip,
     const SAMRAI::hier::Index &jp,
     const double &dx,
     const double &dy);

    double v_operator_3D(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::EdgeData<double> &edge_moduli,
                         const SAMRAI::pdat::CellIndex &cell,
                         const SAMRAI::pdat::EdgeIndex &edge_y,
                         const SAMRAI::pdat::EdgeIndex &edge_z,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::pdat::SideIndex &y,
                         const SAMRAI::pdat::SideIndex &z,
                         const SAMRAI::hier::Index &ip,
                         const SAMRAI::hier::Index &jp,
                         const SAMRAI::hier::Index &kp,
                         const double &dx,
                         const double &dy,
                         const double &dz)
    {
      return aligned_terms(v,cell_moduli,cell,x,ip,dx)
        +lame_mixed(v,cell_moduli,cell,y,ip,jp,dx,dy)
        +lame_mixed(v,cell_moduli,cell,z,ip,kp,dx,dz)
        +shear_noncell(v,edge_moduli,x,y,edge_z,ip,jp,dx,dy)
        +shear_noncell(v,edge_moduli,x,z,edge_y,ip,kp,dx,dz);
    }

    void xeqScheduleProlongation(int v_dst, int v_src, int v_scr, int dest_ln);
    void xeqScheduleURestriction(int v_dst, int v_src, int dest_ln);
    void xeqScheduleRRestriction(int v_dst, int v_src, int dest_ln);
    void xeqScheduleFluxCoarsen(int dst_id, int src_id, int dest_ln);
    void xeqScheduleGhostFill(int v_id, int dest_ln);
    void xeqScheduleGhostFillNoCoarse(int v_id, int dest_ln);

    static void finalizeCallback();

    const SAMRAI::tbox::Dimension d_dim;
    std::string d_object_name;

    // FIXME: We should not need this variable, because it should
    // always be initialized.
    bool initialized;
    int d_ln_min;
    int d_ln_max;

    std::vector<boost::shared_ptr<SAMRAI::hier::CoarseFineBoundary> >
    d_cf_boundary;
    std::string d_cf_discretization;
    std::string v_prolongation_method;
    double d_coarse_solver_tolerance;
    int d_coarse_solver_max_iterations;
    double d_residual_tolerance_during_smoothing;
    int cell_moduli_id, edge_moduli_id, dv_diagonal_id, dv_mixed_id,
      level_set_id;

    static boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    s_cell_scratch_var[SAMRAI::MAX_DIM_VAL];

    static boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;

    int d_side_scratch_id;
    boost::shared_ptr<SAMRAI::hier::RefineOperator>
    v_prolongation_refine_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_prolongation_refine_schedules;

    boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    v_urestriction_coarsen_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    v_urestriction_coarsen_schedules;

    boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    v_rrestriction_coarsen_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    v_rrestriction_coarsen_schedules;

    boost::shared_ptr<SAMRAI::hier::RefineOperator>
    v_ghostfill_refine_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_ghostfill_refine_schedules;

    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_nocoarse_refine_schedules;

    V_Refine_Patch_Strategy v_refine_patch_strategy;
    V_Coarsen_Patch_Strategy v_coarsen_patch_strategy;

    bool d_enable_logging;
    const SAMRAI::solv::FACPreconditioner* d_preconditioner;
    const Boundary_Conditions &d_boundary_conditions;
        
    boost::shared_ptr<SAMRAI::tbox::Timer> t_restrict_solution;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_restrict_residual;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_prolong;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_smooth_error;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_solve_coarsest;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_compute_composite_residual;
    boost::shared_ptr<SAMRAI::tbox::Timer> t_compute_residual_norm;

    static SAMRAI::tbox::StartupShutdownManager::Handler
    s_finalize_handler;
  };

}



