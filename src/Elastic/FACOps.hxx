#pragma once

/// Copyright: (c) 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright: (c) 2013-2016 California Institute of Technology

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
    void enableLogging(bool enable_logging);
    void setCoarsestLevelSolverTolerance(double tol);
    void setCoarsestLevelSolverMaxIterations(int max_iterations);
    void setCoarseFineDiscretization(const std::string& coarsefine_method);
    void set_V_ProlongationMethod(const std::string& prolongation_method);

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

    void setPreconditioner(const SAMRAI::solv::FACPreconditioner*
                           preconditioner);

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

    virtual void
    initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                            const SAMRAI::solv::SAMRAIVectorReal<double>& rhs);

    virtual void
    deallocateOperatorState();

    virtual void
    postprocessOneCycle(int fac_cycle_num,
                        const SAMRAI::solv::SAMRAIVectorReal<double>& current_soln,
                        const SAMRAI::solv::SAMRAIVectorReal<double>& residual);

    void set_boundaries(const int &v_id, const int &l)
    {
      set_boundaries(v_id,l,true);
    }
    void set_boundaries(const int &v_id, const int &l, const bool &rhs)
    {
      boost::shared_ptr<SAMRAI::hier::PatchLevel> level =
        d_hierarchy->getPatchLevel(l);
      set_boundaries(v_id,level,rhs);
    }
    void set_boundaries(const int &v_id,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> &level)
    {
      set_boundaries(v_id,level,true);
    }
    void set_boundaries(const int &v_id, 
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
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

    //@{ @name For executing, caching and resetting communication schedules.

    /*!
     * @brief Execute a refinement schedule
     * for prolonging cell data.
     *
     * General notes regarding internal objects for communication:
     * We maintain objects to support caching schedules to improve
     * efficiency.  Communication is needed in 5 distinct tasks.
     *   -# Prolongation
     *   -# Restriction
     *   -# Flux coarsening.  Changing the coarse grid flux to the
     *      composite grid flux by coarsening the fine grid flux
     *      at the coarse-fine boundaries.
     *   -# Fill boundary data from other patches in the same level
     *      and physical boundary condition.
     *   -# Fill boundary data from same level, coarser levels
     *      and physical boundary condition.
     *
     * For each task, we maintain a refine or coarsen operator,
     * and a array of communication schedules (one for each
     * destination level).
     *
     * The 5 member functions named @c xeqSchedule... execute
     * communication schedules appropriate for five specific tasks.
     * They use a cached schedule if possible or create and cache
     * a new schedule if needed.  These functions and the data
     * they manipulate are as follows:
     * <ol>
     *   <li> xeqScheduleProlongation():
     *        prolongation_refine_operator
     *        prolongation_refine_schedules
     *   <li> xeqScheduleURestriction():
     *        d_restriction_coarsen_operator,
     *        urestriction_coarsen_schedules.
     *   <li> xeqScheduleRRestriction():
     *        restriction_coarsen_operator,
     *        rrestriction_coarsen_schedules.
     *   <li> xeqScheduleFluxCoarsen():
     *        d_flux_coarsen_operator,
     *        d_flux_coarsen_schedules.
     *   <li> xeqScheduleGhostFill():
     *        ghostfill_refine_operator,
     *        ghostfill_refine_schedules.
     *   <li> xeqScheduleGhostFillNoCoarse():
     *        ghostfill_nocoarse_refine_operator,
     *        ghostfill_nocoarse_refine_schedules.
     * </ol>
     *
     * @return refinement schedule for prolongation
     */
    void
    xeqScheduleProlongation(int v_dst, int v_src, int v_scr,
                            int dest_ln);

    /*!
     * @brief Execute schedule for restricting solution to the specified
     * level or reregister an existing one.
     *
     * See general notes for xeqScheduleProlongation().
     *
     * @return coarsening schedule for restriction
     */
    void
    xeqScheduleURestriction(int v_dst, int v_src, int dest_ln);
                                

    /*!
     * @brief Execute schedule for restricting residual to the specified
     * level or reregister an existing one.
     *
     * See general notes for xeqScheduleProlongation().
     *
     * @return coarsening schedule for restriction
     */
    void
    xeqScheduleRRestriction(int v_dst, int v_src, int dest_ln);

    /*!
     * @brief Execute schedule for coarsening flux to the specified
     * level or reregister an existing one.
     *
     * See general notes for xeqScheduleProlongation().
     *
     * @return coarsening schedule for setting composite grid flux at
     * coarse-fine boundaries.
     */
    void
    xeqScheduleFluxCoarsen(
                           int dst_id,
                           int src_id,
                           int dest_ln);

    /*!
     * @brief Execute schedule for filling ghosts on the specified
     * level or reregister an existing one.
     *
     * See general notes for xeqScheduleProlongation().
     *
     * @return refine schedule for filling ghost data from coarser level
     * and physical bc.
     */
    void
    xeqScheduleGhostFill(int v_id, int dest_ln);

    /*!
     * @brief Execute schedule for filling ghosts on the specified
     * level or reregister an existing one.
     * This version does not get data from coarser levels.
     *
     * See general notes for xeqScheduleProlongation().
     *
     * This function is used for the bottom solve level, since it does
     * not access data from any coarser level.  (Ghost data obtained
     * from coarser level must have been placed there before solve begins!)
     *
     * @return refine schedule for filling ghost data from same level
     * and physical bc.
     */
    void
    xeqScheduleGhostFillNoCoarse(int v_id, int dest_ln);

    //@}

    //! @brief Return the patch data index for cell scratch data.
    int
    registerCellScratch() const;
    //! @brief Return the patch data index for flux scratch data.
    int
    registerFluxScratch() const;
    //! @brief Return the patch data index for outerflux scratch data.
    int
    registerOfluxScratch() const;

    //! @brief Free static variables at shutdown time.
    static void
    finalizeCallback();

    /*!
     * @brief Object dimension.
     */
    const SAMRAI::tbox::Dimension d_dim;

    /*!
     * @brief Object name.
     */
    std::string d_object_name;

    //@{ @name Hierarchy-dependent objects.

    /*!
     * @brief Reference hierarchy
     *
     * This variable is non-null between the initializeOperatorState()
     * and deallocateOperatorState() calls.  It is not truly needed,
     * because the hierarchy is obtainable through variables in most
     * function argument lists.  We use it to enforce working on one
     * hierarchy at a time.
     */
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

    /*!
     * @brief Coarsest level for solve.
     */
    int d_ln_min;

    /*!
     * @brief Finest level for solve.
     */
    int d_ln_max;

    /*!
     * @brief Description of coarse-fine boundaries.
     *
     * There is one coarse-fine boundary object for each level.
     * d_coarse_fine_boundary[i] is the description of
     * the coarse-fine boundary between level i and level i-1.
     * The coarse-fine boundary does not exist at the coarsest level,
     * although the SAMRAI::hier::CoarseFineBoundary object still exists (it
     * should not contain any boxes).
     *
     * This array is initialized in initializeOperatorState() and
     * deallocated in deallocateOperatorState().  When allocated,
     * it is allocated for the index range [0,d_ln_max], though
     * the range [0,d_ln_min-1] is not used.  This is okay because
     * SAMRAI::hier::CoarseFineBoundary is a light object before
     * it is set for a level.
     */
    std::vector<boost::shared_ptr<SAMRAI::hier::CoarseFineBoundary> >
    d_cf_boundary;

    //@}

    //@{
    /*!
     * @name Private state variables for solution process.
     */

    /*!
     * @brief Coarse-fine discretization method.
     * @see setCoarseFineDiscretization().
     */
    std::string d_cf_discretization;

    /*!
     * @brief Coarse-fine discretization method.
     *
     * The name of the refinement operator used to prolong the
     * coarse grid correction.
     *
     * @see setProlongationMethod()
     */
    std::string v_prolongation_method;

    /*!
     * @brief Tolerance specified to coarse solver
     * @see setCoarsestLevelSolverTolerance()
     */
    double d_coarse_solver_tolerance;

    /*!
     * @brief Coarse level solver iteration limit.
     * @see setCoarsestLevelSolverMaxIterations()
     */
    int d_coarse_solver_max_iterations;

    /*!
     * @brief Residual tolerance to govern smoothing.
     *
     * When we use one of the internal error smoothing functions
     * and want to terminate the smoothing sweeps at a certain
     * level of residual, this will be set to > 0.  If it is
     * < 0, the smoothing function effectively ignores it.
     *
     * This variable is needed because some coarse-level solver
     * simply runs the smoothing function until convergence.
     * It sets this variable to > 0, calls the smoothing function,
     * then resets it to < 0.
     */
    double d_residual_tolerance_during_smoothing;

    /*!
     * @brief Id of extra terms.
     *
     * @see set_extra_ids.
     */
    int cell_moduli_id, edge_moduli_id, dv_diagonal_id, dv_mixed_id,
      level_set_id;

    /*!
     * @brief Externally provided physical boundary condition object.
     *
     * see setPhysicalBcCoefObject()
     */
    // const RobinBcCoefStrategy* d_physical_bc_coef;

    //@}

    //@{ @name Internal context and scratch data

    static boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
    s_cell_scratch_var[SAMRAI::MAX_DIM_VAL];

    static boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
    s_side_scratch_var[SAMRAI::MAX_DIM_VAL];

    /*!
     * @brief Default context of internally maintained hierarchy data.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;

    /*!
     * @brief ID of the solution-like scratch data.
     *
     * Set in constructor and never changed.
     * Corresponds to a SAMRAI::pdat::CellVariable<double> named
     * @c d_object_name+"::cell_scratch".
     * Scratch data is allocated and removed as needed
     * to reduce memory usage.
     */
    int d_side_scratch_id;

    //@}

    //@{
    /*!
     * @name Various refine and coarsen objects used internally.
     */

    //! @brief Error prolongation (refinement) operator.
    boost::shared_ptr<SAMRAI::hier::RefineOperator>
    v_prolongation_refine_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_prolongation_refine_schedules;

    //! @brief Solution restriction (coarsening) operator.
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    v_urestriction_coarsen_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    v_urestriction_coarsen_schedules;

    //! @brief Residual restriction (coarsening) operator.
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator>
    v_rrestriction_coarsen_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> >
    v_rrestriction_coarsen_schedules;

    //! @brief Refine operator for data from coarser level.
    boost::shared_ptr<SAMRAI::hier::RefineOperator>
    v_ghostfill_refine_operator;
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_ghostfill_refine_schedules;

    //! @brief Refine operator for data from same level.
    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> >
    v_nocoarse_refine_schedules;

    //@}

    /*!
     * @brief Utility object employed in setting ghost cells and providing
     * SAMRAI::xfer::RefinePatchStrategy implementation.
     *
     * Since this class deals only in scalar variables having
     * Robin boundary conditions, we take advantage of the corresponding
     * implementation in CartesianRobinBcHelper.  Whenever
     * we need an implementation of SAMRAI::xfer::RefinePatchStrategy,
     * this object is used.  Note that in the code, before we
     * use this object to set ghost cell values, directly or
     * indirectly by calling SAMRAI::xfer::RefineSchedule::fillData(),
     * we must tell the patch strategies the patch data index we want
     * to set and whether we are setting data with homogeneous
     * boundary condition.
     */
    V_Refine_Patch_Strategy v_refine_patch_strategy;
    V_Coarsen_Patch_Strategy v_coarsen_patch_strategy;

    //@{
    /*!
     * @name Non-essential objects used in outputs and debugging.
     */

    /*!
     * @brief Logging flag.
     */
    bool d_enable_logging;

    /*!
     * @brief Preconditioner using this object.
     *
     * See setPreconditioner().
     */
    const SAMRAI::solv::FACPreconditioner* d_preconditioner;

    const Boundary_Conditions &d_boundary_conditions;
        
    /*!
     * @brief Hierarchy side operator used in debugging.
     */
    boost::shared_ptr<SAMRAI::math::HierarchySideDataOpsReal<double> >
    d_hopsside;

    /*!
     * @brief Timers for performance measurement.
     */
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

#include "Elastic/FACOps.I"


