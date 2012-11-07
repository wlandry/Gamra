/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Operator class for cell-centered scalar Elastic using FAC 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_FACOPS_H
#define GAMRA_ELASTIC_FACOPS_H

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/FACOperatorStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Timer.h"
#include "Elastic/V_Refine_Patch_Strategy.h"
#include "Elastic/V_Coarsen_Patch_Strategy.h"

#include <string>

namespace Elastic {

  /*!
   * @brief FAC operator class to solve Elastic's equation on a SAMR grid,
   * using cell-centered, second-order finite-volume method, with Robin
   * boundary conditions.
   *
   * This class provides operators that are used by the FAC
   * preconditioner FACPreconditioner.
   * It is used to solve the scalar Elastic's equation using a cell-centered
   * second-order finite-volume discretization.
   * It is designed to provide all operations specific to
   * the scalar Elastic's equation,
   * @f[ \nabla \cdot D \nabla u + C u = f @f]
   * (see Elastic::Specifications) where
   * - C, D and f are indpendent of u
   * - C is a cell-centered scalar field
   * - D is the @em diffusion @em coefficients, stored on faces
   * - f is a cell-centered scalar function
   *
   * You are left to provide the source function, initial guess, etc.,
   * by specifying them in specific forms.
   *
   * This class provides:
   * -# 5-point (second order), cell-centered stencil operations
   *    for the discrete Laplacian.
   * -# Red-black Gauss-Seidel smoothing.
   * -# Provisions for working Robin boundary conditions
   *    (see RobinBcCoefStrategy).
   *
   * This class is meant to provide the Elastic-specific operator
   * used by the FAC preconditioner, FACPreconditioner.
   * To use the preconditioner with this class, you will have to provide:
   * -# The solution vector SAMRAIVectorReal,
   *    with appropriate norm weighting for the cell-centered AMR mesh.
   *    This class provides the function computeVectorWeights()
   *    to help with computing the appropriate weights.
   *    Since this is for a scalar equation, only the first depth
   *    of the first component of the vectors are used.
   *    All other parts are ignored.
   * -# The source vector SAMRAIVectorReal for f.
   * -# A Elastic::Specifications objects to specify
   *    the cell-centered scalar field C and the side-centered
   *    diffusion coefficients D
   * -# The boundary condition specifications in terms of the coefficients
   *    @f$ \alpha @f$, @f$ \beta @f$ and @f$ \gamma @f$ in the
   *    Robin formula @f$  \alpha u + \beta u_n = \gamma @f$ applied on the
   *    boundary faces.  See RobinBcCoefStrategy.
   *
   * This class allocates and deallocates only its own scratch data.
   * Other data that it manipuates are passed in as function arguments.
   * Hence, it owns none of the solution vectors, error vectors,
   * diffusion coefficient data, or any such things.
   *
   * Input Examples
   * @verbatim
   * coarse_solver_choice = "hypre"    // see setCoarsestLevelSolverChoice()
   * coarse_solver_tolerance = 1e-14   // see setCoarsestLevelSolverTolerance()
   * coarse_solver_max_iterations = 10 // see setCoarsestLevelSolverMaxIterations()
   * smoothing_choice = "Tackley"     // see setSmoothingChoice()
   * cf_discretization = "Ewing"       // see setCoarseFineDiscretization()
   * prolongation_method = "P_REFINE" // see setProlongationMethod()
   * hypre_solver = { ... }            // tbox::Database for initializing Hypre solver
   * @endverbatim
   */
  class FACOps:
    public SAMRAI::solv::FACOperatorStrategy
  {

  public:
    /*!
     * @brief Constructor.
     *
     * If you want standard output and logging,
     * pass in valid pointers for those streams.
     * @param object_name Ojbect name
     * @param database Input database
     */
    FACOps(const SAMRAI::tbox::Dimension& dim,
           const std::string& object_name,
           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database,
           Boundary_Conditions &bc);

    /*!
     * @brief Destructor.
     *
     * Deallocate internal data.
     */
    ~FACOps(void) {}

    /*!
     * @brief Enable logging.
     *
     * By default, logging is disabled.  The logging flag is
     * propagated to the major components used by this class.
     */
    void
    enableLogging(bool enable_logging);

    //@{
    /*!
     * @name Functions for setting solver mathematic algorithm controls
     */

    /*!
     * @brief Set the choice of smoothing algorithms.
     *
     * Current smoothing choices are:
     * - "Tackley"
     * - "Gerya"
     */
    void
    setSmoothingChoice(const std::string& smoothing_choice);

    /*!
     * @brief Set coarse level solver.
     *
     * Select from these:
     * - @c "Tackley" (red-black smoothing until convergence--very slow!)
     * - @c "Gerya" (red-black smoothing until convergence--very slow!)
     * - @c "hypre" (only if the HYPRE library is available).
     */
    void
    setCoarsestLevelSolverChoice(const std::string& choice);

    /*!
     * @brief Set tolerance for coarse level solve.
     *
     * If the coarse level solver requires a tolerance (currently, they all do),
     * the specified value is used.
     */
    void
    setCoarsestLevelSolverTolerance(double tol);

    /*!
     * @brief Set max iterations for coarse level solve.
     *
     * If the coarse level solver requires a max iteration limit
     * (currently, they all do), the specified value is used.
     */
    void
    setCoarsestLevelSolverMaxIterations(int max_iterations);

    /*!
     * @brief Set the coarse-fine boundary discretization method.
     *
     * Specify the @c op_name std::string which will be passed to
     * xfer::Geometry::lookupRefineOperator() to get the operator
     * for setting fine grid ghost cells from the coarse grid.
     * Note that chosing this operator implicitly choses the
     * discretization method at the coarse-fine boundary.
     *
     * There is one important instance where this std::string is
     * @em not passed to xfer::Geometry::lookupRefineOperator.
     * If this variable is set to "Ewing", Ewing's coarse-fine
     * discretization is used (a constant refinement is performed,
     * and the flux is later corrected to result in Ewing's scheme).
     * For a reference to Ewing's discretization method, see
     * "Local Refinement Techniques for Elliptic Problems on Cell-Centered
     * Grids, I. Error Analysis", Mathematics of Computation, Vol. 56, No. 194,
     * April 1991, pp. 437-461.
     *
     * @param coarsefine_method String selecting the coarse-fine discretization method.
     */
    void
    setCoarseFineDiscretization(const std::string& coarsefine_method);

    /*!
     * @brief Set the name of the prolongation method.
     *
     * Specify the @c op_name std::string which will be passed to
     * xfer::Geometry::lookupRefineOperator() to get the operator
     * for prolonging the coarse-grid correction.
     *
     * By default, "CONSTANT_REFINE" is used.  "LINEAR_REFINE" seems to
     * to lead to faster convergence, but it does NOT satisfy the Galerkin
     * condition.
     *
     * Prolonging using linear refinement requires a Robin bc
     * coefficient implementation that is capable of delivering
     * coefficients for non-hierarchy data, because linear refinement
     * requires boundary conditions to be set on temporary levels.
     *
     * @param prolongation_method String selecting the coarse-fine
     *        discretization method.
     */
    void
    set_V_ProlongationMethod(const std::string& prolongation_method);

    //@}

    //@{
    /*!
     * @name Functions for setting patch data indices and coefficients
     */

    /*!
     * @brief Set the scratch patch data index for the flux.
     *
     * The use of this function is optional.
     * The patch data index should be a pdat::SideData<DIM> type of variable.
     * If the flux id is -1 (the default initial value), scratch space
     * for the flux is allocated as needed and immediately deallocated
     * afterward, level by level.  If you have space preallocated for
     * flux and you would like that to be used, set flux id to the
     * patch data index of that space.
     */
    void set_moduli_id(const int &cell_moduli, const int &edge_moduli)
    {
      cell_moduli_id=cell_moduli;
      edge_moduli_id=edge_moduli;
    }
    //@}

    /*!
     * @brief Provide an implementation for getting the
     * physical bc coefficients
     *
     * If your solution is fixed at the physical boundary
     * ghost cell centers AND those cells have the correct
     * values before entering solveSystem(), you may use a
     * GhostCellRobinBcCoefs object.
     *
     * If your solution is @b not fixed at the ghost cell centers,
     * the ghost cell values will change as the interior
     * cell values change.  In those cases, the flexible
     * Robin boundary conditions are applied.  You must
     * call this function to provide the implementation for
     * determining the boundary condition coefficients.
     *
     * @param physical_bc_coef tbox::Pointer to an object that can
     *        set the Robin bc coefficients.
     */
    // void
    // setPhysicalBcCoefObject(
    //    const RobinBcCoefStrategy* physical_bc_coef);

    /*!
     * @brief Set weight appropriate for computing vector norms.
     *
     * If you this function to set the weights used when you
     * SAMRAIVectorReal::addComponent, you can use the
     * vector norm functions of SAMRAIVectorReal, and
     * the weights will be used to blank out coarse grid
     * regions under fine grids.
     *
     * The weights computed are specific to the cell-centered
     * discretization used by this class.  The weight is equal
     * to the cell volume if the cell has not been refined,
     * and zero if it has.
     *
     * This function is state-independent.  All inputs are in
     * the argument list.
     *
     * @param hierarchy Hierarchy configuration to compute weights for
     * @param weight_id hier::Patch data index of the weight
     * @param coarsest_ln Coarsest level number.  Must be included
     *        in hierarchy.  Must not be greater than @c finest_ln.
     *        Default to 0.
     * @param finest_ln Finest level number.  Must be included
     *        in hierarchy.  Must not be less than @c coarsest_ln.
     *        Default to finest level in @c hierarchy.
     */
    void
    computeVectorWeights(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy> hierarchy,
                         int weight_id,
                         int coarsest_ln = -1,
                         int finest_ln = -1) const;

    /*!
     * @brief Set the FAC preconditioner that will be using this object.
     *
     * The FAC preconditioner is accessed to get convergence data during
     * the cycle postprocessing step.  It is optional.
     */
    void
    setPreconditioner(const SAMRAI::solv::FACPreconditioner* preconditioner);

    //@{ @name FACOperatorStrategy virtuals

    virtual void
    restrictSolution(const SAMRAI::solv::SAMRAIVectorReal<double>& source,
                     SAMRAI::solv::SAMRAIVectorReal<double>& dest,
                     int dest_ln);
    virtual void
    restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<double>& source,
                     SAMRAI::solv::SAMRAIVectorReal<double>& dest,
                     int dest_ln);

    virtual void
    prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<double>& source,
                           SAMRAI::solv::SAMRAIVectorReal<double>& dest,
                           int dest_ln);

    virtual void
    smoothError(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                int ln,
                int num_sweeps);

    virtual int
    solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<double>& error,
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
      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> level =
        d_hierarchy->getPatchLevel(l);
      set_boundaries(v_id,level,rhs);
    }
    void set_boundaries(const int &v_id,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> &level)
    {
      set_boundaries(v_id,level,true);
    }
    void set_boundaries(const int &v_id, 
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> &level,
                        const bool &rhs);

    //@}

  private:
    //@{
    /*!
     * @name Private workhorse functions.
     */

    /*!
     * @brief Red-black Gauss-Seidel error smoothing on a level.
     *
     * Smoothes on the residual equation @f$ Ae=r @f$ on a level.
     *
     * @param error error vector
     * @param residual residual vector
     * @param ln level number
     * @param num_sweeps number of sweeps
     * @param residual_tolerance the maximum residual considered to be
     *        converged
     */
    void smooth_Tackley_2D(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                           const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                           int ln,
                           int num_sweeps,
                           double residual_tolerance = -1.0);

    void smooth_Tackley_3D
    (SAMRAI::solv::SAMRAIVectorReal<double>& solution,
     const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
     int ln,
     int num_sweeps,
     double residual_tolerance = -1.0);

    void smooth_V_2D
    (const int &axis,
     const SAMRAI::hier::Box &pbox,
     const SAMRAI::pdat::CellIndex &center,
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

    void smooth_V_3D
    (const int &ix,
     const SAMRAI::hier::Box &pbox,
     SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_rhs,
     SAMRAI::pdat::CellData<double> &cell_moduli,
     SAMRAI::pdat::EdgeData<double> &edge_moduli,
     const SAMRAI::pdat::CellIndex &center,
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
       perspective of vx, but pass in different values for center_x
       etc. to get vy. */

    double aligned_terms(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::CellIndex &center,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::hier::Index &ip,
                         const double &dx)
    {
      return (( v(x+ip)-v(x   ))*(cell_moduli(center   ,0)+2*cell_moduli(center   ,1))
              -(v(x   )-v(x-ip))*(cell_moduli(center-ip,0)+2*cell_moduli(center-ip,1)))/(dx*dx);
    }

    double lame_mixed(const SAMRAI::pdat::SideData<double> &v,
                      const SAMRAI::pdat::CellData<double> &cell_moduli,
                      const SAMRAI::pdat::CellIndex &center,
                      const SAMRAI::pdat::SideIndex &y,
                      const SAMRAI::hier::Index &ip,
                      const SAMRAI::hier::Index &jp,
                      const double &dx,
                      const double &dy)
    {
      return (+ cell_moduli(center   ,0)*(v(y+jp   )-v(y   ))/dy
              - cell_moduli(center-ip,0)*(v(y+jp-ip)-v(y-ip))/dy)/dx;
    }

    double v_operator_2D(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::NodeData<double> &edge_moduli,
                         const SAMRAI::pdat::CellIndex &center,
                         const SAMRAI::pdat::NodeIndex &edge,
                         const SAMRAI::pdat::SideIndex &x,
                         const SAMRAI::pdat::SideIndex &y,
                         const SAMRAI::hier::Index &ip,
                         const SAMRAI::hier::Index &jp,
                         const double &dx,
                         const double &dy)
    {
      return aligned_terms(v,cell_moduli,center,x,ip,dx)
        +lame_mixed(v,cell_moduli,center,y,ip,jp,dx,dy)
        +shear_noncell(v,edge_moduli,x,y,edge,ip,jp,dx,dy);
    }

    double v_operator_3D(const SAMRAI::pdat::SideData<double> &v,
                         const SAMRAI::pdat::CellData<double> &cell_moduli,
                         const SAMRAI::pdat::EdgeData<double> &edge_moduli,
                         const SAMRAI::pdat::CellIndex &center,
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
      return aligned_terms(v,cell_moduli,center,x,ip,dx)
        +lame_mixed(v,cell_moduli,center,y,ip,jp,dx,dy)
        +lame_mixed(v,cell_moduli,center,z,ip,kp,dx,dz)
        +shear_noncell(v,edge_moduli,x,y,edge_y,ip,jp,dx,dy)
        +shear_noncell(v,edge_moduli,x,z,edge_z,ip,kp,dx,dz);
    }

    /*!
     * @brief Solve the coarsest level using HYPRE
     */
    int
    solveCoarsestLevel_HYPRE(SAMRAI::solv::SAMRAIVectorReal<double>& error,
                             const SAMRAI::solv::SAMRAIVectorReal<double>& residual,
                             int ln);

    /*!
     * @brief AMR-unaware function to red or black smoothing on a single patch,
     * for variable diffusion coefficient and variable scalar field.
     *
     * @param patch patch
     * @param flux_data side-centered flux data
     * @param rhs_data cell-centered rhs data
     * @param scalar_field_data
     *        cell-centered scalar field data
     * @param soln_data cell-centered solution data
     * @param red_or_black red-black switch.  Set to 'r' or 'b'.
     * @param p_maxres max residual output.  Set to NULL to avoid computing.
     */
    void
    redOrBlackSmoothingOnPatch(const SAMRAI::hier::Patch& patch,
                               const SAMRAI::pdat::SideData<double>& flux_data,
                               const SAMRAI::pdat::CellData<double>& rhs_data,
                               SAMRAI::pdat::CellData<double>& soln_data,
                               char red_or_black,
                               double* p_maxres = NULL) const;

    //@}

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
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy> d_hierarchy;

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
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary> >
    d_cf_boundary;

    //@}

    //@{
    /*!
     * @name Private state variables for solution process.
     */

    /*!
     * @brief Smoothing choice.
     * @see setSmoothingChoice.
     */
    std::string d_smoothing_choice;

    /*!
     * @brief Coarse level solver.
     * @see setCoarsestLevelSolverChoice
     */
    std::string d_coarse_solver_choice;

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
     * @brief Id of moduli and dp.
     *
     * @see set_moduli_id.
     */
    int cell_moduli_id, edge_moduli_id;

    /*!
     * @brief Externally provided physical boundary condition object.
     *
     * see setPhysicalBcCoefObject()
     */
    // const RobinBcCoefStrategy* d_physical_bc_coef;

    //@}

    //@{ @name Internal context and scratch data

    static SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<double> >
    s_cell_scratch_var[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    static SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<double> >
    s_side_scratch_var[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

    /*!
     * @brief Default context of internally maintained hierarchy data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

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
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
    v_prolongation_refine_operator;
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule> >
    v_prolongation_refine_schedules;

    //! @brief Solution restriction (coarsening) operator.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
    v_urestriction_coarsen_operator;
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule> >
    v_urestriction_coarsen_schedules;

    //! @brief Residual restriction (coarsening) operator.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator>
    v_rrestriction_coarsen_operator;
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule> >
    v_rrestriction_coarsen_schedules;

    //! @brief Refine operator for data from coarser level.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator>
    v_ghostfill_refine_operator;
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule> >
    v_ghostfill_refine_schedules;

    //! @brief Refine operator for data from same level.
    SAMRAI::tbox::Array<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule> >
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

    Boundary_Conditions &d_boundary_conditions;
        
    /*!
     * @brief Hierarchy side operator used in debugging.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<double> >
    d_hopsside;

    /*!
     * @brief Timers for performance measurement.
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_restrict_solution;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_restrict_residual;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_prolong;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_smooth_error;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_solve_coarsest;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_composite_residual;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_residual_norm;

    static SAMRAI::tbox::StartupShutdownManager::Handler
    s_finalize_handler;
  };

}

#ifdef SAMRAI_INLINE
#include "Elastic/FACOps.I"
#endif

#endif // included_solv_ElasticFACOps
