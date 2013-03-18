/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar Elastic equation. 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_FACSolver_H
#define GAMRA_ELASTIC_FACSolver_H

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/solv/FACPreconditioner.h"
#include "Elastic/FACOps.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/tbox/Database.h"
#include "Elastic/Boundary_Conditions.h"

namespace Elastic {
  /*!
   * @brief Class for solving vector Elastic's equation on SAMR grid,
   * wrapping up lower-level components (FAC cycling, Elastic equation
   * operations and boundary conditions) in a single high-level interface.
   *
   * Note: this class provides a backward-compatible interface to
   * the soon-to-be obsolete Elastic::HierarchySolver<DIM> class.
   * Although this class hides the lower-level components (FAC cycling,
   * Elastic equation operations and boundary conditions), it is
   * perfectly acceptable to use those lower-level components directly.
   *
   * This class is a wrapper, providing a single class that
   * coordinates three major components: the FAC solver, the
   * cell-centered Elastic FAC operator and a default bc
   * implelemtation.  It is perfectly acceptable to use those
   * classes outside of this class.
   *
   * The underlying solver is an FAC solver using staggered
   * discretization.  The difference scheme is second-order
   * central-difference.
   *
   * Typical use of this class is:
   * -# Construct a Elastic::FACSolver object, providing it
   *    the hierarchy and range of levels participating in the solve.
   * -# Call initializeSolverState() to set up information
   *    internal to the solver.  This is step is not required
   *    but will save setup costs if you are making multiple
   *    solves.  This commits the object to the current hierarchy state
   *    and the specific @em types of boundary conditions you selected,
   *    It does NOT commit to the specific @em values of the boundary
   *    condition.  A hierarchy change (through adaption or other means)
   *    invalidates the state, thus you must reinitialize or
   *    deallocateSolverState() the state before another solve.
   * -# Solve the equation with solveSystem().  You provide the
   *    patch data indices for the solution u and the right hand
   *    side f.  u must have at least one ghost cell and where
   *    a Dirichlet boundary condition applies, those cells
   *    must be set to the value on the boundary.  If only Neumann
   *    boundary conditions are used, the ghost cell values
   *    do not matter.
   * -# Call deallocateSolverState() to free up internal resources,
   *    if initializeSolverState() was called before the solve.
   *
   * After the solve, information on the solve can be obtained
   * by calling one of these functions:
   * - getNumberOfIterations() gives the number of FAC cycles used.
   * - getConvergenceFactors() gives the average and final convergence
   *   factors for the solve.
   * - getResidualNorm() gives the final residual
   *
   * Finer solver controls can be set using the functions in this class.
   *
   * Object of this class can be set using input databases.
   * The following parameters can be set.  Each is shown with its
   * default value in the case where hypre is used.
   * @verbatim
   * enable_logging = TRUE // Bool flag to switch logging on/off
   * max_cycles = 10       // Integer number of max FAC cycles to use
   * residual_tol = 1.e-6  // Residual tolerance to solve for
   * num_pre_sweeps = 1    // Number of presmoothing sweeps to use
   * num_post_sweeps = 1   // Number of postsmoothing sweeps to use
   * coarse_fine_discretization = "Ewing" // Name of coarse-fine discretization
   * prolongation_method = "CONSTANT_REFINE" // Name of prolongation method
   * coarse_solver_choice = "hypre"  // Name of coarse level solver
   * coarse_solver_tolerance = 1e-10 // Coarse level tolerance
   * coarse_solver_max_iterations = 20 // Coarse level max iterations
   * use_smg = "FALSE"     // Whether to use hypre's smg solver
   *                       // (alternative is the pfmg solver)
   * @endverbatim
   *
   */
  class FACSolver
  {

  public:
    /*!
     * @brief Construct a solver.
     *
     * If the database is not NULL, initial settings will be set
     * using the database.
     * The solver is uninitialized until initializeSolverState()
     * is called.
     *
     * @param object_name Name of object used in outputs
     * @param database tbox::Database for initialization (may be NULL)
     */
    FACSolver(const SAMRAI::tbox::Dimension& dim,
              const std::string& object_name,
              boost::shared_ptr<SAMRAI::tbox::Database> database,
              Boundary_Conditions &bc);

    /*!
     * @brief Destructor.
     */
    ~FACSolver(void);

    /*!
     * @brief Enable logging.
     *
     * To disable, pass in @c false.
     */
    void
    enableLogging(bool logging);

    /*!
     * @brief Solve Elastic's equation, assuming an uninitialized
     * solver state.
     *
     * Here, u is the "solution" patch data index and f is the
     * right hand side patch data index.
     * The return value is true if the solver converged and false otherwise.
     *
     * This function is a wrapper.
     * It simply initializes the solver state, call the
     * solveSystem(const int,const int) for the initialized solver then
     * deallocates the solver state.
     *
     * Upon return from this function,
     * solution will contain the result of the solve.
     *
     * See initializeSolverState() for opportunities to save overhead
     * when using multiple consecutive solves.
     *
     * @see solveSystem(const int,const int)
     *
     * @param solution hier::Patch data index for solution u
     * @param rhs hier::Patch data index for right hand side f
     * @param hierarchy The patch hierarchy to solve on
     * @param coarse_ln The coarsest level in the solve.
     * @param fine_ln The finest level in the solve.
     *
     * @return whether solver converged to specified level
     *
     * @see initializeSolverState
     */
    bool
    solveSystem(const int cell_moduli_id,
                const int edge_moduli_id,
                const int dv_diagonal_id,
                const int dv_mixed_id, 
                const int level_set_id,
                const int v_id,
                const int v_rhs_id,
                boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
                int coarse_ln = -1,
                int fine_ln = -1);

    /*!
     * @brief Solve Elastic's equation using the current solver state
     * set by initializeSolverState().
     *
     * When the solver state has been initialized, this function may
     * be called repeadedly with different values on the rhs.
     * There is some cost savings for multiple solves when this
     * is done.
     *
     * Before calling this function, the solution and
     * right-hand-side quantities should be set properly by the user
     * on all patch interiors on the range of levels covered by the
     * FAC iteration.  All data for these patch data index should be allocated.
     * Thus, the user is responsible for managing the
     * storage for the solution and right-hand-side.
     *
     * @return whether solver converged to specified level
     *
     * @see solveSystem( const int, const int, boost::shared_ptr< SAMRAI::hier::PatchHierarchy >, int, int);
     */
    bool
    solveSystem(const int v, const int v_rhs);

    //@{ @name Functions for setting solver mathematic algorithm controls

    /*!
     * @brief Set tolerance for coarse level solve.
     *
     * If the coarse level solver requires a tolerance
     * (currently, they all do), the specified value is used.
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
     * @em not passed to xfer::Geometry::lookupRefineOperator().
     * If this variable is set to "Ewing", a constant refinement
     * method is used along with Ewing's correction.
     * For a reference to the correction method, see
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
     * @param prolongation_method String selecting the coarse-fine discretization method.
     */
    void
    set_V_ProlongationMethod(const std::string& prolongation_method);

    /*!
     * @brief Set the number of pre-smoothing sweeps during
     * FAC iteration process.
     *
     * Presmoothing is applied during the fine-to-coarse phase of the
     * iteration.  The default is to use one sweep.
     *
     * @param num_pre_sweeps Number of presmoothing sweeps
     */
    void
    setPresmoothingSweeps(int num_pre_sweeps);

    /*!
     * @brief Set the number of post-smoothing sweeps during
     * FAC iteration process.
     *
     * Postsmoothing is applied during the coarse-to-fine phase of the
     * iteration.  The default is to use one sweep.
     *
     * @param num_post_sweeps Number of postsmoothing sweeps
     */
    void
    setPostsmoothingSweeps(int num_post_sweeps);

    /*!
     * @brief Set the max number of iterations (cycles) to use per solve.
     */
    void
    setMaxCycles(int max_cycles);

    /*!
     * @brief Set the residual tolerance for stopping.
     *
     * If you want the prescribed maximum number of cycles to always be taken,
     * set the residual tolerance to a negative number.
     */
    void
    setResidualTolerance(double residual_tol);

    //@}

    /*!
     * @brief Prepare the solver's internal state for solving
     *
     * In the interest of efficiency, this class may prepare and
     * cache some hierarchy-dependent objects.  Though it is not required,
     * initializing the solver state makes for greater efficiency
     * when you are doing multiple solves on the same system of
     * equation.  If you do not initialize the state, it is initialized
     * and deallocated each time you call solveSystem(const int, const int).
     * The state must be reinitialized if the hierarchy or a boundary
     * condition type changes.
     *
     * To unset the data set in this function,
     * see deallocateSolverState().
     *
     * The @c solution and @c rhs patch data indices in the argument
     * list are used to determine the @em form of the data you
     * plan to use in the solve.  They need not be the same data
     * you solve on, but they should be similar.  Both must represent
     * cell-centered double data.  The solution must have at least one
     * ghost cell width, though this is not checked in the initialize
     * phase, because data is not required yet.
     *
     * @param solution solution patch data index for u
     * @param rhs right hand side patch data index for f
     * @param hierarchy The patch hierarchy to solve on
     * @param coarse_level The coarsest level in the solve
     * @param fine_level The finest level in the solve
     */
    void
    initializeSolverState(const int cell_moduli_id,
                          const int edge_moduli_id,
                          const int dv_diagonal_id,
                          const int dv_mixed_id, 
                          const int level_set_id,
                          const int v_id,
                          const int v_rhs_id,
                          boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
                          const int coarse_level = -1,
                          const int fine_level = -1);

    /*!
     * @brief Remove the solver's internal state data
     *
     * Remove all hierarchy-dependent data set by initializeSolverState.
     * It is safe to call deallocateSolverState() even state is already
     * deallocated, but nothing is done in that case.
     *
     * @see initializeSolverState()
     */
    void
    deallocateSolverState();

    //@{
    //! @name Functions to get data on last solve.

    /*!
     * @brief Return FAC iteration count from last (or current
     * if there is one) FAC iteration process.
     */
    int
    getNumberOfIterations() const;

    /*!
     * @brief Get average convergance rate and convergence rate of
     * the last (or current if there is one) FAC solve.
     *
     * @param avg_factor average convergence factor over current FAC cycles
     * @param final_factor convergence factor of the last FAC cycle
     */
    void
    getConvergenceFactors(double& avg_factor,
                          double& final_factor) const;

    /*!
     * @brief Return residual norm from the just-completed FAC iteration.
     *
     * The norm return value is computed as the maximum norm over all
     * patch levels involved in the solve.  The value corresponds to the
     * norm applied in the user-defined residual computation.
     *
     * The latest computed norm is the one returned.
     */
    double
    getResidualNorm() const;

    void set_boundaries(const int &v_id,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> &level,
                        const bool &homogeneous)
    {
      d_fac_ops.set_boundaries(v_id,level,homogeneous);
    }


    //@}

  private:
    /*!
     * @brief Set state using database
     *
     * See the class description for the parameters that can be set
     * from a database.
     *
     * @param database Input database.  If a NULL pointer is given,
     * nothing is done.
     */
    void
    getFromInput(boost::shared_ptr<SAMRAI::tbox::Database> database);

    /*
     * @brief Set @c d_uv and @c d_fv to vectors wrapping the data
     * specified by patch data indices u and f.
     */
    void
    createVectorWrappers(int v, int v_rhs);

    /*
     * @brief Destroy vector wrappers referenced to by @c d_uv and @c d_fv.
     */
    void
    destroyVectorWrappers();

    /*
     * @brief Initialize static members
     */
    static void
    initializeStatics();

    const SAMRAI::tbox::Dimension d_dim;

    /*!
     * @brief Object name.
     */
    std::string d_object_name;

    Boundary_Conditions &d_boundary_conditions;
    /*!
     * @brief FAC operator implementation corresponding to staggered
     * Elastic discretization.
     */
    FACOps d_fac_ops;

    /*!
     * @brief FAC preconditioner algorithm.
     */
    SAMRAI::solv::FACPreconditioner d_fac_precond;

    /*!
     * @brief Robin bc object in use.
     */
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;
    int d_ln_min;
    int d_ln_max;

    /*!
     * @brief Context for all internally maintained data.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    /*
     * @brief Vector wrapper for solution.
     * @see createVectorWrappers(), destroyVectorWrappers()
     */
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_uv;
    /*
     * @brief Vector wrapper for source.
     * @see createVectorWrappers(), destroyVectorWrappers()
     */
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_fv;

    bool d_solver_is_initialized;
    bool d_enable_logging;

    static bool s_initialized;
    static int s_weight_id[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
    static int s_instance_counter[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
  };

}

#include "Elastic/FACSolver.I"

#endif  // included_solv_ElasticFACSolver
