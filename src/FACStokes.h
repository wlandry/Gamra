/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#ifndef included_FACStokes
#define included_FACStokes

#include "StokesFACSolver.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

namespace SAMRAI {

  /*!
   * @brief Class to solve a sample Stokes equation on a SAMR grid.
   */
  class FACStokes:
    public mesh::StandardTagAndInitStrategy,
    public appu::VisDerivedDataStrategy
  {

  public:
    /*!
     * @brief Constructor.
     *
     * If you want standard output and logging,
     * pass in valid pointers for those streams.
     *
     * @param object_name Ojbect name
     * @param database Input database (may be NULL)
     */
    FACStokes(const std::string& object_name,
              const tbox::Dimension& dim,
              tbox::Pointer<tbox::Database> database =
              tbox::Pointer<tbox::Database>(NULL));

    virtual ~FACStokes() {}

    //@{ @name mesh::StandardTagAndInitStrategy virtuals

    /*!
     * @brief Allocate and initialize data for a new level
     * in the patch hierarchy.
     *
     * This is where you implement the code for initialize data on
     * the grid.  All the information needed to initialize the grid
     * are in the arguments.
     *
     * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
     */
    virtual void
    initializeLevelData(const tbox::Pointer<hier::BasePatchHierarchy> hierarchy,
                        const int level_number,
                        const double init_data_time,
                        const bool can_be_refined,
                        const bool initial_time,
                        const tbox::Pointer<hier::BasePatchLevel> old_level,
                        const bool allocate_data);

    /*!
     * @brief Reset any internal hierarchy-dependent information.
     */
    virtual void
    resetHierarchyConfiguration(tbox::Pointer<hier::BasePatchHierarchy> new_hierarchy,
                                int coarsest_level,
                                int finest_level);

    //@}

    virtual void
    applyGradientDetector(const tbox::Pointer<hier::BasePatchHierarchy> hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const int tag_index,
                          const bool initial_time,
                          const bool uses_richardson_extrapolation);

    void computeAdaptionEstimate(pdat::CellData<double>& estimate_data,
                                 const pdat::CellData<double>& soln_cell_data)
      const;

    //@{ @name appu::VisDerivedDataStrategy virtuals

    virtual bool
    packDerivedDataIntoDoubleBuffer(double* buffer,
                                    const hier::Patch& patch,
                                    const hier::Box& region,
                                    const std::string& variable_name,
                                    int depth_id) const;

    //@}

    /*!
     * @brief Solve using HYPRE Stokes solver
     *
     * Set up the linear algebra problem and use a
     * solv::StokesFACSolver object to solve it.
     * -# Set initial guess
     * -# Set boundary conditions
     * -# Specify Stokes equation parameters
     * -# Call solver
     */
    int
    solveStokes();

#ifdef HAVE_HDF5
    /*!
     * @brief Set up external plotter to plot internal
     * data from this class.
     *
     * After calling this function, the external
     * data writer may be used to write the
     * viz file for this object.
     *
     * The internal hierarchy is used and must be
     * established before calling this function.
     * (This is commonly done by building a hierarchy
     * with the mesh::StandardTagAndInitStrategy virtual
     * functions implemented by this class.)
     *
     * @param viz_writer VisIt writer
     */
    int
    setupPlotter(appu::VisItDataWriter& plotter) const;
#endif

  private:
    void fix_viscosity();
    std::string d_object_name;

    const tbox::Dimension d_dim;

    tbox::Pointer<hier::PatchHierarchy> d_hierarchy;

    //@{
    /*!
     * @name Major algorithm objects.
     */

    /*!
     * @brief FAC stokes solver.
     */
    solv::StokesFACSolver d_stokes_fac_solver;

    /*!
     * @brief Boundary condition coefficient implementation.
     */
    solv::LocationIndexRobinBcCoefs d_bc_coefs;

    //@}

    //@{

    /*!
     * @name Private state variables for solution.
     */

    /*!
     * @brief Context owned by this object.
     */
    tbox::Pointer<hier::VariableContext> d_context;
  
    /*!
     * @brief Descriptor indices of internal data.
     *
     * These are initialized in the constructor and never change.
     */

    double d_adaption_threshold;
    int min_full_refinement_level;
  public:
    int p_id, cell_viscosity_id, edge_viscosity_id, dp_id, p_exact_id,
      p_rhs_id, v_id, v_rhs_id;

    tbox::Array<double> viscosity, viscosity_xyz_max, viscosity_xyz_min;
    tbox::Array<int> viscosity_ijk;

    tbox::Array<double> v_rhs, v_rhs_xyz_max, v_rhs_xyz_min;
    tbox::Array<int> v_rhs_ijk;
    //@}

  };

}

#endif  // included_FACStokes
