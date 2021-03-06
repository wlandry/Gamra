/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Stokes solver 
 *
 ************************************************************************/
#pragma once

#include "Stokes/FACSolver.hxx"
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/hier/Box.h>
#include <SAMRAI/solv/LocationIndexRobinBcCoefs.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/appu/VisDerivedDataStrategy.h>
#include <SAMRAI/appu/VisItDataWriter.h>

namespace Stokes
{
  class FAC:
    public SAMRAI::mesh::StandardTagAndInitStrategy,
    public SAMRAI::appu::VisDerivedDataStrategy
  {
  public:
    FAC(const SAMRAI::tbox::Dimension& dim, SAMRAI::tbox::Database &database);
    virtual ~FAC() {}

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
    initializeLevelData(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                        const int level_number,
                        const double init_data_time,
                        const bool can_be_refined,
                        const bool initial_time,
                        const boost::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
                        const bool allocate_data);

    /*!
     * @brief Reset any internal hierarchy-dependent information.
     */
    virtual void
    resetHierarchyConfiguration(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& new_hierarchy,
                                const int coarsest_level,
                                const int finest_level);

    //@}

    virtual void
    applyGradientDetector(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy> &hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const int tag_index,
                          const bool initial_time,
                          const bool uses_richardson_extrapolation);

    //@{ @name appu::VisDerivedDataStrategy virtuals

    virtual bool
    packDerivedDataIntoDoubleBuffer(double* buffer,
                                    const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& region,
                                    const std::string& variable_name,
                                    int depth_id,
                                    double simulation_time = 0.0) const;

    //@}

    /*!
     * Set up the linear algebra problem and use a
     * Stokes::FACSolver object to solve it.
     * -# Set initial guess
     * -# Set boundary conditions
     * -# Specify Stokes equation parameters
     * -# Call solver
     */
    bool solve();

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
    void
    setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const;
#endif

  private:
    void fix_viscosity();
    const SAMRAI::tbox::Dimension d_dim;

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

    //@{
    /*!
     * @name Major algorithm objects.
     */

    /*!
     * @brief FAC stokes solver.
     */
    Stokes::FACSolver d_stokes_fac_solver;

    /*!
     * @brief Boundary condition coefficient implementation.
     */
    SAMRAI::solv::LocationIndexRobinBcCoefs d_bc_coefs;

    //@}

    //@{

    /*!
     * @name Private state variables for solution.
     */

    /*!
     * @brief Context owned by this object.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
  
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

    std::vector<double> viscosity, viscosity_xyz_max, viscosity_xyz_min;
    std::vector<int> viscosity_ijk;

    std::vector<double> v_rhs, v_rhs_xyz_max, v_rhs_xyz_min;
    std::vector<int> v_rhs_ijk;

    std::vector<double> p_initial, p_initial_xyz_max, p_initial_xyz_min;
    std::vector<int> p_initial_ijk;
    //@}

  };
}


