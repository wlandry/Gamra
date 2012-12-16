/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_FAC_H
#define GAMRA_ELASTIC_FAC_H

#include "Elastic/FACSolver.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "Elastic/Boundary_Conditions.h"
#include "Input_Expression.h"

namespace Elastic {
  /*!
   * @brief Class to solve a sample Elastic equation on a SAMR grid.
   */
  class FAC:
    public SAMRAI::mesh::StandardTagAndInitStrategy,
    public SAMRAI::appu::VisDerivedDataStrategy
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
    FAC(const std::string& object_name,
        const SAMRAI::tbox::Dimension& dim,
        boost::shared_ptr<SAMRAI::tbox::Database> database =
        boost::shared_ptr<SAMRAI::tbox::Database>());

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
    initializeLevelData(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>&
                        hierarchy,
                        const int level_number,
                        const double init_data_time,
                        const bool can_be_refined,
                        const bool initial_time,
                        const boost::shared_ptr<SAMRAI::hier::PatchLevel>&
                        old_level,
                        const bool allocate_data);

    /*!
     * @brief Reset any internal hierarchy-dependent information.
     */
    virtual void
    resetHierarchyConfiguration(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>&
                                new_hierarchy,
                                int coarsest_level,
                                int finest_level);

    //@}

    virtual void
    applyGradientDetector(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const int tag_index,
                          const bool initial_time,
                          const bool uses_richardson_extrapolation);

    void
    computeAdaptionEstimate(SAMRAI::pdat::CellData<double>& estimate_data,
                            const SAMRAI::pdat::CellData<double>& soln_cell_data)
      const;

    //@{ @name appu::VisDerivedDataStrategy virtuals

    virtual bool
    packDerivedDataIntoDoubleBuffer(double* buffer,
                                    const SAMRAI::hier::Patch& patch,
                                    const SAMRAI::hier::Box& region,
                                    const std::string& variable_name,
                                    int depth_id) const;

    //@}

    /*!
     * @brief Solve using HYPRE Elastic solver
     *
     * Set up the linear algebra problem and use a
     * solv::Elastic::FACSolver object to solve it.
     * -# Set initial guess
     * -# Set boundary conditions
     * -# Specify Elastic equation parameters
     * -# Call solver
     */
    int
    solve();

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
    setupPlotter(SAMRAI::appu::VisItDataWriter& plotter) const;
#endif

  private:
    void fix_moduli();
    std::string d_object_name;

    const SAMRAI::tbox::Dimension d_dim;

    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

    //@{
    /*!
     * @name Major algorithm objects.
     */
  public:
    Boundary_Conditions d_boundary_conditions;
  private:
    /*!
     * @brief FAC Elastic solver.
     */
    Elastic::FACSolver d_elastic_fac_solver;

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
    int cell_moduli_id, edge_moduli_id, v_id, v_rhs_id;

    Input_Expression lambda, mu, v_rhs;

    SAMRAI::tbox::Array<double> faults;
    //@}

  };

}

#endif  // included_FACElastic
