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

#include "SAMRAI/solv/CellPoissonFACSolver.h"
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
 *
 * This class demonstrates how use the FAC Stokes solver
 * class to solve Stokes's equation on a single level
 * within a hierarchy.
 *
 * We set up and solve the following problem:
 *
 *   2d: div(grad(u)) = -2 (pi^2) sin(pi x) sin(pi y)
 *
 *   3d: div(grad(u)) = -3 (pi^2) sin(pi x) sin(pi y) sin(pi z)
 *
 * which has the exact solution
 *
 *   2d: u = sin(pi x) sin(pi y)
 *
 *   3d: u = sin(pi x) sin(pi y) sin(pi z)
 *
 * This class inherits and implements virtual functions from
 * - mesh::StandardTagAndInitStrategy to initialize data
 *   on the SAMR grid.
 * - appu::VisDerivedDataStrategy to write out certain data
 *   in a vis file, such as the error of the solution.
 *
 * Inputs:  The only input parameter for this class is
 * "fac_stokes", the input database for the solv::CellStokesFACSolver
 * object.  See the documentation for solv::CellStokesFACSolver
 * for its input parameters.
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
   FACStokes(
      const std::string& object_name,
      const tbox::Dimension& dim,
      tbox::Pointer<tbox::Database> database =
         tbox::Pointer<tbox::Database>(NULL));

   virtual ~FACStokes();

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
   initializeLevelData(
      const tbox::Pointer<hier::BasePatchHierarchy> hierarchy,
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
   resetHierarchyConfiguration(
      tbox::Pointer<hier::BasePatchHierarchy> new_hierarchy,
      int coarsest_level,
      int finest_level);

   //@}

   //@{ @name appu::VisDerivedDataStrategy virtuals

   virtual bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_id) const;

   //@}

   /*!
    * @brief Solve using HYPRE Stokes solver
    *
    * Set up the linear algebra problem and use a
    * solv::CellStokesFACSolver object to solve it.
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
   setupPlotter(
      appu::VisItDataWriter& plotter) const;
#endif

private:
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
   solv::CellPoissonFACSolver d_stokes_fac_solver;

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
   int d_comp_soln_id, d_exact_id, d_rhs_id;

   //@}

};

}

#endif  // included_FACStokes
