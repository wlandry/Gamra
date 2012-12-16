/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition support on cartesian grids. 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_V_REFINE_PATCH_STRATEGY_H
#define GAMRA_ELASTIC_V_REFINE_PATCH_STRATEGY_H

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "Constants.h"
#include "Boundary_Conditions.h"

namespace Elastic {

/*!
 * @brief Helper utility for setting boundary conditions on V.
 *
 * This class inherits and implements virtual functions from
 * xfer::RefinePatchStrategy so it may be used to help create
 * communication schedules if desired.
 */
class V_Refine_Patch_Strategy:
   public SAMRAI::xfer::RefinePatchStrategy
{

public:
   /*!
    * @brief Constructor using.
    *
    * @param object_name Name of the object, for general referencing.
    * @param coef_strategy Coefficients strategy being helped.
    */
   V_Refine_Patch_Strategy(
      const SAMRAI::tbox::Dimension& dim,
      std::string object_name, Boundary_Conditions &bc):
     SAMRAI::xfer::RefinePatchStrategy(dim), d_dim(dim), d_object_name(object_name),
     d_boundary_conditions(bc) {}

   /*!
    * @brief Destructor.
    */
   virtual ~V_Refine_Patch_Strategy(void) {}

   //@{ @name SAMRAI::xfer::RefinePatchStrategy virtuals

   virtual void
   setPhysicalBoundaryConditions(
      SAMRAI::hier::Patch& patch,
      const double ,
      const SAMRAI::hier::IntVector& )
  {
    /* We only set Dirichlet boundaries once, before the solve. */
    d_boundary_conditions.set_boundary(patch,v_id,d_homogeneous_bc);
  }
   SAMRAI::hier::IntVector
   getRefineOpStencilWidth() const
  {
    return SAMRAI::hier::IntVector::getOne(d_dim);
  }
  virtual void
  preprocessRefine(SAMRAI::hier::Patch& ,
                   const SAMRAI::hier::Patch& coarse,
                   const SAMRAI::hier::Box& ,
                   const SAMRAI::hier::IntVector& )
  {
    d_boundary_conditions.set_boundary(coarse,v_id,d_homogeneous_bc);
  }

  virtual void
  postprocessRefineBoxes(SAMRAI::hier::Patch& ,
                         const SAMRAI::hier::Patch& ,
                         const SAMRAI::hier::BoxContainer& ,
                         const SAMRAI::hier::IntVector& ) {}
  virtual void
  postprocessRefine(SAMRAI::hier::Patch& ,
                    const SAMRAI::hier::Patch& ,
                    const SAMRAI::hier::Box& ,
                    const SAMRAI::hier::IntVector& ) {}

   //@}

   /*!
    * @brief Set the data id that should be filled when setting
    * physical boundary conditions.
    *
    * When setPhysicalBoundaryConditions is called, the data
    * specified will be set.  This information is required because
    * the it is not passed in through the argument list of
    * setPhysicalBounaryConditions.
    */
  void setTargetDataId(int id)
  {
    v_id=id;
  }

   /*!
    * @brief Set whether boundary filling should assume homogeneous
    * conditions.
    *
    * In certain circumstances, only the value of a is needed, while
    * the value of g is temporarily not required and taken to be zero.
    * (An example is in setting the boundary condition for error
    * value in an iterative method.)  In such cases, use this function
    * to set a flag that will cause a null pointer to be given to
    * setBcCoefs() to indicate that fact.
    */
  void
  setHomogeneousBc(bool homogeneous_bc)
  {
    d_homogeneous_bc=homogeneous_bc;
  }

   //@}

private:
   /*!
    * @brief Trim a boundary box so that it does not stick out
    * past a given box.
    *
    * Certain boundary-related operations occur on patch data that
    * do not or cannot extend past the edgr or corner of a patch.
    * This function is used to trim down the parts of the boundary box
    * that extend past those points so that a suitable index range
    * is achieved.
    *
    * The boundary box trimmed must be of type 1 or 2.
    *
    * @param boundary_box Boundary box to be trimmed.
    * @param limit_box SAMRAI::hier::Box to not stick past
    *
    * @return New trimmed boundary box.
    */
   // SAMRAI::hier::BoundaryBox
   // trimBoundaryBox(
   //    const SAMRAI::hier::BoundaryBox& boundary_box,
   //    const SAMRAI::hier::Box& limit_box) const;

   /*!
    * @brief Return box describing the index space of boundary nodes
    * defined by a boundary box.
    *
    * Define a box describing the indices of the nodes corresponding
    * to the input boundary box.  These nodes lie on the boundary
    * itself.
    *
    * The input boundary_box must be of type 1
    * (see SAMRAI::hier::BoundaryBox::getBoundaryType()).
    *
    * @param boundary_box input boundary box
    * @return a box to define the node indices corresponding to
    *   boundary_box
    */
   // SAMRAI::hier::Box
   // makeNodeBoundaryBox(
   //    const SAMRAI::hier::BoundaryBox& boundary_box) const;

   /*!
    * @brief Return box describing the index space of faces
    * defined by a boundary box.
    *
    * Define a box describing the indices of the codimension 1
    * surface corresponding to the input boundary box.
    *
    * The input boundary_box must be of type 1
    * (see SAMRAI::hier::BoundaryBox::getBoundaryType()).
    *
    * This is a utility function for working with the
    * indices coresponding to a boundary box but coincide
    * with the patch boundary.
    *
    * @param boundary_box input boundary box
    * @return a box to define the face indices corresponding to
    *    boundary_box
    */
   // SAMRAI::hier::Box
   // makeFaceBoundaryBox(
   //    const SAMRAI::hier::BoundaryBox& boundary_box) const;

   const SAMRAI::tbox::Dimension d_dim;

   std::string d_object_name;

   /*!
    * @brief Coefficient strategy giving a way to get to
    * Robin bc coefficients.
    */
   // const RobinBcCoefStrategy* d_coef_strategy;

   /*!
    * @brief SAMRAI::hier::Index of target patch data when filling ghosts.
    */
   int v_id;

  Boundary_Conditions &d_boundary_conditions;

   /*!
    * @brief Whether to assumg g=0 when filling ghosts.
    */
   bool d_homogeneous_bc;

   /*!
    * @brief Timers for performance measurement.
    */
   // boost::shared_ptr<tbox::Timer> t_set_boundary_values_in_cells;
   // boost::shared_ptr<tbox::Timer> t_use_set_bc_coefs;
};

}
#endif  // included_solv_V_Refine_Patch_Strategy
