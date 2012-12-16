/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef GAMRA_ELASTIC_V_REFINE_H
#define GAMRA_ELASTIC_V_REFINE_H

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "FTensor.hpp"
#include <string>

namespace Elastic {

  /**
   * Class V_Refine implements linear
   * interpolation for side-centered double patch data defined over a Cartesian
   * mesh.  It is derived from the hier::RefineOperator base class.
   * CartesianSideDoubleConservativeLinearRefine does not handle the 
   * boundary for the velocity correctly, so use this instead for
   * velocity.
   *
   * The findRefineOperator() operator function returns true if the input
   * variable is side-centered double, and the std::string is "V_REFINE".
   *
   * @see hier::RefineOperator
   */

  class V_Refine:
    public SAMRAI::hier::RefineOperator
  {
  public:
    /**
     * Uninteresting default constructor.
     */
    explicit V_Refine(const SAMRAI::tbox::Dimension& dim):
      SAMRAI::hier::RefineOperator(dim, "V_REFINE")
    {
      d_name_id = "V_REFINE";
    }


    /**
     * Uninteresting virtual destructor.
     */
    virtual ~V_Refine(){}

    /**
     * Return true if the variable and name std::string match side-centered
     * double linear interpolation; otherwise, return false.
     */
    bool
    findRefineOperator(const boost::shared_ptr<SAMRAI::hier::Variable>& var,
                       const std::string& op_name) const
    {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

      const boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
        cast_var(boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >(var));
      if (cast_var && (op_name == d_name_id)) {
        return true;
      } else {
        return false;
      }
    }
    /**
     * Return name std::string identifier of this refinement operator.
     */
    const std::string& getOperatorName() const
    {
      return d_name_id;
    }

    /**
     * The priority of side-centered double linear interpolation is 0.
     * It will be performed before any user-defined interpolation operations.
     */
    int getOperatorPriority() const
    {
      return 0;
    }

    /**
     * The stencil width of the linear interpolation operator is the vector
     * of ones.  That is, its stencil extends one side outside the fine box.
     */
    SAMRAI::hier::IntVector getStencilWidth() const
    {
      return SAMRAI::hier::IntVector::getOne(getDim());
    }

    /**
     * Refine the source component on the coarse patch to the destination
     * component on the fine patch using the side-centered double linear
     * interpolation operator.  Interpolation is performed on the intersection
     * of the destination patch and the boxes contained in fine_overlap
     * It is assumed that the coarse patch contains sufficient data for the
     * stencil width of the refinement operator.
     */
    void refine(SAMRAI::hier::Patch& fine,
                const SAMRAI::hier::Patch& coarse,
                const int dst_component,
                const int src_component,
                const SAMRAI::hier::BoxOverlap& fine_overlap,
                const SAMRAI::hier::IntVector& ratio) const;

    /**
     * Refine the source component on the coarse patch to the destination
     * component on the fine patch using the side-centered double linear
     * interpolation operator.  Interpolation is performed on the intersection
     * of the destination patch and the fine box.   It is assumed that the
     * coarse patch contains sufficient data for the stencil width of the
     * refinement operator.  This differs from the above refine() method
     * only in that it operates on a single fine box instead of a BoxOverlap.
     */
    void refine(SAMRAI::hier::Patch& fine,
                const SAMRAI::hier::Patch& coarse,
                const int dst_component,
                const int src_component,
                const SAMRAI::hier::Box& fine_box,
                const SAMRAI::hier::IntVector& ratio,
                const int &axis) const;

    double refine_along_line
    (SAMRAI::pdat::SideData<double> &v,
     const int &axis,
     const int &dim,
     const SAMRAI::hier::Index pp[],
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::pdat::SideIndex &coarse,
     const SAMRAI::hier::Box &coarse_box,
     const SAMRAI::geom::CartesianPatchGeometry &coarse_geom) const;

  private:
    std::string d_name_id;

  };

}
#endif
