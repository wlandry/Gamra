/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for cell-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef included_geom_P_Boundary_Refine
#define included_geom_P_Boundary_Refine

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Class P_Boundary_Refine implements the special interpolation needed
 * for the boundary elements on the coarse-fine interface.
 *
 * The findRefineOperator() operator function returns true if the input
 * variable is cell-centered double, and the std::string is "P_BOUNDARY_REFINE".
 *
 * @see xfer::RefineOperator
 */

class P_Boundary_Refine:
  public xfer::RefineOperator
{
public:
  /**
   * Uninteresting default constructor.
   */
  explicit P_Boundary_Refine(const tbox::Dimension& dim):
    xfer::RefineOperator(dim, "P_BOUNDARY_REFINE")
  {
    d_name_id = "P_BOUNDARY_REFINE";
  }


  /**
   * Uninteresting virtual destructor.
   */
  virtual ~P_Boundary_Refine(){}

  /**
   * Return true if the variable and name std::string match cell-centered
   * double linear interpolation; otherwise, return false.
   */
  bool findRefineOperator(const tbox::Pointer<hier::Variable>& var,
                          const std::string& op_name) const
  {
    TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

    const tbox::Pointer<pdat::CellVariable<double> > cast_var(var);
    if (!cast_var.isNull() && (op_name == d_name_id)) {
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
   * The priority of cell-centered double linear interpolation is 0.
   * It will be performed before any user-defined interpolation operations.
   */
  int getOperatorPriority() const
  {
    return 0;
  }

  /**
   * The stencil width of the linear interpolation operator is the vector
   * of ones.  That is, its stencil extends one cell outside the fine box.
   */
  hier::IntVector getStencilWidth() const
  {
    return hier::IntVector::getOne(getDim());
  }

  /**
   * Refine the source component on the coarse patch to the destination
   * component on the fine patch using the cell-centered double linear
   * interpolation operator.  Interpolation is performed on the intersection
   * of the destination patch and the boxes contained in fine_overlap
   * It is assumed that the coarse patch contains sufficient data for the
   * stencil width of the refinement operator.
   */
  virtual void refine(hier::Patch& fine,
                      const hier::Patch& coarse,
                      const int dst_component,
                      const int src_component,
                      const hier::BoxOverlap& fine_overlap,
                      const hier::IntVector& ratio) const;


private:
  std::string d_name_id;

  void Update_P(const pdat::CellIndex &fine,
                const hier::Index &ip, const hier::Index &jp,
                int &j,
                tbox::Pointer<SAMRAI::pdat::CellData<double> > &p,
                tbox::Pointer<SAMRAI::pdat::CellData<double> > &p_fine) const;
};

}
}
#endif
