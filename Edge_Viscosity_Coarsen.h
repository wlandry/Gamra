/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef included_geom_Edge_Viscosity_Coarsen
#define included_geom_Edge_Viscosity_Coarsen

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"

#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Class Edge_Viscosity_Coarsen implements averaging 
 * for the edge-centered component of the viscosity.
 *
 * The name of this method when specifying it in the input file is
 * "EDGE_VISCOSITY_COARSEN"
 *
 * @see xfer::CoarsenOperator
 */

class Edge_Viscosity_Coarsen:
   public xfer::CoarsenOperator
{
public:
  /**
   * Uninteresting default constructor.
   */
  explicit Edge_Viscosity_Coarsen(const tbox::Dimension& dim,
                                  const int &cell_viscosity):
    xfer::CoarsenOperator(dim, "EDGE_VISCOSITY_COARSEN"),
    cell_viscosity_id(cell_viscosity)
  {
    d_name_id = "EDGE_VISCOSITY_COARSEN";
  }

  /**
   * Uninteresting virtual destructor.
   */
  virtual ~Edge_Viscosity_Coarsen(){}

  /**
   * Return true if the variable and name std::string match the edge-centered
   * double weighted averaging; otherwise, return false.
   */
  
  bool findCoarsenOperator(const tbox::Pointer<hier::Variable>& var,
                           const std::string& op_name) const
  {
    TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

    const tbox::Pointer<pdat::NodeVariable<double> > node_var(var);
    const tbox::Pointer<pdat::EdgeVariable<double> > edge_var(var);
    if ((!node_var.isNull() || !edge_var.isNull()) && (op_name == d_name_id)) {
      return true;
    } else {
      return false;
    }
  }

  /**
   * Return name std::string identifier of this coarsening operator.
   */
  const std::string& getOperatorName() const
  {
    return d_name_id;
  }

  /**
   * The priority of edge-centered double weighted averaging is 0.
   * It will be performed before any user-defined coarsen operations.
   */
  int getOperatorPriority() const
  {
    return 0;
  }

  /**
   * The stencil width of the weighted averaging operator is the vector of
   * zeros.  That is, its stencil does not extend outside the fine box.
   */
  hier::IntVector getStencilWidth() const
  {
    return hier::IntVector::getZero(getDim());
  }

  /**
   * Coarsen the source component on the fine patch to the destination
   * component on the coarse patch using the edge-centered double weighted
   * averaging operator.  Coarsening is performed on the intersection of
   * the destination patch and the coarse box.  It is assumed that the
   * fine patch contains sufficient data for the stencil width of the
   * coarsening operator.
   */
  void
  coarsen(
          hier::Patch& coarse,
          const hier::Patch& fine,
          const int dst_component,
          const int src_component,
          const hier::Box& coarse_box,
          const hier::IntVector& ratio) const;

private:
  std::string d_name_id;
  int cell_viscosity_id;
};

}
}
#endif
