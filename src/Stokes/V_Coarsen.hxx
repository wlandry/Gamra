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

#pragma once

#include <SAMRAI/SAMRAI_config.h>

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/pdat/SideVariable.h>

#include <string>

namespace Stokes {

  /**
   * Class V_Coarsen implements averaging 
   * for side-centered double patch data defined over
   * a Cartesian mesh.  It is derived from the SAMRAI::hier::CoarsenOperator base class.
   * The numerical operations for theaveraging use FORTRAN numerical routines.
   *
   * CartesianSideDoubleWeightedAverage averages over the nearest two
   * cells, which is not what we want for multigrid.  This averages over
   * the nearest six points (in 2D).  The findCoarsenOperator() operator
   * function returns true if the input variable is side-centered
   * double, and the std::string is "V_COARSEN".
   *
   * @see SAMRAI::hier::CoarsenOperator
   */

  class V_Coarsen:
    public SAMRAI::hier::CoarsenOperator
  {
  public:
    /**
     * Uninteresting default constructor.
     */
    explicit V_Coarsen():
      SAMRAI::hier::CoarsenOperator("V_COARSEN")
    {
      d_name_id = "V_COARSEN";
    }

    /**
     * Uninteresting virtual destructor.
     */
    virtual ~V_Coarsen(){}

    /**
     * Return true if the variable and name std::string match the side-centered
     * double weighted averaging; otherwise, return false.
     */
  
    bool findCoarsenOperator
    (const boost::shared_ptr<SAMRAI::hier::Variable>& var,
     const std::string& op_name) const
    {
      const boost::shared_ptr<SAMRAI::pdat::SideVariable<double> >
        cast_var(boost::dynamic_pointer_cast<SAMRAI::pdat::SideVariable<double> >
                 (var));
      if (cast_var && (op_name == d_name_id)) {
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
     * The priority of side-centered double weighted averaging is 0.
     * It will be performed before any user-defined coarsen operations.
     */
    int getOperatorPriority() const
    {
      return 0;
    }

    SAMRAI::hier::IntVector getStencilWidth(const SAMRAI::tbox::Dimension& dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }

    /**
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch using the side-centered double weighted
     * averaging operator.  Coarsening is performed on the intersection of
     * the destination patch and the coarse box.  It is assumed that the
     * fine patch contains sufficient data for the stencil width of the
     * coarsening operator.
     */
    void
    coarsen(SAMRAI::hier::Patch& coarse,
            const SAMRAI::hier::Patch& fine,
            const int dst_component,
            const int src_component,
            const SAMRAI::hier::Box& coarse_box,
            const SAMRAI::hier::IntVector& ratio) const
    {
      if(fine.getDim().getValue()==2)
        coarsen_2D(coarse,fine,dst_component,src_component,coarse_box,ratio);
      else
        coarsen_3D(coarse,fine,dst_component,src_component,coarse_box,ratio);
    }

    void
    coarsen_2D(SAMRAI::hier::Patch& coarse,
               const SAMRAI::hier::Patch& fine,
               const int dst_component,
               const int src_component,
               const SAMRAI::hier::Box& coarse_box,
               const SAMRAI::hier::IntVector& ratio) const;

    void
    coarsen_3D(SAMRAI::hier::Patch& coarse,
               const SAMRAI::hier::Patch& fine,
               const int dst_component,
               const int src_component,
               const SAMRAI::hier::Box& coarse_box,
               const SAMRAI::hier::IntVector& ratio) const;

  private:
    std::string d_name_id;

  };

}

