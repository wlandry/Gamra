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

#ifndef GAMRA_ELASTIC_V_BOUNDARY_REFINE_H
#define GAMRA_ELASTIC_V_BOUNDARY_REFINE_H

#include <SAMRAI/SAMRAI_config.h>

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include "Boundary_Conditions.hxx"
#include "quad_offset_interpolate.hxx"

#include <string>

namespace Elastic {

  /**
   * Class V_Boundary_Refine implements the special interpolation needed
   * for the boundary elements on the coarse-fine interface.
   *
   * The findRefineOperator() operator function returns true if the input
   * variable is side-centered double, and the std::string is "V_BOUNDARY_REFINE".
   *
   * @see hier::RefineOperator
   */

  class V_Boundary_Refine:
    public SAMRAI::hier::RefineOperator
  {
  public:

    static bool is_residual;
    static int dv_diagonal_id, dv_mixed_id, level_set_id;
    /**
     * Uninteresting default constructor.
     */
    explicit V_Boundary_Refine():
      SAMRAI::hier::RefineOperator("V_BOUNDARY_REFINE")
    {
      d_name_id = "V_BOUNDARY_REFINE";
    }

    bool have_faults() const
    {
      return dv_diagonal_id!=invalid_id;
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }
    /**
     * Uninteresting virtual destructor.
     */
    virtual ~V_Boundary_Refine(){}

    /**
     * Return true if the variable and name std::string match side-centered
     * double linear interpolation; otherwise, return false.
     */
    bool findRefineOperator
    (const boost::shared_ptr<SAMRAI::hier::Variable>& var,
     const std::string& op_name) const
    {
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
    SAMRAI::hier::IntVector getStencilWidth(const SAMRAI::tbox::Dimension &dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }

    /**
     * Refine the source component on the coarse patch to the destination
     * component on the fine patch using the side-centered double linear
     * interpolation operator.  Interpolation is performed on the intersection
     * of the destination patch and the boxes contained in fine_overlap
     * It is assumed that the coarse patch contains sufficient data for the
     * stencil width of the refinement operator.
     */
    void refine(SAMRAI::hier::Patch &fine_patch,
                const SAMRAI::hier::Patch &coarse_patch,
                const int dst_component,
                const int src_component,
                const SAMRAI::hier::BoxOverlap &fine_overlap,
                const SAMRAI::hier::IntVector &ratio) const;

    /**
     * Refine the source component on the coarse patch to the destination
     * component on the fine patch using the side-centered double linear
     * interpolation operator.  Interpolation is performed on the intersection
     * of the destination patch and the fine box.   It is assumed that the
     * coarse patch contains sufficient data for the stencil width of the
     * refinement operator.  This differs from the above refine() method
     * only in that it operates on a single fine box instead of a BoxOverlap.
     */
    void refine(SAMRAI::hier::Patch &fine_patch,
                const SAMRAI::hier::Patch &coarse_patch,
                const int dst_component,
                const int src_component,
                const SAMRAI::hier::Box &fine_box,
                const SAMRAI::hier::IntVector &ratio,
                const Gamra::Dir &axis) const;

  private:
    std::string d_name_id;

    void Update_V_embedded_2D
    (const int &axis,
     const int &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector &ip,
     const SAMRAI::hier::IntVector &jp,
     const int &i,
     const int &j,
     const SAMRAI::pdat::SideData<double> &level_set,
     const SAMRAI::pdat::SideData<double> &level_set_fine,
     const SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_fine) const;

    void Update_V_2D
    (const int &axis,
     const int &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector &ip,
     const SAMRAI::hier::IntVector &jp,
     const int &i,
     const int &j,
     const SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_fine) const;

    void Correction_2D
    (const int &axis,
     const int &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector &ip,
     const SAMRAI::hier::IntVector &jp,
     const int &i,
     const int &j,
     const int &i_min,
     const int &i_max, 
     const SAMRAI::pdat::CellData<double> &dv_diagonal,
     const SAMRAI::pdat::CellData<double> &dv_diagonal_fine,
     const SAMRAI::pdat::SideData<double> &dv_mixed,
     const SAMRAI::pdat::SideData<double> &dv_mixed_fine,
     SAMRAI::pdat::SideData<double> &v_fine) const;

    void Update_V_3D
    (const Gamra::Dir &ix,
     const Gamra::Dir &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector unit[],
     const SAMRAI::hier::Index &ijk,
     const SAMRAI::hier::Box &coarse_box,
     const SAMRAI::hier::Index &fine_min,
     const SAMRAI::hier::Index &fine_max,
     const SAMRAI::geom::CartesianPatchGeometry &geom,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> > &dv_diagonal,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> > &dv_diagonal_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > &dv_mixed_fine,
     const SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_fine) const;

    void quad_offset_correction(const double &dv_pm, const double &dv_p,
                                const double &dv_m, const double &dv_mp,
                                double &correction_p, double &correction_m) const
    {
      quad_offset_interpolate(dv_pm-dv_p,0,dv_mp-dv_m,correction_p,correction_m);
    }

    double quad_offset_correction(const double &dv_pm, const double &dv_p,
                                  const double &dv_m, const double &dv_mp) const
    {
      return quad_offset_interpolate(dv_pm-dv_p,0,dv_mp-dv_m);
    }
  };

}
#endif
