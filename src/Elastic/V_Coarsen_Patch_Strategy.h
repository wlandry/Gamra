/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Robin boundary condition support on cartesian grids. 
 *
 ************************************************************************/
#ifndef GAMRA_ELASTIC_V_COARSEN_PATCH_STRATEGY_H
#define GAMRA_ELASTIC_V_COARSEN_PATCH_STRATEGY_H

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "Constants.h"
#include "Boundary_Conditions.h"

namespace Elastic {

  /*!
   * @brief Helper utility for setting boundary conditions on V.
   *
   * This class inherits and implements virtual functions from
   * xfer::CoarsenPatchStrategy so it may be used to help create
   * communication schedules if desired.
   */
  class V_Coarsen_Patch_Strategy:
    public SAMRAI::xfer::CoarsenPatchStrategy
  {

  public:
    /*!
     * @brief Constructor using.
     *
     * @param object_name Name of the object, for general referencing.
     * @param coef_strategy Coefficients strategy being helped.
     */
    V_Coarsen_Patch_Strategy(const SAMRAI::tbox::Dimension& dim,
                             std::string object_name,
                             Boundary_Conditions &bc):
      SAMRAI::xfer::CoarsenPatchStrategy(dim),
      d_object_name(object_name),
      dv_diagonal_id(invalid_id),
      dv_mixed_id(invalid_id),
      level_set_id(invalid_id),
      d_boundary_conditions(bc){}

    /*!
     * @brief Destructor.
     */
    virtual ~V_Coarsen_Patch_Strategy(void) {}

    //@{ @name SAMRAI::xfer::CoarsenPatchStrategy virtuals

    virtual SAMRAI::hier::IntVector
    getCoarsenOpStencilWidth() const
    { return SAMRAI::hier::IntVector::getOne(getDim()); }

    virtual void
    preprocessCoarsen(SAMRAI::hier::Patch& ,
                      const SAMRAI::hier::Patch& ,
                      const SAMRAI::hier::Box& ,
                      const SAMRAI::hier::IntVector& ) {}
    virtual void
    postprocessCoarsen(SAMRAI::hier::Patch& coarse,
                       const SAMRAI::hier::Patch& fine,
                       const SAMRAI::hier::Box& coarse_box,
                       const SAMRAI::hier::IntVector& ratio);

    void
    postprocessCoarsen_2D(SAMRAI::hier::Patch& coarse,
                          const SAMRAI::hier::Patch& fine,
                          const SAMRAI::hier::Box& coarse_box,
                          const SAMRAI::hier::IntVector& ratio);

    void
    postprocessCoarsen_3D(SAMRAI::hier::Patch& coarse,
                          const SAMRAI::hier::Patch& fine,
                          const SAMRAI::hier::Box& coarse_box,
                          const SAMRAI::hier::IntVector& ratio);


    void
    coarsen_2D
    (SAMRAI::pdat::SideData<double>& v,
     const SAMRAI::pdat::SideData<double>& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal,
     const SAMRAI::geom::CartesianPatchGeometry& coarse_geom,
     const SAMRAI::hier::Box& coarse_box) const;

    void
    fix_boundary_elements_2D
    (SAMRAI::pdat::SideData<double>& v,
     const SAMRAI::pdat::SideData<double>& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
     const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>& boundaries) const;

    void
    coarsen_3D
    (SAMRAI::pdat::SideData<double>& v,
     const SAMRAI::pdat::SideData<double>& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal,
     const SAMRAI::geom::CartesianPatchGeometry& coarse_geom,
     const SAMRAI::hier::Box& coarse_box) const;

    void
    fix_boundary_elements_3D
    (SAMRAI::pdat::SideData<double>& v,
     const SAMRAI::pdat::SideData<double>& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
     const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>& boundaries) const;


    double coarsen_plane(const SAMRAI::pdat::SideData<double>& v_fine,
                         const SAMRAI::pdat::SideIndex &fine,
                         const SAMRAI::hier::Index jp,
                         const SAMRAI::hier::Index kp) const
    {
      return (v_fine(fine) + v_fine(fine+jp)
              + v_fine(fine+kp) + v_fine(fine+jp+kp))/4;
    }

    double coarsen_plane_correction
    (const SAMRAI::pdat::SideData<double>& dv_mixed,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::Index jp,
     const SAMRAI::hier::Index kp) const
    {
      /* The numbering here (4,7,5,6) is determined by
         the numbering used in FAC::add_faults */
      if(have_faults() && !is_residual)
        return (dv_mixed(fine,4) + dv_mixed(fine+jp,7)
                + dv_mixed(fine+kp,5) + dv_mixed(fine+jp+kp,6))/4;
      return 0;
    }

    double coarsen_point_3D(const SAMRAI::pdat::SideData<double>& v_fine,
                            const boost::shared_ptr<SAMRAI::pdat::SideData<double> > dv_mixed,
                            const boost::shared_ptr<SAMRAI::pdat::CellData<double> > dv_diagonal,
                            const SAMRAI::pdat::SideIndex& fine,
                            const int& axis,
                            const SAMRAI::hier::Index& ip,
                            const SAMRAI::hier::Index& jp,
                            const SAMRAI::hier::Index& kp) const
    {
      double result=coarsen_plane(v_fine,fine,jp,kp)/2
        + (coarsen_plane(v_fine,fine+ip,jp,kp)
           + coarsen_plane(v_fine,fine-ip,jp,kp))/4;
                                   
      if(have_faults() && !is_residual)
        {
          SAMRAI::pdat::CellIndex cell(fine);
          result+=coarsen_plane_correction(*dv_mixed,fine,jp,kp)
            + ((*dv_diagonal)(cell-ip,axis) - (*dv_diagonal)(cell,axis)
               + (*dv_diagonal)(cell-ip+jp,axis) - (*dv_diagonal)(cell+jp,axis)
               + (*dv_diagonal)(cell-ip+kp,axis) - (*dv_diagonal)(cell+kp,axis)
               + (*dv_diagonal)(cell-ip+jp+kp,axis)
               - (*dv_diagonal)(cell+jp+kp,axis))/16;
        }
      return result;
    }


    //@}

    //@{
    /*!
     * @brief Set the data id that should be filled when setting
     * physical boundary conditions.
     *
     * When setPhysicalBoundaryConditions is called, the data
     * specified will be set.  This information is required because
     * the it is not passed in through the argument list of
     * setPhysicalBounaryConditions.
     */

    int data_id;
    bool is_residual;

    void set_extra_ids(const int& Dv_diagonal_id, const int& Dv_mixed_id,
                       const int& Level_set_id)
    {
      dv_diagonal_id=Dv_diagonal_id;
      dv_mixed_id=Dv_mixed_id;
      level_set_id=Level_set_id;
    }

    bool have_faults() const
    {
      return dv_diagonal_id!=invalid_id;
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }

    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::hier::CoarseFineBoundary> >
    coarse_fine;

  private:

    std::string d_object_name;

    /*!
     * @brief SAMRAI::hier::Index of target patch data when filling ghosts.
     */
    int dv_diagonal_id, dv_mixed_id, level_set_id;

    Boundary_Conditions& d_boundary_conditions;

    /*!
     * @brief Timers for performance measurement.
     */
    // boost::shared_ptr<tbox::Timer> t_set_boundary_values_in_cells;
    // boost::shared_ptr<tbox::Timer> t_use_set_bc_coefs;
  };

}
#endif  // included_solv_V_Coarsen_Patch_Strategy
