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
      d_boundary_conditions(bc) {}

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
                       const SAMRAI::hier::IntVector& ratio)
    {
      if(getDim().getValue()==2)
        postprocessCoarsen_2D(coarse,fine,coarse_box,ratio);
      else
        postprocessCoarsen_3D(coarse,fine,coarse_box,ratio);
    }

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
    (boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& dv_mixed,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> >& dv_diagonal,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry>& coarse_geom,
     const SAMRAI::hier::Box& coarse_box) const;

    void
    fix_boundary_elements_2D
    (boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& dv_mixed,
     const SAMRAI::hier::Box& coarse_box,
     const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox> &boundaries) const;

    void
    coarsen_3D
    (boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& dv_mixed,
     const boost::shared_ptr<SAMRAI::pdat::CellData<double> >& dv_diagonal,
     const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry>& coarse_geom,
     const SAMRAI::hier::Box& coarse_box) const;

    void
    fix_boundary_elements_3D
    (boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& v_fine,
     const boost::shared_ptr<SAMRAI::pdat::SideData<double> >& dv_mixed,
     const SAMRAI::hier::Box& coarse_box,
     const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox>& boundaries) const;

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

    void set_extra_ids(const int& cell_moduli, const int& edge_moduli,
                       const int& dv_diagonal, const int& dv_mixed)
    {
      cell_moduli_id=cell_moduli;
      edge_moduli_id=edge_moduli;
      dv_diagonal_id=dv_diagonal;
      dv_mixed_id=dv_mixed;
    }

    SAMRAI::tbox::Array<boost::shared_ptr<SAMRAI::hier::CoarseFineBoundary> >
    coarse_fine;

  private:

    std::string d_object_name;

    /*!
     * @brief SAMRAI::hier::Index of target patch data when filling ghosts.
     */
    int cell_moduli_id, edge_moduli_id, dv_diagonal_id, dv_mixed_id;

    Boundary_Conditions& d_boundary_conditions;

    /*!
     * @brief Timers for performance measurement.
     */
    // boost::shared_ptr<tbox::Timer> t_set_boundary_values_in_cells;
    // boost::shared_ptr<tbox::Timer> t_use_set_bc_coefs;
  };

}
#endif  // included_solv_V_Coarsen_Patch_Strategy
