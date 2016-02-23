#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>

#include "quad_offset_interpolate.hxx"
#include "Constants.hxx"

namespace Elastic {
  class Coarse_Fine_Boundary_Refine: public SAMRAI::hier::RefineOperator
  {
  public:

    static bool is_residual;
    static int dv_diagonal_id, dv_mixed_id, level_set_id;
    explicit Coarse_Fine_Boundary_Refine():
      SAMRAI::hier::RefineOperator("COARSE_FINE_BOUNDARY_REFINE") { }

    bool have_faults() const
    {
      return dv_diagonal_id!=invalid_id;
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }
    virtual ~Coarse_Fine_Boundary_Refine(){}

    virtual int getOperatorPriority() const
    {
      return 0;
    }

    virtual SAMRAI::hier::IntVector getStencilWidth
    (const SAMRAI::tbox::Dimension &dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }

    /// It is assumed that the coarse patch contains sufficient data
    /// for the stencil width of the refinement operator.

    virtual void refine(SAMRAI::hier::Patch &fine_patch,
                        const SAMRAI::hier::Patch &coarse_patch,
                        const int dst_component,
                        const int src_component,
                        const SAMRAI::hier::BoxOverlap &fine_overlap,
                        const SAMRAI::hier::IntVector &ratio) const;

    void refine_box(SAMRAI::hier::Patch &fine_patch,
                    const SAMRAI::hier::Patch &coarse_patch,
                    const int dst_component,
                    const int src_component,
                    const SAMRAI::hier::Box &fine_box,
                    const SAMRAI::hier::IntVector &ratio,
                    const Gamra::Dir &axis) const;

  private:
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

