#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/pdat/SideVariable.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include "Dir.hxx"

namespace Stokes
{
  class V_Boundary_Refine: public SAMRAI::hier::RefineOperator
  {
  public:
    explicit V_Boundary_Refine():
      SAMRAI::hier::RefineOperator("V_BOUNDARY_REFINE") { }
    virtual ~V_Boundary_Refine(){}
    virtual int getOperatorPriority() const
    {
      return 0;
    }

    virtual SAMRAI::hier::IntVector getStencilWidth
    (const SAMRAI::tbox::Dimension& dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }

    void refine(SAMRAI::hier::Patch& fine,
                const SAMRAI::hier::Patch& coarse,
                const int dst_component,
                const int src_component,
                const SAMRAI::hier::BoxOverlap& fine_overlap,
                const SAMRAI::hier::IntVector& ratio) const;

    void refine_box(SAMRAI::hier::Patch& fine,
                    const SAMRAI::hier::Patch& coarse,
                    const int dst_component,
                    const int src_component,
                    const SAMRAI::hier::Box& fine_box,
                    const SAMRAI::hier::IntVector& ratio,
                    const Gamra::Dir &axis) const;

  private:
    void Update_V_2D
    (const int &axis,
     const int &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector &ip, const SAMRAI::hier::IntVector &jp,
     int &i, int &j,
     const int &i_max,
     const int &j_min,
     const int &j_max,
     SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_fine) const;

    void Update_V_3D
    (const Gamra::Dir &axis,
     const Gamra::Dir &boundary_direction,
     const bool &boundary_positive,
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::hier::IntVector pp[],
     const SAMRAI::hier::Index &ijk,
     const SAMRAI::hier::Index &p_min, const SAMRAI::hier::Index &p_max,
     SAMRAI::pdat::SideData<double> &v,
     SAMRAI::pdat::SideData<double> &v_fine) const;
  };
}

