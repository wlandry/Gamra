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

#include <FTensor.hpp>
#include <string>

#include "Constants.hxx"

namespace Elastic {
  class V_Refine: public SAMRAI::hier::RefineOperator
  {
  public:
    explicit V_Refine(): SAMRAI::hier::RefineOperator("V_REFINE") { }
    virtual ~V_Refine(){}
    
    virtual int getOperatorPriority() const { return 0; }

    virtual SAMRAI::hier::IntVector getStencilWidth
    (const SAMRAI::tbox::Dimension &dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }

    bool have_embedded_boundary() const
    {
      return level_set_id!=invalid_id;
    }

    virtual void refine(SAMRAI::hier::Patch& fine,
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

    double refine_along_line
    (SAMRAI::pdat::SideData<double> &v,
     const Gamra::Dir &axis,
     const Gamra::Dir &dim,
     const SAMRAI::hier::Index pp[],
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::pdat::SideIndex &coarse,
     const SAMRAI::hier::Box &coarse_box,
     const SAMRAI::geom::CartesianPatchGeometry &coarse_geom) const;

    double refine_along_line
    (SAMRAI::pdat::SideData<double> &v,
     const Gamra::Dir &ix,
     const Gamra::Dir &dim,
     const SAMRAI::hier::Index unit[],
     const SAMRAI::pdat::SideIndex &fine,
     const SAMRAI::pdat::SideIndex &coarse,
     const SAMRAI::pdat::SideData<double> &level_set_coarse) const;

    // FIXME: Get rid of this global
    static int level_set_id;
  };

}

