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
  class V_Refine:
    public SAMRAI::hier::RefineOperator
  {
  public:
    explicit V_Refine(): SAMRAI::hier::RefineOperator("V_REFINE") { }
    virtual ~V_Refine(){}
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
  };
}

