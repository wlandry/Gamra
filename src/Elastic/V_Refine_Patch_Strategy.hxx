#pragma once

/// Copyright © 1997-2010 Lawrence Livermore National Security, LLC
/// Copyright © 2013-2016 California Institute of Technology
/// Copyright © 2013-2016 Nanyang Technical University

#include <SAMRAI/xfer/RefinePatchStrategy.h>
#include "Boundary_Conditions.hxx"

namespace Elastic {
  class V_Refine_Patch_Strategy: public SAMRAI::xfer::RefinePatchStrategy
  {
  public:
    V_Refine_Patch_Strategy(std::string object_name,
                            const Boundary_Conditions &bc):
      SAMRAI::xfer::RefinePatchStrategy(), d_object_name(object_name),
      boundary_conditions(bc) {}

    virtual ~V_Refine_Patch_Strategy() {}

    virtual void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch,
                                               const double ,
                                               const SAMRAI::hier::IntVector& )
    {
      /// Do not apply normal stress bc's, since those external points
      /// are not used anyway.
      boundary_conditions.set_physical_boundary(patch,data_id,is_residual,
                                                false);
    }
    SAMRAI::hier::IntVector getRefineOpStencilWidth
    (const SAMRAI::tbox::Dimension& dim) const
    {
      return SAMRAI::hier::IntVector::getOne(dim);
    }
    virtual void preprocessRefine(SAMRAI::hier::Patch& ,
                                  const SAMRAI::hier::Patch& coarse,
                                  const SAMRAI::hier::Box& ,
                                  const SAMRAI::hier::IntVector& )
    {
      /// Do not apply normal stress bc's, since those external points
      /// are not used anyway.
      boundary_conditions.set_physical_boundary(coarse,data_id,is_residual,
                                                false);
    }

    virtual void postprocessRefineBoxes(SAMRAI::hier::Patch& ,
                                        const SAMRAI::hier::Patch& ,
                                        const SAMRAI::hier::BoxContainer& ,
                                        const SAMRAI::hier::IntVector& ) {}
    virtual void postprocessRefine(SAMRAI::hier::Patch& ,
                                   const SAMRAI::hier::Patch& ,
                                   const SAMRAI::hier::Box& ,
                                   const SAMRAI::hier::IntVector& ) {}
    int data_id;
    bool is_residual;
  private:
    std::string d_object_name;
    const Boundary_Conditions &boundary_conditions;
  };

}

