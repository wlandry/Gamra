#ifndef GAMRA_ELASTIC_Resid_Coarsen_H
#define GAMRA_ELASTIC_Resid_Coarsen_H

#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include <string>

namespace Elastic {
  /**
   * Coarsens using the moduli as weights.  So in 2D
   resid_coarse = (resid(i,j)*moduli(i,j)
   + resid(i,j+1)*moduli(i,j+1)
   + resid(i+1,j)*moduli(i+1,j)
   + resid(i+1,j+1)*moduli(i+1,j+1))/(4*moduli_coarse)
   * @see xfer::CoarsenOperator
   */

  class Resid_Coarsen:
    public SAMRAI::xfer::CoarsenOperator
  {
  public:
    explicit Resid_Coarsen(const SAMRAI::tbox::Dimension& dim,
                           const int &cell_moduli):
      SAMRAI::xfer::CoarsenOperator(dim, "RESID_COARSEN"),
      cell_moduli_id(cell_moduli)
    {
      d_name_id = "RESID_COARSEN";
    }

    virtual ~Resid_Coarsen(){}

    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable>& var,
                             const std::string& op_name) const
    {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

      const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<double> >
        cast_var(var);
      if (!cast_var.isNull() && (op_name == d_name_id)) {
        return true;
      } else {
        return false;
      }
    }

    const std::string& getOperatorName() const
    {
      return d_name_id;
    }

    int getOperatorPriority() const
    {
      return 0;
    }

    SAMRAI::hier::IntVector getStencilWidth() const
    {
      return SAMRAI::hier::IntVector::getZero(getDim());
    }

    void coarsen(SAMRAI::hier::Patch& coarse,
                 const SAMRAI::hier::Patch& fine,
                 const int dst_component,
                 const int src_component,
                 const SAMRAI::hier::Box& coarse_box,
                 const SAMRAI::hier::IntVector& ratio) const;

  private:
    std::string d_name_id;
    const int cell_moduli_id;
  };

}
#endif
