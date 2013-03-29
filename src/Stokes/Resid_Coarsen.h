#ifndef GAMRA_STOKES_STOKES_RESID_COARSEN_H
#define GAMRA_STOKES_STOKES_RESID_COARSEN_H

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include <string>

namespace Stokes {
  /**
   * Coarsens using the viscosities as weights.  So in 2D
   resid_coarse = (resid(i,j)*viscosity(i,j)
   + resid(i,j+1)*viscosity(i,j+1)
   + resid(i+1,j)*viscosity(i+1,j)
   + resid(i+1,j+1)*viscosity(i+1,j+1))/(4*viscosity_coarse)
   * @see SAMRAI::hier::CoarsenOperator
   */

  class Resid_Coarsen:
    public SAMRAI::hier::CoarsenOperator
  {
  public:
    explicit Resid_Coarsen(const SAMRAI::tbox::Dimension& dim,
                           const int &cell_viscosity):
      SAMRAI::hier::CoarsenOperator(dim, "RESID_COARSEN"),
      cell_viscosity_id(cell_viscosity)
    {
      d_name_id = "RESID_COARSEN";
    }

    virtual ~Resid_Coarsen(){}

    bool findCoarsenOperator(const boost::shared_ptr<SAMRAI::hier::Variable>& var,
                             const std::string& op_name) const
    {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

      const boost::shared_ptr<SAMRAI::pdat::CellVariable<double> >
        cast_var(boost::dynamic_pointer_cast<SAMRAI::pdat::CellVariable<double> >
                 (var));
      if (cast_var && (op_name == d_name_id)) {
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
    const int cell_viscosity_id;
  };

}
#endif
