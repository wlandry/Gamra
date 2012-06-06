#ifndef included_geom_Resid_Coarsen
#define included_geom_Resid_Coarsen

#include "SAMRAI/xfer/CoarsenOperator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Coarsens using the moduli as weights.  So in 2D
   resid_coarse = (resid(i,j)*moduli(i,j)
                   + resid(i,j+1)*moduli(i,j+1)
                   + resid(i+1,j)*moduli(i+1,j)
                   + resid(i+1,j+1)*moduli(i+1,j+1))/(4*moduli_coarse)
 * @see xfer::CoarsenOperator
 */

class Resid_Coarsen:
   public xfer::CoarsenOperator
{
public:
  explicit Resid_Coarsen(const tbox::Dimension& dim,
                         const int &cell_moduli):
    xfer::CoarsenOperator(dim, "RESID_COARSEN"),
    cell_moduli_id(cell_moduli)
  {
    d_name_id = "RESID_COARSEN";
  }

  virtual ~Resid_Coarsen(){}

  bool findCoarsenOperator(const tbox::Pointer<hier::Variable>& var,
                           const std::string& op_name) const
  {
    TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

    const tbox::Pointer<pdat::CellVariable<double> > cast_var(var);
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

  hier::IntVector getStencilWidth() const
  {
    return hier::IntVector::getZero(getDim());
  }

  void coarsen(hier::Patch& coarse,
               const hier::Patch& fine,
               const int dst_component,
               const int src_component,
               const hier::Box& coarse_box,
               const hier::IntVector& ratio) const;

private:
  std::string d_name_id;
  const int cell_moduli_id;
};

}
}
#endif
