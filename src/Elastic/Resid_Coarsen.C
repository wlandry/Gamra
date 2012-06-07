#include "Resid_Coarsen.h"

/**
 * Coarsens using the moduli as weights.  So in 2D
   resid_coarse = (resid(i,j)*moduli(i,j)
                   + resid(i,j+1)*moduli(i,j+1)
                   + resid(i+1,j)*moduli(i+1,j)
                   + resid(i+1,j+1)*moduli(i+1,j+1))/(4*moduli_coarse)
 */

void SAMRAI::geom::Elastic::Resid_Coarsen::coarsen
(hier::Patch& coarse,
 const hier::Patch& fine,
 const int dst_component,
 const int src_component,
 const hier::Box& coarse_box,
 const hier::IntVector& ratio) const
{
  const tbox::Dimension& dimension(getDim());
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, coarse, fine, coarse_box, ratio);
  
  tbox::Pointer<pdat::CellData<double> >
    r_fine_ptr = fine.getPatchData(src_component);
  pdat::CellData<double> &r_fine(*r_fine_ptr);
  tbox::Pointer<pdat::CellData<double> >
    r_ptr = coarse.getPatchData(dst_component);
  pdat::CellData<double> &r(*r_ptr);
  tbox::Pointer<pdat::CellData<double> >
    cell_moduli_fine_ptr = fine.getPatchData(cell_moduli_id);
  pdat::CellData<double> &cell_moduli_fine(*cell_moduli_fine_ptr);

  TBOX_ASSERT(!r_ptr.isNull());
  TBOX_ASSERT(!r_fine_ptr.isNull());
  TBOX_ASSERT(r_fine.getDepth() == r.getDepth());
  TBOX_ASSERT(r.getDepth() == 1);

  hier::Box cell_box(hier::Index::getZeroIndex(dimension),
                     hier::Index::getOneIndex(dimension));

  for(pdat::CellIterator ci(coarse.getBox()); ci; ci++)
    {
      pdat::CellIndex coarse(*ci);
      pdat::CellIndex fine(coarse*2);
      double temp(0), moduli_sum(0);

      for(pdat::CellIterator ii(cell_box); ii; ii++)
        {
          pdat::CellIndex i(*ii);
          temp+=r_fine(fine+i)*(cell_moduli_fine(fine+i,0)+cell_moduli_fine(fine+i,1));
          moduli_sum+=(cell_moduli_fine(fine+i,0)+cell_moduli_fine(fine+i,1));
        }
      r(coarse)=temp/moduli_sum;
    }
}
