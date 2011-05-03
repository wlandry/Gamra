#include "Resid_Coarsen.h"

/**
 * Coarsens using the viscosities as weights.  So in 2D
   resid_coarse = (resid(i,j)*viscosity(i,j)
                   + resid(i,j+1)*viscosity(i,j+1)
                   + resid(i+1,j)*viscosity(i+1,j)
                   + resid(i+1,j+1)*viscosity(i+1,j+1))/(4*viscosity_coarse)
 */

void SAMRAI::geom::Resid_Coarsen::coarsen(hier::Patch& coarse,
                                          const hier::Patch& fine,
                                          const int dst_component,
                                          const int src_component,
                                          const hier::Box& coarse_box,
                                          const hier::IntVector& ratio) const
{
  const tbox::Dimension& dimension(getDim());
  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, coarse, fine, coarse_box, ratio);
  const int dim(dimension.getValue());
  
  tbox::Pointer<pdat::CellData<double> >
    r_fine_ptr = fine.getPatchData(src_component);
  pdat::CellData<double> &r_fine(*r_fine_ptr);
  tbox::Pointer<pdat::CellData<double> >
    r_ptr = coarse.getPatchData(dst_component);
  pdat::CellData<double> &r(*r_ptr);
  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_ptr = coarse.getPatchData(dst_component);
  pdat::CellData<double> &cell_viscosity(*cell_viscosity_ptr);
  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_fine_ptr = fine.getPatchData(dst_component);
  pdat::CellData<double> &cell_viscosity_fine(*cell_viscosity_fine_ptr);

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
      double temp(0);
      for(pdat::CellIterator ii(cell_box); ii; ii++)
        {
          pdat::CellIndex i(*ii);
          temp+=r_fine(fine+i)*cell_viscosity_fine(fine+i);
        }
      r(coarse)=temp/(4*(dim-1)*cell_viscosity(coarse));
    }
}
