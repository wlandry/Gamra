#include "Stokes/Resid_Coarsen.h"

/**
 * Coarsens using the viscosities as weights.  So in 2D
   resid_coarse = (resid(i,j)*viscosity(i,j)
                   + resid(i,j+1)*viscosity(i,j+1)
                   + resid(i+1,j)*viscosity(i+1,j)
                   + resid(i+1,j+1)*viscosity(i+1,j+1))/(4*viscosity_coarse)
 */

void SAMRAI::geom::Stokes::Resid_Coarsen::coarsen(hier::Patch& coarse,
                                                  const hier::Patch& fine,
                                                  const int dst_component,
                                                  const int src_component,
                                                  const hier::Box&,
                                                  const hier::IntVector&) const
{
  const tbox::Dimension& dimension(getDim());
  
  boost::shared_ptr<pdat::CellData<double> > r_fine_ptr =
    boost::dynamic_pointer_cast<pdat::CellData<double> >
    (fine.getPatchData(src_component));
  pdat::CellData<double> &r_fine(*r_fine_ptr);
  boost::shared_ptr<pdat::CellData<double> > r_ptr =
    boost::dynamic_pointer_cast<pdat::CellData<double> >
    (coarse.getPatchData(dst_component));
  pdat::CellData<double> &r(*r_ptr);
  boost::shared_ptr<pdat::CellData<double> > cell_viscosity_fine_ptr = 
    boost::dynamic_pointer_cast<pdat::CellData<double> >
    (fine.getPatchData(cell_viscosity_id));
  pdat::CellData<double> &cell_viscosity_fine(*cell_viscosity_fine_ptr);

  TBOX_ASSERT(r_ptr);
  TBOX_ASSERT(r_fine_ptr);
  TBOX_ASSERT(r_fine.getDepth() == r.getDepth());
  TBOX_ASSERT(r.getDepth() == 1);

  pdat::CellIterator cend(coarse.getBox(),false);
  for(pdat::CellIterator ci(coarse.getBox(),true); ci!=cend; ci++)
    {
      pdat::CellIndex coarse(*ci);
      pdat::CellIndex fine(coarse*2);
      double temp(0), viscosity_sum(0);

      const int dim(dimension.getValue());
      for(int i=0;i<(2 << dim); ++i)
        {
          hier::Index j(i%2, (i/2)%2, i/4);
          temp+=r_fine(fine+j)*cell_viscosity_fine(fine+j);
          viscosity_sum+=cell_viscosity_fine(fine+j);
        }
      r(coarse)=temp/viscosity_sum;
    }
}
