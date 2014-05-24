#include "Elastic/FAC.h"

void
Elastic::FAC::pack_v_initial(double* buffer,
                             const SAMRAI::hier::Patch& patch,
                             const SAMRAI::hier::Box& region,
                             const int &depth) const
{
  const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry>
    &geom(boost::dynamic_pointer_cast
          <SAMRAI::geom::CartesianPatchGeometry>
          (patch.getPatchGeometry()));
  const SAMRAI::hier::Box &pbox(patch.getBox());
  const double *dx=geom->getDx();

  double offset[]={0.5,0.5,0.5};
  if(offset_vector_on_output)
    offset[depth]=0;

  SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
  for(SAMRAI::pdat::CellIterator
        icell(SAMRAI::pdat::CellGeometry::begin(region));
      icell!=iend; ++icell)
    {
      const SAMRAI::pdat::CellIndex &center(*icell);
	
      double coord[3];
      for(int d=0;d<d_dim.getValue();++d)
        coord[d]=geom->getXLower()[d]
          + dx[d]*(center[d]-pbox.lower()[d]+offset[d]);

      *buffer = (v_initial[depth].is_valid ? v_initial[depth].eval(coord) : 0);
      ++buffer;
    }
}
