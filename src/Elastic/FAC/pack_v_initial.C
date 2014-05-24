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

  if(d_dim.getValue()==2)
    {
      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell) {

        const SAMRAI::pdat::CellIndex &center(*icell);
        const SAMRAI::pdat::SideIndex
          x(center,0,SAMRAI::pdat::SideIndex::Lower),
          y(center,1,SAMRAI::pdat::SideIndex::Lower);
	
        double coord[3];
        for(int d=0;d<d_dim.getValue();++d)
          coord[d]=geom->getXLower()[d] + dx[d]*(x[d]-pbox.lower()[d]+0.5);

        if (0==depth)
          {
            *buffer = v_initial[0].eval(coord);
          }
        else
          {
            *buffer = v_initial[1].eval(coord);
          }
        buffer = buffer + 1;
      }
    }
  else
    {
      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell) {

        const SAMRAI::pdat::CellIndex &center(*icell);
	
        double coord[3];
        for(int d=0;d<d_dim.getValue();++d)
          coord[d]=geom->getXLower()[d] + dx[d]*(center[d]-pbox.lower()[d]+0.5);

        if (0==depth)
          {
            *buffer = v_initial[0].eval(coord);
          }
        else if (1==depth)
          {
            *buffer = v_initial[1].eval(coord);
          }
        else
          {
            *buffer = v_initial[2].eval(coord);
          }
        buffer = buffer + 1;
      }
    }
}
