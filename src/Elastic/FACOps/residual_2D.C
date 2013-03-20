#include "Elastic/FACOps.h"
#include "Constants.h"

void Elastic::FACOps::residual_2D
(SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::CellData<double> &cell_moduli,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::SideData<double> &v_resid,
 SAMRAI::hier::Patch &patch,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::geom::CartesianPatchGeometry &geom)
{
  const SAMRAI::hier::Index unit[]={SAMRAI::hier::Index(1,0),
                                    SAMRAI::hier::Index(0,1)};

  boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_moduli_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
    (patch.getPatchData(edge_moduli_id));
  SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];

  SAMRAI::pdat::CellIterator cend(pbox,false);
  for(SAMRAI::pdat::CellIterator ci(pbox,true); ci!=cend; ci++)
    {
      SAMRAI::pdat::CellIndex center(*ci);

      const SAMRAI::pdat::NodeIndex
        edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

      const int dim(2);
      for(int ix=0;ix<dim;++ix)
        {
          const int iy((ix+1)%dim);
          const SAMRAI::pdat::SideIndex
            x(center,ix,SAMRAI::pdat::SideIndex::Lower),
            y(center,iy,SAMRAI::pdat::SideIndex::Lower);
          const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);

          if(center[iy]!=pbox.upper(iy))
            {
              /* If x==0 */
              if((center[ix]==pbox.lower(ix) && v(x-ip)==boundary_value)
                 || (center[ix]==pbox.upper(ix) && v(x+ip)==boundary_value))
                {
                  v_resid(x)=0;
                }
              else
                {
                  v_resid(x)=v_rhs(x)
                    - v_operator_2D(v,cell_moduli,edge_moduli,center,
                                    edge,x,y,ip,jp,dx,dy);
                }
            }
        }
    }
}

