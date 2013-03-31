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

  if(have_embedded_boundary())
    {
      boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr =
        boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
        (patch.getPatchData(level_set_id));
      SAMRAI::pdat::SideData<double> &level_set(*level_set_ptr);
      for(SAMRAI::pdat::CellIterator ci(pbox,true); ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex &cell(*ci);

          const SAMRAI::pdat::NodeIndex
            edge(cell,SAMRAI::pdat::NodeIndex::LowerLeft);

          const int dim(2);
          for(int ix=0;ix<dim;++ix)
            {
              const int iy((ix+1)%dim);
              const SAMRAI::pdat::SideIndex
                x(cell,ix,SAMRAI::pdat::SideIndex::Lower),
                y(cell,iy,SAMRAI::pdat::SideIndex::Lower);
              const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);

              if(cell[iy]!=pbox.upper(iy))
                {
                  if(level_set(x)<0)
                    {
                      v_resid(x)=0;
                    }
                  else if(level_set(x)>1)
                    {
                      v_resid(x)=v_rhs(x)
                        - v_operator_2D(v,cell_moduli,edge_moduli,cell,
                                        edge,x,y,ip,jp,dx,dy);
                    }
                  else
                    {
                      v_resid(x)=v_rhs(x)
                        - v_level_set_operator_2D(level_set,v,cell_moduli,
                                                  edge_moduli,cell,
                                                  edge,x,y,ip,jp,dx,dy);
                    }
                }
            }
        }
    }
  else
    {
      for(SAMRAI::pdat::CellIterator ci(pbox,true); ci!=cend; ++ci)
        {
          const SAMRAI::pdat::CellIndex &cell(*ci);

          const SAMRAI::pdat::NodeIndex
            edge(cell,SAMRAI::pdat::NodeIndex::LowerLeft);

          const int dim(2);
          for(int ix=0;ix<dim;++ix)
            {
              const int iy((ix+1)%dim);
              const SAMRAI::pdat::SideIndex
                x(cell,ix,SAMRAI::pdat::SideIndex::Lower),
                y(cell,iy,SAMRAI::pdat::SideIndex::Lower);
              const SAMRAI::hier::Index ip(unit[ix]), jp(unit[iy]);

              if(cell[iy]!=pbox.upper(iy))
                {
                  /* If x==0 */
                  if((cell[ix]==pbox.lower(ix) && v(x-ip)==boundary_value)
                     || (cell[ix]==pbox.upper(ix) && v(x+ip)==boundary_value))
                    {
                      v_resid(x)=0;
                    }
                  else
                    {
                      v_resid(x)=v_rhs(x)
                        - v_operator_2D(v,cell_moduli,edge_moduli,cell,
                                        edge,x,y,ip,jp,dx,dy);
                    }
                }
            }
        }
    }
}
