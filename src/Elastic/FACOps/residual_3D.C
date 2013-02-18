#include "Elastic/FACOps.h"
#include "Constants.h"

void Elastic::FACOps::residual_3D
(SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::CellData<double> &cell_moduli,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::SideData<double> &v_resid,
 SAMRAI::hier::Patch &patch,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::geom::CartesianPatchGeometry &geom)
{
  boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > edge_moduli_ptr =
    boost::dynamic_pointer_cast<SAMRAI::pdat::EdgeData<double> >
    (patch.getPatchData(edge_moduli_id));
  SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);

  const double *Dx = geom.getDx();
  const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  const SAMRAI::hier::Index unit[]={ip,jp,kp};

  SAMRAI::pdat::CellIterator cend(pbox,false);
  for(SAMRAI::pdat::CellIterator ci(pbox,true); ci!=cend; ci++)
    {
      SAMRAI::pdat::CellIndex center(*ci);
      SAMRAI::pdat::CellIndex up(center), down(center), right(center),
        left(center), front(center), back(center);

      ++right[0];
      --left[0];
      ++up[1];
      --down[1];
      ++front[2];
      --back[2];

      for(int ix=0;ix<3;++ix)
        {
          const int iy((ix+1)%3), iz((ix+2)%3);
          const SAMRAI::pdat::SideIndex
            x(center,ix,SAMRAI::pdat::SideIndex::Lower),
            y(center,iy,SAMRAI::pdat::SideIndex::Lower),
            z(center,iz,SAMRAI::pdat::SideIndex::Lower);
          const SAMRAI::pdat::EdgeIndex
            edge_y(center,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
            edge_z(center,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);

          if(center[iy]!=pbox.upper(iy) && center[iz]!=pbox.upper(iz))
            {
              if((center[ix]==pbox.lower(ix) && v(x-unit[ix])==boundary_value)
                 || (center[ix]==pbox.upper(ix)
                     && v(x+unit[ix])==boundary_value))
                {
                  v_resid(x)=0;
                }
              else
                {
                  v_resid(x)=v_rhs(x)
                    - v_operator_3D(v,cell_moduli,edge_moduli,
                                    center,edge_y,edge_z,x,y,z,
                                    unit[ix],unit[iy],unit[iz],
                                    Dx[ix],Dx[iy],Dx[iz]);
                }
            }
        }          
    }
}

