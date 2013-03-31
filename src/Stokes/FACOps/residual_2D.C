#include "Stokes/FACOps.h"
#include "Constants.h"

void Stokes::FACOps::residual_2D
(SAMRAI::pdat::CellData<double> &p,
 SAMRAI::pdat::SideData<double> &v,
 SAMRAI::pdat::CellData<double> &cell_viscosity,
 SAMRAI::pdat::CellData<double> &p_rhs,
 SAMRAI::pdat::SideData<double> &v_rhs,
 SAMRAI::pdat::CellData<double> &p_resid,
 SAMRAI::pdat::SideData<double> &v_resid,
 SAMRAI::hier::Patch &patch,
 const SAMRAI::hier::Box &pbox,
 const SAMRAI::geom::CartesianPatchGeometry &geom)
{
  const SAMRAI::hier::Index ip(1,0), jp(0,1);

  boost::shared_ptr<SAMRAI::pdat::NodeData<double> > edge_viscosity_ptr = 
    boost::dynamic_pointer_cast<SAMRAI::pdat::NodeData<double> >
    (patch.getPatchData(edge_viscosity_id));
  SAMRAI::pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];

  SAMRAI::pdat::CellIterator cend(pbox,false);
  for(SAMRAI::pdat::CellIterator ci(pbox,true); ci!=cend; ++ci)
    {
      SAMRAI::pdat::CellIndex center(*ci);
      SAMRAI::pdat::CellIndex up(center), down(center), right(center),
        left(center);

      ++up[1];
      --down[1];
      ++right[0];
      --left[0];

      const SAMRAI::pdat::SideIndex
        x(center,0,SAMRAI::pdat::SideIndex::Lower),
        y(center,1,SAMRAI::pdat::SideIndex::Lower);
      const SAMRAI::pdat::NodeIndex
        edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

      /* p */
      if(center[0]!=pbox.upper(0) && center[1]!=pbox.upper(1))
        {
          double dvx_dx=(v(x+ip) - v(x))/dx;
          double dvy_dy=(v(y+jp) - v(y))/dy;
          p_resid(center)=p_rhs(center) - dvx_dx - dvy_dy;
        }

      /* vx */
      if(center[1]!=pbox.upper(1))
        {
          /* If x==0 */
          if((center[0]==pbox.lower(0) && v(x-ip)==boundary_value)
             || (center[0]==pbox.upper(0) && v(x+ip)==boundary_value))
            {
              v_resid(x)=0;
            }
          else
            {
              v_resid(x)=v_rhs(x)
                - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                                edge,x,y,ip,jp,dx,dy);
            }
        }

      /* vy */
      if(center[0]!=pbox.upper(0))
        {
          /* If y==0 */
          if((center[1]==pbox.lower(1) && v(y-jp)==boundary_value)
             || (center[1]==pbox.upper(1) && v(y+jp)==boundary_value))
            {
              v_resid(y)=0;
            }
          else
            {
              v_resid(y)=v_rhs(y)
                - v_operator_2D(v,p,cell_viscosity,edge_viscosity,center,
                                edge,y,x,jp,ip,dy,dx);
            }
        }
    }
}

