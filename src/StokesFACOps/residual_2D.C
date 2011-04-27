#include "StokesFACOps.h"
#include "Constants.h"

void SAMRAI::solv::StokesFACOps::residual_2D
(pdat::CellData<double> &p,
 pdat::SideData<double> &v,
 pdat::CellData<double> &cell_viscosity,
 pdat::CellData<double> &p_rhs,
 pdat::SideData<double> &v_rhs,
 pdat::CellData<double> &p_resid,
 pdat::SideData<double> &v_resid,
 hier::Patch &patch,
 const hier::Box &pbox,
 const geom::CartesianPatchGeometry &geom)
{
  const hier::Index ip(1,0), jp(0,1);

  tbox::Pointer<pdat::NodeData<double> >
    edge_viscosity_ptr = patch.getPatchData(edge_viscosity_id);
  pdat::NodeData<double> &edge_viscosity(*edge_viscosity_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];

  for(pdat::CellIterator ci(pbox); ci; ci++)
    {
      pdat::CellIndex center(*ci);
      pdat::CellIndex up(center), down(center), right(center),
        left(center);

      ++up[1];
      --down[1];
      ++right[0];
      --left[0];

      const pdat::SideIndex
        x(center,0,pdat::SideIndex::Lower),
        y(center,1,pdat::SideIndex::Lower);
      const pdat::NodeIndex
        edge(center,pdat::NodeIndex::LowerLeft);

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

