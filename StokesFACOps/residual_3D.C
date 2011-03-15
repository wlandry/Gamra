#include "StokesFACOps.h"
#include "Boundary.h"

void SAMRAI::solv::StokesFACOps::residual_3D
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
  const hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);

  tbox::Pointer<pdat::EdgeData<double> >
    edge_viscosity_ptr = patch.getPatchData(edge_viscosity_id);
  pdat::EdgeData<double> &edge_viscosity(*edge_viscosity_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];
  double dz = geom.getDx()[2];

  for(pdat::CellIterator ci(pbox); ci; ci++)
    {
      pdat::CellIndex center(*ci);
      pdat::CellIndex up(center), down(center), right(center),
        left(center), front(center), back(center);

      ++right[0];
      --left[0];
      ++up[1];
      --down[1];
      ++front[2];
      --back[2];

      const pdat::SideIndex
        x(center,0,pdat::SideIndex::Lower),
        y(center,1,pdat::SideIndex::Lower),
        z(center,2,pdat::SideIndex::Lower);
      const pdat::EdgeIndex
        edge_x(center,0,pdat::EdgeIndex::LowerLeft),
        edge_y(center,1,pdat::EdgeIndex::LowerLeft),
        edge_z(center,2,pdat::EdgeIndex::LowerLeft);

      /* p */
      if(center[0]!=pbox.upper(0) && center[1]!=pbox.upper(1)
         && center[2]!=pbox.upper(2))
        {
          double dvx_dx=(v(x+ip) - v(x))/dx;
          double dvy_dy=(v(y+jp) - v(y))/dy;
          double dvz_dz=(v(z+kp) - v(z))/dy;
          p_resid(center)=p_rhs(center) - dvx_dx - dvy_dy - dvz_dz;
        }

      /* vx */
      if(center[1]!=pbox.upper(1) && center[2]!=pbox.upper(2))
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
                - v_operator_3D(v,p,cell_viscosity,edge_viscosity,
                                center,edge_y,edge_z,x,y,z,ip,jp,kp,dx,dy,dz);
            }
        }

      /* vy */
      if(center[0]!=pbox.upper(0) && center[2]!=pbox.upper(2))
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
                - v_operator_3D(v,p,cell_viscosity,edge_viscosity,
                                center,edge_z,edge_x,y,z,x,jp,kp,ip,dy,dz,dx);
            }
        }

      /* vz */
      if(center[1]!=pbox.upper(0) && center[2]!=pbox.upper(2))
        {
          /* If z==0 */
          if((center[2]==pbox.lower(1) && v(z-kp)==boundary_value)
             || (center[2]==pbox.upper(1) && v(z+kp)==boundary_value))
            {
              v_resid(z)=0;
            }
          else
            {
              v_resid(z)=v_rhs(z)
                - v_operator_3D(v,p,cell_viscosity,edge_viscosity,
                                center,edge_x,edge_y,z,x,y,kp,ip,jp,dz,dx,dy);
            }
        }
    }
}

