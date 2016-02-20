#include "Elastic/FACOps.hxx"

double Elastic::FACOps::v_level_set_operator_2D
(const SAMRAI::pdat::SideData<double> &level_set,
 const SAMRAI::pdat::SideData<double> &v,
 const SAMRAI::pdat::CellData<double> &cell_moduli,
 const SAMRAI::pdat::NodeData<double> &edge_moduli,
 const SAMRAI::pdat::CellIndex &cell,
 const SAMRAI::pdat::NodeIndex &edge,
 const SAMRAI::pdat::SideIndex &x,
 const SAMRAI::pdat::SideIndex &y,
 const SAMRAI::hier::Index &ip,
 const SAMRAI::hier::Index &jp,
 const double &dx,
 const double &dy)
{
  double v_p(v(x+ip)), v_m(v(x-ip));
  if(level_set(x+ip)<0)
    {
      v_p=0;
    }
  if(level_set(x-ip)<0)
    {
      v_m=0;
    }
  double dx_xx=
    ((v_p-v(x))*(cell_moduli(cell,0)
                 + 2*cell_moduli(cell,1))
     - (v(x)-v_m)*(cell_moduli(cell-ip,0)
                   + 2*cell_moduli(cell-ip,1)))
    /(dx*dx);
              
  double dx_yp(v(x+jp)-v(x)),
    dx_ym(v(x)-v(x-jp));
  if(level_set(x+jp)<0)
    {
      dx_yp=0;
    }
  if(level_set(x-jp)<0)
    {
      dx_ym=0;
    }
  double dx_yy=
    (edge_moduli(edge+jp,1)*dx_yp
     - edge_moduli(edge,1)*dx_ym)/(dy*dy);

  double dvy_xp(v(y+jp)-v(y+jp-ip)),
    dvy_x(v(y)-v(y-ip));
  if(level_set(y+jp)<0 || level_set(y+jp-ip)<0)
    {
      dvy_xp=0;
    }
  if(level_set(y)<0 || level_set(y-ip)<0)
    {
      dvy_x=0;
    }
  double dy_xy=
    (edge_moduli(edge+jp,1)*dvy_xp
     - edge_moduli(edge,1)*dvy_x)/(dx*dy);

  /* FIXME: This needs to be
   * interpolated to the dirichlet
   * conditions */
  double vy_pp(level_set(y+jp)<0 ? 0 : v(y+jp)),
    vy_pm(level_set(y)<0 ? 0 : v(y)),
    vy_mp(level_set(y+jp-ip)<0 ? 0 : v(y+jp-ip)),
    vy_mm(level_set(y-ip)<0 ? 0 : v(y-ip));

  double dvy_yp(vy_pp-vy_pm), dvy_ym(vy_mp-vy_mm);
  double dy_yx=
    (cell_moduli(cell,0)*dvy_yp
     - cell_moduli(cell-ip,0)*dvy_ym)/(dx*dy);

  return dx_xx + dx_yy + dy_xy + dy_yx;
}
