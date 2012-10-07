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
  const SAMRAI::hier::Index ip(1,0), jp(0,1);

  SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<double> >
    edge_moduli_ptr = patch.getPatchData(edge_moduli_id);
  SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);

  double dx = geom.getDx()[0];
  double dy = geom.getDx()[1];

  for(SAMRAI::pdat::CellIterator ci(pbox); ci; ci++)
    {
      SAMRAI::pdat::CellIndex center(*ci);

      const SAMRAI::pdat::SideIndex
        x(center,0,SAMRAI::pdat::SideIndex::Lower),
        y(center,1,SAMRAI::pdat::SideIndex::Lower);
      const SAMRAI::pdat::NodeIndex
        edge(center,SAMRAI::pdat::NodeIndex::LowerLeft);

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
                - v_operator_2D(v,cell_moduli,edge_moduli,center,
                                edge,x,y,ip,jp,dx,dy);

          SAMRAI::tbox::plog << "vx resid "
                             << x << " "
                             << v_resid(x) << " "
                             << "\n";
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
                - v_operator_2D(v,cell_moduli,edge_moduli,center,
                                edge,y,x,jp,ip,dy,dx);
            }
        }
    }
}

