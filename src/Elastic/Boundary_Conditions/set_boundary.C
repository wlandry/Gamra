#include "Constants.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "Elastic/Boundary_Conditions.h"

void Elastic::Boundary_Conditions::set_boundary
(const SAMRAI::hier::Patch& patch, const int &v_id, const bool &homogeneous)
{
  try
    {
      /* We need to set dirichlet boundaries before traction
         boundaries, and tangential traction boudaries before normal
         traction boundaries. */

      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<double> >
        v_ptr = patch.getPatchData(v_id);
      SAMRAI::pdat::SideData<double> &v(*v_ptr);
      const SAMRAI::tbox::Dimension Dim(patch.getDim());
      const int dim(Dim.getValue());
      const SAMRAI::hier::Index zero(SAMRAI::hier::Index::getZeroIndex(Dim));

      SAMRAI::hier::Index pp[]={zero,zero,zero};
      for(int i=0;i<dim;++i)
        pp[i][i]=1;

      const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry>
        geom = patch.getPatchGeometry();
      const double *dx=geom->getDx();

      const SAMRAI::hier::Box pbox=patch.getBox();
      const SAMRAI::hier::Box gbox=v.getGhostBox();

      set_dirichlet(v,pp,dim,pbox,gbox,geom,dx,homogeneous);
      set_shear_derivs(v,pp,dim,pbox,gbox,geom,dx,homogeneous);

      if(dim==2)
        {
          SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<double> >
            edge_moduli_ptr= patch.getPatchData(edge_moduli_id);
          if(!edge_moduli_ptr.isNull())
            {
              SAMRAI::pdat::NodeData<double> &edge_moduli(*edge_moduli_ptr);
              set_normal_stress(v,edge_moduli,pp,dim,pbox,gbox,geom,dx,
                                homogeneous);
            }
        }
      else
        {
          SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<double> >
            edge_moduli_ptr= patch.getPatchData(edge_moduli_id);
          if(!edge_moduli_ptr.isNull())
            {
              SAMRAI::pdat::EdgeData<double> &edge_moduli(*edge_moduli_ptr);
              set_normal_stress(v,edge_moduli,pp,dim,pbox,gbox,geom,dx,
                                homogeneous);
            }
        }
    }
  catch(mu::Parser::exception_type &e)
    {
      TBOX_ERROR("Error in input formula\n"
                 << "Message:  " << e.GetMsg() << "\n"
                 << "Formula:  " << e.GetExpr() << "\n"
                 << "Token:    " << e.GetToken() << "\n"
                 << "Position: " << e.GetPos() << "\n"
                 << "Errc:     " << e.GetCode() << "\n");
    }
}
