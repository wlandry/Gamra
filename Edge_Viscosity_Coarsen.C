/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for side-centered double data on
 *                a Cartesian mesh. 
 *
 ************************************************************************/

#ifndef included_geom_Edge_Viscosity_Coarsen_C
#define included_geom_Edge_Viscosity_Coarsen_C

#include "Edge_Viscosity_Coarsen.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "viscosity_coarsen.h"

using namespace SAMRAI;

void SAMRAI::geom::Edge_Viscosity_Coarsen::coarsen(hier::Patch& coarse,
                                      const hier::Patch& fine,
                                      const int dst_component,
                                      const int src_component,
                                      const hier::Box& coarse_box,
                                      const hier::IntVector& ratio) const
{
  const tbox::Dimension& dimension(getDim());
  const int dim(dimension.getValue());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, coarse, fine, coarse_box, ratio);

  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_fine_ptr = fine.getPatchData(cell_viscosity_id);
  pdat::CellData<double> &cell_viscosity_fine(*cell_viscosity_fine_ptr);

  tbox::Pointer<pdat::NodeData<double> > edge_viscosity_fine_2D_ptr;
  tbox::Pointer<pdat::EdgeData<double> > edge_viscosity_fine_3D_ptr;
  if(dim==2)
    {
      edge_viscosity_fine_2D_ptr = fine.getPatchData(src_component);
      TBOX_ASSERT(!edge_viscosity_fine_2D_ptr.isNull());
    }
  else
    {
      edge_viscosity_fine_3D_ptr = fine.getPatchData(src_component);
      TBOX_ASSERT(!edge_viscosity_fine_3D_ptr.isNull());
    }

  tbox::Pointer<pdat::NodeData<double> > edge_viscosity_coarse_2D_ptr;
  tbox::Pointer<pdat::EdgeData<double> > edge_viscosity_coarse_3D_ptr;
  if(dim==2)
    {
      edge_viscosity_coarse_2D_ptr = coarse.getPatchData(dst_component);
      TBOX_ASSERT(!edge_viscosity_coarse_2D_ptr.isNull());
      TBOX_ASSERT(edge_viscosity_fine_2D_ptr->getDepth()
                  == edge_viscosity_coarse_2D_ptr->getDepth());
      TBOX_ASSERT(edge_viscosity_coarse_2D_ptr->getDepth() == 1);
    }
  else
    {
      edge_viscosity_coarse_3D_ptr = coarse.getPatchData(dst_component);
      TBOX_ASSERT(!edge_viscosity_coarse_3D_ptr.isNull());
      TBOX_ASSERT(edge_viscosity_fine_3D_ptr->getDepth()
                  == edge_viscosity_coarse_3D_ptr->getDepth());
      TBOX_ASSERT(edge_viscosity_coarse_3D_ptr->getDepth() == 1);
    }


  const tbox::Pointer<CartesianPatchGeometry> cgeom =
    coarse.getPatchGeometry();

  hier::Index ip(hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(dim>2)
    kp[2]=1;
  hier::Index pp[]={ip,jp,kp};

  if(dim==2)
    {
      for(pdat::NodeIterator ni(coarse_box); ni; ni++)
        {
          pdat::NodeIndex coarse_edge(*ni);

          (*edge_viscosity_coarse_2D_ptr)(coarse_edge)=
            viscosity_coarsen_2D(cell_viscosity_fine,*edge_viscosity_fine_2D_ptr,
                                 coarse_edge*2);
        }
    }
  else
    {
      for(int axis=0;axis<3;++axis)
        {
          const int axis2((axis+1)%3), axis3((axis+2)%3);
          for(pdat::EdgeIterator ni(coarse_box,axis); ni; ni++)
            {
              pdat::EdgeIndex coarse_edge(*ni);

              (*edge_viscosity_coarse_3D_ptr)(coarse_edge)=
                viscosity_coarsen_3D(cell_viscosity_fine,
                                     *edge_viscosity_fine_3D_ptr,
                                     coarse_edge*2-pp[axis2]-pp[axis3]);
            }
        }
    }
}
#endif
