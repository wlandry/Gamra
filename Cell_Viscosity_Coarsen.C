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

#ifndef included_geom_Cell_Viscosity_Coarsen_C
#define included_geom_Cell_Viscosity_Coarsen_C

#include "Cell_Viscosity_Coarsen.h"

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

void SAMRAI::geom::Cell_Viscosity_Coarsen::coarsen(hier::Patch& coarse,
                                                   const hier::Patch& fine,
                                                   const int dst_component,
                                                   const int src_component,
                                                   const hier::Box& coarse_box,
                                                   const hier::IntVector& ratio)
  const
{
  const tbox::Dimension& dimension(getDim());
  const int dim(dimension.getValue());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dimension, coarse, fine, coarse_box, ratio);

  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_fine_ptr = fine.getPatchData(src_component);
  pdat::CellData<double> &cell_viscosity_fine(*cell_viscosity_fine_ptr);
  tbox::Pointer<pdat::NodeData<double> > edge_viscosity_fine_2D_ptr;
  tbox::Pointer<pdat::EdgeData<double> > edge_viscosity_fine_3D_ptr;
  if(dim==2)
    edge_viscosity_fine_2D_ptr = fine.getPatchData(edge_viscosity_id);
  else
    edge_viscosity_fine_3D_ptr = fine.getPatchData(edge_viscosity_id);

  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_coarse_ptr = coarse.getPatchData(dst_component);
  pdat::CellData<double> &cell_viscosity_coarse(*cell_viscosity_coarse_ptr);

  TBOX_ASSERT(!cell_viscosity_coarse_ptr.isNull());
  TBOX_ASSERT(!cell_viscosity_fine_ptr.isNull());
  TBOX_ASSERT(cell_viscosity_fine.getDepth() == cell_viscosity_coarse.getDepth());
  TBOX_ASSERT(cell_viscosity_coarse.getDepth() == 1);

  const tbox::Pointer<CartesianPatchGeometry> cgeom =
    coarse.getPatchGeometry();

  hier::Index ip(hier::Index::getZeroIndex(dimension)), jp(ip), kp(ip);
  ip[0]=1;
  jp[1]=1;
  if(dim>2)
    kp[2]=1;

  for(pdat::CellIterator ci(coarse_box); ci; ci++)
    {
      pdat::CellIndex coarse_cell(*ci);

      if(dim==2)
        {
          cell_viscosity_coarse(coarse_cell)=
            viscosity_coarsen_2D(cell_viscosity_fine,*edge_viscosity_fine_2D_ptr,
                                 coarse_cell*2+ip+jp);
        }
      else
        {
          cell_viscosity_coarse(coarse_cell)=
            viscosity_coarsen_3D(cell_viscosity_fine,*edge_viscosity_fine_3D_ptr,
                                 coarse_cell*2);
        }
    }
}

#endif
