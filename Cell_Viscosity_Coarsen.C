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
  const tbox::Dimension& dim(getDim());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_fine = fine.getPatchData(src_component);
  tbox::Pointer<pdat::NodeData<double> >
    edge_viscosity_fine = fine.getPatchData(edge_viscosity_id);
  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_coarse = coarse.getPatchData(dst_component);

  TBOX_ASSERT(!cell_viscosity_coarse.isNull());
  TBOX_ASSERT(!cell_viscosity_fine.isNull());
  TBOX_ASSERT(cell_viscosity_fine->getDepth() == cell_viscosity_coarse->getDepth());
  TBOX_ASSERT(cell_viscosity_coarse->getDepth() == 1);

  const tbox::Pointer<CartesianPatchGeometry> cgeom =
    coarse.getPatchGeometry();

  hier::Index ip(1,0), jp(0,1);
   for(int j=coarse_box.lower(1); j<=coarse_box.upper(1); ++j)
     for(int i=coarse_box.lower(0); i<=coarse_box.upper(0); ++i)
       {
         pdat::CellIndex coarse_cell(hier::Index(i,j));

         (*cell_viscosity_coarse)(coarse_cell)=
           viscosity_coarsen(*cell_viscosity_fine,*edge_viscosity_fine,
                             coarse_cell*2+ip+jp);

         tbox::plog << "Cell "
                    << coarse_box.lower(0) << " "
                    << coarse_box.upper(0) << " "
                    << coarse_box.lower(1) << " "
                    << coarse_box.upper(1) << " "
                    << i << " "
                    << j << " "
                    << (*cell_viscosity_coarse)(coarse_cell) << " "
                    << "\n";

       }
}

#endif
