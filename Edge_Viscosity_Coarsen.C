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
  const tbox::Dimension& dim(getDim());

  TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

  tbox::Pointer<pdat::CellData<double> >
    cell_viscosity_fine = fine.getPatchData(cell_viscosity_id);
  tbox::Pointer<pdat::NodeData<double> >
    edge_viscosity_fine = fine.getPatchData(src_component);
  tbox::Pointer<pdat::NodeData<double> >
    edge_viscosity_coarse = coarse.getPatchData(dst_component);

  TBOX_ASSERT(!edge_viscosity_coarse.isNull());
  TBOX_ASSERT(!edge_viscosity_fine.isNull());
  TBOX_ASSERT(edge_viscosity_fine->getDepth() == edge_viscosity_coarse->getDepth());
  TBOX_ASSERT(edge_viscosity_coarse->getDepth() == 1);

  const tbox::Pointer<CartesianPatchGeometry> cgeom =
    coarse.getPatchGeometry();
  hier::Index ip(1,0), jp(0,1);
   for(int j=coarse_box.lower(1); j<=coarse_box.upper(1)+1; ++j)
     for(int i=coarse_box.lower(0); i<=coarse_box.upper(0)+1; ++i)
       {
         pdat::NodeIndex coarse_edge(hier::Index(i,j),
                                     pdat::NodeIndex::LowerLeft);

         (*edge_viscosity_coarse)(coarse_edge)=
           viscosity_coarsen(*cell_viscosity_fine,*edge_viscosity_fine,
                             coarse_edge*2);
       }
}

#endif
