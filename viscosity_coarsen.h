#ifndef GAMR_VISCOSITY_COARSEN_H
#define GAMR_VISCOSITY_COARSEN_H

/* This uses both the cell-centered and node-centered viscosities to
   coarsen the viscosity.  In 2D, if you rotate the grid 45 degrees,
   then it becomes a regular node-based grid.  So we use a standard
   full-weighted coarsening.  There is no refining needed, since
   this is only for viscosity.

   Note that coarsening for the node and for the cell are the same
   except for an offset.

   Ee------e-------Ee------e-------Ee
   |       |       |       |       |
   |   c   |   c   |   c   |   c   |
   |       |       |       |       |
   e-------Ce------e-------Ce------e
   |       |       |       |       |
   |   c   |   c   |   c   |   c   |
   |       |       |       |       |
   Ee------e-------Ee------e-------Ee
   |       |       |       |       |
   |   c   |   c   |   c   |   c   |
   |       |       |       |       |
   e-------Ce------e-------Ce------e
   |       |       |       |       |
   |   c   |   c   |   c   |   c   |
   |       |       |       |       |
   Ee------e-------Ee------e-------Ee

*/

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Index.h"

inline double viscosity_coarsen(const SAMRAI::pdat::CellData<double> &cell,
                                const SAMRAI::pdat::NodeData<double> &edge,
                                const SAMRAI::hier::Index &fine)
{
  SAMRAI::hier::Index ip(1,0), jp(0,1);
  SAMRAI::pdat::CellIndex fine_cell(fine);
  SAMRAI::pdat::NodeIndex fine_edge(fine,SAMRAI::pdat::NodeIndex::LowerLeft);
  return (cell(fine_cell) + cell(fine_cell-ip) + cell(fine_cell-jp)
          + cell(fine_cell-ip-jp))/8
    + edge(fine_edge)/4
    + (edge(fine_edge+ip) + edge(fine_edge+jp) + edge(fine_edge-ip)
       + edge(fine_edge-jp))/16;
}

#endif
