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
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/hier/Index.h"

inline double viscosity_coarsen_2D(const SAMRAI::pdat::CellData<double> &cell,
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

inline double viscosity_coarsen_3D(const SAMRAI::pdat::CellData<double> &cell,
                                   const SAMRAI::pdat::EdgeData<double> &edge,
                                   const SAMRAI::hier::Index &fine)
{
  SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
  SAMRAI::pdat::CellIndex fine_cell(fine);
  SAMRAI::pdat::EdgeIndex
    fine_edge_x(fine_cell,SAMRAI::pdat::EdgeIndex::X,
                SAMRAI::pdat::EdgeIndex::LowerLeft),
    fine_edge_y(fine_cell,SAMRAI::pdat::EdgeIndex::Y,
                SAMRAI::pdat::EdgeIndex::LowerLeft),
    fine_edge_z(fine_cell,SAMRAI::pdat::EdgeIndex::Z,
                SAMRAI::pdat::EdgeIndex::LowerLeft);
  return (cell(fine_cell)
          + cell(fine_cell+ip)
          + cell(fine_cell+jp)
          + cell(fine_cell+kp)
          + cell(fine_cell+ip+jp)
          + cell(fine_cell+ip+kp)
          + cell(fine_cell+jp+kp)
          + cell(fine_cell+ip+jp+kp))/32
    + (edge(fine_edge_x+jp+kp) + edge(fine_edge_x+ip+jp+kp)
       + edge(fine_edge_y+ip+kp) + edge(fine_edge_y+ip+jp+kp)
       + edge(fine_edge_z+ip+jp) + edge(fine_edge_z+ip+jp+kp))/8;
}

#endif
