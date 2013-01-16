/* Helper templates that make it easy to refer to node/edge data
   without having to duplicate code for 2D and 3D. */

#ifndef GAMRA_EDGE_NODE_EVAL_H
#define GAMRA_EDGE_NODE_EVAL_H

#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/SideIndex.h"


inline double &edge_node_eval(SAMRAI::pdat::EdgeData<double> &edge,
                              const SAMRAI::pdat::SideIndex &s,
                              const int &ix,
                              const int &depth)
{
  SAMRAI::pdat::EdgeIndex e(s,ix,SAMRAI::pdat::EdgeIndex::LowerLeft);
  return edge(e,depth);  
}

inline double &edge_node_eval(SAMRAI::pdat::NodeData<double> &edge,
                              const SAMRAI::pdat::SideIndex &s,
                              const int &,
                              const int &depth)
{
  SAMRAI::pdat::NodeIndex e(s,SAMRAI::pdat::NodeIndex::LowerLeft);
  return edge(e,depth);  
}

inline double edge_node_eval(const SAMRAI::pdat::EdgeData<double> &edge,
                             const SAMRAI::pdat::SideIndex &s,
                             const int &ix,
                             const int &depth)
{
  SAMRAI::pdat::EdgeIndex e(s,ix,SAMRAI::pdat::EdgeIndex::LowerLeft);
  return edge(e,depth);  
}

inline double edge_node_eval(const SAMRAI::pdat::NodeData<double> &edge,
                             const SAMRAI::pdat::SideIndex &s,
                             const int &,
                             const int &depth)
{
  SAMRAI::pdat::NodeIndex e(s,SAMRAI::pdat::NodeIndex::LowerLeft);
  return edge(e,depth);  
}

inline double edge_node_average(const SAMRAI::pdat::NodeData<double> &edge,
                                const SAMRAI::pdat::SideIndex &s, const int &ix,
                                const SAMRAI::hier::Index pp[],
                                const int &depth)
{
  const int dim(2);
  const int iy((ix+1)%dim);
  SAMRAI::pdat::NodeIndex n(s,SAMRAI::pdat::NodeIndex::LowerLeft);
  return (edge(n,depth) + edge(n+pp[iy],depth))/2;
}

inline double edge_node_average(const SAMRAI::pdat::EdgeData<double> &edge,
                                const SAMRAI::pdat::SideIndex &s, const int &ix,
                                const SAMRAI::hier::Index pp[],
                                const int &depth)
{
  const int dim(3);
  const int iy((ix+1)%dim);
  const int iz((iy+1)%dim);
  SAMRAI::pdat::EdgeIndex ey(s,iy,SAMRAI::pdat::EdgeIndex::LowerLeft),
    ez(s,iz,SAMRAI::pdat::EdgeIndex::LowerLeft);
  return (edge(ey,depth) + edge(ey+pp[iy],depth)
          + edge(ez,depth) + edge(ez+pp[iz],depth))/4;
}

#endif
