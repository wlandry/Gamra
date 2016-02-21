/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Numerical routines for example FAC Elastic solver 
 *
 ************************************************************************/
#include "Elastic/FAC.hxx"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/solv/SimpleCellRobinBcCoefs.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/math/HierarchyCellDataOpsReal.h>
#include <SAMRAI/pdat/SideData.h>
#include <SAMRAI/tbox/Utilities.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/hier/VariableDatabase.h>

/*
*************************************************************************
* Write derived data to the given stream.                               *
*************************************************************************
*/
void
Elastic::FAC::pack_level_set(double* buffer,
                             const SAMRAI::hier::Patch& patch,
                             const SAMRAI::hier::Box& region) const
{
  boost::shared_ptr<SAMRAI::pdat::SideData<double> > level_set_ptr;
  level_set_ptr = boost::dynamic_pointer_cast<SAMRAI::pdat::SideData<double> >
    (patch.getPatchData(level_set_id));

  SAMRAI::pdat::SideData<double>& level_set = *level_set_ptr;
  if(d_dim.getValue()==2)
    {
      const SAMRAI::hier::Index ip(1,0), jp(0,1);

      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell) {

        const SAMRAI::pdat::CellIndex &center(*icell);
        const SAMRAI::pdat::SideIndex
          x(center,0,SAMRAI::pdat::SideIndex::Lower),
          y(center,1,SAMRAI::pdat::SideIndex::Lower);
	
        double vx=(level_set(x+ip) + level_set(x))/2.;
        double vy=(level_set(y+jp) + level_set(y))/2.;

        *buffer = (vx+vy)/2;
        buffer = buffer + 1;
      }
    }
  else
    {
      const SAMRAI::hier::Index ip(1,0,0), jp(0,1,0), kp(0,0,1);
      SAMRAI::pdat::CellIterator iend(SAMRAI::pdat::CellGeometry::end(region));
      for (SAMRAI::pdat::CellIterator
             icell(SAMRAI::pdat::CellGeometry::begin(region));
           icell!=iend; ++icell) {

        const SAMRAI::pdat::CellIndex &center(*icell);
        const SAMRAI::pdat::SideIndex
          x(center,0,SAMRAI::pdat::SideIndex::Lower),
          y(center,1,SAMRAI::pdat::SideIndex::Lower),
          z(center,2,SAMRAI::pdat::SideIndex::Lower);
	
        double vx=(level_set(x+ip) + level_set(x))/2.;
        double vy=(level_set(y+jp) + level_set(y))/2.;
        double vz=(level_set(z+kp) + level_set(z))/2.;

        *buffer = (vx+vy+vz)/3;
        buffer = buffer + 1;
      }
    }
}